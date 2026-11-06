#!/usr/bin/env python3
"""Bounded EB MOL time refinement against an identified native held-EOS diagnostic.

Uses a completed frozen-case capture as initial data. This study cannot qualify
the unmodified pinned algorithm. It leaves the frozen reference checker intact.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline', 'policy-record', 'probe', 'native-held', 'exporter', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    baseline = args.baseline.resolve()
    policy = args.policy_record.resolve()
    probe, native, exporter = (p.resolve() for p in (args.probe, args.native_held, args.exporter))
    contract = json.loads((root / 'cases/tgreact_eb_amr_reference/acceptance.json').read_text())
    assert contract['root_shape'] == [8, 8, 8]
    assert contract['fine_coarse_box'] == [[1, 1, 1], [4, 4, 4]]
    chemistry = json.loads((root / 'mechanisms/lidryer.json').read_text())
    seed = json.loads((baseline / 'checks.json').read_text())
    provenance = json.loads((policy / 'identity.json').read_text())
    policy_check = json.loads((policy / 'verification.json').read_text())
    assert policy_check['pinned_reference_qualified'] is False
    assert policy_check['diagnostic_field_limit_exceedances'] == 0
    assert digest(probe) == seed['commands'][-1]['binary_sha256'], 'baseline executable changed'
    assert digest(baseline / 'initial.dat') == seed['import_sha256'], 'initial data changed'
    assert digest(root / 'mechanisms/lidryer.json') == seed['bundle_sha256']
    assert digest(root / 'cases/tgreact_eb_amr_reference/acceptance.json') == seed['case_contract_sha256']
    for path in (native, exporter):
        assert digest(path) == provenance[str(path)], ('diagnostic executable changed', path)
    for path, sha in provenance.items():
        assert digest(Path(path)) == sha, ('policy source identity changed', path)
    case = root / 'cases/tgreact_eb_amr_reference'
    for name, sha in contract['input_sha256'].items():
        assert digest(case / name) == sha, ('frozen case changed', name)

    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='temporal_', dir=args.work.resolve()))
    shutil.copyfile(__file__, work / 'checker.py')
    print(work, flush=True)
    names = chemistry['species']
    fields = ['density', 'xmom', 'ymom', 'zmom', 'rho_E'] + ['rho_' + n for n in names] + ['rho_e', 'Temp']
    factors = [1000., 10., 10., 10., .1] + [1000.] * len(names) + [.1, 1.]
    spacing = [length / n for length, n in zip(contract['domain_m'], contract['root_shape'])]
    limits = contract['acceptance']
    keys = {(level, i, j, k) for level in (0, 1) for i in range(1, 9)
            for j in range(1, 9) for k in range(1, 9)}

    def fraction(key):
        dz = spacing[2] / 2**key[0]
        return min(1., max(0., (key[3] * dz - contract['plane_z_m']) / dz))

    cells = sorted(k for k in keys if fraction(k) > 0 and (k[0] == 1 or any(v > 4 for v in k[1:])))
    volumes = {k: math.prod(spacing) / 8**k[0] * fraction(k) for k in cells}
    assert len(cells) == 544

    def load(path, is_native, physical_time):
        data = {}
        for row in rows(path):
            assert all(math.isfinite(float(v)) for v in row.values()), ('nonfinite export', path)
            assert abs(float(row['time' if is_native else 'time_s']) - physical_time) < 1.e-18
            level = int(row['level'])
            index = tuple(round(float(row[a]) / float(row['d' + a]) - .5) + 1 for a in 'xyz') if is_native else tuple(int(row[a]) for a in 'ijk')
            key = (level, *index)
            assert key in keys and key not in data
            assert abs(float(row['vfrac']) - fraction(key)) < 1.e-13
            if is_native:
                assert all(abs(float(row['d' + a]) * .01 / (spacing[d] / 2**level) - 1) < 1.e-12 for d, a in enumerate('xyz'))
            else:
                assert int(row['node']) == level + 1
            data[key] = {f: float(row[f]) * scale for f, scale in zip(fields, factors if is_native else [1.] * len(fields))}
        assert data.keys() == keys
        return data

    def inventory(data):
        inv = {f: math.fsum(data[k][f] * volumes[k] for k in cells) for f in fields[:-2]}
        for element in ('H', 'O', 'N'):
            inv['element_' + element] = math.fsum(inv['rho_' + s['name']] * s['composition'].get(element, 0) / s['molecular_weight'] for s in chemistry['thermo'])
        return inv

    initial = load(baseline / 'fortran_000.csv', False, 0.)
    initial_inventory = inventory(initial)
    env = {**os.environ, 'OMP_NUM_THREADS': '1', 'HWLOC_COMPONENTS': '-gl,-x86'}
    schedules = [('half', 5.e-8, 6), ('quarter', 2.5e-8, 12), ('long', 1.e-7, 12)]
    report = {'scope': 'temporal diagnostic with native held spatial EOS; fixed mesh',
              'qualification': False, 'pinned_reference_status': seed['status'],
              'pinned_reference_violations': seed['violations'], 'status': 'running',
              'baseline': str(baseline), 'policy_record': str(policy),
              'policy_identity': provenance,
              'executables': {str(p): digest(p) for p in (probe, native, exporter)},
              'checker_sha256': digest(Path(__file__)), 'initial_sha256': seed['import_sha256'],
              'case_contract_sha256': seed['case_contract_sha256'], 'schedules': schedules}
    save(work / 'checks.json', report)

    def run_case(schedule):
        name, dt, steps = schedule
        directory = work / name
        directory.mkdir()
        native_dir = directory / 'native'
        native_dir.mkdir()
        result = {'dt_s': dt, 'steps_requested': steps, 'status': 'running', 'commands': [], 'snapshots': [], 'diagnostic_limit_exceedances': []}

        def command(argv, cwd, log_name):
            start = time.monotonic()
            with (directory / log_name).open('w') as log:
                p = subprocess.run([str(a) for a in argv], cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT)
            result['commands'].append({'command': [str(a) for a in argv], 'cwd': str(cwd), 'returncode': p.returncode, 'seconds': time.monotonic() - start, 'log': log_name})
            save(directory / 'checks.json', result)
            assert p.returncode == 0, ('command failed', log_name, p.returncode)

        try:
            for filename in contract['input_sha256']:
                shutil.copyfile(case / filename, native_dir / filename)
            command([native, 'inputs', f'max_step={steps}', f'stop_time={steps*dt:.17e}', f'pelec.fixed_dt={dt:.17e}'], native_dir, 'native.log')
            for step in range(steps + 1):
                command([exporter, native_dir / f'plt{step:05}', directory / f'native_{step:03}.csv', '--levels'], directory, f'export-{step}.log')
            lines = (baseline / 'initial.dat').read_text().splitlines()
            lines[1] = f'{steps} 2 1 {dt:.17e}'
            (directory / 'initial.dat').write_text('\n'.join(lines) + '\n')
            result['input_sha256'] = digest(directory / 'initial.dat')
            command([probe, directory / 'initial.dat', directory / 'fortran'], directory, 'fortran.log')
            for step in range(steps + 1):
                f = load(directory / f'fortran_{step:03}.csv', False, step * dt)
                n = load(directory / f'native_{step:03}.csv', True, step * dt)
                if step == 0:
                    assert f == initial and n == initial, 'initial gas/UEINT/T must match exactly'
                snap = {'step': step, 'time_s': step * dt, 'fields': {}}
                for field in fields:
                    absolute = max(abs(f[k][field] - n[k][field]) for k in cells)
                    scaled = absolute / max(1.e-30, max(abs(n[k][field]) for k in cells))
                    item = {'absolute': absolute, 'scaled': scaled}
                    snap['fields'][field] = item
                    if field.startswith('rho_') and field not in ('rho_E', 'rho_e'):
                        value = max(abs(f[k][field] / f[k]['density'] - n[k][field] / n[k]['density']) for k in cells)
                        item['mass_fraction_absolute'] = value
                        limit = limits['species_mass_fraction_absolute']
                    elif field == 'Temp':
                        value, limit = absolute, limits['temperature_absolute_K']
                    elif field in ('xmom', 'ymom', 'zmom'):
                        value, limit = absolute, limits['momentum_absolute_kg_m2_s']
                    else:
                        value, limit = scaled, limits['density_scaled' if field == 'density' else 'energy_scaled']
                    if value > limit:
                        result['diagnostic_limit_exceedances'].append({'step': step, 'field': field, 'value': value, 'limit': limit})
                for label, data in (('fortran', f), ('native_held', n)):
                    inv = inventory(data)
                    snap[label + '_inventory'] = inv
                    assert all(data[k]['density'] > 0 and all(data[k]['rho_' + name] >= 0 for name in names) for k in cells)
                    closure = max(abs(sum(data[k]['rho_' + name] for name in names) / data[k]['density'] - 1) for k in cells)
                    drift = max(abs(inv[field] / initial_inventory[field] - 1) for field in ('density', 'element_H', 'element_O', 'element_N'))
                    snap[label + '_inventory_drift'] = drift
                    assert max(closure, drift) <= limits['inventory_relative'], (label, step, closure, drift)
                budget = rows(directory / f'fortran_{step:03}.budget.csv')
                assert len(budget) == 1 and int(budget[0]['step']) == step
                assert abs(float(budget[0]['time_s']) - step * dt) < 1.e-18
                residual = (snap['fortran_inventory']['rho_E'] - initial_inventory['rho_E'] - float(budget[0]['numerical_energy_J'])) / abs(initial_inventory['rho_E'])
                snap['energy_budget_relative_residual'] = residual
                assert abs(residual) <= limits['energy_budget_relative'], (step, residual)
                result['snapshots'].append(snap)
                save(directory / 'checks.json', result)
            consumed = {label: 1 - result['snapshots'][-1][label + '_inventory']['rho_H2'] / initial_inventory['rho_H2'] for label in ('fortran', 'native_held')}
            assert all(x >= limits['minimum_H2_consumption_fraction'] for x in consumed.values())
            result['H2_consumed_fraction'] = consumed
            result['H2_consumption_relative_error'] = abs(consumed['fortran'] / consumed['native_held'] - 1)
            if result['H2_consumption_relative_error'] > limits['H2_consumption_relative']:
                result['diagnostic_limit_exceedances'].append({
                    'step': steps, 'field': 'H2_consumption_relative_error',
                    'value': result['H2_consumption_relative_error'],
                    'limit': limits['H2_consumption_relative']})
            result['status'] = 'completed'
        except Exception as error:
            result['status'] = 'failed'
            result['error'] = str(error)
        save(directory / 'checks.json', result)
        return name, result

    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(run_case, s) for s in schedules]
        for future in concurrent.futures.as_completed(futures):
            name, result = future.result()
            results[name] = result
            report['runs'] = results
            save(work / 'checks.json', report)
            print(name, result['status'], flush=True)
    if all(r['status'] == 'completed' for r in results.values()):
        states = [load(baseline / 'fortran_003.csv', False, 3.e-7),
                  load(work / 'half/fortran_006.csv', False, 3.e-7),
                  load(work / 'quarter/fortran_012.csv', False, 3.e-7)]
        convergence = {}
        volume = math.fsum(volumes.values())
        for field in fields:
            norms = [math.sqrt(math.fsum((a[k][field] - b[k][field])**2 * volumes[k] for k in cells) / volume) for a, b in zip(states, states[1:])]
            scale = max(abs(s[k][field]) for s in states for k in cells)
            resolved = min(norms) > 128 * math.ulp(max(scale, 1.e-30))
            convergence[field] = {'volume_weighted_rms_dt_vs_half': norms[0], 'volume_weighted_rms_half_vs_quarter': norms[1], 'observed_order': math.log2(norms[0] / norms[1]) if resolved else None}
        report['temporal_comparison_at_300ns'] = convergence
        report['status'] = 'completed'
    else:
        report['status'] = 'failed'
    save(work / 'checks.json', report)
    if report['status'] != 'completed':
        raise SystemExit('Temporal study incomplete; inspect saved per-run results')


if __name__ == '__main__':
    main()
