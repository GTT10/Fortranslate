#!/usr/bin/env python3
"""Compare uniform LiDryer reacting planar-EB trajectories with pinned PeleC.

This is a single-level application comparison, not spatial-gradient, AMR,
transport or physical ignition validation. Acceptance is fixed in the case.
"""
import argparse
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


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('pelec', 'fortran', 'exporter', 'reference-source', 'work'):
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--reuse-run', type=Path,
                        help='Recompare completed simulation outputs after verifying input and executable hashes')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    case = root / 'cases/pelec_eb_reacting'
    contract = json.loads((case / 'acceptance.json').read_text())
    tolerances = contract['acceptance']
    bundle = root / 'mechanisms/lidryer.json'
    chemistry = json.loads(bundle.read_text())
    assert digest(root / 'mechanisms/lidryer.yaml') == chemistry['source']['sha256']
    provenance = json.loads((root / 'references/pelec_eb_reacting.json').read_text())
    baseline = json.loads((root / 'references/pelec_baseline.json').read_text())
    assert provenance['pelec_commit'] == baseline['pelec']['commit']
    assert json.loads((args.reference_source / 'pelef-reference-provenance.json').read_text()) == baseline
    for path, expected in provenance['source_sha256'].items():
        assert digest(args.reference_source / path) == expected, ('reference source drift', path)
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='matched_eb_', dir=args.work.resolve()))
    print(work, flush=True)
    report = {'scope': __doc__, 'status': 'running', 'contract': contract,
              'reference_commit': provenance['pelec_commit'],
              'bundle_sha256': digest(bundle), 'source_contract_sha256': digest(root / 'references/pelec_eb_reacting.json'),
              'checker_sha256': digest(Path(__file__)), 'commands': [], 'steps': []}
    env = os.environ.copy()
    env.update(HWLOC_COMPONENTS='-gl,-x86', OMP_NUM_THREADS='1')

    def save():
        (work / 'checks.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')

    def run(binary, arguments, directory, log):
        binary = binary.resolve(strict=True)
        command = [str(binary), *map(str, arguments)]
        entry = {'command': command, 'cwd': str(directory), 'binary_sha256': digest(binary)}
        report['commands'].append(entry)
        save()
        start = time.monotonic()
        try:
            with (directory / log).open('w') as stream:
                result = subprocess.run(command, cwd=directory, env=env, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=120)
            entry['returncode'] = result.returncode
            assert result.returncode == 0, directory / log
        finally:
            entry['seconds'] = time.monotonic() - start
            save()

    try:
        previous = None
        if args.reuse_run:
            previous = json.loads((args.reuse_run / 'checks.json').read_text())
            for key in ('bundle_sha256', 'source_contract_sha256', 'contract'):
                assert previous[key] == report[key], ('reuse contract mismatch', key)
        for name, source, binary in (('reference', 'pelec.inp', args.pelec), ('fortran', 'pelef.nml', args.fortran)):
            directory = work / name
            if previous is not None:
                original = args.reuse_run / name
                assert digest(original / source) == digest(case / source) == previous[name + '_input_sha256']
                command = [str(binary.resolve(strict=True)), source]
                matched = [item for item in previous['commands'] if item['command'] == command]
                assert len(matched) == 1 and matched[0]['returncode'] == 0
                assert matched[0]['binary_sha256'] == digest(binary)
                shutil.copytree(original, directory)
                report['commands'].append({**matched[0], 'reused_from': str(original.resolve())})
            else:
                directory.mkdir()
                shutil.copyfile(case / source, directory / source)
            report[name + '_input_sha256'] = digest(directory / source)
            if previous is None:
                run(binary, [source], directory, 'run.log')
        names = [s['name'] for s in chemistry['thermo']]
        molecular_weights = {s['name']: s['molecular_weight'] for s in chemistry['thermo']}
        specific_gas_constant = 8314.46261815324 * sum(
            fraction / molecular_weights[name] for name, fraction in contract['mass_fractions'].items())
        initial_density = contract['pressure_Pa'] / (specific_gas_constant * contract['temperature_K'])
        fields = [('density', 'rho', 1000.), ('xmom', 'rhou', 10.), ('ymom', 'rhov', 10.),
                  ('zmom', 'rhow', 10.), ('rho_E', 'rhoE', .1),
                  ('rho_e', 'rho_internal_energy_J_m3', .1), ('Temp', 'temperature', 1.)]
        fields += [('rho_' + name, 'rhoY_' + name, 1000.) for name in names]
        initial_h2 = None
        final_reference_h2 = final_fortran_h2 = None
        for step in range(contract['steps'] + 1):
            exported = work / 'reference' / f'cells_{step}.csv'
            run(args.exporter, [work / 'reference' / f'plt{step:05}', exported, '--levels'],
                work / 'reference', f'export_{step}.log')
            reference = {}
            volumes = {}
            target_time = step * contract['dt_s']
            for row in rows(exported):
                assert abs(float(row['time']) - target_time) <= 1.e-12 * contract['dt_s']
                assert int(row['level']) == 0
                assert math.isfinite(float(row['vfrac']))
                if float(row['vfrac']) == 0:
                    continue
                assert all(math.isfinite(float(value)) for value in row.values())
                assert 0 < float(row['vfrac']) <= 1
                key = tuple(round(float(row[axis]) / float(row['d' + axis]) - .5) + 1 for axis in 'xyz')
                for index, axis in enumerate('xyz'):
                    expected_dx = contract['domain_m'][index] / contract['root_shape'][index]
                    assert abs(float(row['d' + axis]) * .01 / expected_dx - 1) <= 1.e-12
                    assert 1 <= key[index] <= contract['root_shape'][index]
                dz = contract['domain_m'][2] / contract['root_shape'][2]
                expected_fraction = min(1., max(0., (key[2] * dz - contract['plane_position_m']) / dz))
                assert abs(float(row['vfrac']) - expected_fraction) <= 1.e-12
                assert key not in reference
                reference[key] = {f: float(row[r]) * scale for r, f, scale in fields}
                volumes[key] = float(row['vfrac']) * math.prod(float(row['d' + axis]) * .01 for axis in 'xyz')
            filename = f'result.csv.step{step}.csv'
            stamps = rows(work / 'fortran' / (filename + '.time.csv'))
            assert len(stamps) == 1 and int(stamps[0]['step']) == step
            assert abs(float(stamps[0]['time_s']) - target_time) <= 1.e-12 * contract['dt_s']
            assert int(stamps[0]['nodes']) == 1 and int(stamps[0]['regrid_stage']) == 0
            actual = {}
            for row in rows(work / 'fortran' / filename):
                assert all(math.isfinite(float(value)) for value in row.values())
                assert int(row['level']) == 0 and int(row['node']) == 1
                key = tuple(int(row[axis]) for axis in 'ijk')
                assert key not in actual
                actual[key] = {name: float(value) for name, value in row.items()}
            nx, ny, nz = contract['root_shape']
            covered_layers = math.floor(contract['plane_position_m'] / (contract['domain_m'][2] / nz))
            assert set(actual) == set(reference) and len(actual) == nx * ny * (nz - covered_layers), \
                'matching exposed cells required'
            for key in actual:
                assert abs(actual[key]['fluid_volume'] / volumes[key] - 1) <= 1.e-12
                assert all(actual[key]['rhoY_' + name] >= 0 for name in names)
                assert abs(sum(actual[key]['rhoY_' + name] for name in names) / actual[key]['rho'] - 1) <= 1.e-12
                if step == 0:
                    assert abs(actual[key]['rho'] / initial_density - 1) <= tolerances['initial_scaled']
                    assert abs(actual[key]['temperature'] / contract['temperature_K'] - 1) <= tolerances['initial_scaled']
                    for name in names:
                        assert abs(actual[key]['rhoY_' + name] / actual[key]['rho'] -
                                   contract['mass_fractions'].get(name, 0.)) <= tolerances['initial_scaled']
            errors = {}
            report['steps'].append({'step': step, 'time_s': target_time, 'cells': len(actual), 'errors': errors})
            for _, field, _ in fields:
                absolute = max(abs(actual[key][field] - reference[key][field]) for key in actual)
                scale = max(1.e-30, max(abs(reference[key][field]) for key in actual))
                errors[field] = {'absolute': absolute, 'relative_to_field_max': absolute / scale}
                if step == 0:
                    assert absolute <= tolerances['initial_scaled'] * max(1., scale), ('initial', field, absolute)
                if field.startswith('rhoY_'):
                    y_error = max(abs(actual[key][field] / actual[key]['rho'] -
                                      reference[key][field] / reference[key]['rho']) for key in actual)
                    errors[field]['mass_fraction_absolute'] = y_error
                    assert y_error <= tolerances['species_mass_fraction_absolute'], (step, field, y_error)
                elif field == 'temperature':
                    assert absolute <= tolerances['temperature_absolute_K'], (step, field, absolute)
                elif field in ('rhou', 'rhov', 'rhow'):
                    assert absolute <= tolerances['momentum_absolute_kg_m2_s'], (step, field, absolute)
                else:
                    limit = tolerances['density_scaled'] if field == 'rho' else tolerances['rhoE_scaled']
                    assert absolute / scale <= limit, (step, field, absolute / scale)
            # Independently integrate every conserved component and account for
            # chemistry, walls, boundaries and numerical/forcing contributions.
            budget = rows(work / 'fortran' / (filename + '.budget.csv'))
            conservative = ['rho', 'rhou', 'rhov', 'rhow', 'rhoE'] + ['rhoY_' + name for name in names]
            assert len(budget) == len(conservative)
            for field, entry in zip(conservative, budget):
                inventory = math.fsum(actual[key][field] * volumes[key] for key in actual)
                final, initial = float(entry['final']), float(entry['initial'])
                changes = [float(entry[label]) for label in ('boundary', 'wall', 'chemistry', 'forcing', 'numerical_energy')]
                scale = max(abs(initial), abs(final), math.fsum(abs(v) for v in changes), 1.e-30)
                assert abs(inventory - final) <= 1.e-12 * scale, (step, field, 'CSV inventory')
                assert abs(final - initial - math.fsum(changes)) <= 5.e-10 * scale, (step, field, 'budget')
            final_reference_h2 = math.fsum(reference[key]['rhoY_H2'] * volumes[key] for key in actual)
            final_fortran_h2 = math.fsum(actual[key]['rhoY_H2'] * volumes[key] for key in actual)
            if step == 0:
                initial_h2 = final_reference_h2
            save()
        for name, final in (('reference', final_reference_h2), ('fortran', final_fortran_h2)):
            fraction = (initial_h2 - final) / initial_h2
            assert fraction > tolerances['minimum_H2_consumption_fraction'], ('resolved reaction required', name, fraction)
            report[name + '_hydrogen_consumed_fraction'] = fraction
        consumption_error = abs(final_reference_h2 - final_fortran_h2) / (initial_h2 - final_reference_h2)
        report['hydrogen_consumption_relative_error'] = consumption_error
        assert consumption_error <= tolerances['H2_consumption_relative'], ('chemical change parity', consumption_error)
        report['status'] = 'passed'
        print('PASS: matched LiDryer reacting planar-EB trajectories at all four physical times', flush=True)
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = repr(error)
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
