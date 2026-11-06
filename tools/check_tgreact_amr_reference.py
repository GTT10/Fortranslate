#!/usr/bin/env python3
"""Compare the production LiDryer AMR MOL integrator with unchanged pinned TGReact.

Matched two-level reacting control with periodic, adiabatic or isothermal x walls,
optional Simple transport and serial/MPI probes, no EB. This does not qualify
the full application, regridding, convergence or experiment.
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
from pelephysics_transport_import import import_simple_transport


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('probe', 'exporter', 'work'):
        parser.add_argument('--' + name, required=True, type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--pelec', type=Path)
    mode.add_argument('--reference', type=Path, help='Captured native reference with its SHA manifest')
    parser.add_argument('--reference-source', type=Path)
    parser.add_argument('--x-boundary', choices=('periodic', 'slip_wall', 'no_slip_wall', 'outflow'), default='periodic')
    parser.add_argument('--wall-temperature', type=float,
                        help='Positive common x-wall temperature in K; requires a wall and transport')
    parser.add_argument('--mpi-ranks', type=int, choices=(0, 1, 2, 3), default=0)
    parser.add_argument('--mpiexec', type=Path, default=Path('/usr/bin/mpiexec'))
    parser.add_argument('--reaction-off', action='store_true', help='Matched no-reaction diagnostic with zero initial history')
    parser.add_argument('--initial-velocity', choices=('tg', 'zero'), default='tg',
                        help='Explicit initial-flow control; zero preserves the frozen thermal/composition profile')
    parser.add_argument('--mol-iterations', type=int, choices=(1, 2), default=2)
    parser.add_argument('--cells', type=int, choices=(4, 8, 16, 32, 64), default=4,
                        help='Root cells per direction; preserve the half-domain refinement region')
    parser.add_argument('--transport', action='store_true',
                        help='Enable all molecular transport with pinned Simple fits; requires --reference-source')
    parser.add_argument('--subcycling', type=int, choices=(0, 1), default=1,
                        help='0: common level timestep; 1: recursive level timesteps (default)')
    parser.add_argument('--initialize-history', action='store_true',
                        help='Independently compute initial reaction history with the production initializer')
    parser.add_argument('--time-refinement', type=int, choices=(1, 2, 4), default=1,
                        help='Divide the frozen timestep and multiply steps; keep final time fixed')
    parser.add_argument('--command-timeout', type=float, default=120.,
                        help='Wall-time limit in seconds per subprocess (default: 120)')
    args = parser.parse_args()
    if not math.isfinite(args.command_timeout) or args.command_timeout <= 0:
        parser.error('--command-timeout must be finite and positive')
    assert not (args.reaction_off and args.initialize_history), 'Disabled chemistry cannot initialize reaction history'
    if args.wall_temperature is not None:
        if not math.isfinite(args.wall_temperature) or args.wall_temperature <= 0:
            parser.error('--wall-temperature must be finite and positive')
        if args.x_boundary not in ('slip_wall', 'no_slip_wall') or not args.transport:
            parser.error('--wall-temperature requires --x-boundary slip_wall/no_slip_wall and --transport')
    root = Path(__file__).resolve().parents[1]
    case = root / 'cases/tgreact_amr_reference'
    contract = json.loads((case / 'acceptance.json').read_text())
    contract['root_shape'] = [args.cells] * 3
    contract['dt_s'] /= args.time_refinement
    contract['steps'] *= args.time_refinement
    limits = contract['acceptance']
    chemistry = json.loads((root / 'mechanisms/lidryer.json').read_text())
    assert digest(root / 'mechanisms/lidryer.yaml') == chemistry['source']['sha256']
    provenance_path = root / 'references/pelec_tgreact_amr.json'
    provenance = json.loads(provenance_path.read_text())
    baseline = json.loads((root / 'references/pelec_baseline.json').read_text())
    assert provenance['pelec_commit'] == baseline['pelec']['commit']
    for name, expected in contract['input_sha256'].items():
        assert digest(case / name) == expected, ('case input drift', name)
    fits = None
    if args.transport:
        assert args.reference_source is not None, '--reference-source required with --transport'
        transport_manifest = json.loads((root / 'references/pelec_tgreact_eb_amr.json').read_text())
        entry = transport_manifest['simple_transport_header']
        assert entry['sha256'] == provenance['source_sha256'][entry['path']]
        fits = import_simple_transport(args.reference_source / entry['path'], chemistry['species'], transport_manifest)
    if args.pelec:
        assert args.reference_source is not None, '--reference-source required with --pelec'
        assert json.loads((args.reference_source / 'pelef-reference-provenance.json').read_text()) == baseline
        for name, expected in provenance['source_sha256'].items():
            assert digest(args.reference_source / name) == expected, ('reference source drift', name)
        args.pelec = args.pelec.resolve(strict=True)
    else:
        args.reference = args.reference.resolve(strict=True)
        captured = json.loads((args.reference / 'reference_manifest.json').read_text())
        assert captured.get('wall_temperature_K') == args.wall_temperature, 'reference wall-temperature mismatch'
        assert captured.get('x_boundary', 'periodic') == args.x_boundary, 'reference boundary mismatch'
        assert captured.get('initial_velocity', 'tg') == args.initial_velocity, 'reference initial velocity mismatch'
        assert captured.get('time_refinement', 1) == args.time_refinement, 'reference temporal refinement mismatch'
        assert captured.get('subcycling', 1) == args.subcycling, 'reference timestep mode mismatch'
        assert captured.get('root_cells', 4) == args.cells, 'reference mesh mismatch'
        assert captured.get('mol_iterations', 2) == args.mol_iterations, 'reference MOL iteration mismatch'
        assert captured.get('reaction_enabled', True) == (not args.reaction_off), 'reference reaction mode mismatch'
        assert captured.get('transport_enabled', False) == args.transport, 'reference transport mode mismatch'
        if fits is not None:
            assert captured['simple_transport_source'] == fits['source'], 'reference transport source mismatch'
        assert captured['source_contract_sha256'] == digest(provenance_path)
        assert captured['case_contract_sha256'] == digest(case / 'acceptance.json')
        for name, expected in captured['files'].items():
            assert digest(args.reference / name) == expected, ('reference artifact drift', name)
    args.probe = args.probe.resolve(strict=True)
    args.exporter = args.exporter.resolve(strict=True)
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='matched_tgreact_', dir=args.work.resolve()))
    print(work, flush=True)
    report = {'scope': __doc__, 'status': 'running', 'contract': contract,
              'source_contract_sha256': digest(provenance_path),
              'case_contract_sha256': digest(case / 'acceptance.json'),
              'bundle_sha256': digest(root / 'mechanisms/lidryer.json'),
              'checker_sha256': digest(Path(__file__)), 'commands': [], 'steps': [], 'violations': []}
    report['initial_reaction_history'] = 'computed' if args.initialize_history else 'imported'
    report['x_boundary'] = args.x_boundary
    report['wall_temperature_K'] = args.wall_temperature
    report['initial_velocity'] = args.initial_velocity
    report['mpi_ranks'] = args.mpi_ranks
    report['probe_sha256'] = digest(args.probe)
    report['time_refinement'] = args.time_refinement
    report['subcycling'] = args.subcycling
    report['root_cells'] = args.cells
    report['reaction_enabled'] = not args.reaction_off
    report['mol_iterations'] = args.mol_iterations
    report['transport_enabled'] = args.transport
    if fits is not None:
        report['simple_transport_source'] = fits['source']
    env = os.environ.copy()
    env.update(HWLOC_COMPONENTS='-gl,-x86', OMP_NUM_THREADS='1')

    def save():
        (work / 'checks.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')

    def run(command, directory, log):
        entry = {'command': list(map(str, command)), 'cwd': str(directory),
                 'binary_sha256': digest(Path(command[0])), 'timeout_s': args.command_timeout}
        report['commands'].append(entry)
        save()
        start = time.monotonic()
        try:
            with log.open('w') as stream:
                result = subprocess.run(entry['command'], cwd=directory, env=env, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=args.command_timeout)
            entry['returncode'] = result.returncode
            assert result.returncode == 0, log
        except subprocess.TimeoutExpired:
            entry['timed_out'] = True
            raise
        finally:
            entry['seconds'] = time.monotonic() - start
            save()

    def check(value, limit, label, step):
        assert math.isfinite(value), (step, label, value)
        if value > limit:
            report['violations'].append({'step': step, 'quantity': label, 'value': value, 'limit': limit})

    names = chemistry['species']
    fields = ['density', 'xmom', 'ymom', 'zmom', 'rho_E'] + ['rho_' + name for name in names] + ['rho_e', 'Temp']
    factors = [1000., 10., 10., 10., .1] + [1000.] * len(names) + [.1, 1.]
    history_fields = ['rho_omega_' + name for name in names] + ['rhoe_dot']
    imported_fields = fields + history_fields
    imported_factors = factors + [1000.] * len(names) + [.1]
    spacing = [length / n for length, n in zip(contract['domain_m'], contract['root_shape'])]
    n = args.cells
    all_keys = {(0, i, j, k) for i in range(1, n+1) for j in range(1, n+1) for k in range(1, n+1)}
    all_keys.update((1, i, j, k) for i in range(1, n+1) for j in range(1, 2*n+1) for k in range(1, 2*n+1))
    keys = sorted(key for key in all_keys if key[0] == 1 or key[1] > n//2)
    volumes = {key: math.prod(spacing) / (8 if key[0] else 1) for key in keys}

    def load(path, reference, step):
        result = {}
        for row in read_rows(path):
            if reference and args.reaction_off:
                # Native omits reaction diagnostics entirely when chemistry
                # is disabled. Only these known inactive fields may be absent.
                for field in history_fields:
                    row.setdefault(field, '0')
                    assert float(row[field]) == 0., ('inactive native reaction history', path, field)
            assert all(math.isfinite(float(value)) for value in row.values()), path
            assert abs(float(row['time' if reference else 'time_s']) - step * contract['dt_s']) <= 1.e-19, path
            level = int(row['level'])
            assert level in (0, 1)
            if reference:
                assert all(abs(float(row['d' + axis]) * .01 / (spacing[n] / 2**level) - 1) <= 1.e-12
                           for n, axis in enumerate('xyz'))
                index = tuple(round(float(row[axis]) / float(row['d' + axis]) - .5) + 1 for axis in 'xyz')
            else:
                assert int(row['node']) == level + 1
                index = tuple(int(row[axis]) for axis in 'ijk')
            key = (level, *index)
            assert key not in result
            result[key] = {field: float(row[field]) * factor for field, factor in
                           zip(imported_fields, imported_factors if reference else [1.] * len(imported_fields))}
        assert set(result) == all_keys, ('matching hierarchy required', path)
        return result

    def inventories(data):
        values = {field: math.fsum(data[key][field] * volumes[key] for key in keys) for field in fields[:-2]}
        for element in ('H', 'O', 'N'):
            values['element_' + element] = math.fsum(
                values['rho_' + sp['name']] * sp['composition'].get(element, 0.) / sp['molecular_weight']
                for sp in chemistry['thermo'])
        return values

    try:
        if args.pelec:
            args.reference = work / 'reference'
            args.reference.mkdir()
            for name in contract['input_sha256']:
                shutil.copyfile(case / name, args.reference / name)
            if args.initial_velocity == 'zero':
                inputs = args.reference / 'inputs'
                text = inputs.read_text()
                assert text.count('prob.u0 = 1000') == 1
                inputs.write_text(text.replace('prob.u0 = 1000', 'prob.u0 = 0'))
            if args.x_boundary != 'periodic':
                native_boundary = {'slip_wall': 'SlipWall', 'no_slip_wall': 'NoSlipWall',
                                   'outflow': 'FOExtrap'}[args.x_boundary]
                inputs = args.reference / 'inputs'
                text = inputs.read_text()
                for old, new in [('geometry.is_periodic = 1 1 1', 'geometry.is_periodic = 0 1 1'),
                                 ('pelec.lo_bc = Interior Interior Interior', f'pelec.lo_bc = {native_boundary} Interior Interior'),
                                 ('pelec.hi_bc = Interior Interior Interior', f'pelec.hi_bc = {native_boundary} Interior Interior')]:
                    assert text.count(old) == 1
                    text = text.replace(old, new)
                inputs.write_text(text)
            if args.wall_temperature is not None:
                with (args.reference / 'inputs').open('a') as stream:
                    stream.write(f'\npelec.do_isothermal_walls = 1\n'
                                 f'pelec.domlo_isothermal_temp = {args.wall_temperature:.17g} -1 -1\n'
                                 f'pelec.domhi_isothermal_temp = {args.wall_temperature:.17g} -1 -1\n')
            if args.time_refinement != 1:
                inputs = args.reference / 'inputs'
                text = inputs.read_text()
                for old, new in [('max_step = 3', f"max_step = {contract['steps']}"),
                                 ('pelec.fixed_dt = 1.e-7', f"pelec.fixed_dt = {contract['dt_s']:.17e}")]:
                    assert text.count(old) == 1
                    text = text.replace(old, new)
                inputs.write_text(text)
            if args.mol_iterations != 2:
                inputs = args.reference / 'inputs'
                text = inputs.read_text()
                assert text.count('pelec.mol_iters = 2') == 1
                inputs.write_text(text.replace('pelec.mol_iters = 2', f'pelec.mol_iters = {args.mol_iterations}'))
            if args.reaction_off:
                inputs = args.reference / 'inputs'
                text = inputs.read_text()
                assert text.count('pelec.do_react = 1') == 1
                inputs.write_text(text.replace('pelec.do_react = 1', 'pelec.do_react = 0'))
            if args.cells != 4:
                inputs = args.reference / 'inputs'
                text = inputs.read_text()
                for old, new in [('amr.n_cell = 4 4 4', f'amr.n_cell = {n} {n} {n}'),
                                 ('amr.max_grid_size = 8', f'amr.max_grid_size = {2*n}')]:
                    assert text.count(old) == 1
                    text = text.replace(old, new)
                inputs.write_text(text)
                (args.reference / 'grids.dat').write_text(f'1\n1\n((0,0,0) ({n//2-1},{n-1},{n-1}) (0,0,0))\n')
            if args.subcycling == 0:
                inputs = args.reference / 'inputs'
                text = inputs.read_text()
                assert text.count('amr.subcycling_mode = Auto') == 1
                inputs.write_text(text.replace('amr.subcycling_mode = Auto', 'amr.subcycling_mode = None'))
            if args.transport:
                inputs = args.reference / 'inputs'
                text = inputs.read_text()
                for name in ('vel', 'temp', 'spec', 'enth'):
                    old = 'pelec.diffuse_' + name + ' = 0'
                    assert text.count(old) == 1
                    text = text.replace(old, old[:-1] + '1')
                inputs.write_text(text)
            run([args.pelec, 'inputs'], args.reference, args.reference / 'run.log')
            manifest = {key: report[key] for key in ('source_contract_sha256', 'case_contract_sha256')}
            manifest['x_boundary'] = args.x_boundary
            manifest['wall_temperature_K'] = args.wall_temperature
            manifest['initial_velocity'] = args.initial_velocity
            manifest['time_refinement'] = args.time_refinement
            manifest['subcycling'] = args.subcycling
            manifest['root_cells'] = args.cells
            manifest['reaction_enabled'] = not args.reaction_off
            manifest['mol_iterations'] = args.mol_iterations
            manifest['transport_enabled'] = args.transport
            if fits is not None:
                manifest['simple_transport_source'] = fits['source']
            manifest['pelec_sha256'] = digest(args.pelec)
            manifest['files'] = {str(p.relative_to(args.reference)): digest(p)
                                 for p in args.reference.rglob('*') if p.is_file()}
            (args.reference / 'reference_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        report['reference'] = str(args.reference)
        reference = []
        for step in range(contract['steps'] + 1):
            exported = work / f'reference_{step:03}.csv'
            run([args.exporter, args.reference / f'plt{step:05}', exported, '--levels'], work,
                work / f'export_{step:03}.log')
            reference.append(load(exported, True, step))
        initial = reference[0]
        temperature_span = max(initial[k]['Temp'] for k in keys) - min(initial[k]['Temp'] for k in keys)
        speed = max(math.sqrt(sum(initial[k][f]**2 for f in ('xmom', 'ymom', 'zmom'))) / initial[k]['density'] for k in keys)
        assert temperature_span >= limits['minimum_initial_temperature_span_K']
        if args.initial_velocity == 'zero':
            assert speed == 0., 'Zero-flow control requires exactly zero initial momentum'
        else:
            assert speed >= limits['minimum_initial_speed_m_s']
        report['initial_temperature_span_K'] = temperature_span
        report['initial_maximum_speed_m_s'] = speed
        lines = ['PELEF_REACTING_AMR_REFERENCE_V2', f"{contract['steps']} {args.mol_iterations} {args.subcycling} {contract['dt_s']:.17e}",
                 ' '.join(f'{v:.17e}' for v in spacing), str(len(all_keys))]
        if fits is not None:
            lines[0] = 'PELEF_REACTING_AMR_REFERENCE_V3'
            rows = fits['viscosity'] + fits['conductivity'] + fits['binary_diffusion']
            lines[3:3] = [' '.join(f'{v:.17e}' for v in row) for row in rows]
        if args.cells != 4:
            lines[0] = 'PELEF_REACTING_AMR_REFERENCE_V4'
            lines.insert(1, f'{args.cells} {int(args.transport)}')
        for (level, i, j, k), values in sorted(initial.items()):
            lines.append(f'{level+1} {i} {j} {k} ' + ' '.join(f'{values[f]:.17e}' for f in imported_fields))
        imported = work / 'initial.dat'
        imported.write_text('\n'.join(lines) + '\n')
        report['import_sha256'] = digest(imported)
        launcher = [args.mpiexec.resolve(strict=True), '-n', str(args.mpi_ranks)] if args.mpi_ranks else []
        run(launcher + [args.probe, imported, work / 'fortran'] + (['--initialize-history'] if args.initialize_history else []) +
            (['--reaction-off'] if args.reaction_off else []) +
            (['--slip-x'] if args.x_boundary == 'slip_wall' else []) +
            (['--no-slip-x'] if args.x_boundary == 'no_slip_wall' else []) +
            (['--outflow-x'] if args.x_boundary == 'outflow' else []) +
            ([f'--wall-temperature={args.wall_temperature:.17g}'] if args.wall_temperature is not None else []),
            work, work / 'probe.log')
        initial_inventory = inventories(initial)
        final_inventory = {}
        for step in range(contract['steps'] + 1):
            actual = load(work / f'fortran_{step:03}.csv', False, step)
            native = reference[step]
            if step == 0:
                identity_fields = fields if args.initialize_history else imported_fields
                assert all(actual[k][f] == native[k][f] for k in all_keys for f in identity_fields), \
                    'identical initial gas, UEINT, temperature (and imported history) required'
                if args.initialize_history:
                    report['initial_history_errors'] = {
                        f: {'absolute': max(abs(actual[k][f] - native[k][f]) for k in all_keys),
                            'reference_maximum': max(abs(native[k][f]) for k in all_keys)} for f in history_fields}
            errors = {}
            record = {'step': step, 'time_s': step * contract['dt_s'], 'fields': errors}
            report['steps'].append(record)
            for field in fields:
                absolute = max(abs(actual[k][field] - native[k][field]) for k in keys)
                scale = max(1.e-30, max(abs(native[k][field]) for k in keys))
                errors[field] = {'absolute': absolute, 'scaled': absolute / scale,
                                 'worst_cell': max(keys, key=lambda k: abs(actual[k][field] - native[k][field]))}
                if step == 0:
                    check(absolute / max(1., scale), limits['initial_scaled'], 'initial:' + field, step)
                if field.startswith('rho_') and field not in ('rho_E', 'rho_e'):
                    y_error = max(abs(actual[k][field] / actual[k]['density'] - native[k][field] / native[k]['density']) for k in keys)
                    errors[field]['mass_fraction_absolute'] = y_error
                    check(y_error, limits['species_mass_fraction_absolute'], field, step)
                elif field == 'Temp':
                    check(absolute, limits['temperature_absolute_K'], field, step)
                elif field in ('xmom', 'ymom', 'zmom'):
                    check(absolute, limits['momentum_absolute_kg_m2_s'], field, step)
                else:
                    check(absolute / scale, limits['density_scaled' if field == 'density' else 'energy_scaled'], field, step)
            for label, data in (('reference', native), ('fortran', actual)):
                assert all(data[k]['rho_' + n] >= 0 for k in keys for n in names), (label, step, 'negative species')
                closure = max(abs(sum(data[k]['rho_' + n] for n in names) / data[k]['density'] - 1) for k in keys)
                check(closure, limits['inventory_relative'], label + ':species_sum', step)
                inv = inventories(data)
                record[label + '_inventory'] = inv
                final_inventory[label] = inv
                for field in (() if args.x_boundary == 'outflow' else
                              ('density', 'element_H', 'element_O', 'element_N')):
                    drift = abs(inv[field] / initial_inventory[field] - 1)
                    check(drift, limits['inventory_relative'], label + ':conservation:' + field, step)
            budget = read_rows(work / f'fortran_{step:03}.budget.csv')
            assert len(budget) == 1 and int(budget[0]['step']) == step
            assert abs(float(budget[0]['time_s']) - step * contract['dt_s']) <= 1.e-19
            correction = float(budget[0]['numerical_energy_J'])
            ledger = {f: {'boundary': float(budget[0]['boundary_' + f]),
                          'source': float(budget[0]['source_' + f])} for f in fields[:5]}
            record['boundary_source_ledger'] = ledger
            if args.x_boundary == 'outflow':
                # Open-domain inventories need boundary supply. Native face
                # ledgers are not exported: native checks here are inventory
                # comparisons, not independent native conservation claims.
                species_ledger = {f: {kind: float(budget[0][kind + '_' + f])
                                     for kind in ('boundary', 'source')}
                                  for f in fields[5:-2]}
                record['species_boundary_source_ledger'] = species_ledger
                record['element_budget_relative_residual'] = {}
                for element in ('H', 'O', 'N'):
                    field = 'element_' + element
                    supply = math.fsum(
                        sum(species_ledger['rho_' + sp['name']].values()) *
                        sp['composition'].get(element, 0.) / sp['molecular_weight']
                        for sp in chemistry['thermo'])
                    residual = abs(record['fortran_inventory'][field] - initial_inventory[field] - supply)
                    residual /= max(abs(initial_inventory[field]), abs(supply), 1.e-30)
                    record['element_budget_relative_residual'][field] = residual
                    check(residual, limits['inventory_relative'], 'fortran:open_budget:' + field, step)
                for field in ('density', 'element_H', 'element_O', 'element_N'):
                    difference = abs(record['fortran_inventory'][field] - record['reference_inventory'][field])
                    difference /= max(abs(initial_inventory[field]), 1.e-30)
                    check(difference, limits['inventory_relative'], 'open_inventory_reference:' + field, step)
            record['gas_budget_relative_residual'] = {}
            for field in fields[:5]:
                initial_value = initial_inventory[field]
                final_value = record['fortran_inventory'][field]
                delta = ledger[field]['boundary'] + ledger[field]['source']
                if field == 'rho_E':
                    delta += correction
                absolute_inventory = math.fsum(abs(actual[k][field]) * volumes[k] for k in keys)
                scale = max(abs(initial_value), abs(final_value), abs(delta), absolute_inventory, 1.e-30)
                gas_residual = abs(final_value - initial_value - delta) / scale
                record['gas_budget_relative_residual'][field] = gas_residual
                check(gas_residual, limits['inventory_relative'], 'fortran:gas_budget:' + field, step)
            residual = (record['fortran_inventory']['rho_E'] - initial_inventory['rho_E'] - correction -
                        ledger['rho_E']['boundary'] - ledger['rho_E']['source']) / abs(initial_inventory['rho_E'])
            record['numerical_energy_J'] = correction
            record['energy_budget_relative_residual'] = residual
            check(abs(residual), limits['energy_budget_relative'], 'fortran:energy_budget', step)
            save()
        consumed = {}
        for label, inv in final_inventory.items():
            supply = (species_ledger['rho_H2']['boundary']
                      if args.x_boundary == 'outflow' else 0.)
            consumed[label] = (initial_inventory['rho_H2'] + supply - inv['rho_H2']) / initial_inventory['rho_H2']
            if args.reaction_off:
                check(abs(consumed[label]), limits['inventory_relative'], label + ':inactive_H2_conservation', contract['steps'])
            else:
                check(limits['minimum_H2_consumption_fraction'] - consumed[label], 0., label + ':minimum_H2_consumption', contract['steps'])
        report['H2_consumed_fraction'] = consumed
        if args.x_boundary == 'outflow':
            report['H2_consumption_basis'] = 'Inventory change adjusted by candidate boundary H2 supply for both solvers; native boundary flux is not independently exported'
        if not args.reaction_off and consumed['reference'] > 0:
            relative_error = abs(consumed['reference'] - consumed['fortran']) / consumed['reference']
            report['H2_consumption_relative_error'] = relative_error
            check(relative_error, limits['H2_consumption_relative'], 'H2_consumption_relative_error', contract['steps'])
        assert not report['violations'], report['violations']
        report['status'] = 'passed'
        print('PASS: reacting two-level spatial-gradient MOL trajectory', flush=True)
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = repr(error)
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
