#!/usr/bin/env python3
"""Compare nonuniform reacting EB-AMR MOL with unchanged pinned TGReact.

Frozen two-level planar EB, Simple molecular transport and independent initial
reaction history. Not application initialization, convergence or ignition validation.
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


def reaction_controls(study):
    """Scalar CVODE tolerances with explicit density units.

Native evolves CGS densities plus temperature; PeleF evolves reduced SI
densities and recovers temperature through its EOS. Density tolerance
matching does not make those reactor formulations identical.
"""
    controls = {'relative': 1.e-10, 'fortran_absolute_mass_fraction': 1.e-14,
                'fortran_atol_basis': 'mass_fraction_times_local_density',
                'native_absolute_density_g_cm3': 1.e-18,
                'native_absolute_temperature_K': 1.e-18}
    if study is not None and 'reaction_tolerances' in study:
        settings = study['reaction_tolerances']
        assert set(settings) == {'relative', 'absolute_density_kg_m3'}
        rtol, atol = settings['relative'], settings['absolute_density_kg_m3']
        assert math.isfinite(rtol) and 0 < rtol <= 1.e-10
        assert math.isfinite(atol) and 0 < atol <= 1.e-15, 'Must not loosen native density tolerance'
        native_atol = atol * .001  # kg/m3 -> g/cm3
        assert native_atol > 0, 'Native absolute tolerance underflow'
        del controls['fortran_absolute_mass_fraction']
        controls.update(relative=rtol, fortran_absolute_density_kg_m3=atol,
                        fortran_atol_basis='density_kg_m3',
                        native_absolute_density_g_cm3=native_atol,
                        native_absolute_temperature_K=native_atol)
    return controls


def eos_controls(study):
    """Record the explicit PeleF intervention separately from the native policy."""
    convergence = False
    if study is not None and 'fortran_eos_convergence' in study:
        convergence = study['fortran_eos_convergence']
        assert type(convergence) is bool, 'fortran_eos_convergence must be a boolean'
        assert 'reaction_tolerances' in study, 'V4 requires explicit density reaction tolerances'
    return {'fortran_converged_inverse': convergence, 'native_converged_inverse': False}


def velocity_controls(study):
    """Bounded initial-flow controls, separate from the frozen case."""
    amplitude = 1000.
    if study is not None and 'initial_velocity_cm_s' in study:
        amplitude = study['initial_velocity_cm_s']
        assert type(amplitude) in (int, float) and amplitude in (0., 1., 1000.), 'Only zero/weak/default TG velocity supported'
    return {'tg_amplitude_cm_s': float(amplitude)}


def energy_consistency(data, keys, volumes):
    """Report stored versus total-derived internal energy in SI units.

    This is diagnostic only: independent UEINT can intentionally differ from
    total-minus-kinetic energy, especially after reflux or under cancellation.
    Callers supply only composite fluid cells, excluding covered coarse cells.
    """
    gaps = {}
    for key in keys:
        row = data[key]
        density = row['density']
        assert math.isfinite(density) and density > 0, ('energy diagnostic density', key)
        kinetic = .5 * math.fsum(row[f]**2 for f in ('xmom', 'ymom', 'zmom')) / density
        gap = math.fsum((row['rho_E'], -row['rho_e'], -kinetic))
        assert math.isfinite(gap), ('nonfinite energy diagnostic', key)
        gaps[key] = gap
    result = {}
    for level in sorted({key[0] for key in keys}):
        cells = [key for key in keys if key[0] == level]
        worst = max(cells, key=lambda key: abs(gaps[key]))
        result[str(level)] = {
            'cell_count': len(cells),
            'maximum_absolute_J_m3': abs(gaps[worst]),
            'signed_at_worst_J_m3': gaps[worst],
            'worst_cell': worst,
            'volume_integral_J': math.fsum(gaps[key] * volumes[key] for key in cells),
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('probe', 'exporter', 'work'):
        parser.add_argument('--' + name, required=True, type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--pelec', type=Path)
    mode.add_argument('--reference', type=Path, help='Captured native reference with its SHA manifest')
    parser.add_argument('--reference-source', type=Path, required=True)
    parser.add_argument('--transport-off', action='store_true', help='Matched diagnostic with all molecular transport disabled')
    studies = parser.add_mutually_exclusive_group()
    studies.add_argument('--baseline-mesh-study', type=Path,
                         help='Mesh/time study retaining binaries and mechanism from a passed original comparison')
    studies.add_argument('--mesh-study', type=Path,
                        help='Explicit held-EOS/native-activation diagnostic configuration; never pinned qualification')
    policies = parser.add_mutually_exclusive_group()
    policies.add_argument('--momentum-only', action='store_true', help='Separate-momentum reacting intervention')
    policies.add_argument('--velocity-only', action='store_true', help='Consistent-velocity reacting intervention')
    policies.add_argument('--paired-pressure', action='store_true',
                        help='Compare the opt-in paired pressure policy; never pinned algorithm qualification')
    args = parser.parse_args()
    policy_option = ('--paired-reacting' if args.paired_pressure else '--velocity-reacting'
                     if args.velocity_only else '--momentum-reacting' if args.momentum_only else None)
    root = Path(__file__).resolve().parents[1]
    case = root / 'cases/tgreact_eb_amr_reference'
    contract = json.loads((case / 'acceptance.json').read_text())
    limits = contract['acceptance']
    chemistry = json.loads((root / 'mechanisms/lidryer.json').read_text())
    assert digest(root / 'mechanisms/lidryer.yaml') == chemistry['source']['sha256']
    study = None
    study_path = args.baseline_mesh_study or args.mesh_study
    if study_path:
        assert args.pelec and not args.transport_off, 'Mesh study requires a live, transport-enabled diagnostic'
        study = json.loads(study_path.read_text())
        assert study['schema_version'] == 1
        assert study['cells'] in (8, 16, 32)
        if 'profile_points' in study:
            assert type(study['profile_points']) is int and study['profile_points'] in (129, 257, 513)
            assert (study['profile_points'] - 1) % (4 * study['cells']) == 0, 'Profile must sample every fine cell center'
        assert isinstance(study['steps'], int) and 1 <= study['steps'] <= 100
        assert math.isfinite(study['dt_s']) and 0 < study['dt_s'] <= 1.e-7
        assert 60 <= study['timeout_s'] <= 3600
        contract['root_shape'] = [study['cells']] * 3
        contract['fine_coarse_box'] = [[1] * 3, [study['cells'] // 2] * 3]
        contract['steps'], contract['dt_s'] = study['steps'], study['dt_s']
        if args.baseline_mesh_study:
            assert policy_option is None, 'Baseline study cannot select an experimental halo policy'
            assert set(study) == {'schema_version', 'cells', 'steps', 'dt_s', 'timeout_s', 'baseline_checks'}, \
                'Baseline study only changes mesh/time controls'
            baseline_report_path = Path(study['baseline_checks']).resolve(strict=True)
            baseline_report = json.loads(baseline_report_path.read_text())
            assert baseline_report['status'] == 'passed' and not baseline_report['violations']
            assert 'mesh_study' not in baseline_report and baseline_report['halo_policy_option'] is None
            assert baseline_report['contract'] == json.loads((case / 'acceptance.json').read_text())
            assert baseline_report['bundle_sha256'] == digest(root / 'mechanisms/lidryer.json')
            assert baseline_report['source_contract_sha256'] == digest(root / 'references/pelec_tgreact_eb_amr.json')
            assert baseline_report['case_contract_sha256'] == digest(case / 'acceptance.json')
            for executable in (args.probe, args.exporter):
                matches = [entry for entry in baseline_report['commands']
                           if Path(entry['command'][0]).resolve() == executable.resolve()]
                assert matches and all(entry.get('returncode') == 0 and
                                       entry['binary_sha256'] == digest(executable) for entry in matches), \
                    ('baseline binary drift', executable)
            captured_path = Path(baseline_report['reference']) / 'reference_manifest.json'
            captured_baseline = json.loads(captured_path.read_text())
            assert 'mesh_study_sha256' not in captured_baseline
            assert captured_baseline['pelec_sha256'] == digest(args.pelec), 'Pinned native binary drift'
            for key in ('source_contract_sha256', 'case_contract_sha256'):
                assert captured_baseline[key] == baseline_report[key]
            for name, expected in captured_baseline['files'].items():
                assert digest(captured_path.parent / name) == expected, ('baseline artifact drift', name)
            policy_identity = {str(path.resolve()): digest(path) for path in (args.pelec, args.exporter)}
        else:
            mechanism_record = json.loads(Path(study['mechanism_record']).read_text())
            assert mechanism_record['status'] == 'generated' and mechanism_record['qualification'] is False
            for path, expected in mechanism_record['outputs'].items():
                assert digest(Path(path)) == expected, ('study mechanism drift', path)
            bundle = Path(study['mechanism_bundle']).resolve()
            assert str(bundle) in mechanism_record['outputs']
            variant = json.loads(bundle.read_text())
            assert digest(bundle.parent / variant['source']['file']) == variant['source']['sha256']
            for field in ('species', 'thermo', 'transport'):
                assert variant[field] == chemistry[field], ('study mechanism changed', field)
            from prepare_lidryer_native_activation import verify_activation_only
            verify_activation_only(chemistry, variant, mechanism_record['activation_energy_ratio'])
            chemistry = variant
            policy_identity = json.loads(Path(study['native_policy_identity']).read_text())
            for path in (args.pelec, args.exporter):
                assert digest(path) == policy_identity[str(path.resolve())], ('diagnostic binary drift', path)
            probe_record = json.loads(Path(study['probe_build_record']).read_text())
            assert digest(args.probe) == probe_record['native_activation_binary_sha256']
            assert digest(bundle) == probe_record['mechanism_bundle_sha256']
    reaction = reaction_controls(study)
    eos = eos_controls(study)
    velocity = velocity_controls(study)
    provenance_path = root / 'references/pelec_tgreact_eb_amr.json'
    provenance = json.loads(provenance_path.read_text())
    baseline = json.loads((root / 'references/pelec_baseline.json').read_text())
    assert provenance['pelec_commit'] == baseline['pelec']['commit']
    for name, expected in contract['input_sha256'].items():
        assert digest(case / name) == expected, ('case input drift', name)
    assert json.loads((args.reference_source / 'pelef-reference-provenance.json').read_text()) == baseline
    for name, expected in provenance['source_sha256'].items():
        assert digest(args.reference_source / name) == expected, ('reference source drift', name)
    fits = import_simple_transport(args.reference_source / provenance['simple_transport_header']['path'],
                                   chemistry['species'], provenance)
    fit_rows = fits['viscosity'] + fits['conductivity'] + fits['binary_diffusion']
    if args.pelec:
        args.pelec = args.pelec.resolve(strict=True)
    else:
        args.reference = args.reference.resolve(strict=True)
        captured = json.loads((args.reference / 'reference_manifest.json').read_text())
        assert 'mesh_study_sha256' not in captured, 'Diagnostic mesh captures cannot qualify the pinned reference'
        assert captured['transport_enabled'] == (not args.transport_off)
        assert captured['source_contract_sha256'] == digest(provenance_path)
        assert captured['case_contract_sha256'] == digest(case / 'acceptance.json')
        for name, expected in captured['files'].items():
            assert digest(args.reference / name) == expected, ('reference artifact drift', name)
    args.probe = args.probe.resolve(strict=True)
    args.exporter = args.exporter.resolve(strict=True)
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='matched_eb_amr_', dir=args.work.resolve()))
    print(work, flush=True)
    report = {'scope': __doc__, 'status': 'running', 'contract': contract,
              'source_contract_sha256': digest(provenance_path),
              'case_contract_sha256': digest(case / 'acceptance.json'),
              'bundle_sha256': digest(root / 'mechanisms/lidryer.json'),
              'checker_sha256': digest(Path(__file__)), 'commands': [], 'steps': [], 'violations': []}
    report['initial_reaction_history'] = 'computed'
    report['reaction_controls'] = reaction
    report['eos_controls'] = eos
    report['initial_profile_controls'] = {'kind': 'frozen_table', 'points': 33}
    if args.baseline_mesh_study:
        report.update(scope='Mesh/time study of unchanged pinned native binary and original mechanism',
                      qualification=False, mesh_study=study, mesh_study_sha256=digest(study_path),
                      baseline_checks_sha256=digest(baseline_report_path),
                      baseline_reference_manifest_sha256=digest(captured_path),
                      native_policy_identity=policy_identity, mechanism_basis='original_pinned')
    elif study is not None:
        report.update(scope='Mesh diagnostic: native held spatial EOS and explicit native activation basis',
                      qualification=False, mesh_study=study, mesh_study_sha256=digest(study_path),
                      bundle_sha256=digest(Path(study['mechanism_bundle'])),
                      probe_build_record=probe_record, mechanism_record=mechanism_record,
                      native_policy_identity=policy_identity)
    report['paired_pressure'] = args.paired_pressure
    report['halo_policy_option'] = policy_option
    if policy_option:
        report.update(qualification=False,
                      scope='Halo-policy intervention compared against unchanged pinned reference')
    report['transport_enabled'] = not args.transport_off
    report['simple_transport_source'] = fits['source']
    report['energy_consistency_definition'] = {
        'quantity': 'rho_E - rho_e - (xmom^2 + ymom^2 + zmom^2)/(2*density)',
        'units': 'J/m3', 'cell_selection': 'composite fluid cells by level',
        'acceptance_gate': False,
    }
    env = os.environ.copy()
    report['initial_velocity_controls'] = velocity
    env.update(HWLOC_COMPONENTS='-gl,-x86', OMP_NUM_THREADS='1')

    def save():
        (work / 'checks.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')

    def run(command, directory, log):
        entry = {'command': list(map(str, command)), 'cwd': str(directory),
                 'binary_sha256': digest(Path(command[0]))}
        report['commands'].append(entry)
        save()
        start = time.monotonic()
        try:
            with log.open('w') as stream:
                result = subprocess.run(entry['command'], cwd=directory, env=env, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=study['timeout_s'] if study else 240)
            entry['returncode'] = result.returncode
            assert result.returncode == 0, log
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
    root_size = contract['root_shape'][0]
    all_keys = {(level, i, j, k) for level in (0, 1) for i in range(1, root_size + 1)
                for j in range(1, root_size + 1) for k in range(1, root_size + 1)}
    def fluid_fraction(key):
        dz = spacing[2] / 2**key[0]
        return min(1., max(0., (key[3] * dz - contract['plane_z_m']) / dz))

    keys = sorted(key for key in all_keys if (key[0] == 1 or any(v > root_size // 2 for v in key[1:])) and fluid_fraction(key) > 0)
    volumes = {key: math.prod(spacing) / (8 if key[0] else 1) * fluid_fraction(key) for key in keys}

    def load(path, reference, step):
        result = {}
        for row in read_rows(path):
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
            assert abs(float(row['vfrac']) - fluid_fraction(key)) <= 1.e-13, ('EB geometry mismatch', path, key)
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
            if study is not None and 'profile_points' in study:
                from prepare_tgreact_profile import prepare
                profile_manifest = prepare(case / 'profile.csv', args.reference / 'profile.csv', study['profile_points'])
                report['initial_profile_manifest'] = profile_manifest
                report['initial_profile_generator_sha256'] = digest(root / 'tools/prepare_tgreact_profile.py')
                report['initial_profile_controls'] = {'kind': 'cosine_sampling_diagnostic', 'points': study['profile_points']}
            report['initial_profile_sha256'] = digest(args.reference / 'profile.csv')
            if args.transport_off:
                input_file = args.reference / 'inputs'
                text = input_file.read_text()
                for field in ('vel', 'temp', 'spec', 'enth'):
                    text = text.replace(f'pelec.diffuse_{field} = 1', f'pelec.diffuse_{field} = 0')
                input_file.write_text(text)
            native_command = [args.pelec, 'inputs']
            if study is not None and 'profile_points' in study:
                native_command.append(f'prob.nx={study["profile_points"]}')
            if study is not None:
                half = root_size // 2 - 1
                (args.reference / 'grids.dat').write_text(f'1\n1\n((0,0,0) ({half},{half},{half}) (0,0,0))\n')
                native_command += [f'amr.n_cell={root_size} {root_size} {root_size}',
                                   f'amr.max_grid_size={root_size}', f'max_step={contract["steps"]}',
                                   f'stop_time={contract["steps"] * contract["dt_s"]:.17e}',
                                   f'pelec.fixed_dt={contract["dt_s"]:.17e}']
                if 'initial_velocity_cm_s' in study:
                    native_command.append(f'prob.u0={velocity["tg_amplitude_cm_s"]:.17e}')
                if 'reaction_tolerances' in study:
                    native_command += [f'ode.rtol={reaction["relative"]:.17e}',
                                       f'ode.atol={reaction["native_absolute_density_g_cm3"]:.17e}',
                                       'pelec.use_typ_vals_chem=0', 'pelec.use_typ_vals_chem_usr=0']
            run(native_command, args.reference, args.reference / 'run.log')
            manifest = {key: report[key] for key in ('source_contract_sha256', 'case_contract_sha256')}
            manifest['transport_enabled'] = not args.transport_off
            manifest['pelec_sha256'] = digest(args.pelec)
            if study is not None:
                manifest.update(mesh_study_sha256=digest(study_path), qualification=False)
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
        if velocity['tg_amplitude_cm_s'] == 0:
            assert study is not None and speed == 0., 'Zero-velocity intervention was not applied'
        else:
            assert speed >= limits['minimum_initial_speed_m_s'] * velocity['tg_amplitude_cm_s'] / 1000.
        report['initial_temperature_span_K'] = temperature_span
        report['initial_maximum_speed_m_s'] = speed
        lines = ['PELEF_REACTING_EB_REFERENCE_V1', f"{contract['steps']} 2 1 {contract['dt_s']:.17e}",
                 ' '.join(f'{v:.17e}' for v in spacing),
                 f"{contract['plane_z_m']:.17e} {int(not args.transport_off)}"]
        if study is not None:
            lines[0] = 'PELEF_REACTING_EB_REFERENCE_V2'
            lines.insert(1, str(root_size))
            if 'reaction_tolerances' in study:
                lines[0] = 'PELEF_REACTING_EB_REFERENCE_V3'
                lines.insert(2, f'{reaction["relative"]:.17e} {reaction["fortran_absolute_density_kg_m3"]:.17e}')
            if 'fortran_eos_convergence' in study:
                lines[0] = 'PELEF_REACTING_EB_REFERENCE_V4'
                lines.insert(3, str(int(eos['fortran_converged_inverse'])))
        lines += [' '.join(f'{v:.17e}' for v in row) for row in fit_rows] + [str(len(all_keys))]
        for (level, i, j, k), values in sorted(initial.items()):
            lines.append(f'{level+1} {i} {j} {k} ' + ' '.join(f'{values[f]:.17e}' for f in imported_fields))
        imported = work / 'initial.dat'
        imported.write_text('\n'.join(lines) + '\n')
        report['import_sha256'] = digest(imported)
        run([args.probe, imported, work / 'fortran'] +
            ([policy_option] if policy_option else []), work, work / 'probe.log')
        initial_inventory = inventories(initial)
        final_inventory = {}
        for step in range(contract['steps'] + 1):
            actual = load(work / f'fortran_{step:03}.csv', False, step)
            native = reference[step]
            if step == 0:
                identity_fields = fields
                assert all(actual[k][f] == native[k][f] for k in all_keys for f in identity_fields), \
                    'identical initial gas, UEINT and temperature required'
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
                record.setdefault('energy_consistency', {})[label] = energy_consistency(data, keys, volumes)
                final_inventory[label] = inv
                for field in ('density', 'element_H', 'element_O', 'element_N'):
                    drift = abs(inv[field] / initial_inventory[field] - 1)
                    check(drift, limits['inventory_relative'], label + ':conservation:' + field, step)
            budget = read_rows(work / f'fortran_{step:03}.budget.csv')
            assert len(budget) == 1 and int(budget[0]['step']) == step
            assert abs(float(budget[0]['time_s']) - step * contract['dt_s']) <= 1.e-19
            correction = float(budget[0]['numerical_energy_J'])
            residual = (record['fortran_inventory']['rho_E'] - initial_inventory['rho_E'] - correction) / abs(initial_inventory['rho_E'])
            record['numerical_energy_J'] = correction
            record['energy_budget_relative_residual'] = residual
            check(abs(residual), limits['energy_budget_relative'], 'fortran:energy_budget', step)
            save()
        consumed = {}
        for label, inv in final_inventory.items():
            consumed[label] = (initial_inventory['rho_H2'] - inv['rho_H2']) / initial_inventory['rho_H2']
            check(limits['minimum_H2_consumption_fraction'] - consumed[label], 0., label + ':minimum_H2_consumption', contract['steps'])
        report['H2_consumed_fraction'] = consumed
        if consumed['reference'] > 0:
            relative_error = abs(consumed['reference'] - consumed['fortran']) / consumed['reference']
            report['H2_consumption_relative_error'] = relative_error
            check(relative_error, limits['H2_consumption_relative'], 'H2_consumption_relative_error', contract['steps'])
        assert not report['violations'], report['violations']
        report['status'] = 'diagnostic_passed' if study is not None or policy_option is not None else 'passed'
        print('PASS: diagnostic comparison (not pinned qualification)' if study is not None or policy_option is not None else
              'PASS: reacting spatial-gradient EB-AMR MOL trajectory', flush=True)
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = repr(error)
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
