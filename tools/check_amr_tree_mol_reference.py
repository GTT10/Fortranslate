"""Qualify imported two-level periodic HIT states against pinned MOL runs.

This is a frozen gamma-law operator/AMR control, not combustion validation.
The reference folders must contain the complete plotfile trajectory and inputs.
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('exporter', 'probe', 'work'):
        parser.add_argument('--' + key, type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--reference', type=Path, help='Reuse a captured reference with its SHA manifest')
    mode.add_argument('--pelec', type=Path, help='Run the unchanged pinned PeleC-HIT application')
    parser.add_argument('--reference-source', type=Path)
    parser.add_argument('--case', choices=('sub0_it1', 'sub0_it2', 'sub1_it1', 'sub1_it2'),
                        help='Limit execution to one clock/iteration combination')
    parser.add_argument('--keep-going', action='store_true',
                        help='Continue after field/inventory gate failures; still exit nonzero. '
                             'Execution, provenance and malformed-data failures remain fatal.')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    contract_path = root / 'references/pelec_amr_mol.json'
    contract = json.loads(contract_path.read_text())
    limits = contract['acceptance']
    baseline = json.loads((root / 'references/pelec_baseline.json').read_text())
    assert contract['pelec_commit'] == baseline['pelec']['commit']
    if args.pelec:
        assert args.reference_source is not None, '--reference-source is required with --pelec'
        args.pelec = args.pelec.resolve(strict=True)
        assert json.loads((args.reference_source / 'pelef-reference-provenance.json').read_text()) == baseline
        for name, expected in contract['source_sha256'].items():
            assert digest(args.reference_source / name) == expected, ('reference source drift', name)
    else:
        args.reference = args.reference.resolve(strict=True)
        captured = json.loads((args.reference / 'reference_manifest.json').read_text())
        assert captured['source_contract_sha256'] == digest(contract_path)
        for name, expected in captured['files'].items():
            assert digest(args.reference / name) == expected, ('reference artifact drift', name)
    args.exporter = args.exporter.resolve(strict=True)
    args.probe = args.probe.resolve(strict=True)
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='matched_amr_mol_', dir=args.work.resolve()))
    if args.pelec:
        args.reference = work / 'reference'
        args.reference.mkdir()
    print(work, flush=True)
    fields = ('density', 'xmom', 'ymom', 'zmom', 'rho_E', 'rho_AIR', 'Temp')
    factors = (1000, 10, 10, 10, .1, 1000, 1)
    checks = {'scope': contract['scope'], 'status': 'running', 'acceptance': limits,
              'keep_going': args.keep_going, 'matrix_complete': False,
              'requested_cases': [args.case] if args.case else
                                 ['sub0_it1', 'sub0_it2', 'sub1_it1', 'sub1_it2'],
              'reference': str(args.reference), 'probe_sha256': digest(args.probe),
              'exporter_sha256': digest(args.exporter), 'checker_sha256': digest(Path(__file__)),
              'source_contract_sha256': digest(contract_path), 'commands': [], 'cases': []}
    env = os.environ.copy()
    env.update(HWLOC_COMPONENTS='-gl,-x86', OMP_NUM_THREADS='1')

    def save():
        (work/'checks.json').write_text(json.dumps(checks, indent=2, allow_nan=False)+'\n')

    def run(command, directory, log):
        entry = {'command': list(map(str, command)), 'cwd': str(directory), 'binary_sha256': digest(Path(command[0]))}
        checks['commands'].append(entry)
        save()
        start = time.monotonic()
        try:
            with log.open('w') as stream:
                result = subprocess.run(entry['command'], cwd=directory, env=env,
                                        stdout=stream, stderr=subprocess.STDOUT, timeout=120)
            entry['returncode'] = result.returncode
            assert result.returncode == 0, log
        finally:
            entry['seconds'] = time.monotonic()-start
            save()

    def load(path, reference, expected_time):
        rows = {}
        with path.open() as stream:
            for row in csv.DictReader(stream):
                assert all(math.isfinite(float(value)) for value in row.values()), path
                observed_time = float(row['time' if reference else 'time_s'])
                assert abs(observed_time-expected_time) <= 1.e-19, ('physical time', path, observed_time, expected_time)
                level = int(row['level'])
                assert level in (0, 1)
                if reference:
                    dx = math.pi*.005 / 2**level
                    assert all(abs(float(row['d'+axis])*.01/dx-1) <= 1.e-13 for axis in 'xyz')
                index = tuple(int(round(float(row[axis])/float(row['d'+axis])-.5))+1
                              for axis in 'xyz') if reference else tuple(int(row[axis]) for axis in 'ijk')
                key = (level, *index)
                assert key not in rows
                rows[key] = tuple(float(row[field])*factor for field, factor in
                                  zip(fields, factors if reference else (1,)*7))
        assert len(rows) == 320
        expected = {(0, i, j, k) for i in range(1, 5) for j in range(1, 5) for k in range(1, 5)}
        expected.update((1, i, j, k) for i in range(1, 5) for j in range(1, 9) for k in range(1, 9))
        assert set(rows) == expected, 'matching hierarchy required'
        return rows

    def composite(key):
        level, i, _, _ = key
        return level == 1 or i > 2

    try:
        for recursive in (False, True):
            for iterations in (1, 2):
                name = f'sub{int(recursive)}_it{iterations}'
                if args.case and name != args.case:
                    continue
                folder = work / name; folder.mkdir()
                reference = args.reference / name
                if args.pelec:
                    reference.mkdir()
                    text = (root/'cases/amr_mol_reference/inputs').read_text()
                    text = text.replace('pelec.mol_iters = 2', f'pelec.mol_iters = {iterations}')
                    text = text.replace('amr.subcycling_mode = Auto', 'amr.subcycling_mode = ' + ('Auto' if recursive else 'None'))
                    (reference/'inputs').write_text(text)
                    shutil.copyfile(root/'cases/amr_mol_reference/grids.dat', reference/'grids.dat')
                    shutil.copyfile(args.reference_source/'Exec/RegTests/HIT/hit_ic_4_32.dat', reference/'hit_ic_4_32.dat')
                    run([args.pelec, 'inputs'], reference, reference/'run.log')
                    manifest = {'source_contract_sha256': digest(contract_path), 'pelec_sha256': digest(args.pelec),
                                'files': {str(p.relative_to(args.reference)): digest(p)
                                          for p in args.reference.rglob('*') if p.is_file() and p.name != 'reference_manifest.json'}}
                    (args.reference/'reference_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
                native = []
                for step in range(4):
                    output = folder / f'reference_{step:03}.csv'
                    run([args.exporter, reference/f'plt{step:05}', output, '--levels'], folder,
                        folder/f'export_{step:03}.log')
                    native.append(load(output, True, step*1.e-7))
                spacing = [math.pi*.005]*3  # CGS domain 2*pi cm, four coarse cells.
                lines = ['PELEF_AMR_MOL_REFERENCE_V1', f'3 {iterations} {int(recursive)} 1e-7',
                         ' '.join(f'{v:.17e}' for v in spacing), '320']
                for (level, i, j, k), values in sorted(native[0].items()):
                    lines.append(f'{level+1} {i} {j} {k} ' + ' '.join(f'{v:.17e}' for v in values))
                input_file = folder / 'initial.dat'; input_file.write_text('\n'.join(lines)+'\n')
                run([args.probe, input_file, folder/'fortran'], folder, folder/'probe.log')
                result = {'case': name, 'initial_sha256': hashlib.sha256(input_file.read_bytes()).hexdigest(),
                          'steps': [], 'violations': []}
                checks['cases'].append(result)
                for step in range(4):
                    actual = load(folder/f'fortran_{step:03}.csv', False, step*1.e-7)
                    assert set(actual) == set(native[step])
                    keys = [key for key in actual if composite(key)]
                    assert len(keys) == 288
                    error = {}
                    for component, field in enumerate(fields):
                        scale = max(1e-12, max(abs(native[step][key][component]) for key in keys))
                        absolute = max(abs(actual[key][component]-native[step][key][component]) for key in keys)
                        worst = max(keys, key=lambda key: abs(actual[key][component]-native[step][key][component]))
                        error[field] = {'absolute': absolute, 'scaled': absolute/scale, 'worst_cell': worst}
                    volumes = {key: math.prod(spacing)/(8 if key[0] == 1 else 1) for key in keys}
                    reference_mass = math.fsum(native[step][key][0]*volumes[key] for key in keys)
                    actual_mass = math.fsum(actual[key][0]*volumes[key] for key in keys)
                    initial_mass = math.fsum(native[0][key][0]*volumes[key] for key in keys)
                    result['steps'].append({'step': step, 'fields': error,
                                           'reference_mass_relative_change': (reference_mass-initial_mass)/initial_mass,
                                           'fortran_mass_relative_change': (actual_mass-initial_mass)/initial_mass})
                    field_limit = limits['initial_scaled' if step == 0 else 'field_scaled']
                    for field, value in error.items():
                        if value['scaled'] > field_limit:
                            result['violations'].append({'step': step, 'field': field,
                                                         'scaled': value['scaled'], 'limit': field_limit})
                    mass_difference = abs(actual_mass-reference_mass)/initial_mass
                    if mass_difference > limits['inventory_difference_scaled']:
                        result['violations'].append({'step': step, 'field': 'composite_mass',
                                                     'scaled': mass_difference,
                                                     'limit': limits['inventory_difference_scaled']})
                evolution = max(abs(native[-1][key][component]-native[0][key][component]) /
                                max(1e-12, max(abs(native[0][k][component]) for k in keys))
                                for component in range(7) for key in keys)
                assert evolution >= limits['minimum_evolution_scaled'], 'nontrivial spatial evolution required'
                result['maximum_evolution_scaled'] = evolution
                result['status'] = 'failed' if result['violations'] else 'passed'
                save()
                if result['violations'] and not args.keep_going:
                    raise AssertionError((name, result['violations']))
                print(name, {field: value['scaled'] for field, value in result['steps'][-1]['fields'].items()}, flush=True)
        checks['matrix_complete'] = True
        checks['failed_cases'] = [case['case'] for case in checks['cases'] if case['violations']]
        if checks['failed_cases']:
            raise AssertionError(('comparison gates failed', checks['failed_cases']))
        checks['status'] = 'passed'
    except Exception as error:
        checks['status'] = 'failed'
        checks['error'] = repr(error)
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
