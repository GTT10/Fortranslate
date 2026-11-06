#!/usr/bin/env python3
"""Reacting H2/O2/N2 EB-AMR lifecycle, inventories and changed-rank restart.

This is a bounded application gate, not fixed-PeleC or ignition validation.
Unlike the pure-N2 held-energy fixture it requires measurable chemistry.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
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
    for name in ('serial', 'mpi', 'mpiexec', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--subcycling', action='store_true')
    parser.add_argument('--smoke-only', action='store_true')
    parser.add_argument('--full-cvode', action='store_true',
                        help='Exercise full-coordinate CVODE and reactor restart identity')
    parser.add_argument('--hot-cvode', action='store_true',
                        help='Use the high-temperature CVODE case and only the essential restart path')
    parser.add_argument('--pinned-eos-convergence', action='store_true',
                        help='Exercise the opt-in converged caloric inverse and checkpoint identity')
    parser.add_argument('--separate-momentum', action='store_true',
                        help='Exercise opt-in momentum limiter and checkpoint identity')
    parser.add_argument('--consistent-velocity', action='store_true')
    parser.add_argument('--paired-pressure', action='store_true')
    parser.add_argument('--connected-directions', action='store_true')
    args = parser.parse_args()
    if args.full_cvode:
        args.hot_cvode = True
    if args.connected_directions:
        args.separate_momentum = True
    if args.paired_pressure:
        args.consistent_velocity = True
    if args.consistent_velocity:
        args.separate_momentum = True
    momentum_key = 'consistent_velocity' if args.consistent_velocity else 'separate_momentum'
    mismatch_message = 'Restart velocity policy mismatch' if args.consistent_velocity else 'Restart momentum policy mismatch'
    if args.paired_pressure:
        momentum_key = 'paired_pressure'
        mismatch_message = 'Restart pressure policy mismatch'
    if args.connected_directions:
        momentum_key = 'connected_directions'
        mismatch_message = 'Restart connected-direction policy mismatch'
    assert not args.pinned_eos_convergence or args.hot_cvode, 'Converged EOS check uses the hot CVODE case'
    assert not args.separate_momentum or args.hot_cvode, 'Momentum policy check uses hot CVODE case'
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='reacting_restart_', dir=args.work.resolve()))
    print(work, flush=True)
    bundle = Path(__file__).resolve().parents[2] / 'mechanisms/h2o2_full.json'
    thermo = json.loads(bundle.read_text())['thermo']
    report = {'scope': __doc__, 'status': 'running', 'subcycling': args.subcycling,
              'bundle_sha256': digest(bundle), 'runs': [], 'checks': {}}
    report['pinned_eos_convergence'] = args.pinned_eos_convergence
    report['separate_momentum'] = args.separate_momentum
    report['consistent_velocity'] = args.consistent_velocity
    report['paired_pressure'] = args.paired_pressure
    report['connected_directions'] = args.connected_directions
    report['full_cvode'] = args.full_cvode
    env = os.environ.copy()
    env['OMPI_MCA_rmaps_base_oversubscribe'] = '1'
    dt, end_time = 1.e-7, 3.e-7
    case_name = 'reacting_stiff.nml' if args.hot_cvode else 'reacting_restart.nml'
    case = Path(__file__).resolve().parents[2] / 'cases/reactive_eb_tree_3d' / case_name
    template = case.read_text().rstrip().removesuffix('/')
    report['case_sha256'] = digest(case)
    if args.subcycling:
        template = template.replace('subcycling=.false.', 'subcycling=.true.')
    if args.pinned_eos_convergence:
        template += '\n pinned_eos_extrapolation=.true.,pinned_eos_convergence=.true.,'
    if args.separate_momentum:
        template += '\n separate_momentum=.true.,'
    if args.consistent_velocity:
        template += '\n consistent_velocity=.true.,'
    if args.paired_pressure:
        template += '\n paired_pressure=.true.,\n'
    if args.connected_directions:
        template += '\n connected_directions=.true.,\n'

    def save():
        (work / 'checks.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')

    def run(name, ranks=0, extra='', chemistry=True, stop=False, restart=None,
            final_time=end_time, success=True, output_interval=1, eos_mode=None, momentum_mode=None,
            reactor_mode=None):
        directory = work / name
        directory.mkdir()
        text = template.replace('final_time=3.e-7', f'final_time={final_time:.17e}')
        if eos_mode is not None:
            assert args.pinned_eos_convergence
            text = text.replace('pinned_eos_convergence=.true.',
                                'pinned_eos_convergence=' + ('.true.' if eos_mode else '.false.'))
        if momentum_mode is not None:
            assert args.separate_momentum
            text = text.replace(momentum_key + '=.true.',
                                momentum_key + '=' + ('.true.' if momentum_mode else '.false.'))
        if not chemistry:
            text = text.replace('chemistry_enabled=.true.', 'chemistry_enabled=.false.')
            # The application intentionally rejects CVODE with reactions OFF.
            # In this control no reactor is called, so select native dispatch.
            text = text.replace("sdc_reactor='cvode'", "sdc_reactor='native'")
        if chemistry and args.full_cvode:
            text = text.replace("sdc_reactor='cvode'", "sdc_reactor='cvode_full'")
        if reactor_mode is not None:
            text = text.replace("sdc_reactor='cvode_full'", f"sdc_reactor='{reactor_mode}'")
        if stop:
            extra += "checkpoint_file='saved.chk',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,"
        if restart:
            extra += f"restart_file='{restart}',"
        if args.hot_cvode and chemistry:
            text += f'output_interval_steps={output_interval},'
        (directory / 'input.nml').write_text(text + extra + '\n/\n')
        binary = (args.mpi if ranks else args.serial).resolve()
        command = ([str(args.mpiexec.resolve()), '-n', str(ranks)] if ranks else [])
        command += [str(binary), 'input.nml']
        item = {'name': name, 'command': command, 'executable_sha256': digest(binary),
                'input_sha256': digest(directory / 'input.nml'), 'expected_success': success}
        report['runs'].append(item)
        save()
        start = time.monotonic()
        with (directory / 'run.log').open('w') as stream:
            try:
                completed = subprocess.run(command, cwd=directory, env=env, stdout=stream,
                                           stderr=subprocess.STDOUT, timeout=240)
                item['returncode'] = completed.returncode
            except subprocess.TimeoutExpired:
                item['timeout'] = True
                raise
            finally:
                item['seconds'] = time.monotonic() - start
                save()
        assert (completed.returncode == 0) == success, directory / 'run.log'
        if success:
            cells = rows(directory / 'result.csv')
            assert cells, directory
            for cell in cells:
                assert all(math.isfinite(float(value)) for value in cell.values()), directory
            log = (directory / 'run.log').read_text()
            expected_initializations = int(chemistry and restart is None)
            assert log.count('Initialized MOL reaction history for first accepted step') == expected_initializations
            item['initial_history_count'] = expected_initializations
            times = [float(value.replace('D', 'E')) for value in
                     re.findall(r'^Completed steps:\s*\d+\s+time:\s*(\S+)', log, re.MULTILINE)]
            target = dt if stop else final_time
            assert len(times) == 1 and abs(times[0]-target) <= 1.e-12*target, (name, times, target)
            if stop:
                assert 'Regrid stage: 1 nodes: 1' in (directory / 'run.log').read_text()
            elif final_time == end_time:
                assert {int(cell['level']) for cell in cells} == {0, 1}, name
                assert 'Regrid stage: 2 nodes: 2' in (directory / 'run.log').read_text()
        return directory

    def audit(directory, filename='result.csv'):
        cells = rows(directory / filename)
        budget = rows(directory / (filename + '.budget.csv'))
        fields = ['rho', 'rhou', 'rhov', 'rhow', 'rhoE'] + [f'rhoY_{s["name"]}' for s in thermo]
        assert len(fields) == len(budget)
        # Do not trust the driver's residual alone: independently integrate the
        # exposed-cell CSV and reconstruct each signed boundary/source balance.
        worst = 0.
        for field, row in zip(fields, budget):
            values = {key: float(value) for key, value in row.items() if key != 'component'}
            assert all(math.isfinite(v) for v in values.values()), (directory, field)
            initial, final = values['initial'], values['final']
            inventory = math.fsum(float(c[field]) * float(c['fluid_volume']) for c in cells)
            scale = max(1., abs(initial), abs(final))
            assert abs(inventory-final) <= 3.e-14 * scale, (directory, field, 'CSV inventory')
            changes = [values[key] for key in ('boundary', 'wall', 'chemistry', 'forcing', 'numerical_energy')]
            scale = max(scale, math.fsum(abs(v) for v in changes))
            residual = final - initial - math.fsum(changes)
            assert abs(residual) <= 5.e-10 * scale, (directory, field, residual/scale)
            assert abs(residual-values['residual']) <= 3.e-14 * scale, (directory, field)
            if field.startswith('rhoY_'):
                # A unit-sized floor would hide loss of a small accumulated
                # reaction source. Scale each species by its actual inventory
                # and exchanged amounts, including newly produced trace species.
                species_scale = max(abs(initial), abs(final), math.fsum(abs(v) for v in changes))
                assert abs(residual) <= 5.e-10 * species_scale, (directory, field, 'species balance')
            if field != 'rhoE':
                assert values['numerical_energy'] == 0
            worst = max(worst, abs(residual)/scale)
        for cell in cells:
            density = float(cell['rho'])
            t = float(cell['temperature'])
            assert density > 0 and t > 0
            assert all(float(cell[field]) >= 0 for field in fields[5:]), directory
            assert abs(math.fsum(float(cell[f]) for f in fields[5:])-density) <= 3.e-12*density
            energy = []
            for s in thermo:
                a = s['low_coefficients'] if t < s['temperature_mid'] else s['high_coefficients']
                e = 8314.46261815324/s['molecular_weight'] * (
                    (a[0]-1)*t + a[1]*t*t/2 + a[2]*t**3/3 + a[3]*t**4/4 + a[4]*t**5/5 + a[5])
                energy.append(float(cell[f'rhoY_{s["name"]}'])*e)
            actual, expected = float(cell['rho_internal_energy_J_m3']), math.fsum(energy)
            assert abs(actual-expected) <= 3.e-10*max(1., abs(actual), abs(expected)), (directory, 'NASA energy')
        species_budget = dict(zip((s['name'] for s in thermo), budget[5:]))
        element_errors = {}
        for element in ('H', 'O', 'N'):
            initial = math.fsum(float(species_budget[s['name']]['initial']) *
                                s['composition'].get(element, 0)/s['molecular_weight'] for s in thermo)
            changes = [float(species_budget[s['name']]['chemistry']) *
                       s['composition'].get(element, 0)/s['molecular_weight'] for s in thermo]
            error = abs(math.fsum(changes))/initial
            assert error <= 3.e-12, (directory, element, error)
            element_errors[element] = error
        h2 = species_budget['H2']
        reaction_fraction = -float(h2['chemistry'])/float(h2['initial'])
        return {'maximum_budget_residual': worst, 'chemical_element_residuals': element_errors,
                'hydrogen_reacted_fraction': reaction_fraction}

    try:
        serial = run('serial')
        report['checks']['serial'] = audit(serial)
        minimum_reaction = 1.e-5 if args.hot_cvode else 1.e-9
        assert report['checks']['serial']['hydrogen_reacted_fraction'] > minimum_reaction, 'chemistry must materially evolve'
        if args.smoke_only:
            report['status'] = 'smoke_passed_not_qualified'
            save()
            print('PASS: reacting serial smoke only', flush=True)
            return
        off = run('reaction_off', chemistry=False)
        report['checks']['reaction_off'] = audit(off)
        assert report['checks']['reaction_off']['hydrogen_reacted_fraction'] == 0
        on_cells, off_cells = rows(serial / 'result.csv'), rows(off / 'result.csv')
        assert len(on_cells) == len(off_cells)
        effect = max(abs(float(a['rhoY_H2'])-float(b['rhoY_H2']))/float(a['rho'])
                     for a, b in zip(on_cells, off_cells))
        assert effect > 1.e-11, 'reaction switch must change the exported gas field'
        report['checks']['reaction_off_intervention'] = effect
        parallel = run('mpi2', 2, extra='rebalance_interval_steps=1,rebalance_minimum_gain=0.,', output_interval=2)
        third = parallel if args.hot_cvode else run('mpi3', 3)
        stopped = run('stopped2', 2, stop=True, output_interval=2)
        restarted = run('resumed3', 3, restart=stopped / 'saved.chk')
        if args.full_cvode:
            checkpoint = stopped / 'saved.chk'
            original_hash = digest(checkpoint)
            for ranks in (0, 3):
                for reactor in ('native', 'cvode'):
                    refused = run(f'changed_reactor_{ranks}_{reactor}', ranks, restart=checkpoint,
                                  reactor_mode=reactor, success=False)
                    expected_error = ('Invalid distributed tree checkpoint record' if ranks else
                                      'Invalid tree checkpoint record')
                    assert expected_error in (refused / 'run.log').read_text()
                    assert not list(refused.glob('result.csv*'))
                    assert digest(checkpoint) == original_hash
            serial_resume = run('full_mpi_to_serial', restart=checkpoint)
            assert (serial_resume / 'result.csv').read_bytes() == (serial / 'result.csv').read_bytes()
            legacy = run('reduced_checkpoint', 2, stop=True, reactor_mode='cvode')
            legacy_checkpoint = legacy / 'saved.chk'
            legacy_hash = digest(legacy_checkpoint)
            for ranks in (0, 3):
                refused = run(f'reduced_to_full_{ranks}', ranks, restart=legacy_checkpoint, success=False)
                expected_error = ('Invalid distributed tree checkpoint record' if ranks else
                                  'Invalid tree checkpoint record')
                assert expected_error in (refused / 'run.log').read_text()
                assert not list(refused.glob('result.csv*'))
                assert digest(legacy_checkpoint) == legacy_hash
            report['checks']['full_reactor_restart_mismatch_serial_mpi'] = True
            report['checks']['full_reactor_restart_reverse_mismatch_serial_mpi'] = True
            report['checks']['full_mpi_to_serial'] = audit(serial_resume)
        if args.hot_cvode:
            # One changed-rank continuation is enough for this new stiff path;
            # serial interchange and malformed metadata belong to the existing
            # lifecycle tests rather than another full Cartesian test matrix.
            for directory in (parallel, restarted):
                assert (directory / 'result.csv').read_bytes() == (serial / 'result.csv').read_bytes(), directory
                report['checks'][directory.name] = audit(directory)
            for actual, saved in zip(rows(restarted / 'result.csv.budget.csv'),
                                     rows(stopped / 'result.csv.budget.csv')):
                assert actual['initial'] == saved['initial']
            # Intermediate fields expose accepted synchronized states, with
            # physical time and absolute steps across changed-rank restart.
            snapshot_checks = {}
            for directory, steps in ((serial, (0, 1, 2, 3)), (parallel, (0, 2, 3)),
                                     (stopped, (0, 1)), (restarted, (1, 2, 3))):
                expected = {f'result.csv.step{step}.csv' for step in steps}
                assert {p.name for p in directory.glob('result.csv.step*.csv')
                        if re.fullmatch(r'result\.csv\.step[0-9]+\.csv', p.name)} == expected
                for step in steps:
                    filename = f'result.csv.step{step}.csv'
                    times = rows(directory / (filename + '.time.csv'))
                    assert len(times) == 1
                    stamp = times[0]
                    assert int(stamp['step']) == step
                    assert abs(float(stamp['time_s']) - step * dt) <= 1.e-12 * dt
                    assert int(stamp['nodes']) == (1 if step == 1 else 2)
                    assert int(stamp['regrid_stage']) == min(step, 2)
                    assert (directory / filename).read_bytes() == (serial / filename).read_bytes()
                    snapshot_checks[f'{directory.name}/{step}'] = audit(directory, filename)
                final_snapshot = directory / f'result.csv.step{steps[-1]}.csv'
                assert final_snapshot.read_bytes() == (directory / 'result.csv').read_bytes()
                assert Path(str(final_snapshot) + '.budget.csv').read_bytes() == (
                    directory / 'result.csv.budget.csv').read_bytes()
            assert not list(off.glob('result.csv.step*')), 'default output must remain final-only'
            report['checks']['snapshots'] = snapshot_checks
            invalid = run('invalid_output_interval', output_interval=-1, success=False)
            assert 'Invalid output policy' in (invalid / 'run.log').read_text()
            assert not list(invalid.glob('result.csv*'))
            collision = run('snapshot_checkpoint_alias', extra="checkpoint_file='./result.csv.step0.csv',", success=False)
            assert 'Snapshot aliases input or checkpoint' in (collision / 'run.log').read_text()
            assert not list(collision.glob('result.csv*'))
            if args.separate_momentum:
                checkpoint = stopped / 'saved.chk'
                token = ('PELEF_EB_TREE_CONNECTED_DIRECTIONS_V1' if args.connected_directions else
                         'PELEF_EB_TREE_PAIRED_PRESSURE_V1' if args.paired_pressure else
                         'PELEF_EB_TREE_CONSISTENT_VELOCITY_V1' if args.consistent_velocity else 'PELEF_EB_TREE_SEPARATE_MOMENTUM_V1')
                assert checkpoint.read_text().startswith(token + '\n')
                original_hash = digest(checkpoint)
                refused = run('changed_momentum_policy', 3, restart=checkpoint,
                              momentum_mode=False, success=False)
                assert mismatch_message in (refused / 'run.log').read_text()
                assert not list(refused.glob('result.csv*'))
                assert digest(checkpoint) == original_hash
                legacy = run('legacy_policy_checkpoint', 2, stop=True, momentum_mode=False)
                legacy_checkpoint = legacy / 'saved.chk'
                assert token not in legacy_checkpoint.read_text()
                legacy_hash = digest(legacy_checkpoint)
                refused = run('legacy_to_separate_policy', 3, restart=legacy_checkpoint, success=False)
                assert mismatch_message in (refused / 'run.log').read_text()
                assert not list(refused.glob('result.csv*'))
                assert digest(legacy_checkpoint) == legacy_hash
                report['checks']['momentum_restart_mismatch_both_directions'] = True
                mixed = work / 'mixed_rank_momentum_policy'
                mixed.mkdir()
                for mode in (False, True):
                    text = template.replace(momentum_key + '=.true.',
                                            momentum_key + '=' + ('.true.' if mode else '.false.'))
                    (mixed / f'{int(mode)}.nml').write_text(text + '\n/\n')
                command = [str(args.mpiexec.resolve()), '-n', '1', str(args.mpi.resolve()), '0.nml',
                           ':', '-n', '1', str(args.mpi.resolve()), '1.nml']
                start = time.monotonic()
                with (mixed / 'run.log').open('w') as stream:
                    rejected = subprocess.run(command, cwd=mixed, env=env, stdout=stream,
                                              stderr=subprocess.STDOUT, timeout=30)
                report['runs'].append({'name': mixed.name, 'command': command, 'expected_success': False,
                                       'returncode': rejected.returncode, 'seconds': time.monotonic()-start,
                                       'executable_sha256': digest(args.mpi.resolve()),
                                       'input_sha256': {p.name: digest(p) for p in mixed.glob('*.nml')}})
                assert rejected.returncode != 0
                assert 'MPI EB-tree momentum policy mismatch' in (mixed / 'run.log').read_text()
                assert not list(mixed.glob('result.csv*'))
                report['checks']['mixed_rank_momentum_refused'] = True
            if args.pinned_eos_convergence:
                checkpoint = stopped / 'saved.chk'
                assert 'PELEF_FUEGO_EOS_CONVERGENCE_V1' in checkpoint.read_text()
                original_hash = digest(checkpoint)
                refused = run('changed_eos_precision', 3, restart=checkpoint, eos_mode=False, success=False)
                assert 'Invalid distributed tree checkpoint record' in (refused / 'run.log').read_text()
                assert not list(refused.glob('result.csv*'))
                assert digest(checkpoint) == original_hash
                report['checks']['changed_precision_restart_refused'] = True
                mixed = work / 'mixed_rank_eos_precision'
                mixed.mkdir()
                for mode in (False, True):
                    text = template.replace('pinned_eos_convergence=.true.',
                                            'pinned_eos_convergence=' + ('.true.' if mode else '.false.'))
                    (mixed / f'{int(mode)}.nml').write_text(text + '\n/\n')
                command = [str(args.mpiexec.resolve()), '-n', '1', str(args.mpi.resolve()), '0.nml',
                           ':', '-n', '1', str(args.mpi.resolve()), '1.nml']
                start = time.monotonic()
                with (mixed / 'run.log').open('w') as stream:
                    rejected = subprocess.run(command, cwd=mixed, env=env, stdout=stream,
                                              stderr=subprocess.STDOUT, timeout=30)
                report['runs'].append({'name': mixed.name, 'command': command, 'expected_success': False,
                                       'returncode': rejected.returncode, 'seconds': time.monotonic()-start,
                                       'executable_sha256': digest(args.mpi.resolve()),
                                       'input_sha256': {p.name: digest(p) for p in mixed.glob('*.nml')}})
                assert rejected.returncode != 0
                assert 'MPI EB-tree EOS policy mismatch' in (mixed / 'run.log').read_text()
                assert not list(mixed.glob('result.csv*'))
                report['checks']['mixed_rank_precision_refused'] = True
            report['status'] = 'passed'
            print('PASS: hot CVODE EB-AMR trajectories, reaction-off control and 2-to-3-rank restart', flush=True)
            return
        to_serial = run('mpi_to_serial', restart=stopped / 'saved.chk')
        serial_stop = run('serial_stop', stop=True)
        from_serial = run('serial_to_mpi', 3, restart=serial_stop / 'saved.chk')
        for directory in (parallel, third, restarted, to_serial, from_serial):
            assert (directory / 'result.csv').read_bytes() == (serial / 'result.csv').read_bytes(), directory
            report['checks'][directory.name] = audit(directory)
        # Cross-implementation restart must retain its source's initial
        # inventory, not replace it with the destination's differently grouped
        # serial/MPI reduction. All continuations above independently close
        # every species balance, including the small chemical source prefix.
        for directory, destination, checkpoint_run in (
                (restarted, third, stopped), (to_serial, serial, stopped),
                (from_serial, third, serial_stop)):
            actual = rows(directory / 'result.csv.budget.csv')
            native = rows(destination / 'result.csv.budget.csv')
            saved = rows(checkpoint_run / 'result.csv.budget.csv')
            for a, b, c in zip(actual, native, saved):
                assert a['initial'] == c['initial'], (directory, 'original inventory prefix')
                assert a['numerical_energy'] == b['numerical_energy'], (directory, 'energy correction prefix')
                assert a['internal_energy_pressure_work'] == b['internal_energy_pressure_work']
        restored = run('restore_only3', 3, restart=stopped / 'saved.chk', final_time=dt)
        assert (restored / 'result.csv').read_bytes() == (stopped / 'result.csv').read_bytes()
        assert (restored / 'result.csv.budget.csv').read_bytes() == (stopped / 'result.csv.budget.csv').read_bytes()
        run('changed_chemistry', 3, chemistry=False, restart=stopped / 'saved.chk', success=False)
        assert not (work / 'changed_chemistry/result.csv').exists()
        report['status'] = 'passed'
        print('PASS: reacting H2/O2/N2 EB-AMR, inventories, reaction-off control and changed-rank restart', flush=True)
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = repr(error)
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
