#!/usr/bin/env python3
"""Check rank-neutral tree restart, optional schedules and serial interchange."""
import argparse
import os
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('fixed', 'selected', 'serial-fixed', 'serial-selected', 'case-dir', 'workdir'):
        parser.add_argument('--' + option, type=Path, required=True)
    parser.add_argument('--scheduled', action='store_true')
    parser.add_argument('--subcycling', action='store_true')
    parser.add_argument('--indexed', action='store_true')
    parser.add_argument('--rebalance', action='store_true')
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--numproc-flag', default='-n')
    args = parser.parse_args()
    case = args.case_dir.resolve()
    template = (case / 'hotspot.nml').read_text()
    if not args.scheduled:
        template = template.replace('cfl=0.2', 'cfl=0.002')
    topology = case / 'branched.dat'
    if args.subcycling:
        if args.scheduled:
            # Isolate recursive geometry/event/restart integration; full split is checked separately.
            template = template.replace('chemistry_enabled=.true., transport_enabled=.true.',
                                        'chemistry_enabled=.false., transport_enabled=.false.')
            template = template.replace('wall_velocity(1,3)=2.0, wall_velocity(1,4)=2.0', 'wall_velocity=0.0')
        else:
            args.workdir.mkdir(parents=True, exist_ok=True)
            topology = (args.workdir / 'subcycle_tree.dat').resolve()
            topology.write_text('1\n1 1 1 1 1 1 1 2\n')
        template += '\n&amr_tree_timestepping subcycling=.true. /\n'
    env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')
    for kind in ('fixed', 'selected'):
        config = template
        if kind == 'selected':
            config = config.replace("thermo_model='elementary'", "thermo_model='selected'")
            config = config.replace('  output_file=', "  composition_count=7\n"
                                    "  composition_species='H2','H','O','O2','OH','H2O','N2'\n"
                                    "  composition_mole_fractions=.29570,.00001,.00001,.14784,.00001,0.0,.55643\n"
                                    '  output_file=')
        folder = (args.workdir / kind).resolve()
        folder.mkdir(parents=True, exist_ok=True)
        inputs = folder / 'input.nml'
        checkpoint = folder / 'tree.chk'
        outputs = [folder / name for name in ('tree_hotspot.csv', 'tree_hotspot.csv.budget.csv')]

        def run(label, content, ranks=2, error=None, schedule=None):
            if args.rebalance and ranks and '&amr_tree_distribution' not in content:
                cadence = 0 if label == 'restart_mpi_3' else 1
                content += (f'\n&amr_tree_distribution rebalance_interval_steps={cadence}, '
                            'rebalance_minimum_gain=0.0 /\n')
            inputs.write_text(content)
            binary = getattr(args, kind if ranks else 'serial_' + kind).resolve()
            command = [str(binary), str(inputs), str(topology)]
            if args.scheduled:
                command.append(str(schedule or case / 'regrid.dat'))
            if ranks:
                command = [args.mpiexec, args.numproc_flag, str(ranks), *command]
            result = subprocess.run(command, cwd=folder, env=env, capture_output=True, text=True, timeout=180)
            log = result.stdout + result.stderr
            (folder / (label + '.log')).write_text(log)
            if error is None:
                assert result.returncode == 0, log
            else:
                assert result.returncode != 0 and error in log, log
            return log

        def compare():
            assert [p.read_bytes() for p in outputs] == baseline, f'{kind}: continuation output mismatch'

        log = run('continuous', config)
        steps = int(re.search(r'Steps:\s+(\d+)', log).group(1))
        assert steps > 1, 'Fixture must continue beyond saved step'
        if args.rebalance:
            assert log.count('Rebalance step:') == steps, log
        baseline = [p.read_bytes() for p in outputs]
        if args.scheduled:
            assert 'Regrids: 3' in log
        if args.scheduled or args.subcycling or args.rebalance:
            serial_log = run('continuous_serial', config, 0)
            assert int(re.search(r'Steps:\s+(\d+)', serial_log).group(1)) == steps
            compare()
        save = "\n&amr_tree_checkpoint checkpoint_file='tree.chk', checkpoint_interval_steps=1, stop_after_checkpoint=.true. /\n"
        if args.indexed:
            save = save.replace('stop_after_checkpoint=.true.', 'stop_after_checkpoint=.true., checkpoint_indexed=.true.')
        restart = "\n&amr_tree_checkpoint restart_file='tree.chk' /\n"
        log = run('stop_mpi_2', config + save)
        assert 'Stopped after checkpoint at step 1' in log
        compare()
        saved = checkpoint.read_bytes()
        run('restart_mpi_3', config + restart, 3)
        compare()
        run('restart_serial', config + restart, 0)
        compare()
        run('stop_serial', config + save, 0)
        assert checkpoint.read_bytes() == saved, 'Serial/MPI checkpoint bytes differ'
        run('restart_mpi_1', config + restart, 1)
        compare()
        if args.scheduled:
            # The second event has root-only geometry: two ranks own no cells on reload.
            run('stop_root_only', config + save.replace('interval_steps=1', 'interval_steps=2'))
            run('restart_root_only', config + restart, 3)
            compare()
            checkpoint.write_bytes(saved)
            changed_schedule = folder / 'changed_schedule.dat'
            changed_schedule.write_text("3\n2.0e-10 " + repr(str(case / 'moved_ratio3.dat')) +
                                        "\n5e-10 " + repr(str(case / 'root.dat')) +
                                        "\n1e-9 " + repr(str(topology)) + "\n")
            run('changed_schedule', config + restart, 3, 'Invalid MPI tree restart checkpoint/context', changed_schedule)
        if args.rebalance:
            for option, value, error in [('rebalance_interval_steps', '-1', 'Invalid rebalance interval'),
                                         ('rebalance_minimum_gain', '1.0', 'Invalid rebalance gain')]:
                run('bad_' + option, config + f'\n&amr_tree_distribution {option}={value} /\n', 3, error)
                compare()
        if args.subcycling:
            run('malformed_subcycling', config + '\n&amr_tree_timestepping subcycling=invalid /\n',
                3, 'Invalid tree timestepping configuration')
            compare()
            run('changed_subcycling', config.replace('subcycling=.true.', 'subcycling=.false.') + restart,
                3, 'Invalid MPI tree restart checkpoint/context')
            compare()
        if args.indexed:
            assert saved.startswith(b'PELEF_AMR_TREE_3D_CHECKPOINT_V2\n')
            first, length, rest = saved.split(b'\n', 2)
            context_size = int(length)
            context, records = rest[:context_size], rest[context_size:].splitlines(keepends=True)
            nodes = int(records[0].split()[0])
            records[nodes + 2] = b'1\n'  # First relative node offset must be zero.
            checkpoint.write_bytes(first + b'\n' + length + b'\n' + context + b''.join(records))
            run('wrong_index', config + restart, 3, 'Invalid MPI tree restart checkpoint/context')
            compare()
            checkpoint.write_bytes(saved)
        run('changed_physics', config.replace('initial_temperature=1100.0', 'initial_temperature=1110.0') + restart,
            3, 'Invalid MPI tree restart checkpoint/context')
        checkpoint.write_bytes(saved[:saved.rfind(b'END_AMR_TREE_CHECKPOINT')])
        run('truncated', config + restart, 3, 'Invalid MPI tree restart checkpoint/context')
        # Corrupt only the finite boundary ledger; the physics/context remain unchanged.
        first, length, rest = saved.split(b'\n', 2)
        count = int(length)
        context, records = rest[:count], rest[count:].splitlines(keepends=True)
        ledger = records[2].split()
        ledger[len(ledger) // 2] = b'1.0'
        records[2] = b' '.join(ledger) + b'\n'
        checkpoint.write_bytes(first + b'\n' + length + b'\n' + context + b''.join(records))
        run('wrong_ledger', config + restart, 3, 'Tree restart boundary ledger mismatch')
        checkpoint.write_bytes(saved)
        run('alias', config + save.replace("'tree.chk'", "'./tree_hotspot.csv'"), 3,
            'MPI tree checkpoint aliases input/output')
        run('malformed', config + '\n&amr_tree_checkpoint checkpoint_interval_steps=1\n', 3,
            'Invalid MPI tree checkpoint configuration')
        compare()
    print('PASS: fixed/selected MPI restart, serial interchange, exact output and rejected checkpoint guards; scheduled=', args.scheduled, 'subcycling=', args.subcycling)


if __name__ == '__main__':
    main()
