#!/usr/bin/env python3
"""MPI automatic tree evolution, serial parity and changed-rank continuation."""
import argparse
import csv
import os
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--subcycling', action='store_true')
    parser.add_argument('--fixed', type=Path, required=True)
    parser.add_argument('--selected', type=Path, required=True)
    parser.add_argument('--serial-fixed', type=Path, required=True)
    parser.add_argument('--serial-selected', type=Path, required=True)
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--numproc-flag', default='-n')
    parser.add_argument('--case-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    case = args.case_dir.resolve()
    work = args.work_dir.resolve()
    source = (case / 'hotspot.nml').read_text().replace('cfl=0.2', 'cfl=0.001')
    source = source.replace('chemistry_enabled=.true., transport_enabled=.true.',
                            'chemistry_enabled=.false., transport_enabled=.false.')
    source = source.replace('hotspot_center_x=0.002, hotspot_center_y=0.002, hotspot_center_z=0.002',
                            'hotspot_center_x=0.0005, hotspot_center_y=0.0005, hotspot_center_z=0.0005')
    source = source.replace('wall_velocity(1,3)=2.0, wall_velocity(1,4)=2.0', 'wall_velocity=0.0')
    if args.subcycling:
        source = source.replace("reconstruction='pcm'", "reconstruction='characteristic_plm'")
        source += '\n&amr_tree_timestepping subcycling=.true. /\n'
    adaptation = """
&amr_tree_adaptation
 regrid_interval_steps=1, refinement_levels=2, refinement_ratios=2,2
 temperature_thresholds=1150,1150, tile_shape=1,1,1, tag_buffer=0, maximum_nodes=64
/
"""
    outputs_by_kind = []
    for kind, binary in [('fixed', args.fixed), ('selected', args.selected)]:
        directory = work / kind
        directory.mkdir(parents=True, exist_ok=True)
        config = source
        if kind == 'selected':
            config = config.replace("thermo_model='elementary'", "thermo_model='selected'")
            config = config.replace('  output_file=', "  composition_count=7\n"
                                    "  composition_species='H2','H','O','O2','OH','H2O','N2'\n"
                                    "  composition_mole_fractions=.29570,.00001,.00001,.14784,.00001,0.0,.55643\n"
                                    '  output_file=')
        inputs = directory / 'input.nml'
        output = directory / 'tree_hotspot.csv'
        budget = directory / 'tree_hotspot.csv.budget.csv'

        def run(label, content, success=True, ranks=2):
            inputs.write_text(content)
            executable = binary if ranks else getattr(args, 'serial_' + kind)
            command = [str(executable.resolve()), str(inputs), str(case / 'root.dat')]
            if ranks:
                command = [args.mpiexec, args.numproc_flag, str(ranks), *command]
            result = subprocess.run(command, cwd=directory, text=True, capture_output=True, timeout=180,
                                    env=dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1'))
            log = result.stdout + result.stderr
            (directory / f'{label}.log').write_text(log)
            assert (result.returncode == 0) == success, log
            return log

        log = run('continuous', config + adaptation)
        steps = int(re.search(r'Steps:\s+(\d+)', log)[1])
        assert steps > 2 and log.count('Automatic regrid check at step') == steps, log
        data = list(csv.DictReader(output.open()))
        assert max(int(row['level']) for row in data) == 2
        assert abs(sum(float(row['dx']) * float(row['dy']) * float(row['dz']) for row in data) - .004**3) < 1e-20
        baseline = [output.read_bytes(), budget.read_bytes()]
        serial_log = run('serial_continuous', config + adaptation, ranks=0)
        assert int(re.search(r'Steps:\s+(\d+)', serial_log)[1]) == steps
        assert [output.read_bytes(), budget.read_bytes()] == baseline
        stop = "\n&amr_tree_checkpoint checkpoint_file='auto.chk', checkpoint_interval_steps=1, stop_after_checkpoint=.true. /\n"
        # Exercise both optional-group orders.
        log = run('stop', config + stop + adaptation)
        assert 'Stopped after checkpoint at step 1' in log and 'Automatic regrid check at step 1' in log
        assert [output.read_bytes(), budget.read_bytes()] == baseline
        restart = "\n&amr_tree_checkpoint restart_file='auto.chk' /\n"
        log = run('restart', config + adaptation + restart, ranks=3)
        assert 'Automatic regrid check at step 1 ' not in log
        assert int(re.search(r'Regrids:\s+(\d+)', log)[1]) == steps
        assert [output.read_bytes(), budget.read_bytes()] == baseline
        saved = (directory / 'auto.chk').read_bytes()
        run('serial_restart', config + adaptation + restart, ranks=0)
        assert [output.read_bytes(), budget.read_bytes()] == baseline
        run('serial_stop', config + stop + adaptation, ranks=0)
        assert (directory / 'auto.chk').read_bytes() == saved
        # Corrupt only the adaptation counter, keeping context and state intact.
        first, length, rest = saved.split(b'\n', 2)
        count = int(length)
        context, records = rest[:count], rest[count:].splitlines(keepends=True)
        metadata = records[1].split()
        metadata[-1] = b'0'
        records[1] = b' '.join(metadata) + b'\n'
        (directory / 'auto.chk').write_bytes(first + b'\n' + length + b'\n' + context + b''.join(records))
        run('wrong_adaptation_counter', config + adaptation + restart, False, ranks=3)
        (directory / 'auto.chk').write_bytes(saved)
        if args.subcycling:
            run('changed_subcycling', config.replace('subcycling=.true.', 'subcycling=.false.') + adaptation + restart,
                False, ranks=3)
        run('changed_threshold', config + adaptation.replace('1150,1150', '1160,1160') + restart, False)
        run('changed_cadence', config + adaptation.replace('regrid_interval_steps=1', 'regrid_interval_steps=2') + restart, False)
        run('node_budget', config + adaptation.replace('maximum_nodes=64', 'maximum_nodes=2'), False)
        assert [output.read_bytes(), budget.read_bytes()] == baseline
        outputs_by_kind.append(baseline)
    assert outputs_by_kind[0] == outputs_by_kind[1]
    print('PASS: MPI automatic adaptation; exact serial parity and 2-to-3 restart; context/counter/budget rejection')


if __name__ == '__main__':
    main()
