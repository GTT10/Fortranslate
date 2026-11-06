#!/usr/bin/env python3
"""Exact serial tree continuation across scheduled topology changes."""
import argparse
from pathlib import Path
import subprocess
import re
from amr_tree_checkpoint_test_support import checkpoint_records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fixed', type=Path, required=True)
    parser.add_argument('--selected', type=Path, required=True)
    parser.add_argument('--case-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--subcycling', action='store_true')
    parser.add_argument('--pressure-outlet', action='store_true')
    parser.add_argument('--mach-inlet', action='store_true')
    args = parser.parse_args()
    case = args.case_dir.resolve()
    work = args.work_dir.resolve()
    original = (case / 'hotspot.nml').read_text()
    if args.pressure_outlet:
        original = original.replace("  boundary_condition=",
                                    "  outlet_pressure(2)=102325\n  boundary_condition=")
    if args.mach_inlet:
        original = original.replace("boundary_faces='outflow','outflow'", "boundary_faces='inflow','outflow'")
        original = original.replace('  boundary_condition=',
                                    '  inflow_mach(1)=0.2, inflow_pressure(1)=0, inflow_temperature(1)=900\n'
                                    '  boundary_condition=')
    time_controls = '\n&amr_tree_timestepping subcycling=.true. /\n'
    if args.subcycling:
        original += time_controls
    for kind, binary in [('fixed', args.fixed), ('selected', args.selected)]:
        directory = work / kind
        directory.mkdir(parents=True, exist_ok=True)
        config = original
        if kind == 'selected':
            config = config.replace("thermo_model='elementary'", "thermo_model='selected'")
            config = config.replace('  output_file=', "  composition_count=7\n"
                                    "  composition_species='H2','H','O','O2','OH','H2O','N2'\n"
                                    "  composition_mole_fractions=.29570,.00001,.00001,.14784,.00001,0.0,.55643\n"
                                    '  output_file=')
        inputs = directory / 'input.nml'
        outputs = [directory / 'tree_hotspot.csv', directory / 'tree_hotspot.csv.budget.csv']

        def run(label, content, success=True, schedule=None):
            inputs.write_text(content)
            result = subprocess.run([str(binary.resolve()), str(inputs), str(case / 'branched.dat'),
                                     str(schedule or case / 'regrid.dat')], cwd=directory,
                                    capture_output=True, text=True, timeout=180)
            log = result.stdout + result.stderr
            (directory / f'{label}.log').write_text(log)
            assert (result.returncode == 0) == success, log
            return log

        continuous_log = run('continuous', config)
        continuous_steps = re.search(r'Steps: (\d+)', continuous_log).group(1)
        baseline = [p.read_bytes() for p in outputs]
        for interval in ([1] if args.subcycling or args.pressure_outlet or args.mach_inlet else [1, 2, 3]):
            controls = (f"\n&amr_tree_checkpoint checkpoint_file='tree.chk', checkpoint_interval_steps={interval},"
                        ' stop_after_checkpoint=.true. /\n')
            log = run(f'stop{interval}', config + controls)
            assert f'Stopped after checkpoint at step {interval}' in log, log
            assert [p.read_bytes() for p in outputs] == baseline
            restart = "\n&amr_tree_checkpoint restart_file='tree.chk' /\n"
            log = run(f'restart{interval}', config + restart)
            assert 'Regrids: 3' in log and f'Steps: {continuous_steps} ' in log, log
            assert [p.read_bytes() for p in outputs] == baseline
        saved = (directory / 'tree.chk').read_bytes()
        if args.pressure_outlet:
            run('changed_pressure', config.replace('outlet_pressure(2)=102325',
                                                   'outlet_pressure(2)=103325') + restart, False)
            run('removed_pressure', config.replace('outlet_pressure(2)=102325',
                                                   'outlet_pressure(2)=0') + restart, False)
        if args.mach_inlet:
            for label, old, new in [('mach', 'inflow_mach(1)=0.2', 'inflow_mach(1)=0.3'),
                                    ('inlet_pressure', 'inflow_pressure(1)=0', 'inflow_pressure(1)=120000'),
                                    ('inlet_temperature', 'inflow_temperature(1)=900', 'inflow_temperature(1)=901')]:
                run('changed_' + label, config.replace(old, new) + restart, False)
        if args.subcycling:
            run('changed_timestepping', config.replace(time_controls, '') + restart, False)
        run('changed_physics', config.replace('initial_temperature=1100.0', 'initial_temperature=1110.0') + restart, False)
        altered = directory / 'changed_schedule.dat'
        altered.write_text('3\n2.0e-10 ' + repr(str(case / 'moved_ratio3.dat')) + '\n5e-10 ' +
                           repr(str(case / 'root.dat')) + '\n1e-9 ' + repr(str(case / 'branched.dat')) + '\n')
        run('changed_schedule', config + restart, False, altered)
        (directory / 'tree.chk').write_bytes(saved[:saved.rfind(b'END_AMR_TREE_CHECKPOINT')])
        run('truncated', config + restart, False)
        (directory / 'tree.chk').write_bytes(saved)
        # Locate the scalar metadata after the exact byte-counted context and geometry header.
        prefix, records, ledger_index, _ = checkpoint_records(saved)
        original_metadata = records[1]
        fields = original_metadata.split()
        fields[-1] = b'0'  # final-time topology cannot have zero consumed schedule events
        records[1] = b' '.join(fields) + b'\n'
        (directory / 'tree.chk').write_bytes(prefix + b''.join(records))
        run('wrong_schedule_position', config + restart, False)
        records[1] = original_metadata
        ledger = records[ledger_index].split()
        ledger[len(ledger) // 2] = b'1.0'  # finite but inconsistent mass boundary ledger
        records[ledger_index] = b' '.join(ledger) + b'\n'
        (directory / 'tree.chk').write_bytes(prefix + b''.join(records))
        run('wrong_ledger', config + restart, False)
        (directory / 'tree.chk').write_bytes(saved)
        run('alias', config + "\n&amr_tree_checkpoint checkpoint_file='./tree_hotspot.csv', "
            'checkpoint_interval_steps=1 /\n', False)
        run('malformed', config + "\n&amr_tree_checkpoint checkpoint_file='tree.chk', checkpoint_interval_steps=1\n", False)
        assert [p.read_bytes() for p in outputs] == baseline
    print('PASS: fixed/selected exact tree continuation and input/failure guards; subcycling=', args.subcycling,
          'pressure outlet=', args.pressure_outlet, 'Mach inlet=', args.mach_inlet)


if __name__ == '__main__':
    main()
