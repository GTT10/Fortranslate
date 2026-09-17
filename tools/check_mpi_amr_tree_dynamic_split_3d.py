#!/usr/bin/env python3
"""Exercise recursive chemistry/transport/hydro across regrid, rebalance and V2 restart."""
import argparse
import os
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('fixed', 'selected', 'serial-fixed', 'serial-selected', 'case-dir', 'workdir'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--numproc-flag', default='-n')
    args = parser.parse_args()
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    topology = work / 'branched.dat'
    moved = work / 'moved.dat'
    root = work / 'root.dat'
    schedule = work / 'schedule.dat'
    topology.write_text('3\n1 1 1 1 1 1 1 2\n1 3 3 1 1 1 1 2\n2 1 1 1 1 1 1 2\n')
    moved.write_text('2\n1 2 2 1 1 1 1 3\n2 1 1 1 1 1 1 2\n')
    root.write_text('0\n')
    # Every regrid is followed by a nonzero interval, including the last refinement.
    schedule.write_text(f"3\n2.5e-10 '{moved}'\n5.0e-10 '{root}'\n7.5e-10 '{topology}'\n")
    template = (args.case_dir / 'hotspot.nml').read_text()
    template = template.replace("reconstruction='pcm'", "reconstruction='characteristic_plm'")
    template += '\n&amr_tree_timestepping subcycling=.true. /\n'
    env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')
    reference = None
    for kind in ('fixed', 'selected'):
        config = template
        if kind == 'selected':
            config = config.replace("thermo_model='elementary'", "thermo_model='selected'")
            config = config.replace('  output_file=', "  composition_count=7\n"
                                    "  composition_species='H2','H','O','O2','OH','H2O','N2'\n"
                                    "  composition_mole_fractions=.29570,.00001,.00001,.14784,.00001,0.0,.55643\n"
                                    '  output_file=')
        folder = work / kind
        folder.mkdir(exist_ok=True)
        inputs = folder / 'input.nml'
        checkpoint = folder / 'tree.chk'
        outputs = [folder / 'tree_hotspot.csv', folder / 'tree_hotspot.csv.budget.csv']

        def run(label, content, ranks):
            if ranks:
                policy = 'weighted' if label == 'stop' else 'round_robin'
                content += (f"\n&amr_tree_distribution partitioning='{policy}', "
                            'rebalance_interval_steps=1, rebalance_minimum_gain=0.0 /\n')
            inputs.write_text(content)
            binary = getattr(args, kind if ranks else 'serial_' + kind).resolve()
            command = [str(binary), str(inputs), str(topology), str(schedule)]
            if ranks:
                command = [args.mpiexec, args.numproc_flag, str(ranks), *command]
            result = subprocess.run(command, cwd=folder, env=env, capture_output=True, text=True, timeout=240)
            log = result.stdout + result.stderr
            (folder / (label + '.log')).write_text(log)
            assert result.returncode == 0, log
            return log

        log = run('serial', config, 0)
        assert re.search(r'Steps:\s+4\b', log), log
        assert 'Regrids: 3' in log, log
        baseline = [p.read_bytes() for p in outputs]
        if reference is None:
            reference = baseline
        else:
            assert baseline == reference, 'Fixed/selected serial outputs differ'
        save = ("\n&amr_tree_checkpoint checkpoint_file='tree.chk', checkpoint_indexed=.true., "
                'checkpoint_interval_steps=1, stop_after_checkpoint=.true. /\n')
        log = run('stop', config + save, 2)
        assert 'Stopped after checkpoint at step 1' in log, log
        saved = checkpoint.read_bytes()
        assert saved.startswith(b'PELEF_AMR_TREE_3D_CHECKPOINT_V2\n'), saved[:80]
        log = run('restart', config + "\n&amr_tree_checkpoint restart_file='tree.chk' /\n", 3)
        assert re.search(r'Steps:\s+4\b', log) and 'Regrids: 3' in log, log
        assert log.count('Rebalance step:') == 3, log
        assert [p.read_bytes() for p in outputs] == baseline, f'{kind}: restarted output mismatch'
        run('serial_stop', config + save, 0)
        assert checkpoint.read_bytes() == saved, f'{kind}: indexed checkpoint byte mismatch'
        assert [p.read_bytes() for p in outputs] == baseline, 'Checkpoint stop overwrote final output'
    print('PASS: fixed/selected dynamic recursive full split, rebalancing, V2 2-to-3 restart and exact serial bytes')


if __name__ == '__main__':
    main()
