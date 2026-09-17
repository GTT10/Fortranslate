#!/usr/bin/env python3
"""Check rank-neutral static tree restart, serial interchange and rejected reads."""
import argparse
import os
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('fixed', 'selected', 'serial-fixed', 'serial-selected', 'case-dir', 'workdir'):
        parser.add_argument('--' + option, type=Path, required=True)
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--numproc-flag', default='-n')
    args = parser.parse_args()
    case = args.case_dir.resolve()
    template = (case / 'hotspot.nml').read_text().replace('cfl=0.2', 'cfl=0.002')
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

        def run(label, content, ranks=2, error=None):
            inputs.write_text(content)
            binary = getattr(args, kind if ranks else 'serial_' + kind).resolve()
            command = [str(binary), str(inputs), str(case / 'branched.dat')]
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
        baseline = [p.read_bytes() for p in outputs]
        save = "\n&amr_tree_checkpoint checkpoint_file='tree.chk', checkpoint_interval_steps=1, stop_after_checkpoint=.true. /\n"
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
    print('PASS: fixed/selected MPI 2-to-3 restart, serial interchange, exact output and rejected checkpoint guards')


if __name__ == '__main__':
    main()
