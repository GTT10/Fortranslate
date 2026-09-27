#!/usr/bin/env python3
"""Check MPI LES SDC activation and changed-rank checkpoint continuation."""
import argparse
import os
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--serial', type=Path)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--recursive', action='store_true')
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    template = (args.source / 'cases/amr_tree_reactive_3d/hotspot.nml').read_text()
    template = template.replace("reconstruction='pcm'", "time_integrator='sdc', sdc_iterations=2, reconstruction='characteristic_ppm'")
    template = template.replace("riemann_solver='hllc'", "riemann_solver='pelec'")
    template = template.replace('final_time=1.0e-9, cfl=0.2', 'final_time=3.0e-10, cfl=0.0002')
    template = template.replace('chemistry_enabled=.true.', 'chemistry_enabled=.false.')
    if args.recursive:
        template += '\n&amr_tree_timestepping subcycling=.true. /\n'
    topology = work / 'tree.dat'
    topology.write_text('1\n1 2 2 2 2 2 2 2\n')
    env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')
    les = "\n&amr_tree_les enabled=.true., model='smagorinsky', cs=.17 /\n"

    def run(label, content, ranks=2, error=None):
        inp = work / 'input.nml'
        inp.write_text(content)
        command = [str((args.binary if ranks else args.serial).resolve()), str(inp), str(topology)]
        if ranks:
            command = [args.mpiexec, '-n', str(ranks), *command]
        result = subprocess.run(command,
                                cwd=work, env=env, capture_output=True, text=True, timeout=180)
        log = result.stdout + result.stderr
        (work / (label + '.log')).write_text(log)
        if error:
            assert result.returncode != 0 and error in log, log
        else:
            assert result.returncode == 0, log
        return log

    def outputs():
        return [(work / name).read_bytes() for name in ('tree_hotspot.csv', 'tree_hotspot.csv.budget.csv')]

    run('disabled', template)
    disabled = outputs()
    log = run('continuous', template + les)
    assert int(re.search(r'Steps:\s+(\d+)', log).group(1)) > 1, log
    baseline = outputs()
    assert baseline[0] != disabled[0], 'LES must affect the gas state'
    if args.serial:
        run('serial_continuous', template + les, 0)
        assert outputs() == baseline, 'Serial/MPI application outputs differ'
    save = "\n&amr_tree_checkpoint checkpoint_file='tree.chk', checkpoint_interval_steps=1, stop_after_checkpoint=.true. /\n"
    restart = "\n&amr_tree_checkpoint restart_file='tree.chk' /\n"
    log = run('stop', template + les + save)
    assert 'Stopped after checkpoint at step 1' in log, log
    run('restart_4', template + les + restart, 4)
    assert outputs() == baseline, 'Changed-rank continuation differs'
    if args.serial:
        run('serial_restart', template + les + restart, 0)
        assert outputs() == baseline, 'MPI-to-serial continuation differs'
        run('serial_stop', template + les + save, 0)
        run('serial_to_mpi', template + les + restart, 4)
        assert outputs() == baseline, 'Serial-to-MPI continuation differs'
    run('changed_les', template + les.replace('cs=.17', 'cs=.18') + restart, 4,
        'Invalid MPI tree restart checkpoint/context')
    run('removed_les', template + restart, 4, 'Invalid MPI tree restart checkpoint/context')
    if args.serial:
        print('PASS: serial/MPI outputs and checkpoint interchange are exact')
    print('PASS: active LES changes state; 2-to-4-rank restart is exact; changed/removed LES rejected')


if __name__ == '__main__':
    main()
