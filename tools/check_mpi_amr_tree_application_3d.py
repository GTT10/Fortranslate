#!/usr/bin/env python3
"""Compare fixed-mechanism static MPI tree application output with serial."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--serial', type=Path, required=True)
    parser.add_argument('--mpi', type=Path, required=True)
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--numproc-flag', default='-n')
    parser.add_argument('--case-dir', type=Path, required=True)
    parser.add_argument('--workdir', type=Path, required=True)
    args = parser.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    template = (args.case_dir / 'hotspot.nml').read_text()
    env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')

    def run(name: str, topology: str, ranks: int = 0, extra: str = '', fail: str = '') -> Path:
        folder = (args.workdir / name).resolve()
        folder.mkdir(parents=True, exist_ok=True)
        input_path = folder / 'input.nml'
        input_path.write_text(template + extra)
        if fail:
            (folder / 'tree_hotspot.csv').write_bytes(b'preserve-existing-output\n')
        command = [str(args.serial.resolve()), str(input_path), str((args.case_dir / topology).resolve())]
        if ranks:
            command = [args.mpiexec, args.numproc_flag, str(ranks), str(args.mpi.resolve()), *command[1:]]
        result = subprocess.run(command, cwd=folder, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, timeout=180, check=False)
        (folder / 'run.log').write_text(result.stdout)
        if fail:
            if result.returncode == 0 or fail not in result.stdout:
                raise RuntimeError(f'{name}: expected rejection {fail!r}; see {folder / "run.log"}')
            if (folder / 'tree_hotspot.csv').read_bytes() != b'preserve-existing-output\n':
                raise RuntimeError(f'{name}: rejected input overwrote output')
        else:
            if result.returncode or 'Steps:' not in result.stdout:
                raise RuntimeError(f'{name}: application failed; see {folder / "run.log"}')
        return folder

    def compare(reference: Path, candidate: Path) -> None:
        for filename in ('tree_hotspot.csv', 'tree_hotspot.csv.budget.csv'):
            expected = (reference / filename).read_bytes()
            if not expected or expected != (candidate / filename).read_bytes():
                raise RuntimeError(f'{candidate.name}: {filename} differs from serial')

    reference = run('branched_serial', 'branched.dat')
    for ranks in (2, 3):
        compare(reference, run(f'branched_mpi_{ranks}', 'branched.dat', ranks))
    reference = run('root_serial', 'root.dat')
    compare(reference, run('root_mpi_3', 'root.dat', 3))
    run('unsupported_checkpoint', 'branched.dat', 2,
        extra='\n&amr_tree_checkpoint\n checkpoint_interval_steps=1\n/\n',
        fail='does not yet support checkpoint, adaptation or subcycling groups')
    # Alias a topology to the output name without changing the input template.
    folder = (args.workdir / 'topology_alias').resolve()
    folder.mkdir(parents=True, exist_ok=True)
    input_path = folder / 'input.nml'
    input_path.write_text(template)
    topology = folder / 'tree_hotspot.csv'
    original = (args.case_dir / 'root.dat').read_bytes()
    topology.write_bytes(original)
    result = subprocess.run([args.mpiexec, args.numproc_flag, '2', str(args.mpi.resolve()),
                             str(input_path), str(topology)], cwd=folder, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=30, check=False)
    (folder / 'run.log').write_text(result.stdout)
    if result.returncode == 0 or 'MPI tree output aliases an input' not in result.stdout:
        raise RuntimeError('topology alias was not rejected')
    if topology.read_bytes() != original:
        raise RuntimeError('aliased topology was overwritten')
    print('PASS: static MPI tree application; serial CSV/budget bytes match at 2/3 ranks; empty ranks and input rejection')


if __name__ == '__main__':
    main()
