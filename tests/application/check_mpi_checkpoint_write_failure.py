#!/usr/bin/env python3
"""Inject a late MPI write error into the existing tiny checkpoint unit fixture.

The unit deliberately fails its ordinary successful-write assertion after this
injection. The expected result of this runner is not a successful unit exit:
it is rejection of that write and preservation of the prior public bytes.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--shim', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--expect-preserved', choices=('yes', 'no'), required=True)
    parser.add_argument('--mpiexec', default='/usr/bin/mpiexec')
    parser.add_argument('--mode', choices=('write', 'close'), default='write')
    parser.add_argument('--rank', type=int, choices=(0, 1), default=0)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    assert not list(args.work.glob('*.tmp*')), 'Work directory contains a prior temporary file'
    binary, shim = args.executable.resolve(), args.shim.resolve()
    env = dict(HOME=os.environ.get('HOME', '/tmp'), PATH='/usr/bin:/bin', LANG='C.UTF-8',
               LD_PRELOAD=str(shim),
               PELEF_TEST_FAIL_MPI_WRITE_CALL='3' if args.rank == 0 else '1',
               PELEF_TEST_FAIL_MPI_MODE=args.mode, PELEF_TEST_FAIL_MPI_RANK=str(args.rank),
               OMPI_MCA_rmaps_base_oversubscribe='1')
    command = [args.mpiexec, '-n', '2', str(binary)]
    result = subprocess.run(command, cwd=args.work, env=env, capture_output=True, timeout=60)
    (args.work / 'stdout.log').write_bytes(result.stdout)
    (args.work / 'stderr.log').write_bytes(result.stderr)
    target = args.work / 'mpi_tree_checkpoint_preflight_2.chk'
    actual = target.read_bytes() if target.exists() else b''
    expected = b'preserve\n'
    injected = f'FAULT_SHIM mode={args.mode} rank={args.rank} '.encode() in result.stderr
    rejected = b'owner-local checkpoint write' in result.stdout + result.stderr
    observed = dict(command=command, returncode=result.returncode,
                    binary_sha256=sha(binary), shim_sha256=sha(shim),
                    injection_observed=injected, write_rejection_observed=rejected, fault_mode=args.mode, fault_rank=args.rank, target_exists=target.exists(),
                    target_size=len(actual), target_sha256=hashlib.sha256(actual).hexdigest(),
                    expected_original_sha256=hashlib.sha256(expected).hexdigest(),
                    original_preserved=actual == expected,
                    public_changed_on_rejection=actual != expected and bool(actual),
                    temporary_paths=[p.name for p in args.work.glob('*.tmp*')])
    (args.work / 'result.json').write_text(json.dumps(observed, indent=2) + '\n')
    print(json.dumps(observed, indent=2))
    assert injected and rejected and result.returncode != 0, 'Fault was not observed or unit did not reject it'
    assert observed['original_preserved'] == (args.expect_preserved == 'yes'), observed
    if args.expect_preserved == 'yes':
        assert not observed['temporary_paths'], 'Failed attempt left temporary file(s)'


if __name__ == '__main__':
    main()
