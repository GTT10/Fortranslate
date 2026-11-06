#!/usr/bin/env python3
"""Inject a controlled ENOSPC error in the first tree checkpoint write.

A test-only LD_PRELOAD shim targets tree.chk or tree.chk.tmp after 8192 bytes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import os
import subprocess



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--shim', type=Path, required=True)
    parser.add_argument('--expect-preserved', choices=('yes', 'no'), required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    assert not list(args.work.glob('*.tmp*')), 'Work directory contains a prior temporary file'
    binary = args.executable.resolve()
    target = args.work / 'tree.chk'
    expected = b'preserve\n'
    target.write_bytes(expected)
    command = [str(binary)]
    shim = args.shim.resolve()
    result = subprocess.run(command, cwd=args.work, capture_output=True, timeout=60,
                            env=dict(HOME=os.environ.get('HOME', '/tmp'), PATH='/usr/bin:/bin', LANG='C.UTF-8',
                                     LD_PRELOAD=str(shim)))
    (args.work / 'stdout.log').write_bytes(result.stdout)
    (args.work / 'stderr.log').write_bytes(result.stderr)
    actual = target.read_bytes() if target.exists() else b''
    rejected = b'write snapshot' in result.stdout + result.stderr
    read_rejected = b'read snapshot' in result.stdout + result.stderr
    observed = dict(command=command, returncode=result.returncode,
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    shim_sha256=hashlib.sha256(shim.read_bytes()).hexdigest(),
                    injection_observed=b'FAULT_SHIM serial checkpoint write ENOSPC' in result.stderr,
                    write_rejection_observed=rejected, subsequent_read_rejected=read_rejected,
                    target_exists=target.exists(), target_size=len(actual),
                    target_sha256=hashlib.sha256(actual).hexdigest(),
                    expected_original_sha256=hashlib.sha256(expected).hexdigest(),
                    original_preserved=actual == expected,
                    temporary_paths=[p.name for p in args.work.glob('*.tmp*')])
    (args.work / 'result.json').write_text(json.dumps(observed, indent=2) + '\n')
    print(json.dumps(observed, indent=2))
    assert observed['injection_observed'] and (rejected or read_rejected) and result.returncode != 0, 'Fixture did not reject its checkpoint write'
    assert observed['original_preserved'] == (args.expect_preserved == 'yes'), observed
    if args.expect_preserved == 'yes':
        assert rejected, 'Writer did not report failure'
        assert not observed['temporary_paths'], 'Failed attempt left temporary file(s)'


if __name__ == '__main__':
    main()
