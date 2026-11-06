#!/usr/bin/env python3
"""Check spatial translation invariance of the bounded 4^3 evaporation fixture.

This is a metamorphic stencil regression, not a fixed-PeleC or physical gate.
The momentum tolerance is fixed at 1e-12 kg m^-2 s^-1. Outputs are retained.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--case-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mpiexec', default='/usr/bin/mpiexec.openmpi')
    parser.add_argument('--ranks', type=int, choices=(1, 2), default=1)
    args = parser.parse_args()
    exe, fixture, output = (p.resolve() for p in (args.executable, args.case_dir, args.output))
    tree = (fixture / 'tree.nml').read_text()
    spray = (fixture / 'spray.nml').read_text()
    for token in ('root_shape=4,4,4', 'child_count=0', "boundary_faces=6*'periodic'",
                  'time_step=1e-7', 'final_time=2e-7', 'maximum_steps=2',
                  'domain_upper=.004,.004,.004'):
        if token not in tree:
            raise ValueError(f'Expected bounded fixture token: {token}')
    center = 'center=.0015,.0015,.0015'
    if spray.count(center) != 1:
        raise ValueError('Expected exactly one centered four-cell-fixture source')
    output.mkdir(parents=True, exist_ok=False)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    record = {'qualification': False, 'status': 'running', 'threshold_kg_m2_s': 1e-12,
              'binary': str(exe), 'binary_sha256': sha(exe), 'ranks': args.ranks,
              'input_sha256': {name: sha(fixture / name) for name in
                               ('tree.nml', 'spray.nml', 'level_set.dat')}, 'commands': []}

    def save():
        (output / 'result.json').write_text(json.dumps(record, indent=2) + '\n')

    save()
    env = {k: os.environ[k] for k in ('HOME', 'USER', 'LANG', 'TMPDIR') if k in os.environ}
    env.update(PATH='/usr/bin:/bin', OMP_NUM_THREADS='1', HWLOC_COMPONENTS='-gl,-x86')
    states = []
    try:
        for shift in (0, 1):
            work = output / f'shift-{shift}'
            work.mkdir()
            for name in ('tree.nml', 'spray.nml', 'level_set.dat'):
                shutil.copyfile(fixture / name, work / name)
            if shift:
                (work / 'spray.nml').write_text(spray.replace(center, 'center=.0025,.0015,.0015'))
            command = [args.mpiexec, '-n', str(args.ranks), str(exe), 'tree.nml', 'spray.nml']
            with (work / 'run.log').open('w') as log:
                result = subprocess.run(command, cwd=work, env=env, stdout=log,
                                        stderr=subprocess.STDOUT, timeout=120)
            record['commands'].append({'command': command, 'cwd': str(work),
                                       'returncode': result.returncode})
            save()
            if result.returncode or 'Completed steps: 2' not in (work / 'run.log').read_text():
                raise RuntimeError(f'Incomplete shifted case: {shift}')
            with (work / 'fortran.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            state = {((int(row['i']) - 1 - shift) % 4 + 1, int(row['j']), int(row['k'])):
                     row for row in rows}
            if len(rows) != 64 or len(state) != 64:
                raise ValueError('Expected 64 unique cell rows')
            if not all(math.isfinite(float(row[field])) for row in rows
                       for field in ('rhou', 'rhov', 'rhow')):
                raise ValueError('Nonfinite momentum')
            states.append(state)
        if states[0].keys() != states[1].keys():
            raise ValueError('Mismatched shifted cell keys')
        errors = {field: max(abs(float(states[0][key][field]) - float(states[1][key][field]))
                             for key in states[0]) for field in ('rhou', 'rhov', 'rhow')}
        # NaN comparisons cannot pass this all-components check.
        passed = all(value <= record['threshold_kg_m2_s'] for value in errors.values())
        if sha(exe) != record['binary_sha256']:
            raise RuntimeError('Executable changed during check')
        record.update(status='passed' if passed else 'failed', maximum_shift_difference=errors)
        record['output_sha256'] = {str(p.relative_to(output)): sha(p)
                                   for p in output.glob('shift-*/*.csv')}
        save()
        print(json.dumps(record, indent=2))
        return 0 if passed else 1
    except Exception as error:
        record.update(status='error', error=repr(error))
        save()
        raise


if __name__ == '__main__':
    raise SystemExit(main())
