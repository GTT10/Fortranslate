#!/usr/bin/env python3
"""Run frozen cut-cell pressure/velocity diagnostic controls and check budgets.

The imported input must be V4 with an EB plane and complete cell records.
No physical-validation or production held-MOL qualification is implied.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from audit_eb_initial_interface import sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe', type=Path, required=True)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise ValueError('Use a new output directory to preserve frozen runs')
    lines = a.input.read_text().splitlines()
    if lines[0] != 'PELEF_REACTING_EB_REFERENCE_V4' or lines[3].strip() != '1':
        raise ValueError('Expected converged V4 input')
    a.output.mkdir(parents=True)
    out = a.output.resolve()
    probe = out/'frozen-probe'
    shutil.copy2(a.probe, probe)
    inp = out/'input.dat'
    shutil.copy2(a.input, inp)
    scripts = Path(__file__).resolve().parent
    report = {'status': 'running', 'qualification': False, 'commands': [],
              'binary_sha256': sha(probe), 'input_sha256': sha(inp), 'script_sha256': sha(__file__)}

    def save():
        (out/'record.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    def run(cmd, log):
        start = time.monotonic()
        with log.open('w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, timeout=300)
        report['commands'].append({'argv': list(map(str, cmd)), 'returncode': result.returncode,
                                   'seconds': time.monotonic()-start})
        save()
        result.check_returncode()

    try:
        for profile, option in [('smooth_cut', '--smooth-cut-contact'), ('imported_cut', '--imported-state')]:
            for policy in ('pressure', 'velocity'):
                prefix = out/(profile+'_'+policy)
                cmd = [str(probe), str(inp), str(prefix), '--audit-initial', '--separate-momentum',
                       '--pressure-evolve-reset', option]
                if policy == 'velocity':
                    cmd.append('--consistent-velocity')
                run(cmd, Path(str(prefix)+'.log'))
                for script, suffix in [('check_eb_pressure_evolution.py', 'checks'),
                                       ('localize_eb_pressure_evolution.py', 'location'),
                                       ('check_eb_pressure_sync.py', 'sync')]:
                    cmd = [sys.executable, str(scripts/script), '--prefix', str(prefix),
                           '--output', str(prefix)+'-'+suffix+'.json']
                    if suffix == 'checks':
                        cmd += ['--input', str(inp)]
                    run(cmd, Path(str(prefix)+'-'+suffix+'.log'))
            # The extra candidate flag must leave the baseline mode unchanged.
            before = out/(profile+'_pressure_baseline_evolution.csv')
            after = out/(profile+'_velocity_baseline_evolution.csv')
            if before.read_bytes() != after.read_bytes():
                raise AssertionError('Candidate changed baseline history')
        with (out/'smooth_cut_velocity_paired_evolution.csv').open() as f:
            maximum = max(float(r['max_pressure_offset_Pa']) for r in csv.DictReader(f))
        report['smooth_cut_max_pressure_offset_Pa'] = maximum
        if maximum > 1e-8:
            raise AssertionError('Smooth cut contact equilibrium failed')
        if sha(probe) != report['binary_sha256'] or sha(inp) != report['input_sha256']:
            raise AssertionError('Frozen run identity changed')
        report['status'] = 'completed'
    except Exception as exc:
        report['status'] = 'failed'
        report['error'] = str(exc)
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
