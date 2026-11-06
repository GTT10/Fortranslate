#!/usr/bin/env python3
"""Scalar continuity stress test; not tensor, EOS or flow qualification."""
import argparse
import itertools
import json
import math
from pathlib import Path
import subprocess
from check_eb_extremum_slope import check


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--work-dir', type=Path, required=True)
    args = ap.parse_args()
    report = check(args.work_dir, continuous_quartic=True)
    cases, pairs = [], []
    for q in itertools.product((-1., 0., 1.), repeat=5):
        for axis in range(5):
            lo, hi = list(q), list(q)
            lo[axis] -= 1e-10
            hi[axis] += 1e-10
            pairs.append((len(cases), len(cases)+1, 2e-10))
            cases.extend((lo, hi))
    base = [.7352585883390887, 2.0938987768722894, 3.1337620549887606,
            3.696533849157872, 3.696533849157872]
    for reverse in (False, True):
        previous = None
        for ulps in range(-32, 33):
            q = base.copy()
            for _ in range(abs(ulps)):
                q[4] = math.nextafter(q[4], math.inf if ulps > 0 else -math.inf)
            if reverse:
                q.reverse()
            if previous is not None:
                pairs.append((previous, len(cases), max(abs(a-b) for a,b in zip(cases[previous], q))))
            previous = len(cases)
            cases.append(q)
    stdin = ''.join(' '.join(format(x, '.17g') for x in q)+'\n' for q in cases)
    out = subprocess.check_output([str(args.work_dir/'probe')], input=stdin, text=True)
    rows = [list(map(float, line.split())) for line in out.splitlines()]
    assert len(rows) == len(cases)
    max_gain = 0.
    for a, b, delta in pairs:
        gain = abs(rows[a][0]-rows[b][0])/delta
        max_gain = max(max_gain, gain)
        assert gain < 100, (cases[a], cases[b], rows[a], rows[b], gain)
    report.update(perturbation_cases=len(cases), perturbation_pairs=len(pairs),
                  max_observed_slope_gain=max_gain,
                  limitation='Finite perturbation sampling is not a global Lipschitz proof')
    (args.work_dir/'continuity-inputs.txt').write_text(stdin)
    (args.work_dir/'continuity-outputs.txt').write_text(out)
    (args.work_dir/'continuity-report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(f'PASS: {len(cases)} perturbation samples, {len(pairs)} pairs, max gain {max_gain:g}')


if __name__ == '__main__':
    main()
