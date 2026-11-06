#!/usr/bin/env python3
"""Compiled bounded-slope regression for a nearly flat outer stencil pair."""
import argparse
import json
import math
from pathlib import Path
import subprocess
from check_eb_extremum_slope import check


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir
    report = check(work, bounded_quartic=True)
    # Saved second-step donor line, with the nearly equal outer pair restored
    # exactly equal before controlled representable perturbations.
    base = [.7352585883390887, 2.0938987768722894, 3.1337620549887606,
            3.696533849157872, 3.696533849157872]
    cases = []
    for scale in (1e-8, 1., 1e8):
        for shift in (-10., 0., 10.):
            for sign in (-1., 1.):
                for ulps in range(-4, 5):
                    q = [scale*(sign*v+shift) for v in base]
                    for _ in range(abs(ulps)):
                        q[4] = math.nextafter(q[4], math.inf if ulps > 0 else -math.inf)
                    high = (22*(q[3]-q[1])-3*(q[4]-q[0]))/32
                    cases.append((q, high, 'roundoff'))
    # A resolved opposite outer slope must still take the extremum bound.
    for scale in (1e-8, 1., 1e8):
        q = [scale*v for v in base]
        q[4] -= scale*1e-6
        dm, dp = q[2]-q[1], q[3]-q[2]
        d2m = dm-(q[1]-q[0]); d2c = dp-dm; d2p = (q[4]-q[3])-dp
        expected = min(1.875*min(abs(d2c), -d2m, -d2p), 2*abs(dm))
        cases.append((q, expected, 'resolved'))
    stdin = ''.join(' '.join(format(v, '.17g') for v in q)+'\n' for q, _, _ in cases)
    output = subprocess.check_output([str(work/'probe')], input=stdin, text=True)
    rows = [list(map(float, line.split())) for line in output.splitlines()]
    assert len(rows) == len(cases)
    for (q, expected, kind), (slope, low, high) in zip(cases, rows):
        tol = 64*math.ulp(1.)*max(map(abs, q))
        assert abs(slope-expected) <= tol, (kind, q, slope, expected)
        assert abs((low+high)/2-q[2]) <= tol
    (work/'roundoff-inputs.txt').write_text(stdin)
    (work/'roundoff-outputs.txt').write_text(output)
    report.update(roundoff_cases=len(cases),
                  limitation='Scale and finite-offset regression, not arbitrary-offset invariance or CFD qualification')
    (work/'roundoff-report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(f'PASS: {len(cases)} roundoff/resolved controls plus {report["cases"]} scalar controls')


if __name__ == '__main__':
    main()
