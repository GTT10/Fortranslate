#!/usr/bin/env python3
"""Isolated scalar EP slope: analytic averages and discontinuity bounds.

Source: Sekora & Colella arXiv:0903.4200v2 section 2.1, C_VL=1.25.
This is not an EB family, EOS, evolution, or fixed-PeleC parity test.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess


def check(work, bounded_quartic=False, continuous_quartic=False):
    root = Path(__file__).resolve().parents[1]
    work.mkdir(parents=True, exist_ok=True)
    sources = [root/'tools/reference_probes'/name for name in
               ('eb_extremum_slope_mod.F90',
                'eb_continuous_quartic_probe.F90' if continuous_quartic else
                'eb_bounded_quartic_probe.F90' if bounded_quartic else 'eb_extremum_slope_probe.F90')]
    binary = work/'probe'
    command = ['gfortran', '-std=f2018', '-O0', '-Wall', '-Wextra',
               '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
               '-J'+str(work), '-I'+str(work), *map(str, sources), '-o', str(binary)]
    subprocess.run(command, check=True, capture_output=True, text=True)
    cases = []
    def add(name, q, exact=None, bounded=False):
        cases.append(dict(name=name, q=q, exact=exact, bounded=bounded))
    for scale in (1e-8, 1., 1e8):
        for shift in (-.5, -.25, 0., .25, .5):
            # Exact cell and half-cell averages of scale*(x-shift)^2.
            q = [scale*((j-shift)**2+1/12) for j in range(-2, 3)]
            exact = [scale*((x-shift)**2+1/48) for x in (-.25, .25)]
            add(f'quadratic_{scale}_{shift}', q, exact)
        add(f'constant_{scale}', [scale]*5, [scale]*2)
        add(f'linear_{scale}', [scale*(j+3) for j in range(-2, 3)],
            [scale*2.75, scale*3.25])
    for offset in (k/8 for k in range(-20, 21)):
        q = [max(0., min(1., j+.5-offset)) for j in range(-2, 3)]
        add(f'step_{offset}', q, bounded=True)
    def logcosh(x):
        return abs(x)+math.log1p(math.exp(-2*abs(x)))-math.log(2)
    for width in (.125, .25, .5, 1., 2., 4.):
        for offset in (k/8 for k in range(-8, 9)):
            q = [width*(logcosh((j+.5-offset)/width)-
                        logcosh((j-.5-offset)/width)) for j in range(-2, 3)]
            add(f'tanh_{width}_{offset}', q, bounded=True)
    for n in (8, 16, 32, 64):
        h = 2*math.pi/n
        for j in range(1, n+1):
            q = [math.sin((j+k-.5)*h)*math.sin(h/2)/(h/2) for k in range(-2, 3)]
            exact = [math.sin((j-.5+off)*h)*math.sin(h/4)/(h/4) for off in (-.25, .25)]
            add(f'sine_{n}_{j}', q, exact)
    stdin = ''.join(' '.join(format(v, '.17g') for v in c['q'])+'\n' for c in cases)
    run = subprocess.run([str(binary)], input=stdin, text=True, capture_output=True, check=True)
    (work/'inputs.txt').write_text(stdin)
    (work/'outputs.txt').write_text(run.stdout)
    rows = [list(map(float, line.split())) for line in run.stdout.splitlines()]
    assert len(rows) == len(cases)
    sine = {}
    conventional = {}
    for c, row in zip(cases, rows):
        slope, low, high = row
        assert all(map(math.isfinite, row))
        tol = 2e-14*max(1., *map(abs, c['q']))
        assert abs((low+high)/2-c['q'][2]) <= tol, c
        if c['bounded']:
            assert min(c['q'][1:4])-tol <= min(low, high), c
            assert max(low, high) <= max(c['q'][1:4])+tol, c
        if c['exact']:
            err = max(abs(x-y) for x,y in zip((low,high),c['exact']))
            if c['name'].startswith('sine_'):
                n = c['name'].split('_')[1]
                sine[n] = max(sine.get(n,0.), err)
                q = c['q']
                dm, dp = q[2]-q[1], q[3]-q[2]
                dc = (q[3]-q[1])/2
                mc = math.copysign(min(abs(dc),2*abs(dm),2*abs(dp)),dc) if dm*dp>0 else 0.
                mcerr = max(abs(x-y) for x,y in zip((q[2]-mc/4,q[2]+mc/4),c['exact']))
                conventional[n] = max(conventional.get(n,0.),mcerr)
            else:
                assert err <= tol, (c, row, err)
    result = dict(status='passed', qualification=False, cases=len(cases),
                  variant='experimental_continuous_quartic' if continuous_quartic else
                  'experimental_bounded_quartic' if bounded_quartic else 'published_EP_slope',
                  command=command, source='https://arxiv.org/html/0903.4200v2',
                  compiler=subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
                  source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                  binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  sine_max_child_error=sine,
                  conventional_sine_max_child_error=conventional)
    (work/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--work-dir', type=Path, required=True)
    ap.add_argument('--bounded-quartic', action='store_true')
    ap.add_argument('--continuous-quartic', action='store_true')
    args = ap.parse_args()
    print(json.dumps(check(args.work_dir, args.bounded_quartic, args.continuous_quartic),indent=2))
