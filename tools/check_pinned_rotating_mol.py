#!/usr/bin/env python3
"""Check rotating MOL energy fluxes against compiled pinned Riemann.H bodies."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import random

from check_pinned_ppm_kernels import function, run
from check_pinned_riemann_kernel import PREFIX

DRIVER = r'''
int main() {
  for(auto& row:thermo)for(double& v:row)if(!(std::cin>>v))return 2;
  double l[12],r[12],c[2],phi[2],f[12],us,eig,ei,iu,iv,iw,p,gamma;
  int lax;
  std::cout<<std::setprecision(17);
  while(std::cin>>l[0]) {
    for(int n=1;n<12;n++)std::cin>>l[n];
    for(double& v:r)std::cin>>v;
    if(!(std::cin>>c[0]>>c[1]>>phi[0]>>phi[1]>>lax))return 3;
    double rl=std::sqrt(2*phi[0]),rr=std::sqrt(2*phi[1]);
    if(lax) {
      laxfriedrich_flux(l[0],l[1],l[2],l[3],l[4],l+5,r[0],r[1],r[2],r[3],r[4],r+5,
        1,.5*(c[0]+c[1]),us,eig,f[0],f+5,f[1],f[2],f[3],f[4],ei,1,rl,1,rr);
    } else {
      riemann(l[0],l[1],l[2],l[3],l[4],l+5,r[0],r[1],r[2],r[3],r[4],r+5,
        1,.5*(c[0]+c[1]),us,f[0],f+5,f[1],f[2],f[3],f[4],ei,iu,iv,iw,p,gamma,1,rl,1,rr);
    }
    for(double v:f)std::cout<<v<<' ';
    std::cout<<'\n';
  }
}
'''


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--cxx', default='c++')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    pin = json.loads((root/'references/pelec_baseline.json').read_text())['pelec']['commit']
    header = run(['git', '-C', str(args.reference), 'show', f'{pin}:Source/Riemann.H'])
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    (work/'Riemann.H').write_text(header)
    cpp = work/'reference.cpp'
    cpp.write_text(PREFIX + function(header, 'riemann') + function(header, 'laxfriedrich_flux') + DRIVER)
    command = [args.cxx, '-std=c++17', '-O0', '-ffp-contract=off', str(cpp), '-o', str(work/'reference')]
    run(command)
    exe = str(args.executable.resolve())
    thermo = run([exe, '--thermo-stream'])
    (work/'thermo.txt').write_text(thermo)
    mw = [float(row.split()[0]) for row in thermo.splitlines()]
    assert len(mw) == 7, 'reference adapter requires the seven-species test mechanism'
    rng = random.Random(20260927)

    def state(velocity: float, vary: bool) -> list[float]:
        y = [.01, .001, .001, .2, .001, .05, .737]
        if vary:
            y[0] += rng.uniform(0, .02)
            y[-1] = 1-sum(y[:-1])
        pressure = 1e5*(rng.uniform(.8, 1.2) if vary else 1)
        temperature = rng.uniform(600, 1800) if vary else 900
        gas_constant = sum(v*8314.46261815324/m for v, m in zip(y, mw))
        return [pressure/(gas_constant*temperature), velocity, 3., -2., pressure]+y

    rows = []
    for velocity in (-2500., -100., 0., 100., 2500.):
        for potentials in ((0., 0.), (200., 700.), (4e5, 1e5)):
            for vary in (False, True):
                left = state(velocity, vary)
                right = state(velocity, vary)
                for lax in (0, 1):
                    rows.append(left+right+[700., 700., *potentials, lax])
    inputs = ''.join(' '.join(format(v, '.17g') for v in row)+'\n' for row in rows)
    (work/'inputs.txt').write_text(inputs)
    ref = run([str(work/'reference')], input=thermo+inputs)
    actual = run([exe, '--rotating-stream'], input=inputs)
    (work/'reference.txt').write_text(ref)
    (work/'actual.txt').write_text(actual)
    expected = [[float(v) for v in row.split()] for row in ref.splitlines()]
    observed = [[float(v) for v in row.split()] for row in actual.splitlines()]
    assert len(expected) == len(observed) == len(rows)
    worst = 0.
    for i, (a, b) in enumerate(zip(expected, observed)):
        assert len(a) == len(b) == 12
        assert all(math.isfinite(v) for v in a+b)
        error = max(abs(x-y)/max(1., abs(x), abs(y)) for x, y in zip(a, b))
        assert error <= 3e-12, (i, error, a, b)
        worst = max(worst, error)
    report = dict(reference_commit=pin, header_sha256=hashlib.sha256(header.encode()).hexdigest(),
                  executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
                  cases=len(rows), values=12*len(rows), max_scaled_error=worst, tolerance=3e-12,
                  compile_command=command, eos='NASA7 SI adapter using exported PeleF coefficients',
                  status='PASS')
    (work/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
