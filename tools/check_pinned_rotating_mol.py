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
    std::cout<<ei<<'\n';
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
    # A pinned source archive need not contain .git; never resolve its parent checkout.
    header_bytes = (args.reference/'Source/Riemann.H').read_bytes()
    assert hashlib.sha256(header_bytes).hexdigest() == '073ae14b2f7925de2858ef312d02fc50f033c2a344dccc5b24ecde1f77ed8253'
    header = header_bytes.decode()
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
    # Positive transient face temperatures may leave the accepted-cell fit
    # range after MOL reconstruction. Check both reference flux branches.
    y = [0., 0., 0., .3, 0., 0., .7]
    gas_constant = sum(v*8314.46261815324/m for v, m in zip(y, mw))
    transient_face_cases = 0
    for velocity in (-100., 0., 100.):
        for potentials in ((0., 0.), (400., 1600.)):
            left = [.6, velocity, 0., 0., .6*gas_constant*180.]+y
            right = [1.4, velocity, 0., 0., 1.4*gas_constant*301.]+y
            for lax in (0, 1):
                rows.append(left+right+[700., 700., *potentials, lax])
                transient_face_cases += 1
    # At7000K O2's extrapolated Cv is negative, while this mixture's Cv
    # remains positive. Fuego qualifies the mixture sound speed.
    signed_species_cv_cases = 0
    for velocity in (-100., 0., 100.):
        q = [1., velocity, 3., -2., gas_constant*7000.]+y
        for potentials in ((0., 0.), (400., 1600.)):
            for lax in (0, 1):
                rows.append(q+q+[700., 700., *potentials, lax])
                signed_species_cv_cases += 1
    inputs = ''.join(' '.join(format(v, '.17g') for v in row)+'\n' for row in rows)
    (work/'inputs.txt').write_text(inputs)
    ref = run([str(work/'reference')], input=thermo+inputs)
    actual = run([exe, '--internal-energy-stream'], input=inputs)
    (work/'reference.txt').write_text(ref)
    (work/'actual.txt').write_text(actual)
    expected = [[float(v) for v in row.split()] for row in ref.splitlines()]
    observed = [[float(v) for v in row.split()] for row in actual.splitlines()]
    assert len(expected) == len(observed) == len(rows)
    worst = 0.
    for i, (a, b) in enumerate(zip(expected, observed)):
        assert len(a) == len(b) == 13
        assert all(math.isfinite(v) for v in a+b)
        error = max(abs(x-y)/max(1., abs(x), abs(y)) for x, y in zip(a, b))
        assert error <= 3e-12, (i, error, a, b)
        worst = max(worst, error)
    report = dict(reference_commit=pin, header_sha256=hashlib.sha256(header.encode()).hexdigest(),
                  executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
                  cases=len(rows), transient_face_cases=transient_face_cases,
                  signed_species_cv_cases=signed_species_cv_cases,
                  values=13*len(rows), max_scaled_error=worst, tolerance=3e-12,
                  compile_command=command, eos='NASA7 SI adapter using exported PeleF coefficients',
                  components='mass, momentum, total energy, species, auxiliary internal-energy flux',
                  status='PASS')
    (work/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
