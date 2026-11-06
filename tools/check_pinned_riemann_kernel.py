#!/usr/bin/env python3
"""Compare pinned PeleC Riemann arithmetic with a NASA7 EOS adapter in SI units."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import random
from check_pinned_ppm_kernels import function, run


PREFIX = r'''
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
namespace amrex { using Real=double; using std::min; using std::max; }
constexpr int NUM_SPECIES=7;
double thermo[7][16];
namespace constants {
 double smallu(){return 1e-12;} double small_num(){return 1e-8;}
 double very_small_num(){return std::numeric_limits<double>::epsilon()*1e-100;}
}
struct EOS {
 void properties(double rho,const double* y,double p,double& e,double& c) {
  double r=0,cp=0,h=0;
  for(int n=0;n<7;n++) r+=y[n]*8314.46261815324/thermo[n][0];
  double t=p/(rho*r);
  for(int n=0;n<7;n++) {
   const double* a=&thermo[n][t<=thermo[n][1]?2:9];
   double rn=8314.46261815324/thermo[n][0];
   cp+=y[n]*rn*(a[0]+a[1]*t+a[2]*t*t+a[3]*t*t*t+a[4]*t*t*t*t);
   h+=y[n]*rn*t*(a[0]+a[1]*t/2+a[2]*t*t/3+a[3]*t*t*t/4+a[4]*t*t*t*t/5+a[5]/t);
  }
  e=h-r*t; c=std::sqrt(cp/(cp-r)*p/rho);
 }
 void RPY2Cs(double rho,double p,const double* y,double& c){double e;properties(rho,y,p,e,c);}
 void RYP2E(double rho,const double* y,double p,double& e){double c;properties(rho,y,p,e,c);}
};
namespace pele { namespace physics { struct PhysicsType { static EOS eos(){return EOS{};} }; } }
'''
DRIVER = r'''
int main(){
 for(auto& row:thermo) for(double& v:row) if(!(std::cin>>v))return 2;
 double l[12],r[12],cav,f[12],us,ei,iu,iv,iw,p,gamma;
 std::cout<<std::setprecision(17);
 while(std::cin>>l[0]){
  for(int n=1;n<12;n++)std::cin>>l[n];
  for(double& v:r)std::cin>>v;
  if(!(std::cin>>cav))return 3;
  riemann(l[0],l[1],l[2],l[3],l[4],l+5,r[0],r[1],r[2],r[3],r[4],r+5,1,cav,
    us,f[0],f+5,f[1],f[2],f[3],f[4],ei,iu,iv,iw,p,gamma);
  for(double v:f)std::cout<<v<<' ';
  std::cout<<iu<<' '<<p<<'\n';
 }
}
'''


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--cxx', default='c++')
    parser.add_argument('--face-flux', action='store_true', help='include pc_cmpflx passive upwinding')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    pin = json.loads((root/'references/pelec_baseline.json').read_text())['pelec']['commit']
    header = run(['git','-C',str(args.reference),'show',f'{pin}:Source/Riemann.H'])
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    (work/'Riemann.H').write_text(header)
    cpp = work/'reference.cpp'
    driver = DRIVER
    passive = ''
    if args.face_flux:
        godunov = run(['git','-C',str(args.reference),'show',f'{pin}:Source/Godunov.H'])
        (work/'Godunov.H').write_text(godunov)
        passive = function(godunov,'pc_cmpflx_passive')
        driver = driver.replace('  for(double v:f)',
            '  for(int n=0;n<7;n++)pc_cmpflx_passive(us,f[0],l[5+n],r[5+n],f[5+n]);\n  for(double v:f)')
    cpp.write_text(PREFIX + function(header,'riemann') + passive + driver)
    command = [args.cxx,'-std=c++17','-O0','-ffp-contract=off',str(cpp),'-o',str(work/'reference')]
    run(command)
    exe = str(args.executable.resolve())
    thermo = run([exe,'--thermo-stream'])
    (work/'thermo.txt').write_text(thermo)
    mw = [float(row.split()[0]) for row in thermo.splitlines()]
    rng = random.Random(20260925)
    def state(u: float, randomize: bool) -> list[float]:
        y = [.01, .001, .001, .2, .001, .05, .737]
        if randomize:
            y[0] += rng.uniform(0,.02)
            y[-1] = 1-sum(y[:-1])
        p = 1e5*(rng.uniform(.8,1.2) if randomize else 1)
        t = rng.uniform(600,1800) if randomize else 900
        r = sum(v*8314.46261815324/m for v,m in zip(y,mw))
        return [p/(r*t),u,3.,-2.,p]+y
    rows=[]
    for u in (-2500.,-100.,-1e-8,0.,1e-8,100.,2500.):
        q=state(u,False)
        rows.append(q+q+[700.])
        # Material contact: equal pressure/velocity with different composition and density.
        contact=state(u,True)
        contact[0]*=q[4]/contact[4]
        contact[4]=q[4]
        rows.append(q+contact+[700.])
        for _ in range(6):
            rows.append(state(u+rng.uniform(-50,50),True)+state(u+rng.uniform(-50,50),True)+[700.])
    # Signed CTU passive weights are consumed by Fuego without normalization.
    # Exercise both upwind directions and stationary/moving interfaces.
    for weight in (-1e-24, -.001):
        for velocity in (-2500., -100., 0., 2500.):
            q = state(velocity, False)
            q[5] = weight
            q[-1] = 1-sum(q[5:-1])
            gas_constant = sum(v*8314.46261815324/m for v,m in zip(q[5:],mw))
            q[0] = q[4]/(gas_constant*900)
            right = state(velocity, False)
            rows.append(q+right+[700.])
    if args.face_flux:
        for weight in (.8, 1.2):
            for velocity in (-2500., -100., 0., 100., 2500.):
                left=state(velocity,False);right=state(velocity,False)
                left[5:]=[v*weight for v in left[5:]]
                right[5:]=[v*(2-weight) for v in right[5:]]
                rows.append(left+right+[700.])
    # Both cells are within the NASA fit interval; the acoustic star state
    # expands below it. Fuego evaluates that positive transient temperature.
    y = [0., 0., 0., .3, 0., 0., .7]
    r = sum(v*8314.46261815324/m for v,m in zip(y,mw))
    pressure_expansion_cases = 0
    for ratio in (10., 100., 1000.):
        left = [1e6/(r*300.), 0., 0., 0., 1e6]+y
        right = [left[0]/ratio, 0., 0., 0., left[4]/ratio]+y
        rows.extend((left+right+[700.], right+left+[700.]))
        pressure_expansion_cases += 2
    inputs=''.join(' '.join(format(v,'.17g') for v in row)+'\n' for row in rows)
    (work/'inputs.txt').write_text(inputs)
    ref=run([str(work/'reference')],input=thermo+inputs)
    actual=run([exe,'--face-stream' if args.face_flux else '--reference-stream'],input=inputs)
    (work/'reference.txt').write_text(ref)
    (work/'actual.txt').write_text(actual)
    a=[[float(v) for v in row.split()] for row in ref.splitlines()]
    b=[[float(v) for v in row.split()] for row in actual.splitlines()]
    assert len(a)==len(b)==len(rows)
    worst=0.
    for i,(expected,observed) in enumerate(zip(a,b)):
        assert len(expected)==len(observed)==14
        assert all(math.isfinite(v) for v in expected+observed)
        err=max(abs(x-y)/max(1.,abs(x),abs(y)) for x,y in zip(expected,observed))
        assert err<=3e-12,(i,err,expected,observed)
        worst=max(worst,err)
    report=dict(reference_commit=pin,header_sha256=hashlib.sha256(header.encode()).hexdigest(),
                godunov_sha256=hashlib.sha256(godunov.encode()).hexdigest() if args.face_flux else None,
                executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
                face_flux=args.face_flux,cases=len(rows),pressure_expansion_cases=pressure_expansion_cases,
                values=14*len(rows),max_scaled_error=worst,tolerance=3e-12,
                compile_command=command,eos='NASA7 SI adapter using exported PeleF coefficients',status='PASS')
    (work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
