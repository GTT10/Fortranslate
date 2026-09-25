#!/usr/bin/env python3
"""Compare normal PPM reconstruction/projection against pinned PeleC C++ bodies.

Uses a NASA7 SI adapter, no AMReX execution, WENO, flattening calculation or CTU.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import random
from check_pinned_ppm_kernels import function, run
from check_pinned_riemann_kernel import PREFIX

SHIM = r'''
constexpr int QRHO=0,QUN=1,QUT=2,QUTT=3,QPRES=4,QFS=5,QREINT=12,QVAR=13;
struct Face { double q[13]{}; double& operator()(int,int n){return q[n];} };
void project(double Im[13][3],double Ip[13][3],double un,double cc,Face& qp,Face& qm){
  EOS eos;
  int idir=0,i=0,j=0,k=0,iv=0,ivp1=1,vlo[3]={0,0,0},vhi[3]={0,0,0};
'''
DRIVER = r'''
int main(){
 for(auto& row:thermo)for(double& v:row)if(!(std::cin>>v))return 2;
 double q[7][12],flat,ratio;
 std::cout<<std::setprecision(17);
 while(std::cin>>q[0][0]){
  for(int j=0;j<7;j++)for(int n=0;n<12;n++)if(j||n)std::cin>>q[j][n];
  if(!(std::cin>>flat>>ratio))return 3;
  EOS eos; double cc; eos.RPY2Cs(q[3][0],q[3][4],q[3]+5,cc);
  // Traced internal energy does not enter the published rho/u/p projection.
  // The pinned body recovers final internal energy from rho/Y/p through EOS.
  double Ip[13][3]{},Im[13][3]{};
  for(int n=0;n<12;n++){
   double s[5],sm,sp;
   for(int j=0;j<5;j++)s[j]=q[j+1][n];
   ppm_reconstruct(s,flat,sm,sp);
   ppm_int_profile(sm,sp,s[2],q[3][1],cc,ratio,Ip[n],Im[n]);
  }
  Face qp,qm; project(Im,Ip,q[3][1],cc,qp,qm);
  for(int n=0;n<12;n++)std::cout<<qp.q[n]<<' ';
  for(int n=0;n<12;n++)std::cout<<qm.q[n]<<' ';
  std::cout<<'\n';
 }
}
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--cxx',default='c++')
    a=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    pin=json.loads((root/'references/pelec_baseline.json').read_text())['pelec']['commit']
    sources={name:run(['git','-C',str(a.reference),'show',f'{pin}:Source/{name}']) for name in ('PPM.H','PPM.cpp')}
    work=a.work_dir.resolve();work.mkdir(parents=True,exist_ok=True)
    for name,source in sources.items(): (work/name).write_text(source)
    # Keep the complete species and left/right projection block unchanged.
    source=sources['PPM.cpp']
    projection=source[source.index('    for (int n = QFS;'):source.rindex('  });')]
    cpp=work/'reference.cpp'
    cpp.write_text(PREFIX+'constexpr int im2=0,im1=1,i0=2,ip1=3,ip2=4;\n'+function(sources['PPM.H'],'ppm_reconstruct')+'\n'+
                   function(sources['PPM.H'],'ppm_int_profile')+SHIM+projection+'\n}\n'+DRIVER)
    command=[a.cxx,'-std=c++17','-O0','-ffp-contract=off',str(cpp),'-o',str(work/'reference')]
    run(command)
    exe=str(a.executable.resolve());thermo=run([exe,'--thermo-stream'])
    (work/'thermo.txt').write_text(thermo)
    mw=[float(line.split()[0]) for line in thermo.splitlines()]
    rng=random.Random(20260926);rows=[]
    for velocity in (-2500.,-100.,-1e-8,0.,1e-8,100.,2500.):
        for flat in (0.,.37,1.):
            for shape in ('constant','linear','contact','random'):
                stencil=[]
                for j in range(-3,4):
                    factor=(j*.02 if shape in ('linear','contact') else rng.uniform(-.05,.05) if shape=='random' else 0.)
                    y=[.01+factor*.01,.001,.001,.2,.001,.05,0.]
                    y[-1]=1-sum(y[:-1]);p=1e5*(1+factor if shape!='contact' else 1.)
                    t=900*(1+factor)
                    r=sum(v*8314.46261815324/m for v,m in zip(y,mw))
                    u=velocity+(20*factor if shape!='contact' else 0.)
                    stencil.extend([p/(r*t),u,3+factor,-2+factor,p]+y)
                for ratio in (0.,1e-5,8e-5):rows.append(stencil+[flat,ratio])
    inputs=''.join(' '.join(format(v,'.17g') for v in row)+'\n' for row in rows)
    (work/'inputs.txt').write_text(inputs)
    ref=run([str(work/'reference')],input=thermo+inputs)
    actual=run([exe,'--reference-stream'],input=inputs)
    (work/'reference.txt').write_text(ref);(work/'actual.txt').write_text(actual)
    expected=[[float(v) for v in line.split()] for line in ref.splitlines()]
    observed=[[float(v) for v in line.split()] for line in actual.splitlines()]
    assert len(expected)==len(observed)==len(rows)
    worst=0.
    for i,(x,y) in enumerate(zip(expected,observed)):
        assert len(x)==len(y)==24 and all(math.isfinite(v) for v in x+y)
        error=max(abs(v-w)/max(1.,abs(v),abs(w)) for v,w in zip(x,y))
        assert error<=3e-12,(i,error,x,y)
        worst=max(worst,error)
    report=dict(reference_commit=pin,source_sha256={k:hashlib.sha256(v.encode()).hexdigest() for k,v in sources.items()},
                executable_sha256=hashlib.sha256(a.executable.read_bytes()).hexdigest(),compile_command=command,
                compiler=run([a.cxx,'--version']).splitlines()[0],
                scope='normal primitive states; internal energy, WENO and CTU excluded',
                cases=len(rows),values=24*len(rows),max_scaled_error=worst,tolerance=3e-12,status='PASS')
    (work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
