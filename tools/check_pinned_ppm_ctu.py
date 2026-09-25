#!/usr/bin/env python3
"""Compare unchanged pinned twelve-solve CTU assembly on small periodic grids.

The serial array shim executes reference loops; EOS uses exported NASA7 data.
WENO, forcing and physical boundaries are off. Optional shock/viscosity cases.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path
from check_pinned_ppm_kernels import run
from check_pinned_ppm_transverse import extract
from check_pinned_riemann_kernel import PREFIX

DRIVER=r'''
int main(){
 for(auto& row:thermo)for(double& v:row)if(!(std::cin>>v))return 2;
 double dt,difmag;int flattening;std::cout<<std::setprecision(17);
 while(std::cin>>dt>>flattening>>difmag){
  double cells[4][4][4][12];
  for(auto& plane:cells)for(auto& row:plane)for(auto& cell:row)for(double& v:cell)if(!(std::cin>>v))return 3;
  using namespace amrex;
  Box bx{{0,0,0},{3,3,3}},wide=grow(bx,4),fb{{0,0,0},{4,4,4}};
  FArrayBox prim(wide,QVAR),aux(wide,2),src(wide,QVAR),vol(bx,1),pd(bx,1),cons(wide,NVAR),divergence(fb,1);
  FArrayBox fx(fb,NVAR),fy(fb,NVAR),fz(fb,NVAR),qx(fb,NGDNV),qy(fb,NGDNV),qz(fb,NGDNV),area(fb,1);
  auto q=prim.array(),qa=aux.array(),vv=vol.array(),aa=area.array(),u=cons.array(),div=divergence.array();EOS eos;
  ParallelFor(wide,[&](int i,int j,int k){
   auto& cell=cells[(k%4+4)%4][(j%4+4)%4][(i%4+4)%4];
   for(int n=0;n<12;n++)q(i,j,k,n)=cell[n];
   double e,c;eos.properties(cell[0],cell+5,cell[4],e,c);
   u(i,j,k,URHO)=cell[0];double kin=0;
   for(int n=1;n<4;n++){u(i,j,k,n)=cell[0]*cell[n];kin+=cell[n]*cell[n];}
   u(i,j,k,UEDEN)=cell[0]*(e+.5*kin);u(i,j,k,UEINT)=cell[0]*e;
   for(int n=5;n<12;n++)u(i,j,k,n)=cell[0]*cell[n];
   q(i,j,k,QREINT)=cell[0]*e;qa(i,j,k,QC)=c;qa(i,j,k,QGAMC)=c*c*cell[0]/cell[4];
  });
  ParallelFor(bx,[&](int i,int j,int k){vv(i,j,k)=1.;});
  ParallelFor(fb,[&](int i,int j,int k){aa(i,j,k)=1.;});
  GpuArray<const Array4<double>,3> flx{fx.array(),fy.array(),fz.array()},qec{qx.array(),qy.array(),qz.array()};
  GpuArray<const Array4<const double>,3> areas{aa,aa,aa};
  int bc[3]={0,0,0},lo[3]={0,0,0},hi[3]={3,3,3};
  pc_umeth_3D(bx,bc,bc,lo,hi,q,qa,src.array(),flx,qec,areas,pd.array(),vv,{.25,.25,.25},dt,1,4,flattening!=0,false,0);
  ParallelFor(fb,[&](int i,int j,int k){pc_divu(i,j,k,q,.25,.25,.25,div);});
  for(int d=0;d<3;d++){
   auto out=flx[d];Box faces=surroundingNodes(bx,d);
   ParallelFor(faces,[&](int i,int j,int k){
    pc_artif_visc(i,j,k,out,div,u,.25,difmag,d,0,3,0,0);
    pc_norm_spec_flx(i,j,k,out);
    for(int n=0;n<12;n++)std::cout<<out(i,j,k,n)<<' ';
    std::cout<<'\n';
   });
  }
 }
}
'''


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--work-dir',type=Path,required=True)
    p.add_argument('--cxx',default='c++')
    p.add_argument('--shock-viscosity',action='store_true')
    a=p.parse_args();root=Path(__file__).resolve().parents[1]
    pin=json.loads((root/'references/pelec_baseline.json').read_text())['pelec']['commit']
    names=('PPM.H','PPM.cpp','Riemann.H','Godunov.H','Godunov.cpp','Hydro.H')
    sources={n:run(['git','-C',str(a.reference),'show',f'{pin}:Source/{n}']) for n in names}
    work=a.work_dir.resolve();work.mkdir(parents=True,exist_ok=True)
    for n,s in sources.items():(work/n).write_text(s)
    shim=(root/'tools/reference_shims/ctu_array.hpp').read_text()
    flatten=extract(sources['Godunov.H'],'flatten').replace('void flatten(', 'amrex::Real flatten(',1)
    bodies=flatten+''.join(extract(sources[n],f) for n,f in (
        ('PPM.H','ppm_reconstruct'),('PPM.H','ppm_int_profile'),
        ('Riemann.H','riemann'),('PPM.cpp','trace_ppm'),
        ('Godunov.H','pc_cmpflx_passive'),('Godunov.H','pc_cmpflx'),
        ('Godunov.H','pc_transdo_passive'),('Godunov.H','pc_transdd_passive'),
        ('Godunov.H','pc_transdo'),('Godunov.H','pc_transdd'),('Godunov.H','pc_pdivu'),
        ('Godunov.cpp','pc_umeth_3D'),('Hydro.H','pc_divu'),
        ('Godunov.H','pc_artif_visc'),('Hydro.H','pc_norm_spec_flx')))
    cpp=work/'reference.cpp';cpp.write_text(PREFIX+shim+bodies+DRIVER)
    command=[a.cxx,'-std=c++17','-O0','-ffp-contract=off',str(cpp),'-o',str(work/'reference')]
    # Preserve compiler diagnostics, including failures.
    compiled=subprocess.run(command,text=True,capture_output=True)
    (work/'compile.log').write_text(compiled.stdout+compiled.stderr)
    compiled.check_returncode()
    exe=str(a.executable.resolve());thermo=run([exe,'--thermo-stream'])
    (work/'thermo.txt').write_text(thermo)
    mw=[float(line.split()[0]) for line in thermo.splitlines()]
    cases=[]
    configs=[(mode,dt,0,0.) for mode in range(3) for dt in (0.,1e-5)]
    if a.shock_viscosity:configs += [(3,1e-5,flat,visc) for visc in (0.,.1) for flat in (0,1)]
    for mode,dt,flat,visc in configs:
        rows=[]
        for k in range(4):
            for j in range(4):
                for i in range(4):
                    phase=2*math.pi*(i+2*j+3*k)/4
                    y=[.01,.001,.001,.2,.001,.05,0.]
                    if mode==2:
                        y[0]+=.003*math.cos(phase)
                        y[3]+=.02*math.sin(phase+.3)
                    y[-1]=1-sum(y[:-1])
                    temp=900*(1+.02*mode*math.sin(phase))
                    pressure=1e5*(1+.01*mode*math.cos(phase))
                    r=sum(v*8314.46261815324/m for v,m in zip(y,mw))
                    velocity=[20+mode*math.sin(phase),-10+mode*math.cos(phase),5-.5*mode*math.sin(phase)]
                    if mode==3:
                        pressure=1e5 if i<2 else 8e5
                        temp=900
                        velocity=[40. if i<2 else -40.,10*math.sin(phase),5*math.cos(phase)]
                    rows.append([pressure/(r*temp)]+velocity+[pressure]+y)
        cases.append(f'{dt} {flat} {visc}\n'+''.join(' '.join(format(v,'.17g') for v in row)+'\n' for row in rows))
    inputs=''.join(cases);(work/'inputs.txt').write_text(inputs)
    ref=run([str(work/'reference')],input=thermo+inputs)
    actual=run([exe,'--ctu-stream'],input=inputs)
    (work/'reference.txt').write_text(ref);(work/'actual.txt').write_text(actual)
    expected=[[float(v) for v in line.split()] for line in ref.splitlines()]
    observed=[[float(v) for v in line.split()] for line in actual.splitlines()]
    assert len(expected)==len(observed)==len(cases)*240
    worst=0.
    for case in range(len(cases)):
        x=expected[case*240:(case+1)*240];y=observed[case*240:(case+1)*240]
        assert all(len(row)==12 and all(math.isfinite(v) for v in row) for row in x+y)
        scales=[max(1.,max(abs(row[n]) for row in x+y)) for n in range(12)]
        error=max(abs(v-w)/scales[n] for left,right in zip(x,y) for n,(v,w) in enumerate(zip(left,right)))
        assert error<=3e-12,(case,error)
        worst=max(worst,error)
    if a.shock_viscosity:
        for table in (expected,observed):
            assert table[6*240:7*240]!=table[7*240:8*240], 'flattening inactive'
            assert table[6*240:7*240]!=table[8*240:9*240], 'viscosity inactive'
    report=dict(shock_viscosity=a.shock_viscosity,reference_commit=pin,source_sha256={n:hashlib.sha256(s.encode()).hexdigest() for n,s in sources.items()},
                shim_sha256=hashlib.sha256(shim.encode()).hexdigest(),
                eos_adapter_sha256=hashlib.sha256(PREFIX.encode()).hexdigest(),
                compiler=run([a.cxx,'--version']).splitlines()[0],
                scaling='maximum absolute component flux across all faces per case, floored at one',
                executable_sha256=hashlib.sha256(a.executable.read_bytes()).hexdigest(),compile_command=command,
                cases=len(cases),values=len(expected)*12,max_scaled_error=worst,tolerance=3e-12,status='PASS')
    (work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
