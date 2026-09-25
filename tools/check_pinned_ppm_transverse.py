#!/usr/bin/env python3
"""Compare regular CTU transverse updates to pinned C++ bodies."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import random
from check_pinned_ppm_kernels import run

PREFIX=r'''
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#define AMREX_D_DECL(a,b,c) a,b,c
#define NUM_ADV 0
#define NUM_LIN 0
constexpr int NUM_SPECIES=2;
constexpr int QRHO=0,QU=1,QV=2,QW=3,QPRES=4,QFS=5,QREINT=7,QGAMC=8;
constexpr int URHO=0,UMX=1,UMY=2,UMZ=3,UEDEN=4,UFS=5;
constexpr int GDU=1,GDV=2,GDW=3,GDPRES=4;
namespace constants { double very_small_num(){return std::numeric_limits<double>::epsilon()*1e-100;} }
namespace amrex {
using Real=double; using std::max; using std::min;
template<class T,int N> using GpuArray=std::array<T,N>;
struct IntVect {
 int x,y,z;
 static IntVect TheDimensionVector(int d){return {d==0,d==1,d==2};}
 IntVect operator+(IntVect b) const {return {x+b.x,y+b.y,z+b.z};}
};
template<class T> struct Array4 {
 T* data=nullptr;
 T& operator()(IntVect v,int n=0) const {return data[(((v.z*2)+v.y)*2+v.x)*9+n];}
 explicit operator bool() const {return data!=nullptr;}
};
}
'''
DRIVER=r'''
int main(){
 int nd,dir; double forcing[7],source_e,hdt,normal[7],p,gamma,flux[2][7],pl[2],ph[2],vl[2],vh[2],scale[2];
 std::cout<<std::setprecision(17);
 while(std::cin>>nd>>dir){
  for(double& v:normal)std::cin>>v;
  std::cin>>p>>gamma;
  for(auto& row:flux)for(double& v:row)std::cin>>v;
  for(double& v:pl)std::cin>>v;for(double& v:ph)std::cin>>v;
  for(double& v:vl)std::cin>>v;for(double& v:vh)std::cin>>v;
  for(double& v:scale)std::cin>>v;
  for(double& v:forcing)std::cin>>v;std::cin>>source_e>>hdt;
  if(!std::cin)return 2;
  double normals[72]{},f0[72]{},f1[72]{},q0[72]{},q1[72]{},qa[72]{},src[72]{},outp[72]{},outm[72]{};
  using amrex::Array4;using amrex::IntVect;
  Array4<double> nn{normals},ff0{f0},ff1{f1},qq0{q0},qq1{q1},aa{qa},pp{outp},mm{outm};
  IntVect iv{0,0,0},ivn=IntVect::TheDimensionVector(dir);
  int t0=dir==0?1:0,t1=dir==2?1:2;
  double rho=normal[0],kin=.5*(normal[1]*normal[1]+normal[2]*normal[2]+normal[3]*normal[3])/rho;
  for(IntVect v:{iv,ivn}){
   nn(v,QRHO)=rho;for(int n=1;n<4;n++)nn(v,n)=normal[n]/rho;
   nn(v,QPRES)=p;nn(v,QREINT)=normal[4]-kin;
   for(int n=5;n<7;n++)nn(v,n)=normal[n]/rho;
  }
  for(int n=0;n<7;n++)src[n]=forcing[n];src[QREINT]=source_e;
  aa(iv,QGAMC)=gamma;
  for(int t=0;t<2;t++){
   int axis=t==0?t0:t1;IntVect high=IntVect::TheDimensionVector(axis);
   auto ff=t==0?ff0:ff1;auto qq=t==0?qq0:qq1;
   for(int n=0;n<7;n++)ff(high,n)=flux[t][n];
   qq(iv,GDPRES)=pl[t];qq(high,GDPRES)=ph[t];
   qq(iv,axis+1)=vl[t];qq(high,axis+1)=vh[t];
  }
  if(nd==1)pc_transdo(0,0,0,dir,t0,mm,pp,{normals},{normals},{f0},{qa},{q0},scale[0]);
  else pc_transdd(0,0,0,dir,mm,pp,{normals},{normals},{f0},{f1},{q0},{q1},{qa},{src},hdt,scale[0],scale[1]);
  for(int n=0;n<8;n++)std::cout<<pp(iv,n)<<' ';
  for(int n=0;n<8;n++)std::cout<<mm(ivn,n)<<' ';
  std::cout<<'\n';
 }
}
'''


def extract(source: str,name: str)->str:
    start=source.index('\n'+name+'(')+1
    opening=source.index('\n{',start)+1
    depth=1;end=opening+1
    while depth:
        depth+=(source[end]=='{')-(source[end]=='}');end+=1
    return 'void '+source[start:end]+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--cxx',default='c++')
    a=parser.parse_args();root=Path(__file__).resolve().parents[1]
    pin=json.loads((root/'references/pelec_baseline.json').read_text())['pelec']['commit']
    source=run(['git','-C',str(a.reference),'show',f'{pin}:Source/Godunov.H'])
    work=a.work_dir.resolve();work.mkdir(parents=True,exist_ok=True)
    (work/'Godunov.H').write_text(source)
    cpp=work/'reference.cpp'
    cpp.write_text(PREFIX+''.join(extract(source,n) for n in
        ('pc_transdo_passive','pc_transdd_passive','pc_transdo','pc_transdd'))+DRIVER)
    command=[a.cxx,'-std=c++17','-O0','-ffp-contract=off',str(cpp),'-o',str(work/'reference')]
    run(command)
    rng=random.Random(20260926);rows=[]
    for nd in (1,2):
        for axis in range(3):
            for mode in ('regular','signed','floor'):
                for _ in range(8):
                    rho=rng.uniform(.5,2);vel=[rng.uniform(-30,30) for _ in range(3)]
                    normal=[rho]+[rho*v for v in vel]+[1e5+.5*rho*sum(v*v for v in vel)]+[.3*rho,.7*rho]
                    flux=[[rng.uniform(-.1,.1) for _ in range(7)] for _ in range(2)]
                    if mode=='signed':flux[0][5]=100;flux[0][6]=-100
                    pl=[rng.uniform(8e4,1e5) for _ in range(2)];ph=[rng.uniform(8e4,1e5) for _ in range(2)]
                    vl=[rng.uniform(-3,3) for _ in range(2)];vh=[rng.uniform(-3,3) for _ in range(2)]
                    if mode=='floor':vl=[-1e4,-1e4];vh=[1e4,1e4]
                    rows.append([nd,axis]+normal+[1e5,1.4]+sum(flux,[])+pl+ph+vl+vh+[.01,.02])
    rows=[row+[0.]*9 for row in rows]
    # Repeat all double-transverse fixtures with finite primitive source rates.
    rows += [row[:-9]+[.2,3.,-4.,2.,-2e4,.1,-.1,5e3,.025] for row in rows if row[0]==2]
    inputs=''.join(' '.join(format(v,'.17g') for v in row)+'\n' for row in rows)
    (work/'inputs.txt').write_text(inputs)
    reference=run([str(work/'reference')],input=inputs)
    actual=run([str(a.executable.resolve()),'--transverse-stream'],input=inputs)
    (work/'reference.txt').write_text(reference);(work/'actual.txt').write_text(actual)
    expected=[[float(v) for v in line.split()] for line in reference.splitlines()]
    observed=[[float(v) for v in line.split()] for line in actual.splitlines()]
    assert len(expected)==len(observed)==len(rows)
    worst=0.
    for i,(x,y) in enumerate(zip(expected,observed)):
        assert len(x)==len(y)==16 and all(math.isfinite(v) for v in x+y)
        error=max(abs(v-w)/max(1.,abs(v),abs(w)) for v,w in zip(x,y))
        assert error<=3e-12,(i,error,x,y)
        # Pressure floor is tiny: require exact agreement, not a unit-scaled pass.
        if i%24>=16:assert x[4]==y[4] and x[12]==y[12],(i,x,y)
        worst=max(worst,error)
    report=dict(reference_commit=pin,source_sha256=hashlib.sha256(source.encode()).hexdigest(),
                executable_sha256=hashlib.sha256(a.executable.read_bytes()).hexdigest(),compile_command=command,
                compiler=run([a.cxx,'--version']).splitlines()[0],cases=len(rows),values=16*len(rows),
                max_scaled_error=worst,tolerance=3e-12,status='PASS')
    (work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
