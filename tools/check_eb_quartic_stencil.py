#!/usr/bin/env python3
"""Isolated interpolation comparison, not an EB halo or CFD qualification."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import subprocess

PINNED_HEADER_SHA256 = '10cd49e2a8e880585cbd0d8956e27736d2620d01de5f063804cf84b7ef00296c'

PREFIX = r'''
#include <iostream>
#include <iomanip>
using Real=double;
constexpr double operator"" _rt(long double x){return double(x);}
template<class T,int L,int U> struct Array1D {T data[U-L+1]; constexpr T operator()(int i)const{return data[i-L];} T& operator()(int i){return data[i-L];}};
template<class T,int L,int U,int V,int W> struct Array2D {T data[(U-L+1)*(W-V+1)]; T& operator()(int i,int j){return data[i-L+(U-L+1)*(j-V)];}};
template<class T> struct Array4 {T* data; int lower,side; T& operator()(int i,int j,int k,int)const{return data[i-lower+side*(j-lower+side*(k-lower))];}};
namespace amrex {int coarsen(int i,int ratio){return i/ratio;}}
'''
DRIVER = r'''
int main(){double coarse[125],fine[8]; std::cout<<std::setprecision(17);
while(std::cin>>coarse[0]){for(int n=1;n<125;n++)if(!(std::cin>>coarse[n]))return 2;
Array4<const double> q{coarse,-2,5};Array4<double> f{fine,0,2};
for(int k=0;k<2;k++)for(int j=0;j<2;j++)for(int i=0;i<2;i++)ccquartic_interp(i,j,k,0,q,f);
for(double v:fine)std::cout<<v<<' ';std::cout<<'\n';}}
'''


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--header',type=Path,required=True)
    ap.add_argument('--work-dir',type=Path,required=True)
    args=ap.parse_args();b=args.work_dir;b.mkdir(parents=True,exist_ok=True)
    assert hashlib.sha256(args.header.read_bytes()).hexdigest()==PINNED_HEADER_SHA256, 'Pinned header mismatch'
    text=args.header.read_text();start=text.index('void ccquartic_interp (');opening=text.index('{',start)
    end=opening+1;depth=1
    while depth:
        depth+=(text[end]=='{')-(text[end]=='}');end+=1
    cpp=b/'reference.cpp';cpp.write_text(PREFIX+text[start:end]+DRIVER)
    root=Path(__file__).resolve().parents[1];source=root/'tools/reference_probes/eb_quartic_stencil_probe.F90'
    kernel=root/'tools/reference_probes/eb_quartic_kernel_mod.F90'
    commands=[['c++','-std=c++17','-O0','-ffp-contract=off',str(cpp),'-o',str(b/'reference')],
              ['gfortran','-std=f2018','-O0','-fcheck=all','-ffpe-trap=invalid,zero,overflow','-J'+str(b),'-I'+str(b),str(kernel),str(source),'-o',str(b/'probe')]]
    for command in commands:subprocess.run(command,check=True,capture_output=True,text=True)
    coordinates=list(itertools.product(range(-2,3),repeat=3)) # z,y,x; x contiguous
    children=list(itertools.product(range(2),repeat=3))
    cases=[]
    def avg(power,lo,hi):return (hi**(power+1)-lo**(power+1))/((power+1)*(hi-lo))
    for powers in [(p,0,0) for p in range(5)]+[(0,p,0) for p in range(1,5)]+[(0,0,p) for p in range(1,5)]+[(2,2,2),(3,1,4)]:
        values=[math.prod(avg(p,c-.5,c+.5) for p,c in zip(powers,xyz)) for xyz in coordinates]
        exact=[math.prod(avg(p,(c-1)*.5,c*.5) for p,c in zip(powers,xyz)) for xyz in children]
        cases.append({'name':f'polynomial_{powers}','values':values,'exact':exact})
    rng=random.Random(421)
    for n in range(100):cases.append({'name':f'random_{n}','values':[rng.uniform(-2,3) for _ in coordinates]})
    for n in [8,16,32]:
        width=2*math.pi/n;j=n//4
        cases.append({'name':f'shear_{n}','values':[10*math.sin((j+y-.5)*width)*math.sin(width/2)/(width/2) for z,y,x in coordinates]})
    cases.append({'name':'step','values':[float(y>0) for z,y,x in coordinates]})
    payload=''.join(' '.join(format(v,'.17g') for v in c['values'])+'\n' for c in cases)
    outputs=[]
    for exe in ['reference','probe']:
        result=subprocess.run([str((b/exe).resolve())],input=payload,capture_output=True,text=True,check=True)
        (b/(exe+'.txt')).write_text(result.stdout)
        outputs.append([[float(v) for v in line.split()] for line in result.stdout.splitlines()])
    assert len(outputs[0])==len(outputs[1])==len(cases)
    max_difference=max_conservation=max_polynomial=0.0;shear=[]
    for case,ref,observed in zip(cases,*outputs):
        assert len(ref)==len(observed)==8 and all(math.isfinite(v) for v in ref+observed)
        delta=max(abs(a-b) for a,b in zip(ref,observed));max_difference=max(max_difference,delta)
        assert delta<1e-12*max(1,max(abs(v) for v in ref)),case['name']
        residual=abs(sum(observed)/8-case['values'][62]);max_conservation=max(max_conservation,residual)
        assert residual<1e-12*max(1,abs(case['values'][62])),case['name']
        if 'exact' in case:
            error=max(abs(a-b) for a,b in zip(observed,case['exact']));max_polynomial=max(max_polynomial,error)
            assert error<1e-12,case['name']
        if case['name'].startswith('shear_'):
            n=int(case['name'].split('_')[1]);width=2*math.pi/n
            exact=10*math.sin((n/2-.5)*width/2)*math.sin(width/4)/(width/4)
            q=case['values'];center=q[62];minus=q[57];plus=q[67]
            centered=center+.125*(plus-minus)
            shear.append({'root_cells':n,'exact_fine_average':exact,'limited_error':abs(center-exact),
                          'centered_error':abs(centered-exact),'quartic_error':abs(observed[2]-exact)})
    step=outputs[1][-1];assert min(step)<0 and max(step)>0
    report={'status':'passed_arithmetic_checks','qualification':False,'production_candidate_admitted':False,
            'toolchain':{c:subprocess.check_output([c,'--version'],text=True).splitlines()[0] for c in ['c++','gfortran']},
            'head':subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),
            'cases':len(cases),'commands':commands,'max_reference_difference':max_difference,
            'max_child_mean_residual':max_conservation,'max_polynomial_error':max_polynomial,'shear':shear,
            'step_child_min':min(step),'step_child_max':max(step),
            'limitation':'Unlimited quartic creates negative children at a nonnegative jump; needs admissibility and discontinuity handling',
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),args.header,kernel,source,cpp,b/'reference',b/'probe']}}
    (b/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['commands','sha256']},indent=2))


if __name__=='__main__':main()
