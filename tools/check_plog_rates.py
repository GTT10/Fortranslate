#!/usr/bin/env python3
"""Compare imported/generated PLOG rates and concentration derivatives with Cantera."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    import cantera as ct
    import numpy as np
    parser=argparse.ArgumentParser()
    parser.add_argument('--build-dir',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path,required=True)
    args=parser.parse_args();args.work_dir.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='plog-',dir=args.work_dir.resolve()))
    root=Path(__file__).resolve().parents[1]
    parent=ct.Solution('h2o2.yaml')
    entries=[{'P':p,'A':a,'b':b,'Ea':e} for p,a,b,e in
             [(1e4,1e9,-.2,2e6),(1e5,2e9,.1,1e6),(1e5,-2e8,.1,1e6),(1e6,5e10,.3,3e6)]]
    reaction=ct.Reaction.from_dict({'equation':'H + OH <=> H2O','type':'pressure-dependent-Arrhenius',
                                    'rate-constants':entries},parent)
    gas=ct.Solution(thermo='ideal-gas',kinetics='gas',species=parent.species(),
                    reactions=[reaction],transport_model='mixture-averaged')
    gas.name='plog';gas.write_yaml(work/'plog.yaml')
    subprocess.run([sys.executable,str(root/'tools/ingest_cantera_mechanism.py'),
        '--input',str(work/'plog.yaml'),'--output',str(work/'plog.json'),'--phase','plog',
        '--module-name','plog_bundle_mod','--loader-name','load_plog'],check=True)
    subprocess.run([sys.executable,str(root/'tools/generate_elementary_mechanism.py'),
        '--input',str(work/'plog.json'),'--output',str(work/'plog.F90')],check=True)
    (work/'probe.F90').write_text('''program probe
 use precision_mod, only: dp
 use elementary_kinetics_mod
 use plog_bundle_mod, only: load_plog
 implicit none
 type(elementary_reaction),allocatable :: r(:)
 real(dp) :: t,c(10),k,dk(10)
 integer :: ios
 logical :: ok
 call load_plog(r,ok)
 if (.not.ok) error stop 'load'
 do
   read(*,*,iostat=ios) t,c
   if(ios/=0) exit
   call reaction_effective_forward_rate(r(1),t,c,k,dk,ok)
   if(.not.ok) error stop 'rate'
   write(*,'(*(es26.17e3,1x))') k,dk
 end do
 r(1)%plog(2,2:3)=-abs(r(1)%plog(2,2:3))
 call reaction_effective_forward_rate(r(1),1000.0_dp,c,k,dk,ok)
 if(ok) error stop 'nonpositive group accepted'
end program
''')
    subprocess.run(['gfortran','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
        '-I'+str(args.build_dir.resolve()/'modules'),str(work/'plog.F90'),str(work/'probe.F90'),
        str(args.build_dir.resolve()/'libpelef_core.a'),'-o',str(work/'probe')],cwd=work,check=True)
    inputs=[];expected=[];derivatives=[]
    for t in [500.,1000.,2000.]:
        for p in [1e3,1e4,2e4,7e4,1e5,2e5,7e5,1e6,1e7]:
            gas.TPX=t,p,'N2:1';m=sum(gas.concentrations)
            inputs.append(' '.join(map(str,[t,*gas.concentrations])))
            rate=gas.forward_rate_constants[0];expected.append(rate)
            if p in [1e4,1e5,1e6]:
                derivatives.append(None);continue  # pressure knots need not be differentiable
            samples=[];h=1e-3
            for shift in [-2*h,-h,h,2*h]:
                gas.TD=t,m*(1+shift)*gas.mean_molecular_weight
                samples.append(gas.forward_rate_constants[0])
            derivatives.append((samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*h*m))
    run=subprocess.run([str(work/'probe')],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    (work/'rates.txt').write_text(run.stdout)
    rows=np.array([list(map(float,line.split())) for line in run.stdout.splitlines()])
    rate_error=float(np.max(abs(rows[:,0]/expected-1)));derivative_error=0.
    for idx,(row,ref) in enumerate(zip(rows,derivatives)):
        if ref is None:continue
        if idx in [0,8,9,17,18,26]:continue
        derivative_error=max(derivative_error,abs(row[1]/ref-1))
    for idx in [0,8,9,17,18,26]: assert np.all(rows[idx,1:]==0)
    assert rate_error<1e-12,rate_error
    assert derivative_error<1e-6,derivative_error
    result={'states':len(inputs),'rate_relative_error':rate_error,
            'derivative_relative_error':float(derivative_error),'artifacts':str(work)}
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
