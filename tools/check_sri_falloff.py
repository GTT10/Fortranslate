#!/usr/bin/env python3
"""Cantera pressure-dependent rate and derivative oracle through generated Fortran."""
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
    parser.add_argument('--activated',action='store_true')
    parser.add_argument('--model',choices=['sri','troe','troe-lu','lindemann','tsang'],default='sri')
    parser.add_argument('--build-dir',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path,required=True)
    args=parser.parse_args()
    args.work_dir.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='sri-',dir=args.work_dir.resolve()))
    root=Path(__file__).resolve().parents[1]
    parent=ct.Solution('h2o2.yaml')
    data={'equation':'H + OH (+M) <=> H2O (+M)',
        'type':'chemically-activated' if args.activated else 'falloff',
        'low-P-rate-constant':{'A':2e12,'b':-1.,'Ea':1e6},
        'high-P-rate-constant':{'A':3e9,'b':.2,'Ea':2e6},
        'SRI':{'A':1.2,'B':500.,'C':1000.,'D':.8,'E':.1}}
    if args.model != 'sri': del data['SRI']
    if args.model == 'troe-lu': data['Troe']={'A':1.9816,'T3':5383.7,'T1':4.2932,'T2':-.0795}
    if args.model == 'troe': data['Troe']={'A':.6,'T3':100.,'T1':1000.,'T2':5000.}
    if args.model == 'tsang': data['Tsang']={'A':.4,'B':1e-4}
    reaction=ct.Reaction.from_dict(data,parent)
    gas=ct.Solution(thermo='ideal-gas',kinetics='gas',species=parent.species(),reactions=[reaction],transport_model='mixture-averaged')
    gas.name='sri';gas.write_yaml(work/'sri.yaml')
    subprocess.run([sys.executable,str(root/'tools/ingest_cantera_mechanism.py'),
        '--input',str(work/'sri.yaml'),'--output',str(work/'sri.json'),'--phase','sri',
        '--module-name','sri_bundle_mod','--loader-name','load_sri'],check=True)
    subprocess.run([sys.executable,str(root/'tools/generate_elementary_mechanism.py'),
        '--input',str(work/'sri.json'),'--output',str(work/'sri.F90')],check=True)
    (work/'probe.F90').write_text('''program probe
 use precision_mod, only: dp
 use elementary_kinetics_mod
 use sri_bundle_mod, only: load_sri
 implicit none
 type(elementary_reaction),allocatable :: r(:)
 real(dp) :: t,c(10),k,dk(10)
 integer :: ios
 logical :: ok
 call load_sri(r,ok)
 if (.not.ok) error stop 'load'
 do
   read(*,*,iostat=ios) t,c
   if(ios/=0) exit
   call reaction_effective_forward_rate(r(1),t,c,k,dk,ok)
   if(.not.ok) error stop 'rate'
   write(*,'(*(es26.17e3,1x))') k,dk
 end do
end program
''')
    if args.model == 'tsang':
        probe = work/'probe.F90'
        probe.write_text(probe.read_text().replace('end program', """r(1)%tsang=[.4_dp,-.001_dp]
 call reaction_effective_forward_rate(r(1),1000.0_dp,c,k,dk,ok)
 if(ok) error stop 'nonpositive Tsang center accepted'
end program"""))
    subprocess.run(['gfortran','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
        '-I'+str(args.build_dir.resolve()/'modules'),str(work/'sri.F90'),str(work/'probe.F90'),
        str(args.build_dir.resolve()/'libpelef_core.a'),'-o',str(work/'probe')],cwd=work,check=True)
    inputs=[];expected=[]
    n2=gas.species_index('N2')
    for t in [500.,1000.,2000.]:
        for pr in [1e-6,1e-3,1.,1e3,1e6]:
            m=pr*reaction.rate.high_rate(t)/reaction.rate.low_rate(t)
            gas.TPX=t,ct.one_atm,'N2:1';gas.TD=t,m*gas.mean_molecular_weight
            inputs.append(' '.join(map(str,[t,*gas.concentrations])))
            rate=gas.forward_rate_constants[0]
            # Fourth-order stencil permits a larger interval at near-flat low-pressure limits.
            h=1e-3; samples=[]
            for shift in [-2*h,-h,h,2*h]:
                gas.TD=t,m*(1+shift)*gas.mean_molecular_weight
                samples.append(gas.forward_rate_constants[0])
            expected.append((rate,(samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*h*m)))
    result=subprocess.run([str(work/'probe')],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    (work/'rates.txt').write_text(result.stdout)
    rows=np.array([list(map(float,line.split())) for line in result.stdout.splitlines()])
    expected=np.array(expected)
    rate_error=float(np.max(abs(rows[:,0]/expected[:,0]-1)))
    derivative_error=float(np.max(abs(rows[:,n2+1]/expected[:,1]-1)))
    assert rate_error<1e-12,(rate_error,derivative_error)
    assert derivative_error<1e-6,(rate_error,derivative_error)
    summary={'model':args.model,'activated':args.activated,'states':len(inputs),'rate_relative_error':rate_error,'derivative_relative_error':derivative_error,
             'artifacts':str(work)}
    (work/'result.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
