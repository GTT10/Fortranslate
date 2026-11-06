#!/usr/bin/env python3
"""Compare ghost-only pressure reconstruction with the unchanged initial operator."""
import argparse
import json
from pathlib import Path
from audit_eb_initial_interface import summarize, sha


def compare(directory):
    original=json.loads((directory/'report.json').read_text())
    if original.get('status')!='completed':
        raise ValueError('Baseline diagnostic is incomplete')
    cases={}
    for name, baseline in original['cases'].items():
        # The audit input header supplies the grid; no implicit eight-cell mesh.
        input_path=directory/('_'.join(name.split('_')[:-1])+'.dat')
        n=int(input_path.read_text().splitlines()[1])
        prefix=directory/(name+'_pressure_candidate')
        candidate=summarize(prefix,n)
        if candidate['max_area_fraction_mismatch']>1e-12:
            raise AssertionError('Candidate open areas disagree')
        if '_uniform_' in name:
            if candidate['max_interface_mass_flux_kg_m2_s']>1e-10:
                raise AssertionError('Candidate uniform mass flux failed')
            if candidate['max_area_weighted_flux_mismatch']['energy_W_m2']>1e-5:
                raise AssertionError('Candidate uniform energy flux failed')
        keys=('max_interface_pressure_jump_Pa','max_interface_mass_flux_kg_m2_s')
        cases[name]={'baseline':{k:baseline[k] for k in keys},
                     'candidate':{k:candidate[k] for k in keys},
                     'baseline_flux_mismatch':baseline['max_area_weighted_flux_mismatch'],
                     'candidate_flux_mismatch':candidate['max_area_weighted_flux_mismatch'],
                     'files_sha256':{str(Path(str(prefix)+suffix)):sha(str(prefix)+suffix)
                                     for suffix in ('_initial_pressure.csv','_initial_flux.csv')}}
    return {'status':'completed','qualification':False,
            'scope':'ghost-only pressure interpolation; valid-state immutability asserted by probe; no reflux or time integration qualification',
            'baseline_report_sha256':sha(directory/'report.json'),'cases':cases}


def compare_record(record_path):
    record=json.loads(record_path.read_text())
    if record.get('status')!='completed':
        raise ValueError('Probe run is incomplete')
    cases={}
    for command in record['commands']:
        if command['returncode']!=0:
            raise ValueError('Failed probe run')
        binary, source, prefix=map(Path,command['argv'][:3])
        if sha(binary)!=record['binary_sha256'] or sha(source)!=command['input_sha256']:
            raise ValueError('Probe or input identity changed')
        n=int(source.read_text().splitlines()[1])
        variants={}
        for label,suffix in (('baseline',''),('fine_only','_pressure_candidate'),('paired','_pressure_pair_candidate')):
            result=summarize(Path(str(prefix)+suffix),n)
            if result['max_area_fraction_mismatch']>1e-12:
                raise AssertionError('Coarse/fine areas disagree')
            variants[label]={k:v for k,v in result.items() if k!='faces'}
        paired=variants['paired']
        # Acceptance for these stationary isobaric controls, not arbitrary flows.
        if paired['max_interface_pressure_jump_Pa']>1e-8:
            raise AssertionError('Paired pressure equilibrium failed')
        if paired['max_interface_mass_flux_kg_m2_s']>1e-10:
            raise AssertionError('Paired stationary mass flux failed')
        limits={'mass_kg_m2_s':1e-10,'normal_momentum_Pa':1e-8,'energy_W_m2':1e-5}
        if any(paired['max_area_weighted_flux_mismatch'][k]>v for k,v in limits.items()):
            raise AssertionError('Paired interface flux equilibrium failed')
        cases[source.stem]=variants
    return {'status':'completed','qualification':False,'cases':cases,
            'scope':'stationary isobaric initial spatial operator; no evolved conservation/convergence qualification',
            'run_record_sha256':sha(record_path)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--directory',type=Path)
    source.add_argument('--run-record',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=compare_record(args.run_record) if args.run_record else compare(args.directory)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
