#!/usr/bin/env python3
"""Validate a hydro-only paired-pressure Euler/reflux diagnostic record.

No chemistry, SDC, repeated-step, MPI, or spatial convergence qualification.
"""
import argparse
import csv
import json
import math
from pathlib import Path
from audit_eb_initial_interface import sha


def inventory_check(rows):
    values={r['component']:(float(r['initial']),float(r['final'])) for r in rows}
    if len(values)!=len(rows) or not {'1','5','held'}.issubset(values):
        raise ValueError('Invalid inventory components')
    if not all(math.isfinite(x) for pair in values.values() for x in pair):
        raise ValueError('Nonfinite inventory')
    mass=values['1'][0]
    if mass<=0:
        raise ValueError('Nonpositive initial mass')
    residuals={}
    for component,(initial,final) in values.items():
        if component in ('2','3','4'):
            continue  # Momentum requires external wall forces, not zero budget.
        scale=max(abs(initial),mass*1e-14)
        residuals[component]=abs(final-initial)/scale
    if max(residuals.values())>1e-12:
        raise AssertionError('Composite gas/species/held inventory failed')
    species=[pair for key,pair in values.items() if key!='held' and int(key)>=6]
    if not species:
        raise ValueError('Missing species inventories')
    for time in (0,1):
        if abs(math.fsum(pair[time] for pair in species)-values['1'][time])/mass>1e-12:
            raise AssertionError('Species sum does not close mass')
    return residuals


def check(record_path):
    record=json.loads(record_path.read_text())
    if record.get('status')!='completed' or not record.get('commands'):
        raise ValueError('Incomplete run record')
    result={}
    for command in record['commands']:
        binary,source,prefix=map(Path,command['argv'][:3])
        if command['returncode'] or sha(binary)!=record['binary_sha256'] or sha(source)!=command['input_sha256']:
            raise ValueError('Failed run or changed identity')
        cases={}
        for variant in ('baseline','paired'):
            path=Path(str(prefix)+f'_{variant}_pressure_step.csv')
            with path.open() as f:rows=list(csv.DictReader(f))
            if len(rows)!=1:raise ValueError('Expected one step')
            values={k:float(v) for k,v in rows[0].items()}
            if not all(math.isfinite(v) for v in values.values()) or values['dt_s']<=0:
                raise ValueError('Invalid step diagnostics')
            inventory=Path(str(prefix)+f'_{variant}_pressure_step_inventory.csv')
            with inventory.open() as f:residuals=inventory_check(list(csv.DictReader(f)))
            cases[variant]={'diagnostics':values,'independent_inventory_errors':residuals,
                            'sha256':{str(path):sha(path),str(inventory):sha(inventory)}}
        paired=cases['paired']['diagnostics']
        if any(command['velocity_m_s']):
            if paired['max_state_change']<=1e-6:
                raise AssertionError('Moving control did not evolve')
        elif paired['max_pressure_offset_Pa']>1e-8:
            raise AssertionError('Stationary paired equilibrium failed')
        if paired['max_E_minus_H_minus_KE_J_m3']>1e-6:
            raise AssertionError('Held/total energy inconsistency')
        result[source.stem]=cases
    return {'status':'passed','qualification':False,'cases':result,'run_sha256':sha(record_path)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(check(args.record),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
