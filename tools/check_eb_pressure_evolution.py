#!/usr/bin/env python3
"""Check composite budgets at every accepted hydro-only SSPRK2 step."""
import argparse
import csv
import json
import math
from pathlib import Path
from check_eb_pressure_step import inventory_check
from audit_eb_initial_interface import sha


def validate(prefix, steps, dt, stationary=False):
    cases={}
    for policy in ('baseline','paired'):
        path=Path(str(prefix)+f'_{policy}_evolution.csv')
        with path.open() as f:rows=list(csv.DictReader(f))
        if len(rows)!=steps+1:
            raise ValueError('Incomplete evolution history')
        for it,row in enumerate(rows):
            if int(row['step'])!=it or not math.isclose(float(row['time_s']),it*dt,rel_tol=1e-13,abs_tol=1e-20):
                raise ValueError('Wrong step clock')
            if not all(math.isfinite(float(v)) for v in row.values()):
                raise ValueError('Nonfinite evolution history')
        invpath=Path(str(prefix)+f'_{policy}_evolution_inventory.csv')
        with invpath.open() as f:inv=list(csv.DictReader(f))
        by_step={}
        for row in inv:
            step=int(row['step']);component=row['component']
            if component in by_step.setdefault(step,{}):raise ValueError('Duplicate inventory')
            by_step[step][component]=row['value']
        if set(by_step)!=set(range(steps+1)):raise ValueError('Incomplete inventory clocks')
        start=by_step[0];worst=0
        for it in range(steps+1):
            if set(by_step[it])!=set(start):raise ValueError('Incomplete components')
            correction=float(rows[it].get('numerical_energy_correction_J',0))
            absolute=float(rows[it].get('absolute_energy_correction_J',abs(correction)))
            if absolute+1e-24<abs(correction):
                raise AssertionError('Invalid signed/absolute correction ledger')
            if it==0 and correction!=0:
                raise AssertionError('Initial correction ledger is nonzero')
            checked=[{'component':c,'initial':v,
                      'final':str(float(by_step[it][c])-correction) if c=='5' else by_step[it][c]}
                     for c,v in start.items()]
            residuals=inventory_check(checked)
            worst=max(worst,max(residuals.values()))
            if 'kinetic_energy_J' in rows[it]:
                energy=float(by_step[it]['5']);held=float(by_step[it]['held'])
                kinetic=float(rows[it]['kinetic_energy_J']);gap=float(rows[it]['signed_E_H_KE_gap_J'])
                if abs(energy-held-kinetic-gap)>1e-12*max(abs(energy),abs(held)):
                    raise AssertionError('Thermal/kinetic decomposition does not close')
        if stationary and policy=='paired' and max(float(r['max_pressure_offset_Pa']) for r in rows)>1e-8:
            raise AssertionError('Stationary equilibrium drift')
        cases[policy]={'endpoint':{k:float(v) for k,v in rows[-1].items()},
                       'worst_inventory_relative_error':worst,
                       'sha256':{str(path):sha(path),str(invpath):sha(invpath)}}
    return {'status':'passed','qualification':False,'steps':steps,'dt':dt,'cases':cases}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prefix',type=Path,required=True);p.add_argument('--input',type=Path,required=True)
    p.add_argument('--stationary',action='store_true');p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();lines=a.input.read_text().splitlines()
    if lines[0]!='PELEF_REACTING_EB_REFERENCE_V4':raise ValueError('Expected V4 input')
    row=lines[4].split();result=validate(a.prefix,int(row[0]),float(row[3]),a.stationary)
    result['input_sha256']=sha(a.input)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
