#!/usr/bin/env python3
"""Compare four native MOL halo policies with analytic steady shear cell averages.

Inputs supply geometry, transport fits and V4 controls; the probe replaces all
fluid states with pure N2 at 900 K, 101325 Pa and u=10*sin(2*pi*y/L).
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe',type=Path,required=True)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--maximum-grid-seconds',type=float,default=3600.)
    args=parser.parse_args()
    out=args.work.resolve();out.mkdir(parents=True,exist_ok=False)
    probe=out/'frozen-probe';shutil.copy2(args.probe,probe)
    record={'status':'preparing','qualification':False,'pid':os.getpid(),'cwd':str(Path.cwd()),
            'binary':str(probe),'binary_sha256':digest(probe),'runner_sha256':digest(Path(__file__)),
            'head':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
            'analytic_solution':'rho,T,p,Y constant; u=10 sin(2 pi y/L); v=w=0; exact cell-average momentum and total energy',
            'field_components':{'1':'density','2':'xmom','3':'ymom','4':'zmom','5':'rho_E',
                                '6-14':'species densities in LiDryer order','15':'held rho_e','16':'temperature','17':'pressure'},
            'dt_s':5e-8,'endpoint_s':2e-7,'runs':[]}
    def save():
        p=out/'record.tmp';p.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');p.replace(out/'record.json')
    save()
    try:
        for n in (8,16,32):
            if n==32:
                estimate=record['runs'][-1]['seconds']*12
                record['grid32_forecast_seconds']=estimate
                if estimate>args.maximum_grid_seconds:
                    record['status']='cost_limit';save();return
            source=args.inputs/f'grid{n}.dat';lines=source.read_text().splitlines()
            if lines[0]!='PELEF_REACTING_EB_REFERENCE_V4' or int(lines[1])!=n or int(lines[3])!=1:
                raise ValueError('Expected converged V4 at the requested resolution')
            lines[4]='4 2 1 5e-8'
            inp=out/f'grid{n}.dat';inp.write_text('\n'.join(lines)+'\n')
            prefix=out/f'grid{n}'
            cmd=[str(probe),str(inp),str(prefix),'--audit-shear']
            item={'cells':n,'argv':cmd,'input_sha256':digest(inp),'source':str(source.resolve()),
                  'source_sha256':digest(source),'status':'running'}
            record['runs'].append(item);record['status']='running';save()
            start=time.monotonic()
            with (out/f'grid{n}.log').open('w') as log:
                process=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT)
                item['pid']=process.pid;save()
                try: item['returncode']=process.wait(timeout=args.maximum_grid_seconds)
                except subprocess.TimeoutExpired:
                    process.kill();process.wait();raise
            item['seconds']=time.monotonic()-start
            if item['returncode']!=0:raise RuntimeError(f'grid{n} failed; inspect its log')
            for policy in ('default','momentum','velocity','paired'):
                with Path(str(prefix)+f'_{policy}_native.csv').open() as f:summary=list(csv.DictReader(f))
                with Path(str(prefix)+f'_{policy}_analytic.csv').open() as f:rows=list(csv.DictReader(f))
                if len(summary)!=5 or len(rows)!=85:raise ValueError('Incomplete trajectory')
                if not all(math.isfinite(float(r[k])) for r in rows for k in ('L1','L2','Linf')):
                    raise ValueError('Nonfinite norm')
                if any(float(r['Linf'])>1e-8 for r in rows if int(r['step'])==0):
                    raise ValueError('Initial analytic state mismatch')
                if abs(float(summary[-1]['time_s'])-2e-7)>1e-20:raise ValueError('Endpoint mismatch')
                budget=max(float(r['max_budget_residual']) for r in summary)
                if not math.isfinite(budget) or budget>1e-12:raise ValueError('Budget failed')
                item[policy]={'maximum_budget':budget,'endpoint':summary[-1],
                              'endpoint_norms':{r['component']:{k:float(r[k]) for k in ('L1','L2','Linf')} for r in rows if int(r['step'])==4}}
            item['status']='completed';save();print(f'grid{n}: completed in {item["seconds"]:.1f} s',flush=True)
        record['observed_L2_orders']={}
        for policy in ('default','momentum','velocity','paired'):
            record['observed_L2_orders'][policy]={}
            for c in range(1,18):
                errors=[r[policy]['endpoint_norms'][str(c)]['L2'] for r in record['runs']]
                record['observed_L2_orders'][policy][str(c)]=[math.log2(a/b) if a>0 and b>0 else None for a,b in zip(errors,errors[1:])]
        if digest(probe)!=record['binary_sha256']:raise ValueError('Binary drift')
        record['status']='completed'
    except Exception as exc:
        record['status']='failed';record['error']=str(exc);raise
    finally:save()


if __name__=='__main__':main()
