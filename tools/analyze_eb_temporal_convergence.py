#!/usr/bin/env python3
"""Compare same-time EB composite fields across two or more time spacings.

Input JSON: {"times_s": [...], "runs": [{"label": str, "dt_s": float,
"application": directory, "reference": directory}, ...]}.
Runs must be ordered coarse to fine. Native reference tables include covered
coarse cells; the already validated application's composite keys select the
matching fluid cells. This analyzer does not replace field/inventory gates.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def norms(errors, volumes):
    total = math.fsum(volumes.values())
    worst = max(errors, key=lambda k: abs(errors[k]))
    return {"Linf": abs(errors[worst]),
            "L1": math.fsum(abs(v)*volumes[k] for k, v in errors.items())/total,
            "L2": math.sqrt(math.fsum(v*v*volumes[k] for k, v in errors.items())/total),
            "worst_cell": worst}


def analyze(spec):
    identities = {}

    def read(path):
        path = Path(path).resolve()
        identities[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open() as stream:
            return list(csv.DictReader(stream))

    runs = spec['runs']
    if len(runs) < 2 or any(not math.isfinite(r['dt_s']) or r['dt_s'] <= 0 for r in runs):
        raise ValueError('At least two finite positive time spacings required')
    if any(a['dt_s'] <= b['dt_s'] for a, b in zip(runs, runs[1:])):
        raise ValueError('Runs must be strictly coarse to fine')
    outputs = []
    for time in spec['times_s']:
        tables = []
        for run in runs:
            step = round(time/run['dt_s'])
            if abs(step*run['dt_s']-time) > 1e-19:
                raise ValueError('Requested time is not a saved step')
            app = Path(run['application'])/f'result.csv.step{step}.csv'
            clock = read(str(app)+'.time.csv')
            if len(clock) != 1 or abs(float(clock[0]['time_s'])-time) > 1e-19:
                raise ValueError('Application clock mismatch')
            ca = read(app)
            species = [k[5:] for k in ca[0] if k.startswith('rhoY_')]
            mapping = [('rho','density',1000),('rhou','xmom',10),('rhov','ymom',10),
                       ('rhow','zmom',10),('rhoE','rho_E',.1),
                       ('rho_internal_energy_J_m3','rho_e',.1),('temperature','Temp',1)]
            c, volumes = {}, {}
            for row in ca:
                key = tuple(int(row[k]) for k in ['level','i','j','k'])
                if key in c:
                    raise ValueError('Duplicate composite cell')
                c[key] = {a:float(row[a]) for a, _, _ in mapping}
                c[key].update({'Y_'+s:float(row['rhoY_'+s])/float(row['rho']) for s in species})
                volumes[key] = float(row['fluid_volume'])
            n = {}
            for row in read(Path(run['reference'])/f'reference_{step:03}.csv'):
                if abs(float(row['time'])-time) > 1e-19:
                    raise ValueError('Native clock mismatch')
                key = (int(row['level']), *(round(float(row[d])/float(row['d'+d])-.5)+1 for d in 'xyz'))
                if key not in c:
                    continue
                if key in n:
                    raise ValueError('Duplicate native cell')
                n[key] = {a:float(row[b])*scale for a,b,scale in mapping}
                n[key].update({'Y_'+s:float(row['rho_'+s])/float(row['density']) for s in species})
            if n.keys() != c.keys() or any(not math.isfinite(v) or v <= 0 for v in volumes.values()):
                raise ValueError('Composite coverage/volume mismatch')
            if any(not math.isfinite(v) for table in (c,n) for row in table.values() for v in row.values()):
                raise ValueError('Nonfinite field')
            tables.append((c,n,volumes))
        pairs = []
        for idx, (a,b) in enumerate(zip(tables,tables[1:])):
            if a[2] != b[2] or a[0].keys() != b[0].keys() or a[0][next(iter(a[0]))].keys() != b[0][next(iter(b[0]))].keys():
                raise ValueError('Mesh, volumes, or species changed between spacings')
            pair = {'coarse':runs[idx]['label'], 'fine':runs[idx+1]['label']}
            for solver, column in [('candidate',0),('native',1)]:
                pair[solver] = {field:norms({k:b[column][k][field]-a[column][k][field] for k in a[column]},a[2])
                                for field in next(iter(a[column].values()))}
            pairs.append(pair)
        orders = []
        for idx, (a,b) in enumerate(zip(pairs,pairs[1:])):
            ratio = runs[idx]['dt_s']/runs[idx+1]['dt_s']
            if not math.isclose(ratio,runs[idx+1]['dt_s']/runs[idx+2]['dt_s'],rel_tol=1e-12):
                raise ValueError('Observed order requires equal refinement ratios')
            orders.append({solver:{field:{norm:(math.log(a[solver][field][norm]/b[solver][field][norm])/math.log(ratio)
                                               if a[solver][field][norm]>0 and b[solver][field][norm]>0 else None)
                                         for norm in ['Linf','L1','L2']} for field in a[solver]} for solver in ['candidate','native']})
        outputs.append({'time_s':time,'pairs':pairs,'observed_orders':orders})
    return {'qualification':False,'scope':'volume-weighted same-time composite temporal self-differences; species are mass fractions; two spacings do not establish order',
            'runs':runs,'rows':outputs,'sha256':identities}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    result = analyze(json.loads(args.spec.read_text()))
    for path in [args.spec,Path(__file__)]:
        result['sha256'][str(path.resolve())] = hashlib.sha256(path.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
