#!/usr/bin/env python3
"""Replay recorded EB coarse-halo stencils and attribute shared-factor clipping.

Requires the probe's _initial_limiter.csv. Evaluates the documented limiter
independently, checks actual kernel outputs, and reports which momentum
component changes the scalar factor. This is an operator audit, not a solver
accuracy or fixed-reference acceptance test.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def analyze(path):
    stencils = {}
    with path.open() as stream:
        for row in csv.DictReader(stream):
            key = tuple(int(row[k]) for k in ('i', 'j', 'k'))
            c = int(row['component'])
            if c in stencils.setdefault(key, {}):
                raise ValueError('Duplicate stencil component')
            stencils[key][c] = {k: float(v) for k, v in row.items() if k not in ('i', 'j', 'k', 'component')}
    if not stencils:
        raise ValueError('No recorded stencils')
    directions = []; max_replay = 0.; changes = {}
    for key, rows in stencils.items():
        if set(rows) != set(range(1, 17)):
            raise ValueError('Expected gas(14), temperature, held-energy layout')
        offsets = [rows[1]['offset_' + a] for a in 'xyz']
        predicted = {c: [r['center'], r['center']] for c, r in rows.items()}
        for a, offset in zip('xyz', offsets):
            dc = {}; limited = {}; factors = {}
            for c, r in rows.items():
                left, center, right = r[a+'minus'], r['center'], r[a+'plus']
                dc[c] = .5 * (right-left)
                forward, backward = 2*(right-center), 2*(center-left)
                bound = min(abs(forward), abs(backward)) if forward*backward >= 0 else 0.
                limited[c] = math.copysign(min(bound, abs(dc[c])), dc[c])
                factors[c] = limited[c]/dc[c] if dc[c] != 0 else 1.
            common = min(factors.values())
            scalar = min(f for c, f in factors.items() if c not in (2, 3, 4))
            if scalar-common > 1e-14:
                drivers = [c for c in (2, 3, 4) if factors[c] == common]
                directions.append({'cell': key, 'axis': a, 'common_factor': common,
                                   'scalar_factor': scalar, 'momentum_components': drivers,
                                   'momentum_centered_slopes': {c: dc[c] for c in drivers}})
            for c in rows:
                predicted[c][0] += offset*common*dc[c]
                predicted[c][1] += offset*(limited[c] if c in (2,3,4) else scalar*dc[c])
        for c,r in rows.items():
            actual = [r['default_value'], r['separate_value']]
            for pred, value in zip(predicted[c], actual):
                error = abs(pred-value)/max(1., abs(value))
                max_replay = max(max_replay, error)
                if not math.isfinite(error) or error > 2e-14:
                    raise ValueError('Recorded kernel differs from limiter formula')
            delta = abs(actual[1]-actual[0])
            if c not in changes or delta > changes[c]['absolute']:
                changes[c] = {'absolute': delta, 'cell': key, 'default': actual[0], 'separate': actual[1]}
    return {'status':'passed', 'qualification':False,
            'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'stencils':len(stencils), 'directions':3*len(stencils),
            'changed_scalar_factor_directions':len(directions),
            'common_zero_scalar_nonzero':sum(d['common_factor']==0 and d['scalar_factor']>0 for d in directions),
            'maximum_scaled_kernel_replay_error':max_replay,
            'maximum_component_changes':changes, 'factor_changes':directions}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=analyze(args.input)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='factor_changes'},indent=2))


if __name__=='__main__': main()
