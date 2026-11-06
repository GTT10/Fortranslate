"""Report particle differences for the four-parcel evaporation fixture.

This is a diagnostic, not a parity acceptance test. Reference coordinates,
velocities and diameter use CGS; candidate CSV fields use SI.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference-gas', type=Path, help='CGS plotfile-export CSV for this periodic cube fixture')
    parser.add_argument('--candidate-gas', type=Path, help='SI EB tree CSV for this fixture')
    parser.add_argument('--expected-time', type=float, default=2e-7,
                        help='Expected reference gas endpoint in seconds (default: 2e-7)')
    parser.add_argument('--cells', type=int, choices=(4, 8, 16), default=4,
                        help='Expected cells per side in the 4 mm cube')
    args = parser.parse_args()
    if not math.isfinite(args.expected_time) or args.expected_time <= 0:
        parser.error('--expected-time must be finite and positive')
    if bool(args.reference_gas) != bool(args.candidate_gas):
        parser.error('Both gas CSV paths are required together')
    lines = args.reference.read_text().splitlines()
    if list(map(int, lines[:5])) != [4, 11, 0, 0, 0]:
        raise ValueError('Expected four-parcel, single-fuel PeleC ASCII format')
    reference = {}
    for line in lines[5:]:
        values = list(map(float, line.split()))
        if len(values) != 16 or not all(map(math.isfinite, values)):
            raise ValueError('Invalid reference particle row')
        key = int(values[14])
        if key in reference or values[14] != key:
            raise ValueError('Duplicate or invalid reference particle ID')
        reference[key] = dict(zip(('x', 'y', 'z', 'u', 'v', 'w'),
                                 (v * .01 for v in values[:6])))
        reference[key].update(temperature=values[6],
                              mass=681.4 * math.pi / 6 * (values[7] * .01)**3,
                              multiplicity=values[9])
    with args.candidate.open() as stream:
        rows = list(csv.DictReader(stream))
    candidate = {int(row['id']): row for row in rows}
    if len(rows) != 4 or set(candidate) != set(reference) or len(reference) != 4:
        raise ValueError('Particle IDs/counts differ')
    initial_mass = 681.4 * math.pi / 6 * (20e-6)**3
    initial = dict(x=.0015, y=.0015, z=.0015, u=0, v=0, w=2,
                   temperature=310, mass=initial_mass,
                   multiplicity=1e-10 / (4 * initial_mass))
    particles = {}
    for key, ref in reference.items():
        result = {}
        for field, value in ref.items():
            actual = float(candidate[key][field])
            if not math.isfinite(actual):
                raise ValueError('Nonfinite candidate value')
            change = value - initial[field]
            result[field] = dict(reference=value, candidate=actual,
                                 difference=actual-value,
                                 reference_change=change,
                                 candidate_change=actual-initial[field],
                                 relative_change_error=(actual-value)/abs(change)
                                 if abs(change) > 1e-25 else None)
        particles[key] = result
    report = dict(
        scope='400 K gas / 310 K liquid fixture; diagnostic only, no parity gate',
        inputs={name: dict(path=str(path.resolve()),
                           sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                for name, path in [('reference', args.reference), ('candidate', args.candidate)]},
        particles=particles)
    if args.reference_gas:
        with args.reference_gas.open() as stream:
            reference_gas = list(csv.DictReader(stream))
        with args.candidate_gas.open() as stream:
            rows = list(csv.DictReader(stream))
        candidate_gas = {(int(r['i']),int(r['j']),int(r['k'])):r for r in rows}
        cells = {(i,j,k) for i in range(1,args.cells+1) for j in range(1,args.cells+1) for k in range(1,args.cells+1)}
        if len(rows)!=args.cells**3 or len(reference_gas)!=args.cells**3 or set(candidate_gas)!=cells:
            raise ValueError('Expected one regular patch with the requested cell count')
        spacing_cm = .4 / args.cells
        mapped = {}
        for row in reference_gas:
            xyz = [float(row[c]) for c in ('x','y','z')]
            key = tuple(round(v/spacing_cm+.5) for v in xyz)
            if key not in cells or key in mapped or any(abs(v-(i-.5)*spacing_cm)>1e-12 for v,i in zip(xyz,key)):
                raise ValueError('Invalid reference cell coordinates')
            if (not math.isfinite(float(row['time'])) or
                abs(float(row['time'])-args.expected_time)>1e-20):
                raise ValueError('Unexpected reference final time')
            mapped[key] = row
        fields = {'rho':('density',1000),'rhou':('xmom',10),'rhov':('ymom',10),
                  'rhow':('zmom',10),'rhoE':('rho_E',.1),'temperature':('Temp',1),
                  **{f'rhoY_{s}':(f'rho_{s}',1000) for s in ('NC7H16','O2','N2')}}
        report['gas'] = {}
        for field,(source,factor) in fields.items():
            expected = [float(mapped[key][source])*factor for key in sorted(cells)]
            actual = [float(candidate_gas[key][field]) for key in sorted(cells)]
            if not all(map(math.isfinite,expected+actual)):
                raise ValueError('Nonfinite gas state')
            difference = max(abs(a-b) for a,b in zip(expected,actual))
            scale = max(map(abs,expected))
            report['gas'][field] = dict(max_absolute_difference=difference,reference_max=scale,
                                       relative_to_field_max=difference/scale if scale else None)
        for name,path in [('reference_gas',args.reference_gas),('candidate_gas',args.candidate_gas)]:
            report['inputs'][name] = dict(path=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(args.output)


if __name__ == '__main__':
    main()
