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
    args = parser.parse_args()
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
        scope='Two-step 400 K gas / 310 K liquid fixture; diagnostic only, no parity gate',
        inputs={name: dict(path=str(path.resolve()),
                           sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                for name, path in [('reference', args.reference), ('candidate', args.candidate)]},
        particles=particles)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(args.output)


if __name__ == '__main__':
    main()
