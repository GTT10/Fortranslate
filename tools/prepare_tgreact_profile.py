#!/usr/bin/env python3
"""Prepare a separate, explicitly sampled TGReact cosine-profile diagnostic."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def cosine_row(x, length=.01):
    wave = math.cos(2 * math.pi * x / length)
    hydrogen = .06 + .002 * wave
    oxygen = .5
    nitrogen = 1 - hydrogen - oxygen
    return [x, hydrogen, oxygen, nitrogen, 2000 + 50 * wave, 101325.,
            hydrogen, oxygen, nitrogen, 0., 0., 0., 0., 0., 0.]


def prepare(source, output, points):
    """Keep the frozen profile untouched; verify its declared analytic basis."""
    if type(points) is not int or points not in (129, 257, 513):
        raise ValueError('Diagnostic profile requires 129, 257 or 513 points')
    source, output = Path(source).resolve(), Path(output).resolve()
    if source == output or (output.exists() and output.samefile(source)):
        raise ValueError('The diagnostic must not replace its frozen source profile')
    with source.open() as stream:
        rows = list(csv.reader(stream))
    expected_header = ['x', 'H2_init', 'O2_init', 'N2_init', 'Tad', 'p',
                       'H2', 'O2', 'N2', 'H2O', 'H', 'O', 'OH', 'HO2', 'H2O2']
    if not rows or rows[0] != expected_header or len(rows) != 34:
        raise ValueError('Expected the frozen 33-point TGReact profile')
    for i, row in enumerate(rows[1:]):
        values = list(map(float, row))
        expected = cosine_row(.01 * i / 32)
        if len(values) != 15 or any(not math.isclose(a, b, rel_tol=2.e-15, abs_tol=2.e-17)
                                    for a, b in zip(values, expected)):
            raise ValueError(f'Frozen row {i} differs from the declared cosine profile')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(expected_header)
        writer.writerows([[f'{v:.17e}' for v in cosine_row(.01 * i / (points - 1))]
                         for i in range(points)])
    return {'qualification': False, 'kind': 'separate_cosine_sampling_diagnostic',
            'points': points, 'domain_m': .01,
            'temperature_K': '2000 + 50*cos(2*pi*x/L)',
            'mass_fractions': 'H2=0.06+0.002*cos(2*pi*x/L), O2=0.5, N2=1-H2-O2',
            'pressure_Pa': 101325., 'trace_species': 'zero',
            'source': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
            'output': str(output), 'output_sha256': hashlib.sha256(output.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1] /
                        'cases/tgreact_eb_amr_reference/profile.csv')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--points', type=int, default=129, choices=(129, 257, 513))
    args = parser.parse_args()
    record = prepare(args.source, args.output, args.points)
    args.output.with_suffix('.manifest.json').write_text(json.dumps(record, indent=2) + '\n')
    print(f'Prepared {args.points} points; separate diagnostic, frozen source preserved.')


if __name__ == '__main__':
    main()
