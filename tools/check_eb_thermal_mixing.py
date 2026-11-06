#!/usr/bin/env python3
"""Check local EOS mixing controls; no AMR trajectory qualification."""
import argparse
import csv
import json
import math
from pathlib import Path
from audit_eb_initial_interface import sha


def summarize(rows):
    groups = {}
    for row in rows:
        if not all(math.isfinite(float(v)) for v in row.values()):
            raise ValueError('Nonfinite mixing record')
        variant = int(row['variant'])
        if variant not in (0, 1, 2):
            raise ValueError('Unknown mixing variant')
        key = (variant, *(int(row[k]) for k in ('i', 'j', 'k')))
        groups.setdefault(key, []).append(row)
        a = float(row['fraction'])
        if not 0 < a < 1:
            raise ValueError('Invalid mixing fraction')
        left, right = (float(row[k]) for k in ('p_left_Pa', 'p_right_Pa'))
        target = float(row['linear_pressure_Pa'])
        if abs(target-((1-a)*right+a*left)) > 1e-8 or abs(left-right) > 1e-8:
            raise AssertionError('Donor pressure equilibrium mismatch')
        if variant == 2 and abs(float(row['temperature_left_K'])-float(row['temperature_right_K'])) > 1e-8:
            raise AssertionError('Temperature control mismatch')
    if not groups:
        raise ValueError('Empty mixing controls')
    coverage = {v: {key[1:] for key in groups if key[0] == v} for v in range(3)}
    if not coverage[0] or not coverage[0] == coverage[1] == coverage[2]:
        raise ValueError('Incomplete variant coverage')
    result = {}
    for key, values in groups.items():
        values.sort(key=lambda r: float(r['fraction']), reverse=True)
        fractions = [float(r['fraction']) for r in values]
        if len(values) != 3 or any(not math.isclose(fractions[0]/2**i, a, rel_tol=1e-13) for i, a in enumerate(fractions)):
            raise ValueError('Incomplete fraction sequence')
    for variant in range(3):
        levels = []
        for level in range(3):
            subset = [v[level] for k, v in groups.items() if k[0] == variant]
            peak = max(subset, key=lambda r: abs(float(r['p_mixed_Pa'])-float(r['linear_pressure_Pa'])))
            error = abs(float(peak['p_mixed_Pa'])-float(peak['linear_pressure_Pa']))
            levels.append({'fraction_scale': 2**-level, 'max_pressure_defect_Pa': error,
                           'max_energy_gap_J_m3': max(abs(float(r['pressure_preserving_energy_gap_J_m3'])) for r in subset),
                           'peak': peak})
        if variant == 2 and max(v['max_pressure_defect_Pa'] for v in levels) > 1e-8:
            raise AssertionError('Isothermal mixture equilibrium failed')
        result[str(variant)] = levels
    return {'status': 'passed', 'qualification': False, 'pairs_per_variant': len(coverage[0]), 'variants': result}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    with a.input.open() as f:
        report = summarize(list(csv.DictReader(f)))
    report['input_sha256'] = sha(a.input)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
