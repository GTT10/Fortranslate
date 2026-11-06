#!/usr/bin/env python3
"""Audit paired interior ratio-two shared-limiter stencils, without changing gates.

Input JSON: labels; samples[{name, reference, candidate, offset}]. Each stencil
contains seven vectors: center, x-, x+, y-, y+, z-, z+. Components 0..3 are
density and momentum in consistent SI units. Samples must already be matched
by stage/time/location; this tool does not infer correspondence. Remaining
components participate in the same scalar limiter. Optional sample context is
copied verbatim. Output is diagnostic, never a trajectory qualification.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def replay(stencil, offset, labels):
    count = len(labels)
    if len(stencil) != 7 or any(len(row) != count for row in stencil):
        raise ValueError('Expected center plus six equal-width neighbor vectors')
    if any(not math.isfinite(v) for row in stencil for v in row):
        raise ValueError('Nonfinite stencil')
    if any(row[0] <= 0 for row in stencil):
        raise ValueError('Nonpositive density')
    if len(offset) != 3 or any(v not in (-.25, .25) for v in offset):
        raise ValueError('Expected interior ratio-two child offsets +/-0.25')
    center = stencil[0]
    result = center.copy()
    axes = []
    for axis in range(3):
        low, high = stencil[1+2*axis:3+2*axis]
        components = []
        for name, left, middle, right in zip(labels, low, center, high):
            dc = .5*(right-left)
            dl, dr = 2*(middle-left), 2*(right-middle)
            bound = min(abs(dl), abs(dr)) if dl*dr >= 0 else 0.
            factor = min(abs(dc), bound)/abs(dc) if dc else 1.
            components.append(dict(name=name, centered_slope=dc,
                                   left_difference=middle-left,
                                   right_difference=right-middle, factor=factor))
        factor = min(x['factor'] for x in components)
        for index, component in enumerate(components):
            result[index] += offset[axis]*factor*component['centered_slope']
        axes.append(dict(axis='xyz'[axis], factor=factor,
                         controls=[x['name'] for x in components if x['factor'] == factor],
                         components=components))
    if result[0] <= 0 or any(not math.isfinite(v) for v in result):
        raise ValueError('Invalid interpolated state')
    return result, axes


def analyze(spec):
    labels = spec['labels']
    if len(labels) < 4 or len(set(labels)) != len(labels):
        raise ValueError('Need unique labels, with density and three momenta first')
    samples = spec['samples']
    if not samples or len({s['name'] for s in samples}) != len(samples):
        raise ValueError('Need nonempty, uniquely named samples')
    output = []
    for sample in samples:
        ref, candidate, offset = (sample[k] for k in ('reference', 'candidate', 'offset'))
        rv, ra = replay(ref, offset, labels)
        cv, ca = replay(candidate, offset, labels)
        velocity = lambda v: [v[i]/v[0] for i in (1, 2, 3)]
        reference_velocity = velocity(rv)
        delta = lambda v: [a-b for a, b in zip(velocity(v), reference_velocity)]
        # Reference-anchored identity: fc*sc-fr*sr = (fc-fr)*sr+fc*(sc-sr).
        # This ordering assigns the cross term to slope, not a causal share.
        center = [a-b for a, b in zip(candidate[0], ref[0])]
        terms = []
        for axis, (r_axis, c_axis) in enumerate(zip(ra, ca)):
            terms.append(dict(axis='xyz'[axis],
                factor_change=[offset[axis]*(c_axis['factor']-r_axis['factor'])*r['centered_slope']
                               for r in r_axis['components']],
                slope_change=[offset[axis]*c_axis['factor']*(c['centered_slope']-r['centered_slope'])
                              for r, c in zip(r_axis['components'], c_axis['components'])]))
        total = [math.fsum([center[i]]+[t[k][i] for t in terms
                             for k in ('factor_change', 'slope_change')])
                 for i in range(len(labels))]
        # du = dm/rho_c - u_r*drho/rho_c, avoiding cancellation of reciprocals.
        momentum_term = [total[i]/cv[0] for i in (1, 2, 3)]
        density_term = [-u*total[0]/cv[0] for u in reference_velocity]
        predicted = [a+b for a, b in zip(momentum_term, density_term)]
        decomposition = dict(
            convention='Reference slopes for factor changes; candidate factors for slope changes; cross term belongs to slope. Algebraic identity, not causal shares.',
            labels=labels, center_change=center, directions=terms,
            reconstructed_state_difference=total,
            state_closure_residual=[c-r-t for c, r, t in zip(cv, rv, total)],
            velocity_momentum_term_m_s=momentum_term,
            velocity_density_term_m_s=density_term,
            reconstructed_velocity_difference_m_s=predicted,
            velocity_closure_residual_m_s=[d-p for d, p in zip(delta(cv), predicted)])
        substitutions = []
        for index, label in enumerate(labels):
            mixed = [row.copy() for row in candidate]
            for pos in range(7):
                mixed[pos][index] = ref[pos][index]
            value, axes = replay(mixed, offset, labels)
            substitutions.append(dict(component=label, velocity_difference_m_s=delta(value),
                                      common_factors=[a['factor'] for a in axes]))
        output.append(dict(name=sample['name'], context=sample.get('context'),
                           reference_axes=ra, candidate_axes=ca,
                           velocity_difference_m_s=delta(cv),
                           difference_decomposition=decomposition,
                           maximum_donor_difference_by_component={
                               label: max(abs(a[index]-b[index]) for a, b in zip(candidate, ref))
                               for index, label in enumerate(labels)},
                           substitutions=substitutions))
    return dict(qualification=False, scope='Paired interior shared-limiter conditioning; '
                'single-component substitutions are diagnostic, not causal shares or CFD changes',
                samples=output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        parser.error('Output must not overwrite input')
    result = analyze(json.loads(args.input.read_text()))
    result['sha256'] = {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in (args.input, Path(__file__))}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(f"Audited {len(result['samples'])} paired stencils; qualification=false")


if __name__ == '__main__':
    main()
