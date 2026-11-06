#!/usr/bin/env python3
"""Replay the experimental vector-scaled inflection screen on QBOUND stencils.

Only applies to the isolated candidate with momentum_roundoff_scale and
inflection_hierarchy enabled. This audits decisions, not physical accuracy.
Component and axis numbers in the report are one-based.
"""
import argparse
import collections
import hashlib
import json
import math
import sys
from pathlib import Path


def read_blocks(path):
    blocks = []
    current = None
    expected = {(i, j, k) for i in range(-2, 3)
                for j in range(-2, 3) for k in range(-2, 3)}
    for line in path.read_text().splitlines():
        fields = [value.strip() for value in line.split(',')]
        if fields[0] == 'QBOUND_SCREEN':
            if current is not None or len(fields) != 10 or fields[7] not in ('T', 'F'):
                raise ValueError('Malformed or nested QBOUND_SCREEN')
            current = {'source': list(map(int, fields[1:4])),
                       'smooth': fields[7] == 'T', 'samples': {}}
        elif fields[0] == 'QBOUND_SAMPLE':
            if current is None:
                raise ValueError('Sample outside block')
            index = tuple(map(int, fields[1:4]))
            values = list(map(float, fields[4:]))
            if index not in expected or index in current['samples']:
                raise ValueError('Invalid or duplicate stencil index')
            if len(values) < 5 or not all(map(math.isfinite, values)):
                raise ValueError('Invalid stencil components')
            current['samples'][index] = values
        elif fields[0] == 'QBOUND_END':
            if current is None or set(current['samples']) != expected:
                raise ValueError('Incomplete stencil')
            if len({len(v) for v in current['samples'].values()}) != 1:
                raise ValueError('Inconsistent component count')
            blocks.append(current)
            current = None
    if current is not None or not blocks:
        raise ValueError('Missing or unterminated stencil blocks')
    return blocks


def failures(q):
 for c in range(len(q[(0,0,0)])):
  for axis in range(3):
   for j in range(-2,3):
    for i in range(-2,3):
     indices=[(v,i,j) if axis==0 else (i,v,j) if axis==1 else (i,j,v) for v in range(-2,3)]
     line=[q[idx][c] for idx in indices];scale=max(map(abs,line))
     if c in [1,2,3]:scale=max(scale,max(math.sqrt(sum(x*x for x in q[idx][1:4])) for idx in indices))
     noise=16*sys.float_info.epsilon*scale;d1=[line[k+1]-line[k] for k in range(4)];d2=[d1[k+1]-d1[k] for k in range(3)]
     if all(abs(x)<=noise for x in d2) or all(x>noise for x in d2) or all(x<-noise for x in d2):continue
     reason=None
     if not (all(x>noise for x in d1) or all(x<-noise for x in d1)):reason='not_strictly_monotone'
     elif max(abs(d2[k+1]-d2[k]) for k in range(2))>min(map(abs,d1))+noise:reason='third_difference_exceeds_first'
     else:
      low=sum(w*x for w,x in zip([-3,22,128,-22,3],line))/128;high=2*line[2]-low
      if not(min(line[1:3])-noise<=low<=max(line[1:3])+noise and min(line[2:4])-noise<=high<=max(line[2:4])+noise):reason='child_outside_adjacent_bounds'
     if reason:yield {'component':c+1,'axis':axis+1,'line_indices':[i,j],'reason':reason,'line':line,'d2':d2,'noise':noise}


def first_failure(q):
    return next(failures(q), None)


def analyze(path, all_failures=False):
    counts = collections.Counter()
    results = []
    previous = {}
    transitions = []
    for block in read_blocks(path):
        key = tuple(block['source'])
        counts[key] += 1
        failure = first_failure(block['samples'])
        if (failure is None) != block['smooth']:
            raise ValueError(f'Runtime decision mismatch: {key}, fill {counts[key]}')
        result = {'source': list(key), 'fill_occurrence': counts[key],
                  'runtime_smooth': block['smooth'], 'first_failure': failure}
        if all_failures:
            rejected = list(failures(block['samples']))
            result['all_failures'] = rejected
            result['rejection_components'] = sorted({f['component'] for f in rejected})
        results.append(result)
        if key in previous and previous[key] != block['smooth']:
            transitions.append(result)
        previous[key] = block['smooth']
    fills = []
    for fill in sorted({r['fill_occurrence'] for r in results}):
        selected = [r for r in results if r['fill_occurrence'] == fill]
        reasons = collections.Counter(
            (r['first_failure']['component'], r['first_failure']['axis'],
             r['first_failure']['reason']) for r in selected if r['first_failure'])
        fills.append({'fill_occurrence': fill, 'sampled': len(selected),
                      'accepted': sum(r['runtime_smooth'] for r in selected),
                      'rejections': [{'component': c, 'axis': a, 'reason': reason,
                                      'count': count}
                                     for (c, a, reason), count in sorted(reasons.items())]})
    return {'status': 'passed', 'qualification': False,
            'policy': 'vector_roundoff_and_inflection_hierarchy',
            'input_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'exact_replayed_decisions': len(results), 'fills': fills,
            'transitions': transitions, 'results': results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--all-failures', action='store_true',
                        help='Record every rejecting line, not only the first')
    args = parser.parse_args()
    result = analyze(args.input, args.all_failures)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('results', 'transitions')}, indent=2))


if __name__ == '__main__':
    main()
