#!/usr/bin/env python3
"""Reduce pre/post gas+held reflux pressure changes on identical fluid leaves.

These are intermediate Euler states, not accepted SSPRK2 clocks or budgets.
The pre-sync state already includes any coarse/fine spatial halo influence.
"""
import argparse
import csv
import json
import math
from pathlib import Path
from audit_eb_initial_interface import sha

KEY = ('node', 'i', 'j', 'k')
PRE = 'pre_sync_pressure_offset_Pa'
POST = 'post_sync_pressure_offset_Pa'


def summarize(rows, expected, stages):
    grouped = {}
    for r in rows:
        if not all(math.isfinite(float(v)) for v in r.values()):
            raise ValueError('Nonfinite stage data')
        stage = int(r['stage'])
        key = tuple(int(r[k]) for k in KEY)
        group = grouped.setdefault(stage, {})
        if key in group:
            raise ValueError('Duplicate stage cell')
        group[key] = r
    if set(grouped) != set(range(1, stages+1)):
        raise ValueError('Incomplete stage clocks')
    result = []
    for stage, cells in sorted(grouped.items()):
        if set(cells) != expected:
            raise ValueError('Stage leaf coverage mismatch')
        levels = {}
        for node in sorted({k[0] for k in expected}):
            subset = {k: r for k, r in cells.items() if k[0] == node}
            peak = max(subset, key=lambda k: abs(float(subset[k][POST])-float(subset[k][PRE])))
            r = subset[peak]
            levels[str(node)] = {'count': len(subset),
                'max_pre_pressure_offset_Pa': max(abs(float(v[PRE])) for v in subset.values()),
                'max_post_pressure_offset_Pa': max(abs(float(v[POST])) for v in subset.values()),
                'max_sync_pressure_change_Pa': abs(float(r[POST])-float(r[PRE])),
                'peak_change_cell': peak, 'peak_pre_Pa': float(r[PRE]), 'peak_post_Pa': float(r[POST])}
        result.append({'stage': stage, 'nodes': levels})
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prefix', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    report = {'qualification': False, 'cases': {}, 'sha256': {}}
    for policy in ('baseline', 'paired'):
        base = str(a.prefix)+'_'+policy
        files = [Path(base+s) for s in ('_evolution_cells.csv', '_evolution.csv', '_stages.csv')]
        with files[0].open() as f:
            expected = {tuple(int(r[k]) for k in KEY) for r in csv.DictReader(f) if int(r['step']) == 0}
        with files[1].open() as f:
            stages = 2*(len(list(csv.DictReader(f)))-1)
        with files[2].open() as f:
            report['cases'][policy] = summarize(list(csv.DictReader(f)), expected, stages)
        report['sha256'].update({str(f): sha(f) for f in files})
    report['status'] = 'passed'
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
