#!/usr/bin/env python3
"""Independently reduce composite fluid-leaf pressure histories by level/geometry."""
import argparse
import csv
import json
import math
from pathlib import Path
from audit_eb_initial_interface import sha


def reduce_cells(path, history):
    with Path(path).open() as f:
        rows = list(csv.DictReader(f))
    clocks = {}
    for r in rows:
        if not all(math.isfinite(float(v)) for v in r.values()):
            raise ValueError('Nonfinite cell record')
        it = int(r['step'])
        key = tuple(int(r[k]) for k in ('node', 'i', 'j', 'k'))
        cells = clocks.setdefault(it, {})
        if key in cells:
            raise ValueError('Duplicate fluid leaf')
        if not 0 < float(r['vfrac']) <= 1 or float(r['volume_m3']) <= 0:
            raise ValueError('Invalid fluid volume')
        cells[key] = r
    if set(clocks) != set(range(len(history))):
        raise ValueError('Incomplete cell clocks')
    initial = clocks[0]
    result = []
    for it, h in enumerate(history):
        cells = clocks[it]
        if cells.keys() != initial.keys():
            raise ValueError('Changed leaf coverage')
        for key, r in cells.items():
            if any(r[k] != initial[key][k] for k in ('vfrac', 'volume_m3')):
                raise ValueError('Changed geometry')
        peak = max(cells.values(), key=lambda r: abs(float(r['pressure_offset_Pa'])))
        maximum = abs(float(peak['pressure_offset_Pa']))
        if not math.isclose(maximum, float(h['max_pressure_offset_Pa']), rel_tol=1e-13, abs_tol=1e-14):
            raise AssertionError('Cell/history pressure mismatch')
        groups = {}
        for r in cells.values():
            name = 'level' + str(int(r['node'])-1) + ('_cut' if float(r['vfrac']) < 1 else '_full')
            groups.setdefault(name, []).append(r)
        summary = {}
        for name, group in groups.items():
            vol = math.fsum(float(r['volume_m3']) for r in group)
            summary[name] = {'count': len(group), 'volume_m3': vol,
                'max_pressure_offset_Pa': max(abs(float(r['pressure_offset_Pa'])) for r in group),
                'volume_mean_abs_pressure_offset_Pa': math.fsum(float(r['volume_m3'])*abs(float(r['pressure_offset_Pa'])) for r in group)/vol}
        result.append({'step': it, 'peak': peak, 'groups': summary})
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prefix', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    report = {'qualification': False, 'cases': {}, 'sha256': {}}
    for policy in ('baseline', 'paired'):
        cells = Path(str(a.prefix)+f'_{policy}_evolution_cells.csv')
        hist = Path(str(a.prefix)+f'_{policy}_evolution.csv')
        with hist.open() as f:
            history = list(csv.DictReader(f))
        report['cases'][policy] = reduce_cells(cells, history)
        report['sha256'].update({str(p): sha(p) for p in (cells, hist)})
    report['status'] = 'passed'
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
