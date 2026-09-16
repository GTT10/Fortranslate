#!/usr/bin/env python3
"""Scheduled tree topology changes: executable-level conservation and timing."""
import argparse
import csv
import math
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--case-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    case = args.case_dir.resolve()
    binary = args.executable.resolve()

    def run(schedule, label, success=True):
        command = [str(binary), str(case / 'hotspot.nml'), str(case / 'branched.dat')]
        if schedule:
            command.append(str(schedule))
        result = subprocess.run(command, cwd=work, text=True, capture_output=True, timeout=180)
        log = result.stdout + result.stderr
        (work / f'{label}.log').write_text(log)
        assert (result.returncode == 0) == success, log
        return log

    def rows(path):
        with path.open() as stream:
            return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]

    log = run(case / 'regrid.dat', 'dynamic')
    times = [float(t) for t in re.findall(r'Regrid time:\s+(\S+)', log)]
    assert times == [2.5e-10, 5e-10, 1e-9], times
    assert re.search(r'Regrids:\s+3\b', log)
    assert re.search(r'Regrid time:.*Nodes: 1 Levels: 1', log)
    data = rows(work / 'tree_hotspot.csv')
    assert len(data) == 127
    assert {row['level'] for row in data} == {0, 1, 2}
    assert all(row['time'] == 1e-9 and row['temperature'] > 0 and row['rho'] > 0 for row in data)
    volume = [row['dx'] * row['dy'] * row['dz'] for row in data]
    assert math.isclose(math.fsum(volume), .004**3, rel_tol=2e-14)
    ledger = rows(work / 'tree_hotspot.csv.budget.csv')
    columns = ['rho', 'rhou', 'rhov', 'rhow', 'rhoE'] + [k for k in data[0] if k.startswith('rhoY_')]
    for col, entry in zip(columns, ledger, strict=True):
        integrated = math.fsum(row[col] * v for row, v in zip(data, volume))
        assert math.isclose(integrated, entry['final'], rel_tol=2e-12, abs_tol=1e-25), col
    for col, entry in zip(columns[:5], ledger[:5], strict=True):
        scale = max(abs(entry['initial']), abs(entry['final']), abs(entry['boundary_change']),
                    math.fsum(abs(row[col]) * v for row, v in zip(data, volume)), 1e-20)
        assert abs(entry['final'] - entry['initial'] - entry['boundary_change']) < 2e-10 * scale, col
    for row in data:
        assert abs(math.fsum(v for k, v in row.items() if k.startswith('Y_')) - 1) < 1e-12
    # Invalid schedule is diagnosed before evolving or overwriting final files.
    saved = (work / 'tree_hotspot.csv').read_bytes()
    saved_budget = (work / 'tree_hotspot.csv.budget.csv').read_bytes()
    invalid = work / 'invalid.dat'
    invalid.write_text(f"2\n5e-10 '{case / 'root.dat'}'\n2e-10 '{case / 'branched.dat'}'\n")
    run(invalid, 'invalid', success=False)
    assert (work / 'tree_hotspot.csv').read_bytes() == saved
    assert (work / 'tree_hotspot.csv.budget.csv').read_bytes() == saved_budget
    # An identity event at final_time changes neither solution nor boundary ledger.
    run(None, 'static')
    reference = rows(work / 'tree_hotspot.csv')
    reference_budget = (work / 'tree_hotspot.csv.budget.csv').read_bytes()
    identity = work / 'identity.dat'
    identity.write_text(f"1\n1e-9 '{case / 'branched.dat'}'\n")
    run(identity, 'identity')
    actual = rows(work / 'tree_hotspot.csv')
    for a, b in zip(actual, reference, strict=True):
        for key in b:
            assert math.isclose(a[key], b[key], rel_tol=2e-12, abs_tol=1e-12), key
    assert (work / 'tree_hotspot.csv.budget.csv').read_bytes() == reference_budget
    print('PASS: scheduled refinement/coarsening, exact event times, composite budgets, identity, rejection')


if __name__ == '__main__':
    main()
