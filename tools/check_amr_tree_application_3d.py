#!/usr/bin/env python3
"""Focused executable checks for the serial regular 3D tree driver."""
import argparse
import csv
import math
import re
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--reference-executable', type=Path, required=True)
    parser.add_argument('--case-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = args.work_dir.resolve()
    case = args.case_dir.resolve()

    def run(executable, inputs, label, success=True):
        result = subprocess.run([str(executable.resolve()), *map(str, inputs)], cwd=work,
                                capture_output=True, text=True, timeout=180)
        (work / f'{label}.log').write_text(result.stdout + result.stderr)
        assert (result.returncode == 0) == success, f'{label}: exit {result.returncode}'

    def rows(path):
        with path.open() as stream:
            return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]

    run(args.executable, [case / 'hotspot.nml', case / 'branched.dat'], 'branched')
    data = rows(work / 'tree_hotspot.csv')
    assert len(data) == 127
    assert {int(row['level']) for row in data} == {0, 1, 2}
    assert all(row['time'] == 1e-9 and row['rho'] > 0 and row['temperature'] > 0 for row in data)
    volume = [row['dx'] * row['dy'] * row['dz'] for row in data]
    assert math.isclose(math.fsum(volume), 0.004**3, rel_tol=2e-14)
    for row in data:
        assert abs(math.fsum(v for k, v in row.items() if k.startswith('Y_')) - 1) < 1e-12
    budget = rows(work / 'tree_hotspot.csv.budget.csv')
    columns = ['rho', 'rhou', 'rhov', 'rhow', 'rhoE']
    columns += [k for k in data[0] if k.startswith('rhoY_')]
    for col, ledger in zip(columns, budget, strict=True):
        inventory = math.fsum(row[col] * v for row, v in zip(data, volume))
        assert math.isclose(inventory, ledger['final'], rel_tol=2e-12, abs_tol=1e-25), col
    for col, ledger in zip(columns[:5], budget[:5], strict=True):
        error = abs(ledger['final'] - ledger['initial'] - ledger['boundary_change'])
        # Use the absolute cell inventory when opposite momenta cancel globally.
        scale = max(abs(ledger['initial']), abs(ledger['final']), abs(ledger['boundary_change']),
                    math.fsum(abs(row[col]) * v for row, v in zip(data, volume)), 1e-20)
        assert error < 2e-10 * scale, ledger
    assert any(abs(row['final'] - row['initial'] - row['boundary_change']) > 1e-20 for row in budget[5:])

    # Root-only tree uses the existing regular-grid hydro implementation as oracle.
    config = (case / 'hotspot.nml').read_text().replace('chemistry_enabled=.true.', 'chemistry_enabled=.false.')
    config = config.replace('transport_enabled=.true.', 'transport_enabled=.false.')
    config = config.replace('wall_velocity(1,3)=2.0, wall_velocity(1,4)=2.0',
                            'wall_velocity(1,3)=0.0, wall_velocity(1,4)=0.0')
    config = config.replace('final_time=1.0e-9', 'final_time=1.0e-7')
    (work / 'root.dat').write_text('0\n')
    (work / 'root.nml').write_text(config.replace('tree_hotspot.csv', 'tree_root.csv'))
    (work / 'reference.nml').write_text(config.replace('tree_hotspot.csv', 'reference.csv'))
    run(args.executable, [work / 'root.nml', work / 'root.dat'], 'root')
    run(args.reference_executable, [work / 'reference.nml'], 'reference')
    steps = re.search(r'Steps:\s+(\d+)', (work / 'root.log').read_text())
    assert steps and int(steps[1]) >= 2
    root = rows(work / 'tree_root.csv')
    reference = rows(work / 'reference.csv')
    assert len(root) == len(reference) == 64
    for actual, expected in zip(root, reference):
        for key, value in expected.items():
            assert math.isclose(actual[key], value, rel_tol=2e-12, abs_tol=1e-12), (key, actual[key], value)
    saved = (work / 'tree_root.csv').read_bytes()
    (work / 'invalid.dat').write_text('1\n2 1 1 1 1 1 1 2\n')
    run(args.executable, [work / 'root.nml', work / 'invalid.dat'], 'invalid', success=False)
    assert (work / 'tree_root.csv').read_bytes() == saved
    print('PASS: branched reacting/transport budgets, leaf CSV, root oracle, invalid topology')


if __name__ == '__main__':
    main()
