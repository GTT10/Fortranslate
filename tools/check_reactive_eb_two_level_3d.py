#!/usr/bin/env python3
"""Check the fixed elementary two-level EB application's public artifacts."""
import argparse
import csv
import math
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    exe = args.exe.resolve()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    source = args.case.read_text()

    def run(name, content, success):
        path = work / name
        path.write_text(content)
        result = subprocess.run([str(exe), name], cwd=work, text=True, capture_output=True, timeout=90)
        (work / (name + '.log')).write_text(result.stdout + result.stderr)
        assert (result.returncode == 0) == success, result.stdout + result.stderr
        return result.stdout

    log = run('reference.nml', source, True)
    assert 'Completed steps: 3' in log, log
    output = work / 'eb_two_level.csv'
    budget_path = work / 'eb_two_level.csv.budget.csv'
    with output.open() as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        fields = reader.fieldnames[9:]
    assert len(rows) == 104, len(rows)
    totals = [0.0] * len(fields)
    initial_mass = 0.0
    levels = set()
    keys = set()
    for row in rows:
        assert all(math.isfinite(float(value)) for value in row.values())
        level = int(row['level'])
        ijk = tuple(int(row[c]) for c in ('i', 'j', 'k'))
        assert (level, ijk) not in keys
        keys.add((level, ijk))
        levels.add(level)
        assert level in (0, 1) and all(1 <= i <= 4 for i in ijk)
        if level == 0:
            assert not all(2 <= i <= 3 for i in ijk), 'refined parent duplicated in composite output'
        kappa = float(row['kappa'])
        assert 0 < kappa <= 1 and float(row['temperature']) > 0 and float(row['rho']) > 0
        volume = kappa / (4 if level == 0 else 8) ** 3
        initial_mass += volume * (0.93 if kappa < 1 - 1e-12 else 0.31)
        species_sum = sum(float(value) for key, value in row.items() if key.startswith('rhoY_'))
        assert math.isclose(species_sum, float(row['rho']), rel_tol=1e-10, abs_tol=1e-12)
        for c, name in enumerate(fields):
            totals[c] += volume * float(row[name])
    assert levels == {0, 1}
    with budget_path.open() as stream:
        budgets = list(csv.DictReader(stream))
    assert len(budgets) == len(fields)
    assert math.isclose(float(budgets[0]['initial']), initial_mass, rel_tol=2e-12, abs_tol=1e-13)
    for c, row in enumerate(budgets):
        values = {key: float(value) for key, value in row.items()}
        assert all(math.isfinite(value) for value in values.values())
        assert int(values['component']) == c + 1
        assert math.isclose(totals[c], values['final'], rel_tol=3e-12, abs_tol=3e-12)
        residual = values['final'] - values['initial'] - values['boundary'] - values['wall'] - values['chemistry']
        assert abs(residual) <= 5e-10 * max(1.0, abs(values['initial']), abs(values['final']))
        assert math.isclose(residual, values['residual'], rel_tol=1e-12, abs_tol=1e-15)
    saved = output.read_bytes(), budget_path.read_bytes()
    run('bad_ratio.nml', source.replace('refinement_ratio=2', 'refinement_ratio=3'), False)
    assert saved == (output.read_bytes(), budget_path.read_bytes())
    alias = source.replace("output_file='eb_two_level.csv'", "output_file='alias.nml'")
    run('alias.nml', alias, False)
    assert (work / 'alias.nml').read_text() == alias
    print('PASS: two-level EB application, composite output, inventories and rejected-input preservation')


if __name__ == '__main__':
    main()
