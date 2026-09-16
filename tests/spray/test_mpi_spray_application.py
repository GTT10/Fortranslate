#!/usr/bin/env python3
"""Focused input-to-output check for the fixed MPI spray executable."""
import argparse
import csv
import math
from pathlib import Path
import subprocess
import tempfile


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith('#')))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--mpiexec', type=Path, required=True)
    parser.add_argument('--ranks', type=int, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='run-', dir=args.work_dir.resolve()))
    command = [str(args.mpiexec.resolve()), '-n', str(args.ranks),
               str(args.executable.resolve()), str(args.input.resolve())]
    result = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=90)
    (work / 'run.log').write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    history = rows(work / 'periodic_mixture.history.csv')
    assert len(history) == 3 and int(history[-1]['step']) == 2
    assert math.isclose(float(history[-1]['time']), 2e-8, rel_tol=1e-14)
    for row in history:
        assert all(math.isfinite(float(value)) for value in row.values())
        assert abs(float(row['relative_mass_error'])) < 2e-11
        assert abs(float(row['relative_energy_error'])) < 2e-11
        assert int(row['parcels']) == 4
    gas = [row for path in sorted(work.glob('*.gas.csv')) for row in rows(path)]
    assert len(gas) == 8 * 4 * 4
    cells = {(int(row['i']), int(row['j']), int(row['k'])) for row in gas}
    assert cells == {(i, j, k) for i in range(1, 9) for j in range(1, 5) for k in range(1, 5)}
    for row in gas:
        assert all(math.isfinite(float(value)) for value in row.values())
        rho = float(row['rho'])
        species = [float(value) for name, value in row.items() if name.startswith('rhoY_')]
        assert rho > 0 and float(row['temperature']) > 0 and min(species) >= 0
        assert math.isclose(sum(species), rho, rel_tol=2e-12)
    particles = [row for path in sorted(work.glob('*.particles.csv')) for row in rows(path)]
    assert sorted(int(row['id']) for row in particles) == [1, 2, 3, 4]
    for row in particles:
        assert all(math.isfinite(float(value)) for value in row.values())
        assert math.isclose(float(row['mass']), float(row['component_mass_1']) + float(row['component_mass_2']), rel_tol=2e-12)
        assert float(row['temperature']) > 330
    # A second invocation must fail before replacing the saved history.
    original = (work / 'periodic_mixture.history.csv').read_bytes()
    repeated = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=90)
    (work / 'repeat.log').write_text(repeated.stdout + repeated.stderr)
    assert repeated.returncode != 0
    assert (work / 'periodic_mixture.history.csv').read_bytes() == original
    print(f'PASS MPI spray application: {args.ranks} ranks; artifacts: {work}')


if __name__ == '__main__':
    main()
