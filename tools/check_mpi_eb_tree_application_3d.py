#!/usr/bin/env python3
"""Serial/MPI EB tree state parity, rank-changing restart and output protection."""
import argparse
import csv
import math
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--serial', type=Path, required=True)
    parser.add_argument('--mpi', type=Path, required=True)
    parser.add_argument('--mpiexec', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    template = args.input.read_text().replace("output_file='eb_tree.csv'", "output_file='result.csv'")
    env = os.environ.copy()
    env['OMPI_MCA_rmaps_base_oversubscribe'] = '1'

    def run(name, ranks=0, extra='', text=None, success=True):
        directory = work / name
        directory.mkdir(exist_ok=True)
        source = template if text is None else text
        (directory / 'input.nml').write_text(source.rsplit('/', 1)[0] + extra + '\n/\n')
        command = [str(args.serial.resolve()), 'input.nml']
        if ranks:
            command = [str(args.mpiexec.resolve()), '-n', str(ranks), str(args.mpi.resolve()), 'input.nml']
        result = subprocess.run(command, cwd=directory, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
        (directory / 'run.log').write_text(result.stdout)
        if (result.returncode == 0) != success:
            raise AssertionError(f'{name}: exit {result.returncode}\n{result.stdout}')
        return directory

    serial = run('serial')
    full = run('mpi2', 2)
    third = run('mpi3', 3)
    stop_policy = "checkpoint_file='saved.chk', checkpoint_interval_steps=1, stop_after_checkpoint=.true.,"
    stopped = run('stopped2', 2, stop_policy)
    resumed = run('resumed2', 2, f"restart_file='{stopped / 'saved.chk'}',")
    changed = run('resumed3', 3, f"restart_file='{stopped / 'saved.chk'}',")
    to_serial = run('mpi_to_serial', 0, f"restart_file='{stopped / 'saved.chk'}',")
    serial_stop = run('serial_stop', 0, stop_policy)
    from_serial = run('serial_to_mpi', 3, f"restart_file='{serial_stop / 'saved.chk'}',")
    expected = (serial / 'result.csv').read_bytes()
    for directory in (full, third, resumed, changed, to_serial, from_serial):
        assert (directory / 'result.csv').read_bytes() == expected, directory
    for filename in ('result.csv', 'result.csv.budget.csv'):
        assert (full / filename).read_bytes() == (resumed / filename).read_bytes(), filename

    reference_budget = list(csv.DictReader((serial / 'result.csv.budget.csv').open()))
    worst = 0.0
    for directory in (serial, full, third, resumed, changed, to_serial, from_serial):
        with (directory / 'result.csv').open() as stream:
            cells = list(csv.DictReader(stream))
        with (directory / 'result.csv.budget.csv').open() as stream:
            budget = list(csv.DictReader(stream))
        fields = ['rho', 'rhou', 'rhov', 'rhow', 'rhoE'] + [key for key in cells[0] if key.startswith('rhoY_')]
        assert len(fields) == len(budget) == len(reference_budget)
        for field, row, reference in zip(fields, budget, reference_budget):
            initial, final = float(row['initial']), float(row['final'])
            inventory = math.fsum(float(cell[field]) * float(cell['fluid_volume']) for cell in cells)
            assert math.isfinite(inventory) and abs(inventory-final) <= 3e-14*max(1, abs(final))
            exchange = sum(float(row[key]) for key in ('boundary', 'wall', 'chemistry'))
            residual = abs(final-initial-exchange)/max(1, abs(initial), abs(final), abs(exchange))
            assert residual <= 5e-10
            worst = max(worst, residual)
            for key in ('initial', 'final', 'boundary', 'wall', 'chemistry'):
                assert abs(float(row[key])-float(reference[key])) <= 3e-14*max(1, abs(initial)), (directory, key)

    broken = work / 'broken.chk'
    broken.write_text((stopped / 'saved.chk').read_text().replace('END_EB_TREE_CHECKPOINT', 'BROKEN_END'))
    for name, extra in (
        ('corrupt', f"restart_file='{broken}',"),
        ('wrong_policy', f"restart_file='{stopped / 'saved.chk'}',time_step=2.e-9,"),
        ('unsupported_auto', 'automatic_timestep=.true.,'),
    ):
        directory = work / name
        directory.mkdir(exist_ok=True)
        paths = [directory / filename for filename in ('result.csv', 'result.csv.budget.csv', 'protected.chk')]
        for path in paths:
            path.write_text('sentinel\n')
        run(name, 2, extra + "checkpoint_file='protected.chk',checkpoint_interval_steps=1,", success=False)
        assert all(path.read_text() == 'sentinel\n' for path in paths), name
    print(f'PASS: serial/2/3-rank exact state, same-rank exact budget, changed-rank and serial restart; residual={worst:.17e}')


if __name__ == '__main__':
    main()
