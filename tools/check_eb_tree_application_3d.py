#!/usr/bin/env python3
"""Focused multilevel EB application inventory and exact restart check."""
import argparse
import csv
import math
import re
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--automatic', action='store_true')
    args = parser.parse_args()
    exe, work = args.exe.resolve(), args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    template = args.input.read_text().replace("output_file='eb_tree.csv'", "output_file='result.csv'")

    if args.automatic:
        template = template.replace('maximum_steps=20', 'maximum_steps=200, automatic_timestep=.true., hydro_cfl=1.e-5')

    def run(name, extra='', text=None, success=True):
        directory = work / name
        directory.mkdir(exist_ok=True)
        source = template if text is None else text
        (directory / 'input.nml').write_text(source.rsplit('/', 1)[0] + extra + '\n/\n')
        result = subprocess.run([str(exe), 'input.nml'], cwd=directory, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
        (directory / 'run.log').write_text(result.stdout)
        if (result.returncode == 0) != success:
            raise AssertionError(f'{name}: exit {result.returncode}\n{result.stdout}')
        return directory

    full = run('full')
    if args.automatic:
        interval = float(re.search(r'Minimum accepted interval:\s+(\S+)', (full / 'run.log').read_text())[1])
        assert 0 < interval < 1e-9, interval
    stopped = run('stopped', "checkpoint_file='saved.chk', checkpoint_interval_steps=1, stop_after_checkpoint=.true.,")
    resumed = run('resumed', f"restart_file='{stopped / 'saved.chk'}',")
    for filename in ('result.csv', 'result.csv.budget.csv'):
        assert (full / filename).read_bytes() == (resumed / filename).read_bytes(), filename
    with (full / 'result.csv').open() as stream:
        cells = list(csv.DictReader(stream))
    assert {int(row['level']) for row in cells} == {0, 1, 2}
    assert all(float(row['fluid_volume']) > 0 and float(row['temperature']) > 0 for row in cells)
    fields = ['rho', 'rhou', 'rhov', 'rhow', 'rhoE'] + [key for key in cells[0] if key.startswith('rhoY_')]
    with (full / 'result.csv.budget.csv').open() as stream:
        budgets = list(csv.DictReader(stream))
    assert len(fields) == len(budgets)
    worst = 0.0
    for field, row in zip(fields, budgets):
        inventory = math.fsum(float(cell[field]) * float(cell['fluid_volume']) for cell in cells)
        final, initial = float(row['final']), float(row['initial'])
        assert math.isfinite(inventory) and abs(inventory-final) <= 3e-14*max(1, abs(final))
        exchange = sum(float(row[key]) for key in ('boundary', 'wall', 'chemistry'))
        residual = abs(final-initial-exchange)/max(1, abs(initial), abs(final), abs(exchange))
        assert residual <= 5e-10
        worst = max(worst, residual)
    # Reject both policy mismatch and a late corrupted record before output publication.
    corrupt = work / 'broken.chk'
    corrupt.write_text((stopped / 'saved.chk').read_text().replace('END_EB_TREE_CHECKPOINT', 'BROKEN_END'))
    for name, restart, source in (
        ('bad_policy', stopped / 'saved.chk', template.replace('time_step=1.e-9', 'time_step=2.e-9')),
        ('bad_record', corrupt, template),
        ('bad_mode', stopped / 'saved.chk', template.replace('automatic_timestep=.true.', 'automatic_timestep=.false.')
         if args.automatic else template.replace('maximum_steps=20', 'maximum_steps=20, automatic_timestep=.true.')),
    ):
        directory = work / name
        directory.mkdir(exist_ok=True)
        for filename in ('result.csv', 'result.csv.budget.csv', 'protected.chk'):
            (directory / filename).write_text('preserve\n')
        run(name, f"restart_file='{restart}', checkpoint_file='protected.chk', checkpoint_interval_steps=1,",
            source, success=False)
        for filename in ('result.csv', 'result.csv.budget.csv', 'protected.chk'):
            assert (directory / filename).read_text() == 'preserve\n'
    print(f'PASS: {len(cells)} exposed cells, exact restart CSV/budget, residual={worst:.17e}; rejection preserves files')


if __name__ == '__main__':
    main()
