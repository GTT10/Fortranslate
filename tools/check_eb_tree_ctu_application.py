#!/usr/bin/env python3
"""Focused serial/MPI CTU tree regrid and changed-rank continuation check."""
import argparse
import csv
import math
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--serial', type=Path, required=True)
    parser.add_argument('--mpi', type=Path, required=True)
    parser.add_argument('--mpiexec', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--selected', action='store_true')
    parser.add_argument('--uniform-density', action='store_true')
    boundaries = parser.add_mutually_exclusive_group()
    boundaries.add_argument('--inflow', action='store_true')
    boundaries.add_argument('--periodic', action='store_true')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='run-', dir=args.work.resolve()))
    template = args.input.read_text().replace("riemann_solver='hllc'", "riemann_solver='pelec', hydro_method='plm_ctu'")
    template = template.replace('&eb_tree', '&eb_tree\n subcycling=.true.,', 1)
    template = template.replace("output_file='eb_tree.csv'", "output_file='result.csv'")
    if args.selected:
        template = template.replace("chemistry_integrator='explicit'", "chemistry_integrator='implicit'")
        template = template.replace('&eb_tree', "&eb_tree\n composition_count=3, "
                                    "composition_species='H2','O2','N2', "
                                    "composition_mole_fractions=.29572,.14785,.55643,", 1)
    if args.uniform_density:
        template = template.replace('&eb_tree', '&eb_tree\n cut_density=.31,', 1)
    if args.inflow:
        template = template.replace('&eb_tree', "&eb_tree\n "
            "boundary_faces='inflow','outflow','outflow','outflow','outflow','outflow', "
            "inflow_pressure(1)=135000, inflow_velocity(1,1)=100,", 1)
    if args.periodic:
        template = template.replace('&eb_tree', "&eb_tree\n "
            "boundary_faces='periodic','periodic','periodic','periodic','outflow','outflow',", 1)
    env = os.environ.copy()
    env.pop('DISPLAY', None)
    env['OMPI_MCA_rmaps_base_oversubscribe'] = '1'

    def run(name, ranks=0, extra='', text=None, success=True):
        directory = work / name
        directory.mkdir()
        source = template if text is None else text
        (directory / 'input.nml').write_text(source.rsplit('/', 1)[0] + extra + '\n/\n')
        command = [str(args.serial.resolve()), 'input.nml']
        if ranks:
            command = [str(args.mpiexec.resolve()), '-n', str(ranks), str(args.mpi.resolve()), 'input.nml']
        result = subprocess.run(command, cwd=directory, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
        (directory / 'run.log').write_text(result.stdout)
        if (result.returncode == 0) != success:
            raise AssertionError(f'{directory}: exit {result.returncode}\n{result.stdout}')
        return directory

    serial = run('serial')
    full = run('mpi2', 2)
    stopped = run('stopped2', 2, "checkpoint_file='saved.chk', checkpoint_interval_steps=2, stop_after_checkpoint=.true.,")
    restart = f"restart_file='{stopped / 'saved.chk'}',"
    resumed = run('resumed3', 3, restart)
    expected = (serial / 'result.csv').read_bytes()
    assert (full / 'result.csv').read_bytes() == expected
    assert (resumed / 'result.csv').read_bytes() == expected
    if not args.selected or args.uniform_density:
        assert 'nodes: 1' in (stopped / 'run.log').read_text()
    worst = 0.0
    for directory in (serial, full, resumed):
        with (directory / 'result.csv').open() as stream:
            cells = list(csv.DictReader(stream))
        assert {int(row['level']) for row in cells} == {0, 1, 2}
        assert all(math.isfinite(float(row['temperature'])) and float(row['temperature']) > 0 for row in cells)
        with (directory / 'result.csv.budget.csv').open() as stream:
            budget = list(csv.DictReader(stream))
        fields = ['rho', 'rhou', 'rhov', 'rhow', 'rhoE'] + [key for key in cells[0] if key.startswith('rhoY_')]
        assert len(fields) == len(budget)
        for field, row in zip(fields, budget):
            initial, final = float(row['initial']), float(row['final'])
            inventory = math.fsum(float(cell[field]) * float(cell['fluid_volume']) for cell in cells)
            assert math.isfinite(inventory) and abs(inventory-final) <= 3e-14*max(1, abs(final))
            exchange = sum(float(row[key]) for key in ('boundary', 'wall', 'chemistry'))
            residual = abs(final-initial-exchange)/max(1, abs(initial), abs(final), abs(exchange))
            assert residual <= 5e-10
            worst = max(worst, residual)
    run('wrong_method', 3, restart, template.replace("hydro_method='plm_ctu'", "hydro_method='pcm'"), False)
    run('wrong_solver', text=template.replace("riemann_solver='pelec'", "riemann_solver='hllc'"), success=False)
    run('unknown_method', text=template.replace("hydro_method='plm_ctu'", "hydro_method='unknown'"), success=False)
    rejections = 3
    if args.inflow:
        run('changed_inlet', 3, restart, template.replace('inflow_velocity(1,1)=100', 'inflow_velocity(1,1)=101'), False)
        run('changed_boundary', 3, restart, template.replace("'inflow'", "'outflow'"), False)
        run('conflicting_mode', text=template.replace('&eb_tree', '&eb_tree\n inflow_mach(1)=.2,', 1), success=False)
        rejections += 3
    if args.periodic:
        run('changed_boundary', 3, restart, template.replace("'periodic'", "'outflow'"), False)
        run('unpaired_boundary', text=template.replace("'periodic'", "'outflow'", 1), success=False)
        run('incompatible_geometry', text=template.replace("'outflow'", "'periodic'"), success=False)
        rejections += 3
    print(f'PASS: 4 successful runs, {rejections} expected rejections; maximum residual={worst:.6e}; artifacts={work}')


if __name__ == '__main__':
    main()
