"""Check LES input, coupled budgets and changed-rank particle restart."""
import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('exe', 'mpiexec', 'source', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--recursive', action='store_true')
    parser.add_argument('--selected', action='store_true', help='Use selected H2O2 composition input')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='les_', dir=args.work.resolve()))
    env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')
    archive = work / 'stopped' / 'saved.step1'
    template = (args.source / 'cases/mpi_eb_spray_3d/tree.nml').read_text()
    template = template.replace("riemann_solver='hllc'", "riemann_solver='pelec',hydro_method='plm_ctu'")
    template = template.replace('chemistry_enabled=.true.', 'chemistry_enabled=.false.')
    if args.selected:
        template = template.replace('&eb_tree', "&eb_tree\n composition_count=1, "
                                    "composition_species(1)='N2', composition_mole_fractions(1)=1.,", 1)
    template = template.replace("output_file='result.csv'", "output_file='result.csv',time_integrator='sdc',"
                                "sdc_iterations=2,boundary_faces=6*'outflow'")
    if args.recursive:
        template = template.replace('&eb_tree', '&eb_tree\n subcycling=.true.,', 1)
    spray = """&eb_spray
 nc=1,number=2,vapor_names(1)='H2O',fractions(1)=1.,particle_outflow=.true.,
 total_mass=1.e-8,center=.0015,.0025,.0015,
 axis=1.,0.,0.,speed=100.,diameter=1.e-5,liquid_temperature=350.
/
"""
    les = "\n&eb_spray_les les_enabled=.true.,les_model='smagorinsky',les_cs=.17 /\n"

    def run(label, controls='', ranks=4, settings=les, cut=False, error=None):
        directory = work / label
        shutil.copytree(args.source / 'cases/mpi_eb_spray_3d', directory)
        if not cut:
            (directory / 'level_set.dat').write_text('8 8 8\n' + '1.0\n' * 9**3)
        (directory / 'tree.nml').write_text(template.replace('&eb_tree', '&eb_tree\n' + controls, 1))
        (directory / 'spray.nml').write_text(spray + settings)
        result = subprocess.run([str(args.mpiexec.resolve()), '-n', str(ranks), str(args.exe.resolve()),
                                 'tree.nml', 'spray.nml'], cwd=directory, env=env, capture_output=True,
                                text=True, timeout=180)
        log = result.stdout + result.stderr
        (directory / 'run.log').write_text(log)
        if error:
            assert result.returncode != 0 and error in log, log
        else:
            assert result.returncode == 0, log
        return directory

    def outputs(directory):
        particles = [row for path in sorted(directory.glob('result.csv.particles.*.csv')) for row in rows(path)]
        particles.sort(key=lambda row: int(row['id']))
        return ((directory / 'result.csv').read_bytes(),
                (directory / 'result.csv.budget.csv').read_bytes(), particles)

    disabled = run('disabled', settings='')
    continuous = run('continuous')
    baseline = outputs(continuous)
    assert baseline[0] != outputs(disabled)[0], 'LES must change the coupled gas CSV'
    budget = rows(continuous / 'result.csv.budget.csv')
    mass = max(abs(float(budget[0][key])) for key in ('initial', 'final', 'injected'))
    energy = max(abs(float(budget[4][key])) for key in ('initial', 'final', 'injected'))
    maximum = 0.0
    for row in budget:
        residual = float(row['final']) - float(row['initial'])
        residual -= sum(float(row[key]) for key in ('boundary', 'wall', 'chemistry', 'injected', 'force'))
        scale = max(abs(float(row[key])) for key in ('initial', 'final', 'injected'))
        if 2 <= int(row['component']) <= 4:
            scale = max(scale, math.sqrt(mass * energy))
        scale = max(scale, 1e-30)
        assert abs(residual) <= 3e-12 * scale, (row, residual)
        maximum = max(maximum, abs(residual) / scale)
    run('stopped', "checkpoint_file='saved',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,", 2)
    assert archive.exists(), 'Checkpoint was not produced'
    restart = f"restart_file='{archive}',"
    continued = run('continued', restart)
    assert outputs(continued) == baseline, 'Changed-rank gas/particle continuation differs'
    run('changed_les', restart, settings=les.replace('les_cs=.17', 'les_cs=.18'), error='EB spray context mismatch')
    run('removed_les', restart, settings='', error='EB spray context mismatch')
    run('cut_rejected', cut=True, error='Invalid particle LES state or non-regular geometry')
    print(f'PASS: LES changes gas; combined residual {maximum:.17g}; exact 2-to-4-rank restart; changed LES/cut cells rejected')
    print(work)


if __name__ == '__main__':
    main()
