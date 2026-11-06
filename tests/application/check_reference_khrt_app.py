"""Reference KHRT children: coupled inventory and changed-rank continuation."""
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
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='khrt_', dir=args.work.resolve()))
    env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')
    template = (args.source / 'cases/mpi_eb_spray_3d/tree.nml').read_text()
    template = template.replace("riemann_solver='hllc'", "riemann_solver='pelec',hydro_method='plm_ctu'")
    template = template.replace('chemistry_enabled=.true.', 'chemistry_enabled=.false.')
    template = template.replace('time_step=1.e-8,final_time=2.e-8', 'time_step=1.e-7,final_time=4.e-7')
    template = template.replace("output_file='result.csv'", "output_file='result.csv',time_integrator='sdc',"
                                "sdc_iterations=2,boundary_faces=6*'outflow'")
    template = template.replace('&eb_tree', '&eb_tree\n regular_density=1.,cut_density=1.,'
                                'initial_pressure=435483.87096774194,', 1)
    if args.recursive:
        template = template.replace('&eb_tree', '&eb_tree\n subcycling=.true.,', 1)
    spray = """&eb_spray
 nc=1,number=2,vapor_names(1)='H2O',fractions(1)=1.,particle_outflow=.true.,
 total_mass=1.e-10,center=.0015,.0025,.0015,
 axis=1.,0.,0.,speed=500.,diameter=1.e-5,liquid_temperature=350.,
 particle_options%breakup_model='khrt',particle_options%khrt_children=.true.,
 particle_options%khrt_reference_diameter_cubed=1.e-25,
 particle_options%surface_tension=.072,particle_options%liquid_viscosity=.001,
 particle_options%drag=.false.,particle_options%heat=.false.,particle_options%evaporation=.false.,
 particle_options%khrt_parcel_factor=.8
/
"""

    def run(label, controls='', ranks=4, factor='.8', error=None):
        directory = work / label
        shutil.copytree(args.source / 'cases/mpi_eb_spray_3d', directory)
        (directory / 'tree.nml').write_text(template.replace('&eb_tree', '&eb_tree\n' + controls, 1))
        (directory / 'spray.nml').write_text(spray.replace('khrt_parcel_factor=.8', 'khrt_parcel_factor=' + factor))
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

    legacy = run('legacy', factor='-1')
    continuous = run('continuous')
    baseline = outputs(continuous)
    assert len(baseline[2]) > 2 and any(int(row['id']) > 2 for row in baseline[2])
    assert baseline[2] != outputs(legacy)[2], 'Reference policy must change the child population'
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
    stopped = run('stopped', "checkpoint_file='saved',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,", 2)
    archive = stopped / 'saved.step1'
    continued = run('continued', f"restart_file='{archive}',")
    assert outputs(continued) == baseline, 'Gas/particle continuation differs after rank change'
    for name, factor in [('changed_factor', '.9'), ('removed_factor', '-1')]:
        run(name, f"restart_file='{archive}',", factor=factor, error='EB spray context mismatch')
    print(f'PASS: {len(baseline[2])} parcels; combined residual {maximum:.17g}; exact 2-to-4 restart; changed policy rejected')
    print(work)


if __name__ == '__main__':
    main()
