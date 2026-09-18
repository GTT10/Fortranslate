"""MPI EB deposition/splash, film evolution and rank-changing application restart."""
import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--exe', type=Path, required=True)
p.add_argument('--mpiexec', type=Path, required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
p.add_argument('--selected', action='store_true')
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)
work = Path(tempfile.mkdtemp(prefix='impact_', dir=a.work.resolve()))
env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def population(directory, kind):
    return sorted((r for f in directory.glob(f'result.csv.{kind}.rank*.csv') for r in rows(f)),
                  key=lambda r: tuple(float(v) for v in r.values()))


def compare(left, right):
    assert len(left) == len(right), (len(left), len(right))
    for first, second in zip(left, right):
        assert first.keys() == second.keys()
        for key in first:
            x, y = float(first[key]), float(second[key])
            assert math.isfinite(x) and math.isfinite(y)
            assert abs(x-y) <= 3e-12 * max(abs(x), abs(y), 1e-30), (key, x, y)


def run(mode, label, ranks, extra='', success=True, altered=False):
    directory = work / f'{mode}_{label}'
    shutil.copytree(a.source / 'cases/mpi_eb_spray_3d', directory)
    tree = (directory / 'tree.nml').read_text().replace('time_step=1.e-8', 'time_step=5.e-9')
    controls = ('regrid_count=3,regrid_times=5.e-9,1.e-8,1.5e-8,regrid_child_count=1,0,1,'
                'regrid_lower(:,1,1)=3,1,1,regrid_upper(:,1,1)=4,4,4,'
                'regrid_lower(:,1,3)=1,1,1,regrid_upper(:,1,3)=2,4,4,'
                'rebalance_interval_steps=1,rebalance_minimum_gain=0.,')
    if a.selected:
        controls += "composition_count=1,composition_species(1)='N2',composition_mole_fractions(1)=1.,"
    (directory / 'tree.nml').write_text(tree.replace('&eb_tree', '&eb_tree\n' + controls + extra, 1))
    spray = (directory / 'spray.nml').read_text()
    spray = spray.replace('&eb_spray', '&eb_spray\n wall_impacts=.true.,'
                          'walls%enabled=6*.true.,walls%temperature=6*400.,walls%evolve=.true.,'
                          'walls%liquid_conductivity=.6,particle_options%liquid_viscosity=.001,'
                          'particle_options%surface_tension=.072,', 1)
    if mode == 'splash':
        spray = spray.replace('speed=20.', 'speed=200.').replace('total_mass=1.e-10', 'total_mass=1.e-16')
    if altered:
        spray = spray.replace('6*400.', '6*401.')
    (directory / 'spray.nml').write_text(spray)
    result = subprocess.run([str(a.mpiexec.resolve()), '-n', str(ranks), str(a.exe.resolve()), 'tree.nml', 'spray.nml'],
                            cwd=directory, env=env, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=60)
    (directory / 'run.log').write_text(result.stdout)
    assert (result.returncode == 0) == success, result.stdout
    return directory


for mode in ('deposit', 'splash'):
    full = run(mode, 'full', 4)
    films = population(full, 'films')
    particles = population(full, 'particles')
    assert len(films) == 2 and len(particles) == (0 if mode == 'deposit' else 8)
    for film in films:
        assert float(film['mass']) > 0 and float(film['temperature']) > 350
        assert float(film['ny']) == 1 and float(film['y']) == .00125
        assert float(film['component_mass_1']) == float(film['mass'])
    if particles:
        assert sorted(int(r['id']) for r in particles) == list(range(3, 11))
        assert all(float(r['v']) > 0 and float(r['y']) > .00125 for r in particles)
    log = (full / 'run.log').read_text()
    assert all(f'Regrid stage: {n}' in log for n in (1, 2, 3)), log
    stopped = run(mode, 'stopped', 2, "checkpoint_file='saved',checkpoint_interval_steps=2,stop_after_checkpoint=.true.,")
    archive = stopped / 'saved.step2'
    assert Path(str(archive) + '.particles').read_text().splitlines()[0] == 'PELEF_EB_PARTICLES_V2'
    resumed = run(mode, 'resumed', 4, f"restart_file='{archive}',")
    compare(rows(full / 'result.csv'), rows(resumed / 'result.csv'))
    compare(particles, population(resumed, 'particles'))
    compare(films, population(resumed, 'films'))
    for x, y in zip(rows(full / 'result.csv.budget.csv'), rows(resumed / 'result.csv.budget.csv')):
        compare([{k: x[k] for k in ('initial', 'final', 'boundary', 'wall', 'chemistry', 'injected')}],
                [{k: y[k] for k in ('initial', 'final', 'boundary', 'wall', 'chemistry', 'injected')}])
    mismatch = run(mode, 'mismatch', 2, f"restart_file='{archive}',", success=False, altered=True)
    assert 'EB spray context mismatch' in (mismatch / 'run.log').read_text()
    assert not (mismatch / 'result.csv').exists()
print(f'PASS: deposition/splash, evolved films, regrid/rebalance and MPI 2-to-4 restart; {work}')
