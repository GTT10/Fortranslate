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
p.add_argument('--khrt', action='store_true')
p.add_argument('--domain-wall', action='store_true')
p.add_argument('--periodic', action='store_true')
p.add_argument('--pcm', action='store_true')
p.add_argument('--split', action='store_true')
p.add_argument('--outflow', action='store_true')
p.add_argument('--specular-outflow', action='store_true')
p.add_argument('--force', action='store_true')
a = p.parse_args()
if a.periodic:
    assert a.khrt and not a.domain_wall, '--periodic requires KH without --domain-wall'
if a.pcm:
    assert a.domain_wall, '--pcm requires --domain-wall'
if a.domain_wall:
    assert a.khrt, '--domain-wall requires --khrt'
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
    if a.domain_wall:
        if not a.pcm:
            tree = tree.replace("riemann_solver='hllc'", "riemann_solver='pelec',hydro_method='plm_ctu'")
        controls += "boundary_faces='no_slip_wall','no_slip_wall','outflow','outflow','outflow','outflow',"
    if a.periodic:
        tree = tree.replace("riemann_solver='hllc'", "riemann_solver='pelec',hydro_method='plm_ctu'")
        controls += "boundary_faces='periodic','periodic','outflow','outflow','outflow','outflow',"
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
    if mode == 'khrt':
        spray = spray.replace('wall_impacts=.true.', 'wall_impacts=.false.').replace('walls%evolve=.true.', 'walls%evolve=.false.')
        spray = spray.replace('speed=20.', 'speed=1000.').replace('total_mass=1.e-10', 'total_mass=1.e-16')
        spray = spray.replace('axis=0.,-1.,0.', 'axis=1.,0.,0.').replace('center=.0015,.001250001,.0015',
                                                                        'center=.0015,.0015,.0015')
        spray = spray.replace('&eb_spray', "&eb_spray\n particle_options%breakup_model='khrt',"
                              'particle_options%khrt_children=.true.,particle_options%khrt_reference_diameter_cubed=1.e-25,'
                              'particle_options%drag=.false.,particle_options%heat=.false.,particle_options%evaporation=.false.,', 1)
    if a.domain_wall or a.periodic:
        spray = spray.replace('&eb_spray', '&eb_spray\n particle_outflow=.true.,', 1)
        spray = spray.replace('center=.0015,.0015,.0015', 'center=.003999999,.0015,.0015')
    if mode == 'split':
        spray = spray.replace('wall_impacts=.true.', 'wall_impacts=.false.').replace('walls%evolve=.true.', 'walls%evolve=.false.')
        spray = spray.replace('particle_options%surface_tension=.072', 'particle_options%surface_tension=1000.')
        spray = spray.replace('axis=0.,-1.,0.', 'axis=1.,0.,0.').replace('center=.0015,.001250001,.0015',
                                                                        'center=.0015,.0015,.0015')
        spray = spray.replace('&eb_spray', "&eb_spray\n particle_options%breakup_model='etab',"
                              'particle_options%maximum_multiplicity=60.,'
                              'particle_options%drag=.false.,particle_options%heat=.false.,particle_options%evaporation=.false.,', 1)
    if mode == 'outflow':
        spray = spray.replace('wall_impacts=.true.', 'wall_impacts=.false.').replace('walls%evolve=.true.', 'walls%evolve=.false.')
        spray = spray.replace('&eb_spray', '&eb_spray\n particle_outflow=.true.,', 1)
        spray = spray.replace('center=.0015,.001250001,.0015', 'center=.003999999,.0015,.0015')
        spray = spray.replace('axis=0.,-1.,0.', 'axis=1.,0.,0.')
    if a.outflow:
        assert mode == 'splash'
        spray = spray.replace('&eb_spray', '&eb_spray\n particle_outflow=.true.,', 1)
        spray = spray.replace('center=.0015,.001250001,.0015', 'center=.0015,.001250001,.003999999')
    if a.force:
        spray = spray.replace('&eb_spray', '&eb_spray\n acceleration=0,0,1e6,', 1)
    if altered:
        spray = spray.replace('6*400.', '6*401.')
    (directory / 'spray.nml').write_text(spray)
    result = subprocess.run([str(a.mpiexec.resolve()), '-n', str(ranks), str(a.exe.resolve()), 'tree.nml', 'spray.nml'],
                            cwd=directory, env=env, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=60)
    (directory / 'run.log').write_text(result.stdout)
    assert (result.returncode == 0) == success, result.stdout
    return directory


for mode in (('deposit',) if a.force else ('outflow',) if a.specular_outflow else ('splash',) if a.outflow else ('split',) if a.split else ('khrt',) if a.khrt else ('deposit', 'splash')):
    full = run(mode, 'full', 4)
    films = population(full, 'films')
    particles = population(full, 'particles')
    if mode == 'outflow':
        assert not films and not particles
    elif mode == 'split':
        assert len(films) == 0 and len(particles) == 4
        assert sorted(int(r['id']) for r in particles) == [1, 2, 3, 4]
        assert all(0 < float(r['multiplicity']) <= 60 for r in particles)
    elif mode == 'khrt':
        assert len(films) == 0 and len(particles) > 2 and (len(particles)-2) % 4 == 0
        assert sorted(int(r['id']) for r in particles) == list(range(1, len(particles)+1))
        if a.domain_wall:
            assert all(float(r['u']) < 0 and 0 < float(r['x']) < .004 for r in particles)
        elif a.periodic:
            assert all(float(r['u']) > 0 and 0 <= float(r['x']) < .001 for r in particles)
        else:
            assert all(float(r['u']) > 0 and float(r['x']) > .0015 for r in particles)
    elif a.outflow:
        assert len(films) == 2 and 0 < len(particles) < 8
    else:
        assert len(films) == 2 and len(particles) == (0 if mode == 'deposit' else 8)
    for film in films:
        assert float(film['mass']) > 0 and float(film['temperature']) > 350
        assert float(film['ny']) == 1 and float(film['y']) == .00125
        assert float(film['component_mass_1']) == float(film['mass'])
    if particles and mode not in ('khrt', 'split'):
        ids = sorted(int(r['id']) for r in particles)
        assert ids == sorted(set(ids)) and all(3 <= i <= 10 for i in ids)
        if not a.outflow:
            assert ids == list(range(3, 11))
        assert all(float(r['v']) > 0 and float(r['y']) > .00125 for r in particles)
    if a.force:
        assert float(rows(full / 'result.csv.budget.csv')[3]['force']) > 0
    log = (full / 'run.log').read_text()
    assert all(f'Regrid stage: {n}' in log for n in (1, 2, 3)), log
    stopped = run(mode, 'stopped', 2, "checkpoint_file='saved',checkpoint_interval_steps=2,stop_after_checkpoint=.true.,")
    archive = stopped / 'saved.step2'
    payload = Path(str(archive) + '.particles').read_text().splitlines()
    assert payload[0] == ('PELEF_EB_PARTICLES_V1' if mode in ('khrt', 'split', 'outflow') else 'PELEF_EB_PARTICLES_V2')
    if mode == 'khrt':
        np = int(payload[2].split()[2])
        stored = [line.split() for line in payload[4:4+np]]
        assert any(int(row[1]) > 0 for row in stored)
        assert any(float(row[15]) > 0 for row in stored)
        assert all(float(row[16]) >= 0 for row in stored)
    resumed = run(mode, 'resumed', 4, f"restart_file='{archive}',")
    compare(rows(full / 'result.csv'), rows(resumed / 'result.csv'))
    compare(particles, population(resumed, 'particles'))
    compare(films, population(resumed, 'films'))
    for x, y in zip(rows(full / 'result.csv.budget.csv'), rows(resumed / 'result.csv.budget.csv')):
        compare([{k: x[k] for k in ('initial', 'final', 'boundary', 'wall', 'chemistry', 'injected', 'force')}],
                [{k: y[k] for k in ('initial', 'final', 'boundary', 'wall', 'chemistry', 'injected', 'force')}])
    mismatch = run(mode, 'mismatch', 2, f"restart_file='{archive}',", success=False, altered=True)
    assert 'EB spray context mismatch' in (mismatch / 'run.log').read_text()
    assert not (mismatch / 'result.csv').exists()
    if mode == 'split':
        gas_archive = Path(str(archive) + '.gas')
        original = gas_archive.read_bytes()
        lines = original.decode().splitlines()
        index = next(i for i, line in enumerate(lines) if 'GENERATION:' in line) + 1
        metadata = lines[index].split()
        assert [int(x) for x in metadata[1:]] == [4, 2]
        metadata[2] = '1'
        lines[index] = ' '.join(metadata)
        gas_archive.write_text('\n'.join(lines) + '\n')
        try:
            invalid = run(mode, 'bad_split_count', 2, f"restart_file='{archive}',", success=False)
            assert 'Invalid restart ID watermark' in (invalid / 'run.log').read_text()
            assert not (invalid / 'result.csv').exists()
        finally:
            gas_archive.write_bytes(original)
print(f'PASS: KH mode={a.khrt}, split={a.split}, particle/film state, regrid routing and MPI 2-to-4 restart; {work}')
