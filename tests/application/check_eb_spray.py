#!/usr/bin/env python3
"""Focused input-to-CSV EB spray application check; retain each unique run."""
import argparse
import csv
import math
from pathlib import Path
import shutil
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--exe', type=Path, required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
p.add_argument('--selected', action='store_true')
p.add_argument('--ctu', action='store_true')
p.add_argument('--inflow', action='store_true')
p.add_argument('--inlet-composition', action='store_true')
p.add_argument('--periodic', action='store_true')
p.add_argument('--wall', action='store_true')
p.add_argument('--scheduled', action='store_true')
p.add_argument('--restart', action='store_true')
p.add_argument('--split', action='store_true')
p.add_argument('--khrt', action='store_true')
p.add_argument('--outflow', action='store_true')
p.add_argument('--force', action='store_true')
p.add_argument('--difmag', type=float)
p.add_argument('--state-redist-order', type=int, choices=(0, 2, 4))
p.add_argument('--impacts', choices=('deposit', 'splash'))
a = p.parse_args()
assert not a.inlet_composition or (a.inflow and a.restart)
if a.inflow or a.periodic or a.wall:
    a.ctu = True
if a.periodic or a.wall:
    assert not (a.periodic and a.wall)
    assert not any((a.inflow, a.scheduled, a.split, a.outflow, a.force))
a.work.mkdir(parents=True, exist_ok=True)
work = Path(tempfile.mkdtemp(prefix='eb_spray_', dir=a.work.resolve()))
for name in ('gas.nml', 'spray.nml', 'level_set.dat'):
    shutil.copy(a.source / 'cases/eb_spray_3d' / name, work / name)
if a.ctu:
    gas = (work / 'gas.nml').read_text().replace('&reactive_eb_3d',
        "&reactive_eb_3d\n hydro_method='plm_ctu', riemann_solver='pelec',", 1)
    (work / 'gas.nml').write_text(gas)
if a.inflow:
    gas = (work / 'gas.nml').read_text().replace('&reactive_eb_3d',
        "&reactive_eb_3d\n boundary_faces='outflow','outflow','inflow','outflow','outflow','outflow',\n"
        " inflow_velocity(2,3)=200,", 1)
    (work / 'gas.nml').write_text(gas)
    spray = (work / 'spray.nml').read_text().replace('&eb_spray', '&eb_spray\n particle_outflow=.true.,', 1)
    (work / 'spray.nml').write_text(spray)
if a.inlet_composition:
    gas = (work / 'gas.nml').read_text().replace('&reactive_eb_3d',
        "&reactive_eb_3d\n inflow_composition_count(3)=2, "
        "inflow_composition_species(1:2,3)='O2','N2', "
        "inflow_composition_mole_fractions(1:2,3)=.2,.8,", 1)
    (work / 'gas.nml').write_text(gas)
if a.periodic or a.wall:
    gas = (work / 'gas.nml').read_text().replace('&reactive_eb_3d',
        "&reactive_eb_3d\n boundary_faces='outflow','outflow','periodic','periodic','periodic','periodic',", 1)
    if a.wall:
        gas = gas.replace("'periodic'", "'isothermal_wall'").replace('&reactive_eb_3d',
            '&reactive_eb_3d\n wall_temperature=1200,1200,1200,1200,1200,1200,', 1)
    (work / 'gas.nml').write_text(gas)
    spray = (work / 'spray.nml').read_text().replace('&eb_spray', '&eb_spray\n particle_outflow=.true.,', 1)
    if a.impacts:
        spray = spray.replace('center=.001250001,.0015,.0015, axis=-1,0,0',
                              'center=.001250001,.0039999999,.002, axis=-1,1,0')
    else:
        spray = spray.replace('center=.001250001,.0015,.0015, axis=-1,0,0',
                              'center=.003,.003999999,.002, axis=0,1,0')
    (work / 'spray.nml').write_text(spray)
if a.selected:
    gas = (work / 'gas.nml').read_text().replace("thermo_model='elementary'", "thermo_model='selected'")
    gas = gas.replace(" output_file=", " composition_count=1, composition_species(1)='N2', composition_mole_fractions(1)=1\n output_file=")
    (work / 'gas.nml').write_text(gas)
if a.scheduled:
    spray = (work / 'spray.nml').read_text().replace('&eb_spray',
        '&eb_spray\n injection_events=5, injection_times=5e-9,5e-9,5e-9,1e-8,2e-8,\n'
        ' injection_scales=1,0,1,.5,2,')
    (work / 'spray.nml').write_text(spray)
if a.impacts:
    spray = (work / 'spray.nml').read_text().replace('&eb_spray',
        '&eb_spray\n wall_impacts=.true., walls%evolve=.true., walls%maximum_films=32,\n'
        ' options%surface_tension=.072, options%liquid_viscosity=.001, options%maximum_substeps=8,')
    if a.impacts == 'splash':
        spray = spray.replace('speed=20', 'speed=200').replace('total_mass=1e-10', 'total_mass=1e-16')
    (work / 'spray.nml').write_text(spray)
if a.split:
    assert not a.impacts
    spray = (work / 'spray.nml').read_text().replace('&eb_spray',
        "&eb_spray\n options%breakup_model='etab', options%maximum_multiplicity=60,\n"
        ' options%surface_tension=1000, options%liquid_viscosity=.001,\n'
        ' options%drag=.false., options%heat=.false., options%evaporation=.false.,')
    (work / 'spray.nml').write_text(spray)
if a.khrt:
    assert not a.impacts and not a.scheduled and not a.split
    gas = (work / 'gas.nml').read_text().replace('regular_density=.31, cut_density=.31', 'regular_density=1, cut_density=1')
    gas = gas.replace('initial_pressure=135000', 'initial_pressure=435483.87096774194')
    (work / 'gas.nml').write_text(gas)
    spray = (work / 'spray.nml').read_text().replace('&eb_spray',
        "&eb_spray\n options%breakup_model='khrt', options%khrt_children=.true.,\n"
        ' options%khrt_reference_diameter_cubed=1e-25, options%surface_tension=.072, options%liquid_viscosity=.001,\n'
        ' options%drag=.false., options%heat=.false., options%evaporation=.false.,')
    spray = spray.replace('center=.001250001,.0015,.0015, axis=-1,0,0', 'center=.002,.0015,.0015, axis=1,0,0')
    spray = spray.replace('speed=20', 'speed=1000').replace('total_mass=1e-10', 'total_mass=1e-14')
    (work / 'spray.nml').write_text(spray)
if a.outflow:
    assert not a.scheduled and not a.khrt and not a.split
    assert a.impacts in (None, "splash")
    spray = (work / 'spray.nml').read_text().replace('&eb_spray', '&eb_spray\n particle_outflow=.true.,')
    if a.impacts:
        spray = spray.replace('center=.001250001,.0015,.0015', 'center=.001250001,.0015,.003999999')
        spray = spray.replace('maximum_dt=1e-8', 'maximum_dt=5e-9')
    else:
        spray = spray.replace('center=.001250001,.0015,.0015, axis=-1,0,0',
                              'center=.003999999,.0015,.0015, axis=1,0,0')
    (work / 'spray.nml').write_text(spray)
if a.force:
    spray = (work / 'spray.nml').read_text().replace('&eb_spray', '&eb_spray\n acceleration=0,0,1e6,', 1)
    (work / 'spray.nml').write_text(spray)
labels = ('Injected inventory:', 'Boundary inventory:', 'Wall inventory:', 'Chemistry inventory:')
if a.force:
    labels += ('Force inventory:',)
exe = str(a.exe.resolve())
if a.difmag is not None:
    gas = (work / 'gas.nml').read_text().replace('&reactive_eb_3d', f'&reactive_eb_3d\n difmag={a.difmag},', 1)
    (work / 'gas.nml').write_text(gas)
if a.state_redist_order is not None:
    gas = (work / 'gas.nml').read_text().replace('&reactive_eb_3d',
        f'&reactive_eb_3d\n state_redist_max_order={a.state_redist_order},', 1)
    (work / 'gas.nml').write_text(gas)
result = subprocess.run([exe, 'gas.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
(work / 'run.log').write_text(result.stdout + result.stderr)
assert result.returncode == 0, (work, result.stderr)
with (work / 'gas.csv.particles.csv').open() as stream:
    rows = list(csv.DictReader(stream))
if a.outflow and not a.impacts:
    assert [int(row['id']) for row in rows] == [1, 2]
    assert all(float(row['mass']) == 0 and float(row['component_mass_1']) == 0 for row in rows)
    boundary_line = next(x for x in result.stdout.splitlines() if x.startswith('Boundary inventory:'))
    assert float(boundary_line.split(':', 1)[1].split()[0]) < 0
elif a.wall and not a.impacts and not a.khrt:
    assert len(rows) == 2
    assert all(.003 < float(row['y']) < .004 and float(row['mass']) > 0 for row in rows)
    assert all(float(row['temperature']) > 350 and float(row['v']) < 0 for row in rows)
    assert all(math.isfinite(float(value)) for row in rows for value in row.values())
elif a.periodic and not a.impacts and not a.khrt:
    assert len(rows) == 2
    assert all(0 <= float(row['y']) < .001 and float(row['mass']) > 0 for row in rows)
    assert all(float(row['temperature']) > 350 and float(row['v']) > 0 for row in rows)
    assert all(math.isfinite(float(value)) for row in rows for value in row.values())
elif a.khrt:
    assert len(rows)>2 and (len(rows)-2)%4==0
    assert [int(row['id']) for row in rows] == list(range(1,len(rows)+1))
    if a.periodic:
        assert all(float(row['v'])>0 and 0 <= float(row['y']) < .001 for row in rows)
    elif a.wall:
        assert all(float(row['v'])<0 and .003 < float(row['y']) < .004 for row in rows)
    else:
        assert all(float(row['u'])>0 and float(row['x'])>.002 for row in rows)
    assert all(float(row['component_mass_1']) == float(row['mass']) for row in rows)
elif a.split:
    assert len(rows) == (12 if a.scheduled else 4)
    assert [int(row['id']) for row in rows] == list(range(1, len(rows)+1))
    evolved = rows[:-2] if a.scheduled else rows
    assert all(0 < float(row['multiplicity']) <= 60 for row in evolved)
    assert all(float(row['u']) > 0 for row in evolved)
    assert all(float(row['component_mass_1']) == float(row['mass']) for row in rows)
    if a.scheduled:
        assert all(float(row['multiplicity']) > 60 and float(row['u']) == -20 for row in rows[-2:])
elif not a.impacts:
    assert len(rows) == (8 if a.scheduled else 2)
    assert [int(row["id"]) for row in rows] == list(range(1, len(rows) + 1))
    initial_mass = math.pi / 6 * 1000 * 1e-15
    for row in rows:
        assert None not in row, row
        assert all(math.isfinite(float(value)) for value in row.values())
        if a.scheduled and int(row['id']) > 6:
            assert float(row['x']) == .001250001 and float(row['u']) == -20
            assert float(row['temperature']) == 350
            assert math.isclose(float(row['mass']), initial_mass, rel_tol=1e-14)
        else:
            assert float(row['x']) > .00125 and float(row['u']) > 0
            assert 0 < float(row['mass']) < initial_mass
            assert float(row['temperature']) > 350
        assert float(row['component_mass_1']) == float(row['mass'])
else:
    with (work / 'gas.csv.films.csv').open() as stream:
        films = list(csv.DictReader(stream))
    if a.impacts == 'deposit':
        assert films, work
    if a.outflow:
        assert films and 0 < len(rows) < 8, work
        assert all(3 <= int(row['id']) <= 10 for row in rows)
        boundary_line = next(x for x in result.stdout.splitlines() if x.startswith('Boundary inventory:'))
        assert float(boundary_line.split(':', 1)[1].split()[0]) < 0
    else:
        assert len(rows) == ((2 if a.scheduled else 0) if a.impacts == 'deposit' else (26 if a.scheduled else 8))
    ids = [int(row['id']) for row in rows]
    assert ids == sorted(set(ids))
    for row in rows + films:
        assert None not in row and all(math.isfinite(float(value)) for value in row.values())
        assert float(row['mass']) > 0 and float(row['component_mass_1']) == float(row['mass'])
    for film in films:
        assert abs(float(film['x']) - .00125) < 1e-15
        assert float(film['nx']) == 1 and float(film['height']) > 0
        if a.periodic:
            assert 0 <= float(film['y']) < .001
        if a.wall:
            assert .003 < float(film['y']) < .0039999999
    if a.scheduled:
        assert [int(row['id']) for row in rows[-2:]] == ([7, 8] if a.impacts == 'deposit' else [31, 32])
        for row in rows[-2:]:
            assert float(row['temperature']) == 350
            assert float(row['u']) == (-20 if a.impacts == 'deposit' else -200)
    if a.impacts == 'splash':
        for row in rows[:-2] if a.scheduled else rows:
            assert float(row['u']) > 0 and float(row['x']) > .00125
with (work / 'gas.csv').open() as stream:
    gas_rows = list(csv.DictReader(stream))
assert len(gas_rows) == 64
assert all(abs(float(row['time']) - 2e-8) < 1e-20 for row in gas_rows)
for label in labels:
    line = next(line for line in result.stdout.splitlines() if line.startswith(label))
    assert all(math.isfinite(float(x)) for x in line.split(':', 1)[1].split())
line = next(line for line in result.stdout.splitlines() if line.startswith('Injected inventory:'))
injected = [float(x) for x in line.split(':', 1)[1].split()]
assert math.isclose(injected[0], (4.5 if a.scheduled else 1) * (1e-14 if a.khrt else 1e-16 if a.impacts == 'splash' else 1e-10), rel_tol=1e-14)
if a.force:
    force = next(line for line in result.stdout.splitlines() if line.startswith('Force inventory:'))
    assert float(force.split(':', 1)[1].split()[3]) > 0
if a.scheduled:
    assert 'Injection events processed: 5' in result.stdout
    # Tied positive events have identical evolution, and the zero event consumes no IDs.
    for left, right in ([] if a.impacts or a.split else zip(rows[:2], rows[2:4])):
        assert all(left[key] == right[key] for key in left if key != 'id')
if a.restart:
    full_gas = (work / 'gas.csv').read_bytes()
    full_particles = (work / 'gas.csv.particles.csv').read_bytes()
    full_films = (work / 'gas.csv.films.csv').read_bytes()
    full_log = result.stdout
    base_gas = (work / 'gas.nml').read_text()
    split = base_gas.replace(" output_file=", " checkpoint_interval_steps=1, stop_after_checkpoint=.true., checkpoint_file='save.chk'\n output_file=")
    if a.impacts or a.split:
        split = split.replace('checkpoint_interval_steps=1', 'checkpoint_interval_steps=2')
    checkpoint_name = 'save.chk.step2' if a.impacts or a.split else 'save.chk.step1'
    (work / 'split.nml').write_text(split)
    run = subprocess.run([exe, 'split.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
    (work / 'save.log').write_text(run.stdout + run.stderr)
    assert run.returncode == 0, (work, run.stderr)
    checkpoint = work / checkpoint_name
    assert checkpoint.exists()
    resume = base_gas.replace(" output_file=", f" restart_file='{checkpoint_name}'\n output_file=")
    (work / 'resume.nml').write_text(resume)
    run = subprocess.run([exe, 'resume.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
    (work / 'resume.log').write_text(run.stdout + run.stderr)
    assert run.returncode == 0, (work, run.stderr)
    assert (work / 'gas.csv').read_bytes() == full_gas
    assert (work / 'gas.csv.particles.csv').read_bytes() == full_particles
    assert (work / 'gas.csv.films.csv').read_bytes() == full_films
    for label in labels:
        before = next(line for line in full_log.splitlines() if line.startswith(label))
        after = next(line for line in run.stdout.splitlines() if line.startswith(label))
        assert before == after, (label, work)
    if a.state_redist_order is not None:
        changed_order = 0 if a.state_redist_order else 2
        altered = resume.replace(f'state_redist_max_order={a.state_redist_order}',
                                 f'state_redist_max_order={changed_order}', 1)
        (work / 'changed_order.nml').write_text(altered)
        rejected = subprocess.run([exe, 'changed_order.nml', 'spray.nml'], cwd=work,
                                  capture_output=True, text=True)
        (work / 'changed_order.log').write_text(rejected.stdout + rejected.stderr)
        assert rejected.returncode != 0, work
        assert 'incompatible spray checkpoint' in (rejected.stdout + rejected.stderr).lower(), work
        assert (work / 'gas.csv').read_bytes() == full_gas
        assert (work / 'gas.csv.particles.csv').read_bytes() == full_particles
        assert (work / 'gas.csv.films.csv').read_bytes() == full_films
    if a.ctu:
        altered = resume.replace("hydro_method='plm_ctu'", "hydro_method='pcm'")
        (work / 'changed_hydro.nml').write_text(altered)
        rejected = subprocess.run([exe, 'changed_hydro.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
        (work / 'changed_hydro.log').write_text(rejected.stdout + rejected.stderr)
        assert rejected.returncode != 0
        assert (work / 'gas.csv').read_bytes() == full_gas
        assert (work / 'gas.csv.particles.csv').read_bytes() == full_particles
        assert (work / 'gas.csv.films.csv').read_bytes() == full_films
    if a.inlet_composition:
        for name, changed in [
            ('changed_composition', resume.replace('=.2,.8', '=.3,.7')),
            ('unknown_species', resume.replace("'O2','N2'", "'missing','N2'")),
            ('negative_fraction', resume.replace('=.2,.8', '=-.2,.8')),
        ]:
            (work / f'{name}.nml').write_text(changed)
            rejected = subprocess.run([exe, f'{name}.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
            (work / f'{name}.log').write_text(rejected.stdout + rejected.stderr)
            assert rejected.returncode != 0, name
        assert 'incompatible spray checkpoint' in (work / 'changed_composition.log').read_text().lower()
    if a.inflow:
        for name, altered in [
            ('changed_inlet', resume.replace('inflow_velocity(2,3)=200', 'inflow_velocity(2,3)=201')),
            ('changed_faces', resume.replace("'inflow'", "'outflow'")),
            ('invalid_boundary', resume.replace("'inflow'", "'invalid'")),
        ]:
            (work / f'{name}.nml').write_text(altered)
            rejected = subprocess.run([exe, f'{name}.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
            (work / f'{name}.log').write_text(rejected.stdout + rejected.stderr)
            assert rejected.returncode != 0, name
            assert (work / 'gas.csv').read_bytes() == full_gas
            assert (work / 'gas.csv.particles.csv').read_bytes() == full_particles
            assert (work / 'gas.csv.films.csv').read_bytes() == full_films
        (work / 'closed_particles.nml').write_text((work / 'spray.nml').read_text().replace(
            'particle_outflow=.true.', 'particle_outflow=.false.'))
        rejected = subprocess.run([exe, 'gas.nml', 'closed_particles.nml'], cwd=work, capture_output=True, text=True)
        (work / 'closed_particles.log').write_text(rejected.stdout + rejected.stderr)
        assert rejected.returncode != 0 and 'require particle_outflow' in rejected.stderr
    if a.periodic:
        for name, altered in [
            ('changed_periodic', resume.replace("'periodic'", "'outflow'")),
            ('unpaired_periodic', resume.replace("'periodic'", "'outflow'", 1)),
            ('incompatible_geometry', resume.replace("'outflow','outflow'", "'periodic','periodic'")),
        ]:
            (work / f'{name}.nml').write_text(altered)
            rejected = subprocess.run([exe, f'{name}.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
            (work / f'{name}.log').write_text(rejected.stdout + rejected.stderr)
            assert rejected.returncode != 0, name
            assert (work / 'gas.csv').read_bytes() == full_gas
            assert (work / 'gas.csv.particles.csv').read_bytes() == full_particles
            assert (work / 'gas.csv.films.csv').read_bytes() == full_films
    if a.wall:
        for name, altered in [
            ('changed_wall', resume.replace('1200', '1201')),
            ('invalid_wall', resume.replace('1200', '-1')),
            ('changed_wall_type', resume.replace("'isothermal_wall'", "'no_slip_wall'")),
        ]:
            (work / f'{name}.nml').write_text(altered)
            rejected = subprocess.run([exe, f'{name}.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
            (work / f'{name}.log').write_text(rejected.stdout + rejected.stderr)
            assert rejected.returncode != 0, name
            assert (work / 'gas.csv').read_bytes() == full_gas
            assert (work / 'gas.csv.particles.csv').read_bytes() == full_particles
            assert (work / 'gas.csv.films.csv').read_bytes() == full_films
    if a.difmag is not None:
        for name, value in [('changed_difmag', a.difmag + .01), ('negative_difmag', -1), ('nan_difmag', 'NaN')]:
            altered = resume.replace(f'difmag={a.difmag}', f'difmag={value}')
            (work / f'{name}.nml').write_text(altered)
            rejected = subprocess.run([exe, f'{name}.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
            (work / f'{name}.log').write_text(rejected.stdout + rejected.stderr)
            assert rejected.returncode != 0, name
            assert (work / 'gas.csv').read_bytes() == full_gas
            assert (work / 'gas.csv.particles.csv').read_bytes() == full_particles
            assert (work / 'gas.csv.films.csv').read_bytes() == full_films
    # Immutable publication, context change and a truncated payload must reject.
    run = subprocess.run([exe, 'split.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
    (work / 'overwrite.log').write_text(run.stdout + run.stderr)
    assert run.returncode != 0
    changed = (work / 'spray.nml').read_text().replace('liquid_temperature=350', 'liquid_temperature=351')
    (work / 'changed.nml').write_text(changed)
    run = subprocess.run([exe, 'resume.nml', 'changed.nml'], cwd=work, capture_output=True, text=True)
    (work / 'context.log').write_text(run.stdout + run.stderr)
    assert run.returncode != 0
    geometry = (work / 'level_set.dat').read_bytes()
    lines = geometry.decode().splitlines()
    lines[1:] = [' '.join(f'{float(value) + 1e-8:.16e}' for value in line.split()) for line in lines[1:]]
    (work / 'level_set.dat').write_text('\n'.join(lines) + '\n')
    run = subprocess.run([exe, 'resume.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
    (work / 'geometry_context.log').write_text(run.stdout + run.stderr)
    (work / 'level_set.dat').write_bytes(geometry)
    assert run.returncode != 0
    payload = checkpoint.read_bytes()
    (work / 'truncated.chk').write_bytes(payload[:-50])
    (work / 'truncated.nml').write_text(resume.replace(checkpoint_name, 'truncated.chk'))
    run = subprocess.run([exe, 'truncated.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
    (work / 'truncated.log').write_text(run.stdout + run.stderr)
    assert run.returncode != 0 and checkpoint.read_bytes() == payload
    if a.khrt:
        lines = payload.decode().splitlines()
        assert int(lines[2].split()[-1]) > 2
        assert (int(lines[4])-2) % 4 == 0
        lines[4] = str(int(lines[4])+1)
        (work / 'bad_ids.chk').write_text('\n'.join(lines)+'\n')
        (work / 'bad_ids.nml').write_text(resume.replace(checkpoint_name,'bad_ids.chk'))
        run = subprocess.run([exe,'bad_ids.nml','spray.nml'],cwd=work,capture_output=True,text=True)
        (work / 'bad_ids.log').write_text(run.stdout+run.stderr)
        assert run.returncode != 0 and 'Invalid restart child IDs' in run.stderr
    if a.split:
        lines = payload.decode().splitlines()
        assert lines[0].endswith('_SCOUNT_V1')
        assert int(lines[5]) == (4 if a.scheduled else 2)
        lines[5] = str(int(lines[5]) - 1)
        (work / 'bad_count.chk').write_text('\n'.join(lines) + '\n')
        (work / 'bad_count.nml').write_text(resume.replace(checkpoint_name, 'bad_count.chk'))
        run = subprocess.run([exe, 'bad_count.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
        (work / 'bad_count.log').write_text(run.stdout + run.stderr)
        assert run.returncode != 0 and 'Invalid restart injection IDs' in run.stderr

if a.impacts == 'splash':
    wet = (work / 'spray.nml').read_text().replace('total_mass=1e-16', 'total_mass=1e-10')
    wet = wet.replace('center=.001250001,.0015,.003999999', 'center=.001250001,.0015,.0015')
    (work / 'invalid_wet.nml').write_text(wet)
    run = subprocess.run([exe, 'gas.nml', 'invalid_wet.nml'], cwd=work,
                         capture_output=True, text=True, timeout=10)
    (work / 'wet_rejection.log').write_text(run.stdout + run.stderr)
    assert run.returncode != 0 and 'EB impact correlation outside admissible' in run.stderr
# An output alias must fail before overwriting the spray input (absolute/relative mix).
original = (work / 'spray.nml').read_bytes()
gas = (work / 'gas.nml').read_text().replace("output_file='gas.csv'", f"output_file='{work / 'spray.nml'}'")
(work / 'bad.nml').write_text(gas)
result = subprocess.run([exe, 'bad.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
(work / 'alias.log').write_text(result.stdout + result.stderr)
assert result.returncode != 0 and (work / 'spray.nml').read_bytes() == original
if a.scheduled:
    invalid = original.decode().replace('5e-9,5e-9,5e-9,1e-8,2e-8', '5e-9,0,5e-9,1e-8,2e-8')
    (work / 'invalid_spray.nml').write_text(invalid)
    result = subprocess.run([exe, 'gas.nml', 'invalid_spray.nml'], cwd=work, capture_output=True, text=True)
    (work / 'invalid_schedule.log').write_text(result.stdout + result.stderr)
    assert result.returncode != 0 and 'nondecreasing' in result.stderr
print(f'PASS: EB spray application, reflection/transfer, final CSVs and alias rejection: {work}')
