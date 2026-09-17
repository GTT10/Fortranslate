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
p.add_argument('--scheduled', action='store_true')
p.add_argument('--restart', action='store_true')
p.add_argument('--impacts', choices=('deposit', 'splash'))
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)
work = Path(tempfile.mkdtemp(prefix='eb_spray_', dir=a.work.resolve()))
for name in ('gas.nml', 'spray.nml', 'level_set.dat'):
    shutil.copy(a.source / 'cases/eb_spray_3d' / name, work / name)
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
exe = str(a.exe.resolve())
result = subprocess.run([exe, 'gas.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
(work / 'run.log').write_text(result.stdout + result.stderr)
assert result.returncode == 0, (work, result.stderr)
with (work / 'gas.csv.particles.csv').open() as stream:
    rows = list(csv.DictReader(stream))
if not a.impacts:
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
    assert len(rows) == ((2 if a.scheduled else 0) if a.impacts == 'deposit' else (26 if a.scheduled else 8))
    ids = [int(row['id']) for row in rows]
    assert ids == sorted(set(ids))
    for row in rows + films:
        assert None not in row and all(math.isfinite(float(value)) for value in row.values())
        assert float(row['mass']) > 0 and float(row['component_mass_1']) == float(row['mass'])
    for film in films:
        assert abs(float(film['x']) - .00125) < 1e-15
        assert float(film['nx']) == 1 and float(film['height']) > 0
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
for label in ('Injected inventory:', 'Boundary inventory:', 'Wall inventory:', 'Chemistry inventory:'):
    line = next(line for line in result.stdout.splitlines() if line.startswith(label))
    assert all(math.isfinite(float(x)) for x in line.split(':', 1)[1].split())
line = next(line for line in result.stdout.splitlines() if line.startswith('Injected inventory:'))
injected = [float(x) for x in line.split(':', 1)[1].split()]
assert math.isclose(injected[0], (4.5 if a.scheduled else 1) * (1e-16 if a.impacts == 'splash' else 1e-10), rel_tol=1e-14)
if a.scheduled:
    assert 'Injection events processed: 5' in result.stdout
    # Tied positive events have identical evolution, and the zero event consumes no IDs.
    for left, right in ([] if a.impacts else zip(rows[:2], rows[2:4])):
        assert all(left[key] == right[key] for key in left if key != 'id')
if a.restart:
    full_gas = (work / 'gas.csv').read_bytes()
    full_particles = (work / 'gas.csv.particles.csv').read_bytes()
    full_films = (work / 'gas.csv.films.csv').read_bytes()
    full_log = result.stdout
    base_gas = (work / 'gas.nml').read_text()
    split = base_gas.replace(" output_file=", " checkpoint_interval_steps=1, stop_after_checkpoint=.true., checkpoint_file='save.chk'\n output_file=")
    if a.impacts:
        split = split.replace('checkpoint_interval_steps=1', 'checkpoint_interval_steps=2')
    checkpoint_name = 'save.chk.step2' if a.impacts else 'save.chk.step1'
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
    for label in ('Injected inventory:', 'Boundary inventory:', 'Wall inventory:', 'Chemistry inventory:'):
        before = next(line for line in full_log.splitlines() if line.startswith(label))
        after = next(line for line in run.stdout.splitlines() if line.startswith(label))
        assert before == after, (label, work)
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
if a.impacts == 'splash':
    wet = (work / 'spray.nml').read_text().replace('total_mass=1e-16', 'total_mass=1e-10')
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
