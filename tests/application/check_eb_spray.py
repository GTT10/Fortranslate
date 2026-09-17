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
exe = str(a.exe.resolve())
result = subprocess.run([exe, 'gas.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
(work / 'run.log').write_text(result.stdout + result.stderr)
assert result.returncode == 0, (work, result.stderr)
with (work / 'gas.csv.particles.csv').open() as stream:
    rows = list(csv.DictReader(stream))
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
with (work / 'gas.csv').open() as stream:
    gas_rows = list(csv.DictReader(stream))
assert len(gas_rows) == 64
assert all(abs(float(row['time']) - 2e-8) < 1e-20 for row in gas_rows)
for label in ('Injected inventory:', 'Boundary inventory:', 'Wall inventory:', 'Chemistry inventory:'):
    line = next(line for line in result.stdout.splitlines() if line.startswith(label))
    assert all(math.isfinite(float(x)) for x in line.split(':', 1)[1].split())
line = next(line for line in result.stdout.splitlines() if line.startswith('Injected inventory:'))
injected = [float(x) for x in line.split(':', 1)[1].split()]
assert math.isclose(injected[0], 4.5e-10 if a.scheduled else 1e-10, rel_tol=1e-14)
if a.scheduled:
    assert 'Injection events processed: 5' in result.stdout
    # Tied positive events have identical evolution, and the zero event consumes no IDs.
    for left, right in zip(rows[:2], rows[2:4]):
        assert all(left[key] == right[key] for key in left if key != 'id')
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
