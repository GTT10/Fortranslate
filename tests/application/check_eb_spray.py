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
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)
work = Path(tempfile.mkdtemp(prefix='eb_spray_', dir=a.work.resolve()))
for name in ('gas.nml', 'spray.nml', 'level_set.dat'):
    shutil.copy(a.source / 'cases/eb_spray_3d' / name, work / name)
if a.selected:
    gas = (work / 'gas.nml').read_text().replace("thermo_model='elementary'", "thermo_model='selected'")
    gas = gas.replace(" output_file=", " composition_count=1, composition_species(1)='N2', composition_mole_fractions(1)=1\n output_file=")
    (work / 'gas.nml').write_text(gas)
exe = str(a.exe.resolve())
result = subprocess.run([exe, 'gas.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
(work / 'run.log').write_text(result.stdout + result.stderr)
assert result.returncode == 0, (work, result.stderr)
with (work / 'gas.csv.particles.csv').open() as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 2
initial_mass = math.pi / 6 * 1000 * 1e-15
for row in rows:
    assert None not in row, row
    assert all(math.isfinite(float(value)) for value in row.values())
    assert float(row['x']) > .00125 and float(row['u']) > 0
    assert 0 < float(row['mass']) < initial_mass
    assert float(row['temperature']) > 350
    assert float(row['component_mass_1']) == float(row['mass'])
with (work / 'gas.csv').open() as stream:
    gas_rows = list(csv.DictReader(stream))
assert len(gas_rows) == 64
assert all(abs(float(row['time']) - 2e-8) < 1e-20 for row in gas_rows)
for label in ('Boundary inventory:', 'Wall inventory:', 'Chemistry inventory:'):
    line = next(line for line in result.stdout.splitlines() if line.startswith(label))
    assert all(math.isfinite(float(x)) for x in line.split(':', 1)[1].split())
# An output alias must fail before overwriting the spray input (absolute/relative mix).
original = (work / 'spray.nml').read_bytes()
gas = (work / 'gas.nml').read_text().replace("output_file='gas.csv'", f"output_file='{work / 'spray.nml'}'")
(work / 'bad.nml').write_text(gas)
result = subprocess.run([exe, 'bad.nml', 'spray.nml'], cwd=work, capture_output=True, text=True)
(work / 'alias.log').write_text(result.stdout + result.stderr)
assert result.returncode != 0 and (work / 'spray.nml').read_bytes() == original
print(f'PASS: EB spray application, reflection/transfer, final CSVs and alias rejection: {work}')
