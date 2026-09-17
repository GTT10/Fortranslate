"""Bounded CLI exercise: injection, coupled steps, regrid and rank parity."""
import argparse
import csv
import math
from pathlib import Path
import subprocess

p = argparse.ArgumentParser()
p.add_argument('--exe', type=Path, required=True)
p.add_argument('--mpiexec', required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
p.add_argument('--selected', action='store_true')
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)
case = a.source / 'cases/mpi_amr_spray_3d'


def run(name, ranks, empty=False):
    work = a.work / name
    work.mkdir(exist_ok=True)
    gas = (case / 'gas.nml').read_text()
    if a.selected:
        gas = gas.replace("thermo_model='full_h2o2'", "thermo_model='selected'")
        gas = gas.replace(" output_file=", " composition_count=3,composition_species='H2','O2','N2'\n"
                          " composition_mole_fractions=2.,1.,3.76\n output_file=")
    (work / 'gas.nml').write_text(gas)
    spray = (case / 'spray.nml').read_text()
    if empty:
        spray = spray.replace('parcel_count=4,total_mass=1.e-10', 'parcel_count=0,total_mass=0.')
    spray = spray.replace('\n/\n', f"\n regrid_time=1.e-8,regrid_topology='{case / 'root.dat'}'\n/\n")
    (work / 'spray.nml').write_text(spray)
    result = subprocess.run([a.mpiexec, '-n', str(ranks), str(a.exe.resolve()), 'gas.nml',
                             str(case / 'tree.dat'), 'spray.nml'], cwd=work,
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
    (work / 'run.log').write_text(result.stdout)
    assert result.returncode == 0, result.stdout
    assert 'Regridded: T' in result.stdout, result.stdout
    with (work / 'amr_spray.csv').open() as f:
        gas_rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]
    assert len(gas_rows) == 64 and all(row['level'] == 0 for row in gas_rows)
    assert all(row['time'] == 2e-8 for row in gas_rows)
    particles = []
    for rank in range(ranks):
        with (work / f'amr_spray.csv.particles.{rank}.csv').open() as f:
            rows = list(csv.reader(f))[1:]
            particles.extend([[float(x) for x in row] for row in rows])
    particles.sort(key=lambda row: row[0])
    assert len(particles) == (0 if empty else 4)
    for row in particles:
        assert all(math.isfinite(x) for x in row)
        assert abs(row[10] - 350) > 1e-8 and row[3] > .00075, row
    if not empty:
        assert math.fsum(row[9] * row[11] for row in particles) < 1e-10
    gas_rows.sort(key=lambda row: (row['x'], row['y'], row['z']))
    return gas_rows, particles


def integrate(rows, key):
    return math.fsum(r[key] * r['dx'] * r['dy'] * r['dz'] for r in rows)


baseline, _ = run('empty', 2, True)
gas2, particles2 = run('two', 2)
gas4, particles4 = run('four', 4)
for gas, particles in [(gas2, particles2), (gas4, particles4)]:
    mass = integrate(gas, 'rho') + math.fsum(r[9] * r[11] for r in particles)
    expected = integrate(baseline, 'rho') + 1e-10
    assert abs(mass - expected) <= 3e-12 * expected, (mass, expected)
    for c, (name, fraction) in enumerate([('H2O', .25), ('H2O2', .75)]):
        total = integrate(gas, 'rhoY_' + name) + math.fsum(r[12+c] * r[11] for r in particles)
        expected = integrate(baseline, 'rhoY_' + name) + fraction * 1e-10
        assert abs(total - expected) <= 3e-12 * expected, (name, total, expected)
for left, right in zip(gas2, gas4):
    for key in left:
        assert math.isclose(left[key], right[key], rel_tol=2e-12, abs_tol=1e-25), (key, left[key], right[key])
for left, right in zip(particles2, particles4):
    for x, y in zip(left, right):
        assert math.isclose(x, y, rel_tol=2e-12, abs_tol=1e-25), (x, y)
print('PASS: AMR spray CLI, initial injection, regrid, mass/species budget and 2/4-rank parity')
