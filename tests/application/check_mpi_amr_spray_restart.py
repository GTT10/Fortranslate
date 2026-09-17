"""Changed-rank AMR spray continuation and checkpoint rejection checks."""
import argparse
import csv
import math
from pathlib import Path
import shutil
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--exe', type=Path, required=True)
p.add_argument('--mpiexec', required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
p.add_argument('--selected', action='store_true')
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)
work = Path(tempfile.mkdtemp(prefix='restart_', dir=a.work.resolve()))
case = a.source / 'cases/mpi_amr_spray_3d'
gas = (case / 'gas.nml').read_text()
if a.selected:
    gas = gas.replace("thermo_model='full_h2o2'", "thermo_model='selected'")
    gas = gas.replace(' output_file=', " composition_count=3,composition_species='H2','O2','N2'\n"
                      ' composition_mole_fractions=2.,1.,3.76\n output_file=')
spray = (case / 'spray.nml').read_text().replace('\n/\n',
    f"\n regrid_time=1.e-8,regrid_topology='{case / 'root.dat'}'\n"
    ' injection_events=3,injection_times=0.,1.e-8,1.5e-8,event_parcels=4,4,4,event_mass=3.e-11,3.e-11,4.e-11\n/\n')


def run(name, ranks, controls, good=True, changed=False):
    folder = work / name
    folder.mkdir()
    (folder / 'gas.nml').write_text(gas)
    local_spray = spray.replace('speed=20.', 'speed=21.') if changed else spray
    (folder / 'spray.nml').write_text(local_spray + '\n&amr_spray_checkpoint\n' + controls + '\n/\n')
    result = subprocess.run([a.mpiexec, '-n', str(ranks), str(a.exe.resolve()), 'gas.nml',
                             str(case / 'tree.dat'), 'spray.nml'], cwd=folder,
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
    (folder / 'run.log').write_text(result.stdout)
    assert (result.returncode == 0) == good, result.stdout
    if not good:
        assert not (folder / 'amr_spray.csv').exists()
    return folder


def read(folder, ranks):
    with (folder / 'amr_spray.csv').open() as f:
        fields = [[float(x) for x in r] for r in list(csv.reader(f))[1:]]
    particles = []
    for rank in range(ranks):
        with (folder / f'amr_spray.csv.particles.{rank}.csv').open() as f:
            particles += [[float(x) for x in r] for r in list(csv.reader(f))[1:]]
    return fields, sorted(particles)


full = run('full', 2, "checkpoint_file='state',checkpoint_time=1.e-8")
split = run('split', 2, "checkpoint_file='state',checkpoint_time=1.e-8,stop_after_checkpoint=.true.")
assert not (split / 'amr_spray.csv').exists()
checkpoint = split / 'state'
assert checkpoint.exists() and (split / 'state.gas').exists() and (split / 'state.particles').exists()
continued = run('continued', 4, f"restart_file='{checkpoint}'")
for expected, actual in zip(read(full, 2), read(continued, 4)):
    assert len(expected) == len(actual)
    for left, right in zip(expected, actual):
        assert len(left) == len(right)
        for x, y in zip(left, right):
            assert math.isclose(x, y, rel_tol=2e-12, abs_tol=1e-25), (x, y)
assert 'Injection events applied: 3' in (continued / 'run.log').read_text()
assert len(read(continued, 4)[1]) == 12
run('changed', 2, f"restart_file='{checkpoint}'", good=False, changed=True)
# A complete marker with a truncated particle payload must fail before CSV output.
for suffix in ('', '.gas', '.particles'):
    shutil.copyfile(str(checkpoint) + suffix, work / ('broken' + suffix))
particle_file = work / 'broken.particles'
particle_file.write_text(particle_file.read_text().replace('END_AMR_PARTICLES', ''))
run('truncated', 2, f"restart_file='{work / 'broken'}'", good=False)
# Immutable publication: an existing checkpoint must remain byte-for-byte intact.
saved = {s: Path(str(checkpoint)+s).read_bytes() for s in ('', '.gas', '.particles', '.pending')}
run('overwrite', 2, f"checkpoint_file='{checkpoint}',checkpoint_time=1.e-8", good=False)
assert all(Path(str(checkpoint)+s).read_bytes() == value for s, value in saved.items())
print(f'PASS: 2-to-4 rank restart, event continuation, context/truncation/overwrite rejection; artifacts {work}')
