"""Serial EB finite solid/hollow injection, size weighting and restart."""
import argparse
import csv
import math
from pathlib import Path
import shutil
import subprocess
import tempfile

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--exe', type=Path, required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)
work = Path(tempfile.mkdtemp(prefix='nozzle_', dir=a.work.resolve()))


def run(profile, label, extra='', invalid=False):
    directory = work / profile / label
    shutil.copytree(a.source / 'cases/eb_spray_3d', directory)
    gas = directory / 'gas.nml'
    gas.write_text(gas.read_text().replace('&reactive_eb_3d', '&reactive_eb_3d\n' + extra, 1))
    spray = directory / 'spray.nml'
    values = spray.read_text().replace('number=2', 'number=8').replace(
        'center=.001250001,.0015,.0015, axis=-1,0,0', 'center=.002,.0015,.0015, axis=1,0,0')
    controls = ("injection_events=2,injection_times=0.,2.e-8,nozzle_radius=.0001,half_angle=.1,"
                f"nozzle%profile='{profile}',nozzle%swirl=.4,"
                "sizes%minimum=5.e-6,sizes%maximum=2.e-5,")
    if profile == 'solid':
        controls += 'sizes%shape=2.5,'
    else:
        controls += "nozzle%spread=.02,sizes%distribution='lognormal',sizes%weighting='number',sizes%std=2.e-6,"
    if invalid:
        controls = controls.replace('nozzle_radius=.0001', 'nozzle_radius=-.0001')
    spray.write_text(values.replace('&eb_spray', '&eb_spray\n' + controls, 1))
    result = subprocess.run([str(a.exe.resolve()), 'gas.nml', 'spray.nml'], cwd=directory,
                            capture_output=True, text=True, timeout=30)
    (directory / 'run.log').write_text(result.stdout + result.stderr)
    assert (result.returncode == 0) != invalid, (directory, result.stderr)
    return directory


for profile in ('solid', 'hollow'):
    full = run(profile, 'full')
    with (full / 'gas.csv.particles.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert [int(r['id']) for r in rows] == list(range(1, 17))
    fresh = rows[8:]
    diameters, doses, counts = [], [], []
    for row in fresh:
        assert all(math.isfinite(float(v)) for v in row.values())
        dy, dz = float(row['y'])-.0015, float(row['z'])-.0015
        radius = math.hypot(dy, dz)
        assert float(row['x']) == .002 and radius > 0
        if profile == 'hollow':
            assert math.isclose(radius, .0001, rel_tol=1e-13)
        else:
            assert radius < .0001
        assert abs(dy*float(row['w'])-dz*float(row['v'])) > 1e-8
        mass, count = float(row['mass']), float(row['multiplicity'])
        diameter = (6*mass/(math.pi*1000))**(1/3)
        assert 5e-6 <= diameter <= 2e-5
        assert float(row['component_mass_1']) == mass
        diameters.append(diameter)
        doses.append(mass*count)
        counts.append(count)
    assert max(diameters) > min(diameters)*1.1
    assert math.isclose(sum(doses), 1e-10, rel_tol=1e-13)
    weights = doses if profile == 'solid' else counts
    assert max(weights)-min(weights) <= 1e-13*max(weights)
    stopped = run(profile, 'stopped', "checkpoint_interval_steps=1,stop_after_checkpoint=.true.,checkpoint_file='saved',")
    archive = stopped / 'saved.step1'
    resumed = run(profile, 'resumed', f"restart_file='{archive}',")
    for filename in ('gas.csv', 'gas.csv.particles.csv', 'gas.csv.films.csv'):
        assert (full / filename).read_bytes() == (resumed / filename).read_bytes()
    for label in ('Injected inventory:', 'Boundary inventory:', 'Wall inventory:', 'Chemistry inventory:'):
        before = next(x for x in (full / 'run.log').read_text().splitlines() if x.startswith(label))
        after = next(x for x in (resumed / 'run.log').read_text().splitlines() if x.startswith(label))
        assert before == after
    bad = run(profile, 'invalid', invalid=True)
    assert 'Invalid scheduled EB injection' in (bad / 'run.log').read_text()
    assert not (bad / 'gas.csv').exists()
print(f'PASS: serial EB finite nozzles, size weights and restart: {work}')
