"""Finite solid/hollow MPI EB ports with swirl, size distributions and restart."""
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
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)
work = Path(tempfile.mkdtemp(prefix='nozzle_', dir=a.work.resolve()))
env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def particles(directory):
    return sorted((r for f in directory.glob('result.csv.particles.rank*.csv') for r in rows(f)),
                  key=lambda r: int(r['id']))


def compare(left, right):
    assert len(left) == len(right)
    for x, y in zip(left, right):
        for key in x:
            u, v = float(x[key]), float(y[key])
            assert math.isfinite(u) and math.isfinite(v)
            assert abs(u-v) <= 3e-12 * max(abs(u), abs(v), 1e-30), (key, u, v)


def run(label, ranks, extra='', invalid=False):
    directory = work / label
    shutil.copytree(a.source / 'cases/mpi_eb_spray_3d', directory)
    tree = directory / 'tree.nml'
    tree.write_text(tree.read_text().replace('&eb_tree', '&eb_tree\n' + extra, 1))
    spray = directory / 'spray.nml'
    values = spray.read_text().replace('number=2', 'number=8').replace('center=.0015,.001250001,.0015',
                                                                       'center=.0015,.0015,.0015')
    controls = ("port_count=2,injection_events=4,event_ports=1,2,1,2,injection_times=0.,0.,2.e-8,2.e-8,"
                "half_angle=.1,port_centers(:,2)=.003,.0015,.0015,port_radii=2*.0001,"
                "port_nozzles(1)%profile='solid',port_nozzles(1)%swirl=.4,"
                "port_nozzles(2)%profile='hollow',port_nozzles(2)%swirl=.4,port_nozzles(2)%spread=.02,"
                "port_sizes(1)%shape=2.5,port_sizes(1)%minimum=5.e-6,port_sizes(1)%maximum=2.e-5,"
                "port_sizes(2)%distribution='lognormal',port_sizes(2)%weighting='number',"
                "port_sizes(2)%std=2.e-6,port_sizes(2)%minimum=5.e-6,port_sizes(2)%maximum=2.e-5,")
    if invalid:
        controls = controls.replace('port_radii=2*.0001', 'port_radii=2*-.0001')
    spray.write_text(values.replace('&eb_spray', '&eb_spray\n' + controls, 1))
    result = subprocess.run([str(a.mpiexec.resolve()), '-n', str(ranks), str(a.exe.resolve()), 'tree.nml', 'spray.nml'],
                            cwd=directory, env=env, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=60)
    (directory / 'run.log').write_text(result.stdout)
    assert (result.returncode == 0) != invalid, result.stdout
    return directory


full = [run(f'full{r}', r) for r in (2, 4)]
for directory in full:
    data = particles(directory)
    assert [int(r['id']) for r in data] == list(range(1, 33))
    for port, start in ((1, 16), (2, 24)):
        fresh = data[start:start+8]
        sizes, counts, doses = [], [], []
        for row in fresh:
            assert int(row['material_id']) == port
            x = float(row['x']) - (.0015 if port == 1 else .003)
            z = float(row['z']) - .0015
            radial = math.hypot(x, z)
            assert 0 < radial <= .0001 * (1 + 1e-12)
            if port == 2:
                assert abs(radial-.0001) <= 1e-12*.0001
            u, w = float(row['u']), float(row['w'])
            assert abs((x*u+z*w)/(radial*math.hypot(u, w))-math.cos(.4)) < 1e-12
            assert float(row['y']) == .0015 and float(row['temperature']) == 350
            mass, count = float(row['mass']), float(row['multiplicity'])
            sizes.append((6*mass/(math.pi*1000))**(1/3))
            counts.append(count)
            doses.append(mass*count)
        assert all(5e-6 <= d <= 2e-5 for d in sizes) and max(sizes) > 1.1*min(sizes)
        values = doses if port == 1 else counts
        assert max(values)-min(values) <= 1e-13*max(values)
        assert abs(sum(doses)-1e-10) <= 1e-23
    assert any(float(r['temperature']) != 350 for r in data[:16])
compare(rows(full[0] / 'result.csv'), rows(full[1] / 'result.csv'))
compare(particles(full[0]), particles(full[1]))
stopped = run('stopped', 2, "checkpoint_file='saved',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,")
resumed = run('resumed', 4, f"restart_file='{stopped / 'saved.step1'}',")
compare(rows(full[1] / 'result.csv'), rows(resumed / 'result.csv'))
compare(particles(full[1]), particles(resumed))
invalid = run('invalid', 2, invalid=True)
assert 'Invalid injection port template' in (invalid / 'run.log').read_text()
assert not (invalid / 'result.csv').exists()
print(f'PASS: solid/hollow ports, swirl, mass/number size distributions, MPI 2/4 parity and restart; {work}')
