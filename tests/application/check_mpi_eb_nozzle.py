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
p.add_argument('--rate', action='store_true')
p.add_argument('--rate-files', action='store_true')
p.add_argument('--properties', action='store_true')
a = p.parse_args()
if a.rate_files and not a.rate:
    p.error('--rate-files requires --rate')
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


def run(label, ranks, extra='', invalid=False, changed=False, inline=False):
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
    if a.rate:
        controls = controls.replace(
            'injection_events=4,event_ports=1,2,1,2,injection_times=0.,0.,2.e-8,2.e-8,',
            'injection_events=5,event_ports=1,2,1,1,2,injection_times=0.,0.,1.e-8,2.e-8,2.e-8,'
            'injection_scales=1.,1.,0.,.5,2.,')
        controls += ('port_velocity_from_rate=2*.true.,port_discharge_coefficients=.8,.5,'
                     'port_profile_points=2,2,port_profile_times(:,1)=0.,4.e-8,62*0.,'
                     'port_profile_times(:,2)=0.,4.e-8,62*0.,port_profile_rates(:,1)=2*.002,62*0.,'
                     'port_profile_rates(:,2)=.001,.003,62*0.,'
                     'injection_end_times=5.e-9,5.e-9,1.5e-8,2.5e-8,2.5e-8,')
    if a.properties:
        values = values.replace('nc=1,', "nc=2,vapor_names(2)='H2',")
        controls += ("port_fractions(:,1)=.25,.75,14*0.,port_fractions(:,2)=.75,.25,14*0.,"
                     "mixture_density_models=2*'polynomial',"
                     "mixture_density_coefficients(:,1)=1100.,-.2,0.,0.,"
                     "mixture_density_coefficients(:,2)=900.,-.1,0.,0.,"
                     "mixture_density_temperature_ranges(:,1)=300.,450.,"
                     "mixture_density_temperature_ranges(:,2)=300.,450.,"
                     "mixture_cp_models=2*'polynomial',mixture_cp_coefficients(:,1)=3800.,1.,0.,0.,"
                     "mixture_cp_coefficients(:,2)=1800.,.5,0.,0.,"
                     "mixture_cp_temperature_ranges(:,1)=300.,450.,mixture_cp_temperature_ranges(:,2)=300.,450.,"
                     "mixture_vapor_pressure_models(1)='antoine',"
                     "mixture_antoine_coefficients(:,1)=8.07131,1730.63,-39.724,133.322,"
                     "mixture_antoine_temperature_ranges(:,1)=300.,450.,")
    file_group = ''
    if a.rate_files and not inline:
        controls = controls.replace(
            'port_profile_points=2,2,port_profile_times(:,1)=0.,4.e-8,62*0.,'
            'port_profile_times(:,2)=0.,4.e-8,62*0.,port_profile_rates(:,1)=2*.002,62*0.,'
            'port_profile_rates(:,2)=.001,.003,62*0.,', '')
        # Renamed files and comments on continuation must not change checkpoint identity.
        first = 'renamed.dat' if label == 'resumed' else 'rate1.dat'
        rate = '-.002' if invalid else '.0021' if changed else '.002'
        (directory / first).write_text(f'# {label}\n0 .002\n4.e-8 {rate} ! end\n')
        (directory / 'rate2.dat').write_text('0\t.001\n\n4.e-8 .003\n')
        file_group = f"\n&eb_spray_rate_files profile_files='{first}','rate2.dat' /\n"
    elif invalid and a.properties:
        controls = controls.replace('1100.,-.2,0.,0.', '-1100.,0.,0.,0.')
    elif invalid:
        controls = controls.replace('port_radii=2*.0001', 'port_radii=2*-.0001')
    spray.write_text(values.replace('&eb_spray', '&eb_spray\n' + controls, 1) + file_group)
    result = subprocess.run([str(a.mpiexec.resolve()), '-n', str(ranks), str(a.exe.resolve()), 'tree.nml', 'spray.nml'],
                            cwd=directory, env=env, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=60)
    (directory / 'run.log').write_text(result.stdout)
    assert (result.returncode == 0) != (invalid or changed), result.stdout
    return directory


full = [run(f'full{r}', r) for r in (2, 4)]
for directory in full:
    data = particles(directory)
    assert [int(r['id']) for r in data] == list(range(1, 33))
    for port, start in ((1, 16), (2, 24)):
        fresh = data[start:start+8]
        sizes, counts, doses = [], [], []
        dose = (5e-12 if port == 1 else 2.125e-11) if a.rate else 1e-10
        fractions = (.25, .75) if port == 1 else (.75, .25)
        density = 1 / (fractions[0]/1030 + fractions[1]/865) if a.properties else 1000
        expected_speed = dose / (5e-9 * density * math.pi * .0001**2 * (.8 if port == 1 else .5)) if a.rate else 20
        for row in fresh:
            assert int(row['material_id']) == port
            x = float(row['x']) - (.0015 if port == 1 else .003)
            z = float(row['z']) - .0015
            radial = math.hypot(x, z)
            assert 0 < radial <= .0001 * (1 + 1e-12)
            if port == 2:
                assert abs(radial-.0001) <= 1e-12*.0001
            u, w = float(row['u']), float(row['w'])
            speed = math.sqrt(u*u + float(row['v'])**2 + w*w)
            assert abs(speed-expected_speed) <= 1e-12*expected_speed
            assert abs((x*u+z*w)/(radial*math.hypot(u, w))-math.cos(.4)) < 1e-12
            assert float(row['y']) == .0015 and float(row['temperature']) == 350
            mass, count = float(row['mass']), float(row['multiplicity'])
            sizes.append((6*mass/(math.pi*density))**(1/3))
            if a.properties:
                for c, fraction in enumerate(fractions, 1):
                    assert abs(float(row[f'component_mass_{c}'])/mass-fraction) < 1e-13
            counts.append(count)
            doses.append(mass*count)
        assert all(5e-6 <= d <= 2e-5 for d in sizes) and max(sizes) > 1.1*min(sizes)
        values = doses if port == 1 else counts
        assert max(values)-min(values) <= 1e-13*max(values)
        assert abs(sum(doses)-dose) <= 1e-13*dose
    assert any(float(r['temperature']) != 350 for r in data[:16])
compare(rows(full[0] / 'result.csv'), rows(full[1] / 'result.csv'))
compare(particles(full[0]), particles(full[1]))
stopped = run('stopped', 2, "checkpoint_file='saved',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,")
resumed = run('resumed', 4, f"restart_file='{stopped / 'saved.step1'}',")
compare(rows(full[1] / 'result.csv'), rows(resumed / 'result.csv'))
compare(particles(full[1]), particles(resumed))
if a.rate_files:
    baseline = run('inline', 2, inline=True)
    compare(rows(full[0] / 'result.csv'), rows(baseline / 'result.csv'))
    compare(particles(full[0]), particles(baseline))
    rejected = run('changed', 2, f"restart_file='{stopped / 'saved.step1'}',", changed=True)
    assert 'EB spray context mismatch' in (rejected / 'run.log').read_text()
    assert not (rejected / 'result.csv').exists()
invalid = run('invalid', 2, invalid=True)
expected_error = ('Invalid EB spray rate file' if a.rate_files else
                  'Invalid injected port mixture' if a.properties else
                  'Rate-driven outlet requires positive radius' if a.rate else 'Invalid injection port template')
assert expected_error in (invalid / 'run.log').read_text()
assert not (invalid / 'result.csv').exists()
print(f'PASS: finite ports, size distributions, rate velocity={a.rate}, property models={a.properties}, MPI 2/4 parity and restart; {work}')
