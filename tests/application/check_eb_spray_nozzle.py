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
p.add_argument('--properties', action='store_true')
p.add_argument('--rate', action='store_true')
p.add_argument('--rate-files', action='store_true')
p.add_argument('--ports', action='store_true')
p.add_argument('--models', action='store_true')
a = p.parse_args()
if a.rate_files and not (a.rate and a.ports):
    p.error('--rate-files requires --rate --ports')
a.work.mkdir(parents=True, exist_ok=True)
work = Path(tempfile.mkdtemp(prefix='nozzle_', dir=a.work.resolve()))


def run(profile, label, extra='', invalid=False, changed=False, inline=False):
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
    if a.ports:
        assert a.properties
        controls = controls.replace('injection_events=2,injection_times=0.,2.e-8,',
            'port_count=2,injection_events=5,injection_times=0.,0.,1.e-8,2.e-8,2.e-8,'
            'injection_scales=1.,1.,0.,1.,1.,event_ports=1,2,2,1,2,'
            'port_numbers(2)=4,port_centers(:,2)=.003,.0015,.0015,port_masses(2)=2.e-10,'
            'port_speeds(2)=30.,port_fractions(:,2)=.75,.25,14*0.,')
    if a.rate:
        if a.ports:
            controls += ('injection_profile_points=2,injection_profile_times=0.,4.e-8,injection_profile_rates=2*.002,'
                         'injection_end_times=5.e-9,5.e-9,1.5e-8,2.5e-8,2.5e-8,'
                         'velocity_from_rate=.true.,discharge_coefficient=.8,'
                         'port_profile_points(2)=2,port_profile_times(1:2,2)=0.,4.e-8,'
                         'port_profile_rates(1:2,2)=.001,.003,port_radii(2)=.0002,port_discharge_coefficients(2)=.5,')
        else:
            controls = controls.replace('injection_events=2,injection_times=0.,2.e-8,',
                'injection_events=3,injection_times=0.,1.e-8,2.e-8,injection_scales=1.,0.,2.,')
            controls += ('injection_profile_points=2,injection_profile_times=0.,4.e-8,'
                         'injection_end_times=5.e-9,1.5e-8,2.5e-8,velocity_from_rate=.true.,discharge_coefficient=.8,')
            controls += 'injection_profile_rates=2*.002,' if profile == 'solid' else 'injection_profile_rates=.001,.003,'
    if a.properties:
        values = values.replace("nc=1, number=8, vapor_names(1)='H2O', fractions(1)=1",
                                "nc=2, number=8, vapor_names='H2O','H2', fractions=.25,.75")
        controls += ("mixture_density_models=2*'polynomial',"
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
    if a.models:
        assert a.ports
        if profile == 'solid':
            controls += ("port_nozzles(2)%profile='hollow',port_nozzles(2)%spread=.02,"
                         "port_sizes(2)%distribution='lognormal',port_sizes(2)%weighting='number',"
                         "port_sizes(2)%std=2.e-6,port_sizes(2)%shape=0.,")
        else:
            controls += ("port_nozzles(2)%profile='solid',port_nozzles(2)%spread=0.,"
                         "port_sizes(2)%distribution='weibull',port_sizes(2)%weighting='mass',"
                         "port_sizes(2)%shape=3.,port_sizes(2)%std=0.,")
    file_group = ''
    if a.rate_files and not inline:
        controls = controls.replace(
            'injection_profile_points=2,injection_profile_times=0.,4.e-8,injection_profile_rates=2*.002,', '')
        controls = controls.replace(
            'port_profile_points(2)=2,port_profile_times(1:2,2)=0.,4.e-8,'
            'port_profile_rates(1:2,2)=.001,.003,', '')
        first = 'renamed.dat' if label == 'resumed' else 'shared.dat'
        rate = '-.002' if invalid else '.0021' if changed else '.002'
        (directory / first).write_text(f'# {label}\n0 .002\n4.e-8 {rate}\n')
        (directory / 'port2.dat').write_text('0\t.001\n\n4.e-8 .003 ! end\n')
        file_group = f"\n&eb_spray_rate_files profile_file='{first}',profile_files(2)='port2.dat' /\n"
    elif invalid and a.models:
        controls += "port_nozzles(2)%profile='invalid',"
    elif invalid and a.properties:
        controls = controls.replace('1100.,-.2,0.,0.', '-1100.,0.,0.,0.')
    elif invalid:
        controls = controls.replace('nozzle_radius=.0001', 'nozzle_radius=-.0001')
    spray.write_text(values.replace('&eb_spray', '&eb_spray\n' + controls, 1) + file_group)
    result = subprocess.run([str(a.exe.resolve()), 'gas.nml', 'spray.nml'], cwd=directory,
                            capture_output=True, text=True, timeout=30)
    (directory / 'run.log').write_text(result.stdout + result.stderr)
    assert (result.returncode == 0) != (invalid or changed), (directory, result.stderr)
    return directory


for profile in ('solid', 'hollow'):
    full = run(profile, 'full')
    with (full / 'gas.csv.particles.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert [int(r['id']) for r in rows] == list(range(1, 25 if a.ports else 17))
    fresh = rows[12:] if a.ports else rows[8:]
    dose = (2e-11 if profile == 'solid' else 2.125e-11) if a.rate else 1e-10
    diameters, doses, counts = [], [], []
    for row in fresh:
        assert all(math.isfinite(float(v)) for v in row.values())
        dy, dz = float(row['y'])-.0015, float(row['z'])-.0015
        radius = math.hypot(dy, dz)
        second = a.ports and int(row['material_id']) == 2
        outlet_radius = .0002 if second and a.rate else .0001
        assert float(row['x']) == (.003 if second else .002) and radius > 0
        hollow = (profile == 'hollow') != (a.models and second)
        if hollow:
            assert math.isclose(radius, outlet_radius, rel_tol=1e-13)
        else:
            assert radius < outlet_radius
        assert abs(dy*float(row['w'])-dz*float(row['v'])) > 1e-8
        mass, count = float(row['mass']), float(row['multiplicity'])
        fraction = .75 if second else .25
        density = 1/(fraction/1030+(1-fraction)/865) if a.properties else 1000
        diameter = (6*mass/(math.pi*density))**(1/3)
        if a.rate:
            event_dose = (1.0625e-11 if second else 1e-11) if a.ports else dose
            coefficient = .5 if second else .8
            expected_speed = event_dose/(5e-9*density*math.pi*outlet_radius**2*coefficient)
            actual_speed = math.sqrt(sum(float(row[k])**2 for k in ('u','v','w')))
            assert math.isclose(actual_speed, expected_speed, rel_tol=1e-12)
        assert 5e-6 <= diameter <= 2e-5
        if a.properties:
            assert math.isclose(float(row['component_mass_1'])/mass, fraction, rel_tol=1e-14)
            assert math.isclose(float(row['component_mass_2'])/mass, 1-fraction, rel_tol=1e-14)
        else:
            assert float(row['component_mass_1']) == mass
        diameters.append(diameter)
        doses.append(mass*count)
        counts.append(count)
    assert max(diameters) > min(diameters)*1.1
    if a.ports:
        assert [int(r['material_id']) for r in rows] == [1]*8+[2]*4+[1]*8+[2]*4
        assert math.isclose(sum(doses[8:]), 1.0625e-11 if a.rate else 2e-10, rel_tol=1e-13)
        assert math.isclose(sum(doses[:8]), 1e-11 if a.rate else 1e-10, rel_tol=1e-13)
        speeds = [math.sqrt(sum(float(r[k])**2 for k in ('u','v','w'))) for r in fresh]
        if not a.rate:
            assert all(math.isclose(v,20 if i<8 else 30,rel_tol=1e-13) for i,v in enumerate(speeds))
    else:
        assert math.isclose(sum(doses), dose, rel_tol=1e-13)
    weights = doses if profile == 'solid' else counts
    if a.models:
        groups = (weights[:8], (counts if profile == 'solid' else doses)[8:])
    else:
        groups = (weights[:8], weights[8:]) if a.ports else (weights,)
    for group in groups:
        assert max(group)-min(group) <= 1e-13*max(group)
    stopped = run(profile, 'stopped', "checkpoint_interval_steps=1,stop_after_checkpoint=.true.,checkpoint_file='saved',")
    archive = stopped / 'saved.step1'
    resumed = run(profile, 'resumed', f"restart_file='{archive}',")
    for filename in ('gas.csv', 'gas.csv.particles.csv', 'gas.csv.films.csv'):
        assert (full / filename).read_bytes() == (resumed / filename).read_bytes()
    if a.rate_files:
        baseline = run(profile, 'inline', inline=True)
        for filename in ('gas.csv', 'gas.csv.particles.csv', 'gas.csv.films.csv'):
            assert (full / filename).read_bytes() == (baseline / filename).read_bytes()
        rejected = run(profile, 'changed', f"restart_file='{archive}',", changed=True)
        assert 'EB spray restart rejected' in (rejected / 'run.log').read_text()
        assert not (rejected / 'gas.csv').exists()
    for label in ('Injected inventory:', 'Boundary inventory:', 'Wall inventory:', 'Chemistry inventory:'):
        before = next(x for x in (full / 'run.log').read_text().splitlines() if x.startswith(label))
        after = next(x for x in (resumed / 'run.log').read_text().splitlines() if x.startswith(label))
        assert before == after
    bad = run(profile, 'invalid', invalid=True)
    assert ('Invalid EB spray rate file' if a.rate_files else 'Invalid scheduled EB injection' if a.models else 'Invalid injected liquid mixture' if a.properties else 'Invalid rate-driven outlet radius' if a.rate else 'Invalid scheduled EB injection') in (bad / 'run.log').read_text()
    assert not (bad / 'gas.csv').exists()
print(f'PASS: serial EB finite nozzles, size weights and restart: {work}')
