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
p.add_argument('--ppm', action='store_true')
p.add_argument('--sdc', action='store_true')
p.add_argument('--mol', action='store_true')
p.add_argument('--mol-lf', action='store_true')
p.add_argument('--rotating', action='store_true')
p.add_argument('--sdc-cvode', action='store_true')
p.add_argument('--gas-force', action='store_true')
p.add_argument('--reacting', action='store_true')
p.add_argument('--difmag', action='store_true')
p.add_argument('--multiport', action='store_true')
p.add_argument('--adaptive', action='store_true')
p.add_argument('--rebalance', action='store_true')
p.add_argument('--physical', action='store_true')
p.add_argument('--inflow', action='store_true')
p.add_argument('--films', action='store_true')
p.add_argument('--les', action='store_true')
p.add_argument('--subcycling', action='store_true')
p.add_argument('--dynamic', action='store_true')
p.add_argument('--dynamic-open', action='store_true')
p.add_argument('--mach-inlet', action='store_true')
p.add_argument('--dynamic-walls', action='store_true')
p.add_argument('--dynamic-pressure', action='store_true')
p.add_argument('--pressure', action='store_true')
p.add_argument('--force', action='store_true')
p.add_argument('--properties', action='store_true')
p.add_argument('--pelemp-cap', action='store_true')
p.add_argument('--sizes', action='store_true')
p.add_argument('--rates', action='store_true')
p.add_argument('--discharge', action='store_true')
p.add_argument('--split', action='store_true')
p.add_argument('--rate-files', action='store_true')
a = p.parse_args()
assert not (a.sdc and a.mol), 'Select one time integrator'
staged = a.sdc or a.mol
assert not a.mol_lf or a.mol
assert not a.rotating or a.mol
a.work.mkdir(parents=True, exist_ok=True)
work = Path(tempfile.mkdtemp(prefix='restart_', dir=a.work.resolve()))
case = a.source.resolve() / 'cases/mpi_amr_spray_3d'
gas = (case / 'gas.nml').read_text()
if a.ppm:
    gas = gas.replace("reconstruction='pcm'", "reconstruction='characteristic_ppm'")
if staged:
    gas = gas.replace("reconstruction='pcm'", "reconstruction='characteristic_ppm'")
    gas = gas.replace("riemann_solver='hllc'", "riemann_solver='pelec',time_integrator='sdc',sdc_iterations=2")
if a.mol:
    gas = gas.replace("reconstruction='characteristic_ppm'", "reconstruction='characteristic_plm'")
    gas = gas.replace("time_integrator='sdc'", "time_integrator='mol'")
    if a.mol_lf:
        gas = gas.replace("time_integrator='mol'", "time_integrator='mol',mol_lax_friedrichs=.true.")
if a.rotating:
    gas = gas.replace("boundary_condition='periodic'", "boundary_condition='physical',"
                      "boundary_faces=4*'slip_wall',2*'periodic'")
    gas += "\n&reactive_rotation enabled=.true.,axis=2,angular_speed=200.,axis_point=.002,.002,0. /\n"
if a.sdc_cvode:
    assert staged
    gas = gas.replace("time_integrator='sdc'", "time_integrator='sdc',sdc_reactor='cvode'")
if a.reacting or a.sdc_cvode:
    gas = gas.replace('chemistry_enabled=.false.', 'chemistry_enabled=.true.')
if a.gas_force:
    assert staged
    gas = gas.replace(' output_file=', ' external_forcing=100.,-200.,300.\n output_file=')
if a.difmag:
    assert a.ppm
    gas = gas.replace("riemann_solver='hllc'", "riemann_solver='pelec',difmag=0.1")
if a.selected:
    gas = gas.replace("thermo_model='full_h2o2'", "thermo_model='selected'")
    gas = gas.replace(' output_file=', " composition_count=3,composition_species='H2','O2','N2'\n"
                      ' composition_mole_fractions=2.,1.,3.76\n output_file=')
spray = (case / 'spray.nml').read_text().replace('\n/\n',
    f"\n regrid_time=1.e-8,regrid_topology='{case / 'root.dat'}'\n"
    ' injection_events=3,injection_times=0.,1.e-8,1.5e-8,event_parcels=4,4,4,event_mass=3.e-11,3.e-11,4.e-11\n/\n')

if a.subcycling:
    assert spray.count('subcycling=.false.') == 1
    spray = spray.replace('subcycling=.false.', 'subcycling=.true.')

if a.multiport:
    spray = spray.replace('injection_events=3,injection_times=0.,1.e-8,1.5e-8,event_parcels=4,4,4,event_mass=3.e-11,3.e-11,4.e-11',
        'injection_events=4,injection_times=0.,1.e-8,1.e-8,1.5e-8,event_parcels=4,4,4,4,\n'
        ' event_mass=2.e-11,3.e-11,1.e-11,4.e-11,event_ports=1,1,2,2,port_count=2,\n'
        ' port_centers(:,1)=.00075,.00075,.00075,port_centers(:,2)=.00325,.00325,.00325,\n'
        ' port_axes(:,1)=1.,0.,0.,port_axes(:,2)=-1.,0.,0.,port_speeds=20.,30.,\n'
        ' port_half_angles=.1,.2,port_diameters=1.e-5,2.e-5,port_temperatures=350.,340.,\n'
        ' port_fractions(1:2,1)=.25,.75,port_fractions(1:2,2)=.75,.25,port_radii=1.e-5,2.e-5,\n'
        " port_nozzles(1)%profile='solid',port_nozzles(2)%profile='hollow',port_nozzles(2)%swirl=.2")

if a.adaptive:
    gas = gas.replace('hotspot_temperature_rise=0.', 'hotspot_temperature_rise=100.')
    spray = spray.replace(f"regrid_time=1.e-8,regrid_topology='{case / 'root.dat'}'",
        'adaptation_interval=1,refinement_levels=1,refinement_ratios=2,\n'
        ' temperature_thresholds=1005.,tile_shape=1,1,1,tag_buffer=0,maximum_nodes=100')

if a.rebalance:
    spray = spray.replace('\n/\n', '\n rebalance_interval=1,rebalance_minimum_gain=0.\n/\n')

if a.physical:
    assert a.multiport, 'physical boundary fixture uses both ports'
    gas = gas.replace("boundary_condition='periodic'",
        "boundary_condition='periodic',boundary_faces='slip_wall','outflow','periodic','periodic','periodic','periodic'")
    spray = spray.replace('port_centers(:,1)=.00075,.00075,.00075,port_centers(:,2)=.00325,.00325,.00325',
        'port_centers(:,1)=.000001,.00075,.00075,port_centers(:,2)=.003999,.00325,.00325')
    spray = spray.replace('port_axes(:,1)=1.,0.,0.,port_axes(:,2)=-1.,0.,0.,port_speeds=20.,30.',
        'port_axes(:,1)=-1.,0.,0.,port_axes(:,2)=1.,0.,0.,port_speeds=1000.,1000.')
    spray = spray.replace('port_half_angles=.1,.2', 'port_half_angles=0.,0.')

if a.inflow:
    assert a.physical
    gas = gas.replace("'slip_wall','outflow'", "'inflow','outflow'")
    gas = gas.replace(" output_file=", " inflow_temperature=1000.,inflow_pressure=101325.,inflow_velocity(:,1)=25.,0.,0.\n output_file=")
    spray = spray.replace('port_centers(:,1)=.000001', 'port_centers(:,1)=.00075')
    spray = spray.replace('port_axes(:,1)=-1.', 'port_axes(:,1)=1.')

if a.pressure:
    assert a.physical
    gas = gas.replace(' output_file=', ' outlet_pressure(2)=100000.\n output_file=')

if a.films:
    assert a.physical and not a.inflow
    spray = spray.replace('port_centers(:,1)=.000001', 'port_centers(:,1)=1.e-9')
    spray = spray.replace('port_speeds=1000.,1000.', 'port_speeds=1.,1000.')
    spray = spray.replace('\n/\n', '\n walls%enabled(1)=.true.,walls%temperature=300.,walls%evolve=.true.,\n'
        ' walls%liquid_conductivity=.6,walls%minimum_height=1.e-10,\n'
        ' options%surface_tension=.072,options%liquid_viscosity=.001\n/\n')

if a.les:
    spray += "\n&amr_spray_les les_cs=.17,les_ci=.1,les_model='vreman' /\n"
if a.dynamic:
    assert a.les and not a.physical and not a.pressure
    spray = spray.replace("les_cs=.17,les_ci=.1,les_model='vreman'", "les_model='dynamic'")
    spray += "\n&amr_spray_dynamic les_test_model='gaussian',les_test_ratio=2,les_grid_ratio=1 /\n"
if a.dynamic_open:
    assert a.dynamic
    gas = gas.replace(' output_file=', " boundary_faces='inflow','outflow',4*'periodic'\n"
                      ' inflow_temperature=1000.,inflow_pressure=101325.,inflow_velocity(:,1)=25.,0.,0.\n output_file=')

if a.mach_inlet:
    assert a.dynamic_open and not a.dynamic_pressure
    gas = gas.replace('inflow_pressure=101325.', 'inflow_pressure=0.,inflow_mach(1)=.2')
    gas = gas.replace('inflow_velocity(:,1)=25.,0.,0.', 'inflow_velocity(:,1)=0.,0.,0.')

if a.dynamic_pressure:
    assert a.dynamic_open
    gas = gas.replace(' output_file=', ' outlet_pressure(2)=100000.\n output_file=')

if a.dynamic_walls:
    assert a.dynamic and not a.dynamic_open
    gas = gas.replace("transport_enabled=.false.", "transport_enabled=.true.,viscosity_enabled=.true.")
    gas = gas.replace(' output_file=', " boundary_faces=2*'slip_wall',2*'no_slip_wall',2*'isothermal_wall'\n"
                      ' wall_temperature=1100.,wall_velocity(1,3:6)=2.\n output_file=')


if a.force:
    spray += '\n&amr_spray_force acceleration=100.,-200.,300. /\n'


if a.properties:
    spray += ("\n&amr_spray_properties\n"
              " mixture_density_models=2*'polynomial',\n"
              " mixture_density_coefficients(:,1)=1100.,-.2,0.,0.,\n"
              " mixture_density_coefficients(:,2)=900.,-.1,0.,0.,\n"
              " mixture_density_temperature_ranges(:,1)=300.,450.,\n"
              " mixture_density_temperature_ranges(:,2)=300.,450.,\n"
              " mixture_cp_models=2*'polynomial',mixture_cp_coefficients(:,1)=3800.,1.,0.,0.,\n"
              " mixture_cp_coefficients(:,2)=1800.,.5,0.,0.,\n"
              " mixture_cp_temperature_ranges(:,1)=300.,450.,mixture_cp_temperature_ranges(:,2)=300.,450.,\n"
              " mixture_vapor_pressure_models(1)='antoine',\n"
              " mixture_antoine_coefficients(:,1)=8.07131,1730.63,-39.724,133.322,\n"
              " mixture_antoine_temperature_ranges(:,1)=300.,450.\n/\n")

if a.pelemp_cap:
    spray = spray.replace('drop_temperature=350.', 'drop_temperature=400.')
    spray += ("\n&amr_spray_properties\n"
              " mixture_boiling_models='pelemp-cap','pelemp-cap',\n"
              " mixture_critical_temperatures=647.,700.,\n"
              " mixture_normal_boiling_temperatures=373.15,450.\n/\n")


if a.sizes:
    assert a.multiport
    spray += ("\n&amr_spray_sizes\n"
              " port_sizes(1)%distribution='normal',port_sizes(1)%std=2.e-6,\n"
              " port_sizes(1)%minimum=2.e-6,port_sizes(1)%maximum=3.e-5,\n"
              " port_sizes(2)%shape=2.5,port_sizes(2)%weighting='number',\n"
              " port_sizes(2)%minimum=2.e-6,port_sizes(2)%maximum=4.e-5\n/\n")


if a.rates:
    assert a.multiport
    spray = spray.replace('event_mass=2.e-11,3.e-11,1.e-11,4.e-11', 'event_mass=4*0.')
    spray += ("\n&amr_spray_rates\n profile_points=2,3,\n"
              " profile_times(:,1)=0.,2.e-8,profile_rates(:,1)=.002,.004,\n"
              " profile_times(:,2)=0.,1.e-8,2.e-8,profile_rates(:,2)=.001,.002,.006,\n"
              " event_end_times=1.e-8,1.5e-8,1.5e-8,2.e-8\n/\n")


if a.discharge:
    assert a.rates and a.multiport
    spray = spray.replace('port_half_angles=.1,.2', 'port_half_angles=0.,0.')
    spray = spray.replace('port_radii=1.e-5,2.e-5', 'port_radii=1.e-4,2.e-4')
    spray = spray.replace('port_nozzles(2)%swirl=.2', 'port_nozzles(2)%swirl=0.')
    spray += "\n&amr_spray_velocity rate_velocity=2*.true.,discharge_coefficients=.8,.9 /\n"


if a.split:
    assert not (a.multiport or a.physical or a.films or a.sizes)
    spray = spray.replace('&amr_spray\n',
        "&amr_spray\n options%breakup_model='etab',options%maximum_multiplicity=5.,\n"
        " options%surface_tension=1000.,options%liquid_viscosity=.001,\n"
        " options%drag=.false.,options%heat=.false.,options%evaporation=.false.,\n", 1)


if a.rate_files:
    assert a.rates
    spray = spray.replace('profile_points=2,3,', '')
    spray = spray.replace('profile_times(:,1)=0.,2.e-8,profile_rates(:,1)=.002,.004,', '')
    spray = spray.replace('profile_times(:,2)=0.,1.e-8,2.e-8,profile_rates(:,2)=.001,.002,.006,', '')
    spray += "\n&amr_spray_rate_files profile_files='rate1.dat','rate2.dat' /\n"


def run(name, ranks, controls, good=True, changed=False, changed_boundary=False, changed_les=False, changed_pressure=False, changed_force=False, changed_properties=False, invalid_properties=False, changed_boiling=False, invalid_boiling=False, surface_cap=False, changed_sizes=False, invalid_sizes=False, changed_rates=False, invalid_rates=False, changed_discharge=False, invalid_discharge=False, changed_split=False, changed_subcycling=False, changed_mach=None, changed_reconstruction=False, coefficient=None, changed_sdc=None, changed_mol_lf=False, changed_rotation=False):
    folder = work / name
    folder.mkdir()
    local_gas = gas.replace("'outflow'", "'slip_wall'") if changed_boundary else gas
    if changed_sdc == 'iterations':
        local_gas = local_gas.replace('sdc_iterations=2', 'sdc_iterations=3')
    elif changed_sdc == 'forcing':
        local_gas = local_gas.replace('external_forcing=100.', 'external_forcing=101.')
    elif changed_sdc == 'reactor':
        local_gas = local_gas.replace("sdc_reactor='cvode'", "sdc_reactor='native'")
    elif changed_sdc == 'integrator':
        local_gas = local_gas.replace(f"time_integrator='{'mol' if a.mol else 'sdc'}'", "time_integrator='split'")
    if changed_mol_lf:
        local_gas = local_gas.replace('mol_lax_friedrichs=.true.', 'mol_lax_friedrichs=.false.')
    if coefficient is not None:
        local_gas = local_gas.replace('difmag=0.1', f'difmag={coefficient}')
    if changed_reconstruction:
        local_gas = local_gas.replace("reconstruction='characteristic_ppm'", "reconstruction='characteristic_plm'")
    if changed_pressure:
        local_gas = local_gas.replace('outlet_pressure(2)=100000.', 'outlet_pressure(2)=100100.')
    if changed_mach:
        old, new = {'mach': ('inflow_mach(1)=.2', 'inflow_mach(1)=.3'),
                    'temperature': ('inflow_temperature=1000.', 'inflow_temperature=1001.'),
                    'pressure': ('inflow_pressure=0.', 'inflow_pressure=120000.')}[changed_mach]
        local_gas = local_gas.replace(old, new)
    if changed_rotation:
        local_gas = local_gas.replace('angular_speed=200.', 'angular_speed=201.')
    (folder / 'gas.nml').write_text(local_gas)
    local_spray = spray
    if changed:
        local_spray = (spray.replace('port_temperatures=350.', 'port_temperatures=351.') if a.multiport
                       else spray.replace('speed=20.', 'speed=21.'))
    if changed_subcycling:
        local_spray = local_spray.replace('subcycling=.true.', 'subcycling=.false.')
    if changed_les:
        local_spray = local_spray.replace('les_ci=.1', 'les_ci=.2')
        if a.dynamic:
            local_spray = local_spray.replace('les_grid_ratio=1', 'les_grid_ratio=2')
    if changed_force:
        local_spray = local_spray.replace('acceleration=100.', 'acceleration=101.')
    if changed_properties:
        local_spray = local_spray.replace('1100.,-.2', '1101.,-.2')
    if invalid_properties:
        local_spray = local_spray.replace('1100.,-.2', '-1100.,0.')
    if changed_boiling:
        local_spray = local_spray.replace('373.15,450.', '374.15,450.')
    if invalid_boiling:
        local_spray = local_spray.replace('373.15,450.', '647.,450.')
    if surface_cap:
        local_spray = local_spray.replace("'pelemp-cap','pelemp-cap'", "'surface-cap','surface-cap'")
        local_spray = local_spray.replace('373.15,450.', '0.,0.')
    if changed_sizes:
        local_spray = local_spray.replace('%std=2.e-6', '%std=3.e-6')
    if invalid_sizes:
        local_spray = local_spray.replace('%std=2.e-6', '%std=-2.e-6')
    if changed_rates:
        local_spray = local_spray.replace('=.002,.004', '=.002,.005')
    if invalid_rates:
        local_spray = local_spray.replace('=.002,.004', '=-.002,.004')
    if changed_discharge:
        local_spray = local_spray.replace('discharge_coefficients=.8,.9', 'discharge_coefficients=.7,.9')
    if invalid_discharge:
        local_spray = local_spray.replace('discharge_coefficients=.8,.9', 'discharge_coefficients=0.,.9')
    if changed_split:
        local_spray = local_spray.replace('maximum_multiplicity=5.', 'maximum_multiplicity=6.')
    if a.rate_files:
        rate1 = '# time [s] rate [kg/s]\n\n0 .002\n2.e-8 .004 ! endpoint\n'
        if changed_rates:
            rate1 = rate1.replace('.004', '.005')
        if invalid_rates:
            rate1 = rate1.replace('0 .002', '0 -.002')
        (folder / 'rate1.dat').write_text(rate1)
        (folder / 'rate2.dat').write_text('0 .001\n1.e-8\t.002\n2.e-8 .006\n')
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
    if a.films:
        films = []
        for rank in range(ranks):
            with (folder / f'amr_spray.csv.films.{rank}.csv').open() as f:
                films += [[float(x) for x in r] for r in list(csv.reader(f))[1:]]
        return fields, sorted(particles), sorted(films)
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
assert f'Injection events applied: {4 if a.multiport else 3}' in (continued / 'run.log').read_text()
if staged:
    assert f"Time integrator: {'mol' if a.mol else 'sdc'}" in (continued / 'run.log').read_text()
    expected_reactor = 'cvode' if a.sdc_cvode else 'native'
    assert f'SDC reactor: {expected_reactor}' in (continued / 'run.log').read_text()
    saved_gas = (split / 'state.gas').read_text().splitlines()
    if a.mol:
        assert saved_gas[0] == 'PELEF_AMR_TREE_INTERNAL_ENERGY_V1'
        assert math.isfinite(float(saved_gas[1]))
        assert saved_gas[2] == 'PELEF_AMR_TREE_3D_CHECKPOINT_V6'
        assert 'SPRAY_ACTUAL_UEINT_V1' in saved_gas
    else:
        assert saved_gas[0] == 'PELEF_AMR_TREE_3D_CHECKPOINT_V6'
    with (continued / 'amr_spray.csv').open() as stream:
        heat = [float(row['heat_release_W_m3']) for row in csv.DictReader(stream)]
    assert heat and all(math.isfinite(value) for value in heat)
    if not (a.reacting or a.sdc_cvode):
        assert all(value == 0 for value in heat), 'Nonreacting stage heat must be zero'


if a.rotating:
    saved_rotation = (split / 'state.gas').read_text()
    assert 'CHECKPOINT_ROTATION_V1' in saved_rotation
    assert 'SPRAY_COORDINATE_ENERGY_V1' in saved_rotation
    particle_lines = (split / 'state.particles').read_text().splitlines()
    counter_index = particle_lines.index('SPRAY_COORDINATE_ENERGY_V1') + 1
    saved_counter = float(particle_lines[counter_index])
    assert math.isfinite(saved_counter) and saved_counter != 0
    def coordinate_counter(folder):
        line = next(line for line in (folder / 'run.log').read_text().splitlines()
                    if line.startswith('Coordinate energy inventory:'))
        return float(line.split(':', 1)[1])
    assert math.isclose(coordinate_counter(full), coordinate_counter(continued), rel_tol=2e-12, abs_tol=1e-25)
    run('changed_rotation', 2, f"restart_file='{checkpoint}'", good=False, changed_rotation=True)
    corrupted = work / 'corrupt_coordinate_checkpoint'
    shutil.copyfile(checkpoint, corrupted)
    shutil.copyfile(split / 'state.gas', Path(str(corrupted) + '.gas'))
    particle_lines[counter_index] = 'nan'
    Path(str(corrupted) + '.particles').write_text('\n'.join(particle_lines) + '\n')
    run('invalid_coordinate', 2, f"restart_file='{corrupted}'", good=False)

particles = read(continued, 4)[1]
assert len(particles) == (52 if a.split else 0 if a.films else 8 if a.physical else 16 if a.multiport else 12)
if a.multiport and not a.films:
    assert [row[2] for row in particles] == [1.]*8 + ([] if a.physical else [2.]*8)
    assert all(row[3] < .001 if row[2] == 1 else row[3] > .003 for row in particles)
    assert all(row[12] < row[13] if row[2] == 1 else row[12] > row[13] for row in particles)
if a.physical:
    def boundary_inventory(folder):
        line = next(line for line in (folder / 'run.log').read_text().splitlines()
                    if line.startswith('Boundary inventory:'))
        return [float(x) for x in line.split(':', 1)[1].split()]
    expected_boundary = boundary_inventory(full)
    actual_boundary = boundary_inventory(continued)
    assert expected_boundary[0] < 0
    if not a.inflow and not a.films:
        assert expected_boundary[1] > 0
    for left, right in zip(expected_boundary, actual_boundary):
        assert math.isclose(left, right, rel_tol=2e-12, abs_tol=1e-25), (left, right)
    assert all(row[6] > 0 for row in particles), 'low-wall parcels must reflect'
if a.films:
    films = read(continued, 4)[2]
    assert len(films) == 8
    assert all(row[0] == 1 and row[1] == 0 and row[4] == 1 for row in films)
    assert all(row[7] > 0 and row[9] > 0 and row[10] > 0 for row in films)
    assert sum(row[7] for row in films) < 5e-11, 'deposited films must evaporate'
    saved_lines = (split / 'state.particles').read_text().splitlines()
    assert saved_lines[0] == 'PELEF_AMR_PARTICLES_V2'
    saved_parcels = int(saved_lines[2].split()[1])
    saved_films = int(saved_lines[4 + saved_parcels])
    assert saved_films == 4, 'first injection must deposit before checkpoint'
    saved_film_mass = sum(float(line.split()[7]) for line in
                          saved_lines[5 + saved_parcels:5 + saved_parcels + saved_films])
    assert sum(sorted(row[7] for row in films)[:4]) < saved_film_mass, 'saved films must continue evaporating'
if a.adaptive:
    fields = read(continued, 4)[0]
    assert {row[1] for row in fields} == {0., 1.}
    assert len({row[0] for row in fields}) > 2
    full_log = (full / 'run.log').read_text()
    resumed_log = (continued / 'run.log').read_text()
    assert 'Adaptation check step: 1 ' in full_log
    assert 'Adaptation check step: 1 ' not in resumed_log
    assert 'Adaptation check step: 2 ' in resumed_log
if a.rebalance:
    assert 'Rebalance step: 2 ' in (continued / 'run.log').read_text()
if a.split:
    assert len({row[0] for row in particles}) == 52 and max(row[0] for row in particles) == 52
    assert sum(row[1] > 0 for row in particles) == 40
    assert all(0 < row[11] <= 5 for row in particles)
    assert math.isclose(sum(row[9]*row[11] for row in particles), 1.e-10, rel_tol=2e-12)
    assert all(math.isclose(row[12]/row[9], .25, rel_tol=2e-12) for row in particles)
    saved_header = (split / 'state.particles').read_text().splitlines()[2].split()
    assert int(saved_header[1]) == 20 and int(saved_header[2]) == 20
    run('changed_split', 2, f"restart_file='{checkpoint}'", good=False, changed_split=True)
if a.discharge:
    line = next(x for x in (continued / 'run.log').read_text().splitlines() if x.startswith('Injected inventory:'))
    impulse = float(line.split(':', 1)[1].split()[1])
    rho1, rho2 = 1/(.25/1000+.75/800), 1/(.75/1000+.25/800)
    expected = sum(m*m/(dt*rho1*math.pi*1e-8*.8) for m, dt in ((2.5e-11, 1e-8), (1.625e-11, 5e-9)))
    expected -= sum(m*m/(5e-9*rho2*math.pi*4e-8*.9) for m in (1.5e-11, 2.5e-11))
    assert math.isclose(impulse, expected, rel_tol=2e-12), (impulse, expected)
    run('changed_discharge', 2, f"restart_file='{checkpoint}'", good=False, changed_discharge=True)
    run('invalid_discharge', 2, '', good=False, invalid_discharge=True)
if a.rates:
    line = next(x for x in (continued / 'run.log').read_text().splitlines() if x.startswith('Injected inventory:'))
    injected_mass = float(line.split(':', 1)[1].split()[0])
    assert math.isclose(injected_mass, 8.125e-11, rel_tol=2e-12)
    run('changed_rates', 2, f"restart_file='{checkpoint}'", good=False, changed_rates=True)
    run('invalid_rates', 2, '', good=False, invalid_rates=True)
if a.sizes:
    for port in (1., 2.):
        drop_masses = [row[9] for row in particles if row[2] == port]
        assert max(drop_masses) > 1.1 * min(drop_masses), 'size distribution must vary drop mass'
    run('changed_sizes', 2, f"restart_file='{checkpoint}'", good=False, changed_sizes=True)
    run('invalid_sizes', 2, '', good=False, invalid_sizes=True)
if a.properties:
    run('changed_properties', 2, f"restart_file='{checkpoint}'", good=False, changed_properties=True)
    run('invalid_properties', 2, '', good=False, invalid_properties=True)
if a.pelemp_cap:
    run('changed_boiling', 2, f"restart_file='{checkpoint}'", good=False, changed_boiling=True)
    run('invalid_boiling', 2, '', good=False, invalid_boiling=True)
    surface_cap_case = run('surface_cap', 2, '', surface_cap=True)
    pelemp_particles = read(full, 2)[1]
    surface_particles = read(surface_cap_case, 2)[1]
    assert max(abs(left[9] - right[9]) for left, right in zip(pelemp_particles, surface_particles)) > 1e-17, \
        'PeleMP vapor state must change parcel mass'
if a.force:
    def force_inventory(folder):
        line = next(x for x in (folder / 'run.log').read_text().splitlines() if x.startswith('Source inventory:'))
        return [float(x) for x in line.split(':', 1)[1].split()]
    left, right = force_inventory(full), force_inventory(continued)
    assert left[1] > 0 and left[2] < 0 and left[3] > 0
    assert all(math.isclose(x, y, rel_tol=2e-12, abs_tol=1e-25) for x, y in zip(left, right))
    assert 'SPRAY_SOURCE_V1' in (split / 'state.particles').read_text()
    run('changed_force', 2, f"restart_file='{checkpoint}'", good=False, changed_force=True)
if a.pressure:
    run('changed_pressure', 2, f"restart_file='{checkpoint}'", good=False, changed_pressure=True)
if a.difmag:
    equivalent = run('auto_coefficient', 4, f"restart_file='{checkpoint}'", coefficient=-1)
    assert read(equivalent, 4) == read(continued, 4)
    rejected = run('changed_coefficient', 2, f"restart_file='{checkpoint}'", good=False, coefficient=.2)
    assert 'AMR spray checkpoint operation failed' in (rejected / 'run.log').read_text()
if a.ppm:
    run('changed_reconstruction', 2, f"restart_file='{checkpoint}'", good=False, changed_reconstruction=True)
if a.subcycling:
    run('changed_subcycling', 2, f"restart_file='{checkpoint}'", good=False, changed_subcycling=True)
if a.les:
    run('changed_les', 2, f"restart_file='{checkpoint}'", good=False, changed_les=True)
if a.mach_inlet:
    for control in ('mach', 'temperature', 'pressure'):
        run('changed_inlet_' + control, 2, f"restart_file='{checkpoint}'", good=False, changed_mach=control)
run('changed', 2, f"restart_file='{checkpoint}'", good=False, changed=True)
if a.physical:
    run('changed_boundary', 2, f"restart_file='{checkpoint}'", good=False, changed_boundary=True)
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

if staged:
    run('changed_sdc_iterations', 2, f"restart_file='{checkpoint}'", good=False, changed_sdc='iterations')
    if not (a.sdc_cvode or a.gas_force):
        run('changed_sdc_integrator', 2, f"restart_file='{checkpoint}'", good=False, changed_sdc='integrator')
if a.sdc_cvode:
    run('changed_sdc_reactor', 2, f"restart_file='{checkpoint}'", good=False, changed_sdc='reactor')
if a.gas_force:
    run('changed_gas_forcing', 2, f"restart_file='{checkpoint}'", good=False, changed_sdc='forcing')
if a.mol_lf:
    run('changed_mol_lf', 2, f"restart_file='{checkpoint}'", good=False, changed_mol_lf=True)
print('PASS: requested SDC/MOL restart controls checked')
