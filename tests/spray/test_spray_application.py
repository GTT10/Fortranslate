"""End-to-end spray CLI, conservation, injection and strict restart tests."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from statistics import NormalDist
import subprocess
import tempfile
import unittest

PARSER = argparse.ArgumentParser()
PARSER.add_argument('--exe', required=True, type=Path)
ARGS, REST = PARSER.parse_known_args()
EXE = ARGS.exe.resolve()


class SprayApplication(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix='pelef-spray-cli-')
        self.root = Path(self.directory.name)

    def tearDown(self):
        self.directory.cleanup()

    def run_case(self, name, extra='', fail=False, common=''):
        path = self.root / (name + '.nml')
        path.write_text(f'''&spray
 cells=2,2,2, hydro_reconstruction='pcm', chemistry_enabled=.false.,
 final_time=1e-6, maximum_dt=1e-7, parcel_count=2,
 output_prefix='{self.root / name}',
 {common}
 {extra}
/\n''')
        result = subprocess.run([str(EXE), str(path)], cwd=self.root,
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
        (self.root / (name+'.log')).write_text(result.stdout)
        if fail:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        else:
            self.assertEqual(result.returncode, 0, result.stdout)
        return result

    def history(self, name):
        with (self.root / (name+'.history.csv')).open() as stream:
            return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]

    def check_balance(self, name):
        rows = self.history(name)
        self.assertGreater(len(rows), 1)
        for row in rows:
            self.assertTrue(all(math.isfinite(x) for x in row.values()))
            self.assertLess(abs(row['mass_error']), 2e-12*row['system_mass'])
            self.assertLess(abs(row['energy_error']), 2e-12*max(abs(row['total_energy']), 1e-12))
            self.assertLess(max(abs(row[k]) for k in ['px_error','py_error','pz_error']), 1e-18)
            self.assertGreater(row['tmin'], 0)
        self.assertEqual(rows[-1]['time'], 1e-6)
        return rows

    def test_pelephysics_drag_restart(self):
        settings = ("drag_model='pelephysics',breakup_model='etab',surface_tension=.072,liquid_viscosity=.001,"
                    "injection_speed=100,heat_enabled=.false.,evaporation_enabled=.false.,")
        checkpoint = self.root/'reference_drag.chk'
        self.run_case('reference_drag', common=settings)
        self.check_balance('reference_drag')
        self.run_case('reference_drag_stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
        self.run_case('reference_drag_resumed', f"restart_file='{checkpoint}',", common=settings)
        for suffix in ['gas.csv', 'parcels.csv']:
            self.assertEqual((self.root/f'reference_drag.{suffix}').read_bytes(),
                             (self.root/f'reference_drag_resumed.{suffix}').read_bytes())
        self.assertEqual(self.history('reference_drag')[-1], self.history('reference_drag_resumed')[-1])
        result = self.run_case('drag_changed', f"restart_file='{checkpoint}',", fail=True,
                               common=settings.replace("drag_model='pelephysics'", "drag_model='schiller-naumann'"))
        self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('drag_unknown', common="drag_model='unknown',", fail=True)

    def test_film_evolution_restart(self):
        settings = ("boundary_faces='periodic','periodic','periodic','periodic','no_slip_wall','no_slip_wall',"
                    "wall_impact_faces=.false.,.false.,.false.,.false.,.true.,.false.,"
                    "wall_film_evolution=.true.,wall_liquid_conductivity=.6,wall_temperature=340,340,340,340,340,340,"
                    "surface_tension=.072,liquid_viscosity=.001,injection_speed=10,"
                    "injection_center=.002,.002,5e-7,injection_axis=0,0,-1,cone_angle=0,parcel_count=1,"
                    "mole_names='N2','','',mole_amounts=1,0,0,gas_temperature=900,"
                    "drag_enabled=.false.,transport_enabled=.false.,")
        checkpoint = self.root/'evolving.chk'
        self.run_case('evolving', common=settings)
        history = self.check_balance('evolving')
        self.assertGreater(history[-1]['film_mass'], 0)
        self.assertLess(history[-1]['film_mass'], max(r['film_mass'] for r in history))
        with (self.root/'evolving.gas.csv').open() as stream:
            self.assertGreater(sum(float(r['rhoY_H2O']) for r in csv.DictReader(stream)), 0)
        self.run_case('evolving_stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
        self.run_case('evolving_resumed', f"restart_file='{checkpoint}',", common=settings)
        for suffix in ['gas.csv', 'parcels.csv', 'films.csv']:
            self.assertEqual((self.root/f'evolving.{suffix}').read_bytes(),
                             (self.root/f'evolving_resumed.{suffix}').read_bytes())
        self.assertEqual(history[-1], self.history('evolving_resumed')[-1])
        old_checkpoint = self.root/'evolving_old_skin.chk'
        old_checkpoint.write_bytes(checkpoint.read_bytes().replace(
            b'wall-film-evolution-v2-half-skin', b'wall-film-evolution-v1'))
        result = self.run_case('evolving_old_skin', f"restart_file='{old_checkpoint}',", fail=True, common=settings)
        self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        result = self.run_case('evolving_changed', f"restart_file='{checkpoint}',", fail=True,
                               common=settings.replace('conductivity=.6', 'conductivity=.7'))
        self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('evolving_bad_lambda', common=settings.replace('conductivity=.6', 'conductivity=0'), fail=True)
        vanishing = settings+'minimum_film_height=1e-3,'
        empty_checkpoint = self.root/'vanishing.chk'
        self.run_case('vanishing', common=vanishing)
        self.assertEqual(self.check_balance('vanishing')[-1]['film_mass'], 0)
        self.run_case('vanishing_stop', f"stop_after_steps=3,checkpoint_file='{empty_checkpoint}',", common=vanishing)
        self.run_case('vanishing_resumed', f"restart_file='{empty_checkpoint}',", common=vanishing)
        for suffix in ['gas.csv', 'parcels.csv', 'films.csv']:
            self.assertEqual((self.root/f'vanishing.{suffix}').read_bytes(),
                             (self.root/f'vanishing_resumed.{suffix}').read_bytes())

    def test_film_conductivity_models(self):
        common = ("boundary_faces='periodic','periodic','periodic','periodic','no_slip_wall','no_slip_wall',"
                  "wall_impact_faces=.false.,.false.,.false.,.false.,.true.,.false.,"
                  "wall_film_evolution=.true.,wall_liquid_conductivity=.6,wall_temperature=340,340,340,340,340,340,"
                  "surface_tension=.072,liquid_viscosity=.001,injection_speed=10,cone_angle=0,parcel_count=1,"
                  "mole_names='N2','','',mole_amounts=1,0,0,gas_temperature=900,"
                  "drag_enabled=.false.,transport_enabled=.false.,evaporation_enabled=.false.,")
        cases = [common+"injection_center=.002,.002,5e-7,injection_axis=0,0,-1,"
                 "wall_conductivity_model='polynomial',wall_conductivity_coefficients=.1,.001,0,0,"
                 "wall_conductivity_temperature_range=300,600,",
                 common+"injection_port_count=2,injection_port_centers=.001,.001,5e-7,.003,.001,5e-7,"
                 "injection_port_axes=0,0,-1,0,0,-1,injection_port_fractions=.5,.5,"
                 "port_wall_conductivity_models='polynomial','constant',port_wall_conductivities(2)=.1,"
                 "port_wall_conductivity_coefficients(:,1)=.1,.001,0,0,"
                 "port_wall_conductivity_temperature_ranges(:,1)=300,600,"]
        for index, settings in enumerate(cases):
            name = f'lambda{index}'
            checkpoint = self.root/f'{name}.chk'
            self.run_case(name, common=settings)
            self.check_balance(name)
            self.run_case(name+'stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
            self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=settings)
            for suffix in ['gas.csv', 'parcels.csv', 'films.csv']:
                self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                 (self.root/f'{name}resumed.{suffix}').read_bytes())
            self.assertEqual(self.history(name)[-1], self.history(name+'resumed')[-1])
            result = self.run_case(name+'changed', f"restart_file='{checkpoint}',", fail=True,
                                   common=settings.replace('=.1,.001,0,0', '=.2,.001,0,0'))
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        with (self.root/'lambda1.films.csv').open() as stream:
            films = {int(r['material_id']): r for r in csv.DictReader(stream)}
        self.assertEqual(set(films), {1, 2})
        self.assertGreater(float(films[1]['temperature']), float(films[2]['temperature']))
        self.run_case('lambda_outside', common=cases[0].replace('300,600', '400,600'), fail=True)
        self.run_case('lambda_negative_interior', common=cases[0].replace('=.1,.001,0,0', '=15.5,-.08,.0001,0'),
                      fail=True)
        self.run_case('lambda_unused_port', common=cases[0]+'port_wall_conductivities(3)=.5,', fail=True)

    def test_wall_film_restart(self):
        common = ("boundary_faces='periodic','periodic','periodic','periodic','no_slip_wall','no_slip_wall',"
                  "wall_impact_faces=.false.,.false.,.false.,.false.,.true.,.false.,"
                  "surface_tension=.072,liquid_viscosity=.001,"
                  "injection_center=.002,.002,5e-7,injection_axis=0,0,-1,cone_angle=0,parcel_count=1,"
                  "mole_names='N2','','',mole_amounts=1,0,0,gas_temperature=900,"
                  "drag_enabled=.false.,heat_enabled=.false.,evaporation_enabled=.false.,transport_enabled=.false.,")
        for speed in [10, 300]:
            name = f'wall{speed}'
            settings = common+f'injection_speed={speed},'
            checkpoint = self.root/f'{name}.chk'
            self.run_case(name, common=settings)
            history = self.check_balance(name)
            self.assertGreater(history[-1]['film_mass'], 0)
            with (self.root/f'{name}.films.csv').open() as stream:
                films = list(csv.DictReader(stream))
            self.assertEqual(len(films), 1)
            self.assertAlmostEqual(sum(float(f['mass']) for f in films)/history[-1]['film_mass'], 1, delta=1e-13)
            with (self.root/f'{name}.parcels.csv').open() as stream:
                parcels = list(csv.DictReader(stream))
            self.assertEqual(len(parcels), 1 if speed == 10 else 5)
            if speed == 300:
                self.assertEqual([int(p['parent_id']) for p in parcels], [0, 1, 1, 1, 1])
            self.run_case(name+'stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
            self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=settings)
            for suffix in ['gas.csv', 'parcels.csv', 'films.csv']:
                self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                 (self.root/f'{name}resumed.{suffix}').read_bytes())
            self.assertEqual(history[-1], self.history(name+'resumed')[-1])
            result = self.run_case(name+'changed', f"restart_file='{checkpoint}',wall_contact_angle=1.1,",
                                   common=settings, fail=True)
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
            lines = checkpoint.read_text().splitlines()
            values = lines[-2].split()
            values[-1] = '-1'
            lines[-2] = ' '.join(values)
            bad = self.root/f'{name}_bad.chk'
            bad.write_text('\n'.join(lines)+'\n')
            self.run_case(name+'bad', f"restart_file='{bad}',", common=settings, fail=True)
        self.run_case('impact_bad_face', common='wall_impact_faces=.true.,', fail=True)
        self.run_case('impact_bad_sigma', common=common.replace('surface_tension=.072', 'surface_tension=0'), fail=True)

    def test_heat_capacity_models_restart(self):
        cases = ["liquid_cp_model='polynomial',liquid_cp_coefficients=4000,4,0,0,liquid_cp_temperature_range=300,600,",
                 "injection_port_count=2,injection_port_centers=.001,.001,.001,.002,.001,.001,"
                 "injection_port_axes=1,0,0,1,0,0,injection_port_fractions=.5,.5,"
                 "port_cp_models='polynomial','constant',port_cp_coefficients(:,1)=4000,4,0,0,"
                 "port_cp_temperature_ranges(:,1)=300,600,"]
        for index, settings in enumerate(cases):
            name = f'cp{index}'
            checkpoint = self.root/f'{name}.chk'
            self.run_case(name, common=settings)
            self.check_balance(name)
            self.run_case(name+'stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
            self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=settings)
            for suffix in ['gas.csv', 'parcels.csv']:
                self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                 (self.root/f'{name}resumed.{suffix}').read_bytes())
            self.assertEqual(self.history(name)[-1], self.history(name+'resumed')[-1])
            result = self.run_case(name+'changed', f"restart_file='{checkpoint}',", fail=True,
                                   common=settings.replace('4000,4', '4100,4'))
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('cp_bad_range', common=cases[0].replace('300,600', '350,600'), fail=True)
        self.run_case('cp_bad_shape', common=cases[0].replace('4000,4,0,0', '202499,-900,1,0'), fail=True)
        self.run_case('cp_unused', common='port_cp_coefficients(1,1)=4000,', fail=True)

    def test_port_breakup_properties_restart(self):
        common = ("injection_port_count=2,injection_port_centers=.001,.001,.001,.002,.001,.001,"
                  "injection_port_axes=1,0,0,1,0,0,injection_port_fractions=.5,.5,"
                  "injection_speed=500,cone_angle=0,gas_pressure=1e7,"
                  "surface_tension=.072,liquid_viscosity=.001,"
                  "port_surface_tensions=.072,.035,port_liquid_viscosities=.001,.002,"
                  "port_viscosity_models='constant','inverse-temperature',"
                  "port_viscosity_coefficients(:,2)=.001,.33,0,0,"
                  "port_viscosity_temperature_ranges(:,2)=300,400,"
                  "drag_enabled=.false.,heat_enabled=.false.,evaporation_enabled=.false.,transport_enabled=.false.,")
        for model in ['etab', 'khrt']:
            settings = common+f"breakup_model='{model}',"
            name = 'port_'+model
            checkpoint = self.root/f'{name}.chk'
            self.run_case(name, common=settings)
            self.check_balance(name)
            with (self.root/f'{name}.parcels.csv').open() as stream:
                parcels = list(csv.DictReader(stream))
            self.assertEqual([int(p['material_id']) for p in parcels], [1, 1, 2, 2])
            # Different properties must alter the actual breakup state, not only metadata.
            self.assertNotEqual((parcels[0]['diameter'], parcels[0]['deformation']),
                                (parcels[2]['diameter'], parcels[2]['deformation']))
            self.run_case(name+'stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
            self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=settings)
            for suffix in ['gas.csv', 'parcels.csv']:
                self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                 (self.root/f'{name}resumed.{suffix}').read_bytes())
            self.assertEqual(self.history(name)[-1], self.history(name+'resumed')[-1])
            result = self.run_case(name+'changed', f"restart_file='{checkpoint}',", fail=True,
                                   common=settings.replace('.072,.035', '.072,.04'))
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('bad_port_sigma', common=common+"breakup_model='etab',port_surface_tensions=.072,0,", fail=True)
        self.run_case('bad_port_mu_range', common=(common+"breakup_model='etab',").replace('300,400', '350,400'),
                      fail=True)
        self.run_case('unused_port_mu', common='port_liquid_viscosities=.001,', fail=True)

    def test_port_property_models_restart(self):
        settings = ("injection_port_count=2,injection_port_centers=.001,.001,.001,.002,.001,.001,"
                    "injection_port_axes=1,0,0,-1,0,0,injection_port_fractions=.5,.5,"
                    "injection_port_temperatures=290,400,"
                    "port_density_models='polynomial','constant',"
                    "port_density_coefficients(:,1)=1150,-.5,0,0,port_density_temperature_ranges(:,1)=280,500,"
                    "port_vapor_pressure_models='antoine','clausius',"
                    "port_antoine_coefficients(:,1)=5.40221,1838.675,-31.737,1e5,"
                    "port_antoine_temperature_ranges(:,1)=273,303,"
                    "port_boiling_models='reject','surface-cap',port_critical_temperatures=0,647,"
                    "mass_flow_rate=1e-4,injection_start=1e-7,injection_end=4e-7,injection_interval=2e-7,")
        checkpoint = self.root/'fits.chk'
        self.run_case('fits', common=settings)
        self.check_balance('fits')
        with (self.root/'fits.parcels.csv').open() as stream:
            parcels = list(csv.DictReader(stream))
        self.assertEqual({int(p['material_id']) for p in parcels}, {1, 2})
        for parcel in parcels:
            density = 1150-.5*float(parcel['temperature']) if int(parcel['material_id']) == 1 else 1000
            expected_mass = density*math.pi/6*float(parcel['diameter'])**3
            self.assertAlmostEqual(float(parcel['mass'])/expected_mass, 1, delta=1e-12)
        self.run_case('fits_stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
        self.run_case('fits_resumed', f"restart_file='{checkpoint}',", common=settings)
        for suffix in ['gas.csv', 'parcels.csv']:
            self.assertEqual((self.root/f'fits.{suffix}').read_bytes(),
                             (self.root/f'fits_resumed.{suffix}').read_bytes())
        self.assertEqual(self.history('fits')[-1], self.history('fits_resumed')[-1])
        for index, (old, new) in enumerate([('1150,-.5', '1151,-.5'), ('5.40221', '5.40222'), ('0,647', '0,648')]):
            result = self.run_case(f'fits_changed{index}', f"restart_file='{checkpoint}',", fail=True,
                                   common=settings.replace(old, new))
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('fits_bad_range', common=settings.replace('273,303', '273,280'), fail=True)
        self.run_case('fits_bad_cap', common=settings.replace("'reject','surface-cap'", "'surface-cap','surface-cap'"),
                      fail=True)
        self.run_case('fits_unused', common='port_density_coefficients(1,1)=1000,', fail=True)

    def test_port_materials_restart(self):
        settings = ("injection_port_count=2,injection_port_centers=.001,.001,.001,.002,.001,.001,"
                    "injection_port_axes=1,0,0,-1,0,0,injection_port_fractions=.25,.75,"
                    "port_vapor_names='H2O','H2O2',port_liquid_densities=1000,1400,"
                    "port_liquid_cps=4184,2600,port_reference_temperatures=373.15,420,"
                    "port_latent_heats=2.257e6,1.5e6,")
        schedules = ["mass_flow_rate=1e-4,injection_start=1e-7,injection_end=4e-7,injection_interval=2e-7,",
                     "independent_port_schedules=.true.,port_mass_flow_rates=1e-4,2e-4,"
                     "port_injection_starts=1e-7,1e-7,port_injection_ends=5e-7,6e-7,"
                     "port_injection_intervals=2e-7,3e-7,"]
        for index, schedule in enumerate(schedules):
            with self.subTest(schedule=index):
                name = f'materials{index}'
                checkpoint = self.root/f'{name}.chk'
                common = settings+schedule
                self.run_case(name, common=common)
                self.check_balance(name)
                with (self.root/f'{name}.parcels.csv').open() as stream:
                    parcels = list(csv.DictReader(stream))
                self.assertEqual({int(p['material_id']) for p in parcels}, {1, 2})
                for parcel in parcels:
                    density = 1000 if int(parcel['material_id']) == 1 else 1400
                    expected_mass = density*math.pi/6*float(parcel['diameter'])**3
                    self.assertAlmostEqual(float(parcel['mass'])/expected_mass, 1, delta=1e-12)
                with (self.root/f'{name}.gas.csv').open() as stream:
                    gas = list(csv.DictReader(stream))
                for species in ['H2O', 'H2O2']:
                    self.assertGreater(sum(float(g[f'rhoY_{species}']) for g in gas), 0)
                self.run_case(name+'stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=common)
                self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=common)
                for suffix in ['gas.csv', 'parcels.csv']:
                    self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                     (self.root/f'{name}resumed.{suffix}').read_bytes())
                self.assertEqual(self.history(name)[-1], self.history(name+'resumed')[-1])
                result = self.run_case(name+'changed', f"restart_file='{checkpoint}',", fail=True,
                                       common=common.replace('cps=4184,2600', 'cps=4184,2700'))
                self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('invalid_material', common=settings.replace("'H2O2'", "'unknown'"), fail=True)
        self.run_case('invalid_material_density', common=settings.replace('1000,1400', '1000,0'), fail=True)
        self.run_case('unused_material', common="port_vapor_names='H2O',", fail=True)

    def test_port_state_overrides_restart(self):
        settings = ("injection_port_count=2,injection_port_centers=.001,.001,.001,.002,.001,.001,"
                    "injection_port_axes=1,0,0,-1,0,0,injection_port_fractions=.25,.75,"
                    "injection_port_speeds=10,20,injection_port_diameters=1e-5,2e-5,"
                    "injection_port_temperatures=300,350,injection_port_cone_angles=0,0,"
                    "mass_flow_rate=1e-4,injection_start=1e-7,injection_end=4e-7,injection_interval=2e-7,"
                    "drag_enabled=.false.,heat_enabled=.false.,evaporation_enabled=.false.,")
        checkpoint = self.root/'port-states.chk'
        self.run_case('states', common=settings)
        rows = self.check_balance('states')
        self.assertEqual(rows[-1]['live_parcels'], 12)
        with (self.root/'states.parcels.csv').open() as stream:
            parcels = list(csv.DictReader(stream))
        for index, parcel in enumerate(parcels):
            first = index % 4 < 2
            self.assertEqual(float(parcel['u']), 10 if first else -20)
            self.assertEqual(float(parcel['temperature']), 300 if first else 350)
            self.assertAlmostEqual(float(parcel['diameter'])/(1e-5 if first else 2e-5), 1, delta=1e-13)
        self.run_case('stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
        self.run_case('resumed', f"restart_file='{checkpoint}',", common=settings)
        for suffix in ['gas.csv', 'parcels.csv']:
            self.assertEqual((self.root/f'states.{suffix}').read_bytes(),
                             (self.root/f'resumed.{suffix}').read_bytes())
        result = self.run_case('changed', f"restart_file='{checkpoint}',",
                               common=settings.replace('speeds=10,20', 'speeds=10,21'), fail=True)
        self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('invalid_diameter', common=settings.replace('diameters=1e-5,2e-5', 'diameters=1e-5,0'), fail=True)

    def test_independent_port_schedule_restart(self):
        settings = ("injection_port_count=2,injection_port_centers=.001,.001,.001,.002,.001,.001,"
                    "injection_port_axes=1,0,0,-1,0,0,injection_port_fractions=.5,.5,"
                    "injection_port_speeds=10,20,independent_port_schedules=.true.,port_mass_flow_rates=1e-4,2e-4,"
                    "port_injection_starts=1e-7,1e-7,port_injection_ends=5e-7,6e-7,"
                    "port_injection_intervals=2e-7,3e-7,initial_liquid_mass=0,cone_angle=0,"
                    "drag_enabled=.false.,heat_enabled=.false.,evaporation_enabled=.false.,")
        self.run_case('independent', common=settings)
        rows = self.check_balance('independent')
        self.assertEqual(rows[-1]['live_parcels'], 8)
        self.assertAlmostEqual(rows[-1]['injected_mass']/1.4e-10, 1., delta=1e-13)
        with (self.root/'independent.parcels.csv').open() as stream:
            parcels = list(csv.DictReader(stream))
        self.assertEqual([int(p['id']) for p in parcels], list(range(1,9)))
        self.assertEqual([float(p['u']) for p in parcels], [10,10,-20,-20,10,10,-20,-20])
        for stop in [1,3]:
            checkpoint = self.root/f'independent{stop}.chk'
            self.run_case(f'stop{stop}', f"stop_after_steps={stop},checkpoint_file='{checkpoint}',", common=settings)
            self.run_case(f'resumed{stop}', f"restart_file='{checkpoint}',", common=settings)
            for suffix in ['gas.csv', 'parcels.csv']:
                self.assertEqual((self.root/f'independent.{suffix}').read_bytes(),
                                 (self.root/f'resumed{stop}.{suffix}').read_bytes())
        result = self.run_case('changed', f"restart_file='{checkpoint}',",
                               common=settings.replace('1e-4,2e-4', '1e-4,3e-4'), fail=True)
        self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('bad_interval', common=settings.replace('2e-7,3e-7', '0,3e-7'), fail=True)
        self.run_case('mixed_schedule', common=settings+'mass_flow_rate=1e-4,', fail=True)

    def test_multiport_schedule_restart(self):
        settings = ("domain_length=.001,.001,.001,injection_port_count=2,"
                    "injection_port_centers=.0002,.0005,.0005,.0008,.0005,.0005,"
                    "injection_port_axes=1,0,0,-1,0,0,injection_port_fractions=.25,.75,"
                    "mass_flow_rate=1e-4,injection_start=1e-7,injection_end=4.5e-7,injection_interval=1.5e-7,")
        checkpoint = self.root/'ports.chk'
        self.run_case('ports', common=settings)
        rows = self.check_balance('ports')
        self.assertEqual(rows[-1]['live_parcels'], 16)
        self.assertAlmostEqual(rows[-1]['injected_mass']/3.5e-11, 1., delta=1e-13)
        self.run_case('stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
        self.run_case('resumed', f"restart_file='{checkpoint}',", common=settings)
        for suffix in ['gas.csv', 'parcels.csv']:
            self.assertEqual((self.root/f'ports.{suffix}').read_bytes(),
                             (self.root/f'resumed.{suffix}').read_bytes())
        result = self.run_case('changed', f"restart_file='{checkpoint}',",
                               common=settings.replace('.25,.75', '.5,.5'), fail=True)
        self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('bad_shares', common=settings.replace('.25,.75', '.25,.5'), fail=True)
        self.run_case('storage', common=settings+'maximum_parcels=3,', fail=True)

    def test_boiling_surface_cap_restart(self):
        settings = "boiling_model='surface-cap',critical_temperature=647,drop_temperature=400,"
        checkpoint = self.root/'surface-cap.chk'
        self.run_case('boiling', common=settings)
        rows = self.check_balance('boiling')
        self.assertLess(rows[-1]['liquid_mass'], rows[0]['liquid_mass'])
        self.run_case('stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
        self.run_case('resumed', f"restart_file='{checkpoint}',", common=settings)
        for suffix in ['gas.csv', 'parcels.csv']:
            self.assertEqual((self.root/f'boiling.{suffix}').read_bytes(),
                             (self.root/f'resumed.{suffix}').read_bytes())
        result = self.run_case('changed', f"restart_file='{checkpoint}',",
                               common=settings.replace('647', '648'), fail=True)
        self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('default_boiling', common='drop_temperature=400,', fail=True)

    def test_viscosity_temperature_restart(self):
        settings = ("breakup_model='khrt',surface_tension=.072,"
                    "liquid_viscosity_model='inverse-temperature',liquid_viscosity_coefficients=.0002,.2,10,100,"
                    "liquid_viscosity_temperature_range=250,600,injection_speed=200,")
        checkpoint = self.root/'viscosity.chk'
        self.run_case('viscosity', common=settings)
        self.check_balance('viscosity')
        self.run_case('stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
        self.run_case('resumed', f"restart_file='{checkpoint}',", common=settings)
        for suffix in ['gas.csv', 'parcels.csv']:
            self.assertEqual((self.root/f'viscosity.{suffix}').read_bytes(),
                             (self.root/f'resumed.{suffix}').read_bytes())
        result = self.run_case('changed', f"restart_file='{checkpoint}',",
                               common=settings.replace('.0002', '.0003'), fail=True)
        self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('invalid_viscosity', common=settings.replace('250,600', '250,300'), fail=True)

    def test_density_polynomial_restart(self):
        settings = ("liquid_density_model='polynomial',liquid_density_coefficients=1200,-.5,0,0,"
                    "liquid_density_temperature_range=250,600,")
        checkpoint = self.root/'density.chk'
        self.run_case('density', common=settings)
        rows = self.check_balance('density')
        self.assertLess(rows[-1]['liquid_mass'], rows[0]['liquid_mass'])
        self.run_case('stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
        self.run_case('restart', f"restart_file='{checkpoint}',", common=settings)
        for suffix in ['gas.csv', 'parcels.csv']:
            self.assertEqual((self.root/f'density.{suffix}').read_bytes(),
                             (self.root/f'restart.{suffix}').read_bytes())
        result = self.run_case('changed', f"restart_file='{checkpoint}',",
                               common=settings.replace('1200', '1201'), fail=True)
        self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('invalid_density', common=settings.replace('1200', '1'), fail=True)

    def test_antoine_evaporation_restart_and_coefficients(self):
        settings = ("vapor_pressure_model='antoine',antoine_coefficients=5.40221,1838.675,-31.737,1e5,"
                    "antoine_temperature_range=273,303,drop_temperature=290,")
        checkpoint = self.root/'antoine.chk'
        self.run_case('antoine', common=settings)
        rows = self.check_balance('antoine')
        self.assertLess(rows[-1]['liquid_mass'], rows[0]['liquid_mass'])
        self.run_case('stopped', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
        self.run_case('resumed', f"restart_file='{checkpoint}',", common=settings)
        for suffix in ['gas.csv', 'parcels.csv']:
            self.assertEqual((self.root/f'antoine.{suffix}').read_bytes(),
                             (self.root/f'resumed.{suffix}').read_bytes())
        for field, changed in [('5.40221', '5.40222'), ('273,303', '272,303')]:
            result = self.run_case('mismatch'+field[:3], f"restart_file='{checkpoint}',",
                                   common=settings.replace(field, changed), fail=True)
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('range', common=settings.replace('drop_temperature=290', 'drop_temperature=310'), fail=True)

    def test_evaporation_feedback_and_conservation(self):
        self.run_case('evaporate')
        rows = self.check_balance('evaporate')
        self.assertLess(rows[-1]['liquid_mass'], rows[0]['liquid_mass'])
        self.assertGreater(rows[-1]['gas_mass'], rows[0]['gas_mass'])

    def test_scheduled_mass_and_exact_restart(self):
        schedule = 'mass_flow_rate=1e-4,injection_start=1e-7,injection_end=4.5e-7,injection_interval=1.5e-7,'
        checkpoint = self.root / 'split.chk'
        self.run_case('full', common=schedule)
        self.run_case('first', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=schedule)
        self.run_case('second', f"restart_file='{checkpoint}',", common=schedule)
        for suffix in ['gas.csv', 'parcels.csv']:
            self.assertEqual((self.root/f'full.{suffix}').read_bytes(), (self.root/f'second.{suffix}').read_bytes())
        whole = self.history('full'); tail = self.history('second')
        self.assertEqual([r for r in whole if r['steps']>=tail[0]['steps']], tail)
        rows = self.check_balance('full')
        self.assertAlmostEqual(rows[-1]['injected_mass']/3.5e-11, 1., delta=1e-13)
        self.assertEqual(rows[-1]['live_parcels'], 8)

    def test_size_distribution_and_restart(self):
        settings='parcel_count=16,diameter_shape=3,diameter_min=1e-5,diameter_max=4e-5,'
        checkpoint=self.root/'sizes.chk'
        self.run_case('sizes_full',common=settings)
        self.run_case('sizes_first',f"stop_after_steps=3,checkpoint_file='{checkpoint}',",common=settings)
        self.run_case('sizes_after',f"restart_file='{checkpoint}',",common=settings)
        self.check_balance('sizes_full')
        for suffix in ['gas.csv','parcels.csv']:
            self.assertEqual((self.root/f'sizes_full.{suffix}').read_bytes(),
                             (self.root/f'sizes_after.{suffix}').read_bytes())
        with (self.root/'sizes_full.parcels.csv').open() as stream:
            data=list(csv.DictReader(stream))
        self.assertEqual(len({r['diameter'] for r in data}),16)

    def test_rate_derived_injection_velocity(self):
        base = ("initial_liquid_mass=0,velocity_from_rate=.true.,discharge_coefficient=.5,nozzle_radius=1e-4,"
                "injection_axis=1,0,0,cone_angle=0,injection_start=0,injection_end=8e-7,injection_interval=1e-7,"
                "drag_enabled=.false.,heat_enabled=.false.,evaporation_enabled=.false.,")
        shared = ("injection_profile_points=2,injection_profile_times=0,8e-7,injection_profile_rates=1e-4,2e-4,")
        ports = ("injection_port_count=2,injection_port_centers=.001,.002,.002,.003,.002,.002,"
                 "injection_port_axes=1,0,0,-1,0,0,injection_port_fractions=.3,.7,"
                 "port_nozzle_radii=1e-4,2e-4,port_discharge_coefficients=.5,.8,port_liquid_densities=1000,800,")
        independent = (ports+"independent_port_schedules=.true.,port_injection_starts=0,0,"
                       "port_injection_ends=8e-7,8e-7,port_injection_intervals=1e-7,2e-7,"
                       "port_profile_points=2,0,port_profile_times(:,1)=0,8e-7,"
                       "port_profile_rates(:,1)=1e-4,2e-4,port_mass_flow_rates=0,1e-4,")
        for index, schedule in enumerate([shared, shared+ports, independent]):
            name = f'rate_velocity{index}'
            settings = base+schedule
            checkpoint = self.root/f'{name}.chk'
            self.run_case(name, common=settings)
            self.check_balance(name)
            with (self.root/f'{name}.parcels.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            speeds = []
            for row in rows:
                port2 = index > 0 and int(row['material_id']) == 2
                radius, density, cd = (2e-4, 800, .8) if port2 else (1e-4, 1000, .5)
                duration = 2e-7 if index == 2 and port2 else 1e-7
                rate = 2*float(row['mass'])*float(row['multiplicity'])/duration
                expected = rate/(density*math.pi*radius**2*cd)
                speed = math.sqrt(sum(float(row[k])**2 for k in ['u','v','w']))
                self.assertAlmostEqual(speed/expected, 1, delta=2e-13)
                speeds.append(speed)
            self.assertGreater(max(speeds), min(speeds))
            self.run_case(name+'stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
            self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=settings)
            for suffix in ['gas.csv', 'parcels.csv']:
                self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                 (self.root/f'{name}resumed.{suffix}').read_bytes())
            self.assertEqual(self.history(name)[-1], self.history(name+'resumed')[-1])
            result = self.run_case(name+'changed', f"restart_file='{checkpoint}',", fail=True,
                                   common=settings.replace('discharge_coefficient=.5', 'discharge_coefficient=.6'))
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('rate_velocity_point', common=base+shared+'nozzle_radius=0,', fail=True)
        self.run_case('rate_velocity_cd', common=base+shared+'discharge_coefficient=0,', fail=True)
        self.run_case('rate_velocity_initial', common=base+shared+'initial_liquid_mass=1e-10,', fail=True)

    def test_independent_port_rate_profiles(self):
        base = ("initial_liquid_mass=0,injection_port_count=2,independent_port_schedules=.true.,"
                "injection_port_centers=.001,.001,.001,.003,.001,.001,injection_port_axes=1,0,0,-1,0,0,"
                "injection_port_fractions=.5,.5,injection_speed=0,port_injection_starts=0,0,"
                "port_injection_ends=8e-7,8e-7,port_injection_intervals=1e-7,2e-7,"
                "port_profile_points=5,0,port_profile_times(:,1)=0,2e-7,4e-7,6e-7,8e-7,"
                "port_profile_rates(:,1)=0,1e-4,0,0,2e-4,"
                "drag_enabled=.false.,heat_enabled=.false.,evaporation_enabled=.false.,")
        variants = [base+'port_mass_flow_rates=0,1e-4,',
                    base+"port_profile_points(2)=3,port_profile_times(:,2)=0,4e-7,8e-7,"
                         "port_profile_rates(:,2)=0,2e-4,0,"]
        for index, settings in enumerate(variants):
            name = f'port_rates{index}'
            checkpoint = self.root/f'{name}.chk'
            self.run_case(name, common=settings)
            history = self.check_balance(name)
            with (self.root/f'{name}.parcels.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 20)
            masses = [sum(float(r['mass'])*float(r['multiplicity']) for r in rows
                          if (float(r['x']) < .002) == (port == 0)) for port in range(2)]
            self.assertAlmostEqual(masses[0]/4e-11, 1, delta=2e-13)
            self.assertAlmostEqual(masses[1]/8e-11, 1, delta=2e-13)
            self.run_case(name+'stop', f"stop_after_steps=6,checkpoint_file='{checkpoint}',", common=settings)
            self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=settings)
            for suffix in ['gas.csv', 'parcels.csv']:
                self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                 (self.root/f'{name}resumed.{suffix}').read_bytes())
            self.assertEqual(history[-1], self.history(name+'resumed')[-1])
            result = self.run_case(name+'changed', f"restart_file='{checkpoint}',", fail=True,
                                   common=settings.replace('rates(:,1)=0,1e-4', 'rates(:,1)=0,1.1e-4'))
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('port_rate_conflict', common=variants[0]+'port_mass_flow_rates(1)=1e-4,', fail=True)
        self.run_case('port_rate_unused', common=variants[0]+'port_profile_points(3)=2,', fail=True)
        self.run_case('port_rate_bad_times', common=variants[0].replace('0,2e-7,4e-7', '0,0,4e-7'), fail=True)

    def test_piecewise_injection_profile(self):
        base = ("initial_liquid_mass=0,injection_profile_points=5,"
                "injection_profile_times=0,2e-7,4e-7,6e-7,8e-7,injection_profile_rates=0,1e-4,0,0,2e-4,"
                "injection_start=0,injection_end=8e-7,injection_interval=1e-7,"
                "drag_enabled=.false.,heat_enabled=.false.,evaporation_enabled=.false.,")
        ports = ("injection_port_count=2,injection_port_centers=.001,.001,.001,.003,.001,.001,"
                 "injection_port_axes=1,0,0,-1,0,0,injection_port_fractions=.25,.75,")
        for index, settings in enumerate([base, base+ports]):
            name = f'rate_profile{index}'
            checkpoint = self.root/f'{name}.chk'
            self.run_case(name, common=settings)
            history = self.check_balance(name)
            with (self.root/f'{name}.parcels.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 12*(index+1))
            mass = sum(float(r['mass'])*float(r['multiplicity']) for r in rows)
            self.assertAlmostEqual(mass/4e-11, 1, delta=2e-13)
            self.run_case(name+'stop', f"stop_after_steps=6,checkpoint_file='{checkpoint}',", common=settings)
            self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=settings)
            for suffix in ['gas.csv', 'parcels.csv']:
                self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                 (self.root/f'{name}resumed.{suffix}').read_bytes())
            self.assertEqual(history[-1], self.history(name+'resumed')[-1])
            result = self.run_case(name+'changed', f"restart_file='{checkpoint}',", fail=True,
                                   common=settings.replace('rates=0,1e-4', 'rates=0,1.1e-4'))
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('profile_bad_times', common=base.replace('0,2e-7,4e-7', '0,0,4e-7'), fail=True)
        self.run_case('profile_bad_rates', common=base.replace('rates=0,1e-4', 'rates=0,-1e-4'), fail=True)
        self.run_case('profile_conflict', common=base+'mass_flow_rate=1e-4,', fail=True)

    def test_nozzle_profiles_and_swirl(self):
        base = ("parcel_count=8,nozzle_radius=3e-4,nozzle_profile='hollow',swirl_angle=.3,hollow_spread=.1,"
                "injection_speed=10,cone_angle=.2,injection_center=.002,.002,.002,injection_axis=0,0,1,"
                "drag_enabled=.false.,heat_enabled=.false.,evaporation_enabled=.false.,")
        ports = ("injection_port_count=2,injection_port_centers=.001,.002,.002,.003,.002,.002,"
                 "injection_port_axes=0,0,1,0,0,1,injection_port_fractions=.5,.5,"
                 "port_nozzle_profiles='solid','hollow',port_swirl_angles=-.4,.4,port_hollow_spreads=0,.1,")
        independent = (ports+"initial_liquid_mass=0,independent_port_schedules=.true.,port_mass_flow_rates=1e-4,2e-4,"
                       "port_injection_starts=1e-7,1e-7,port_injection_ends=5e-7,6e-7,"
                       "port_injection_intervals=2e-7,3e-7,")
        for index, schedule in enumerate(['', ports, independent]):
            name = f'profile{index}'
            settings = base+schedule
            checkpoint = self.root/f'{name}.chk'
            self.run_case(name, common=settings)
            self.check_balance(name)
            with (self.root/f'{name}.parcels.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            for start in range(0, len(rows), 8):
                batch = rows[start:start+8]
                solid = index > 0 and float(batch[0]['x']) < .002
                center_x = .002 if index == 0 else (.001 if solid else .003)
                swirl = .3 if index == 0 else (-.4 if solid else .4)
                for q, row in enumerate(batch):
                    u,v,w = [float(row[k]) for k in ['u','v','w']]
                    age = (float(row['z'])-.002)/w
                    dx,dy = float(row['x'])-age*u-center_x, float(row['y'])-age*v-.002
                    radius = math.hypot(dx,dy)
                    theta = .2*(q+.5)/8 if solid else .2+.1*((q+.5)/8-.5)
                    self.assertAlmostEqual(radius/3e-4, (q+.5)/8 if solid else 1, delta=5e-14)
                    self.assertAlmostEqual(w/10, math.cos(theta), delta=2e-14)
                    self.assertAlmostEqual(math.sqrt(u*u+v*v+w*w), 10, delta=2e-13)
                    self.assertAlmostEqual((dx*u+dy*v)/radius, 10*math.sin(theta)*math.cos(swirl), delta=2e-12)
                    self.assertAlmostEqual((dx*v-dy*u)/radius, 10*math.sin(theta)*math.sin(swirl), delta=2e-12)
            self.run_case(name+'stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
            self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=settings)
            for suffix in ['gas.csv', 'parcels.csv']:
                self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                 (self.root/f'{name}resumed.{suffix}').read_bytes())
            result = self.run_case(name+'changed', f"restart_file='{checkpoint}',", fail=True,
                                   common=settings.replace('swirl_angle=.3', 'swirl_angle=.35'))
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('profile_zero_radius', common=base+'nozzle_radius=0,', fail=True)
        self.run_case('profile_bad_angle', common=base+'swirl_angle=4,', fail=True)
        self.run_case('profile_unused', common=base+"port_nozzle_profiles(1)='solid',", fail=True)

    def test_finite_nozzle_restart(self):
        base = ("parcel_count=8,nozzle_radius=3e-4,injection_speed=0,drag_enabled=.false.,"
                "heat_enabled=.false.,evaporation_enabled=.false.,")
        ports = ("injection_port_count=2,injection_port_centers=.001,.002,.002,.003,.002,.002,"
                 "injection_port_axes=0,0,1,1,0,0,injection_port_fractions=.5,.5,port_nozzle_radii=2e-4,4e-4,")
        independent = (ports+"initial_liquid_mass=0,independent_port_schedules=.true.,port_mass_flow_rates=1e-4,2e-4,"
                       "port_injection_starts=1e-7,1e-7,port_injection_ends=5e-7,6e-7,"
                       "port_injection_intervals=2e-7,3e-7,")
        for index, schedule in enumerate(["injection_center=.002,.002,.002,injection_axis=1,2,3,", ports, independent]):
            name = f'nozzle{index}'
            settings = base+schedule
            checkpoint = self.root/f'{name}.chk'
            self.run_case(name, common=settings)
            self.check_balance(name)
            with (self.root/f'{name}.parcels.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            for start in range(0, len(rows), 8):
                batch = rows[start:start+8]
                if index == 0:
                    center, axis, radius = [.002]*3, [1,2,3], 3e-4
                elif float(batch[0]['x']) < .002:
                    center, axis, radius = [.001,.002,.002], [0,0,1], 2e-4
                else:
                    center, axis, radius = [.003,.002,.002], [1,0,0], 4e-4
                for q, row in enumerate(batch):
                    delta = [float(row[k])-center[j] for j, k in enumerate(['x','y','z'])]
                    self.assertLess(abs(sum(d*a for d,a in zip(delta,axis))), 1e-17)
                    self.assertAlmostEqual(math.sqrt(sum(d*d for d in delta))/radius, (q+.5)/8, delta=2e-14)
            self.run_case(name+'stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
            self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=settings)
            for suffix in ['gas.csv', 'parcels.csv']:
                self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                 (self.root/f'{name}resumed.{suffix}').read_bytes())
            self.assertEqual(self.history(name)[-1], self.history(name+'resumed')[-1])
            changed = settings.replace('radius=3e-4', 'radius=3.1e-4') if index == 0 else settings.replace('radii=2e-4', 'radii=2.1e-4')
            result = self.run_case(name+'changed', f"restart_file='{checkpoint}',", fail=True, common=changed)
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('nozzle_negative', common='nozzle_radius=-1,', fail=True)
        self.run_case('nozzle_outside', common='nozzle_radius=.003,', fail=True)
        self.run_case('nozzle_unused', common='port_nozzle_radii(1)=1e-4,', fail=True)

    def test_port_size_distributions(self):
        base = ("parcel_count=8,injection_port_count=2,injection_port_centers=.001,.001,.001,.003,.001,.001,"
                "injection_port_axes=1,0,0,-1,0,0,injection_port_fractions=.4,.6,injection_speed=0,"
                "injection_port_diameters=2e-5,3e-5,port_diameter_distributions='normal','',"
                "port_diameter_weightings='number','',port_diameter_stds=3e-6,-1,port_diameter_shapes=0,4,"
                "port_diameter_minima=1e-5,2e-5,port_diameter_maxima=3e-5,6e-5,"
                "drag_enabled=.false.,heat_enabled=.false.,evaporation_enabled=.false.,")
        schedules = ["mass_flow_rate=1e-4,injection_start=1e-7,injection_end=4e-7,injection_interval=2e-7,",
                     "initial_liquid_mass=0,independent_port_schedules=.true.,port_mass_flow_rates=1e-4,2e-4,"
                     "port_injection_starts=1e-7,1e-7,port_injection_ends=5e-7,6e-7,"
                     "port_injection_intervals=2e-7,3e-7,"]
        for index, schedule in enumerate(schedules):
            name = f'port_sizes{index}'
            settings = base+schedule
            checkpoint = self.root/f'{name}.chk'
            self.run_case(name, common=settings)
            self.check_balance(name)
            with (self.root/f'{name}.parcels.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            for start in range(0, len(rows), 8):
                batch = rows[start:start+8]
                normal = float(batch[0]['x']) < .002
                if normal:
                    self.assertEqual(len({r['multiplicity'] for r in batch}), 1)
                else:
                    masses = [float(r['mass'])*float(r['multiplicity']) for r in batch]
                    self.assertLess((max(masses)-min(masses))/max(masses), 2e-14)
                for q, row in enumerate(batch):
                    quantile = (q+.5)/8
                    if normal:
                        oracle = NormalDist(2e-5, 3e-6)
                        expected = oracle.inv_cdf(oracle.cdf(1e-5)+quantile*(oracle.cdf(3e-5)-oracle.cdf(1e-5)))
                    else:
                        low, high = math.exp(-(2/3)**4), math.exp(-2**4)
                        expected = 3e-5*(-math.log(low-quantile*(low-high)))**.25
                    self.assertAlmostEqual(float(row['diameter'])/expected, 1, delta=2e-13)
            self.run_case(name+'stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
            self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=settings)
            for suffix in ['gas.csv', 'parcels.csv']:
                self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                 (self.root/f'{name}resumed.{suffix}').read_bytes())
            self.assertEqual(self.history(name)[-1], self.history(name+'resumed')[-1])
            result = self.run_case(name+'changed', f"restart_file='{checkpoint}',", fail=True,
                                   common=settings.replace('stds=3e-6', 'stds=4e-6'))
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('port_sizes_unused', common=base+"port_diameter_distributions(3)='normal',", fail=True)
        self.run_case('port_sizes_invalid', common=base.replace('stds=3e-6', 'stds=-2'), fail=True)

    def test_chi_squared_sizes(self):
        settings = ("parcel_count=8,diameter_distribution='chi-squared',diameter_weighting='number',"
                    "injection_port_count=2,injection_port_centers=.001,.001,.001,.002,.001,.001,"
                    "injection_port_axes=1,0,0,-1,0,0,injection_port_fractions=.5,.5,"
                    "injection_port_diameters=2e-5,4e-5,drag_enabled=.false.,heat_enabled=.false.,"
                    "evaporation_enabled=.false.,")
        checkpoint = self.root/'chi.chk'
        self.run_case('chi', common=settings)
        self.check_balance('chi')
        with (self.root/'chi.parcels.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        # Pinned 100-bin table, conditioned on its positive-diameter entries.
        bins = [12, 18, 23, 28, 33, 38, 46, 61]
        for port, d32 in enumerate([2e-5, 4e-5]):
            batch = rows[8*port:8*(port+1)]
            self.assertEqual(len({r['multiplicity'] for r in batch}), 1)
            for row, bin_index in zip(batch, bins):
                self.assertAlmostEqual(float(row['diameter'])/(bin_index*.12*d32/3), 1, delta=2e-14)
        self.run_case('chi_stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
        self.run_case('chi_resumed', f"restart_file='{checkpoint}',", common=settings)
        for suffix in ['gas.csv', 'parcels.csv']:
            self.assertEqual((self.root/f'chi.{suffix}').read_bytes(),
                             (self.root/f'chi_resumed.{suffix}').read_bytes())
        self.assertEqual(self.history('chi')[-1], self.history('chi_resumed')[-1])
        result = self.run_case('chi_changed', f"restart_file='{checkpoint}',", fail=True,
                               common=settings.replace("distribution='chi-squared'", "distribution='weibull'"))
        self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('chi_truncated', common=settings+'diameter_min=1e-5,diameter_max=2e-5,')
        self.check_balance('chi_truncated')
        with (self.root/'chi_truncated.parcels.csv').open() as stream:
            for row in csv.DictReader(stream):
                # CSV reconstructs d from mass via a cube root; allow a few ulps.
                self.assertTrue(1e-5*(1-4e-15) <= float(row['diameter']) <= 2e-5*(1+4e-15))
        self.run_case('chi_empty', common=settings+'diameter_min=1e-9,diameter_max=2e-9,', fail=True)
        self.run_case('chi_std', common=settings+'diameter_std=1e-6,', fail=True)

    def test_normal_lognormal_sizes(self):
        base = ("parcel_count=8,drop_diameter=2e-5,diameter_std=7e-6,diameter_min=1e-5,diameter_max=4e-5,"
                "diameter_weighting='number',drag_enabled=.false.,heat_enabled=.false.,evaporation_enabled=.false.,")
        ports = ("injection_port_count=2,injection_port_centers=.001,.001,.001,.002,.001,.001,"
                 "injection_port_axes=1,0,0,-1,0,0,injection_port_fractions=.5,.5,")
        independent = (ports+"initial_liquid_mass=0,independent_port_schedules=.true.,port_mass_flow_rates=1e-4,2e-4,"
                       "port_injection_starts=1e-7,1e-7,port_injection_ends=5e-7,6e-7,"
                       "port_injection_intervals=2e-7,3e-7,")
        for model in ['normal', 'lognormal']:
            width = math.sqrt(math.log1p((7e-6/2e-5)**2))
            oracle = NormalDist(2e-5, 7e-6) if model == 'normal' else NormalDist(math.log(2e-5)-width**2/2, width)
            transform = (lambda d: d) if model == 'normal' else math.log
            low, high = oracle.cdf(transform(1e-5)), oracle.cdf(transform(4e-5))
            for index, schedule in enumerate(['', ports, independent]):
                name = f'{model}{index}'
                settings = base+f"diameter_distribution='{model}',"+schedule
                checkpoint = self.root/f'{name}.chk'
                self.run_case(name, common=settings)
                self.check_balance(name)
                with (self.root/f'{name}.parcels.csv').open() as stream:
                    rows = list(csv.DictReader(stream))
                for start in range(0, len(rows), 8):
                    batch = rows[start:start+8]
                    self.assertEqual(len({r['multiplicity'] for r in batch}), 1)
                    for q, row in enumerate(batch):
                        expected = oracle.inv_cdf(low+(q+.5)/8*(high-low))
                        if model == 'lognormal':
                            expected = math.exp(expected)
                        self.assertAlmostEqual(float(row['diameter'])/expected, 1, delta=2e-13)
                self.run_case(name+'stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
                self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=settings)
                for suffix in ['gas.csv', 'parcels.csv']:
                    self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                     (self.root/f'{name}resumed.{suffix}').read_bytes())
                result = self.run_case(name+'changed', f"restart_file='{checkpoint}',", fail=True,
                                       common=settings.replace('diameter_std=7e-6', 'diameter_std=8e-6'))
                self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        for extra in ["diameter_distribution='bad',", "diameter_distribution='normal',diameter_std=-1,",
                      "diameter_distribution='normal',diameter_shape=3,", "diameter_std=1e-6,"]:
            self.run_case('bad_distribution'+str(len(list(self.root.glob('bad_distribution*.nml')))), common=extra, fail=True)

    def test_number_weighted_sizes_restart(self):
        base = ("parcel_count=8,drop_diameter=2e-5,diameter_shape=3,diameter_min=1e-5,diameter_max=4e-5,"
                "diameter_weighting='number',drag_enabled=.false.,heat_enabled=.false.,evaporation_enabled=.false.,")
        ports = ("injection_port_count=2,injection_port_centers=.001,.001,.001,.002,.001,.001,"
                 "injection_port_axes=1,0,0,-1,0,0,injection_port_fractions=.5,.5,")
        schedules = ['', ports, ports+"initial_liquid_mass=0,independent_port_schedules=.true.,"
                     "port_mass_flow_rates=1e-4,2e-4,port_injection_starts=1e-7,1e-7,"
                     "port_injection_ends=5e-7,6e-7,port_injection_intervals=2e-7,3e-7,"]
        for index, schedule in enumerate(schedules):
            name = f'number_sizes{index}'
            settings = base+schedule
            checkpoint = self.root/f'{name}.chk'
            self.run_case(name, common=settings)
            self.check_balance(name)
            self.run_case(name+'stop', f"stop_after_steps=3,checkpoint_file='{checkpoint}',", common=settings)
            self.run_case(name+'resumed', f"restart_file='{checkpoint}',", common=settings)
            for suffix in ['gas.csv', 'parcels.csv']:
                self.assertEqual((self.root/f'{name}.{suffix}').read_bytes(),
                                 (self.root/f'{name}resumed.{suffix}').read_bytes())
            self.assertEqual(self.history(name)[-1], self.history(name+'resumed')[-1])
            with (self.root/f'{name}.parcels.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            for start in range(0, len(rows), 8):
                batch = rows[start:start+8]
                self.assertEqual(len({r['multiplicity'] for r in batch}), 1)
                self.assertEqual(len({r['mass'] for r in batch}), 8)
            result = self.run_case(name+'changed', f"restart_file='{checkpoint}',", fail=True,
                                   common=settings.replace("weighting='number'", "weighting='mass'"))
            self.assertIn('incompatible spray checkpoint', result.stdout.lower())
        self.run_case('invalid_weighting', common="diameter_weighting='invalid',", fail=True)

    def test_les_switch_is_active_and_conservative(self):
        common = 'cells=4,4,4,hydro_reconstruction="characteristic_plm",shear_velocity=50,initial_liquid_mass=0,'
        self.run_case('sgs_off', common=common)
        self.run_case('sgs_on', 'smagorinsky_constant=.16,', common=common)
        self.check_balance('sgs_on')
        self.assertNotEqual((self.root/'sgs_off.gas.csv').read_bytes(), (self.root/'sgs_on.gas.csv').read_bytes())

    def test_existing_outputs_and_checkpoint_are_not_overwritten(self):
        self.run_case('original')
        old = (self.root/'original.gas.csv').read_bytes()
        self.run_case('original', fail=True)
        self.assertEqual((self.root/'original.gas.csv').read_bytes(), old)
        checkpoint = self.root/'existing.chk'; checkpoint.write_bytes(b'owner data\n')
        self.run_case('another', f"checkpoint_file='{checkpoint}',", fail=True)
        self.assertEqual(checkpoint.read_bytes(), b'owner data\n')
        self.assertFalse((self.root/'another.history.csv').exists())

    def test_invalid_checkpoint_records_rejected_before_output(self):
        checkpoint = self.root/'good.chk'
        self.run_case('first', f"stop_after_steps=2,checkpoint_file='{checkpoint}',")
        original = checkpoint.read_text().splitlines()
        parcel_start = 6+8
        mutations = {
            'magic': (0, 'wrong magic'), 'context': (1, 'wrong context'),
            'shape': (2, '15 3 2 2 2'), 'negative_count': (2, '15 2 2 2 -1'),
            'nan': (6, 'NaN '+' '.join(original[6].split()[1:])),
            'null': (6, '/'), 'repeat': (6, '16*1'), 'comma': (6, '1,,2'),
            'extra': (6, original[6]+' 1'),
            'trailing': (len(original), 'not an allowed trailing record'),
            'negative_temperature': (6, '-100 '+' '.join(original[6].split()[1:])),
            'event': (3, original[3].rsplit(' ',1)[0]+' 100'),
            'ledger': (4, '0 '+' '.join(original[4].split()[1:])),
            'duplicate_id': (parcel_start+1, original[parcel_start]),
        }
        for name,(line,value) in mutations.items():
            with self.subTest(name=name):
                changed = original.copy()
                if line==len(changed): changed.append(value)
                else: changed[line] = value
                bad = self.root/(name+'.chk'); bad.write_text('\n'.join(changed)+'\n')
                self.run_case(name, f"restart_file='{bad}',", fail=True)
                self.assertFalse((self.root/(name+'.history.csv')).exists())
        truncated = self.root/'truncated.chk'; truncated.write_text('\n'.join(original[:10])+'\n')
        self.run_case('truncated', f"restart_file='{truncated}',", fail=True)
        self.assertFalse((self.root/'truncated.history.csv').exists())

    def test_changed_physics_rejected_on_restart(self):
        checkpoint=self.root/'good.chk'
        self.run_case('first', f"stop_after_steps=2,checkpoint_file='{checkpoint}',")
        for i,change in enumerate(['liquid_cp=4000,','smagorinsky_constant=.1,','maximum_dt=5e-8,',
                                    'evaporation_enabled=.false.,','mass_flow_rate=1e-6,','diameter_shape=3,']):
            self.run_case(f'changed{i}', change+f"restart_file='{checkpoint}',", fail=True)
            self.assertFalse((self.root/f'changed{i}.history.csv').exists())

    def test_invalid_inputs_rejected(self):
        for i,change in enumerate(['mole_names="UNKNOWN",','cells=1,2,2,','liquid_cp=-1,',
                                   'injection_axis=0,0,0,','maximum_dt=0,','mole_amounts=-1,',
                                   'smagorinsky_constant=2,','injection_center=-1,0,0,',
                                   'diameter_shape=-1,','diameter_min=3e-5,diameter_max=1e-5,']):
            self.run_case(f'invalid{i}', change, fail=True)
            self.assertFalse((self.root/f'invalid{i}.history.csv').exists())


if __name__=='__main__':
    unittest.main(argv=[__file__]+REST)
