"""End-to-end spray CLI, conservation, injection and strict restart tests."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
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
