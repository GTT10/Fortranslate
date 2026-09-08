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
