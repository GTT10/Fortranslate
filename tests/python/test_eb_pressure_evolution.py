from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from check_eb_pressure_evolution import validate


class EvolutionHistoryTest(unittest.TestCase):
    def fixture(self,prefix,times=('0','1e-8'),leak=False):
        for policy in ('baseline','paired'):
            Path(str(prefix)+f'_{policy}_evolution.csv').write_text(
                'step,time_s,max_pressure_offset_Pa\n'+''.join(f'{i},{t},0\n' for i,t in enumerate(times)))
            rows=['step,component,value']
            for it in range(len(times)):
                for c,v in [('1',1),('5',100),('6',1),('held',99)]:
                    rows.append(f'{it},{c},{v+(0.1 if leak and it==1 and c=="5" else 0)}')
            Path(str(prefix)+f'_{policy}_evolution_inventory.csv').write_text('\n'.join(rows)+'\n')

    def test_closed_history(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'case';self.fixture(p)
            self.assertEqual(validate(p,1,1e-8,True)['status'],'passed')

    def test_missing_endpoint_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'case';self.fixture(p,('0',))
            with self.assertRaisesRegex(ValueError,'Incomplete'):validate(p,1,1e-8)

    def test_bad_clock_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'case';self.fixture(p,('0','2e-8'))
            with self.assertRaisesRegex(ValueError,'clock'):validate(p,1,1e-8)

    def test_explicit_energy_correction_closes_budget(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'case';self.fixture(p,leak=True)
            for policy in ('baseline','paired'):
                path=Path(str(p)+f'_{policy}_evolution.csv')
                lines=path.read_text().splitlines()
                lines[0]+=',numerical_energy_correction_J,absolute_energy_correction_J'
                lines[1]+=',0,0';lines[2]+=',0.1,0.1'
                path.write_text('\n'.join(lines)+'\n')
            self.assertEqual(validate(p,1,1e-8)['status'],'passed')

    def test_impossible_correction_ledger_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'case';self.fixture(p,leak=True)
            path=Path(str(p)+'_baseline_evolution.csv');lines=path.read_text().splitlines()
            lines[0]+=',numerical_energy_correction_J,absolute_energy_correction_J'
            lines[1]+=',0,0';lines[2]+=',0.1,0.01'
            path.write_text('\n'.join(lines)+'\n')
            with self.assertRaisesRegex(AssertionError,'ledger'):validate(p,1,1e-8)

    def test_kinetic_decomposition_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'case';self.fixture(p)
            for policy in ('baseline','paired'):
                path=Path(str(p)+f'_{policy}_evolution.csv')
                lines=path.read_text().splitlines()
                lines[0]+=',kinetic_energy_J,signed_E_H_KE_gap_J'
                lines[1]+=',1,0';lines[2]+=',2,0'
                path.write_text('\n'.join(lines)+'\n')
            with self.assertRaisesRegex(AssertionError,'decomposition'):validate(p,1,1e-8)

    def test_evolving_energy_leak_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'case';self.fixture(p,leak=True)
            with self.assertRaisesRegex(AssertionError,'inventory'):validate(p,1,1e-8)


if __name__=='__main__':unittest.main()
