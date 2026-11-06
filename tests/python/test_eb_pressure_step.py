from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from check_eb_pressure_step import inventory_check


class InventoryTest(unittest.TestCase):
    def rows(self):
        return [{'component':c,'initial':v,'final':v} for c,v in
                [('1','1'),('2','0'),('3','0'),('4','0'),('5','100'),('6','0.25'),('7','0.75'),('held','99')]]

    def test_closed(self):
        self.assertTrue(all(v==0 for v in inventory_check(self.rows()).values()))

    def test_energy_leak(self):
        rows=self.rows();rows[4]['final']='100.001'
        with self.assertRaisesRegex(AssertionError,'inventory'):inventory_check(rows)

    def test_mass_species_inconsistent(self):
        rows=self.rows();rows[5].update(initial='0.2',final='0.2')
        with self.assertRaisesRegex(AssertionError,'Species sum'):inventory_check(rows)

    def test_nonfinite(self):
        rows=self.rows();rows[7]['final']='nan'
        with self.assertRaisesRegex(ValueError,'Nonfinite'):inventory_check(rows)


if __name__=='__main__':unittest.main()
