from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from localize_eb_pressure_evolution import reduce_cells


class LocationTest(unittest.TestCase):
    def run_case(self, rows, history):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'cells.csv'
            p.write_text('step,node,i,j,k,vfrac,volume_m3,pressure_offset_Pa,E_minus_H_minus_KE_J_m3\n'+rows)
            return reduce_cells(p, history)

    def test_volume_weight_and_signed_peak(self):
        result = self.run_case('0,1,1,1,1,1,1,-4,0\n0,1,2,1,1,1,3,2,0\n',
                               [{'max_pressure_offset_Pa': 4}])[0]
        self.assertEqual(result['peak']['i'], '1')
        self.assertEqual(result['groups']['level0_full']['volume_mean_abs_pressure_offset_Pa'], 2.5)

    def test_duplicate_refused(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            self.run_case('0,1,1,1,1,1,1,0,0\n'*2, [{'max_pressure_offset_Pa': 0}])

    def test_missing_leaf_refused(self):
        with self.assertRaisesRegex(ValueError, 'coverage'):
            self.run_case('0,1,1,1,1,1,1,0,0\n0,2,1,1,1,1,1,0,0\n1,1,1,1,1,1,1,0,0\n',
                          [{'max_pressure_offset_Pa': 0}]*2)

    def test_history_mismatch_refused(self):
        with self.assertRaisesRegex(AssertionError, 'mismatch'):
            self.run_case('0,1,1,1,1,1,1,2,0\n', [{'max_pressure_offset_Pa': 0}])


if __name__ == '__main__':
    unittest.main()
