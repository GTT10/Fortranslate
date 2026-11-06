from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from check_eb_thermal_mixing import summarize


class MixingTest(unittest.TestCase):
    def rows(self):
        return [dict(variant=v, i=6, j=2, k=3, fraction=.01/2**l,
                     p_left_Pa=101325, p_right_Pa=101325, p_mixed_Pa=101325,
                     linear_pressure_Pa=101325, temperature_left_K=900,
                     temperature_right_K=900, temperature_mixed_K=900,
                     pressure_preserving_energy_gap_J_m3=0)
                for v in range(3) for l in range(3)]

    def test_equilibrium(self):
        self.assertEqual(summarize(self.rows())['pairs_per_variant'], 1)

    def test_missing_variant(self):
        with self.assertRaisesRegex(ValueError, 'coverage'):
            summarize(self.rows()[:6])

    def test_bad_fraction(self):
        r = self.rows();r[1]['fraction'] = .009
        with self.assertRaisesRegex(ValueError, 'sequence'):
            summarize(r)

    def test_isothermal_failure(self):
        r = self.rows();r[-1]['p_mixed_Pa'] += .01
        with self.assertRaisesRegex(AssertionError, 'Isothermal'):
            summarize(r)

    def test_donor_failure(self):
        r = self.rows();r[0]['p_left_Pa'] += 1
        with self.assertRaisesRegex(AssertionError, 'Donor'):
            summarize(r)


if __name__ == '__main__':
    unittest.main()
