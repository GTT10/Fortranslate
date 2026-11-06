"""Near-uniform trace species can control well-resolved momentum slopes."""
import copy
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    'conditioning', ROOT/'tools/analyze_eb_limiter_conditioning.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ConditioningTest(unittest.TestCase):
    def fixture(self):
        reference = [[1., 10., 0., 0., 1e-8] for _ in range(7)]
        reference[1][1], reference[2][1] = 8., 12.
        reference[1][4], reference[2][4] = 1e-8-1e-12, 1e-8+9e-12
        candidate = copy.deepcopy(reference)
        candidate[1][4], candidate[2][4] = 1e-8-2e-12, 1e-8+8e-12
        return dict(labels=['rho', 'mx', 'my', 'mz', 'trace'], samples=[dict(
            name='known_trace_control', reference=reference, candidate=candidate,
            offset=[.25, -.25, -.25])])

    def test_trace_controls_velocity_and_substitution(self):
        out = MODULE.analyze(self.fixture())
        self.assertFalse(out['qualification'])
        sample = out['samples'][0]
        self.assertEqual(sample['reference_axes'][0]['controls'], ['trace'])
        self.assertAlmostEqual(sample['reference_axes'][0]['factor'], .4, places=10)
        self.assertAlmostEqual(sample['candidate_axes'][0]['factor'], .8, places=10)
        self.assertAlmostEqual(sample['velocity_difference_m_s'][0], .2, places=10)
        trace = sample['substitutions'][4]
        self.assertEqual(trace['velocity_difference_m_s'], [0., 0., 0.])
        self.assertEqual(sample['substitutions'][1]['velocity_difference_m_s'],
                         sample['velocity_difference_m_s'])
        decomposition = sample['difference_decomposition']
        self.assertAlmostEqual(decomposition['directions'][0]['factor_change'][1], .2, places=10)
        self.assertEqual(decomposition['directions'][0]['slope_change'][1], 0.)
        self.assertLess(max(map(abs, decomposition['velocity_closure_residual_m_s'])), 1e-14)

    def test_density_center_and_slope_close_together(self):
        spec = self.fixture()
        # Changed uniform density, central momentum and its slope simultaneously.
        candidate = spec['samples'][0]['candidate']
        for row in candidate:
            row[0] = 2.
            row[1] += 3.
        candidate[1][1] -= 1.
        candidate[2][1] += 1.
        sample = MODULE.analyze(spec)['samples'][0]
        d = sample['difference_decomposition']
        self.assertEqual(d['center_change'][:2], [1., 3.])
        self.assertAlmostEqual(d['directions'][0]['slope_change'][1], .2, places=10)
        self.assertAlmostEqual(d['velocity_momentum_term_m_s'][0], 1.7, places=10)
        self.assertAlmostEqual(d['velocity_density_term_m_s'][0], -5.1, places=10)
        self.assertAlmostEqual(sample['velocity_difference_m_s'][0], -3.4, places=10)
        self.assertLess(max(map(abs, d['state_closure_residual'])), 1e-14)
        self.assertLess(max(map(abs, d['velocity_closure_residual_m_s'])), 1e-14)

    def test_refusals(self):
        for mode in ['nan', 'density', 'shape', 'offset', 'labels']:
            with self.subTest(mode=mode):
                spec = self.fixture()
                sample = spec['samples'][0]
                if mode == 'nan': sample['candidate'][0][4] = float('nan')
                if mode == 'density': sample['reference'][0][0] = 0.
                if mode == 'shape': sample['candidate'].pop()
                if mode == 'offset': sample['offset'][0] = .5
                if mode == 'labels': spec['labels'][1] = 'rho'
                with self.assertRaises(ValueError): MODULE.analyze(spec)


if __name__ == '__main__':
    unittest.main()
