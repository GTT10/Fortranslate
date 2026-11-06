#!/usr/bin/env python3
"""Reject native activation mappings with the right values on wrong reactions."""
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
import prepare_lidryer_native_activation as prepare


class ActivationMappingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = (ROOT / 'mechanisms/lidryer.yaml').read_text()
        cls.rc = 1.9872155831739962
        reactions = re.findall(r'^- equation:.*?(?=^- equation:|\Z)', cls.source,
                               re.MULTILINE | re.DOTALL)
        lines = [f'ruc = {cls.rc};', 'comp_qfqr(', '// Evaluate the kfs']
        for index, reaction in enumerate(reactions):
            equation = reaction.splitlines()[0].split(':', 1)[1].split('#', 1)[0]
            lines.append(f'// reaction {index}: {equation}')
            rates = dict(re.findall(r'^  ([\w-]*rate-constant):[^\n]*Ea:\s*([-+0-9.eE]+)',
                                    reaction, re.MULTILINE))
            for branch, value in rates.items():
                variable = 'redP' if branch == 'low-P-rate-constant' else 'k_f'
                theta = float(value) / cls.rc
                expression = f'exp(-({theta:.17e}) * invT)' if theta else '1.0'
                lines.append(f'{variable} = {expression};')
        lines.append('productionRate(')
        cls.header = '\n'.join(lines)

    def test_maps_all_reactions_and_pressure_branches(self):
        result = prepare.activation_mapping(self.source, self.header, self.rc)
        self.assertEqual(len(result), 23)
        self.assertEqual({r['source_index'] for r in result}, set(range(1, 22)))
        self.assertEqual(sum(r['activation_temperature_K'] == 0 for r in result), 6)
        self.assertEqual(sum(r['branch'] == 'low-P-rate-constant' for r in result), 2)
        changed, rc, temperatures = prepare.prepare_source(self.source, self.header)
        self.assertEqual(rc, self.rc)
        self.assertEqual(len(temperatures), 23)
        self.assertIn('activation-energy: K', changed)
        self.assertIn('activation-energy: cal/mol', self.source)

    def test_rejects_swapped_values_with_unchanged_multiset(self):
        values = re.findall(r'\(([-+0-9.eE]+)\) \* invT', self.header)
        a, b = values[:2]
        changed = self.header.replace(a, 'SWAP', 1).replace(b, a, 1).replace('SWAP', b, 1)
        with self.assertRaisesRegex(ValueError, 'reaction 1 rate-constant'):
            prepare.prepare_source(self.source, changed)

    def test_rejects_swapped_falloff_branches(self):
        changed = self.header.replace('k_f =', 'TEMP =', 1)  # first reaction must have k_f
        with self.assertRaisesRegex(ValueError, 'Invalid native k_f'):
            prepare.prepare_source(self.source, changed)
        block = self.header.split('// reaction 15:')[1].split('// reaction 16:')[0]
        swapped = block.replace('k_f =', 'TEMP =').replace('redP =', 'k_f =').replace('TEMP =', 'redP =')
        with self.assertRaisesRegex(ValueError, 'reaction 16'):
            prepare.prepare_source(self.source, self.header.replace(block, swapped))

    def test_rejects_nonzero_activation_on_zero_branch(self):
        changed = self.header.replace('k_f = 1.0;', 'k_f = exp(-(1.0) * invT);', 1)
        with self.assertRaisesRegex(ValueError, 'activation mismatch'):
            prepare.prepare_source(self.source, changed)

    def test_rejects_unknown_inverse_temperature_expression(self):
        changed = self.header.replace(' * invT', ' * invT * invT', 1)
        with self.assertRaisesRegex(ValueError, 'Unsupported native activation'):
            prepare.prepare_source(self.source, changed)

    def test_rejects_wrong_equation_and_duplicate_index(self):
        changed = self.header.replace('H + O2 <=> O + OH', 'H + O2 <=> O + H', 1)
        with self.assertRaisesRegex(ValueError, 'equation mismatch'):
            prepare.prepare_source(self.source, changed)
        with self.assertRaisesRegex(ValueError, 'reaction mapping'):
            prepare.prepare_source(self.source, self.header.replace('// reaction 1:', '// reaction 0:', 1))

    def test_stoichiometry_ignores_order_but_not_coefficients(self):
        self.assertEqual(prepare.equation_key('H2 + O <=> OH + H'),
                         prepare.equation_key('O + H2 <=> H + OH'))
        self.assertEqual(prepare.equation_key('H2O2 (+M) <=> OH + OH (+M)'),
                         prepare.equation_key('H2O2 (+M) <=> 2 OH (+M)'))
        self.assertNotEqual(prepare.equation_key('H2 <=> 2 H'),
                            prepare.equation_key('H2 <=> H'))


if __name__ == '__main__':
    unittest.main()
