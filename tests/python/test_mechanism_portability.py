"""Regression coverage for real chemical labels and supported NASA7/falloff forms."""
from __future__ import annotations

import copy
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
import generate_elementary_mechanism as generator
import ingest_cantera_mechanism as ingest


class ChemicalLabels(unittest.TestCase):
    def bundle(self):
        return {'module_name': 'portable_mod', 'loader_name': 'load_portable',
                'kernel_name': 'portable_rates', 'symbol_prefix': 'portable',
                'species': ['CH2(S)', 'CH2', 'CH2-S'],
                'reactions': [{'equation': 'CH2(S) <=> CH2',
                               'reactants': {'CH2(S)': 1}, 'products': {'CH2': 1},
                               'reversible': True, 'arrhenius': {'A': 1, 'b': 0, 'Ea': 0}}]}

    def test_real_labels_have_distinct_deterministic_symbols(self):
        data = self.bundle()
        ingest._validate_species_names(data['species'])
        generator.validate(data)
        text = generator.generate(data)
        for name in data['species']:
            self.assertIn('portable_' + generator.species_symbol_stem(name), text)
        self.assertEqual(len({generator.species_symbol_stem(n) for n in data['species']}), 3)
        self.assertEqual(text, generator.generate(copy.deepcopy(data)))
        self.assertIn('CH2(S) <=> CH2', text)

    def test_legacy_symbols_unchanged(self):
        for name in ['H2O2', 'N2', 'h2o', 'A_1', 'NC16H34']:
            self.assertEqual(generator.species_symbol_stem(name), name.lower())

    def test_invalid_or_injectable_labels_rejected(self):
        for name in ['bad name', 'A"', "A'", 'A/B', 'A\\B', 'A\nB', '', 'A'*25]:
            with self.subTest(name=name):
                data = self.bundle(); data['species'][2] = name
                with self.assertRaises(ValueError):
                    generator.validate(data)
                with self.assertRaises(ValueError):
                    ingest._validate_species_names(data['species'])

    def test_case_and_generated_symbol_collisions_rejected(self):
        for label in ['ch2(s)', generator.species_symbol_stem('CH2(S)')]:
            data = self.bundle(); data['species'][2] = label
            with self.assertRaises(ValueError):
                generator.validate(data)


class SingleIntervalNASA7(unittest.TestCase):
    coefficients = [2.5, 0., 0., 0., 0., 25473.66, -0.44668285]

    def species(self, ranges, rows):
        return SimpleNamespace(name='H', molecular_weight=1.008,
                               thermo=SimpleNamespace(reference_pressure=101325),
                               input_data={'composition': {'H': 1},
                                           'thermo': {'model': 'NASA7',
                                                      'temperature-ranges': ranges, 'data': rows}})

    def test_one_interval_and_cantera_export_are_same_function(self):
        c = self.coefficients
        for ranges, rows in [([200, 6000], [c]), ([200, 6000, 6000], [c, c]),
                             ([200, 200, 6000], [c, c])]:
            record = ingest._thermo_record(self.species(ranges, rows), 0)
            self.assertEqual(record['low_coefficients'], c)
            self.assertEqual(record['high_coefficients'], c)
            self.assertEqual([record['temperature_min'], record['temperature_mid'], record['temperature_max']],
                             [200, 3100, 6000])

    def test_degenerate_unequal_or_invalid_ranges_rejected(self):
        c = self.coefficients; other = c.copy(); other[0] = 3.
        for ranges, rows in [([200, 6000, 6000], [c, other]), ([6000, 200], [c]),
                             ([200, 200], [c]), ([200, 6000], [c, c]),
                             ([200, 6000, 5000], [c, c])]:
            with self.subTest(ranges=ranges, rows=rows):
                with self.assertRaises(ValueError):
                    ingest._thermo_record(self.species(ranges, rows), 0)


class Lindemann(unittest.TestCase):
    def reaction(self):
        return SimpleNamespace(equation='A (+M) <=> B (+M)', reactants={'A': 1}, products={'B': 1},
                               reaction_type='falloff-Lindemann', rate=SimpleNamespace(sub_type='Lindemann'),
                               reversible=True, duplicate=False,
                               third_body=SimpleNamespace(default_efficiency=1., efficiencies={'A': 2.}),
                               input_data={'type': 'falloff',
                                           'low-P-rate-constant': {'A': 2., 'b': 0., 'Ea': 0.},
                                           'high-P-rate-constant': {'A': 3., 'b': 0., 'Ea': 0.}})

    def test_lindemann_uses_existing_no_troe_runtime_path(self):
        record = ingest._reaction_record(self.reaction(), ['A', 'B'], 0)
        self.assertEqual(record['type'], 'falloff')
        self.assertNotIn('troe', record)
        self.assertEqual(record['efficiencies'], {'A': 2.})
        data = {'module_name': 'lin_mod', 'loader_name': 'load_lin', 'kernel_name': 'lin_rates',
                'symbol_prefix': 'lin', 'species': ['A', 'B'], 'reactions': [record]}
        generator.validate(data)
        text = generator.generate(data)
        self.assertIn('reaction_kind_falloff', text)
        self.assertNotIn('%troe%enabled = .true.', text)

    def test_unknown_or_malformed_falloff_not_silently_lindemann(self):
        for subtype, kind in [('Troe', 'falloff-Troe'), ('', 'falloff-Lindemann'),
                              ('Lindemann', 'falloff-unknown'), ('SRI', 'falloff-SRI'),
                              ('Tsang', 'falloff-Tsang')]:
            r = self.reaction(); r.rate.sub_type = subtype; r.reaction_type = kind
            with self.subTest(subtype=subtype, kind=kind):
                with self.assertRaises(ValueError):
                    ingest._reaction_record(r, ['A', 'B'], 0)


if __name__ == '__main__':
    unittest.main()
