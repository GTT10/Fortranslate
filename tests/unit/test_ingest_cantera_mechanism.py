#!/usr/bin/env python3
"""Unit tests for the deterministic Cantera-to-JSON mechanism ingestion."""

from __future__ import annotations

import hashlib
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ingest_cantera_mechanism as ingest  # noqa: E402


try:
    import cantera as _cantera  # noqa: F401

    HAVE_CANTERA = True
except ImportError:
    HAVE_CANTERA = False


class IngestHelpersTest(unittest.TestCase):
    def test_requires_pinned_cantera_runtime_and_phase_models(self) -> None:
        with self.assertRaisesRegex(ValueError, "expected 3.2.0"):
            ingest._validate_cantera_runtime(SimpleNamespace(
                __version__="3.1.0", __git_commit__="old"
            ))
        with self.assertRaisesRegex(ValueError, "git commit"):
            ingest._validate_cantera_runtime(SimpleNamespace(__version__="3.2.0"))
        with self.assertRaisesRegex(ValueError, "expected 4a8358e"):
            ingest._validate_cantera_runtime(SimpleNamespace(
                __version__="3.2.0", __git_commit__="different"
            ))

        valid = {
            "name": "gas",
            "thermo_model": "ideal-gas",
            "kinetics_model": "bulk",
            "transport_model": "mixture-averaged",
        }
        for attribute, value, message in (
            ("thermo_model", "", "thermo model"),
            ("kinetics_model", "surface", "kinetics model"),
            ("transport_model", "multicomponent", "transport model"),
        ):
            phase = SimpleNamespace(**{**valid, attribute: value})
            with self.subTest(attribute=attribute):
                with self.assertRaisesRegex(ValueError, message):
                    ingest._validate_phase_models(phase)

    def test_reads_inline_and_block_source_units(self) -> None:
        inline = (
            b"cantera-version: 3.2.0\n"
            b"units: {length: cm, quantity: mol}\n"
        )
        block = (
            b"cantera-version: 3.2.0\nunits:\n"
            b"  length: cm\n  quantity: mol # comment\nphases: []\n"
        )
        self.assertEqual(
            ingest._yaml_metadata(inline),
            ("3.2.0", {"length": "cm", "quantity": "mol"}),
        )
        self.assertEqual(
            ingest._yaml_metadata(block),
            ("3.2.0", {"length": "cm", "quantity": "mol"}),
        )

    def test_rejects_unsupported_rate_families(self) -> None:
        base = {
            "equation": "A <=> B",
            "reactants": {"A": 1.0},
            "products": {"B": 1.0},
            "reversible": True,
        }
        unsupported = (
            ({"type": "pressure-dependent-Arrhenius", "rate-constants": []}, "PLOG"),
            ({"type": "Chebyshev", "data": []}, "Chebyshev"),
            ({"type": "falloff", "SRI": {}}, "SRI"),
            ({"type": "falloff", "Tsang": {}}, "Tsang"),
            ({"type": "chemically-activated"}, "chemically-activated"),
            ({"type": "elementary", "orders": {"A": 0.5}}, "custom reaction"),
            ({"type": "elementary", "reverse-rate-constant": {}}, "reverse rates"),
        )
        for raw, label in unsupported:
            reaction = SimpleNamespace(
                input_data={**base, **raw},
                equation=base["equation"],
                rate=SimpleNamespace(),
                reaction_type=str(raw["type"]),
                third_body=None,
            )
            with self.subTest(label=label):
                with self.assertRaisesRegex(ValueError, label):
                    ingest._reaction_kind(reaction, reaction.input_data, base["equation"])

    def test_rejects_lindemann_falloff(self) -> None:
        reaction = SimpleNamespace(
            input_data={
                "type": "falloff",
                "low-P-rate-constant": {"A": 1.0, "b": 0.0, "Ea": 0.0},
                "high-P-rate-constant": {"A": 1.0, "b": 0.0, "Ea": 0.0},
            },
            equation="A + M <=> B + M",
            reactants={"A": 1.0},
            products={"B": 1.0},
            reversible=True,
            duplicate=False,
            rate=SimpleNamespace(),
            reaction_type="falloff-Lindemann",
            third_body=SimpleNamespace(default_efficiency=1.0, efficiencies={}),
        )
        with self.assertRaisesRegex(ValueError, "Lindemann"):
            ingest._reaction_record(reaction, ["A", "B"], 0)

    def test_rejects_api_level_custom_orders_and_chemically_activated_rate(self) -> None:
        raw = {"type": "elementary"}
        ordered = SimpleNamespace(
            rate=SimpleNamespace(sub_type=""),
            reaction_type="elementary-Arrhenius",
            third_body=None,
            orders={"A": 0.5},
            allow_negative_orders=False,
            allow_nonreactant_orders=False,
        )
        with self.assertRaisesRegex(ValueError, "custom reaction orders"):
            ingest._reaction_kind(ordered, raw, "A <=> B")

        activated = SimpleNamespace(
            rate=SimpleNamespace(sub_type="", chemically_activated=True),
            reaction_type="falloff",
            third_body=SimpleNamespace(),
            orders={},
            allow_negative_orders=False,
            allow_nonreactant_orders=False,
        )
        with self.assertRaisesRegex(ValueError, "chemically-activated"):
            ingest._reaction_kind(activated, {"type": "falloff"}, "A <=> B")

    def test_rejects_nasa9_and_missing_transport(self) -> None:
        nasa9 = SimpleNamespace(
            name="A",
            molecular_weight=1.0,
            input_data={
                "thermo": {
                    "model": "NASA9",
                    "temperature-ranges": [200.0, 1000.0, 3000.0],
                    "data": [],
                }
            },
        )
        with self.assertRaisesRegex(ValueError, "NASA9"):
            ingest._thermo_record(nasa9, 0)

        missing_transport = SimpleNamespace(
            name="A",
            input_data={
                "thermo": {
                    "model": "NASA7",
                    "temperature-ranges": [200.0, 1000.0, 3000.0],
                    "data": [[2.5] * 7, [2.5] * 7],
                }
            },
        )
        with self.assertRaisesRegex(ValueError, "missing transport"):
            ingest._transport_record(missing_transport, 0)

    def test_rejects_unbalanced_elementary_reaction(self) -> None:
        reaction = SimpleNamespace(
            input_data={
                "type": "elementary",
                "rate-constant": {"A": 1.0, "b": 0.0, "Ea": 0.0},
            },
            equation="A <=> B",
            reactants={"A": 1.0},
            products={"B": 1.0},
            reversible=True,
            duplicate=False,
            rate=SimpleNamespace(),
            third_body=None,
        )
        with self.assertRaisesRegex(ValueError, "element balance"):
            ingest._reaction_record(
                reaction,
                ["A", "B"],
                0,
                {"A": {"H": 2.0}, "B": {"H": 1.0}},
            )

    def test_atomic_publication_preserves_existing_output_on_replace_error(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "bundle.json"
            output.write_text("sentinel\n", encoding="utf-8")
            with mock.patch.object(ingest.os, "replace", side_effect=OSError("boom")):
                with self.assertRaisesRegex(OSError, "boom"):
                    ingest.write_bundle({"value": 1}, output)
            self.assertEqual(output.read_text(encoding="utf-8"), "sentinel\n")
            self.assertEqual([path.name for path in Path(directory).iterdir()], [
                "bundle.json"
            ])


@unittest.skipUnless(HAVE_CANTERA, "Cantera Python package is not installed")
class PinnedCanteraBundleTest(unittest.TestCase):
    yaml_path = ROOT / "mechanisms" / "h2o2_cantera.yaml"

    def test_bundle_is_complete_ordered_and_generator_inputtable(self) -> None:
        metadata = {
            "phase": "ohmech",
            "module_name": "h2o2_full_mechanism_mod",
            "loader_name": "load_h2o2_full_mechanism",
            "kernel_name": "h2o2_full_production_rates",
            "jacobian_name": "h2o2_full_mass_fraction_jacobian",
            "thermo_loader_name": "load_h2o2_full_thermo",
            "transport_loader_name": "load_h2o2_full_transport",
            "symbol_prefix": "h2o2_full",
        }
        bundle = ingest.build_bundle(self.yaml_path, **metadata)
        self.assertEqual(bundle["species"], [
            "H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"
        ])
        self.assertEqual(len(bundle["thermo"]), 10)
        self.assertEqual(len(bundle["transport"]), 10)
        self.assertEqual(len(bundle["reactions"]), 29)
        self.assertEqual(bundle["thermo"][0]["low_coefficients"][0], 2.34433112)
        self.assertEqual(bundle["thermo"][0]["composition"], {"H": 2.0})
        self.assertEqual(bundle["thermo"][0]["reference_pressure"], 101325.0)
        self.assertEqual(bundle["transport"][5]["geometry"], "nonlinear")

        phase = ingest.load_cantera_phase(self.yaml_path, "ohmech")
        expected_duplicate = [bool(phase.reaction(i).duplicate) for i in range(phase.n_reactions)]
        self.assertEqual(
            [item["duplicate"] for item in bundle["reactions"]], expected_duplicate
        )
        self.assertEqual(
            [item["source_index"] for item in bundle["reactions"]],
            list(range(1, 30)),
        )
        self.assertEqual(sum(expected_duplicate), 6)
        self.assertEqual(bundle["reactions"][21]["type"], "falloff")
        self.assertEqual(bundle["reactions"][21]["troe"]["T2"], 5182.0)
        self.assertEqual(bundle["reactions"][0]["default_efficiency"], 1.0)

        digest = hashlib.sha256(self.yaml_path.read_bytes()).hexdigest()
        self.assertEqual(bundle["source"]["sha256"], digest)
        self.assertEqual(bundle["source"]["file"], self.yaml_path.name)
        self.assertEqual(bundle["source"]["cantera_version"], "2.5.0")
        self.assertEqual(bundle["source"]["runtime_cantera_version"], "3.2.0")
        self.assertEqual(bundle["source"]["runtime_cantera_git_commit"], "4a8358e")
        self.assertEqual(bundle["source"]["source_units"]["quantity"], "mol")
        self.assertEqual(bundle["source"]["target_units"]["quantity"], "kmol")
        self.assertEqual(bundle["source"]["phase"], "ohmech")
        self.assertEqual(bundle["schema_version"], 1)
        self.assertEqual(bundle["chemistry_integrator"], "implicit")
        for index, collider in (
            (6, "O2"), (7, "H2O"), (8, "N2"), (9, "AR"),
            (12, "H2"), (13, "H2O"),
        ):
            reaction = bundle["reactions"][index]
            self.assertEqual(reaction["type"], "three-body")
            self.assertEqual(reaction["default_efficiency"], 0.0)
            self.assertEqual(reaction["efficiencies"], {collider: 1.0})

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            first = directory_path / "first.json"
            second = directory_path / "second.json"
            generated = directory_path / "generated.F90"
            ingest.write_bundle(bundle, first)
            ingest.write_bundle(bundle, second)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "tools" / "generate_elementary_mechanism.py"),
                    "--input",
                    str(first),
                    "--output",
                    str(generated),
                ],
                check=True,
            )
            self.assertIn("load_h2o2_full_thermo", generated.read_text(encoding="utf-8"))
            self.assertIn("load_h2o2_full_transport", generated.read_text(encoding="utf-8"))

    def test_cli_failure_does_not_create_partial_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "rejected.json"
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "tools" / "ingest_cantera_mechanism.py"),
                    "--input",
                    str(self.yaml_path),
                    "--phase",
                    "does-not-exist",
                    "--output",
                    str(output),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(output.exists())
            self.assertIn("cannot load Cantera YAML", result.stderr)

    def test_rejects_nonideal_phase(self) -> None:
        with self.assertRaisesRegex(ValueError, "only ideal-gas is supported"):
            ingest.build_bundle(self.yaml_path, phase="ohmech-RK")


if __name__ == "__main__":
    unittest.main()
