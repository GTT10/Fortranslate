"""Dependency-free regressions for the pre-configure provenance check."""
import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
import check_mechanism_provenance as checker


class ProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "gas.yaml"
        self.source.write_bytes(b"description: test\n")
        self.path = self.root / "bundle.json"
        self.bundle = {"source": {"file": "gas.yaml", "phase": "gas",
            "sha256": hashlib.sha256(self.source.read_bytes()).hexdigest()},
            "reactions": [{"A": 1.0}]}
        self.save()
        self.out = io.StringIO()
        self.enterContext(contextlib.redirect_stdout(self.out))
        self.enterContext(contextlib.redirect_stderr(self.out))

    def save(self):
        self.path.write_text(json.dumps(self.bundle), encoding="utf-8")

    def test_matching_source(self):
        self.assertTrue(checker.check(self.path))

    def test_changed_source_rejected_without_writing_inputs(self):
        self.source.write_bytes(b"description: changed\n")
        before = (self.path.read_bytes(), self.source.read_bytes())
        self.assertFalse(checker.check(self.path))
        self.assertEqual(before, (self.path.read_bytes(), self.source.read_bytes()))

    def test_line_endings_are_not_normalized(self):
        self.source.write_bytes(b"description: test\r\n")
        self.assertFalse(checker.check(self.path))

    def test_missing_source_is_error(self):
        self.source.unlink()
        self.assertEqual(checker.main(["--bundle", str(self.path)]), 2)

    def test_malformed_hash_rejected(self):
        for value in (None, "", "a" * 63, "G" * 64):
            with self.subTest(value=value):
                self.bundle["source"]["sha256"] = value
                self.save()
                self.assertEqual(checker.main(["--bundle", str(self.path)]), 2)

    def test_duplicate_keys_rejected(self):
        self.path.write_text('{"source": {}, "source": {}}', encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            checker.read_bundle(self.path)

    def test_nonfinite_json_rejected(self):
        self.bundle["reactions"][0]["A"] = float("nan")
        self.save()
        with self.assertRaisesRegex(ValueError, "non-finite"):
            checker.read_bundle(self.path)

    def test_implicit_path_traversal_rejected(self):
        for name in ("../gas.yaml", "/gas.yaml", "..\\gas.yaml"):
            with self.subTest(name=name):
                self.bundle["source"]["file"] = name
                self.save()
                self.assertEqual(checker.main(["--bundle", str(self.path)]), 2)

    def test_explicit_source_allowed(self):
        self.bundle["source"]["file"] = "external/gas.yaml"
        self.save()
        self.assertTrue(checker.check(self.path, self.source))

    def test_output_cannot_overwrite_inputs(self):
        for output in (self.path, self.source):
            with self.assertRaisesRegex(ValueError, "overwrite"):
                checker.check(self.path, regenerate_output=True, output=output)

    def test_output_requires_regeneration(self):
        with self.assertRaisesRegex(ValueError, "requires"):
            checker.check(self.path, output=self.root / "candidate.json")

    def test_matching_regeneration(self):
        with patch.object(checker, "regenerate", return_value=self.bundle):
            self.assertTrue(checker.check(self.path, regenerate_output=True))

    def test_regeneration_detects_chemistry_changes(self):
        candidate = copy.deepcopy(self.bundle)
        candidate["reactions"][0]["A"] = 2.0
        before = self.path.read_bytes()
        output = self.root / "evidence" / "candidate.json"
        with patch.object(checker, "regenerate", return_value=candidate):
            self.assertFalse(checker.check(self.path, regenerate_output=True,
                                           output=output))
        self.assertEqual(json.loads(output.read_text()), candidate)
        self.assertEqual(self.path.read_bytes(), before)
        self.assertIn('"A": 2.0', self.out.getvalue())

    def test_runtime_failure_is_error(self):
        with patch.object(checker, "regenerate", side_effect=RuntimeError("pin")):
            self.assertEqual(checker.main(["--bundle", str(self.path),
                                           "--regenerate"]), 2)


if __name__ == "__main__":
    unittest.main()
