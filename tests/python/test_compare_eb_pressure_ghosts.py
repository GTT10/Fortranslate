import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from compare_eb_pressure_ghosts import compare_record


class PressureGhostEvidenceTest(unittest.TestCase):
    def test_partial_run_is_not_a_success(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'run.json';path.write_text(json.dumps({'status':'running'}))
            with self.assertRaisesRegex(ValueError,'incomplete'):
                compare_record(path)

    def test_changed_binary_refused_before_reading_fields(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);binary=root/'probe';binary.write_text('changed')
            record={'status':'completed','binary_sha256':'old','commands':[
                {'returncode':0,'argv':[str(binary),str(root/'input'),str(root/'fields')],'input_sha256':'old'}]}
            path=root/'run.json';path.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError,'identity'):
                compare_record(path)


if __name__=='__main__':
    unittest.main()
