from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from check_eb_pressure_sync import summarize, PRE, POST


class SyncTest(unittest.TestCase):
    def row(self, pre=-2, post=3):
        return dict(stage=1, node=1, i=1, j=1, k=1, **{PRE: pre, POST: post})

    def test_signed_change(self):
        d = summarize([self.row()], {(1,1,1,1)}, 1)[0]['nodes']['1']
        self.assertEqual(d['max_sync_pressure_change_Pa'], 5)
        self.assertEqual(d['max_pre_pressure_offset_Pa'], 2)

    def test_coverage(self):
        with self.assertRaisesRegex(ValueError, 'coverage'):
            summarize([self.row()], {(1,2,1,1)}, 1)

    def test_missing_stage(self):
        with self.assertRaisesRegex(ValueError, 'clocks'):
            summarize([self.row()], {(1,1,1,1)}, 2)

    def test_duplicate(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            summarize([self.row()]*2, {(1,1,1,1)}, 1)


if __name__ == '__main__':
    unittest.main()
