"""Check physical face pairing and rejection independently of the solver."""
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('audit',Path(__file__).resolve().parents[2]/'tools/audit_eb_initial_interface.py')
audit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class InterfaceAuditTest(unittest.TestCase):
    def fixture(self, prefix, bad=False):
        with Path(str(prefix)+'_initial_pressure.csv').open('w') as f:
            writer=csv.writer(f);writer.writerow(['node','i','j','k','pressure_Pa'])
            for node in (1,2):
                for i in range(3):
                    for j in (1,2):
                        for k in (1,2):
                            writer.writerow([node,i,j,k,'nan' if bad else 101325.])
        with Path(str(prefix)+'_initial_flux.csv').open('w') as f:
            writer=csv.writer(f);writer.writerow(['node','axis','i','j','k','aperture',*audit.FIELDS])
            # Fine lower x face maps to periodic coarse upper x face (i=2).
            writer.writerow([1,1,2,1,1,.5,3,101325,6])
            for j in (1,2):
                for k in (1,2):
                    aperture=1 if j==1 else 0
                    writer.writerow([2,1,0,j,k,aperture,3,101325,6])

    def test_periodic_pair_and_open_area(self):
        with tempfile.TemporaryDirectory() as d:
            prefix=Path(d)/'test';self.fixture(prefix)
            result=audit.summarize(prefix,2)
            self.assertEqual(result['interface_coarse_faces'],1)
            self.assertEqual(result['interface_fine_faces'],2)
            self.assertEqual(result['max_area_fraction_mismatch'],0)
            self.assertTrue(all(v==0 for v in result['max_area_weighted_flux_mismatch'].values()))

    def test_regular_restriction_uses_all_children(self):
        rows=['1 1 1 1 99 999', '1 2 2 2 123 456']
        for i in (1,2):
            for j in (1,2):
                for k in (1,2):
                    value=i+2*j+4*k
                    rows.append(f'2 {i} {j} {k} {value} {value*3}')
        result=audit.restrict_regular_records(rows,2)
        self.assertEqual(list(map(float,result[0].split()[4:])),[10.5,31.5])
        self.assertEqual(result[1:],rows[1:])

    def test_replay_conservation_failure_refused(self):
        with tempfile.TemporaryDirectory() as d:
            prefix=Path(d)/'test'
            Path(str(prefix)+'_initial_stencils.csv').write_text(
                'halo_scaled_residual,child_mean_scaled_residual,pressure_envelope_excess_Pa\n0,0.01,0\n')
            Path(str(prefix)+'_initial_restriction.csv').write_text(
                'gas_scaled_residual,held_scaled_residual\n0,0\n')
            with self.assertRaisesRegex(AssertionError,'child-average'):
                audit.summarize_stencils(prefix)

    def test_nonfinite_refused(self):
        with tempfile.TemporaryDirectory() as d:
            prefix=Path(d)/'test';self.fixture(prefix,bad=True)
            with self.assertRaisesRegex(ValueError,'Nonfinite'):
                audit.summarize(prefix,2)


if __name__=='__main__':
    unittest.main()
