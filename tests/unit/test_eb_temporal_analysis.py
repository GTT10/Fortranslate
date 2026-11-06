"""Known linear-time fields verify units, composite masks, species and orders."""
import csv
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('temporal', ROOT/'tools/analyze_eb_temporal_convergence.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


class TemporalAnalysisTest(unittest.TestCase):
    def test_linear_fields_and_clock_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            runs = []
            for dt in [1., .5, .25]:
                directory = Path(tmp)/str(dt)
                directory.mkdir()
                step = round(1/dt)
                candidate, native = [], []
                for i, slope, volume in [(1, 1., 1.), (2, -3., 3.)]:
                    rho = 10+dt*slope
                    y = .1+.01*dt*slope
                    value = 100+dt*slope
                    candidate.append(dict(level=0,i=i,j=1,k=1,fluid_volume=volume,rho=rho,
                                          rhou=value,rhov=value,rhow=value,rhoE=value,
                                          rho_internal_energy_J_m3=value,temperature=value,rhoY_H2=rho*y))
                    native.append(dict(level=0,x=i-.5,y=.5,z=.5,dx=1,dy=1,dz=1,time=1,
                                       density=rho/1000,xmom=value/10,ymom=value/10,zmom=value/10,
                                       rho_E=value/.1,rho_e=value/.1,Temp=value,rho_H2=rho*y/1000))
                # Covered/non-composite native entries must not affect norms.
                native.append(dict(native[0],x=98.5,density=1.e12))
                app = directory/f'result.csv.step{step}.csv'
                write(app,candidate)
                write(Path(str(app)+'.time.csv'),[dict(step=step,time_s=1)])
                write(directory/f'reference_{step:03}.csv',native)
                runs.append(dict(label=str(dt),dt_s=dt,application=str(directory),reference=str(directory)))
            result = MODULE.analyze(dict(times_s=[1.],runs=runs))
            row = result['rows'][0]
            for solver in ['candidate','native']:
                metrics = row['pairs'][0][solver]['rho']
                self.assertAlmostEqual(metrics['L1'],1.25)
                self.assertAlmostEqual(metrics['L2'],math.sqrt(1.75))
                self.assertAlmostEqual(metrics['Linf'],1.5)
                self.assertAlmostEqual(row['pairs'][0][solver]['Y_H2']['L1'],.0125)
                for field in row['observed_orders'][0][solver].values():
                    for value in field.values():
                        self.assertAlmostEqual(value,1.,places=10)
            clock = Path(runs[0]['application'])/'result.csv.step1.csv.time.csv'
            write(clock,[dict(step=1,time_s=2)])
            with self.assertRaisesRegex(ValueError,'Application clock mismatch'):
                MODULE.analyze(dict(times_s=[1.],runs=runs))


if __name__ == '__main__':
    unittest.main()
