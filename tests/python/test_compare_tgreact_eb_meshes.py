"""Analytic cut-cell projection and scientific-control checks."""
import copy
import importlib.util
import math
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

PATH = Path(__file__).resolve().parents[2] / 'tools/compare_tgreact_eb_meshes.py'
spec = importlib.util.spec_from_file_location('mesh_comparison', PATH)
mesh = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mesh)


class MeshComparisonTests(unittest.TestCase):
    def test_bounded_initial_velocity_controls(self):
        sys.path.insert(0, str(PATH.parent))
        try:
            from check_tgreact_eb_amr_reference import velocity_controls
            self.assertEqual(velocity_controls(None), {'tg_amplitude_cm_s': 1000.})
            for amplitude in (0, 1, 1000):
                self.assertEqual(velocity_controls({'initial_velocity_cm_s': amplitude}),
                                 {'tg_amplitude_cm_s': float(amplitude)})
            for invalid in (True, -1, 2, 1001, float('nan'), '0'):
                with self.assertRaises(AssertionError):
                    velocity_controls({'initial_velocity_cm_s': invalid})
        finally:
            sys.path.pop(0)

    def test_temporal_order_and_control_refusal(self):
        first = {'status': 'diagnostic_passed', 'qualification': False, 'violations': [],
            'steps': [{'fields': {'density': {}}}] * 3, 'mesh_study': {'cells': 8},
            'contract': {'steps': 2, 'dt_s': .5e-7, 'domain_m': [.01] * 3, 'plane_z_m': .0033,
                'root_shape': [8] * 3, 'fine_coarse_box': [[1] * 3, [4] * 3], 'input_sha256': {}, 'acceptance': {}},
            'bundle_sha256': 'same', 'native_policy_identity': {}, 'simple_transport_source': {},
            'initial_reaction_history': 'computed'}
        reports = []
        for factor in (1, 2, 4):
            report = copy.deepcopy(first)
            report['contract']['dt_s'] /= factor
            report['contract']['steps'] *= factor
            report['steps'] = [first['steps'][0]] * (2 * factor + 1)
            reports.append(report)
        self.assertEqual(mesh.validate_reports(reports, temporal=True), [8] * 3)
        for kind in ('endpoint', 'ratio', 'mesh', 'mechanism'):
            bad = copy.deepcopy(reports)
            if kind == 'endpoint':
                bad[1]['contract']['steps'] += 1
            elif kind == 'ratio':
                bad[1]['contract']['dt_s'] *= 1.5
            elif kind == 'mesh':
                bad[1]['mesh_study']['cells'] = 16
            else:
                bad[1]['bundle_sha256'] = 'different'
            with self.assertRaises(ValueError):
                mesh.validate_reports(bad, temporal=True)
        # Manufactured endpoint error proportional to dt**2 must give order 2.
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = []
            for factor, report in zip((1, 2, 4), reports):
                path = root / str(factor)
                path.mkdir()
                paths.append(path)
                (path / 'checks.json').write_text(json.dumps(report))
                steps = report['contract']['steps']
                for native in (False, True):
                    for step in (0, steps):
                        filename = f'{"reference" if native else "fortran"}_{step:03}.csv'
                        fields = ['level', 'vfrac', 'density'] + (['time', 'x', 'y', 'z', 'dx', 'dy', 'dz']
                            if native else ['time_s', 'node', 'i', 'j', 'k'])
                        with (path / filename).open('w') as stream:
                            writer = csv.DictWriter(stream, fieldnames=fields)
                            writer.writeheader()
                            for level in (0, 1):
                                dx = .01 / 8 / 2**level
                                for i, j, k in mesh.itertools.product(range(1, 9), repeat=3):
                                    row = {'level': level, 'vfrac': min(1., max(0., (k * dx - .0033) / dx)),
                                           'density': (2. if step == 0 else 3. + .01 / factor**2) / (1000 if native else 1)}
                                    if native:
                                        row.update(time=step * report['contract']['dt_s'])
                                        for axis, index in zip('xyz', (i, j, k)):
                                            row[axis] = (index - .5) * dx * 100
                                            row['d' + axis] = dx * 100
                                    else:
                                        row.update(time_s=step * report['contract']['dt_s'], node=level + 1, i=i, j=j, k=k)
                                    writer.writerow(row)
            output = root / 'comparison.json'
            subprocess.run([sys.executable, str(PATH), '--temporal', '--runs', *map(str, paths),
                            '--output', str(output)], check=True, capture_output=True)
            result = json.loads(output.read_text())
            for solver in result['solvers'].values():
                self.assertAlmostEqual(solver['observed_orders']['endpoint']['density'], 2., places=10)
                self.assertIsNone(solver['observed_orders']['initial']['density'])

    def test_cut_cell_constant_and_linear_volume_means(self):
        # The exact fluid centroid differs from the Cartesian cell centre in
        # the cut row: unweighted averaging would fail this linear-field case.
        domain, plane = [.01] * 3, .0033
        coarse = mesh.fluid_cells(8, domain, plane)
        fine = mesh.fluid_cells(16, domain, plane)

        def exact(key, size):
            dz = domain[2] / size / 2**key[0]
            z = .5 * (max((key[3] - 1) * dz, plane) + key[3] * dz)
            return {'density': 2., 'xmom': z - .5 * (plane + domain[2])}

        data = {k: exact(k, 16) for k in fine}
        projected, residual = mesh.restrict(data, fine, coarse, ('density', 'xmom'))
        self.assertAlmostEqual(math.fsum(coarse.values()), .01**2 * (.01 - plane), delta=1.e-21)
        self.assertEqual(len(coarse), 544)
        self.assertEqual(len(fine), 4160)
        for key in coarse:
            self.assertAlmostEqual(projected[key]['density'], 2., delta=2.e-14)
            self.assertAlmostEqual(projected[key]['xmom'], exact(key, 8)['xmom'], delta=2.e-17)
        self.assertLess(max(residual.values()), 1.e-14)
        del data[next(iter(data))]
        with self.assertRaises((KeyError, ValueError)):
            mesh.restrict(data, fine, coarse, ('density', 'xmom'))

    def test_refuse_mixed_time_step_and_mechanism(self):
        first = {'status': 'diagnostic_passed', 'qualification': False, 'violations': [],
            'steps': [{}, {}, {}], 'mesh_study': {'cells': 8},
            'contract': {'steps': 2, 'dt_s': 5.e-8, 'domain_m': [.01] * 3, 'plane_z_m': .0033,
                'root_shape': [8] * 3, 'fine_coarse_box': [[1] * 3, [4] * 3], 'input_sha256': {}, 'acceptance': {}},
            'bundle_sha256': 'same', 'native_policy_identity': {}, 'simple_transport_source': {},
            'initial_reaction_history': 'computed'}
        second = copy.deepcopy(first)
        second['mesh_study']['cells'] = 16
        second['contract'].update(root_shape=[16] * 3, fine_coarse_box=[[1] * 3, [8] * 3])
        self.assertEqual(mesh.validate_reports([first, second]), [8, 16])
        for kind in ('time', 'mechanism', 'reaction', 'eos', 'profile', 'profile_hash', 'velocity'):
            bad = copy.deepcopy(second)
            if kind == 'time':
                bad['contract']['dt_s'] *= 2
            elif kind == 'mechanism':
                bad['bundle_sha256'] = 'different'
            elif kind == 'reaction':
                bad['reaction_controls'] = {'relative': 1.e-12}
            elif kind == 'eos':
                bad['eos_controls'] = {'fortran_converged_inverse': True, 'native_converged_inverse': False}
            elif kind == 'profile':
                bad['initial_profile_controls'] = {'kind': 'cosine_sampling_diagnostic', 'points': 129}
            elif kind == 'velocity':
                bad['initial_velocity_controls'] = {'tg_amplitude_cm_s': 0.}
            else:
                bad['initial_profile_sha256'] = 'different'
            with self.assertRaises(ValueError):
                mesh.validate_reports([first, bad])


if __name__ == '__main__':
    unittest.main()
