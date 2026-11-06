#!/usr/bin/env python3
"""Connect saved initial shear halos to observed fine spatial-stage fluxes.

Read-only diagnosis of the 16-root-cell, L=0.01 m, ratio-two shear fixture.
The analytic mesh table is a local interpolation estimate, not CFD convergence.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial', type=Path, required=True)
    parser.add_argument('--stages', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    identities = {}

    def table(path):
        identities[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        return read(path)

    observations = []
    for policy in ('default', 'momentum', 'velocity', 'paired'):
        initial = args.initial / f'cut_{policy}_initial_pressure.csv'
        pressure = {tuple(int(r[k]) for k in ('node', 'i', 'j', 'k')): r for r in table(initial)}
        flux = {tuple(int(r[k]) for k in ('node', 'axis', 'i', 'j', 'k')): r
                for r in table(args.initial / f'cut_{policy}_initial_flux.csv')}
        stages = table(args.stages / f'cut_{policy}_stage_flux.csv')
        # Later probes also record coarse donors and neighboring fine cells.
        # This audit only links the two original fine-cell intersections to
        # the saved initial halo stencil; do not count unrelated probes.
        first = [r for r in stages if int(r['stage']) == 0
                 and float(r['time_start_s']) == 0 and int(r['probe']) in (1, 2)]
        assert len(first) == 12, 'Requires two six-face fine probes at time zero'
        for probe, cell in ((1, (1, 8, 11)), (2, (1, 9, 16))):
            subset = [r for r in first if int(r['probe']) == probe]
            assert len(subset) == 6
            assert {(int(r['axis']), int(r['side'])) for r in subset} == {
                (axis, side) for axis in (1, 2, 3) for side in (0, 1)}
            maximum = 0.0
            for r in subset:
                assert tuple(int(r[k]) for k in ('i', 'j', 'k')) == cell
                face = list(cell)
                axis = int(r['axis'])
                face[axis - 1] -= 1 - int(r['side'])
                old = flux[(2, axis, *face)]
                for c, name in (('mx', 'x_momentum_Pa'), ('my', 'y_momentum_Pa'), ('mz', 'z_momentum_Pa')):
                    old_value, new_value = float(old[name]), float(r['centroid_' + c])
                    relative = abs(old_value - new_value) / max(1, abs(old_value), abs(new_value))
                    assert relative < 1e-12, (policy, probe, axis, c, relative)
                    maximum = max(maximum, relative)
            ghost = pressure[(2, 0, cell[1], cell[2])]
            valid = pressure[(2, *cell)]
            coarse_j, coarse_k = (cell[1] + 1) // 2, (cell[2] + 1) // 2
            coarse = [float(pressure[(1, 16, coarse_j + j, coarse_k)]['u_m_s']) for j in (-1, 0, 1)]
            dm, dp = coarse[1] - coarse[0], coarse[2] - coarse[1]
            centered = 0.5 * (dm + dp)
            limited = math.copysign(min(abs(centered), 2 * abs(dm), 2 * abs(dp)), centered) if dm * dp >= 0 else 0.0
            offset = 0.25 if cell[1] % 2 == 0 else -0.25
            predicted = coarse[1] + offset * limited
            assert abs(float(ghost['u_m_s']) - predicted) < 1e-12
            observations.append({'policy': policy, 'probe': probe, 'coarse_velocity_stencil_m_s': coarse,
                                 'centered_slope_m_s': centered, 'limited_slope_m_s': limited,
                                 'ghost_velocity_m_s': float(ghost['u_m_s']),
                                 'fine_velocity_m_s': float(valid['u_m_s']),
                                 'velocity_jump_m_s': float(valid['u_m_s']) - float(ghost['u_m_s']),
                                 'pressure_jump_Pa': float(valid['pressure_Pa']) - float(ghost['pressure_Pa']),
                                 'maximum_initial_face_relative_difference': maximum})
    analytic = []
    for cells in (8, 16, 32):
        width = 2 * math.pi / cells
        j = cells // 4
        coarse = 10 * math.sin((j - 0.5) * width) * math.sin(width / 2) / (width / 2)
        fine = 10 * math.sin((cells / 2 - 0.5) * width / 2) * math.sin(width / 4) / (width / 4)
        analytic.append({'root_cells': cells, 'coarse_average_m_s': coarse,
                         'fine_average_m_s': fine, 'constant_donor_error_m_s': fine - coarse})
    output = {'status': 'passed', 'qualification': False, 'files_sha256': identities,
              'observations': observations, 'analytic_cell_average_estimate_only': analytic,
              'scope': 'Local initial halo attribution; no numerical policy change or new mesh trajectory'}
    args.output.write_text(json.dumps(output, indent=2) + '\n')
    print('PASS: 8 initial halo probes; 48 faces match current stages; limited donor velocities reproduced')


if __name__ == '__main__':
    main()
