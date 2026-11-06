#!/usr/bin/env python3
"""Compare completed diagnostic meshes by conservative fluid-volume restriction.

Reports initial sampling, endpoint and evolution-increment differences separately.
Three grids yield an observed spatial order, not physical or asymptotic validation.
All grids must use identical physical times, time steps and numerical conventions.
"""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fluid_cells(size, domain, plane):
    spacing = [length / size for length in domain]
    volumes = {}
    for level in (0, 1):
        dx = [value / 2**level for value in spacing]
        for i, j, k in itertools.product(range(1, size + 1), repeat=3):
            if level == 0 and max(i, j, k) <= size // 2:
                continue
            fraction = min(1., max(0., (k * dx[2] - plane) / dx[2]))
            if fraction > 0:
                volumes[level, i, j, k] = math.prod(dx) * fraction
    return volumes


def restrict(fine, fine_volumes, coarse_volumes, fields):
    """Project each same-level block of eight fluid children, excluding solids."""
    projected = {}
    used = set()
    for key, volume in coarse_volumes.items():
        level, i, j, k = key
        children = [(level, a, b, c) for a, b, c in itertools.product(
            (2 * i - 1, 2 * i), (2 * j - 1, 2 * j), (2 * k - 1, 2 * k))
            if (level, a, b, c) in fine_volumes]
        child_volume = math.fsum(fine_volumes[c] for c in children)
        if not math.isclose(child_volume, volume, rel_tol=64 * math.ulp(1.), abs_tol=0.):
            raise ValueError(f'Fluid-volume partition mismatch at {key}')
        if used.intersection(children):
            raise ValueError('Overlapping restriction blocks')
        used.update(children)
        projected[key] = {f: math.fsum(fine[c][f] * fine_volumes[c] for c in children) / volume for f in fields}
    if used != fine_volumes.keys() or fine.keys() != fine_volumes.keys():
        raise ValueError('Incomplete composite-fluid restriction')
    residuals = {}
    for field in fields:
        old = math.fsum(fine[k][field] * v for k, v in fine_volumes.items())
        new = math.fsum(projected[k][field] * v for k, v in coarse_volumes.items())
        scale = math.fsum(abs(fine[k][field]) * v for k, v in fine_volumes.items())
        residual = abs(new - old) / max(scale, 1.e-300)
        if not math.isfinite(residual) or residual > 256 * math.ulp(1.):
            raise ValueError(f'Nonconservative restriction: {field}')
        residuals[field] = residual
    return projected, residuals


def admit_original_reports(paths, reports):
    """Reuse an immutable passed baseline explicitly bound by a pinned study.

    Only comparison metadata is derived; input reports and tables stay intact.
    """
    admitted = []
    for path, report in zip(paths, reports):
        if 'mesh_study' in report:
            admitted.append(report)
            continue
        receipt = digest(path / 'checks.json')
        peers = [r for r in reports if r.get('mechanism_basis') == 'original_pinned'
                 and r.get('baseline_checks_sha256') == receipt
                 and r.get('status') == 'diagnostic_passed' and not r['violations']]
        if not peers or report['status'] != 'passed' or report['violations']:
            raise ValueError('Original capture requires a passed study binding its exact report hash')
        if report.get('halo_policy_option') is not None or not report['transport_enabled']:
            raise ValueError('Original capture has nonbaseline controls')
        peer = peers[0]
        for key in ('bundle_sha256', 'source_contract_sha256', 'case_contract_sha256'):
            if report[key] != peer[key]:
                raise ValueError(f'Original capture provenance mismatch: {key}')
        # The final command is the Fortran probe in both checker workflows.
        for entry in (report['commands'][-1], peer['commands'][-1]):
            if entry.get('returncode') != 0:
                raise ValueError('Unfinished original or study probe')
        if report['commands'][-1]['binary_sha256'] != peer['commands'][-1]['binary_sha256']:
            raise ValueError('Original and study probe binaries differ')
        manifest_path = Path(report['reference']) / 'reference_manifest.json'
        if digest(manifest_path) != peer['baseline_reference_manifest_sha256']:
            raise ValueError('Original reference manifest drift')
        manifest = json.loads(manifest_path.read_text())
        if 'mesh_study_sha256' in manifest or manifest['pelec_sha256'] not in peer['native_policy_identity'].values():
            raise ValueError('Original native reference identity mismatch')
        shape = report['contract']['root_shape']
        if len(shape) != 3 or len(set(shape)) != 1:
            raise ValueError('Original comparison requires a cubic root grid')
        admitted.append(dict(report, mesh_study={'cells': shape[0]}, qualification=False,
                             mechanism_basis='original_pinned',
                             native_policy_identity=peer['native_policy_identity'],
                             comparison_admission={'kind': 'reused_passed_original',
                                                   'source_checks_sha256': receipt}))
    return admitted


def validate_reports(reports, temporal=False):
    sizes = [r['mesh_study']['cells'] for r in reports]
    if len(sizes) not in (2, 3) or any(n not in (8, 16, 32) for n in sizes):
        raise ValueError('Need two or three supported grids')
    if temporal:
        times = [r['contract']['dt_s'] for r in reports]
        if len(set(sizes)) != 1 or any(not math.isfinite(t) or t <= 0 for t in times):
            raise ValueError('Temporal comparison requires one fixed mesh and positive time steps')
        if any(not math.isclose(a, 2 * b, rel_tol=1.e-14) for a, b in zip(times, times[1:])):
            raise ValueError('Temporal comparison requires factor-two time refinement')
        endpoints = [r['contract']['steps'] * t for r, t in zip(reports, times)]
        if any(not math.isclose(t, endpoints[0], rel_tol=1.e-14) for t in endpoints):
            raise ValueError('Temporal comparison requires identical endpoint times')
    elif any(b != 2 * a for a, b in zip(sizes, sizes[1:])):
        raise ValueError('Need two or three consecutive factor-two grids')
    first = reports[0]
    for report in reports:
        reused = report.get('comparison_admission', {}).get('kind') == 'reused_passed_original'
        if report['status'] != ('passed' if reused else 'diagnostic_passed') or report['qualification'] is not False:
            raise ValueError('Only completed, passing diagnostic captures are comparable')
        if report.get('mechanism_basis', 'diagnostic_activation') != first.get('mechanism_basis', 'diagnostic_activation'):
            raise ValueError('Mixed original and diagnostic mechanism bases')
        size = report['mesh_study']['cells']
        if report['contract']['root_shape'] != [size] * 3 or report['contract']['fine_coarse_box'] != [[1] * 3, [size // 2] * 3]:
            raise ValueError('Inconsistent hierarchy metadata')
        for key in ('dt_s', 'steps', 'domain_m', 'plane_z_m', 'input_sha256', 'acceptance'):
            if temporal and key in ('dt_s', 'steps'):
                continue
            if report['contract'][key] != first['contract'][key]:
                raise ValueError(f'Mixed physical/time controls: {key}')
        for key in ('bundle_sha256', 'native_policy_identity', 'simple_transport_source', 'initial_reaction_history'):
            if report[key] != first[key]:
                raise ValueError(f'Mixed numerical conventions: {key}')
        # Old reports predate explicit tolerance metadata and used V1/V2 defaults.
        legacy_reaction = {'relative': 1.e-10, 'fortran_absolute_mass_fraction': 1.e-14,
                           'fortran_atol_basis': 'mass_fraction_times_local_density',
                           'native_absolute_density_g_cm3': 1.e-18,
                           'native_absolute_temperature_K': 1.e-18}
        if report.get('reaction_controls', legacy_reaction) != first.get('reaction_controls', legacy_reaction):
            raise ValueError('Mixed chemical integration tolerances')
        legacy_eos = {'fortran_converged_inverse': False, 'native_converged_inverse': False}
        if report.get('eos_controls', legacy_eos) != first.get('eos_controls', legacy_eos):
            raise ValueError('Mixed caloric inversion precision')
        legacy_velocity = {'tg_amplitude_cm_s': 1000.}
        if report.get('initial_velocity_controls', legacy_velocity) != first.get('initial_velocity_controls', legacy_velocity):
            raise ValueError('Mixed initial velocity controls')
        legacy_profile = {'kind': 'frozen_table', 'points': 33}
        if report.get('initial_profile_controls', legacy_profile) != first.get('initial_profile_controls', legacy_profile):
            raise ValueError('Mixed initial profile sampling')
        profile_hash = report.get('initial_profile_sha256', report['contract']['input_sha256'].get('profile.csv'))
        first_hash = first.get('initial_profile_sha256', first['contract']['input_sha256'].get('profile.csv'))
        if profile_hash != first_hash:
            raise ValueError('Mixed initial profile data')
        if len(report['steps']) != report['contract']['steps'] + 1 or report['violations']:
            raise ValueError('Incomplete or failing field capture')
    return sizes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', nargs='+', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--temporal', action='store_true',
                        help='Compare fixed-mesh, factor-two time refinement at one physical endpoint')
    args = parser.parse_args()
    def sort_key(path):
        report = json.loads((path / 'checks.json').read_text())
        return -report['contract']['dt_s'] if args.temporal else report['contract']['root_shape'][0]
    paths = sorted((p.resolve() for p in args.runs), key=sort_key)
    reports = [json.loads((p / 'checks.json').read_text()) for p in paths]
    reports = admit_original_reports(paths, reports)
    sizes = validate_reports(reports, temporal=args.temporal)
    first = reports[0]
    fields = list(first['steps'][0]['fields'])
    domain, plane = first['contract']['domain_m'], first['contract']['plane_z_m']
    steps, dt = first['contract']['steps'], first['contract']['dt_s']
    volumes = [fluid_cells(n, domain, plane) for n in sizes]
    input_hashes = {str(p / 'checks.json'): digest(p / 'checks.json') for p in paths}

    def load(directory, size, cells, step, native, timestep):
        path = directory / f'{"reference" if native else "fortran"}_{step:03}.csv'
        input_hashes[str(path)] = digest(path)
        data, seen = {}, set()
        with path.open() as stream:
            for row in csv.DictReader(stream):
                if not all(math.isfinite(float(v)) for v in row.values()):
                    raise ValueError(f'Nonfinite snapshot: {path}')
                if abs(float(row['time' if native else 'time_s']) - step * timestep) > 1.e-19:
                    raise ValueError(f'Incorrect physical time: {path}')
                level = int(row['level'])
                index = tuple(round(float(row[a]) / float(row['d' + a]) - .5) + 1 for a in 'xyz') if native else tuple(int(row[a]) for a in 'ijk')
                key = (level, *index)
                if level not in (0, 1) or any(i < 1 or i > size for i in index) or key in seen:
                    raise ValueError(f'Invalid or duplicate cell: {key}')
                seen.add(key)
                dx = [length / size / 2**level for length in domain]
                fraction = min(1., max(0., (index[2] * dx[2] - plane) / dx[2]))
                if abs(float(row['vfrac']) - fraction) > 1.e-13:
                    raise ValueError(f'Incorrect EB fraction: {key}')
                if native and any(abs(float(row['d' + a]) * .01 / dx[d] - 1) > 1.e-12 for d, a in enumerate('xyz')):
                    raise ValueError('Incorrect native cell spacing')
                if not native and int(row['node']) != level + 1:
                    raise ValueError('Incorrect Fortran node')
                if key in cells:
                    data[key] = {f: float(row[f]) * (1. if not native or f == 'Temp' else
                        .1 if f in ('rho_E', 'rho_e') else 10. if f in ('xmom', 'ymom', 'zmom') else 1000.) for f in fields}
        if len(seen) != 2 * size**3 or data.keys() != cells.keys():
            raise ValueError('Incomplete mesh snapshot')
        return data

    result = {'qualification': False, 'scope': 'common-dt fixed-time spatial differences; observed order only',
              'sizes': sizes, 'dt_s': dt, 'endpoint_s': steps * dt,
              'temperature_restriction': 'fluid-volume mean of stored temperature, not EOS reconstruction',
              'checker_sha256': digest(Path(__file__)), 'solvers': {}}
    if args.temporal:
        result.update(scope='fixed-mesh fixed-endpoint temporal differences; observed order only',
                      timesteps_s=[r['contract']['dt_s'] for r in reports])
    native_label = 'native_pinned' if first.get('mechanism_basis') == 'original_pinned' else 'native_held'
    result['mechanism_basis'] = first.get('mechanism_basis', 'diagnostic_activation')
    for solver, native in (('fortran', False), (native_label, True)):
        states = [(load(p, n, v, 0, native, r['contract']['dt_s']),
                   load(p, n, v, r['contract']['steps'], native, r['contract']['dt_s']))
                  for p, n, v, r in zip(paths, sizes, volumes, reports)]
        if args.temporal and any(s[0] != states[0][0] for s in states[1:]):
            raise ValueError('Temporal comparison requires identical initial fields')
        pairs = []
        for index in range(len(sizes) - 1):
            coarse, fine = states[index], states[index + 1]
            coarse_volumes = volumes[index]
            volume = math.fsum(coarse_volumes.values())
            projections = ([(s, {f: 0. for f in fields}) for s in fine] if args.temporal else
                           [restrict(s, volumes[index + 1], coarse_volumes, fields) for s in fine])
            pair = {'coarse_cells': sizes[index], 'fine_cells': sizes[index + 1],
                    'restriction_scaled_residual': [p[1] for p in projections], 'norms': {}}
            if args.temporal:
                pair['timesteps_s'] = result['timesteps_s'][index:index + 2]
            for stage in ('initial', 'endpoint', 'change_from_initial'):
                norms = {}
                for field in fields:
                    if stage == 'change_from_initial':
                        errors = {k: (coarse[1][k][field] - coarse[0][k][field]) -
                                  (projections[1][0][k][field] - projections[0][0][k][field]) for k in coarse_volumes}
                    else:
                        step = int(stage == 'endpoint')
                        errors = {k: coarse[step][k][field] - projections[step][0][k][field] for k in coarse_volumes}
                    norms[field] = {'volume_rms': math.sqrt(math.fsum(e**2 * coarse_volumes[k] for k, e in errors.items()) / volume),
                                    'maximum_absolute': max(abs(e) for e in errors.values())}
                pair['norms'][stage] = norms
            pairs.append(pair)
        orders = {}
        maximum_orders = {}
        if len(pairs) == 2:
            for stage in pairs[0]['norms']:
                orders[stage] = {}
                maximum_orders[stage] = {}
                for field in fields:
                    a, b = [p['norms'][stage][field]['volume_rms'] for p in pairs]
                    scale = max(abs(s[k][field]) for levels in states for s in levels for k in s)
                    orders[stage][field] = math.log2(a / b) if min(a, b) > 128 * math.ulp(max(scale, 1.e-30)) else None
                    amax, bmax = [p['norms'][stage][field]['maximum_absolute'] for p in pairs]
                    maximum_orders[stage][field] = (math.log2(amax / bmax)
                        if min(amax, bmax) > 128 * math.ulp(max(scale, 1.e-30)) else None)
        result['solvers'][solver] = {'pairs': pairs, 'observed_orders': orders,
                                     'observed_maximum_orders': maximum_orders}
    result['reused_original_captures'] = {str(p): r['comparison_admission']
        for p, r in zip(paths, reports) if 'comparison_admission' in r}
    result['input_sha256'] = input_hashes
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(args.output)


if __name__ == '__main__':
    main()
