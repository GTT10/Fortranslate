"""Recursive SDC: regrid, pending injection, and 2-to-4-rank restart."""
import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def heat_snapshot(data):
    lines = data.decode().splitlines()
    header = lines.index('PELEF_EB_TREE_CHECKPOINT_V6')
    count = int(lines[header + 1])
    shapes = [tuple(map(int, lines[header + 2 + node].split()[2:5])) for node in range(count)]
    markers = [index for index, line in enumerate(lines) if line == 'SDC_HEAT_RELEASE']
    assert len(markers) == count
    blocks = {}
    for node, (shape, index) in enumerate(zip(shapes, markers), 1):
        cells = math.prod(shape)
        values = tuple(float(value) for value in lines[index + 1:index + 1 + cells])
        assert len(values) == cells and all(math.isfinite(value) for value in values)
        blocks[node] = (shape, values)
    return blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('exe', 'mpiexec', 'source', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--selected', action='store_true')
    parser.add_argument('--rotating', action='store_true', help='Exercise MOL rotating coordinate energy history')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='recursive_restart_', dir=args.work.resolve()))
    template = args.source / 'cases/mpi_eb_spray_3d'
    env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')

    def run(name, ranks, controls='', recursive=True, success=True):
        directory = work / name
        shutil.copytree(template, directory)
        text = (directory / 'tree.nml').read_text().replace("riemann_solver='hllc'", "riemann_solver='pelec'")
        config = ("time_integrator='sdc',sdc_iterations=2,hydro_method='plm_ctu',"
                  f"subcycling={'.true.' if recursive else '.false.'},"
                  "boundary_faces=6*'outflow',external_forcing=100.,-200.,300.,"
                  'regrid_count=1,regrid_times(1)=1.e-8,regrid_child_count(1)=1,'
                  'regrid_lower(:,1,1)=3,1,1,regrid_upper(:,1,1)=4,4,4,')
        if args.rotating:
            config = config.replace("hydro_method='plm_ctu'",
                                    "hydro_method='plm_ctu',mol_mode=.true.,rotating_frame_enabled=.true.,"
                                    "rotation_axis=2,rotation_angular_speed=1000.,rotation_axis_point=.001,.001,.001")
        if args.selected:
            config += ("composition_count=1,composition_species(1)='N2',"
                       "composition_mole_fractions(1)=1.,")
        (directory / 'tree.nml').write_text(text.replace('&eb_tree', '&eb_tree\n' + config + controls, 1))
        spray = directory / 'spray.nml'
        spray.write_text(spray.read_text().replace('&eb_spray', '&eb_spray\n'
            'particle_outflow=.true.,injection_events=2,injection_times=0.,1.e-8,injection_scales=1.,.5,', 1))
        command = [str(args.mpiexec.resolve()), '-n', str(ranks), str(args.exe.resolve()), 'tree.nml', 'spray.nml']
        try:
            result = subprocess.run(command, cwd=directory, env=env, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    timeout=180 if args.rotating else 90)
        except subprocess.TimeoutExpired as failure:
            output = failure.stdout or ''
            if isinstance(output, bytes):
                output = output.decode(errors='replace')
            (directory / 'run.log').write_text(output)
            raise
        (directory / 'run.log').write_text(result.stdout)
        assert (result.returncode == 0) == success, result.stdout
        return directory

    final_snapshot = "checkpoint_file='final',checkpoint_interval_steps=2,"
    uninterrupted = run('continuous4', 4, final_snapshot)
    stopped = run('checkpoint2', 2, "checkpoint_file='saved',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,")
    archive = stopped / 'saved.step1'
    payload = {suffix: Path(str(archive) + suffix).read_bytes() for suffix in ('', '.gas', '.particles')}
    assert all(b'RECURSIVE_SDC_PHYSICAL_SYNC_V1' in payload[suffix] for suffix in ('.gas', '.particles'))
    assert payload['.gas'].count(b'SDC_REACTION_SOURCE') == 2
    if args.rotating:
        assert all(b'ROTATING_PARTICLE_ENERGY_V1' in payload[suffix] for suffix in ('.gas', '.particles'))

    heat_snapshot(payload['.gas'])
    continued = run('continued4', 4, f"restart_file='{archive}'," + final_snapshot)
    rejected = run('common_mode_rejected', 2, f"restart_file='{archive}',", recursive=False, success=False)
    assert 'Restart subcycling mode mismatch' in (rejected / 'run.log').read_text()
    assert not (rejected / 'result.csv').exists()
    assert all(Path(str(archive) + suffix).read_bytes() == data for suffix, data in payload.items())
    for directory in (uninterrupted, stopped):
        assert 'Regrid stage: 1 nodes: 2' in (directory / 'run.log').read_text()
    assert 'Processed injection events: 1' in (stopped / 'run.log').read_text()
    for directory in (uninterrupted, continued):
        log = (directory / 'run.log').read_text()
        assert 'Processed injection events: 2' in log and 'Completed steps: 2' in log

    def particles(directory):
        return sorted((row for path in directory.glob('result.csv.particles.rank*.csv') for row in rows(path)),
                      key=lambda row: int(row['id']))

    actual, expected = particles(continued), particles(uninterrupted)
    assert [int(row['id']) for row in actual] == [1, 2, 3, 4]
    assert actual == expected, 'particle payload differs after changed-rank restart'
    actual_heat = heat_snapshot((continued / 'final.step2.gas').read_bytes())
    expected_heat = heat_snapshot((uninterrupted / 'final.step2.gas').read_bytes())
    assert actual_heat == expected_heat, 'accepted stage heat differs after changed-rank restart'
    assert any(value != 0 for shape, values in actual_heat.values() for value in values)
    gas = rows(continued / 'result.csv')
    for row in gas:
        node, i, j, k = (int(row[key]) for key in ('node', 'i', 'j', 'k'))
        heat = float(row['heat_release_W_m3'])
        assert math.isfinite(heat)
        shape, values = actual_heat[node]
        nx, ny, nz = shape
        assert 1 <= i <= nx and 1 <= j <= ny and 1 <= k <= nz
        assert heat == values[(k - 1) * nx * ny + (j - 1) * nx + i - 1], row

    assert gas == rows(uninterrupted / 'result.csv'), 'gas/temperature/source-history continuation differs'
    assert any(int(row['level']) == 1 and int(row['node']) == 2 for row in gas)
    left, right = rows(continued / 'result.csv.budget.csv'), rows(uninterrupted / 'result.csv.budget.csv')
    assert len(left) == len(right)
    for first, second in zip(left, right):
        assert first['component'] == second['component']
        keys = ('initial', 'final', 'boundary', 'wall', 'chemistry', 'injected', 'force')
        if args.rotating:
            keys += ('coordinate_energy',)
        for key in keys:
            x, y = float(first[key]), float(second[key])
            assert math.isfinite(x) and math.isfinite(y)
            assert abs(x - y) <= 3e-12 * max(abs(x), abs(y), 1e-30), (key, x, y)
    if args.rotating:
        assert any(float(row['coordinate_energy']) != 0 for row in left)
        assert all(float(row['coordinate_energy']) == 0 for row in left if int(row['component']) != 5)
        for row in left + right:
            expected_residual = (float(row['final']) - float(row['initial']) - float(row['boundary'])
                                 - float(row['wall']) - float(row['chemistry']) - float(row['injected'])
                                 - float(row['force']) - float(row['coordinate_energy']))
            assert abs(float(row['residual']) - expected_residual) <= 1e-24, row
        # The new history tag is required for both members of the transaction.
        # A pre-ledger archive cannot silently acquire a zero cumulative value.
        stale_base = stopped / 'pre_ledger.step1'
        for suffix, data in payload.items():
            Path(str(stale_base) + suffix).write_bytes(data.replace(b'ROTATING_PARTICLE_ENERGY_V1|', b''))
        stale = run('pre_ledger_rejected', 2, f"restart_file='{stale_base}',", success=False)
        assert 'EB spray context mismatch' in (stale / 'run.log').read_text()
        assert not (stale / 'result.csv').exists()
    print(f'PASS: recursive regrid, pending injection, 2-to-4 restart, common-mode rejection; {work}')


if __name__ == '__main__':
    main()
