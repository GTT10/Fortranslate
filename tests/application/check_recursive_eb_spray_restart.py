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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('exe', 'mpiexec', 'source', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
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
        (directory / 'tree.nml').write_text(text.replace('&eb_tree', '&eb_tree\n' + config + controls, 1))
        spray = directory / 'spray.nml'
        spray.write_text(spray.read_text().replace('&eb_spray', '&eb_spray\n'
            'particle_outflow=.true.,injection_events=2,injection_times=0.,1.e-8,injection_scales=1.,.5,', 1))
        command = [str(args.mpiexec.resolve()), '-n', str(ranks), str(args.exe.resolve()), 'tree.nml', 'spray.nml']
        result = subprocess.run(command, cwd=directory, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
        (directory / 'run.log').write_text(result.stdout)
        assert (result.returncode == 0) == success, result.stdout
        return directory

    uninterrupted = run('continuous4', 4)
    stopped = run('checkpoint2', 2, "checkpoint_file='saved',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,")
    archive = stopped / 'saved.step1'
    payload = {suffix: Path(str(archive) + suffix).read_bytes() for suffix in ('', '.gas', '.particles')}
    assert all(b'RECURSIVE_SDC_PHYSICAL_SYNC_V1' in payload[suffix] for suffix in ('.gas', '.particles'))
    assert payload['.gas'].count(b'SDC_REACTION_SOURCE') == 2
    continued = run('continued4', 4, f"restart_file='{archive}',")
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
    gas = rows(continued / 'result.csv')
    assert gas == rows(uninterrupted / 'result.csv'), 'gas/temperature/source-history continuation differs'
    assert any(int(row['level']) == 1 and int(row['node']) == 2 for row in gas)
    left, right = rows(continued / 'result.csv.budget.csv'), rows(uninterrupted / 'result.csv.budget.csv')
    assert len(left) == len(right)
    for first, second in zip(left, right):
        assert first['component'] == second['component']
        for key in ('initial', 'final', 'boundary', 'wall', 'chemistry', 'injected', 'force'):
            x, y = float(first[key]), float(second[key])
            assert math.isfinite(x) and math.isfinite(y)
            assert abs(x - y) <= 3e-12 * max(abs(x), abs(y), 1e-30), (key, x, y)
    print(f'PASS: recursive regrid, pending injection, 2-to-4 restart, common-mode rejection; {work}')


if __name__ == '__main__':
    main()
