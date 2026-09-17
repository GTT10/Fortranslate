"""Focused MPI EB spray application rank parity and unsupported restart rejection."""
import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe', type=Path, required=True)
    p.add_argument('--mpiexec', type=Path, required=True)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--selected', action='store_true')
    p.add_argument('--scheduled', action='store_true')
    p.add_argument('--restart', action='store_true')
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='eb_spray_', dir=a.work.resolve()))
    template = a.source / 'cases/mpi_eb_spray_3d'
    env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')

    def run(ranks, label, extra='', success=True, decreasing=False, altered=False):
        directory = work / label
        shutil.copytree(template, directory)
        text = (directory / 'tree.nml').read_text()
        if a.scheduled:
            spray = directory / 'spray.nml'
            times = '5.e-9,5.e-9,5.e-9,1.e-8,2.e-8,3.e-8'
            if decreasing:
                times = '5.e-9,4.e-9,5.e-9,1.e-8,2.e-8,3.e-8'
            spray.write_text(spray.read_text().replace('&eb_spray',
                f'&eb_spray\n injection_events=6,injection_times={times},'
                'injection_scales=1.,0.,1.,.5,2.,1.,', 1))
        if a.selected:
            extra += "composition_count=1,composition_species(1)='N2',composition_mole_fractions(1)=1.,"
        (directory / 'tree.nml').write_text(text.replace('&eb_tree', '&eb_tree\n' + extra, 1))
        if altered:
            spray = directory / 'spray.nml'
            spray.write_text(spray.read_text().replace('liquid_temperature=350.', 'liquid_temperature=351.'))
        command = [str(a.mpiexec.resolve()), '-n', str(ranks), str(a.exe.resolve()), 'tree.nml', 'spray.nml']
        result = subprocess.run(command, cwd=directory, env=env, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=90)
        (directory / 'run.log').write_text(result.stdout)
        assert (result.returncode == 0) == success, result.stdout
        return directory

    def rows(path):
        with path.open() as stream:
            return list(csv.DictReader(stream))

    results = [run(n, f'mpi{n}') for n in (2, 4)]
    gas = []
    particles = []
    for directory, ranks in zip(results, (2, 4)):
        gas.append(rows(directory / 'result.csv'))
        files = sorted(directory.glob('result.csv.particles.rank*.csv'))
        assert len(files) == ranks
        drops = sorted((row for file in files for row in rows(file)), key=lambda row: int(row['id']))
        assert [int(row['id']) for row in drops] == list(range(1, 9) if a.scheduled else range(1, 3))
        if a.scheduled:
            assert 'Processed injection events: 5' in (directory / 'run.log').read_text()
            budget = rows(directory / 'result.csv.budget.csv')
            assert abs(float(budget[0]['injected']) - 4.5e-10) < 4.5e-24
            for first, second in zip(drops[:2], drops[2:4]):
                assert {k: v for k, v in first.items() if k != 'id'} == {
                    k: v for k, v in second.items() if k != 'id'}

        for row in drops:
            assert all(math.isfinite(float(v)) for v in row.values())
            assert float(row['mass']) > 0
            if a.scheduled and int(row['id']) >= 7:
                assert float(row['temperature']) == 350 and float(row['v']) == -20, row
                assert float(row['y']) == .001250001, row
            else:
                assert float(row['y']) > .00125 and float(row['v']) > 0, row
                assert float(row['temperature']) != 350, row
            assert float(row['mass']) == float(row['component_mass_1'])
        particles.append(drops)
    for left, right in ((gas[0], gas[1]), (particles[0], particles[1])):
        assert len(left) == len(right) and len(left) > 0
        for first, second in zip(left, right):
            assert first.keys() == second.keys()
            for key in first:
                x, y = float(first[key]), float(second[key])
                assert math.isfinite(x) and math.isfinite(y)
                assert abs(x-y) <= 3e-12 * max(abs(x), abs(y), 1e-30), (key, x, y)
    rejected = run(2, 'restart_rejected', "restart_file='missing.chk',", False)
    assert not (rejected / 'result.csv').exists()
    assert 'Cannot open EB spray manifest' in (rejected / 'run.log').read_text()
    if a.restart:
        stopped = run(2, 'stopped', "checkpoint_file='saved',checkpoint_interval_steps=2,stop_after_checkpoint=.true.,")
        archive = stopped / 'saved.step2'
        assert archive.exists() and Path(str(archive) + '.gas').exists() and Path(str(archive) + '.particles').exists()
        before = {suffix: Path(str(archive) + suffix).read_bytes() for suffix in ('', '.gas', '.particles')}
        overwrite = run(2, 'overwrite', f"checkpoint_file='{stopped / 'saved'}',checkpoint_interval_steps=2,"
                        'stop_after_checkpoint=.true.,', success=False)
        assert 'Cannot reserve EB spray checkpoint' in (overwrite / 'run.log').read_text()
        assert all(Path(str(archive) + suffix).read_bytes() == data for suffix, data in before.items())
        continued = run(4, 'continued', f"restart_file='{archive}',")
        for name in ('result.csv',):
            assert rows(continued / name) == rows(results[1] / name), name
        restored = sorted((r for f in continued.glob('result.csv.particles.rank*.csv') for r in rows(f)),
                          key=lambda r: int(r['id']))
        assert restored == particles[1]
        for first, second in zip(rows(continued / 'result.csv.budget.csv'), rows(results[1] / 'result.csv.budget.csv')):
            for key in ('initial', 'final', 'boundary', 'wall', 'chemistry', 'injected'):
                x, y = float(first[key]), float(second[key])
                assert abs(x-y) <= 3e-12 * max(abs(x), abs(y), 1e-30), (key, x, y)
        mismatch = run(2, 'context_mismatch', f"restart_file='{archive}',", success=False, altered=True)
        assert 'EB spray context mismatch' in (mismatch / 'run.log').read_text()
        payload = Path(str(archive) + '.particles')
        original = payload.read_bytes()
        payload.write_bytes(original[:-40])
        try:
            invalid = run(4, 'truncated', f"restart_file='{archive}',", success=False)
            assert 'Particle snapshot failed' in (invalid / 'run.log').read_text()
            assert not (invalid / 'result.csv').exists()
        finally:
            payload.write_bytes(original)
    if a.scheduled:
        invalid = run(2, 'decreasing', success=False, decreasing=True)
        assert 'Injection times must be nondecreasing' in (invalid / 'run.log').read_text()
        assert not (invalid / 'result.csv').exists()
    print(f'PASS: MPI 2/4 EB gas/particle rank parity, wall reflection, restart rejection; {work}')


if __name__ == '__main__':
    main()
