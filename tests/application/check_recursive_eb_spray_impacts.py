"""Two/four-rank recursive EB splash and newborn-film application check."""
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
    parser.add_argument('--restart', action='store_true')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='impact_', dir=args.work.resolve()))
    outputs = []
    env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')
    archive = work / 'stopped' / 'saved.step1'
    checkpoint_bytes = {}
    budgets = []
    runs = [(2, '2', ''), (4, '4', '')]
    if args.restart:
        runs = [(4, 'continuous', ''),
                (2, 'stopped', "checkpoint_file='saved',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,"),
                (4, 'continued', f"restart_file='{archive}',")]
    for ranks, label, controls in runs:
        directory = work / label
        shutil.copytree(args.source / 'cases/mpi_eb_spray_3d', directory)
        tree = directory / 'tree.nml'
        text = tree.read_text()
        if not args.restart:
            text = text.replace('final_time=2.e-8', 'final_time=1.e-8')
        text = text.replace("riemann_solver='hllc'", "riemann_solver='pelec',hydro_method='plm_ctu'")
        text = text.replace("output_file='result.csv'", "output_file='result.csv',time_integrator='sdc',"
                            "sdc_iterations=2,subcycling=.true.,boundary_faces=6*'outflow'")
        tree.write_text(text.replace("&eb_tree", "&eb_tree\n" + controls, 1))
        (directory / 'spray.nml').write_text("""&eb_spray
 nc=1,number=1,vapor_names(1)='H2O',fractions(1)=1.,
 total_mass=1.e-9,center=.0012,.001250001,.0015,
 axis=0.,-1.,0.,speed=100.,diameter=1.e-4,liquid_temperature=350.,
 particle_outflow=.true.,wall_impacts=.true.,
 particle_options%surface_tension=.072,particle_options%liquid_viscosity=.001,
 walls%enabled=6*.true.,walls%temperature=6*350.,walls%evolve=.true.,
 walls%liquid_conductivity=.6,walls%minimum_height=1.e-15
/
""")
        with (directory / 'run.log').open('w') as log:
            subprocess.run([str(args.mpiexec.resolve()), '-n', str(ranks), str(args.exe.resolve()),
                            'tree.nml', 'spray.nml'], cwd=directory, env=env, stdout=log,
                           stderr=subprocess.STDOUT, timeout=120, check=True)
        if label == 'stopped':
            checkpoint_bytes = {suffix: Path(str(archive) + suffix).read_bytes()
                                for suffix in ('', '.gas', '.particles')}
            metadata = checkpoint_bytes['.particles'].decode().splitlines()[2].split()
            assert int(metadata[3]) > 0 and int(metadata[4]) > 1, metadata
            assert 'Completed steps: 1' in (directory / 'run.log').read_text()
            continue
        gas = rows(directory / 'result.csv')
        particles = sorted((row for path in directory.glob('result.csv.particles.rank*.csv')
                            for row in rows(path)), key=lambda row: int(row['id']))
        films = sorted((row for path in directory.glob('result.csv.films.rank*.csv')
                        for row in rows(path)), key=lambda row: tuple(float(v) for v in row.values()))
        assert len(particles) > 1 and films, (len(particles), len(films))
        ids = [int(row['id']) for row in particles]
        assert len(set(ids)) == len(ids) and max(ids) > 1, ids
        for row in particles + films:
            assert all(math.isfinite(float(v)) for v in row.values()), row
            assert float(row['mass']) > 0, row
            mass = sum(float(v) for key, v in row.items() if key.startswith('component_mass_'))
            assert abs(mass - float(row['mass'])) <= 1e-14 * mass, row
        budget = rows(directory / 'result.csv.budget.csv')
        mass_scale = max(abs(float(budget[0][key])) for key in ('initial', 'final', 'injected'))
        energy_scale = max(abs(float(budget[4][key])) for key in ('initial', 'final', 'injected'))
        maximum = 0.0
        for row in budget:
            residual = float(row['final']) - float(row['initial'])
            residual -= sum(float(row[key]) for key in ('boundary', 'wall', 'chemistry', 'injected', 'force'))
            scale = max(abs(float(row[key])) for key in ('initial', 'final', 'injected'))
            if 2 <= int(row['component']) <= 4:
                scale = max(scale, math.sqrt(mass_scale * energy_scale))
            scale = max(scale, 1e-30)
            assert abs(residual) <= 3e-12 * scale, (row, residual)
            maximum = max(maximum, abs(residual) / scale)
        budgets.append(budget)
        print(f'{label}, {ranks} ranks: {len(particles)} particles, {len(films)} films, residual {maximum:.17g}')
        outputs.append((gas, particles, films))
    for left, right in zip(*outputs):
        assert len(left) == len(right)
        for first, second in zip(left, right):
            assert first.keys() == second.keys()
            for key in first:
                x, y = float(first[key]), float(second[key])
                assert math.isfinite(x) and math.isfinite(y)
                assert abs(x - y) <= 3e-12 * max(abs(x), abs(y), 1e-30), (key, x, y)
    if args.restart:
        for left, right in zip(*outputs):
            assert left == right, 'restart gas/particle/film CSVs must exactly match continuous run'
        assert all(Path(str(archive) + suffix).read_bytes() == data for suffix, data in checkpoint_bytes.items())
        for first, second in zip(*budgets):
            for key in ('initial', 'final', 'boundary', 'wall', 'chemistry', 'injected', 'force'):
                x, y = float(first[key]), float(second[key])
                assert abs(x - y) <= 3e-12 * max(abs(x), abs(y), 1e-30), (key, x, y)
        assert 'Completed steps: 2' in (work / 'continued' / 'run.log').read_text()
    print(f'PASS: recursive EB splash, films, inventories, rank comparison; restart={args.restart}; {work}')


if __name__ == '__main__':
    main()
