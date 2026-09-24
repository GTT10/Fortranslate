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
    p.add_argument('--rebalance', action='store_true')
    p.add_argument('--regrid', action='store_true')
    p.add_argument('--ports', action='store_true')
    p.add_argument('--profiles', action='store_true')
    p.add_argument('--mach', action='store_true')
    p.add_argument('--pelemp-cap', action='store_true')
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
        if a.rebalance:
            text = text.replace('child_count=1,parent(1)=1,ratio(1)=2,lower(:,1)=1,1,1,upper(:,1)=2,4,4,',
                'child_count=4,parent=1,1,1,1,ratio=2,2,2,2,\n'
                'lower(:,1)=1,1,1,upper(:,1)=2,2,4,\n'
                'lower(:,2)=1,3,1,upper(:,2)=2,4,4,\n'
                'lower(:,3)=3,1,1,upper(:,3)=4,2,4,\n'
                'lower(:,4)=3,3,1,upper(:,4)=4,4,4,')
            extra += 'rebalance_interval_steps=1,rebalance_minimum_gain=0.,'
        if a.regrid:
            extra += ('regrid_count=3,regrid_times=5.e-9,1.e-8,1.5e-8,regrid_child_count=1,0,1,'
                      'regrid_lower(:,1,1)=3,1,1,regrid_upper(:,1,1)=4,4,4,'
                      'regrid_lower(:,1,3)=1,1,1,regrid_upper(:,1,3)=2,4,4,')
        if a.selected:
            extra += "composition_count=1,composition_species(1)='N2',composition_mole_fractions(1)=1.,"
        if a.mach:
            extra = ('inflow_mach(4)=.1,inflow_temperature(4)=400.,inflow_pressure(4)=120000.,'
                     'automatic_timestep=.true.,subcycling=.true.,' + extra)
        (directory / 'tree.nml').write_text(text.replace('&eb_tree', '&eb_tree\n' + extra, 1))
        if a.ports:
            spray = directory / 'spray.nml'
            values = spray.read_text().replace('nc=1,', "nc=2,vapor_names(2)='H2',")
            values = values.replace('&eb_spray', '&eb_spray\n port_count=2,event_ports=1,2,2,1,2,1,'
                'port_numbers(2)=3,port_centers(:,2)=.003,.001250001,.0015,'
                'port_masses(2)=2.e-10,port_speeds(2)=30.,port_diameters(2)=2.e-5,'
                'port_temperatures(2)=360.,port_fractions(:,2)=.25,.75,14*0.,', 1)
            spray.write_text(values)
        if a.profiles:
            spray = directory / 'spray.nml'
            values = spray.read_text().replace('&eb_spray', '&eb_spray\n'
                'port_profile_points=2,3,port_profile_times(:,1)=0.,4.e-8,62*0.,'
                'port_profile_rates(:,1)=.02,.02,62*0.,'
                'port_profile_times(:,2)=0.,2.e-8,4.e-8,61*0.,'
                'port_profile_rates(:,2)=0.,.08,0.,61*0.,'
                'injection_end_times=1.e-8,5.e-9,1.e-8,1.5e-8,2.5e-8,3.5e-8,', 1)
            values = values.replace('injection_scales=1.,0.,1.,.5,2.,1.,',
                                    'injection_scales=1.,1.,1.,.5,2.,1.,')
            spray.write_text(values)
        if a.pelemp_cap:
            assert not a.ports, 'PeleMP cap fixture uses one liquid component'
            spray = directory / 'spray.nml'
            spray.write_text(spray.read_text().replace('&eb_spray', '&eb_spray\n'
                "reference_temperatures(1)=300.,mixture_boiling_models(1)='pelemp-cap',"
                'mixture_critical_temperatures(1)=647.,mixture_normal_boiling_temperatures(1)=330.,', 1))
        if altered:
            spray = directory / 'spray.nml'
            if a.pelemp_cap:
                spray.write_text(spray.read_text().replace('temperatures(1)=330.', 'temperatures(1)=331.'))
            else:
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
    if a.regrid:
        for directory in results:
            log = (directory / 'run.log').read_text()
            assert all(f'Regrid stage: {n}' in log for n in (1, 2, 3)), log
    if a.rebalance:
        assert any('Rebalanced at step:' in (d / 'run.log').read_text() for d in results)
    gas = []
    particles = []
    for directory, ranks in zip(results, (2, 4)):
        gas.append(rows(directory / 'result.csv'))
        files = sorted(directory.glob('result.csv.particles.rank*.csv'))
        assert len(files) == ranks
        drops = sorted((row for file in files for row in rows(file)), key=lambda row: int(row['id']))
        assert [int(row['id']) for row in drops] == list(range(1, 11) if a.ports else range(1, 9) if a.scheduled else range(1, 3))
        if a.scheduled:
            assert 'Processed injection events: 5' in (directory / 'run.log').read_text()
            budget = rows(directory / 'result.csv.budget.csv')
            expected_mass = 1.e-9 if a.profiles else 7.5e-10 if a.ports else 4.5e-10
            assert abs(float(budget[0]['injected']) - expected_mass) < expected_mass * 1e-14
            for first, second in ([] if a.ports else zip(drops[:2], drops[2:4])):
                assert {k: v for k, v in first.items() if k != 'id'} == {
                    k: v for k, v in second.items() if k != 'id'}

        for row in drops:
            assert all(math.isfinite(float(v)) for v in row.values())
            assert float(row['mass']) > 0
            port = 2 if a.ports and int(row['id']) in (3, 4, 5, 8, 9, 10) else 1
            assert int(row['material_id']) == port
            initial_t = 360 if port == 2 else 350
            if a.scheduled and int(row['id']) >= (8 if a.ports else 7):
                assert float(row['temperature']) == initial_t and float(row['v']) == (-30 if port == 2 else -20), row
                assert float(row['y']) == .001250001, row
                if a.ports:
                    assert float(row['x']) == .003
                    assert abs(float(row['component_mass_1']) / float(row['mass']) - .25) < 1e-14
            else:
                assert float(row['y']) > .00125 and float(row['v']) > 0, row
                assert float(row['temperature']) != initial_t, row
            component_sum = sum(float(v) for k, v in row.items() if k.startswith('component_mass_'))
            assert abs(float(row['mass']) - component_sum) <= 1e-14 * float(row['mass'])
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
        if a.mach:
            mismatch = run(2, 'inlet_mismatch', f"restart_file='{archive}',inflow_mach(4)=.2,", success=False)
            assert 'Restart inlet mismatch' in (mismatch / 'run.log').read_text()
            assert not (mismatch / 'result.csv').exists()
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
