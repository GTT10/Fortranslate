#!/usr/bin/env python3
"""Focused input-to-output check for the fixed MPI spray executable."""
import argparse
import csv
import math
from pathlib import Path
import subprocess
import tempfile


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith('#')))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--mpiexec', type=Path, required=True)
    parser.add_argument('--ranks', type=int, required=True)
    parser.add_argument('--species-count', type=int, default=10)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='run-', dir=args.work_dir.resolve()))
    command = [str(args.mpiexec.resolve()), '-n', str(args.ranks),
               str(args.executable.resolve()), str(args.input.resolve())]
    result = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=90)
    (work / 'run.log').write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    history = rows(work / 'periodic_mixture.history.csv')
    assert len(history) == 3 and int(history[-1]['step']) == 2
    assert math.isclose(float(history[-1]['time']), 2e-8, rel_tol=1e-14)
    for row in history:
        assert all(math.isfinite(float(value)) for value in row.values())
        assert abs(float(row['relative_mass_error'])) < 2e-11
        assert abs(float(row['relative_energy_error'])) < 2e-11
        assert int(row['parcels']) == 4
    gas = [row for path in sorted(work.glob('*.gas.csv')) for row in rows(path)]
    assert len(gas) == 8 * 4 * 4
    cells = {(int(row['i']), int(row['j']), int(row['k'])) for row in gas}
    assert cells == {(i, j, k) for i in range(1, 9) for j in range(1, 5) for k in range(1, 5)}
    for row in gas:
        assert all(math.isfinite(float(value)) for value in row.values())
        rho = float(row['rho'])
        species = [float(value) for name, value in row.items() if name.startswith('rhoY_')]
        assert len(species) == args.species_count
        assert rho > 0 and float(row['temperature']) > 0 and min(species) >= 0
        assert math.isclose(sum(species), rho, rel_tol=2e-12)
    particles = [row for path in sorted(work.glob('*.particles.csv')) for row in rows(path)]
    assert sorted(int(row['id']) for row in particles) == [1, 2, 3, 4]
    for row in particles:
        assert all(math.isfinite(float(value)) for value in row.values())
        assert math.isclose(float(row['mass']), float(row['component_mass_1']) + float(row['component_mass_2']), rel_tol=2e-12)
        assert float(row['temperature']) > 330
    # A second invocation must fail before replacing the saved history.
    original = (work / 'periodic_mixture.history.csv').read_bytes()
    repeated = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=90)
    (work / 'repeat.log').write_text(repeated.stdout + repeated.stderr)
    assert repeated.returncode != 0
    assert (work / 'periodic_mixture.history.csv').read_bytes() == original
    source = args.input.read_text()

    def run_variant(name, rank_count, additions, stop=None, expect_success=True):
        text = source.replace("output_prefix='periodic_mixture'", f"output_prefix='{name}'")
        if stop is not None:
            text = text.replace('final_time=2e-8', f'final_time={stop}')
        text = text.replace('\n/', '\n' + additions + '\n/')
        path = work / f'{name}.nml'
        path.write_text(text)
        invocation = [str(args.mpiexec.resolve()), '-n', str(rank_count),
                      str(args.executable.resolve()), str(path)]
        output = subprocess.run(invocation, cwd=work, capture_output=True, text=True, timeout=90)
        (work / f'{name}.log').write_text(output.stdout + output.stderr)
        assert (output.returncode == 0) == expect_success, output.stdout + output.stderr

    def snapshots(prefix, suffix, keys):
        records = [row for path in sorted(work.glob(f'{prefix}.rank*.{suffix}.csv')) for row in rows(path)]
        return {tuple(int(row[key]) for key in keys): row for row in records}

    run_variant('split', args.ranks, " checkpoint_path='split.chk'", stop='1e-8')
    assert (work / 'split.chk').is_file()
    reference_gas = snapshots('periodic_mixture', 'gas', ['i', 'j', 'k'])
    reference_particles = snapshots('periodic_mixture', 'particles', ['id'])
    for name, count in [('resumed', args.ranks), ('redistributed', 6 - args.ranks)]:
        run_variant(name, count, " restart_path='split.chk'")
        for reference, actual in [
            (reference_gas, snapshots(name, 'gas', ['i', 'j', 'k'])),
            (reference_particles, snapshots(name, 'particles', ['id'])),
        ]:
            assert actual.keys() == reference.keys()
            for key, expected in reference.items():
                for column, value in expected.items():
                    if name == 'resumed':
                        assert actual[key][column] == value, (name, key, column)
                    else:
                        assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (name, key, column)
        resumed_history = rows(work / f'{name}.history.csv')
        assert int(resumed_history[0]['step']) == 1 and int(resumed_history[-1]['step']) == 2
        assert int(resumed_history[-1]['parcels']) == 4
        assert abs(float(resumed_history[-1]['relative_energy_error'])) < 2e-11
    run_variant('incompatible', args.ranks, " restart_path='split.chk'\n liquid_cp=4000,2000", expect_success=False)
    assert not (work / 'incompatible.history.csv').exists()
    schedule = ' injection_events=3\n injection_interval=7.5e-9'
    run_variant('scheduled', args.ranks, schedule)
    schedule_history = rows(work / 'scheduled.history.csv')
    assert [int(row['events']) for row in schedule_history] == [1, 2, 3, 3]
    assert [int(row['parcels']) for row in schedule_history] == [4, 8, 12, 12]
    assert math.isclose(float(schedule_history[-1]['injected_mass']), 3e-10, rel_tol=1e-14)
    assert int(schedule_history[-1]['last_id']) == 12
    assert all(abs(float(row['relative_mass_error'])) < 2e-11 and
               abs(float(row['relative_energy_error'])) < 2e-11 for row in schedule_history)
    run_variant('scheduled_split', args.ranks, schedule + "\n checkpoint_path='scheduled.chk'", stop='7.5e-9')
    split_history = rows(work / 'scheduled_split.history.csv')
    assert int(split_history[-1]['events']) == 2
    run_variant('scheduled_resume', 6 - args.ranks, schedule + "\n restart_path='scheduled.chk'")
    for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
        reference = snapshots('scheduled', suffix, keys)
        actual = snapshots('scheduled_resume', suffix, keys)
        assert actual.keys() == reference.keys()
        for key, expected in reference.items():
            for column, value in expected.items():
                assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (key, column)
    resumed_history = rows(work / 'scheduled_resume.history.csv')
    assert int(resumed_history[0]['events']) == 2 and int(resumed_history[-1]['events']) == 3
    assert int(resumed_history[-1]['last_id']) == 12 and int(resumed_history[-1]['parcels']) == 12
    profile = (schedule + '\n injection_profile_points=4' +
               '\n injection_profile_times=0,7.5e-9,1.5e-8,2.25e-8' +
               '\n injection_profile_rates=0,0,.02,0')
    run_variant('profiled', args.ranks, profile)
    profile_history = rows(work / 'profiled.history.csv')
    assert [int(row['parcels']) for row in profile_history] == [0, 4, 8, 8]
    assert int(profile_history[-1]['events']) == 3 and int(profile_history[-1]['last_id']) == 8
    assert math.isclose(float(profile_history[-1]['injected_mass']), 1.5e-10, rel_tol=1e-14)
    assert all(abs(float(row['relative_mass_error'])) < 2e-11 and
               abs(float(row['relative_energy_error'])) < 2e-11 for row in profile_history)
    run_variant('profiled_split', args.ranks, profile + "\n checkpoint_path='profiled.chk'", stop='7.5e-9')
    run_variant('profiled_resume', 6 - args.ranks, profile + "\n restart_path='profiled.chk'")
    for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
        reference = snapshots('profiled', suffix, keys)
        actual = snapshots('profiled_resume', suffix, keys)
        assert actual.keys() == reference.keys()
        for key, expected in reference.items():
            for column, value in expected.items():
                assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (key, column)
    run_variant('profiled_mismatch', args.ranks,
                profile.replace('rates=0,0,.02,0', 'rates=0,0,.03,0') + "\n restart_path='profiled.chk'", expect_success=False)
    assert not (work / 'profiled_mismatch.history.csv').exists()
    ports = """ port_count=2
 port_events=2,2
 port_start=0,7.5e-9
 port_interval=7.5e-9,7.5e-9
 port_number=2,3
 port_mass=1e-10,2e-10
 port_center(:,1)=.002,.002,.002
 port_center(:,2)=.006,.002,.002
 port_axis(:,2)=-1,0,0
 port_fractions(:,1)=1,0
 port_fractions(:,2)=0,1
 port_diameter=2e-5,3e-5
 port_temperature=330,340"""
    run_variant('ports', args.ranks, ports)
    port_history = rows(work / 'ports.history.csv')
    assert [int(row['events']) for row in port_history] == [1, 3, 4, 4]
    assert [int(row['parcels']) for row in port_history] == [2, 7, 10, 10]
    assert math.isclose(float(port_history[-1]['injected_mass']), 6e-10, rel_tol=1e-14)
    assert all(abs(float(row['relative_mass_error'])) < 2e-11 and
               abs(float(row['relative_energy_error'])) < 2e-11 for row in port_history)
    port_particles = snapshots('ports', 'particles', ['id'])
    assert set(port_particles) == {(i,) for i in range(1, 11)}
    for (ident,), row in port_particles.items():
        # At the tie, port 1 reserves IDs 3/4 before port 2 reserves 5/6/7.
        assert float(row['component_mass_2' if ident <= 4 else 'component_mass_1']) == 0
        assert float(row['x']) < .004 if ident <= 4 else float(row['x']) > .004
        assert float(row['u']) > 0 if ident <= 4 else float(row['u']) < 0
    run_variant('ports_split', args.ranks, ports + "\n checkpoint_path='ports.chk'", stop='7.5e-9')
    assert int(rows(work / 'ports_split.history.csv')[-1]['events']) == 3
    run_variant('ports_resume', 6 - args.ranks, ports + "\n restart_path='ports.chk'")
    for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
        reference = snapshots('ports', suffix, keys)
        actual = snapshots('ports_resume', suffix, keys)
        assert actual.keys() == reference.keys()
        for key, expected in reference.items():
            for column, value in expected.items():
                assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (key, column)
    run_variant('ports_mismatch', args.ranks,
                ports.replace('port_mass=1e-10,2e-10', 'port_mass=1e-10,3e-10') +
                "\n restart_path='ports.chk'", expect_success=False)
    assert not (work / 'ports_mismatch.history.csv').exists()
    run_variant('ports_profile', args.ranks, ports +
                '\n injection_profile_points=4' +
                '\n injection_profile_times=0,7.5e-9,1.5e-8,2.25e-8' +
                '\n injection_profile_rates=0,0,.02,0')
    profiled_ports = rows(work / 'ports_profile.history.csv')
    assert [int(row['parcels']) for row in profiled_ports] == [0, 5, 8, 8]
    assert math.isclose(float(profiled_ports[-1]['injected_mass']), 2.25e-10, rel_tol=1e-14)
    assert int(profiled_ports[-1]['last_id']) == 8
    print(f'PASS MPI spray application, profile, multiport schedule and restart: {args.ranks} ranks; artifacts: {work}')



if __name__ == '__main__':
    main()
