"""Check serial/MPI particle timestep limits and accepted history on restart."""
import argparse
import csv
from pathlib import Path
import subprocess
import sys
import tempfile


def rows(path):
    return list(csv.DictReader(line for line in path.read_text().splitlines() if not line.startswith('#')))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('serial', 'mpi', 'mpiexec', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--policy', choices=('auto', 'fixed', 'combined', 'events'), default='auto')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='particle_timestep_', dir=args.work.resolve()))
    for parallel in (False, True):
        for mode in ('sdc', 'mol'):
            name = ('mpi_' if parallel else 'serial_') + mode
            # Uniform stationary particles remove trajectory differences from the
            # timestep control; their nonzero mass must survive rank-changing I/O.
            if parallel:
                base = """&mpi_spray cells=8,4,4,speed=0,spacing=.001,.001,.001,number=2,total_mass=1e-10,
final_time=1.655e-8,dt_max=1e-8,initial_timestep_factor=.5,maximum_timestep_growth=1.1,
chemistry_enabled=.false.,solver='pelec',time_integrator='MODE',sdc_iterations=2,
hydro_reconstruction='RECON',OUTPUT,CONTROL /"""
            else:
                base = """&spray cells=4,4,4,injection_speed=0,parcel_count=2,chemistry_enabled=.false.,evaporation_enabled=.false.,
final_time=1.655e-8,maximum_dt=1e-8,initial_timestep_factor=.5,maximum_timestep_growth=1.1,
time_integrator='MODE',sdc_iterations=2,hydro_reconstruction='RECON',riemann_solver='pelec',
OUTPUT,CONTROL /"""
            first, last, expected_steps = 5e-9, 6.05e-9, 3
            extra = ''
            if args.policy in ('fixed', 'combined'):
                base = base.replace('dt_max=1e-8', 'dt_max=6e-9').replace('maximum_dt=1e-8', 'maximum_dt=6e-9')
                extra = 'fixed_dt=6e-9,'
                first, last, expected_steps = 3e-9, 1.55e-9, 4
                if args.policy == 'combined':
                    extra += 'initial_dt=8e-9,'
                    first, last, expected_steps = 8e-9, 2.55e-9, 3
            elif args.policy == 'events':
                extra = ('injection_events=3,injection_interval=5e-9,' if parallel else
                         'mass_flow_rate=1e-3,injection_end=1.5e-8,injection_interval=5e-9,')
                last, expected_steps = 1.05e-9, 4
            base = base.replace('OUTPUT,CONTROL', extra + 'OUTPUT,CONTROL')
            base = base.replace('MODE', mode).replace('RECON', 'characteristic_ppm' if mode == 'sdc' else 'characteristic_plm')

            def run(label, controls='', ranks=1, reject=False, changed=''):
                folder = work / (name + '_' + label)
                folder.mkdir()
                content = base.replace('OUTPUT', "output_prefix='out'").replace('CONTROL', controls)
                if changed == 'growth':
                    content = content.replace('maximum_timestep_growth=1.1', 'maximum_timestep_growth=1.2')
                elif changed == 'fixed':
                    content = content.replace('fixed_dt=6e-9', 'fixed_dt=5e-9')
                elif changed == 'initial':
                    content = content.replace('initial_dt=8e-9', 'initial_dt=7e-9')
                elif changed == 'invalid_fixed':
                    content = content.replace('fixed_dt=6e-9', 'fixed_dt=7e-9')
                (folder / 'input.nml').write_text(content + '\n')
                command = ([str(args.mpiexec.resolve()), '-n', str(ranks), str(args.mpi.resolve())]
                           if parallel else [str(args.serial.resolve())]) + ['input.nml']
                result = subprocess.run(command, cwd=folder, text=True, capture_output=True, timeout=120)
                log = result.stdout + result.stderr
                (folder / 'run.log').write_text(log)
                if reject:
                    assert result.returncode != 0, log
                else:
                    assert result.returncode == 0, log
                return folder

            checkpoint_key = 'checkpoint_path' if parallel else 'checkpoint_file'
            restart_key = 'restart_path' if parallel else 'restart_file'
            stop = 'stop_after_steps=1' if not parallel else f'final_time={first:.17e}'
            a = run('full', f"{checkpoint_key}='full.chk'", ranks=2 if parallel else 1)
            stopped = run('stop', f"{checkpoint_key}='stop.chk',{stop}", ranks=2 if parallel else 1)
            restart = f"{restart_key}='{stopped / 'stop.chk'}'"
            b = run('resume', restart + f",{checkpoint_key}='resumed.chk'", ranks=1)
            full = (a / 'full.chk').read_bytes()
            resumed = (b / 'resumed.chk').read_bytes()
            if parallel:
                same = run('resume_same_ranks', restart + f",{checkpoint_key}='resumed.chk'", ranks=2)
                assert full == (same / 'resumed.chk').read_bytes(), name + ' same ranks'
                # Rank changes alter reduction and halo operation ordering.
                # Metadata, particles and source histories remain exact; compare
                # gas and boundary roundoff with a machine-precision bound.
                original = full.decode().splitlines()
                restored = resumed.decode().splitlines()
                assert len(original) == len(restored)
                gas_end = 11 + 8 * 4 * 4
                assert original[:10] == restored[:10], name + ' metadata/history'
                if args.policy != 'events':
                    assert original[gas_end:] == restored[gas_end:], name + ' particles/source history'
                else:
                    # Retain changed-rank particle differences for diagnosis;
                    # the strict equality assertion follows the gas comparison.
                    differences = [(i + 1, x, y) for i, (x, y) in enumerate(zip(original, restored))
                                   if i >= gas_end and x != y]
                    (b / 'changed_rank_particle_differences.txt').write_text(
                        '\n'.join(f'{i}: {x}\n{y}' for i, x, y in differences) + '\n')
                    print(name, 'changed-rank particle/source differing records', len(differences))
                max_scaled = 0.0
                for xrow, yrow in zip(original[11:gas_end], restored[11:gas_end]):
                    xs, ys = list(map(float, xrow.split())), list(map(float, yrow.split()))
                    assert len(xs) == len(ys) == 16
                    for x, y in zip(xs, ys):
                        bound = 128 * sys.float_info.epsilon * max(1e-12, abs(x), abs(y))
                        max_scaled = max(max_scaled, abs(x-y)/bound)
                        assert abs(x-y) <= bound, (name, x, y, bound)
                print(name, 'changed-rank gas maximum roundoff bound fraction', max_scaled)
                scales = [abs(float(x)) for x in original[8].split()]
                for x, y, scale in zip(original[10].split(), restored[10].split(), scales):
                    assert abs(float(x) - float(y)) <= 128 * sys.float_info.epsilon * max(1e-12, scale), (x, y)
                if args.policy == 'events':
                    assert not differences, name + ' changed-rank particle/source mismatch'
            else:
                assert full == resumed, name
            history = rows(a / 'out.history.csv')
            assert int(history[-1]['step' if parallel else 'steps']) == expected_steps, history[-1]
            history = rows(stopped / 'out.history.csv')
            assert int(history[-1]['step' if parallel else 'steps']) == 1 and float(history[-1]['time']) == first, history[-1]
            assert (stopped / 'stop.chk').read_text().startswith('PELEF_SPRAY_TIMESTEP_HISTORY_V1\n')
            assert abs(float((stopped / 'stop.chk').read_text().splitlines()[7]) - first) < 1e-23
            assert abs(float(full.decode().splitlines()[7]) - last) < 1e-23
            if args.policy in ('fixed', 'combined'):
                run('changed_fixed', restart, ranks=2 if parallel else 1, reject=True, changed='fixed')
                run('invalid_fixed', ranks=2 if parallel else 1, reject=True, changed='invalid_fixed')
            if args.policy == 'combined':
                run('changed_initial', restart, ranks=2 if parallel else 1, reject=True, changed='initial')
            run('changed_growth', restart, ranks=2 if parallel else 1, reject=True, changed='growth')
    print('PASS', args.policy, 'serial/MPI SDC/MOL particle interval history and 2-to-1-rank restart:', work)


if __name__ == '__main__':
    main()
