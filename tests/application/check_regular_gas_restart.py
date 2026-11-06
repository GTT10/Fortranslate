"""Regular gas state/IR/heat and cumulative-ledger continuation."""
import argparse
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('exe', 'source', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--selected', action='store_true')
    parser.add_argument('--rotation', action='store_true')
    parser.add_argument('--dual-energy-policy', action='store_true')
    parser.add_argument('--timestep-history', action='store_true')
    parser.add_argument('--prescribed-timestep', choices=('fixed', 'combined'))
    parser.add_argument('--les', action='store_true')
    parser.add_argument('--reactor', choices=('native', 'cvode'), default='native')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='regular_restart_', dir=args.work.resolve()))
    text = (args.source / 'cases/reactive_boundaries_3d/sdc_reacting.nml').read_text()
    text = text.replace('final_time=1e-8', 'final_time=3e-10')
    if args.selected:
        text = text.replace("thermo_model='full_h2o2'", "thermo_model='selected'")
        text = text.replace("x_h2=.29570,x_h=1e-5,x_o=1e-5,x_o2=.14784,x_oh=1e-5,x_n2=.55643",
                            "composition_count=5,composition_species='H2','O2','N2','H','OH',"
                            "composition_mole_fractions=.29570,.14784,.55644,1e-5,1e-5")
    # CVODE is supported by SDC and MOL, while split uses the bundle integrator.
    policy = 'initial_timestep_factor=.5,maximum_timestep_growth=1.1,' if args.timestep_history else ''
    maximum = '1e-10'
    first_interval = 5e-11 if args.timestep_history else 1e-10
    last_interval = 6.795e-11
    if args.prescribed_timestep:
        policy = 'initial_timestep_factor=.5,maximum_timestep_growth=1,fixed_dt=6e-11,'
        maximum = '6e-11'
        first_interval, last_interval = 3e-11, 3e-11
        if args.prescribed_timestep == 'combined':
            policy += 'initial_dt=8e-11,'
            first_interval, last_interval = 8e-11, 4e-11
    text += f'\n&reactive_timestep {policy}maximum_dt={maximum} /\n'

    def run(label, content, restart='', checkpoint='', stop=0, reject=False, rejection='Regular gas restart rejected'):
        directory = work / label
        directory.mkdir()
        inp = directory / 'input.nml'
        inp.write_text(content + f"\n&reactive_restart restart_file='{restart}',"
                       f"checkpoint_file='{checkpoint}',stop_after_steps={stop} /\n")
        result = subprocess.run([str(args.exe.resolve()), str(inp)], cwd=directory,
                                capture_output=True, text=True, timeout=120)
        log = result.stdout + result.stderr
        (directory / 'run.log').write_text(log)
        if reject:
            assert result.returncode != 0 and rejection in log, log
            assert not (directory / 'sdc_reacting.csv').exists()
        else:
            assert result.returncode == 0, log
            residual = float(re.search(r'Maximum normalized Euler conserved change:\s+(\S+)', log).group(1))
            assert residual < 3e-12, residual
        return directory

    modes = ('sdc', 'mol') if args.reactor == 'cvode' else ('split', 'sdc', 'mol')
    rotation = "\n&reactive_rotation enabled=.true.,axis=2,angular_speed=300,axis_point=.02,-.03,.01 /\n"
    if args.rotation:
        modes = ('mol',)
        text += rotation
    if args.dual_energy_policy:
        modes = ('mol',)
        text = text.replace("time_integrator='sdc'", "time_integrator='mol',allow_small_energy=.true.,"
                            "allow_negative_energy=.false.,dual_energy_update_E_from_e=.true.,dual_energy_eta2=1")
    for mode in modes:
        content = text
        if args.reactor == 'cvode':
            content = content.replace("sdc_iterations=2", "sdc_iterations=2,sdc_reactor='cvode'")
        if mode == 'split':
            content = content.replace("time_integrator='sdc'", "time_integrator='split'")
        if mode == 'mol':
            content = content.replace("time_integrator='sdc'", "time_integrator='mol'")
            content = content.replace("reconstruction='characteristic_ppm'", "reconstruction='characteristic_plm'")
        if args.prescribed_timestep:
            run(mode + '_fixed_above_max', content.replace('fixed_dt=6e-11', 'fixed_dt=7e-11'),
                reject=True, rejection='Fixed timestep exceeds maximum timestep')
        if args.prescribed_timestep == 'fixed':
            disabled = work / f'{mode}_negative_disabled.chk'
            inactive = content.replace('fixed_dt=6e-11', 'fixed_dt=6e-11,initial_dt=-1')
            run(mode + '_nonpositive_initial', inactive, checkpoint=disabled, stop=1)
            records = disabled.read_text().splitlines()
            marker = records.index('PELEF_REGULAR_GAS_TIMESTEP_V1')
            assert float(records[marker + 1]) == first_interval, 'Nonpositive initial_dt must be inactive'
        full = work / f'{mode}_full.chk'
        stop = work / f'{mode}_stop.chk'
        resumed = work / f'{mode}_resumed.chk'
        a = run(mode + '_full', content, checkpoint=full)
        run(mode + '_stop', content, checkpoint=stop, stop=1)
        if args.timestep_history or args.prescribed_timestep:
            lines = stop.read_text().splitlines()
            header = lines.index('PELEF_REGULAR_GAS_TIMESTEP_V1')
            assert float(lines[header + 1]) == first_interval, 'First interval policy not applied'
        b = run(mode + '_resume', content, restart=stop, checkpoint=resumed)
        assert (a / 'sdc_reacting.csv').read_bytes() == (b / 'sdc_reacting.csv').read_bytes()
        assert full.read_bytes() == resumed.read_bytes(), 'State, IR, heat, progress and cumulative ledgers differ'
        if args.timestep_history or args.prescribed_timestep:
            lines = full.read_text().splitlines()
            header = lines.index('PELEF_REGULAR_GAS_TIMESTEP_V1')
            assert abs(float(lines[header + 1]) - last_interval) < 1e-24, 'Growth history lost'
            old_growth = 'maximum_timestep_growth=1,' if args.prescribed_timestep else 'maximum_timestep_growth=1.1'
            new_growth = 'maximum_timestep_growth=1.2,' if args.prescribed_timestep else 'maximum_timestep_growth=1.2'
            run(mode + '_changed_growth', content.replace(old_growth, new_growth), restart=stop, reject=True)
        if args.prescribed_timestep:
            identity_end = lines.index('END_REGULAR_GAS_IDENTITY')
            progress = lines[identity_end + 2].split()
            assert int(progress[1]) == (6 if args.prescribed_timestep == 'fixed' else 5), progress
            run(mode + '_changed_fixed', content.replace('fixed_dt=6e-11', 'fixed_dt=5e-11'),
                restart=stop, reject=True)
            if args.prescribed_timestep == 'combined':
                run(mode + '_changed_initial', content.replace('initial_dt=8e-11', 'initial_dt=7e-11'),
                    restart=stop, reject=True)
        repacked = work / f'{mode}_repacked.chk'
        run(mode + '_repack', content.replace('final_time=3e-10', f'final_time={first_interval:.17e}'),
            restart=stop, checkpoint=repacked)
        assert stop.read_bytes() == repacked.read_bytes()
        if args.dual_energy_policy:
            for name, before, after in (
                ('small', 'allow_small_energy=.true.', 'allow_small_energy=.false.'),
                ('negative', 'allow_negative_energy=.false.', 'allow_negative_energy=.true.'),
                ('update', 'dual_energy_update_E_from_e=.true.', 'dual_energy_update_E_from_e=.false.'),
                ('eta', 'dual_energy_eta2=1', 'dual_energy_eta2=.5'),
            ):
                run(mode + '_changed_energy_' + name, content.replace(before, after), restart=stop, reject=True)
            for name, value in (('negative', '-1'), ('nonfinite', 'NaN')):
                run(mode + '_invalid_eta_' + name, content.replace('dual_energy_eta2=1', 'dual_energy_eta2=' + value),
                    reject=True, rejection='Dual energy eta2 must be finite and nonnegative')
        run(mode + '_changed_force', content.replace('external_forcing=20000', 'external_forcing=20001'),
            restart=stop, reject=True)
        damaged = work / f'{mode}_truncated.chk'
        records = stop.read_text().splitlines()
        damaged.write_text('\n'.join(records[:-2]) + '\n')
        run(mode + '_truncated', content, restart=damaged, reject=True)
        bad_ledger = work / f'{mode}_bad_ledger.chk'
        ledger_records = stop.read_text().splitlines()
        index = ledger_records.index('END_REGULAR_GAS_IDENTITY') + 3
        values = ledger_records[index].split()
        values[0] = str(2 * float(values[0]))
        ledger_records[index] = ' '.join(values)
        bad_ledger.write_text('\n'.join(ledger_records) + '\n')
        run(mode + '_finite_bad_ledger', content, restart=bad_ledger, reject=True)
        corrupted = work / f'{mode}_nan.chk'
        values = records[-2].split()
        values[-1] = 'NaN'
        records[-2] = ' '.join(values)
        corrupted.write_text('\n'.join(records) + '\n')
        run(mode + '_nonfinite', content, restart=corrupted, reject=True)
        if args.rotation:
            for name, changed in (
                ('axis', content.replace('axis=2', 'axis=1')),
                ('speed', content.replace('angular_speed=300', 'angular_speed=301')),
                ('point', content.replace('axis_point=.02', 'axis_point=.021')),
                ('disabled', content.replace(rotation, '')),
            ):
                run(mode + '_rotation_' + name, changed, restart=stop, reject=True)
            zero = run(mode + '_zero_rotation', content.replace('angular_speed=300', 'angular_speed=0'))
            plain = run(mode + '_ordinary', content.replace(rotation, ''))
            assert (zero / 'sdc_reacting.csv').read_bytes() == (plain / 'sdc_reacting.csv').read_bytes()
            assert (a / 'sdc_reacting.csv').read_bytes() != (plain / 'sdc_reacting.csv').read_bytes()
            if args.les:
                les_content = content + "\n&reactive_les enabled=.true.,model='dynamic'," \
                    "test_model='gaussian',test_ratio=4,grid_ratio=2 /\n"
                les_full = work / 'mol_les_full.chk'
                les_stop = work / 'mol_les_stop.chk'
                les_resume = work / 'mol_les_resume.chk'
                la = run('mol_les_full', les_content, checkpoint=les_full)
                run('mol_les_stop', les_content, checkpoint=les_stop, stop=1)
                lb = run('mol_les_resume', les_content, restart=les_stop, checkpoint=les_resume)
                assert (la / 'sdc_reacting.csv').read_bytes() == (lb / 'sdc_reacting.csv').read_bytes()
                assert les_full.read_bytes() == les_resume.read_bytes()


    print(f'PASS regular {modes} selected={args.selected} reactor={args.reactor} rotation={args.rotation} les={args.les} '
          f'CSV and checkpoint exact restart, repack and rejection: {work}')


if __name__ == '__main__':
    main()
