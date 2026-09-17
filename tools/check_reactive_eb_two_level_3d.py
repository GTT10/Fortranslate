#!/usr/bin/env python3
"""Check the fixed elementary two-level EB application's public artifacts."""
import argparse
import csv
import math
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    exe = args.exe.resolve()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    source = args.case.read_text()

    def run(name, content, success):
        path = work / name
        path.write_text(content)
        result = subprocess.run([str(exe), name], cwd=work, text=True, capture_output=True, timeout=90)
        (work / (name + '.log')).write_text(result.stdout + result.stderr)
        assert (result.returncode == 0) == success, result.stdout + result.stderr
        return result.stdout

    log = run('reference.nml', source, True)
    assert 'Completed steps: 3' in log, log
    output = work / 'eb_two_level.csv'
    budget_path = work / 'eb_two_level.csv.budget.csv'
    with output.open() as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        fields = reader.fieldnames[9:]
    assert len(rows) == 104, len(rows)
    totals = [0.0] * len(fields)
    initial_mass = 0.0
    levels = set()
    keys = set()
    for row in rows:
        assert all(math.isfinite(float(value)) for value in row.values())
        level = int(row['level'])
        ijk = tuple(int(row[c]) for c in ('i', 'j', 'k'))
        assert (level, ijk) not in keys
        keys.add((level, ijk))
        levels.add(level)
        assert level in (0, 1) and all(1 <= i <= 4 for i in ijk)
        if level == 0:
            assert not all(2 <= i <= 3 for i in ijk), 'refined parent duplicated in composite output'
        kappa = float(row['kappa'])
        assert 0 < kappa <= 1 and float(row['temperature']) > 0 and float(row['rho']) > 0
        volume = kappa / (4 if level == 0 else 8) ** 3
        initial_mass += volume * (0.93 if kappa < 1 - 1e-12 else 0.31)
        species_sum = sum(float(value) for key, value in row.items() if key.startswith('rhoY_'))
        assert math.isclose(species_sum, float(row['rho']), rel_tol=1e-10, abs_tol=1e-12)
        for c, name in enumerate(fields):
            totals[c] += volume * float(row[name])
    assert levels == {0, 1}
    with budget_path.open() as stream:
        budgets = list(csv.DictReader(stream))
    assert len(budgets) == len(fields)
    assert math.isclose(float(budgets[0]['initial']), initial_mass, rel_tol=2e-12, abs_tol=1e-13)
    for c, row in enumerate(budgets):
        values = {key: float(value) for key, value in row.items()}
        assert all(math.isfinite(value) for value in values.values())
        assert int(values['component']) == c + 1
        assert math.isclose(totals[c], values['final'], rel_tol=3e-12, abs_tol=3e-12)
        residual = values['final'] - values['initial'] - values['boundary'] - values['wall'] - values['chemistry']
        assert abs(residual) <= 5e-10 * max(1.0, abs(values['initial']), abs(values['final']))
        assert math.isclose(residual, values['residual'], rel_tol=1e-12, abs_tol=1e-15)
    saved = output.read_bytes(), budget_path.read_bytes()
    run('bad_ratio.nml', source.replace('refinement_ratio=2', 'refinement_ratio=3'), False)
    assert saved == (output.read_bytes(), budget_path.read_bytes())
    alias = source.replace("output_file='eb_two_level.csv'", "output_file='alias.nml'")
    run('alias.nml', alias, False)
    assert (work / 'alias.nml').read_text() == alias
    def options(**values):
        extra = ''.join(f"  {key}={value},\n" for key, value in values.items())
        return source.replace("\n/", "\n" + extra + "/")

    stopped = options(checkpoint_file="'two_level.chk'", checkpoint_interval_steps=1, stop_after_checkpoint='.true.')
    stopped = stopped.replace('final_time=3.0e-8', 'final_time=1.0e-8').replace('maximum_steps=10', 'maximum_steps=1')
    log = run('stopped.nml', stopped, True)
    assert 'Completed steps: 1' in log
    checkpoint = work / 'two_level.chk'
    original_checkpoint = checkpoint.read_bytes()
    log = run('resumed.nml', options(restart_file="'two_level.chk'"), True)
    assert 'Completed steps: 3' in log
    assert saved == (output.read_bytes(), budget_path.read_bytes()), 'restart differs from uninterrupted output'
    assert checkpoint.read_bytes() == original_checkpoint
    # Reject corruption in the second patch, trailing data and incompatible
    # numerical controls before writing either application output.
    records = original_checkpoint.decode()
    second = records.index('PELEF_REACTIVE_EB_3D_CHECKPOINT', records.index('PELEF_REACTIVE_EB_3D_CHECKPOINT') + 1)
    damaged_budget = records.splitlines()
    values = damaged_budget[4].split()
    values[2 * len(fields)] = '1.0e10'
    damaged_budget[4] = ' '.join(values)
    for name, content in [('truncated', records[:second + 100]),
                          ('budget', '\n'.join(damaged_budget) + '\n'),
                          ('trailing', records + 'garbage\n'),
                          ('terminator', records.replace('END_TWO_LEVEL_CHECKPOINT', 'BAD_END'))]:
        (work / (name + '.chk')).write_text(content)
        run(name + '.nml', options(restart_file=f"'{name}.chk'"), False)
        assert saved == (output.read_bytes(), budget_path.read_bytes())
    run('changed_dt.nml', options(restart_file="'two_level.chk'").replace('time_step=1.0e-8', 'time_step=2.0e-8'), False)
    run('changed_physics.nml', options(restart_file="'two_level.chk'").replace('transport_enabled=.true.',
                                                                            'transport_enabled=.false.'), False)
    assert saved == (output.read_bytes(), budget_path.read_bytes())
    alias = options(checkpoint_file="'checkpoint_alias.nml'", checkpoint_interval_steps=1)
    run('checkpoint_alias.nml', alias, False)
    assert (work / 'checkpoint_alias.nml').read_text() == alias
    retry_case = source.replace('chemistry_enabled=.true.', 'chemistry_enabled=.false.').replace(
        'transport_enabled=.true.', 'transport_enabled=.false.').replace('time_step=1.0e-8', 'time_step=1.0e-3').replace(
        'final_time=3.0e-8', 'final_time=1.0e-3').replace('maximum_steps=10', 'maximum_steps=1000')
    log = run('retry.nml', retry_case, True)
    rejected = int(log.split('Rejected trials this invocation: ')[1].splitlines()[0])
    assert rejected > 0, 'case must exercise rejected updates'
    retry_saved = output.read_bytes(), budget_path.read_bytes()
    exhausted = retry_case.replace('\n/', '\n maximum_step_retries=0,\n/')
    run('exhausted.nml', exhausted, False)
    assert retry_saved == (output.read_bytes(), budget_path.read_bytes())
    stopped_retry = retry_case.replace('\n/', "\n checkpoint_file='retry.chk', checkpoint_interval_steps=1, "
                                      "stop_after_checkpoint=.true.,\n/")
    log = run('retry_stop.nml', stopped_retry, True)
    assert 'Completed steps: 1' in log
    first_saved = output.read_bytes(), budget_path.read_bytes()
    accepted_dt = float((work / 'retry.chk').read_text().splitlines()[3].split()[1])
    assert 0 < accepted_dt < 1e-3
    # Replaying only the accepted interval with retries disabled proves the
    # rejected trials did not alter state or cumulative inventory diagnostics.
    direct = exhausted.replace('time_step=1.0e-3', f'time_step={accepted_dt:.17e}').replace(
        'final_time=1.0e-3', f'final_time={accepted_dt:.17e}')
    run('retry_direct.nml', direct, True)
    assert first_saved == (output.read_bytes(), budget_path.read_bytes()), 'rejected trial leaked state'
    resumed_retry = retry_case.replace('\n/', "\n restart_file='retry.chk',\n/")
    run('retry_resume.nml', resumed_retry, True)
    assert retry_saved == (output.read_bytes(), budget_path.read_bytes()), 'retry continuation mismatch'
    invalid_retry = retry_case.replace('\n/', '\n maximum_step_retries=-1,\n/')
    run('invalid_retry.nml', invalid_retry, False)
    assert retry_saved == (output.read_bytes(), budget_path.read_bytes())
    # Automatic mode also visits transport selection in the reference coupled
    # run; its small user cap should retain the original accepted arithmetic.
    run('auto_coupled.nml', options(automatic_timestep='.true.'), True)
    assert saved == (output.read_bytes(), budget_path.read_bytes())
    auto_case = retry_case.replace('final_time=1.0e-3', 'final_time=1.0e-4').replace(
        '\n/', '\n automatic_timestep=.true.,\n/')
    log = run('auto.nml', auto_case, True)
    assert 'Rejected trials this invocation: 0' in log, log
    auto_saved = output.read_bytes(), budget_path.read_bytes()
    auto_stop = auto_case.replace('\n/', "\n checkpoint_file='auto.chk', checkpoint_interval_steps=1, "
                                 "stop_after_checkpoint=.true.,\n/")
    run('auto_stop.nml', auto_stop, True)
    first_auto = output.read_bytes(), budget_path.read_bytes()
    auto_records = (work / 'auto.chk').read_text().splitlines()
    assert auto_records[0].strip() == 'PELEF_EB_TWO_LEVEL_CHECKPOINT_V2'
    auto_dt = float(auto_records[3].split()[1])
    assert 0 < auto_dt < 1e-4
    direct_auto = auto_case.replace('automatic_timestep=.true.', 'automatic_timestep=.false.').replace(
        'time_step=1.0e-3', f'time_step={auto_dt:.17e}').replace('final_time=1.0e-4', f'final_time={auto_dt:.17e}')
    run('auto_direct.nml', direct_auto, True)
    assert first_auto == (output.read_bytes(), budget_path.read_bytes())
    auto_resume = auto_case.replace('\n/', "\n restart_file='auto.chk',\n/")
    run('auto_resume.nml', auto_resume, True)
    assert auto_saved == (output.read_bytes(), budget_path.read_bytes())
    run('auto_wrong_mode.nml', auto_resume.replace('automatic_timestep=.true.', 'automatic_timestep=.false.'), False)
    run('auto_wrong_cfl.nml', auto_resume.replace('\n/', '\n hydro_cfl=0.25,\n/'), False)
    assert auto_saved == (output.read_bytes(), budget_path.read_bytes())
    # Micrometre geometry makes diffusion dominate; halve only transport CFL
    # and require the saved first interval to halve as well.
    diffusion_case = source.replace('plane_position=0.31', 'plane_position=0.31e-6').replace(
        'chemistry_enabled=.true.', 'chemistry_enabled=.false.').replace('time_step=1.0e-8', 'time_step=1.0e-6')
    diffusion_case = diffusion_case.replace('\n/', "\n domain_upper=1e-6,1e-6,1e-6, automatic_timestep=.true., "
                                            "checkpoint_file='diffusion.chk', checkpoint_interval_steps=1, "
                                            "stop_after_checkpoint=.true.,\n/")
    log = run('auto_diffusion.nml', diffusion_case, True)
    assert 'Rejected trials this invocation: 0' in log
    diffusion_dt = float((work / 'diffusion.chk').read_text().splitlines()[3].split()[1])
    log = run('auto_diffusion_half.nml', diffusion_case.replace('\n/', '\n transport_cfl=0.175,\n/'), True)
    assert 'Rejected trials this invocation: 0' in log
    half_dt = float((work / 'diffusion.chk').read_text().splitlines()[3].split()[1])
    assert math.isclose(half_dt, diffusion_dt / 2, rel_tol=1e-12)
    coupled_subcycle = options(subcycling='.true.', automatic_timestep='.true.')
    run('subcycle_coupled.nml', coupled_subcycle, True)
    subcycle_saved = output.read_bytes(), budget_path.read_bytes()
    stop_subcycle = coupled_subcycle.replace('\n/', "\n checkpoint_file='subcycle.chk', "
                                            "checkpoint_interval_steps=1, stop_after_checkpoint=.true.,\n/")
    run('subcycle_stop.nml', stop_subcycle, True)
    assert (work / 'subcycle.chk').read_text().splitlines()[0].strip() == 'PELEF_EB_TWO_LEVEL_CHECKPOINT_V4'
    resume_subcycle = coupled_subcycle.replace('\n/', "\n restart_file='subcycle.chk',\n/")
    run('subcycle_resume.nml', resume_subcycle, True)
    assert subcycle_saved == (output.read_bytes(), budget_path.read_bytes())
    run('subcycle_wrong_mode.nml', resume_subcycle.replace('subcycling=.true.', 'subcycling=.false.'), False)
    assert subcycle_saved == (output.read_bytes(), budget_path.read_bytes())
    # Fine-grid hydro constraint converts to a coarse interval by ratio.
    subcycle_hydro = auto_stop.replace('\n/', '\n subcycling=.true.,\n/')
    run('subcycle_hydro_stop.nml', subcycle_hydro, True)
    subcycle_dt = float((work / 'auto.chk').read_text().splitlines()[3].split()[1])
    assert auto_dt < subcycle_dt <= 2 * auto_dt * (1 + 1e-12)
    # Fine transport uses ratio squared, still bounded by the coarse estimate.
    subcycle_diffusion = diffusion_case.replace('\n/', '\n subcycling=.true.,\n/')
    run('subcycle_diffusion_stop.nml', subcycle_diffusion, True)
    subcycle_diffusion_dt = float((work / 'diffusion.chk').read_text().splitlines()[3].split()[1])
    assert diffusion_dt < subcycle_diffusion_dt <= 4 * diffusion_dt * (1 + 1e-12)
    print('PASS: two-level EB application, budgets, restart, retries, automatic intervals and subcycled coupling')


if __name__ == '__main__':
    main()
