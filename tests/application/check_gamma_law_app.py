"""Calorically perfect HIT/SDC/LES application and independent ideal-gas equations."""
import argparse
import csv
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('exe', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='gamma_', dir=args.work.resolve()))
    velocity = work / 'velocity.csv'
    velocity.write_text('# x,y,z,u,v,w\n' + '\n'.join(
        f'{x},{y},{z},{y},{z},{x}' for z in (.25, .75) for y in (.25, .75) for x in (.25, .75)) + '\n')
    base = f"""&reactive_3d
nx=4,ny=4,nz=4,x_upper=1,y_upper=1,z_upper=1
final_time=1e-8,maximum_steps=100,time_integrator='sdc',sdc_iterations=2
problem='uniform_reactor',thermo_model='gamma_law'
reconstruction='characteristic_ppm',riemann_solver='pelec',boundary_condition='periodic'
chemistry_enabled=.false.,transport_enabled=.true.
species_diffusion_enabled=.false.,barodiffusion_enabled=.false.
initial_temperature=300,initial_pressure=101300
initial_velocity_x=0,initial_velocity_y=0,initial_velocity_z=0
output_file='gamma.csv'
/
&reactive_hit input_file='{velocity}',inres=2,mach_t0=.1,uin_norm=1 /
&reactive_constant_transport enabled=.true.,mode='hit',lambda0=.5,reynolds_lambda0=100,prandtl=.71 /
&reactive_les enabled=.true.,model='wale',cs=.5 /
"""
    env = dict(os.environ)
    for key in ('LD_LIBRARY_PATH', 'LIBRARY_PATH', 'DISPLAY'):
        env.pop(key, None)

    def run(label, text, error=None):
        directory = work / label
        directory.mkdir()
        inp = directory / 'input.nml'
        inp.write_text(text)
        result = subprocess.run([str(args.exe.resolve()), str(inp)], cwd=directory, env=env,
                                capture_output=True, text=True, timeout=120)
        log = result.stdout + result.stderr
        (directory / 'run.log').write_text(log)
        output = directory / 'gamma.csv'
        if error:
            assert result.returncode != 0 and error in log, log
            assert not output.exists()
            return
        assert result.returncode == 0, log
        residual = float(re.search(r'Maximum normalized Euler conserved change:\s+(\S+)', log).group(1))
        assert residual < 3e-12, residual
        return output, log

    maximum = 0
    gas_r = 8314.46261815324 / 28.97
    for gamma, model in ((1.4, 'wale'), (1.67, 'wale'), (1.4, 'smagorinsky'), (1.4, 'vreman'), (1.4, 'dynamic')):
        label = str(gamma) if model == 'wale' else model
        content = base.replace("model='wale'", f"model='{model}'")
        output, log = run(label, content + f'&gamma_law gamma={gamma} /\n')
        cp = gas_r * gamma / (gamma - 1)
        rho0 = 101300 / (gas_r * 300)
        rms = .1 * math.sqrt(gamma * gas_r * 300) / math.sqrt(3)
        mu = rho0 * rms * .5 / 100
        expected = mu, mu * cp / .71
        actual = map(float, re.search(r'Constant viscosity/conductivity:\s+(\S+)\s+(\S+)', log).groups())
        assert max(abs(a-b)/b for a, b in zip(actual, expected)) < 3e-14
        with output.open() as stream:
            for row in csv.DictReader(stream):
                rho, pressure, temperature, total = (float(row[k]) for k in ('rho', 'pressure', 'temperature', 'rhoE'))
                kinetic = .5 * rho * sum(float(row[k])**2 for k in ('u', 'v', 'w'))
                error = max(abs(pressure-rho*gas_r*temperature)/pressure,
                            abs(total-kinetic-pressure/(gamma-1))/abs(total))
                maximum = max(maximum, error)
                assert error < 3e-13, row
                assert float(row['Y_AIR']) == 1
    cgs_velocity = work / 'velocity_cm.csv'
    cgs_velocity.write_text('# x,y,z,u,v,w\n' + '\n'.join(
        f'{100*x},{100*y},{100*z},{y},{z},{x}'
        for z in (.25, .75) for y in (.25, .75) for x in (.25, .75)) + '\n')
    converted, _ = run('coordinate_conversion', base.replace(str(velocity), str(cgs_velocity)).replace(
        'inres=2,mach_t0', 'inres=2,coordinate_scale=.01,mach_t0'))
    assert converted.read_bytes() == (work / '1.4/gamma.csv').read_bytes()
    run('invalid_coordinate_scale', base.replace('inres=2,mach_t0', 'inres=2,coordinate_scale=0,mach_t0'),
        'Invalid HIT velocity input')
    default, _ = run('default', base)
    assert default.read_bytes() == (work / '1.4/gamma.csv').read_bytes()
    run('invalid_gamma', base + '&gamma_law gamma=1 /\n', 'Invalid gamma-law heat-capacity ratio')
    run('chemistry', base.replace('chemistry_enabled=.false.', 'chemistry_enabled=.true.'),
        'Gamma-law requires nonreacting SDC')
    run('diffusion', base.replace('species_diffusion_enabled=.false.', 'species_diffusion_enabled=.true.'),
        'Gamma-law requires species diffusion disabled')
    run('molecular', re.sub(r'&reactive_constant_transport[^\n]*\n', '', base),
        'Gamma-law transport requires explicit constant transport')
    print(f'PASS: gamma-law HIT/SDC/WALE, ideal-gas state and transport scales; maximum EOS residual {maximum:.17g}')
    print(work)


if __name__ == '__main__':
    main()
