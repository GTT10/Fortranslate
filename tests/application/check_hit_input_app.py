"""HIT interpolation oracle, file formats, LES stepping and fixed/selected identity."""
import argparse
import csv
import os
from pathlib import Path
import re
import struct
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('fixed', 'selected', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--constant-transport', action='store_true')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='hit_', dir=args.work.resolve()))
    rows = [(x, y, z, 2*x+y, y+3*z, z+x)
            for z in (.25, .75) for y in (.25, .75) for x in (.25, .75)]
    ascii_file = work / 'velocity.csv'
    ascii_file.write_text('# x,y,z,u,v,w\n' + '\n'.join(','.join(map(str, r)) for r in rows) + '\n')
    binary_file = work / 'velocity.bin'
    binary_file.write_bytes(b''.join(struct.pack('=6d', *r) for r in rows))
    env = dict(os.environ)
    for key in ('LD_LIBRARY_PATH', 'LIBRARY_PATH', 'DISPLAY'):
        env.pop(key, None)

    def run(label, path, binary=False, selected=False, end=1e-20, les=False, size=2, error=False, coefficients=None, boundary=None):
        directory = work / label
        directory.mkdir()
        composition = ("thermo_model='selected',composition_count=1,composition_species='N2',"
                       "composition_mole_fractions=1" if selected else
                       "thermo_model='full_h2o2',x_h2=0,x_o2=0,x_h=0,x_o=0,x_oh=0,x_h2o=0,x_ho2=0,x_h2o2=0,x_ar=0,x_n2=1")
        text = f"""&reactive_3d
nx=8,ny=4,nz=4,x_upper=2,y_upper=1,z_upper=1
final_time={end},maximum_steps=100,time_integrator='sdc',sdc_iterations=2
problem='uniform_reactor',{composition}
reconstruction='characteristic_ppm',riemann_solver='pelec',boundary_condition='periodic'
chemistry_enabled=.false.,transport_enabled=.true.
initial_temperature=400,initial_pressure=101325
initial_velocity_x=3,initial_velocity_y=2,initial_velocity_z=1
output_file='hit.csv'
/
&reactive_hit input_file='{path}',inres={size},binfmt={'.true.' if binary else '.false.'},mach_t0=.1,uin_norm=1 /
"""
        if args.constant_transport:
            values = coefficients if coefficients is not None else ('.001', '.5')
            text += f"&reactive_constant_transport enabled=.true.,viscosity={values[0]},conductivity={values[1]} /\n"
        if boundary:
            text = text.replace("boundary_condition='periodic'", f"boundary_condition='{boundary}'")
        if les:
            text += "&reactive_les enabled=.true.,model='wale' /\n"
        inp = directory / 'input.nml'
        inp.write_text(text)
        exe = args.selected if selected else args.fixed
        result = subprocess.run([str(exe.resolve()), str(inp)], cwd=directory, env=env,
                                capture_output=True, text=True, timeout=120)
        log = result.stdout + result.stderr
        (directory / 'run.log').write_text(log)
        output = directory / 'hit.csv'
        if error:
            expected_error = 'Invalid HIT velocity input' if error is True else error
            assert result.returncode != 0 and expected_error in log, log
            assert not output.exists()
            return
        assert result.returncode == 0, log
        residual = float(re.search(r'Maximum normalized Euler conserved change:\s+(\S+)', log).group(1))
        assert residual < 3e-12, residual
        return output, log

    actual, log = run('ascii', ascii_file)
    binary, _ = run('binary', binary_file, binary=True)
    selected, _ = run('selected', ascii_file, selected=True)
    assert actual.read_bytes() == binary.read_bytes() == selected.read_bytes()
    scale = float(re.search(r'HIT velocity scale:\s+(\S+)', log).group(1))
    maximum = 0.0
    with actual.open() as stream:
        for row in csv.DictReader(stream):
            # Reference locate clamps below the first center; its upper neighbor
            # wraps after the last center. This is deliberately not clipped interpolation.
            def mapped(value):
                q = float(value) % 1
                return q if q < .75 else 1.5 - q
            x, y, z = (mapped(row[key]) for key in ('x', 'y', 'z'))
            expected = (3+scale*(2*x+y), 2+scale*(y+3*z), 1+scale*(z+x))
            error = max(abs(float(row[key])-value) for key, value in zip(('u', 'v', 'w'), expected))
            maximum = max(maximum, error)
            assert error < 1e-10, (row, expected)
            assert abs(float(row['temperature'])-400) < 1e-9
    fixed, _ = run('les', ascii_file, end=1e-8, les=True)
    selected, _ = run('les_selected', ascii_file, selected=True, end=1e-8, les=True)
    assert fixed.read_bytes() == selected.read_bytes()
    truncated = work / 'truncated.bin'
    truncated.write_bytes(binary_file.read_bytes()[:-8])
    run('truncated', truncated, binary=True, error=True)
    unordered = work / 'unordered.csv'
    unordered.write_text(ascii_file.read_text().replace('0.75', '0.25'))
    run('unordered', unordered, error=True)
    if args.constant_transport:
        active, _ = run('constant_only', ascii_file, end=1e-8)
        zero, _ = run('zero_transport', ascii_file, end=1e-8, coefficients=('0', '0'))
        assert active.read_bytes() != zero.read_bytes(), 'Constant transport must change the solution'
        run('negative_transport', ascii_file, coefficients=('-1', '0'),
            error='Negative constant transport coefficient')
        run('nonperiodic_transport', ascii_file, boundary='outflow',
            error='Constant transport requires periodic SDC with transport enabled')
    if args.reference:
        fixed, _ = run('reference', args.reference.resolve(), size=32, end=1e-8, les=True)
        selected, _ = run('reference_selected', args.reference.resolve(), size=32, end=1e-8, les=True, selected=True)
        assert fixed.read_bytes() == selected.read_bytes()
    print(f'PASS: interpolation error {maximum:.17g}; ASCII/binary/selected identity; LES and input rejection')
    print(work)


if __name__ == '__main__':
    main()
