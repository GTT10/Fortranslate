#!/usr/bin/env python3
"""Require ETAB size reduction, inventory conservation and exact restart."""
import argparse
import csv
from pathlib import Path
import subprocess
import tempfile


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='etab-', dir=args.work_dir.resolve()))
    base = (Path(__file__).resolve().parents[1] / 'cases/spray_3d/etab.nml').read_text()
    for mode in ('periodic', 'physical'):
        def run(name, extra=''):
            tag = mode + '_' + name
            text = base.replace("output_prefix='etab'", f"output_prefix='{tag}'")
            if mode == 'physical':
                extra += "boundary_faces=6*'slip_wall',"
            text = text.replace('&spray', '&spray\n' + extra)
            path = work / (tag + '.nml')
            path.write_text(text)
            result = subprocess.run([str(args.executable.resolve()), str(path)], cwd=work,
                                    capture_output=True, text=True, timeout=240)
            (work / (tag + '.log')).write_text(result.stdout + result.stderr)
            if result.returncode:
                raise RuntimeError(result.stdout + result.stderr)
            return tag
        reference = run('reference')
        run('stop', f"stop_after_steps=2,checkpoint_file='{mode}.chk',")
        resumed = run('resumed', f"restart_file='{mode}.chk',")
        for suffix in ('.gas.csv', '.parcels.csv'):
            assert (work / (reference + suffix)).read_bytes() == (work / (resumed + suffix)).read_bytes()
        drops = rows(work / (reference + '.parcels.csv'))
        assert all(float(p['diameter']) < 20e-6 for p in drops)
        assert abs(sum(float(p['mass']) * float(p['multiplicity']) for p in drops) / 1e-10 - 1) < 1e-12
        row = rows(work / (reference + '.history.csv'))[-1]
        assert row == rows(work / (resumed + '.history.csv'))[-1]
        assert abs(float(row['energy_error'])) < 1e-10 * abs(float(row['total_energy']))
        if mode == 'periodic':
            bad = (work / (resumed + '.nml')).read_text().replace('surface_tension=.072', 'surface_tension=.073')
            bad = bad.replace(resumed, 'changed_surface_tension')
            bad_path = work / 'changed_surface_tension.nml'
            bad_path.write_text(bad)
            result = subprocess.run([str(args.executable.resolve()), str(bad_path)], cwd=work,
                                    capture_output=True, text=True, timeout=30)
            (work / 'changed_surface_tension.log').write_text(result.stdout + result.stderr)
            assert result.returncode != 0 and 'incompatible spray checkpoint' in result.stdout + result.stderr
        print(mode, 'PASS: ETAB breakup, inventory, exact restart')
    print('artifacts:', work)


if __name__ == '__main__':
    main()
