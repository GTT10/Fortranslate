#!/usr/bin/env python3
"""Compare production tree inlet halos with compiled pinned AMReX kernels.

Eight supplied components and the auxiliary energy limiter are compared.
This is a coarse interpolation control,
not a full PeleC auxiliary-state, EOS, physical-boundary trajectory validation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('native', 'probe', 'work'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    native = args.native.resolve(strict=True)
    probe = args.probe.resolve(strict=True)
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='amr_boundary_interp_', dir=args.work.resolve()))
    print(work, flush=True)
    with (work/'input.txt').open('w') as stream, (work/'probe.log').open('w') as log:
        subprocess.run([str(probe)], stdout=stream, stderr=log, check=True, timeout=120)
    with (work/'input.txt').open() as stream, (work/'native.txt').open('w') as output:
        subprocess.run([str(native)], stdin=stream, stdout=output, stderr=subprocess.STDOUT,
                       check=True, timeout=120)
    def parse(path):
        cases = []
        for line in path.read_text().splitlines():
            if not line.startswith('CASE '):
                continue
            words = line.split()
            assert len(words) == 18, line
            values = tuple(float(value) for value in words[2:])
            errors = tuple(abs(values[c]-values[c+1])/max(1, abs(values[c+1]))
                           for c in range(0, 16, 2))
            cases.append({'case': int(words[1]), 'actual': values[::2], 'reference': values[1::2],
                          'scaled_errors': errors, 'passed': max(errors) <= 64*sys.float_info.epsilon})
        assert [case['case'] for case in cases] == list(range(1, 49))
        return cases

    cases = parse(work/'native.txt')
    controls = {}
    for name in ('no_aux', 'wrong_time'):
        input_file = work/'input.txt'
        if name == 'no_aux':
            input_file = work/'no_aux_input.txt'
            with input_file.open('w') as stream, (work/'no_aux_probe.log').open('w') as log:
                subprocess.run([str(probe), 'no_aux'], stdout=stream, stderr=log, check=True, timeout=120)
        output_file = work/(name + '.txt')
        command = [str(native)] + (['wrong_time'] if name == 'wrong_time' else [])
        with input_file.open() as stream, output_file.open('w') as output:
            subprocess.run(command, stdin=stream, stdout=output, stderr=subprocess.STDOUT,
                           check=True, timeout=120)
        control_cases = parse(output_file)
        maximum = max(max(case['scaled_errors']) for case in control_cases)
        controls[name] = {'cases': control_cases, 'max_scaled_error': maximum,
                          'detected': maximum > 64*sys.float_info.epsilon}
    checks = {'scope': __doc__, 'native_sha256': hashlib.sha256(native.read_bytes()).hexdigest(),
              'probe_sha256': hashlib.sha256(probe.read_bytes()).hexdigest(),
              'input_sha256': hashlib.sha256((work/'input.txt').read_bytes()).hexdigest(),
              'cases': cases, 'bound': 64*sys.float_info.epsilon,
              'max_scaled_error': max(max(case['scaled_errors']) for case in cases),
              'negative_controls': controls,
              'passed': all(case['passed'] for case in cases) and all(v['detected'] for v in controls.values())}
    (work/'checks.json').write_text(json.dumps(checks, indent=2)+'\n')
    print(f"48 cases: maximum scaled difference {checks['max_scaled_error']:.17g}", flush=True)
    assert checks['passed'], work/'checks.json'


if __name__ == '__main__':
    main()
