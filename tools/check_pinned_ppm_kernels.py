#!/usr/bin/env python3
"""Compare PeleF scalar PPM kernels to compiled, unchanged pinned PeleC bodies.

This checks reconstruction/profile integration only, not EOS or complete CTU.
The reference checkout is read at the manifest commit, independent of its HEAD.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess


def run(command: list[str], **kwargs) -> str:
    return subprocess.run(command, check=True, text=True, capture_output=True,
                          timeout=120, **kwargs).stdout


def function(source: str, name: str) -> str:
    start = source.index('\n' + name + '(') + 1
    opening = source.index('{', start)
    depth = 1
    end = opening + 1
    while depth:
        depth += (source[end] == '{') - (source[end] == '}')
        end += 1
    return 'void ' + source[start:end]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--cxx', default='c++')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    pin = json.loads((root / 'references/pelec_baseline.json').read_text())['pelec']['commit']
    header = run(['git', '-C', str(args.reference), 'show', f'{pin}:Source/PPM.H'])
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    (work / 'PPM.H').write_text(header)
    prefix = '''#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
namespace amrex { using Real=double; using std::min; using std::max; }
constexpr int im2=0, im1=1, i0=2, ip1=3, ip2=4;
'''
    driver = '''
int main() {
  double s[5], flat, u, c, dt, sm, sp, ip[3], im[3];
  std::cout << std::setprecision(17);
  while (std::cin >> s[0]) {
    for (int j=1;j<5;j++) std::cin >> s[j];
    if (!(std::cin >> flat >> u >> c >> dt)) return 2;
    ppm_reconstruct(s,flat,sm,sp);
    ppm_int_profile(sm,sp,s[2],u,c,dt,ip,im);
    std::cout << sm << ' ' << sp;
    for (double x:ip) std::cout << ' ' << x;
    for (double x:im) std::cout << ' ' << x;
    std::cout << '\\n';
  }
}
'''
    cpp = work / 'reference.cpp'
    cpp.write_text(prefix + function(header, 'ppm_reconstruct') + '\n' +
                   function(header, 'ppm_int_profile') + driver)
    compile_command = [args.cxx, '-std=c++17', '-O0', '-ffp-contract=off', str(cpp), '-o', str(work / 'reference')]
    run(compile_command)
    rng = random.Random(20260925)
    rows = []
    for scale in (1e-8, 1., 1e6):
        shapes = [[1., 2., 3., 4., 5.], [5., 4., 3., 2., 1.],
                  [1., 1., 1., 2., 2.], [1., 2., 3., 2., 1.],
                  [2., 2., 2., 2., 2.]]
        shapes += [[rng.uniform(-2, 3) for _ in range(5)] for _ in range(10)]
        for shape in shapes:
            for flat in (0., .37, 1.):
                for u in (-2., -1., 0., 1., 2.):
                    rows.append([v * scale for v in shape] + [flat, u, 1., .2])
    inputs = ''.join(' '.join(format(v, '.17g') for v in row) + '\n' for row in rows)
    (work / 'inputs.txt').write_text(inputs)
    reference = run([str(work / 'reference')], input=inputs)
    actual = run([str(args.executable.resolve()), '--reference-stream'], input=inputs)
    (work / 'reference.txt').write_text(reference)
    (work / 'actual.txt').write_text(actual)
    expected = [[float(v) for v in row.split()] for row in reference.splitlines()]
    observed = [[float(v) for v in row.split()] for row in actual.splitlines()]
    assert len(expected) == len(observed) == len(rows)
    worst = 0.
    # Scale by the input stencil magnitude: an edge can legitimately cross zero.
    for index, (a, b, inputs_row) in enumerate(zip(expected, observed, rows)):
        assert len(a) == len(b) == 8
        assert all(math.isfinite(v) for v in a + b), index
        scale = max(abs(x) for x in inputs_row[:5])
        error = max(abs(x-y) / scale for x, y in zip(a, b))
        assert error <= 3e-14, (index, error, inputs_row, a, b)
        worst = max(worst, error)
    report = dict(reference_commit=pin, header_sha256=hashlib.sha256(header.encode()).hexdigest(),
                  executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
                  cases=len(rows), values=8*len(rows), max_scaled_error=worst,
                  tolerance=3e-14, compile_command=compile_command,
                  compiler=run([args.cxx, '--version']).splitlines()[0], status='PASS')
    (work / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
