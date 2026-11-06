#!/usr/bin/env python3
"""Compare independent dual-energy policy with an unchanged pinned reset body.

The EOS adapter is e=cv*T+formation. This isolates reset policy, including
negative formation energy and frame potential; it does not qualify NASA7,
state transport, restart or full solver trajectories.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from check_pinned_ppm_kernels import function, run

PREFIX = r'''
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
namespace amrex {
 using Real=double; using std::max;
 [[noreturn]] void Error(const char* s) { throw std::runtime_error(s); }
 template<class T> struct Array4 {
  T* v; T& operator()(int,int,int,int c) const { return v[c]; }
 };
}
constexpr int NUM_SPECIES=1,URHO=0,UMX=1,UMY=2,UMZ=3,UEDEN=4,UEINT=5,UTEMP=6,UFS=7;
double cv=10,formation=0;
namespace constants { double very_small_num(){ return std::numeric_limits<double>::epsilon()*1e-100; } }
namespace pele::physics {
 namespace eos { struct SRK {}; }
 struct Ideal { void T2Ei(double t,double* e) { e[0]=cv*t+formation; } };
 struct PhysicsType { using eos_type=Ideal; static Ideal eos(){ return {}; } };
}
'''
MAIN = r'''
int main() {
 std::cout << std::setprecision(17);
 double rho,mx,my,mz,total,internal,temp,eta,phi; int small,negative,update;
 while(std::cin >> rho >> mx >> my >> mz >> total >> internal >> temp >> small >> negative >> update
       >> eta >> phi >> formation >> cv) {
  double v[8]={rho,mx,my,mz,total,internal,temp,rho};
  amrex::Array4<double> state{v};
  // omega=sqrt(2*phi), rad=1 makes the source potential term equal phi.
  pc_rst_int_e(0,0,0,state,small,negative,update,eta,std::sqrt(2*phi),1);
  std::cout << v[4] << ' ' << v[5] << ' ' << v[4]-total << '\n';
 }
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--cxx', default='/usr/bin/g++')
    args = parser.parse_args()
    header = args.reference/'Source/PeleCUtilities.H'
    # Byte-pinned PeleC bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
    expected = '2a77e64ce3cf3cb4eb00716ae6de66b9e470f007a040c0c812878142ea6e940d'
    assert hashlib.sha256(header.read_bytes()).hexdigest() == expected
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='dual_energy_reference_', dir=args.work.resolve()))
    print(work, flush=True)
    source = work/'reference.cpp'
    source.write_text(PREFIX + function(header.read_text(), 'pc_rst_int_e') + MAIN)
    native = work/'reference'
    run([args.cxx, '-std=c++20', '-O2', str(source), '-o', str(native)])
    cases = []
    for flags, phi, formation, eta, pair in itertools.product(
            itertools.product((0, 1), repeat=3), (0, 8), (0, -100), (.01, .1),
            ((200, 20), (102, 20), (100, 20), (100, -1), (0, 20), (-210, -220))):
        total, internal = pair
        # Phi=8 is exact for omega=4 and rad=1 in the pinned expression.
        cases.append([2, 20, 0, 0, total, internal, 3, *flags, eta, phi, formation, 10])
    data = ''.join(' '.join(str(v) for v in row)+'\n' for row in cases)
    (work/'input.txt').write_text(data)
    def calculate(command, name):
        output = run(command, input=data)
        (work/(name+'.txt')).write_text(output)
        rows = [tuple(map(float, line.split())) for line in output.splitlines()]
        assert len(rows) == len(cases) and all(len(row) == 3 for row in rows)
        return rows
    actual = calculate([str(args.probe.resolve(strict=True))], 'fortran')
    reference = calculate([str(native)], 'reference')
    errors = [max(abs(a-b)/max(1, abs(b)) for a, b in zip(u, v)) for u, v in zip(actual, reference)]
    checks = {'scope': __doc__, 'source_sha256': expected, 'cases': len(cases),
              'max_scaled_error': max(errors), 'bound': 64*sys.float_info.epsilon,
              'probe_sha256': hashlib.sha256(args.probe.read_bytes()).hexdigest(),
              'comparisons': [{'input': c, 'actual': a, 'reference': r, 'error': e}
                              for c, a, r, e in zip(cases, actual, reference, errors)]}
    (work/'checks.json').write_text(json.dumps(checks, indent=2)+'\n')
    print(len(cases), 'cases: maximum scaled difference', max(errors), flush=True)
    assert max(errors) <= checks['bound'], work/'checks.json'


if __name__ == '__main__':
    main()
