#!/usr/bin/env python3
"""Compare imported/generated Chebyshev rates and concentration derivatives with Cantera."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    import cantera as ct
    import numpy as np
    parser = argparse.ArgumentParser()
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='chebyshev-', dir=args.work_dir.resolve()))
    root = Path(__file__).resolve().parents[1]
    parent = ct.Solution('h2o2.yaml')
    reaction = ct.Reaction.from_dict({
        'equation': 'H + OH <=> H2O', 'type': 'Chebyshev',
        'temperature-range': [300., 2500.], 'pressure-range': [1e3, 1e7],
        'data': [[9., .4, -.02], [-.3, .08, .005], [.02, -.003, .001]]}, parent)
    gas = ct.Solution(thermo='ideal-gas', kinetics='gas', species=parent.species(),
                      reactions=[reaction], transport_model='mixture-averaged')
    gas.name = 'chebyshev'
    # Export in non-SI units, reload, and compare with the original SI rate.
    gas.write_yaml(work/'chebyshev.yaml', units={'length': 'cm', 'quantity': 'mol', 'pressure': 'atm'})
    reloaded = ct.Solution(work/'chebyshev.yaml')
    assert np.allclose(reloaded.reaction(0).rate.data, reaction.rate.data, rtol=0, atol=1e-14)
    subprocess.run([sys.executable, str(root/'tools/ingest_cantera_mechanism.py'),
        '--input', str(work/'chebyshev.yaml'), '--output', str(work/'chebyshev.json'),
        '--phase', 'chebyshev', '--module-name', 'chebyshev_bundle_mod',
        '--loader-name', 'load_chebyshev'], check=True)
    subprocess.run([sys.executable, str(root/'tools/generate_elementary_mechanism.py'),
        '--input', str(work/'chebyshev.json'), '--output', str(work/'chebyshev.F90')], check=True)
    (work/'probe.F90').write_text('''program probe
 use precision_mod, only: dp
 use elementary_kinetics_mod
 use chebyshev_bundle_mod, only: load_chebyshev
 implicit none
 type(elementary_reaction),allocatable :: r(:)
 real(dp) :: t,c(10),k,dk(10)
 integer :: ios
 logical :: ok
 call load_chebyshev(r,ok)
 if (.not.ok) error stop 'load'
 do
   read(*,*,iostat=ios) t,c
   if(ios/=0) exit
   call reaction_effective_forward_rate(r(1),t,c,k,dk,ok)
   if(.not.ok) error stop 'rate'
   write(*,'(*(es26.17e3,1x))') k,dk
 end do
 c=0
 call reaction_effective_forward_rate(r(1),1000.0_dp,c,k,dk,ok)
 if(ok) error stop 'zero pressure accepted'
 c=0.01_dp
 r(1)%chebyshev_shape=[1,1];r(1)%chebyshev=0;r(1)%chebyshev(1,1)=9
 call reaction_effective_forward_rate(r(1),1000.0_dp,c,k,dk,ok)
 if(.not.ok.or.abs(k/1e9_dp-1)>1e-13_dp.or.any(dk/=0)) error stop 'constant rate'
 r(1)%chebyshev_shape=[9,1]
 call reaction_effective_forward_rate(r(1),1000.0_dp,c,k,dk,ok)
 if(ok) error stop 'oversized table accepted'
 r(1)%chebyshev_shape=[1,1];r(1)%chebyshev_bounds(2)=r(1)%chebyshev_bounds(1)
 call reaction_effective_forward_rate(r(1),1000.0_dp,c,k,dk,ok)
 if(ok) error stop 'degenerate interval accepted'
end program
''')
    subprocess.run(['gfortran', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
        '-I'+str(args.build_dir.resolve()/'modules'), str(work/'chebyshev.F90'), str(work/'probe.F90'),
        str(args.build_dir.resolve()/'libpelef_core.a'), '-o', str(work/'probe')], cwd=work, check=True)
    inputs = []; expected = []; derivatives = []
    for t in [250., 300., 700., 1600., 2500., 3000.]:
        for p in [1e2, 1e3, 3e4, 1e5, 7e6, 1e7, 1e8]:
            gas.TPX = t, p, 'N2:1'
            m = sum(gas.concentrations)
            inputs.append(' '.join(map(str, [t, *gas.concentrations])))
            expected.append(gas.forward_rate_constants[0])
            samples = []; h = 1e-3
            for shift in [-2*h, -h, h, 2*h]:
                gas.TD = t, m*(1+shift)*gas.mean_molecular_weight
                samples.append(gas.forward_rate_constants[0])
            derivatives.append((samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*h*m))
    run = subprocess.run([str(work/'probe')], input='\n'.join(inputs)+'\n',
                         text=True, capture_output=True, check=True)
    (work/'rates.txt').write_text(run.stdout)
    rows = np.array([list(map(float, line.split())) for line in run.stdout.splitlines()])
    assert rows.shape == (len(inputs), 11)
    rate_error = float(np.max(abs(rows[:, 0]/expected-1)))
    derivative_error = float(np.max(abs(rows[:, 1:]/np.array(derivatives)[:, None]-1)))
    assert rate_error < 1e-12, rate_error
    assert derivative_error < 1e-6, derivative_error
    result = {'states': len(inputs), 'rate_relative_error': rate_error,
              'derivative_relative_error': derivative_error, 'artifacts': str(work)}
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
