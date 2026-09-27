"""Compare the Fortran Simple transport kernel with compiled pinned PelePhysics.

Requires the existing heptane Spray-Evaporation reference and candidate builds.
This checks the transport kernel only, not its application integration.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shlex
import subprocess


CPP = r'''
#include "PelePhysics.H"
#include <fstream>
#include <iomanip>
int main() {
  static_assert(NUM_SPECIES == 3, "Expected heptane_3sp");
  using namespace pele::physics;
  transport::TransParm<eos::Fuego, transport::SimpleTransport> p;
  egtransetWT(p.wt.data());
  egtransetCOFETA(p.fitmu.data()); egtransetCOFLAM(p.fitlam.data());
  egtransetCOFD(p.fitdbin.data());
  for(int i=0;i<NUM_SPECIES;++i) p.iwt[i]=1/p.wt[i];
  std::ofstream fits("fits.dat"); fits << std::setprecision(17);
  for(auto v:p.wt) fits << v << '\n';
  for(auto v:p.fitmu) fits << v << '\n';
  for(auto v:p.fitlam) fits << v << '\n';
  for(auto v:p.fitdbin) fits << v << '\n';
  std::ifstream states("states.dat");
  std::ofstream out("reference.dat"); out << std::setprecision(17);
  double T,P,Y[3],D[3],mu,xi,lam;
  while(states >> T >> P >> Y[0] >> Y[1] >> Y[2]) {
    transport::SimpleTransport::transport(false,true,true,true,false,
      T,1.,Y,D,nullptr,mu,xi,lam,&p);
    out << mu*.1 << ' ' << lam*1e-5;
    for(auto v:D) out << ' ' << v*.1;
    out << '\n';
  }
}
'''

FORTRAN = '''program check
  use precision_mod, only: dp
  use mixture_transport_mod, only: pelephysics_simple_transport
  implicit none
  real(dp) :: mw(3),fm(4,3),fl(4,3),fd(4,6),y(3),t,p,mu,k,rd(3),db(3,3)
  integer :: ios,i,j
  logical :: ok
  open(10,file='fits.dat',status='old')
  do i=1,3
    read(10,*) mw(i)
  end do
  do j=1,3
    do i=1,4
      read(10,*) fm(i,j)
    end do
  end do
  do j=1,3
    do i=1,4
      read(10,*) fl(i,j)
    end do
  end do
  do j=1,6
    do i=1,4
      read(10,*) fd(i,j)
    end do
  end do
  close(10)
  open(10,file='states.dat',status='old')
  open(11,file='candidate.dat',status='new')
  do
    read(10,*,iostat=ios) t,p,y
    if(ios<0) exit
    if(ios/=0) error stop 'Invalid state'
    call pelephysics_simple_transport(mw,y,t,p,fm,fl,fd,mu,k,rd,db,ok)
    if(.not.ok) error stop 'Transport rejected'
    write(11,'(5es26.17)') mu,k,rd
    if(any(db<=0).or.maxval(abs(db-transpose(db)))>0) error stop 'Invalid binary diffusion'
  end do
  y=[-1.0_dp,1.0_dp,1.0_dp]
  call pelephysics_simple_transport(mw,y,t,p,fm,fl,fd,mu,k,rd,db,ok)
  if(ok.or.mu/=0.or.k/=0.or.any(rd/=0).or.any(db/=0)) error stop 'Invalid state accepted'
end program
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-build', type=Path, required=True)
    parser.add_argument('--candidate-build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bundle', type=Path, help='Use imported/generated Fortran coefficients')
    args = parser.parse_args()
    ref, candidate, out = (p.resolve() for p in
                           (args.reference_build, args.candidate_build, args.output))
    out.mkdir(parents=True, exist_ok=False)
    flags_file = ref/'Exec/RegTests/Spray-Evaporation/CMakeFiles/PeleC-Spray-Evaporation.dir/flags.make'
    flags = dict(line.split(' = ', 1) for line in flags_file.read_text().splitlines() if ' = ' in line)
    cache = dict(line.split('=', 1) for line in (ref/'CMakeCache.txt').read_text().splitlines()
                 if line.startswith('CMAKE_HOME_DIRECTORY:'))
    source = Path(cache['CMAKE_HOME_DIRECTORY:INTERNAL'])
    identity = json.loads((source/'pelef-reference-provenance.json').read_text())
    physics = source/'Submodules/PelePhysics'
    source_hashes = {str(path.relative_to(source)): hashlib.sha256(path.read_bytes()).hexdigest()
                     for path in [physics/'Source/Transport/Simple.H',
                                  physics/'Mechanisms/heptane_3sp/mechanism.H']}
    (out/'reference.cpp').write_text(CPP)
    fortran = FORTRAN
    generated_sources = []
    bundle_hash = None
    if args.bundle:
        from check_mechanism_provenance import check
        from generate_elementary_mechanism import validate, generate
        if not check(args.bundle):
            raise ValueError('Bundle provenance failed')
        bundle = json.loads(args.bundle.read_text())
        validate(bundle)
        if bundle['species'] != ['NC7H16','O2','N2'] or 'simple_transport' not in bundle:
            raise ValueError('Expected heptane bundle with Simple transport fits')
        (out/'generated.f90').write_text(generate(bundle))
        generated_sources = ['generated.f90']
        use = (f"  use {bundle['module_name']}, only: load_fits => {bundle['transport_loader_name']}_simple, &\n"
               f"    load_thermo => {bundle['thermo_loader_name']}\n"
               '  use nasa7_thermo_mod, only: nasa7_species\n')
        fortran = fortran.replace('  implicit none', use+'  implicit none')
        fortran = fortran.replace('  integer :: ios,i,j',
                                 '  type(nasa7_species), allocatable :: species(:)\n'
                                 '  real(dp) :: coefficients(4,12)\n  integer :: ios,i,j')
        begin = fortran.index("  open(10,file='fits.dat'")
        end = fortran.index('  close(10)',begin)+len('  close(10)')
        fortran = fortran[:begin]+'''  call load_thermo(species,ok)
  if(.not.ok) error stop 'Generated thermo failed'
  mw=species%molecular_weight
  call load_fits(coefficients,ok)
  if(.not.ok) error stop 'Generated fits failed'
  fm=coefficients(:,1:3);fl=coefficients(:,4:6);fd=coefficients(:,7:12)
'''+fortran[end:]
        bundle_hash = hashlib.sha256(args.bundle.read_bytes()).hexdigest()
    (out/'candidate.f90').write_text(fortran)
    states = [(340,101325,.25,.17,.58), (400,101325,0,.233,.767),
              (1200,3e6,.4,.3,.3), (300,101325,1,0,0),
              (500,101325,0,1,0), (800,101325,0,0,1)]
    (out/'states.dat').write_text(''.join(' '.join(map(str,s))+'\n' for s in states))
    commands = [
        ['/usr/bin/g++','-O2','-std=c++20',*shlex.split(flags['CXX_DEFINES']),
         *shlex.split(flags['CXX_INCLUDES']),'reference.cpp','-o','reference'],
        ['./reference'],
        ['/usr/bin/gfortran','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
         '-I'+str(candidate/'selected_reactive_1d_runtime_modules'),
         '-I'+str(candidate/'selected_runtime_modules'),*generated_sources,'candidate.f90',
         str(candidate/'libpelef_selected_reactive_1d_runtime.a'),
         str(candidate/'libpelef_selected_runtime_core.a'),'-o','candidate'],
        ['./candidate']]
    with (out/'commands.log').open('w') as log:
        for command in commands:
            log.write(shlex.join(command)+'\n'); log.flush()
            subprocess.run(command,cwd=out,stdout=log,stderr=log,check=True)
    read = lambda name: [list(map(float,l.split())) for l in (out/name).read_text().splitlines()]
    expected, actual = read('reference.dat'), read('candidate.dat')
    if len(expected)!=len(states) or len(actual)!=len(states):
        raise ValueError('Incomplete output')
    for row in expected+actual:
        if len(row)!=5 or not all(math.isfinite(v) and v>0 for v in row):
            raise ValueError('Invalid transport output')
    errors = [max(abs(a-b)/abs(a) for a,b in zip(r,c)) for r,c in zip(expected,actual)]
    report = dict(states=states, relative_errors=errors, maximum_relative_error=max(errors),
                  tolerance=2e-12, commands=commands,
                  reference_identity=identity, source_sha256=source_hashes,
                  bundle_sha256=bundle_hash,
                  fits_sha256=hashlib.sha256((out/'fits.dat').read_bytes()).hexdigest())
    (out/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    if max(errors)>2e-12:
        raise ValueError(f'Transport mismatch: {max(errors)}')
    print(f'PASS: {len(states)} states, one invalid-state rejection; maximum error {max(errors):.3e}')


if __name__=='__main__':
    main()
