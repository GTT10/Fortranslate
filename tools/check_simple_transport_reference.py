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
  egtransetEPS(p.eps.data()); egtransetZROT(p.zrot.data()); egtransetNLIN(p.nlin.data());
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
    parser.add_argument('--bulk', action='store_true', help='Also check the fitted bulk viscosity; requires --bundle')
    args = parser.parse_args()
    if args.bulk and not args.bundle:
        parser.error('--bulk requires --bundle')
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
    cpp = CPP
    if args.bulk:
        cpp = cpp.replace('transport(false,true,true,true,false,','transport(true,true,true,true,false,')
        cpp = cpp.replace("    out << '\\n';", "    out << ' ' << xi*.1 << '\\n';")
    (out/'reference.cpp').write_text(cpp)
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
               f"    load_thermo => {bundle['thermo_loader_name']}, load_transport => {bundle['transport_loader_name']}\n"
               '  use nasa7_thermo_mod, only: nasa7_species\n'
               '  use gas_transport_mod, only: gas_transport_species,load_gas_transport_data,set_simple_transport_fits\n'
               '  use mixture_transport_mod, only: mixture_transport_coefficients,binary_diffusion_coefficients\n')
        fortran = fortran.replace('  implicit none', use+'  implicit none')
        fortran = fortran.replace('  integer :: ios,i,j',
                                 '  type(nasa7_species), allocatable :: species(:)\n'
                                 '  type(gas_transport_species), allocatable :: records(:)\n'
                                 '  character(len=24) :: names(3)\n'
                                 '  integer :: geometries(3)\n'
                                 '  real(dp) :: coefficients(4,12),values(5,3),rho\n  integer :: ios,i,j')
        begin = fortran.index("  open(10,file='fits.dat'")
        end = fortran.index('  close(10)',begin)+len('  close(10)')
        fortran = fortran[:begin]+'''  call load_thermo(species,ok)
  if(.not.ok) error stop 'Generated thermo failed'
  mw=species%molecular_weight
  call load_fits(coefficients,ok)
  if(.not.ok) error stop 'Generated fits failed'
  fm=coefficients(:,1:3);fl=coefficients(:,4:6);fd=coefficients(:,7:12)
  call load_transport(names,geometries,values,ok)
  if(.not.ok) error stop 'Generated transport failed'
  call load_gas_transport_data(names,geometries,values,records,ok)
  if(.not.ok) error stop 'Transport records failed'
  call set_simple_transport_fits(records,coefficients,ok)
  if(.not.ok) error stop 'Fit attachment failed'
'''+fortran[end:]
        fortran = fortran.replace(
            '    call pelephysics_simple_transport(mw,y,t,p,fm,fl,fd,mu,k,rd,db,ok)',
            '''    call mixture_transport_coefficients(species,records,y,t,p,mu,k,rd,ok)
    if(.not.ok) error stop 'Record transport failed'
    rho=p/(8314.46261815324_dp*t*sum(y/mw))
    rd=rd*rho
    call binary_diffusion_coefficients(species,records,t,p,db,ok)''',1)
        bundle_hash = hashlib.sha256(args.bundle.read_bytes()).hexdigest()
    if args.bulk:
        fortran = fortran.replace('  implicit none',
            '  use mixture_transport_mod, only: simple_bulk_viscosity\n'
            '  use reactive_transport_3d_mod, only: reactive_transport_face_flux_3d,reactive_transport_timestep_3d\n'
            '  use reactive_transport_2d_mod, only: face_transport_data\n'
            '  use gas_transport_mod, only: write_simple_transport_fits,check_simple_transport_fits\n'
            '  use reactive_1d_mod, only: reactive_primitive_to_conserved\n  implicit none')
        fortran = fortran.replace('  integer :: ios,i,j',
            '  real(dp) :: bulk,q(8),gradient(3,3),flux(8),species_energy,dt,dmax,sound\n'
            '  real(dp) :: state(8,2,2,2),temperatures(2,2,2)\n'
            '  real(dp) :: tl,tr,pl,pr,yl(3),yr(3),ql(8),qr(8),fmu,fk,df(3),ol(3),orr(3),yf(3),xl(3),xr(3),hf(3)\n'
            '  integer :: ios,i,j,ii,jj,kk,record_status')
        fortran = fortran.replace("    write(11,'(5es26.17)') mu,k,rd",'''    call simple_bulk_viscosity(species,records,y,t,bulk,ok)
    if(.not.ok) error stop 'Bulk viscosity failed'
    write(11,'(6es26.17)') mu,k,rd,bulk''')
        fortran = fortran.replace('  y=[-1.0_dp,1.0_dp,1.0_dp]', '''  q=[rho,1._dp,2._dp,3._dp,p,y]
  gradient=0
  do i=1,3
    gradient(i,i)=1
  end do
  call reactive_transport_face_flux_3d(species,records,q,q,t,t,1._dp,1,gradient, &
    .true.,.false.,.false.,.false.,flux,species_energy,ok)
  if(.not.ok) error stop 'Bulk face flux failed'
  if(abs(flux(2)+3*bulk)>1e-12_dp*bulk) error stop 'Isotropic compression stress mismatch'
  if(abs(flux(5)+3*bulk)>1e-12_dp*bulk) error stop 'Bulk viscous work mismatch'
  do kk=1,2
    do jj=1,2
      do ii=1,2
        temperatures(ii,jj,kk)=t
        call reactive_primitive_to_conserved(species,q,state(:,ii,jj,kk),temperatures(ii,jj,kk),sound,ok)
        if(.not.ok) error stop 'Timestep state failed'
      end do
    end do
  end do
  call reactive_transport_timestep_3d(species,records,state,temperatures,2,2,2,1._dp,1._dp,1._dp,.5_dp, &
    .true.,.false.,.false.,dt,dmax,ok)
  if(.not.ok.or.dt<=0) error stop 'Bulk timestep failed'
  if(abs(dmax-(4*mu/3+bulk)/rho)>1e-12_dp*dmax) error stop 'Bulk timestep omitted dilatation viscosity'
  open(12,file='states.dat',status='old')
  open(13,file='candidate-faces.dat',status='new')
  do i=1,3
    read(12,*) tl,pl,yl
    read(12,*) tr,pr,yr
    ql=[pl/(8314.46261815324_dp*tl*sum(yl/mw)),0._dp,0._dp,0._dp,pl,yl]
    qr=[pr/(8314.46261815324_dp*tr*sum(yr/mw)),0._dp,0._dp,0._dp,pr,yr]
    call face_transport_data(species,records,ql,qr,tl,tr,fmu,fk,df,ol,orr,yf,xl,xr,hf,ok)
    if(.not.ok) error stop 'Nonuniform face failed'
    write(13,'(5es26.17)') fmu,fk,.5_dp*(ql(1)+qr(1))*df
  end do
  close(12);close(13)
  open(14,status='scratch',form='formatted')
  call write_simple_transport_fits(14,records(1),record_status)
  if(record_status/=0) error stop 'Fit record write failed'
  rewind(14)
  call check_simple_transport_fits(14,records(1),record_status)
  if(record_status/=0) error stop 'Fit record round trip failed'
  rewind(14)
  write(14,'(a)') 'PELEF_SIMPLE_TRANSPORT_V2'
  rewind(14)
  call check_simple_transport_fits(14,records(1),record_status)
  if(record_status==0) error stop 'Old face model accepted'
  close(14)
  y=[-1.0_dp,1.0_dp,1.0_dp]''')
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
    if args.bulk:
        commands[2][1:1] = ['-I'+str(candidate/'selected_reactive_3d_private_modules'),
                            '-I'+str(candidate/'selected_reactive_2d_runtime_modules')]
        index = commands[2].index(str(candidate/'libpelef_selected_reactive_1d_runtime.a'))
        commands[2][index:index] = [str(candidate/'libpelef_selected_reactive_3d_runtime.a'),
                                   str(candidate/'libpelef_selected_reactive_2d_runtime.a')]
    with (out/'commands.log').open('w') as log:
        for command in commands:
            log.write(shlex.join(command)+'\n'); log.flush()
            subprocess.run(command,cwd=out,stdout=log,stderr=log,check=True)
    read = lambda name: [list(map(float,l.split())) for l in (out/name).read_text().splitlines()]
    expected, actual = read('reference.dat'), read('candidate.dat')
    if len(expected)!=len(states) or len(actual)!=len(states):
        raise ValueError('Incomplete output')
    for row in expected+actual:
        if len(row)!=(6 if args.bulk else 5) or not all(math.isfinite(v) and v>0 for v in row):
            raise ValueError('Invalid transport output')
    errors = [max(abs(a-b)/abs(a) for a,b in zip(r,c)) for r,c in zip(expected,actual)]
    face_error = None
    if args.bulk:
        faces = read('candidate-faces.dat')
        if len(faces)!=3 or any(len(row)!=5 for row in faces):
            raise ValueError('Incomplete face output')
        face_errors = []
        for i, row in enumerate(faces):
            if not all(math.isfinite(v) and v>0 for v in row):
                raise ValueError('Invalid face transport')
            ref = [2*a*b/(a+b) for a,b in zip(expected[2*i][:5],expected[2*i+1][:5])]
            face_errors.extend(abs(a-b)/abs(a) for a,b in zip(ref,row))
        face_error = max(face_errors)
        if face_error>2e-12:
            raise ValueError(f'Face transport mismatch: {face_error}')
    report = dict(states=states, relative_errors=errors, maximum_relative_error=max(errors),
                  tolerance=2e-12, commands=commands,
                  reference_identity=identity, source_sha256=source_hashes,
                  bundle_sha256=bundle_hash,
                  bulk_viscosity=args.bulk,
                  bulk_stress_and_timestep_checks=args.bulk,
                  face_maximum_relative_error=face_error,
                  fitted_checkpoint_record_checks=args.bulk,
                  fits_sha256=hashlib.sha256((out/'fits.dat').read_bytes()).hexdigest())
    (out/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    if max(errors)>2e-12:
        raise ValueError(f'Transport mismatch: {max(errors)}')
    print(f'PASS: {len(states)} states, one invalid-state rejection; maximum error {max(errors):.3e}')


if __name__=='__main__':
    main()
