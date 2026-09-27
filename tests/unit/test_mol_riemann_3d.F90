program test_mol_riemann_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved
  use mol_riemann_3d_mod, only: mol_riemann_flux_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: q(:,:),center(:),jump(:),face(:),state(:),flux(:),expected(:),other_state(:)
  real(dp) :: t,c,sound(2),nan,diffusion_speed
  logical :: ok,connected(2,2)
  integer :: axis,order,i,checks
  checks=0
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  allocate(q(reactive_nprim(size(species)),-1:2))
  allocate(center(size(q,1)),jump(size(q,1)),face(size(q,1)))
  allocate(state(reactive_nvar(size(species))))
  allocate(flux(size(state)),expected(size(state)),other_state(size(state)))
  center=0;center(1:5)=[2._dp,10._dp,-20._dp,30._dp,200000._dp]
  do i=1,size(species)
    if(species(i)%name=='N2') center(5+i)=.7_dp
    if(species(i)%name=='O2') center(5+i)=.3_dp
  end do
  nan=ieee_value(0._dp,ieee_quiet_nan)
  connected=.true.
  do axis=1,3
    do order=1,2
      jump=0
      if(order==2) then
        jump(1)=.02_dp;jump(2:4)=[.2_dp,-.3_dp,.4_dp];jump(5)=100
      end if
      do i=-1,2
        q(:,i)=center+real(i,dp)*jump
      end do
      face=center+.5_dp*jump
      call reactive_primitive_to_conserved(species,face,state,t,c,ok)
      call require(ok,'face EOS')
      sound=c
      expected=state*face(axis+1)
      expected(axis+1)=expected(axis+1)+face(5)
      expected(5)=(state(5)+face(5))*face(axis+1)
      call mol_riemann_flux_3d(species,q,sound,axis,connected,order,flux,ok)
      call require(ok,'MOL flux')
      call require(maxval(abs(flux-expected)/max(1._dp,abs(expected)))<2e-13_dp,'analytic physical face flux')
      call require(abs(sum(flux(6:))-flux(1))<1e-12_dp,'species flux closes mass')
      call mol_riemann_flux_3d(species,q,sound,axis,connected,order,flux,ok,lax_friedrichs=.true.)
      call require(ok,'MOL LF flux')
      call require(maxval(abs(flux-expected)/max(1._dp,abs(expected)))<2e-13_dp,'LF analytic physical face flux')
    end do
    ! Stationary equal-pressure contact: the reference contact speed is zero,
    ! and the LF diffusion speed is the arithmetic mean of EOS sound speeds.
    face=center;face(2:4)=0
    q(:,0)=face
    call reactive_primitive_to_conserved(species,face,state,t,sound(1),ok)
    call require(ok,'LF left stationary contact')
    face(1)=1.5_dp; q(:,1)=face
    call reactive_primitive_to_conserved(species,face,other_state,t,sound(2),ok)
    call require(ok,'LF right stationary contact')
    expected=-.25_dp*sum(sound)*(other_state-state)
    expected(axis+1)=expected(axis+1)+face(5)
    q(:,-1)=nan;q(:,2)=nan
    call mol_riemann_flux_3d(species,q,sound,axis,connected,1,flux,ok,lax_friedrichs=.true.)
    call require(ok.and.maxval(abs(flux-expected)/max(1._dp,abs(expected)))<2e-13_dp,'LF stationary contact diffusion')
    call require(abs(sum(flux(6:))-flux(1))<1e-12_dp,'LF species mass closure')
    ! Equal density and zero velocities isolate pressure-driven contact speed.
    face=q(:,0);face(5)=300000; q(:,1)=face
    call reactive_primitive_to_conserved(species,face,other_state,t,sound(2),ok)
    call require(ok,'LF pressure jump EOS')
    diffusion_speed=100000._dp/(face(1)*sum(sound))+.5_dp*sum(sound)
    expected=-.5_dp*diffusion_speed*(other_state-state)
    expected(axis+1)=250000
    call mol_riemann_flux_3d(species,q,sound,axis,connected,1,flux,ok,lax_friedrichs=.true.)
    call require(ok.and.maxval(abs(flux-expected)/max(1._dp,abs(expected)))<2e-13_dp,'LF pressure-driven speed')
    q(:,-1)=nan;q(:,2)=nan;q(:,0)=center;q(:,1)=center
    connected(:,1)=[.false.,.true.];connected(:,2)=[.true.,.false.]
    call mol_riemann_flux_3d(species,q,sound,axis,connected,2,flux,ok)
    call require(ok,'inactive outer donors ignored')
    connected=.true.;flux=19
    call mol_riemann_flux_3d(species,q,sound,axis,connected,2,flux,ok)
    call require(.not.ok.and.all(flux==19),'active NaN rollback')
    call mol_riemann_flux_3d(species,q,sound,axis,connected,2,flux,ok,lax_friedrichs=.true.)
    call require(.not.ok.and.all(flux==19),'LF active NaN rollback')
  end do
  print *, 'PASS MOL Riemann checks:',checks
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL ',label
      error stop 1
    end if
    checks=checks+1
  end subroutine
end program
