program test_sdc_cvode_bridge
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use fixture_mechanism_mod, only: load_fixture_thermo_data,load_fixture_mechanism
  use mixture_thermo_mod, only: mixture_mass_properties
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: advance_reactive_sdc_3d
  use sundials_forced_reactor_mod, only: advance_forced_reactor_cvode
  implicit none
  integer, parameter :: n=4,nv=7
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: state(nv,n,n,n),saved(nv,n,n,n),temperature(n,n,n),ir(3,n,n,n)
  real(dp) :: primitive(nv),cell(nv),y(2),rho_y(2),reference_t,mw,r,cp,cv,gamma,h,energy,s,delta
  real(dp), parameter :: dt=1e-10_dp
  integer :: i,j,k,calls,steps
  logical :: ok,fail
  call load_fixture_thermo_data(species,ok)
  call require(ok,'thermo')
  call load_fixture_mechanism(reactions,ok)
  call require(ok,'reactions')
  y=[.8_dp,.2_dp]
  call mixture_mass_properties(species,y,1000._dp,mw,r,cp,cv,gamma,h,energy,s,ok)
  call require(ok,'initial energy')
  primitive=[1._dp,0._dp,0._dp,0._dp,r*1000._dp,y]
  call reactive_primitive_to_conserved(species,primitive,cell,reference_t,delta,ok)
  call require(ok,'conserved state')
  do k=1,n
  do j=1,n
  do i=1,n
    state(:,i,j,k)=cell
  enddo
  enddo
  enddo
  saved=state;temperature=1000;ir=0;calls=0;fail=.false.
  rho_y=y;reference_t=1000
  call advance_forced_reactor_cvode(species,reactions,energy,[.1_dp,-.02_dp],1e6_dp,dt, &
    1e-10_dp,1e-13_dp,rho_y,reference_t,steps,ok)
  call require(ok,'single-cell reference')
  call run_step()
  call require(ok.and.calls==2*n**3,'CVODE invoked for every cell in both SDC iterations')
  delta=0
  do k=1,n
  do j=1,n
  do i=1,n
    delta=max(delta,maxval(abs(state(6:,i,j,k)-rho_y)))
    call require(abs(temperature(i,j,k)-reference_t)<2e-7_dp,'SDC temperature matches forced reactor')
    call require(abs(state(1,i,j,k)-1-dt*.08_dp)<1e-12_dp,'SDC external mass inventory')
    call require(abs(state(5,i,j,k)-energy-dt*1e6_dp)<1e-8_dp,'SDC external energy inventory')
    call require(maxval(abs(dt*ir(:2,i,j,k)-(state(6:,i,j,k)-y-dt*[.1_dp,-.02_dp])))<1e-15_dp, &
      'reaction-source history excludes external source')
  enddo
  enddo
  enddo
  call require(delta<1e-10_dp,'SDC species match forced reactor')
  call require(maxval(abs(state(6:,:,:,:)-saved(6:,:,:,:)))>1e-9_dp,'chemistry changes gas')
  print *, 'Maximum SDC/single-cell species difference:',delta
  state=saved;temperature=1000;ir=0;calls=0;fail=.true.
  call run_step()
  call require(.not.ok.and.calls==n**3+1,'second-iteration reactor rejection')
  call require(all(state==saved).and.all(temperature==1000).and.all(ir==0),'SDC transaction rollback')
contains
  subroutine run_step()
    call advance_reactive_sdc_3d(species,state,temperature,[.01_dp,.01_dp,.01_dp],dt,0._dp,2,1, &
      sources,ok,reactions=reactions,reaction_source=ir,chemistry_relative_tolerance=1e-10_dp, &
      chemistry_absolute_tolerance=1e-13_dp,evaluate_reaction=reactor)
  end subroutine
  subroutine sources(u,t,time,rates,good)
    real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),time
    real(dp), intent(out) :: rates(:,:,:,:,:)
    logical, intent(out) :: good
    rates=0;rates(1,:,:,:,1)=.08_dp;rates(5,:,:,:,1)=1e6_dp
    rates(6,:,:,:,1)=.1_dp;rates(7,:,:,:,1)=-.02_dp;good=.true.
  end subroutine
  subroutine reactor(sp,rx,e0,external_y,external_e,duration,rtol,atol,densities,temp,count,good,maximum_steps)
    type(nasa7_species), intent(in) :: sp(:)
    type(elementary_reaction), intent(in) :: rx(:)
    real(dp), intent(in) :: e0,external_y(:),external_e,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),temp
    integer, intent(out) :: count
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    calls=calls+1
    if(fail.and.calls==n**3+1) then
      densities=-42;temp=-7;count=0;good=.false.;return
    endif
    call advance_forced_reactor_cvode(sp,rx,e0,external_y,external_e,duration,rtol,atol, &
      densities,temp,count,good,maximum_steps)
  end subroutine
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition) then
      print *, 'FAIL: ',message
      error stop 1
    endif
  end subroutine
end program
