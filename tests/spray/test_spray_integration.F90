program test_spray_integration
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use mixture_thermo_mod, only: mixture_specific_gas_constant, mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use spray_parcel_mod
  use spray_coupling_3d_mod, only: spray_system_totals
  use spray_reactive_3d_mod
  use spray_checkpoint_mod
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_parcel), allocatable :: drops(:),before_drops(:)
  type(spray_liquid) :: liquid
  type(spray_options) :: options
  real(dp) :: state(15,5,5,5),temperature(5,5,5),before(15,5,5,5),before_t(5,5,5)
  real(dp) :: x(10),y(10),primitive(15),r,sound,local_t,totals0(15),totals1(15),injected(15),time
  integer :: i,j,k,iteration,steps,event
  logical :: ok
  character(len=32) :: which,method
  character(len=256) :: message
  call get_command_argument(1,which)
  call load_h2o2_full_thermo(species,ok); call require(ok,'thermo')
  call load_h2o2_full_mechanism(reactions,ok); call require(ok,'reactions')
  call load_h2o2_full_transport(transport,ok); call require(ok,'transport')
  x=0.0_dp; x([1,4,10])=[2.0_dp,1.0_dp,3.76_dp]; x=x/sum(x)
  call mass_fractions_from_mole_fractions(species,x,y,ok); call require(ok,'composition')
  r=mixture_specific_gas_constant(species,y,ok); call require(ok,'R')
  primitive=0.0_dp; primitive(1)=101325.0_dp/(r*1400.0_dp); primitive(5)=101325.0_dp; primitive(6:)=y
  do k=1,5
    do j=1,5
      do i=1,5
        call reactive_primitive_to_conserved(species,primitive,state(:,i,j,k),local_t,sound,ok)
        call require(ok,'initial state')
      end do
    end do
  end do
  temperature=1400.0_dp; before=state; before_t=temperature
  call initialize_spray_liquid(species,'H2O',1000.0_dp,4184.0_dp,373.15_dp,101325.0_dp,2.257e6_dp,liquid,ok)
  call require(ok,'liquid')
  call inject_cone_parcels(liquid,1.0e-10_dp,2,1_int64,[0.002_dp,0.002_dp,0.002_dp], &
    [0.0_dp,0.0_dp,1.0_dp],10.0_dp,0.1_dp,2.0e-5_dp,330.0_dp,drops,ok)
  call require(ok,'injection'); before_drops=drops
  call spray_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,liquid,totals0)
  select case(trim(which))
  case('native','cvode')
    method='implicit'
    if (trim(which)=='cvode') method='cvode'
    ! 125 cells exceed the backend context-slot table. Three successive calls
    ! expose contexts not finalized after each independently changing cell.
    do iteration=1,3
      call advance_spray_chemistry_3d(species,reactions,state,temperature,5.0e-7_dp, &
        1.0e-7_dp,1.0e-12_dp,method,ok,message)
      call require(ok,trim(message))
    end do
    call require(all(state(1:5,:,:,:)==before(1:5,:,:,:)),'chemistry preserves conserved hydro variables')
    call require(maxval(abs(state(6:,:,:,:)-before(6:,:,:,:)))>1.0e-10_dp,'chemistry actually advances')
    call require(all(state(6:,:,:,:)>=0.0_dp),'nonnegative chemistry')
    call require(maxval(abs(sum(state(6:,:,:,:),dim=1)-state(1,:,:,:)))<1.0e-12_dp,'species closure')
    call require(maxval(abs(state(15,:,:,:)-before(15,:,:,:)))<1.0e-12_dp,'inert nitrogen unchanged')
  case('coupled')
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2.0e-7_dp,.true.,.true.,'implicit', &
      1.0e-7_dp,1.0e-12_dp,0.1_dp,0.7_dp,0.7_dp,ok,message)
    call require(ok,trim(message))
    call spray_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,liquid,totals1)
    call require(abs(totals1(1)-totals0(1))<1.0e-12_dp*totals0(1),'coupled mass balance')
    call require(maxval(abs(totals1(2:4)-totals0(2:4)))<1.0e-18_dp,'coupled momentum balance')
    call require(abs(totals1(5)-totals0(5))<1.0e-12_dp*abs(totals0(5)),'coupled energy balance')
    call require(sum(drops%mass)<sum(before_drops%mass),'coupled evaporation active')
    call require(maxval(abs(state(6,:,:,:)-before(6,:,:,:)))>1.0e-12_dp,'coupled chemistry active')
  case('rollback')
    ! The first spray/LES stages run before the invalid backend is encountered.
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2.0e-7_dp,.true.,.true.,'invalid', &
      1.0e-7_dp,1.0e-12_dp,0.1_dp,0.7_dp,0.7_dp,ok,message)
    call require(.not.ok,'late chemistry rejection')
    call require(all(state==before) .and. all(temperature==before_t),'gas whole-step rollback')
    do i=1,size(drops)
      call require(all(parcel_totals(drops(i),liquid)==parcel_totals(before_drops(i),liquid)),'parcel rollback')
      call require(all(drops(i)%position==before_drops(i)%position),'position rollback')
    end do
  case('checkpoint')
    injected=0.0_dp; time=1.0e-6_dp; steps=3; event=0
    call write_spray_checkpoint('roundtrip.chk','integration-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message)
    call require(ok,trim(message))
    call write_spray_checkpoint('roundtrip.chk','integration-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message)
    call require(.not.ok,'existing checkpoint not overwritten')
    state=0.0_dp; temperature=0.0_dp; time=0.0_dp; steps=0; event=0; totals1=0.0_dp
    deallocate(drops); allocate(drops(0))
    call read_spray_checkpoint('roundtrip.chk','wrong-context',species,state,temperature,drops,time,steps,event, &
      totals1,injected,ok,message)
    call require(.not.ok .and. all(state==0.0_dp) .and. size(drops)==0,'invalid restart transaction')
    call read_spray_checkpoint('roundtrip.chk','integration-context',species,state,temperature,drops,time,steps,event, &
      totals1,injected,ok,message)
    call require(ok,trim(message))
    call require(all(state==before) .and. all(temperature==before_t),'restart exact state')
    call require(all(totals0==totals1) .and. size(drops)==2 .and. steps==3,'restart metadata')
    open(newunit=i,file='roundtrip.chk',status='old'); close(i,status='delete')
  case default
    error stop 'unknown integration test'
  end select
  print *,trim(which),' PASS'
contains
  subroutine require(condition,text)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: text
    if (.not.condition) then
      print *,trim(text)
      error stop 1
    end if
  end subroutine
end program test_spray_integration
