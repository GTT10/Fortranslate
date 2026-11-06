program test_eb_wall_heat_source_3d
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use mixture_transport_mod, only: mixture_transport_coefficients
  use mixture_thermo_mod, only: mixture_mass_properties
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_mass_fraction_component,reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_wall_heat_source_3d_mod, only: compute_eb_wall_heat_source_3d,compute_eb_wall_transport_source_3d
  use eb_reactive_transport_3d_mod, only: compute_sdc_eb_transport_3d,advance_reactive_eb_transport_3d, &
    reactive_eb_transport_patch_timestep_3d,reactive_eb_isothermal_timestep_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use state_indices_mod, only: iet
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(eb_geometry_3d) :: g
  type(eb_integrated_flux_3d) :: face_flux
  real(dp), allocatable :: state(:,:,:,:),temp(:,:,:),q(:),y(:),diff(:)
  real(dp) :: rates(6,6,6),saved(6,6,6),power,bound,mu,conductivity,mw,rgas,cp,cv,gamma,enthalpy,energy,entropy
  real(dp) :: distance,expected,inventory,old_power,old_bound,recovered_t,sound
  real(dp), allocatable :: baseline(:,:,:,:),stage(:,:,:,:)
  real(dp) :: theta,wall_energy,stage_frequency,coefficients(4),dt,diffusivity,wall_dt
  real(dp), allocatable :: advanced(:,:,:,:),advanced_temperature(:,:,:),moving(:,:,:,:)
  real(dp), allocatable :: exchange(:)
  real(dp) :: constant_power,constant_bound
  character(len=20) :: physical_faces(6)
  integer :: i,j,k,s
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  if(.not.ok) error stop 'thermo load'
  call load_h2o2_elementary_transport(transport,ok)
  if(.not.ok) error stop 'transport load'
  allocate(state(reactive_nvar(size(species)),6,6,6),temp(6,6,6),q(reactive_nprim(size(species))), &
    y(size(species)),diff(size(species)))
  y=0
  do s=1,size(species)
    if(species(s)%name=='N2') y(s)=1
  end do
  if(sum(y)/=1) error stop 'nitrogen fixture'
  call build_axis_plane_eb_geometry_3d(6,6,6,0._dp,.06_dp,0._dp,.06_dp,0._dp,.06_dp,'x',.023_dp,g,ok)
  if(.not.ok) error stop 'geometry'
  do i=1,6
    distance=(i-.5_dp)*.01_dp-.023_dp
    temp(i,:,:)=400+1000*distance+1e4_dp*distance**2
    call mixture_mass_properties(species,y,temp(i,1,1),mw,rgas,cp,cv,gamma,enthalpy,energy,entropy,ok)
    if(.not.ok) error stop 'properties'
    q=0;q(1)=1;q(5)=rgas*temp(i,1,1)
    do s=1,size(species)
      q(reactive_mass_fraction_component(s))=y(s)
    end do
    do k=1,6
      do j=1,6
        call reactive_primitive_to_conserved(species,q,state(:,i,j,k),recovered_t,sound,ok)
        if(.not.ok) error stop 'state initialization'
      end do
    end do
  end do
  call mixture_transport_coefficients(species,transport,y,temp(3,1,1),rgas*temp(3,1,1),mu,conductivity,diff,ok)
  if(.not.ok) error stop 'coefficient'
  rates=123;power=456;bound=789
  call compute_eb_wall_heat_source_3d(species,transport,state,temp,g,[1,1,1],[6,6,6],400._dp,rates,power,bound,ok)
  if(.not.ok) error stop 'wall source'
  expected=-conductivity*1000*.06_dp**2
  if(abs(power-expected)>1e-11_dp*abs(expected)) error stop 'analytic wall heat'
  inventory=sum(rates*g%volume_fraction)*g%dx*g%dy*g%dz
  if(abs(inventory-power)>1e-13_dp*abs(power)) error stop 'source wall ledger balance'
  if(any(rates(1:2,:,:)/=0).or.any(rates(4:6,:,:)/=0).or.bound<=0) error stop 'source support/frequency'
  allocate(baseline(size(state,1),6,6,6),stage(size(state,1),6,6,6))
  baseline=0;stage=0;theta=1
  call compute_sdc_eb_transport_3d(species,transport,state,temp,g,1e-7_dp,[.false.,.true.,.false.,.false.], &
    baseline,face_flux,theta,ok)
  if(.not.ok) error stop 'baseline SDC transport'
  call compute_sdc_eb_transport_3d(species,transport,state,temp,g,1e-7_dp,[.false.,.true.,.false.,.false.], &
    stage,face_flux,theta,ok,embedded_wall_temperature=400._dp,embedded_wall_energy=wall_energy, &
    embedded_frequency=stage_frequency)
  if(.not.ok) error stop 'isothermal SDC transport'
  stage(iet,:,:,:)=stage(iet,:,:,:)-rates
  if(maxval(abs(stage-baseline))>1e-9_dp) error stop 'SDC source routing'
  if(abs(wall_energy-power*1e-7_dp)>1e-20_dp.or.stage_frequency/=bound) error stop 'SDC wall ledger'
  ! Full constants bypass the molecular database at faces and embedded walls.
  coefficients=[1.e-5_dp,.02_dp,3.e-5_dp,2.e-5_dp]
  transport%diameter=0
  allocate(exchange(size(state,1)))
  call compute_eb_wall_heat_source_3d(species,transport,state,temp,g,[1,1,1],[6,6,6],400._dp, &
    rates,power,bound,ok,constant_coefficients=coefficients)
  if(.not.ok) error stop 'constant embedded heat bypass'
  expected=-coefficients(2)*1000*.06_dp**2
  if(abs(power-expected)>1e-11_dp*abs(expected)) error stop 'constant analytic heat'
  constant_power=power;constant_bound=bound
  call compute_eb_wall_transport_source_3d(species,transport,state,temp,g,[1,1,1],[6,6,6], &
    .true.,.true.,400._dp,stage,exchange,stage_frequency,ok,constant_coefficients=coefficients)
  if(.not.ok.or.any(exchange(2:4)/=0).or.exchange(5)/=power) error stop 'constant no-slip wall'
  moving=state
  do i=1,6
    q=0;q(1)=1;q(5)=rgas*temp(i,1,1);q(2)=2*((i-.5_dp)*.01_dp-.023_dp)
    do s=1,size(species)
      q(reactive_mass_fraction_component(s))=y(s)
    enddo
    do k=1,6
      do j=1,6
        call reactive_primitive_to_conserved(species,q,moving(:,i,j,k),recovered_t,sound,ok)
        if(.not.ok) error stop 'moving wall fixture'
      enddo
    enddo
  enddo
  call compute_eb_wall_transport_source_3d(species,transport,moving,temp,g,[1,1,1],[6,6,6], &
    .false.,.true.,400._dp,stage,exchange,stage_frequency,ok,constant_coefficients=coefficients)
  if(.not.ok) error stop 'constant viscous wall bypass'
  expected=(4._dp/3*coefficients(1)+coefficients(3))*2*.06_dp**2
  if(abs(abs(exchange(2))-expected)>1e-11_dp*expected.or.any(exchange(3:5)/=0)) &
    error stop 'analytic normal viscous wall traction'
  inventory=sum(stage(2,:,:,:)*g%volume_fraction)*g%dx*g%dy*g%dz
  if(abs(inventory-exchange(2))>1e-13_dp*expected) error stop 'wall force volume ledger'
  old_bound=stage_frequency;old_power=exchange(2)
  coefficients(1)=2*coefficients(1);coefficients(3)=2*coefficients(3)
  call compute_eb_wall_transport_source_3d(species,transport,moving,temp,g,[1,1,1],[6,6,6], &
    .false.,.true.,400._dp,stage,exchange,stage_frequency,ok,constant_coefficients=coefficients)
  if(.not.ok.or.abs(exchange(2)-2*old_power)>1e-13_dp*expected.or. &
    abs(stage_frequency-2*old_bound)>1e-13_dp*old_bound) error stop 'bulk wall stress/CFL scaling'
  coefficients(1)=1.e-5_dp;coefficients(3)=3.e-5_dp
  coefficients(2)=2*coefficients(2)
  call compute_eb_wall_heat_source_3d(species,transport,state,temp,g,[1,1,1],[6,6,6],400._dp, &
    rates,power,bound,ok,constant_coefficients=coefficients)
  if(.not.ok.or.abs(power-2*constant_power)>1e-13_dp.or. &
    abs(bound-2*constant_bound)>1e-11_dp*constant_bound) error stop 'constant wall scaling'
  coefficients(2)=.02_dp
  call compute_eb_wall_heat_source_3d(species,transport,state,temp,g,[1,1,1],[6,6,6],400._dp, &
    rates,power,bound,ok,constant_coefficients=coefficients)
  if(.not.ok) error stop 'restore constant wall source'
  call compute_sdc_eb_transport_3d(species,transport,state,temp,g,1e-7_dp,[.true.,.true.,.true.,.true.], &
    baseline,face_flux,theta,ok,constant_coefficients=coefficients)
  if(.not.ok) error stop 'constant baseline SDC'
  call compute_sdc_eb_transport_3d(species,transport,state,temp,g,1e-7_dp,[.true.,.true.,.true.,.true.], &
    stage,face_flux,theta,ok,embedded_wall_temperature=400._dp,embedded_wall_energy=wall_energy, &
    embedded_frequency=stage_frequency,constant_coefficients=coefficients)
  if(.not.ok) error stop 'constant embedded SDC'
  stage(iet,:,:,:)=stage(iet,:,:,:)-rates
  if(maxval(abs(stage-baseline))>1e-9_dp) error stop 'constant SDC wall routing'
  if(abs(wall_energy-power*1e-7_dp)>1e-20_dp.or.stage_frequency/=bound) error stop 'constant SDC wall ledger'
  call reactive_eb_transport_patch_timestep_3d(species,transport,state,temp,g,.2_dp, &
    .true.,.true.,.true.,dt,diffusivity,ok,constant_coefficients=coefficients)
  if(.not.ok.or.dt<=0.or.diffusivity<=0) error stop 'constant cut cell CFL'
  physical_faces='outflow';physical_faces(4)='isothermal_wall'
  call reactive_eb_isothermal_timestep_3d(species,transport,state,temp,g,physical_faces, &
    spread(400._dp,1,6),.2_dp,wall_dt,ok,constant_coefficients=coefficients)
  if(.not.ok.or.wall_dt<=0) error stop 'constant physical wall CFL'
  advanced=state;advanced_temperature=temp
  call advance_reactive_eb_transport_3d(species,transport,advanced,advanced_temperature,g,1e-9_dp, &
    .true.,.true.,.true.,.true.,.5_dp,theta,ok,constant_coefficients=coefficients)
  if(.not.ok) error stop 'constant SSPRK2 forwarding'
  saved=rates;old_power=power;old_bound=bound
  call compute_eb_wall_heat_source_3d(species,transport,state,temp,g,[1,1,1],[6,6,6],400._dp, &
    rates,power,bound,ok,constant_coefficients=[1._dp,2._dp,3._dp])
  if(ok.or.any(rates/=saved).or.power/=old_power.or.bound/=old_bound) error stop 'invalid constants rollback'
  stage=123;theta=321
  call compute_sdc_eb_transport_3d(species,transport,state,temp,g,1e-7_dp,[.false.,.true.,.false.,.false.], &
    stage,face_flux,theta,ok,embedded_wall_temperature=400._dp)
  if(ok.or.any(stage/=123).or.theta/=321) error stop 'SDC missing ledger rollback'
  call check_held_wall()
  saved=rates;old_power=power;old_bound=bound
  ! Active stencil donors include i=4; failure must not publish partial rows.
  state(1,4,6,6)=-1
  call compute_eb_wall_heat_source_3d(species,transport,state,temp,g,[1,1,1],[6,6,6],400._dp,rates,power,bound,ok, &
    constant_coefficients=coefficients)
  if(ok.or.any(rates/=saved).or.power/=old_power.or.bound/=old_bound) error stop 'late donor failure rollback'
  print *, 'PASS: analytic heat, fluid-volume ledger, SDC routing and rollback',expected,inventory
contains
  subroutine check_held_wall()
    real(dp) :: held(6,6,6),thermal(6,6,6),trial(size(state,1),6,6,6),out(size(state,1),6,6,6)
    real(dp) :: exchange_held(size(state,1)),frequency_held,pressure_work,viscous_work
    real(dp) :: saved_out(size(state,1),6,6,6),saved_exchange(size(state,1)),saved_frequency
    integer :: x
    do x=1,6
      call mixture_mass_properties(species,y,temp(x,1,1),mw,rgas,cp,cv,gamma,enthalpy,energy,entropy,ok)
      if(.not.ok) error stop 'held wall NASA fixture'
      held(x,:,:)=energy
    enddo
    ! Independently held thermal energy remains resolved when total energy
    ! cannot support EOS recovery. Covered donors do not contribute.
    held(1:2,:,:)=ieee_value(0._dp,ieee_quiet_nan)
    thermal=temp;trial=moving;trial(iet,:,:,:)=-1.e20_dp
    out=123;exchange_held=456;frequency_held=789
    call compute_eb_wall_transport_source_3d(species,transport,trial,thermal,g,[1,1,1],[6,6,6], &
      .true.,.true.,400._dp,out,exchange_held,frequency_held,ok,constant_coefficients=coefficients, &
      internal_energy_density=held)
    if(.not.ok) error stop 'held quadratic thermal/viscous wall EOS'
    pressure_work=-coefficients(2)*1000*.06_dp**2
    viscous_work=(4._dp/3*coefficients(1)+coefficients(3))*2*.06_dp**2
    if(abs(exchange_held(5)-pressure_work)>1.e-11_dp*abs(pressure_work)) &
      error stop 'held independent quadratic heat flux'
    if(abs(abs(exchange_held(2))-viscous_work)>1.e-11_dp*viscous_work.or. &
      any(exchange_held(3:4)/=0).or.frequency_held<=0) error stop 'held independent wall traction'
    if(abs(sum(out(5,:,:,:)*g%volume_fraction)*g%dx*g%dy*g%dz-exchange_held(5))> &
      1.e-13_dp*abs(pressure_work)) error stop 'held wall thermal ledger'
    if(abs(sum(out(2,:,:,:)*g%volume_fraction)*g%dx*g%dy*g%dz-exchange_held(2))> &
      1.e-13_dp*viscous_work) error stop 'held wall momentum ledger'
    saved_out=out;saved_exchange=exchange_held;saved_frequency=frequency_held
    held(4,6,6)=1.e30_dp
    call compute_eb_wall_transport_source_3d(species,transport,trial,thermal,g,[1,1,1],[6,6,6], &
      .true.,.true.,400._dp,out,exchange_held,frequency_held,ok,constant_coefficients=coefficients, &
      internal_energy_density=held)
    if(ok.or.any(out/=saved_out).or.any(exchange_held/=saved_exchange).or.frequency_held/=saved_frequency) &
      error stop 'late finite held wall donor rollback'
  end subroutine check_held_wall
end program
