program test_eb_wall_heat_source_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use mixture_transport_mod, only: mixture_transport_coefficients
  use mixture_thermo_mod, only: mixture_mass_properties
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_mass_fraction_component,reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_wall_heat_source_3d_mod, only: compute_eb_wall_heat_source_3d
  use eb_reactive_transport_3d_mod, only: compute_sdc_eb_transport_3d
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
  real(dp) :: theta,wall_energy,stage_frequency
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
  stage=123;theta=321
  call compute_sdc_eb_transport_3d(species,transport,state,temp,g,1e-7_dp,[.false.,.true.,.false.,.false.], &
    stage,face_flux,theta,ok,embedded_wall_temperature=400._dp)
  if(ok.or.any(stage/=123).or.theta/=321) error stop 'SDC missing ledger rollback'
  saved=rates;old_power=power;old_bound=bound
  ! Active stencil donors include i=4; failure must not publish partial rows.
  state(1,4,6,6)=-1
  call compute_eb_wall_heat_source_3d(species,transport,state,temp,g,[1,1,1],[6,6,6],400._dp,rates,power,bound,ok)
  if(ok.or.any(rates/=saved).or.power/=old_power.or.bound/=old_bound) error stop 'late donor failure rollback'
  print *, 'PASS: analytic heat, fluid-volume ledger, SDC routing and rollback',expected,inventory
end program
