module spray_eb_reactive_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_eb_step_mod, only: advance_eb_mixture_coupling
  use reactive_eb_3d_driver_mod, only: advance_reactive_eb_full_3d
  implicit none
  private
  public :: advance_eb_spray_reactive
contains
  ! Particle half / existing reactive EB gas step / particle half.
  ! Inventories use the change in the gas-plus-particle system sign convention.
  subroutine advance_eb_spray_reactive(species,reactions,transport,liquids,options,g,phi,state,temperature, &
      parcels,component_mass,dt,solver,redistribution,chemistry,transport_enabled,viscosity,conduction, &
      diffusion,barodiffusion,rtol,atol,particle_steps,minimum_theta,boundary_change,wall_change,chemistry_change, &
      ok,target_volume_fraction,chemistry_integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(eb_geometry_3d), intent(in) :: g
    real(dp), intent(in) :: phi(0:,0:,0:),dt,rtol,atol
    character(len=*), intent(in) :: solver,redistribution
    logical, intent(in) :: chemistry,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),component_mass(:,:)
    type(spray_parcel), intent(inout) :: parcels(:)
    real(dp), intent(inout) :: minimum_theta,boundary_change(:),wall_change(:),chemistry_change(:)
    integer, intent(out) :: particle_steps
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: target_volume_fraction
    character(len=*), optional, intent(in) :: chemistry_integrator
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),next_gas(:,:,:,:),next_temp(:,:,:),masses(:,:)
    type(spray_parcel), allocatable :: drops(:)
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    real(dp) :: chemical(reactive_nvar(size(species))),first_wall(5),second_wall(5),theta
    integer :: first_steps,second_steps,nv
    logical :: valid
    ok=.false.;particle_steps=0;nv=reactive_nvar(size(species))
    if (.not.ieee_is_finite(dt)) return
    if (dt<=0.or.dt>1e6_dp) return
    if (size(boundary_change)/=nv.or.size(wall_change)/=nv.or.size(chemistry_change)/=nv) return
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    first_wall=0;second_wall=0;boundary=0;wall=0;chemical=0
    call advance_eb_mixture_coupling(species,transport,liquids,options,g,phi,gas,temp,drops,masses, &
      .5_dp*dt,first_steps,first_wall,valid)
    if (.not.valid) return
    next_gas=gas;next_temp=temp
    call advance_reactive_eb_full_3d(species,reactions,transport,gas,temp,g,solver,redistribution,dt, &
      chemistry,rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion,next_gas,next_temp, &
      theta,valid,target_volume_fraction,chemistry_integrator,boundary,wall,chemical)
    if (.not.valid) return
    call advance_eb_mixture_coupling(species,transport,liquids,options,g,phi,next_gas,next_temp,drops,masses, &
      .5_dp*dt,second_steps,second_wall,valid)
    if (.not.valid) return
    wall(1:5)=wall(1:5)+first_wall+second_wall
    state=next_gas;temperature=next_temp;parcels=drops;component_mass=masses
    minimum_theta=theta;boundary_change=boundary;wall_change=wall;chemistry_change=chemical
    particle_steps=first_steps+second_steps;ok=.true.
  end subroutine
end module
