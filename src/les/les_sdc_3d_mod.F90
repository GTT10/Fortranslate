module les_sdc_3d_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_3d_mod, only: advance_reactive_sdc_3d,compute_reactive_external_source_3d
  use les_smagorinsky_3d_mod, only: compute_les_rhs_3d
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_boundary_rhs_3d
  implicit none
  private
  public :: advance_les_sdc_3d
contains
  ! LES is an old/new source in the SDC predictor and correction. Its surface
  ! inventory belongs to the boundary ledger, not the volumetric forcing ledger.
  subroutine advance_les_sdc_3d(species,state,temperature,spacing,dt,time,iterations,cs,prt,sct, &
      boundary_change,source_change,maximum_diffusivity,ok,model,ci,external_force,transport,transport_flags, &
      reactions,reaction_source,chemistry_relative_tolerance,chemistry_absolute_tolerance,difmag,boundary_faces, &
      inflow_state,inflow_temperature,inflow_mach,inflow_pressure,outlet_pressure,wall_temperature,wall_velocity, &
      test_filter,test_ratio,grid_ratio,harmonic_mean)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),dt,time,cs,prt,sct
    integer, intent(in) :: iterations
    real(dp), intent(inout) :: boundary_change(:),source_change(:),maximum_diffusivity
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: model,test_filter,boundary_faces(6)
    real(dp), optional, intent(in) :: ci,external_force(3),difmag
    type(gas_transport_species), optional, intent(in) :: transport(:)
    logical, optional, intent(in) :: transport_flags(4),harmonic_mean
    type(elementary_reaction), optional, intent(in) :: reactions(:)
    real(dp), optional, intent(inout) :: reaction_source(:,:,:,:)
    real(dp), optional, intent(in) :: chemistry_relative_tolerance,chemistry_absolute_tolerance
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: outlet_pressure(6),wall_temperature(:),wall_velocity(:,:)
    integer, optional, intent(in) :: test_ratio,grid_ratio
    real(dp), allocatable :: candidate(:,:,:,:),ct(:,:,:),ir(:,:,:,:)
    real(dp) :: boundary(size(state,1)),source(size(state,1)),old_rate(size(state,1)),new_rate(size(state,1))
    real(dp) :: force(3),isotropic,diffusion,les_boundary(size(state,1))
    integer :: calls,ratio
    character(len=32) :: closure,filter
    logical :: valid
    ok=.false.;calls=0;diffusion=0;old_rate=0;new_rate=0
    if(size(boundary_change)/=size(state,1).or.size(source_change)/=size(state,1)) return
    closure='smagorinsky';filter='box';ratio=2;isotropic=0;force=0
    if(present(model)) closure=model
    if(present(test_filter)) filter=test_filter
    if(present(test_ratio)) ratio=test_ratio
    if(present(ci)) isotropic=ci
    if(present(external_force)) force=external_force
    if(.not.all(ieee_is_finite([cs,prt,sct,isotropic,force]))) return
    if(cs<0.or.cs>2.or.min(prt,sct)<1e-6_dp.or.max(prt,sct)>1e6_dp) return
    if(isotropic<0.or.isotropic>1) return
    select case(closure)
    case('smagorinsky','wale','vreman')
      if(present(harmonic_mean)) then
        if(harmonic_mean) return
      end if
    case('dynamic')
      if(isotropic/=0) return
    case default
      return
    end select
    candidate=state;ct=temperature
    if(present(reaction_source)) then
      if(.not.present(reactions)) return
      allocate(ir,source=reaction_source)
    end if
    call advance_reactive_sdc_3d(species,candidate,ct,spacing,dt,time,iterations, &
      2+merge(1,0,present(transport)),evaluate,valid,difmag=difmag, &
      source_change=source,transport=transport,transport_flags=transport_flags,reactions=reactions,reaction_source=ir, &
      chemistry_relative_tolerance=chemistry_relative_tolerance,chemistry_absolute_tolerance=chemistry_absolute_tolerance, &
      boundary_faces=boundary_faces,boundary_change=boundary,outlet_pressure=outlet_pressure, &
      inflow_state=inflow_state,inflow_temperature=inflow_temperature,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure, &
      wall_temperature=wall_temperature,wall_velocity=wall_velocity)
    if(.not.valid) return
    les_boundary=.5_dp*dt*(old_rate+new_rate)
    boundary=boundary+les_boundary;source=source-les_boundary
    if(.not.all(ieee_is_finite([boundary,source,diffusion]))) return
    state=candidate;temperature=ct
    if(present(reaction_source)) reaction_source=ir
    boundary_change=boundary;source_change=source;maximum_diffusivity=diffusion;ok=.true.
  contains
    subroutine evaluate(u,t,stage_time,rates,good)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),stage_time
      real(dp), intent(out) :: rates(:,:,:,:,:)
      logical, intent(out) :: good
      real(dp) :: rate(size(state,1)),viscosity
      rates=0;good=.false.;rate=0;viscosity=0
      if(.not.ieee_is_finite(stage_time).or.size(rates,5)/=2) return
      call compute_reactive_external_source_3d(u,force,rates(:,:,:,:,1),good)
      if(.not.good) return
      if(closure=='dynamic') then
        call compute_dynamic_les_boundary_rhs_3d(species,u,t,spacing,filter,ratio,rates(:,:,:,:,2),viscosity,good, &
          grid_ratio=grid_ratio,harmonic_mean=harmonic_mean,boundary_faces=boundary_faces,boundary_rate=rate, &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_temperature=wall_temperature, &
          wall_velocity=wall_velocity,outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      else
        call compute_les_rhs_3d(species,u,t,spacing,cs,prt,sct,rates(:,:,:,:,2),viscosity,good, &
          boundary_faces=boundary_faces,boundary_rate=rate,model=closure,ci=isotropic, &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_temperature=wall_temperature, &
          wall_velocity=wall_velocity,outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      end if
      if(.not.good) return
      calls=calls+1
      if(calls==1) old_rate=rate
      new_rate=rate;diffusion=max(diffusion,viscosity)
    end subroutine evaluate
  end subroutine advance_les_sdc_3d
end module les_sdc_3d_mod
