module spray_eb_reactive_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options,valid_spray_parcel
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_breakup_mod, only: split_spray_parcels
  use spray_eb_step_mod, only: advance_eb_mixture_coupling_allocated
  use spray_eb_impact_mod, only: spray_eb_wall_options
  use spray_wall_impact_mod, only: spray_wall_film
  use, intrinsic :: iso_fortran_env, only: int64
  use reactive_eb_3d_driver_mod, only: advance_reactive_eb_full_3d
  implicit none
  private
  public :: advance_eb_spray_reactive,advance_eb_spray_reactive_allocated
contains
  subroutine advance_eb_spray_reactive(species,reactions,transport,liquids,options,g,phi,state,temperature, &
      parcels,component_mass,dt,solver,redistribution,chemistry,transport_enabled,viscosity,conduction, &
      diffusion,barodiffusion,rtol,atol,particle_steps,minimum_theta,boundary_change,wall_change,chemistry_change, &
      ok,target_volume_fraction,chemistry_integrator,acceleration,source_change)
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
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(inout) :: source_change(:)
    real(dp), optional, intent(in) :: target_volume_fraction
    character(len=*), optional, intent(in) :: chemistry_integrator
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    ok=.false.;particle_steps=0
    if (options%maximum_multiplicity>0) return
    drops=parcels;masses=component_mass
    call advance_eb_spray_reactive_allocated(species,reactions,transport,liquids,options,g,phi,state,temperature, &
      drops,masses,dt,solver,redistribution,chemistry,transport_enabled,viscosity,conduction, &
      diffusion,barodiffusion,rtol,atol,particle_steps,minimum_theta,boundary_change,wall_change,chemistry_change, &
      ok,target_volume_fraction,chemistry_integrator,acceleration=acceleration,source_change=source_change)
    if (.not.ok) return
    parcels=drops;component_mass=masses
  end subroutine
  ! Particle half / existing reactive EB gas step / particle half.
  ! Inventories use the change in the gas-plus-particle system sign convention.
  subroutine advance_eb_spray_reactive_allocated(species,reactions,transport,liquids,options,g,phi,state,temperature, &
      parcels,component_mass,dt,solver,redistribution,chemistry,transport_enabled,viscosity,conduction, &
      diffusion,barodiffusion,rtol,atol,particle_steps,minimum_theta,boundary_change,wall_change,chemistry_change, &
      ok,target_volume_fraction,chemistry_integrator,wall_options,films,film_mass,last_id,model_rejected, &
      split_created,particle_outflow,acceleration,source_change)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(eb_geometry_3d), intent(in) :: g
    real(dp), intent(in) :: phi(0:,0:,0:),dt,rtol,atol
    character(len=*), intent(in) :: solver,redistribution
    logical, intent(in) :: chemistry,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    type(spray_eb_wall_options), optional, intent(in) :: wall_options
    type(spray_wall_film), allocatable, optional, intent(inout) :: films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    integer(int64), optional, intent(inout) :: last_id
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    integer(int64), optional, intent(out) :: split_created
    integer(int64) :: watermark,before_split
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(inout) :: source_change(:)
    real(dp) :: accel(3),force_ledger(reactive_nvar(size(species)))
    logical, optional, intent(in) :: particle_outflow
    logical :: impacts,rejected,outflow
    character(len=16) :: outer_faces(6)
    real(dp) :: first_exterior(reactive_nvar(size(species))),second_exterior(reactive_nvar(size(species)))
    logical, optional, intent(out) :: model_rejected
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
    if (present(model_rejected)) model_rejected=.false.
    if (present(split_created)) split_created=0
    ok=.false.;particle_steps=0;nv=reactive_nvar(size(species))
    if (.not.ieee_is_finite(dt)) return
    if (dt<=0.or.dt>1e6_dp) return
    if (size(boundary_change)/=nv.or.size(wall_change)/=nv.or.size(chemistry_change)/=nv) return
    if (.not.allocated(parcels).or..not.allocated(component_mass)) return
    accel=0;force_ledger=0
    if(present(acceleration)) accel=acceleration
    if(.not.all(ieee_is_finite(accel))) return
    if(any(abs(accel)>1e6_dp)) return
    if(any(accel/=0).and..not.present(source_change)) return
    if(present(source_change)) then
      if(size(source_change)/=nv) return
    end if
    impacts=present(wall_options);outflow=.false.
    if(present(particle_outflow)) outflow=particle_outflow
    outer_faces='outflow';first_exterior=0;second_exterior=0
    if (present(films).neqv.impacts) return
    if (present(film_mass).neqv.impacts) return
    if (impacts.and..not.present(last_id)) return
    if ((options%maximum_multiplicity>0.or.options%khrt_children).and..not.present(last_id)) return
    watermark=0
    if(outflow) watermark=max(0_int64,maxval(parcels%id))
    if (present(last_id)) watermark=last_id
    if (impacts) then
      if (.not.allocated(films).or..not.allocated(film_mass)) return
      patches=films;patch_mass=film_mass;watermark=last_id
    end if
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    if(any(accel/=0)) then
      call force_kick(gas,temp,drops,valid)
      if(.not.valid) return
    end if
    first_wall=0;second_wall=0;boundary=0;wall=0;chemical=0
    if(outflow.and..not.impacts) then
      call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,gas,temp,drops,masses, &
        .5_dp*dt,first_steps,first_wall,valid,last_id=watermark,outer_faces=outer_faces,outer_change=first_exterior)
    else if (impacts) then
      if(outflow) then
        call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,gas,temp,drops,masses, &
          .5_dp*dt,first_steps,first_wall,valid,wall_options,patches,patch_mass,watermark,rejected, &
          outer_faces,first_exterior)
      else
        call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,gas,temp,drops,masses, &
          .5_dp*dt,first_steps,first_wall,valid,wall_options,patches,patch_mass,watermark,rejected)
      end if
      if (rejected) then
        if (present(model_rejected)) model_rejected=.true.
        return
      end if
    else if(present(last_id)) then
      call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,gas,temp,drops,masses, &
        .5_dp*dt,first_steps,first_wall,valid,last_id=watermark)
    else
      call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,gas,temp,drops,masses, &
        .5_dp*dt,first_steps,first_wall,valid)
    end if
    if (.not.valid) return
    next_gas=gas;next_temp=temp
    call advance_reactive_eb_full_3d(species,reactions,transport,gas,temp,g,solver,redistribution,dt, &
      chemistry,rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion,next_gas,next_temp, &
      theta,valid,target_volume_fraction,chemistry_integrator,boundary,wall,chemical)
    if (.not.valid) return
    if(outflow.and..not.impacts) then
      call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,next_gas,next_temp,drops,masses, &
        .5_dp*dt,second_steps,second_wall,valid,last_id=watermark,outer_faces=outer_faces,outer_change=second_exterior)
    else if (impacts) then
      if(outflow) then
        call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,next_gas,next_temp,drops,masses, &
          .5_dp*dt,second_steps,second_wall,valid,wall_options,patches,patch_mass,watermark,rejected, &
          outer_faces,second_exterior)
      else
        call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,next_gas,next_temp,drops,masses, &
          .5_dp*dt,second_steps,second_wall,valid,wall_options,patches,patch_mass,watermark,rejected)
      end if
      if (rejected) then
        if (present(model_rejected)) model_rejected=.true.
        return
      end if
    else if(present(last_id)) then
      call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,next_gas,next_temp,drops,masses, &
        .5_dp*dt,second_steps,second_wall,valid,last_id=watermark)
    else
      call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,next_gas,next_temp,drops,masses, &
        .5_dp*dt,second_steps,second_wall,valid)
    end if
    if (.not.valid) return
    before_split=watermark
    if (options%maximum_multiplicity>0) then
      call split_spray_parcels(drops,options%maximum_multiplicity,min(125000,options%maximum_parcels), &
        valid,watermark,masses)
      if (.not.valid) return
    end if
    if(any(accel/=0)) then
      call force_kick(next_gas,next_temp,drops,valid)
      if(.not.valid) return
    end if
    boundary=boundary+first_exterior+second_exterior
    wall(1:5)=wall(1:5)+first_wall+second_wall
    state=next_gas;temperature=next_temp;parcels=drops;component_mass=masses
    if (impacts) then
      films=patches;film_mass=patch_mass
    end if
    if (present(last_id)) last_id=watermark
    if (present(split_created)) split_created=watermark-before_split
    minimum_theta=theta;boundary_change=boundary;wall_change=wall;chemistry_change=chemical
    if(present(source_change)) source_change=force_ledger
    particle_steps=first_steps+second_steps;ok=.true.
  contains
    subroutine force_kick(gas_state,gas_temperature,particles,accepted)
      real(dp), intent(inout) :: gas_state(:,:,:,:)
      real(dp), intent(in) :: gas_temperature(:,:,:)
      type(spray_parcel), intent(inout) :: particles(:)
      logical, intent(out) :: accepted
      real(dp) :: dv(3),momentum(3),energy,mass,volume,velocity(3),checked_t,sound,primitive(nv)
      integer :: i,j,k,p
      accepted=.false.;dv=.5_dp*dt*accel
      if(.not.g%is_valid()) return
      if(any(shape(gas_state)/=[nv,g%nx,g%ny,g%nz])) return
      if(any(shape(gas_temperature)/=[g%nx,g%ny,g%nz])) return
      do k=1,g%nz
      do j=1,g%ny
      do i=1,g%nx
        if(g%volume_fraction(i,j,k)<=0) cycle
        mass=gas_state(1,i,j,k)
        if(.not.ieee_is_finite(mass).or.mass<=0) return
        if(maxval(abs(gas_state(2:4,i,j,k)/mass+dv))>1e7_dp) return
        momentum=mass*dv
        energy=dot_product(gas_state(2:4,i,j,k)+.5_dp*momentum,dv)
        gas_state(2:4,i,j,k)=gas_state(2:4,i,j,k)+momentum
        gas_state(5,i,j,k)=gas_state(5,i,j,k)+energy
        volume=g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz
        force_ledger(2:4)=force_ledger(2:4)+volume*momentum
        force_ledger(5)=force_ledger(5)+volume*energy
        call reactive_conserved_to_primitive(species,gas_state(:,i,j,k),gas_temperature(i,j,k), &
          primitive,checked_t,sound,accepted)
        if(.not.accepted) return
        accepted=.false.
        if(abs(checked_t-gas_temperature(i,j,k))>1e-9_dp*max(1._dp,abs(gas_temperature(i,j,k)))) return
      end do
      end do
      end do
      do p=1,size(particles)
        if(particles(p)%mass<=0) cycle
        mass=particles(p)%mass*particles(p)%multiplicity;velocity=particles(p)%velocity
        particles(p)%velocity=velocity+dv
        if(.not.valid_spray_parcel(particles(p))) return
        force_ledger(2:4)=force_ledger(2:4)+mass*dv
        force_ledger(5)=force_ledger(5)+mass*dot_product(velocity+.5_dp*dv,dv)
      end do
      accepted=all(ieee_is_finite(force_ledger))
    end subroutine
  end subroutine
end module
