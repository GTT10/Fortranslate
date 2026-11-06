module spray_mpi_amr_reactive_mod
  use les_filter_3d_mod, only: les_filter_weights
  use les_mpi_amr_tree_mod, only: advance_mpi_amr_tree_les_3d
  use les_mpi_amr_subcycle_mod, only: advance_mpi_amr_tree_subcycle_les_3d
  use spray_mpi_mod, only: split_mpi_spray_parcels
  use mpi_f08
  use precision_mod, only: dp
  use dual_energy_mod, only: dual_energy_policy
  use mpi_amr_patch_tree_split_3d_mod, only: agree_mpi_amr_tree_energy_policy_3d,close_mpi_amr_tree_held_energy_3d
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use spray_mpi_amr_balance_mod, only: sum_mpi_spray_rotating_potential
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use amr_patch_tree_split_3d_mod, only: amr_sdc_diffusion_evaluator
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options,valid_spray_parcel,valid_spray_options
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options
  use spray_mpi_amr_film_step_mod, only: advance_mpi_amr_film_coupling
  use spray_mpi_amr_film_step_mod, only: spray_sdc_film_stage,advance_mpi_amr_recursive_film_stage, &
    append_mpi_amr_film_births,finish_mpi_eb_film_level
  use spray_mpi_amr_film_mod, only: spray_level_film_population,initialize_mpi_amr_film_levels,migrate_mpi_amr_wall_films
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use spray_mpi_amr_step_mod, only: advance_mpi_amr_mixture_coupling
  use mpi_amr_patch_tree_split_3d_mod, only: advance_mpi_amr_tree_split_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use mpi_amr_patch_tree_transport_3d_mod, only: agree_mpi_constant_transport_coefficients
  use spray_mpi_amr_tree_mod, only: spray_level_population,initialize_mpi_amr_sdc_populations,collect_mpi_amr_level_particles
  use spray_mpi_amr_step_mod, only: spray_sdc_population_stage,advance_mpi_amr_recursive_population_stage, &
    convert_mpi_amr_tree_rotating_energy
  use mpi_amr_patch_tree_split_3d_mod, only: advance_mpi_amr_tree_sdc_3d,advance_mpi_amr_tree_sdc_subcycle_3d
  use spray_mpi_amr_tree_mod, only: redistribute_mpi_amr_subtree_particles
  use mpi_amr_patch_tree_ghost_3d_mod, only: fill_mpi_amr_tree_ghosts_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  implicit none
  private
  public :: advance_mpi_amr_spray_reactive,advance_mpi_amr_spray_sdc
contains
  subroutine advance_mpi_amr_spray_sdc(comm,tree,owners,species,reactions,transport,liquids,options,state,temperature, &
      parcels,component_mass,lower,spacing,dt,iterations,chemistry,transport_flags,force,rtol,atol,faces, &
      boundary_change,source_change,gas_particle_change,particle_wall_change,particle_boundary_change, &
      maximum_correction,minimum_theta,last_id,ok,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
      outlet_pressure,inflow_mach,inflow_pressure,difmag,constant_coefficients,particle_acceleration,particle_force_change, &
      subcycling,evaluate_diffusion,evaluate_reaction,wall_films,film_mass,film_settings,film_steps,film_impacts, &
      reaction_source,heat_release,mol_mode,mol_order,mol_lax_friedrichs,frame,frame_energy_change, &
      internal_energy,energy_correction,energy_policy,internal_energy_pressure_work_change)
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    real(dp), optional, intent(inout) :: energy_correction,internal_energy_pressure_work_change
    real(dp), allocatable :: held_work
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    type(amr_tree_state_3d), allocatable :: held
    type(dual_energy_policy), allocatable :: held_policy
    real(dp), allocatable :: held_prefix
    type(dual_energy_policy) :: agreed_policy
    logical :: has_held
    procedure(amr_sdc_diffusion_evaluator), optional :: evaluate_diffusion
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    integer, optional, intent(out) :: film_steps
    logical, optional, intent(in) :: film_impacts
    type(spray_level_film_population), allocatable :: film_physical(:),film_virtuals(:),film_ghosts(:)
    type(spray_sdc_film_stage), allocatable :: film_work(:)
    type(spray_wall_film), allocatable :: final_films(:)
    real(dp), allocatable :: final_film_mass(:,:)
    integer :: film_flags,film_low,film_high,film_count,cursor,film_kicks
    logical :: impacts
    logical, optional, intent(in) :: subcycling,mol_mode,mol_lax_friedrichs
    integer, optional, intent(in) :: mol_order
    integer :: particle_iterations
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(inout) :: frame_energy_change
    real(dp) :: potential_before,potential_after,frame_source,frame_boundary,frame_change
    real(dp), allocatable :: frame_lower(:)
    logical :: mol
    logical :: recursive_step
    real(dp) :: callback_time
    type(amr_tree_state_3d), allocatable :: ancestor_first,ancestor_last
    real(dp), allocatable :: ancestor_fraction(:)
    real(dp), optional, intent(in) :: particle_acceleration(3)
    real(dp), optional, intent(inout) :: particle_force_change(5)
    real(dp) :: accel(3),particle_force(5),accel_low(3),accel_high(3)
    integer :: force_flags,force_low,force_high
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),iterations
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),dt,force(3),rtol,atol
    logical, intent(in) :: chemistry,transport_flags(4)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(inout) :: boundary_change(:),source_change(:),gas_particle_change(:)
    real(dp), intent(inout) :: particle_wall_change(5),particle_boundary_change(:),maximum_correction,minimum_theta
    integer(int64), intent(inout) :: last_id
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6),difmag,constant_coefficients(:)
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source,heat_release
    type(amr_tree_state_3d), allocatable :: rates,heat
    type(amr_tree_state_3d) :: gas,temp,correction,ordinary
    type(spray_level_population), allocatable :: physical(:),virtuals(:),ghosts(:)
    type(spray_sdc_population_stage), allocatable :: work(:)
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    real(dp) :: boundary(size(boundary_change)),external(size(source_change)),exchange(size(gas_particle_change))
    real(dp) :: particle_wall(5),outside(size(particle_boundary_change)),largest,theta,integral(size(gas_particle_change))
    integer(int64) :: watermark
    integer :: rank,status,nv,node,level,n(3),i,j,k
    logical :: valid,agreed,periodic(3)
    ok=.false.
    if(present(film_steps)) film_steps=0
    impacts=.false.;film_kicks=0
    film_flags=merge(1,0,present(wall_films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))
    call MPI_Allreduce(film_flags,film_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flags,film_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.film_low/=film_high) return
    if(film_flags/=0.and.film_flags/=7) return
    valid=.true.
    if(present(wall_films)) then
      valid=allocated(wall_films).and.allocated(film_mass)
      impacts=any(film_settings%enabled)
    endif
    if(present(film_impacts)) impacts=film_impacts
    valid=valid.and.(.not.impacts.or.present(wall_films))
    call agree(valid);if(.not.valid) return
    film_flags=merge(1,0,impacts)
    call MPI_Allreduce(film_flags,film_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flags,film_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.film_low/=film_high) return
    recursive_step=.false.;callback_time=0;mol=.false.
    if(present(mol_mode)) mol=mol_mode
    call agree_mpi_amr_tree_energy_policy_3d(mol,internal_energy,energy_correction,energy_policy,comm, &
      agreed_policy,has_held,valid)
    if(.not.valid) return
    particle_iterations=iterations
    ! Pinned MOL constructs particle/film sources once at each endpoint.
    ! Later reacting corrections reevaluate only hydro/diffusion/reaction.
    ! Keep the requested gas count for pinned flux registration.
    if(mol) particle_iterations=1
    if(present(subcycling)) recursive_step=subcycling
    ok=.false.;nv=5+size(species);accel=0;particle_force=0
    force_flags=merge(1,0,present(particle_acceleration))+2*merge(1,0,present(particle_force_change))+ &
      4*merge(1,0,recursive_step)+8*merge(1,0,mol)+16*merge(1,0,present(frame))+ &
      32*merge(1,0,present(frame_energy_change))+64*merge(1,0,present(internal_energy_pressure_work_change))
    call MPI_Allreduce(force_flags,force_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(force_flags,force_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.force_low/=force_high) return
    if(present(internal_energy_pressure_work_change).and..not.has_held) return
    if(present(internal_energy_pressure_work_change)) allocate(held_work,source=0._dp)
    if(present(particle_acceleration)) accel=particle_acceleration
    valid=all(ieee_is_finite(accel)).and.all(abs(accel)<=1e6_dp)
    valid=valid.and.(.not.any(accel/=0).or.present(particle_force_change))
    call agree(valid);if(.not.valid) return
    call MPI_Allreduce(accel,accel_low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(accel,accel_high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(accel_low/=accel_high)) return
    valid=size(boundary_change)==nv.and.size(source_change)==nv.and.size(gas_particle_change)==nv
    valid=valid.and.size(particle_boundary_change)==nv.and.iterations>=1
    valid=valid.and.allocated(parcels).and.allocated(component_mass)
    call agree(valid);if(.not.valid) return
    call MPI_Comm_rank(comm,rank,status);if(status/=MPI_SUCCESS) return
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,valid)
    if(.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,valid);if(.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,valid);if(.not.valid) return
    if(has_held) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,valid)
      if(.not.valid) return
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        if(.not.all(ieee_is_finite(internal_energy%fields(node)%values))) valid=.false.
      enddo
      call agree(valid);if(.not.valid) return
      allocate(held,source=internal_energy);allocate(held_prefix,source=energy_correction)
      allocate(held_policy,source=agreed_policy)
    endif
    potential_before=0;potential_after=0;frame_source=0;frame_boundary=0;frame_change=0
    if(present(frame)) then
      valid=mol.and.frame%is_valid()
      call agree(valid);if(.not.valid) return
      allocate(frame_lower(3));frame_lower=lower
      ! Validate identical frame/coordinates before any private stage begins.
      ordinary=state
      call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,ordinary,lower,spacing,frame,.false.,valid)
      if(.not.valid) return
    endif
    if(present(frame)) then
      call sum_mpi_spray_rotating_potential(comm,frame,parcels,potential_before,valid,wall_films)
      if(.not.valid) return
    endif
    periodic=[faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']
    call initialize_mpi_amr_sdc_populations(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
      2,options%maximum_parcels,physical,virtuals,ghosts,valid)
    if(.not.valid) return
    allocate(work(tree%level_count()))
    if(present(wall_films)) then
      call initialize_mpi_amr_film_levels(comm,tree,owners,wall_films,film_mass,lower,spacing, &
        film_settings%maximum_films,film_physical,valid,periodic)
      if(.not.valid) return
      allocate(film_virtuals(size(work)),film_ghosts(size(work)),film_work(size(work)))
      do level=1,size(work)
        allocate(film_virtuals(level)%films(0),film_virtuals(level)%component_mass(size(liquids),0), &
          film_virtuals(level)%nodes(0),film_virtuals(level)%probes(3,0))
        film_ghosts(level)=film_virtuals(level)
      enddo
    endif
    if(present(reaction_source)) then
      allocate(rates);rates=reaction_source
    end if
    if(present(heat_release)) then
      allocate(heat);heat=heat_release
    end if
    gas=state;temp=temperature;correction=state
    do node=1,size(owners)
      if(owners(node)==rank) correction%fields(node)%values=0
    enddo
    watermark=last_id;particle_wall=0;outside=0;exchange=0;boundary=0;external=0;largest=0;theta=1
    ! Stored Godunov also advances complete parent intervals before children.
    ! Preserve its old/new Nv+T+H ancestor samples and post-level lifecycle.
    if(recursive_step.or.mol.or.has_held) then
    call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,gas,temp,spacing,dt,faces, &
      iterations,chemistry,transport_flags,force,rtol,atol,boundary,external,largest,theta,comm,valid, &
      inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_temperature=wall_temperature, &
      wall_velocity=wall_velocity,outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure, &
      difmag=difmag,constant_coefficients=constant_coefficients,reaction_source=rates,heat_release=heat, &
      halo_source=recursive_particle_source,auxiliary_change=exchange,post_level=particle_post_level, &
      subcycling=recursive_step, &
      evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
      mol_mode=mol_mode,mol_order=mol_order,mol_lax_friedrichs=mol_lax_friedrichs, &
      frame=frame,domain_lower=frame_lower,internal_energy=held,energy_correction=held_prefix, &
      energy_policy=held_policy,auxiliary_samples_internal_energy=has_held, &
      internal_energy_pressure_work_change=held_work)
    else
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,gas,temp,spacing,dt,faces, &
      iterations,chemistry,transport_flags,force,rtol,atol,boundary,external,largest,theta,comm,valid, &
      inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_temperature=wall_temperature, &
      wall_velocity=wall_velocity,outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure, &
      difmag=difmag,constant_coefficients=constant_coefficients,reaction_source=rates,heat_release=heat, &
      auxiliary_source=particle_source,auxiliary_change=exchange, &
      evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
      mol_mode=mol_mode,mol_order=mol_order,mol_lax_friedrichs=mol_lax_friedrichs, &
      frame=frame,domain_lower=frame_lower,internal_energy=held,energy_correction=held_prefix, &
      energy_policy=held_policy,auxiliary_samples_internal_energy=has_held, &
      internal_energy_pressure_work_change=held_work)
    endif
    if(.not.valid) return
    ! Replace provisional copy-supported leaf exchange by actual physical exchange.
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(covered(tree,node,[i,j,k])) cycle
        gas%fields(node)%values(:,i,j,k)=gas%fields(node)%values(:,i,j,k)+correction%fields(node)%values(:,i,j,k)
      enddo
      enddo
      enddo
    enddo
    call integrate_correction(integral,valid);if(.not.valid) return
    exchange=exchange+integral
    call average_down_mpi_amr_tree_3d(tree,owners,gas,nv,comm,valid);if(.not.valid) return
    if(has_held) then
      call close_mpi_amr_tree_held_energy_3d(species,tree,owners,gas,temp,held,spacing,held_policy,held_prefix, &
        comm,valid,frame=frame,domain_lower=frame_lower)
    else
    if(present(frame)) then
      ordinary=gas
      call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,ordinary,lower,spacing,frame,.false.,valid)
      if(.not.valid) return
      call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,ordinary,temp,comm,valid)
    else
      call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,gas,temp,comm,valid)
    endif
    endif
    if(.not.valid) return
    ! EOS evaluation accepts signed roundoff for intermediate states. Published
    ! states must satisfy the strict species contract used by reactors/restarts.
    ! Reject the whole private interval so the caller can retry its timestep.
    valid=.true.
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      if(any(gas%fields(node)%values(6:,:,:,:)<0)) valid=.false.
    enddo
    call agree(valid);if(.not.valid) return
    call collect_mpi_amr_level_particles(comm,tree,owners,physical,lower,spacing,periodic, &
      options%maximum_parcels,drops,masses,valid)
    if(.not.valid) return
    if(present(wall_films)) then
      film_count=0
      do level=1,size(work)
        valid=size(film_physical(level)%films)<=film_settings%maximum_films-film_count
        call agree(valid);if(.not.valid) return
        film_count=film_count+size(film_physical(level)%films)
      enddo
      allocate(final_films(film_count),final_film_mass(size(liquids),film_count));cursor=0
      do level=1,size(work)
        film_count=size(film_physical(level)%films)
        final_films(cursor+1:cursor+film_count)=film_physical(level)%films
        final_film_mass(:,cursor+1:cursor+film_count)=film_physical(level)%component_mass
        cursor=cursor+film_count
      enddo
      call migrate_mpi_amr_wall_films(comm,tree,owners,final_films,final_film_mass,lower,spacing, &
        film_settings%maximum_films,valid,periodic=periodic)
      if(.not.valid) return
    endif
    if(options%maximum_multiplicity>0.and.options%etab_parcel_factor<0) then
      call split_mpi_spray_parcels(comm,drops,options%maximum_multiplicity,options%maximum_parcels,watermark,valid,masses)
      if(.not.valid) return
    endif
    if(present(frame)) then
      call sum_mpi_spray_rotating_potential(comm,frame,drops,potential_after,valid,final_films)
      if(.not.valid) return
      ! E_rot=E-M*Phi: gas sources use deposition-cell Phi, while liquid
      ! inventory uses particle/film positions and exterior flux uses exit Phi.
      frame_change=frame_source+potential_before-potential_after-frame_boundary
      valid=ieee_is_finite(frame_change);call agree(valid);if(.not.valid) return
    endif
    if(present(reaction_source)) reaction_source=rates
    if(present(heat_release)) heat_release=heat
    if(has_held) then
      internal_energy=held;energy_correction=held_prefix
    endif
    if(present(internal_energy_pressure_work_change)) internal_energy_pressure_work_change=held_work
    state=gas;temperature=temp
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    if(present(wall_films)) then
      call move_alloc(final_films,wall_films);call move_alloc(final_film_mass,film_mass)
      if(present(film_steps)) film_steps=film_kicks
    endif
    boundary_change=boundary;source_change=external;gas_particle_change=exchange
    particle_wall_change=particle_wall;particle_boundary_change=outside
    if(present(particle_force_change)) particle_force_change=particle_force
    if(present(frame_energy_change)) frame_energy_change=frame_change
    maximum_correction=largest;minimum_theta=theta;last_id=watermark;ok=.true.
  contains
    subroutine agree(good)
      logical, intent(inout) :: good
      call MPI_Allreduce(good,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      good=status==MPI_SUCCESS.and.agreed
    end subroutine
    subroutine recursive_particle_source(u,t,interval,stage,depth,time,rate,source_comm,good,first,last,fractions)
      type(amr_tree_state_3d), intent(in) :: u,t,first,last
      type(amr_tree_state_3d), intent(inout) :: rate
      real(dp), intent(in) :: interval,time,fractions(:)
      integer, intent(in) :: stage,depth
      type(MPI_Comm), intent(in) :: source_comm
      logical, intent(out) :: good
      callback_time=time;ancestor_first=first;ancestor_last=last;ancestor_fraction=fractions
      call particle_source(u,t,interval,stage,depth,rate,source_comm,good)
    end subroutine
    subroutine particle_post_level(depth,time,interval,iteration,ncycle,post_comm,good)
      integer, intent(in) :: depth,iteration,ncycle
      real(dp), intent(in) :: time,interval
      type(MPI_Comm), intent(in) :: post_comm
      logical, intent(out) :: good
      integer :: l
      good=.true.
      do l=depth+1,size(work)
        good=good.and.work(l)%next_stage==0.and. &
          abs(work(l)%time+work(l)%interval-time-interval)<=32*epsilon(dt)*max(dt,abs(time+interval))
      enddo
      call agree(good);if(.not.good) return
      if(present(wall_films)) then
        call finish_mpi_eb_film_level(post_comm,depth,time,interval,iteration,ncycle, &
          film_virtuals,film_ghosts,film_work,good)
        if(.not.good) return
      endif
      if(depth==0.or.(iteration<ncycle.and.depth<size(work)-1)) then
        call redistribute_mpi_amr_subtree_particles(post_comm,tree,owners,physical,lower,spacing,periodic, &
          depth,iteration,options%maximum_parcels,good)
        if(.not.good) return
      endif
      virtuals(depth+1)%parcels=virtuals(depth+1)%parcels(:0)
      virtuals(depth+1)%component_mass=virtuals(depth+1)%component_mass(:,:0)
      virtuals(depth+1)%nodes=virtuals(depth+1)%nodes(:0)
      if(iteration==ncycle) then
        ghosts(depth+1)%parcels=ghosts(depth+1)%parcels(:0)
        ghosts(depth+1)%component_mass=ghosts(depth+1)%component_mass(:,:0)
        ghosts(depth+1)%nodes=ghosts(depth+1)%nodes(:0)
      endif
    end subroutine
    subroutine particle_source(u,t,interval,stage,depth,rate,source_comm,good)
      type(amr_tree_state_3d), intent(in) :: u,t
      type(amr_tree_state_3d), intent(inout) :: rate
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(MPI_Comm), intent(in) :: source_comm
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: packed,future,halo,halo_gas,halo_temp,part,remainder,actual
      type(amr_tree_state_3d), allocatable :: halo_held
      character(len=:), allocatable :: held_limiter
      real(dp), allocatable :: inlet(:,:),full_fraction(:)
      type(amr_tree_state_3d) :: film_rate,film_remainder,film_actual
      type(spray_wall_film), allocatable :: impact_films(:)
      real(dp), allocatable :: impact_mass(:,:)
      real(dp) :: wall(5),film_wall(5),exterior(nv)
      integer :: p,l,d(3),old_films,np,face,allocation_status
      good=depth>=-1;call agree(good);if(.not.good) return
      np=nv+1+merge(1,0,has_held)
      if(has_held) then
        allocate(character(len=6) :: held_limiter);held_limiter='linear'
        allocate(halo_held)
      endif
      packed=u
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        d=tree%nodes(p)%shape
        deallocate(packed%fields(p)%values)
        allocate(packed%fields(p)%values(np,d(1),d(2),d(3)))
        packed%fields(p)%values(:nv,:,:,:)=u%fields(p)%values(:nv,:,:,:)
        packed%fields(p)%values(nv+1,:,:,:)=t%fields(p)%values(1,:,:,:)
        if(has_held) packed%fields(p)%values(np,:,:,:)=u%fields(p)%values(nv+1,:,:,:)
      enddo
      if(present(inflow_state)) then
        good=present(inflow_temperature).and.size(inflow_state,1)==nv.and.size(inflow_state,2)==6
        call agree(good);if(.not.good) return
        good=size(inflow_temperature)==6;call agree(good);if(.not.good) return
        allocate(inlet(np,6));inlet(:nv,:)=inflow_state;inlet(nv+1,:)=inflow_temperature
        if(has_held) then
          inlet(np,:)=0
          do face=1,6
            if(inflow_state(1,face)>0) then
              inlet(np,face)=inflow_state(5,face)-.5_dp*sum(inflow_state(2:4,face)**2)/inflow_state(1,face)
            endif
          enddo
          good=all(ieee_is_finite(inlet(np,:)));call agree(good);if(.not.good) return
        endif
      endif
      if(depth>0) then
        do p=1,size(owners)
          if(owners(p)/=rank.or.tree%nodes(p)%level>=depth) cycle
          packed%fields(p)%values=ancestor_first%fields(p)%values
        enddo
        future=packed
        do p=1,size(owners)
          if(owners(p)/=rank.or.tree%nodes(p)%level>=depth) cycle
          future%fields(p)%values=ancestor_last%fields(p)%values
        enddo
        allocate(full_fraction(tree%level_count()-1));full_fraction=0
        full_fraction(:depth)=ancestor_fraction
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,packed,np,faces,halo,source_comm,good, &
        viscous_wall=.true.,inflow_state=inlet,wall_velocity=wall_velocity,species_count=size(species), &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        ghost_width=2,coarse_limiter=held_limiter, &
        limit_internal_energy=has_held,stored_internal_energy=has_held,coarse_future=future,coarse_fractions=full_fraction, &
        frame=frame,domain_lower=frame_lower,root_spacing=spacing)
      else
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,packed,np,faces,halo,source_comm,good, &
        viscous_wall=.true.,inflow_state=inlet,wall_velocity=wall_velocity,species_count=size(species), &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        ghost_width=2,coarse_limiter=held_limiter, &
        limit_internal_energy=has_held,stored_internal_energy=has_held,frame=frame,domain_lower=frame_lower,root_spacing=spacing)
      endif
      if(.not.good) return
      halo_gas=halo;halo_temp=halo
      if(has_held) halo_held=halo
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        if(.not.allocated(halo%fields(p)%values)) cycle
        d=tree%nodes(p)%shape
        deallocate(halo_gas%fields(p)%values,halo_temp%fields(p)%values)
        allocate(halo_gas%fields(p)%values(nv,-1:d(1)+2,-1:d(2)+2,-1:d(3)+2))
        allocate(halo_temp%fields(p)%values(1,-1:d(1)+2,-1:d(2)+2,-1:d(3)+2))
        halo_gas%fields(p)%values=halo%fields(p)%values(:nv,:,:,:)
        halo_temp%fields(p)%values=halo%fields(p)%values(nv+1:nv+1,:,:,:)
        if(has_held) then
          deallocate(halo_held%fields(p)%values)
          allocate(halo_held%fields(p)%values(1,-1:d(1)+2,-1:d(2)+2,-1:d(3)+2))
          halo_held%fields(p)%values=halo%fields(p)%values(np:np,:,:,:)
        endif
      enddo
      rate=u
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        if(has_held) then
          d=tree%nodes(p)%shape
          deallocate(rate%fields(p)%values)
          allocate(rate%fields(p)%values(nv,d(1),d(2),d(3)),stat=allocation_status)
          if(allocation_status/=0) then
            good=.false.;exit
          endif
        endif
        rate%fields(p)%values=0
      enddo
      call agree(good);if(.not.good) return
      do l=1,size(work)
        if(depth>=0.and.l/=depth+1) cycle
        if(impacts) then
          impact_films=film_physical(l)%films;impact_mass=film_physical(l)%component_mass
          old_films=size(impact_films)
        call advance_mpi_amr_recursive_population_stage(source_comm,species,transport,liquids,options,tree,owners, &
          halo_gas,halo_temp,2,l-1,2,options%maximum_parcels,physical,virtuals,ghosts,work,callback_time,interval,stage, &
          particle_iterations,faces,part,remainder,wall,exterior,good,lower,spacing,physical_source=actual,last_id=watermark, &
          acceleration=accel,frame=frame,common_step=.not.(recursive_step.or.mol.or.has_held), &
          wall_films=impact_films,film_mass=impact_mass,film_settings=film_settings,internal_energy_halo=halo_held)
        else
        call advance_mpi_amr_recursive_population_stage(source_comm,species,transport,liquids,options,tree,owners, &
          halo_gas,halo_temp,2,l-1,2,options%maximum_parcels,physical,virtuals,ghosts,work,callback_time,interval,stage, &
          particle_iterations,faces,part,remainder,wall,exterior,good,lower,spacing,physical_source=actual,last_id=watermark, &
          acceleration=accel,frame=frame,common_step=.not.(recursive_step.or.mol.or.has_held), &
            internal_energy_halo=halo_held)
        endif
        if(.not.good) return
        if(present(wall_films)) then
          call advance_mpi_amr_recursive_film_stage(source_comm,species,transport,liquids,options,film_settings, &
            tree,owners,halo_gas,halo_temp,2,l-1,2,film_settings%maximum_films,film_physical,film_virtuals,film_ghosts, &
            film_work,lower,spacing,callback_time,interval,stage,particle_iterations,faces, &
            film_rate,film_remainder,film_wall,good, &
            physical_source=film_actual,frame=frame,common_step=.not.(recursive_step.or.mol.or.has_held), &
            internal_energy_halo=halo_held)
          if(.not.good) return
          if(impacts.and.stage==0) then
            call append_mpi_amr_film_births(source_comm,tree,owners,lower,spacing,2,film_settings%maximum_films, &
              work(l)%midpoint_films(old_films+1:),work(l)%midpoint_film_mass(:,old_films+1:),film_work(l),good,periodic)
            if(.not.good) return
          endif
          call add_rate(part,film_rate,good);if(.not.good) return
          call add_rate(actual,film_actual,good);if(.not.good) return
          good=all(ieee_is_finite(wall)).and.all(ieee_is_finite(film_wall))
          if(good) good=.not.any(wall>huge(1._dp)-max(film_wall,0._dp)).and. &
            .not.any(wall< -huge(1._dp)-min(film_wall,0._dp))
          call agree(good);if(.not.good) return
          wall=wall+film_wall
          if(stage==0.or.stage==particle_iterations) film_kicks=film_kicks+1
        endif
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          if(stage==0.or.stage==particle_iterations) then
            correction%fields(p)%values=correction%fields(p)%values+.5_dp*interval*actual%fields(p)%values
            ! Later reacting MOL passes omit the new auxiliary family; retain
            ! its completed physical half kick in the accepted gas inventory.
            if(tree%nodes(p)%level==l-1.and..not.(mol.and.chemistry.and.iterations>1.and.stage==1)) &
              correction%fields(p)%values=correction%fields(p)%values-.5_dp*interval*part%fields(p)%values
          endif
          if(tree%nodes(p)%level==l-1) rate%fields(p)%values=part%fields(p)%values
        enddo
        if(stage==particle_iterations) then
          particle_wall=particle_wall+wall;outside=outside+exterior
          particle_force=particle_force+work(l)%force_change
          frame_source=frame_source+work(l)%frame_source_change
          frame_boundary=frame_boundary+work(l)%frame_boundary_change
          if(present(wall_films)) frame_source=frame_source+film_work(l)%frame_source_change
        endif
      enddo
      good=.true.
    end subroutine
    subroutine add_rate(destination,addition,good)
      type(amr_tree_state_3d), intent(inout) :: destination
      type(amr_tree_state_3d), intent(in) :: addition
      logical, intent(out) :: good
      integer :: peer
      good=.true.
      do peer=1,size(owners)
        if(owners(peer)/=rank) cycle
        associate(a=>destination%fields(peer)%values,b=>addition%fields(peer)%values)
          good=good.and.all(ieee_is_finite(a)).and.all(ieee_is_finite(b))
          if(good) good=.not.any(a>huge(1._dp)-max(b,0._dp)).and. &
            .not.any(a< -huge(1._dp)-min(b,0._dp))
        end associate
      enddo
      call agree(good);if(.not.good) return
      do peer=1,size(owners)
        if(owners(peer)==rank) destination%fields(peer)%values= &
          destination%fields(peer)%values+addition%fields(peer)%values
      enddo
    end subroutine
    subroutine integrate_correction(result,good)
      real(dp), intent(out) :: result(:)
      logical, intent(out) :: good
      real(dp) :: local(size(result)),volume
      integer :: p,ancestor,x,y,z,d(3)
      local=0
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        volume=product(spacing);ancestor=p
        do while(ancestor>1)
          volume=volume/real(tree%nodes(ancestor)%patch%refinement_ratio,dp)**3
          ancestor=tree%nodes(ancestor)%parent
        enddo
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(covered(tree,p,[x,y,z])) cycle
          local=local+volume*correction%fields(p)%values(:,x,y,z)
        enddo
        enddo
        enddo
      enddo
      call MPI_Allreduce(local,result,size(result),MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
      good=status==MPI_SUCCESS.and.all(ieee_is_finite(result))
    end subroutine
  end subroutine advance_mpi_amr_spray_sdc

  ! S/2 - F/2 - L/2 - (R/2 T/2 H T/2 R/2) - L/2 - F/2 - S/2 with atomic publication.
  ! Common tables/options and globally unique parcel IDs are caller contracts.
  subroutine advance_mpi_amr_spray_reactive(comm,tree,owners,species,reactions,transport,liquids,options, &
      state,temperature,parcels,component_mass,lower,spacing,dt,reconstruction,limiter,solver, &
      chemistry,viscosity,conduction,diffusion,barodiffusion,rtol,atol,particle_steps, &
      maximum_correction,minimum_theta,ok,last_id,chemistry_integrator,subcycling, &
      boundary_faces,boundary_change,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
      wall_films,film_mass,film_settings,film_steps,les_settings,les_model,outlet_pressure, &
      acceleration,source_change,split_created, &
      les_test_model,les_test_ratio,les_grid_ratio,les_harmonic_mean,inflow_mach,inflow_pressure,difmag, &
      constant_coefficients)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    logical :: constant_agreed
    integer(int64), optional, intent(out) :: split_created
    integer(int64) :: before_split
    real(dp) :: low_multiplicity,high_multiplicity
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(inout) :: source_change(:)
    real(dp) :: accel(3),low_accel(3),high_accel(3)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: difmag
    real(dp) :: pressure(6),pressure_min(6),pressure_max(6)
    ! Static LES settings are [eddy coefficient, PrT, ScT, CI]; dynamic is enabled by model.
    real(dp), optional, intent(in) :: les_settings(4)
    character(len=*), optional, intent(in) :: les_model,les_test_model
    integer, optional, intent(in) :: les_test_ratio,les_grid_ratio
    logical, optional, intent(in) :: les_harmonic_mean
    character(len=64) :: filter,reference_filter
    integer :: filter_controls(3),filter_low(3),filter_high(3)
    real(dp), allocatable :: weights(:)
    logical :: filter_ok
    real(dp) :: sgs(4),sgs_min(4),sgs_max(4),les_correction,les_max
    character(len=16) :: model,reference_model
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(nasa7_species), intent(in) :: species(:)
    real(dp) :: les_boundary(reactive_nvar(size(species))),les_ledger(reactive_nvar(size(species)))
    real(dp) :: force_ledger(reactive_nvar(size(species))),constraint_change(5)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),dt,rtol,atol
    character(len=*), intent(in) :: reconstruction,limiter,solver
    logical, intent(in) :: chemistry,viscosity,conduction,diffusion,barodiffusion
    integer, intent(out) :: particle_steps
    real(dp), intent(inout) :: maximum_correction,minimum_theta
    logical, intent(out) :: ok
    integer(int64), optional, intent(inout) :: last_id
    character(len=*), optional, intent(in) :: chemistry_integrator
    logical, optional, intent(in) :: subcycling
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(inout) :: boundary_change(:)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    integer, optional, intent(out) :: film_steps
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    real(dp) :: film_boundary(reactive_nvar(size(species))),film_ledger(reactive_nvar(size(species)))
    integer :: film_first,film_second,film_flag,low_flag,high_flag
    type(amr_tree_state_3d) :: gas,temp
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    real(dp) :: boundary(reactive_nvar(size(species))),correction,theta
    real(dp) :: first_boundary(reactive_nvar(size(species))),second_boundary(reactive_nvar(size(species)))
    integer(int64) :: candidate_id
    integer :: first_steps,second_steps,ierr
    logical :: valid,agreed,recursive
    integer :: recursive_flag,recursive_low,recursive_high
    character(len=32) :: faces(6)
    ok=.false.;particle_steps=0;film_first=0;film_second=0;film_ledger=0
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,constant_agreed)
    if(.not.constant_agreed) return
    if (present(split_created)) split_created=0
    if (present(film_steps)) film_steps=0
    sgs=[0._dp,.7_dp,.7_dp,0._dp];model='smagorinsky';les_ledger=0;les_max=0
    if (present(les_settings)) sgs=les_settings
    if (present(les_model)) model=les_model
    valid=allocated(parcels).and.allocated(component_mass).and.ieee_is_finite(dt)
    if (valid) valid=dt>0.and.dt<=1e6_dp
    valid=valid.and.all(ieee_is_finite(sgs))
    valid=valid.and.sgs(1)>=0.and.sgs(1)<=2.and.sgs(4)>=0.and.sgs(4)<=1
    valid=valid.and.minval(sgs(2:3))>=1e-6_dp.and.maxval(sgs(2:3))<=1e6_dp
    valid=valid.and.(model=='smagorinsky'.or.model=='wale'.or.model=='vreman'.or.model=='dynamic')
    valid=valid.and.valid_spray_options(options)
    valid=valid.and.(options%maximum_multiplicity==0.or.present(last_id))
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    if (present(boundary_change)) valid=valid.and.size(boundary_change)==size(boundary)
    valid=valid.and.(present(wall_films).eqv.present(film_mass))
    valid=valid.and.(present(wall_films).eqv.present(film_settings))
    if (present(wall_films)) valid=valid.and.allocated(wall_films).and.present(last_id)
    if (present(film_mass)) valid=valid.and.allocated(film_mass)
    faces='periodic'
    if (present(boundary_faces)) then
      valid=valid.and.all(len_trim(boundary_faces)<=len(faces))
      faces=boundary_faces
    end if
    accel=0;force_ledger=0
    if (present(acceleration)) accel=acceleration
    valid=valid.and.all(ieee_is_finite(accel)).and.all(abs(accel)<=1e6_dp)
    if (any(accel/=0).or.options%fixed_parts) valid=valid.and.present(source_change)
    if (present(source_change)) valid=valid.and.size(source_change)==size(boundary)
    pressure=0
    if (present(outlet_pressure)) pressure=outlet_pressure
    valid=valid.and.all(ieee_is_finite(pressure))
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    valid=valid.and..not.any(pressure>0.and.faces/='outflow')
    filter='box';filter_controls=[2,1,0]
    if(present(les_test_model)) filter=les_test_model
    if(present(les_test_ratio)) filter_controls(1)=les_test_ratio
    if(present(les_grid_ratio)) filter_controls(2)=les_grid_ratio
    if(present(les_harmonic_mean)) filter_controls(3)=merge(1,0,les_harmonic_mean)
    if(model=='dynamic') then
      call les_filter_weights(filter,filter_controls(1),weights,filter_ok)
      valid=valid.and.filter_ok.and.filter_controls(2)>=1.and.filter_controls(2)<=64
      valid=valid.and.all(faces=='periodic'.or.faces=='inflow'.or.faces=='outflow'.or. &
        faces=='slip_wall'.or.faces=='no_slip_wall'.or.faces=='isothermal_wall')
    else
      if(present(les_test_model).or.present(les_test_ratio).or. &
        present(les_grid_ratio).or.present(les_harmonic_mean)) valid=.false.
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(options%maximum_multiplicity,low_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(options%maximum_multiplicity,high_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_multiplicity/=high_multiplicity) return
    call MPI_Allreduce(accel,low_accel,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(accel,high_accel,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.any(low_accel/=high_accel)) return
    call MPI_Allreduce(pressure,pressure_min,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(pressure,pressure_max,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(pressure_min/=pressure_max)) return
    call MPI_Allreduce(sgs,sgs_min,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(sgs,sgs_max,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(sgs_min/=sgs_max)) return
    call MPI_Allreduce(filter_controls,filter_low,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(filter_controls,filter_high,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(filter_low/=filter_high)) return
    reference_filter=filter
    call MPI_Bcast(reference_filter,len(reference_filter),MPI_CHARACTER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=filter==reference_filter
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    reference_model=model
    call MPI_Bcast(reference_model,len(reference_model),MPI_CHARACTER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=model==reference_model
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    film_flag=merge(1,0,present(wall_films))
    call MPI_Allreduce(film_flag,low_flag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,high_flag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_flag/=high_flag) return
    recursive=.false.
    if (present(subcycling)) recursive=subcycling
    recursive_flag=merge(1,0,recursive)
    call MPI_Allreduce(recursive_flag,recursive_low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(recursive_flag,recursive_high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.recursive_low/=recursive_high) return
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    if (present(wall_films)) then
      patches=wall_films;patch_mass=film_mass
      call advance_mpi_amr_film_coupling(comm,species,transport,liquids,options,film_settings,gas,temp, &
        tree,owners,patches,patch_mass,lower,spacing,.5_dp*dt,faces,film_first,film_boundary,valid)
      if (.not.valid) return
      film_ledger=film_boundary
    end if
    candidate_id=0
    if (present(last_id)) candidate_id=last_id
    if (options%etab_parcel_factor>=0) before_split=candidate_id
    if (present(wall_films)) then
      call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
        drops,masses,lower,spacing,.5_dp*dt,first_steps,valid,candidate_id,faces,first_boundary, &
        patches,patch_mass,film_settings, &
        force_change=constraint_change)
    else
      call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
        drops,masses,lower,spacing,.5_dp*dt,first_steps,valid,candidate_id,faces,first_boundary, &
        force_change=constraint_change)
    end if
    if (.not.valid) return
    force_ledger(:5)=force_ledger(:5)+constraint_change
    call force_kick()
    if (.not.valid) return
    call les_half_step()
    if (.not.valid) return
    boundary=0;correction=0;theta=1
    call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,gas,temp,spacing,dt,faces, &
      reconstruction,limiter,solver,chemistry,viscosity,conduction,diffusion,barodiffusion,rtol,atol, &
      boundary,correction,theta,comm,valid,chemistry_integrator=chemistry_integrator,subcycling=subcycling, &
      inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_temperature=wall_temperature, &
      wall_velocity=wall_velocity,outlet_pressure=pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure, &
      difmag=difmag, &
      constant_coefficients=constant_coefficients)
    if (.not.valid) return
    call les_half_step()
    if (.not.valid) return
    call force_kick()
    if (.not.valid) return
    if (present(wall_films)) then
      call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
        drops,masses,lower,spacing,.5_dp*dt,second_steps,valid,candidate_id,faces,second_boundary, &
        patches,patch_mass,film_settings, &
        force_change=constraint_change)
    else
      call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
        drops,masses,lower,spacing,.5_dp*dt,second_steps,valid,candidate_id,faces,second_boundary, &
        force_change=constraint_change)
    end if
    if (.not.valid) return
    force_ledger(:5)=force_ledger(:5)+constraint_change
    if (present(wall_films)) then
      call advance_mpi_amr_film_coupling(comm,species,transport,liquids,options,film_settings,gas,temp, &
        tree,owners,patches,patch_mass,lower,spacing,.5_dp*dt,faces,film_second,film_boundary,valid)
      if (.not.valid) return
      film_ledger=film_ledger+film_boundary
    end if
    if (options%etab_parcel_factor<0) before_split=candidate_id
    if (options%maximum_multiplicity>0.and.options%etab_parcel_factor<0) then
      call split_mpi_spray_parcels(comm,drops,options%maximum_multiplicity,options%maximum_parcels, &
        candidate_id,valid,masses)
      if (.not.valid) return
    end if
    boundary=boundary+first_boundary+second_boundary+film_ledger+les_ledger
    valid=all(ieee_is_finite(boundary)).and.all(ieee_is_finite(force_ledger))
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    if (present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    if (present(film_steps)) film_steps=film_first+film_second
    if (present(split_created)) split_created=candidate_id-before_split
    if (present(last_id)) last_id=candidate_id
    if (present(boundary_change)) boundary_change=boundary
    if (present(source_change)) source_change=force_ledger
    maximum_correction=max(correction,les_max);minimum_theta=theta
    particle_steps=first_steps+second_steps;ok=.true.
  contains
    subroutine force_kick()
      real(dp) :: local(size(boundary)),global(size(boundary))
      valid=.true.
      if (all(accel==0)) return
      call validate_mpi_amr_tree_local_3d(tree,owners,gas,size(boundary),comm,valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temp,1,comm,valid)
      if(.not.valid) return
      call local_force_kick(local)
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      valid=ierr==MPI_SUCCESS.and.agreed
      if(.not.valid) return
      call MPI_Allreduce(local,global,size(boundary),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      valid=ierr==MPI_SUCCESS
      if(valid) force_ledger=force_ledger+global
    end subroutine
    subroutine local_force_kick(local)
      real(dp), intent(out) :: local(size(boundary))
      real(dp) :: dv(3),mass,momentum(3),energy,volume,velocity(3),prim(size(boundary)),checked_t,sound
      real(dp) :: h(3,size(tree%nodes))
      integer :: rank,status,node,i,j,k,p,n(3)
      valid=.false.;local=0;dv=.5_dp*dt*accel
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
      h(:,1)=spacing
      do node=2,size(tree%nodes)
        h(:,node)=h(:,tree%nodes(node)%parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      end do
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          mass=gas%fields(node)%values(1,i,j,k)
          if(.not.ieee_is_finite(mass).or.mass<=0) return
          velocity=gas%fields(node)%values(2:4,i,j,k)/mass
          if(any(abs(velocity+dv)>1e7_dp)) return
          momentum=mass*dv
          energy=dot_product(gas%fields(node)%values(2:4,i,j,k)+.5_dp*momentum,dv)
          gas%fields(node)%values(2:4,i,j,k)=gas%fields(node)%values(2:4,i,j,k)+momentum
          gas%fields(node)%values(5,i,j,k)=gas%fields(node)%values(5,i,j,k)+energy
          call reactive_conserved_to_primitive(species,gas%fields(node)%values(:,i,j,k), &
            temp%fields(node)%values(1,i,j,k),prim,checked_t,sound,valid)
          if(.not.valid) return
          valid=.false.
          if(abs(checked_t-temp%fields(node)%values(1,i,j,k))> &
            1e-9_dp*max(1._dp,abs(temp%fields(node)%values(1,i,j,k)))) return
          if(covered(tree,node,[i,j,k])) cycle
          volume=product(h(:,node))
          local(2:4)=local(2:4)+volume*momentum;local(5)=local(5)+volume*energy
        end do
        end do
        end do
      end do
      do p=1,size(drops)
        if(options%fixed_parts) cycle
        if(drops(p)%mass<=0) cycle
        mass=drops(p)%mass*drops(p)%multiplicity;velocity=drops(p)%velocity
        drops(p)%velocity=velocity+dv
        if(.not.valid_spray_parcel(drops(p))) return
        local(2:4)=local(2:4)+mass*dv
        local(5)=local(5)+mass*dot_product(velocity+.5_dp*dv,dv)
      end do
      valid=all(ieee_is_finite(local))
    end subroutine
    subroutine les_half_step()
      real(dp) :: maximum_diffusivity
      valid=.true.
      if (sgs(1)==0.and.sgs(4)==0.and.model/='dynamic') return
      if (recursive) then
        maximum_diffusivity=0
        call advance_mpi_amr_tree_subcycle_les_3d(tree,owners,species,gas,temp,spacing,.5_dp*dt,faces, &
          sgs(1),sgs(2),sgs(3),les_boundary,les_correction,maximum_diffusivity,comm,valid, &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature, &
          wall_temperature=wall_temperature,wall_velocity=wall_velocity,model=model,ci=sgs(4),outlet_pressure=pressure, &
          test_model=les_test_model,test_ratio=les_test_ratio,grid_ratio=les_grid_ratio,harmonic_mean=les_harmonic_mean, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      else
        call advance_mpi_amr_tree_les_3d(tree,owners,species,gas,temp,spacing,.5_dp*dt,faces, &
          sgs(1),sgs(2),sgs(3),les_boundary,les_correction,comm,valid, &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature, &
          wall_temperature=wall_temperature,wall_velocity=wall_velocity,model=model,ci=sgs(4),outlet_pressure=pressure, &
          test_model=les_test_model,test_ratio=les_test_ratio,grid_ratio=les_grid_ratio,harmonic_mean=les_harmonic_mean, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      end if
      if (.not.valid) return
      les_ledger=les_ledger+les_boundary;les_max=max(les_max,les_correction)
    end subroutine
  end subroutine
end module
