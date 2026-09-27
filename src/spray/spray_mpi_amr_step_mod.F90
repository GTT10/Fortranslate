module spray_mpi_amr_step_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory,liquid_mixture_properties, &
    advance_mixture_frozen_gas,mixture_transfer_rates,mixture_trace_fraction,build_mixture_film
  use spray_breakup_mod, only: advance_spray_breakup
  use spray_coupling_3d_mod, only: valid_spray_boundary_faces,apply_parcel_boundaries
  use spray_mpi_mod, only: advance_mpi_khrt_shedding,impact_mpi_planar_walls
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,valid_wall_options
  use spray_mpi_amr_film_mod, only: migrate_mpi_amr_wall_films
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use spray_mpi_amr_field_mod, only: sample_mpi_amr_tree_spray_field,deposit_mpi_amr_tree_spray_gas, &
    deposit_mpi_amr_tree_spray_field, &
    sample_mpi_amr_particle_halo,deposit_mpi_amr_particle_halo,sum_mpi_amr_particle_source_halo
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels,migrate_eb_tree_spray_parcels, &
    redistribute_mpi_amr_ghost_particles,redistribute_mpi_amr_level_particles,spray_level_population, &
    refresh_mpi_amr_level_copies
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_eb_drift_mod, only: drift_spray_level_set
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  use spray_eb_impact_mod, only: spray_eb_wall_options
  use spray_mpi_eb_impact_mod, only: advance_mpi_level_set_impacts
  implicit none
  private
  public :: advance_mpi_amr_mixture_coupling,kick_mpi_amr_mixture_halo,advance_mpi_amr_ghost_halo
  public :: advance_mpi_eb_ghost_stage,advance_mpi_amr_level_halo,advance_mpi_eb_population_stage
  public :: advance_mpi_eb_population_sdc_stage,advance_mpi_eb_recursive_population_stage
  type, public :: spray_sdc_population_stage
    integer :: next_stage=0,iterations=0,depth=-1
    real(dp) :: time=0,interval=0
    type(spray_level_population) :: midpoint(3)
    type(amr_tree_state_3d) :: old_source,old_remainder,old_physical_source
    logical :: exports_physical=.false.
    real(dp) :: old_wall(5)=0,force_change(5)=0,acceleration(3)=0
    real(dp), allocatable :: old_boundary(:)
  end type spray_sdc_population_stage
contains
  ! Advance a caller-owned hierarchy one recursive SDC callback at a time.
  ! The complete hierarchy must stay private until its outer gas step succeeds.
  subroutine advance_mpi_eb_recursive_population_stage(comm,species,transport,components,options,tree,owners,geometry, &
      state,temperature,width,target_level,grow,maximum_parcels,physical,virtuals,ghosts,work, &
      time,interval,stage,iterations,faces, &
      source,remainder,wall_change,boundary_change,ok,nodal_phi,centroids,target_volume_fraction, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,maximum_received,physical_source,acceleration,last_id)
    integer(int64), optional, intent(inout) :: last_id
    integer(int64), allocatable :: watermark
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: owners(:),width,target_level,grow,maximum_parcels
    type(spray_level_population), intent(inout) :: physical(:),virtuals(:),ghosts(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature,centroids,target_volume_fraction
    real(dp), intent(in) :: time,interval,nodal_phi(0:,0:,0:)
    integer, intent(in) :: stage,iterations
    type(spray_sdc_population_stage), intent(inout) :: work(:)
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5),boundary_change(:)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:),inflow_state(:,:),inflow_mach(6)
    real(dp), optional, intent(in) :: inflow_temperature(6),inflow_pressure(6)
    integer, optional, intent(in) :: maximum_received
    real(dp), optional, intent(in) :: acceleration(3)
    type(amr_tree_state_3d), optional, intent(inout) :: physical_source
    type(amr_tree_state_3d), allocatable :: physical_delta
    type(spray_sdc_population_stage) :: candidate
    type(spray_level_population) :: group(3),new_virtual,next_ghost
    type(amr_tree_state_3d) :: change,leftover
    real(dp) :: wall(5),exterior(reactive_nvar(size(species))),lower(3),spacing(3),tolerance,finish
    integer :: levels,level,p,status,flags(4),low(4),high(4)
    real(dp) :: times(2),low_time(2),high_time(2)
    logical :: valid,agreed,periodic(3)
    if(present(last_id)) then
      allocate(watermark);watermark=last_id
    end if
    ok=.false.;valid=tree%is_valid()
    call agree();if(.not.agreed) return
    levels=tree%level_count();level=target_level+1
    valid=all([size(physical),size(virtuals),size(ghosts),size(work)]==levels)
    valid=valid.and.level>=1.and.level<=levels.and.all(ieee_is_finite([time,interval]))
    valid=valid.and.time>=0.and.interval>0.and.size(boundary_change)==size(exterior)
    call agree();if(.not.agreed) return
    flags=[levels,target_level,stage,iterations];times=[time,interval]
    call MPI_Allreduce(flags,low,4,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,4,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    call MPI_Allreduce(times,low_time,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(times,high_time,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low_time/=high_time)) return
    tolerance=32*epsilon(1._dp)*max(abs(time),interval,tiny(1._dp))
    if(stage==0) then
      ! Current level and descendants must start together; parents have already
      ! completed their provisional interval before recursive descent.
      do p=level,levels
        finish=work(p)%time+work(p)%interval
        valid=valid.and.work(p)%next_stage==0.and.abs(finish-time)<=tolerance
      end do
      if(level>1) then
        valid=valid.and.work(level-1)%iterations>0.and.work(level-1)%next_stage==0
        valid=valid.and.time>=work(level-1)%time-tolerance
        valid=valid.and.interval<=work(level-1)%time+work(level-1)%interval-time+tolerance
      end if
    end if
    call agree();if(.not.agreed) return
    lower=[geometry(1)%x_lower,geometry(1)%y_lower,geometry(1)%z_lower]
    spacing=[geometry(1)%dx,geometry(1)%dy,geometry(1)%dz]
    periodic=[faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']
    if(present(physical_source)) allocate(physical_delta)
    candidate=work(level);group(1)=physical(level);group(2)=virtuals(level);group(3)=ghosts(level)
    if(stage==0) then
      call refresh_mpi_amr_level_copies(comm,tree,owners,physical,lower,spacing,periodic,target_level,grow, &
        maximum_parcels,new_virtual,next_ghost,valid,nodal_phi)
      if(.not.valid) return
      group(2)=new_virtual
    end if
    call advance_mpi_eb_population_sdc_stage(comm,species,transport,components,options,tree,owners,geometry, &
      state,temperature,width,target_level,grow,maximum_parcels,group,candidate,time,interval,stage,iterations,faces, &
      change,leftover,wall,exterior,valid,nodal_phi,centroids,target_volume_fraction, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,maximum_received,physical_delta,acceleration,watermark)
    if(.not.valid) return
    physical(level)=group(1);virtuals(level)=group(2);ghosts(level)=group(3);work(level)=candidate
    if(stage==0.and.level<levels) ghosts(level+1)=next_ghost
    if(present(physical_source)) call move_alloc(physical_delta%fields,physical_source%fields)
    call move_alloc(change%fields,source%fields);call move_alloc(leftover%fields,remainder%fields)
    if(present(last_id)) last_id=watermark
    wall_change=wall;boundary_change=exterior;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
  end subroutine advance_mpi_eb_recursive_population_stage

  ! Old half kick + full drift, then a new half kick only at the final
  ! iteration, matching pinned Particle.cpp. Sources are stage rates;
  ! stage>0 ledgers contain old plus current-new physical boundary changes.
  subroutine advance_mpi_eb_population_sdc_stage(comm,species,transport,components,options,tree,owners,geometry, &
      state,temperature,width,target_level,grow,maximum_parcels,populations,work,time,interval,stage,iterations,faces, &
      source,remainder,wall_change,boundary_change,ok,nodal_phi,centroids,target_volume_fraction, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,maximum_received,physical_source,acceleration,last_id)
    integer(int64), optional, intent(inout) :: last_id
    integer(int64), allocatable :: watermark
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: owners(:),width,target_level,grow,maximum_parcels
    type(spray_level_population), intent(inout) :: populations(3)
    type(amr_tree_state_3d), intent(in) :: state,temperature,centroids,target_volume_fraction
    real(dp), intent(in) :: time,interval,nodal_phi(0:,0:,0:)
    integer, intent(in) :: stage,iterations
    type(spray_sdc_population_stage), intent(inout) :: work
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5),boundary_change(:)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:),inflow_state(:,:),inflow_mach(6)
    real(dp), optional, intent(in) :: inflow_temperature(6),inflow_pressure(6)
    integer, optional, intent(in) :: maximum_received
    real(dp), optional, intent(in) :: acceleration(3)
    type(amr_tree_state_3d), optional, intent(inout) :: physical_source
    type(amr_tree_state_3d), allocatable :: physical_delta
    type(spray_sdc_population_stage) :: candidate
    type(spray_level_population) :: trial(3)
    type(amr_tree_state_3d) :: change,leftover
    real(dp) :: half_interval,drift,wall(5),exterior(reactive_nvar(size(species)))
    real(dp) :: times(5),low_time(5),high_time(5),accel(3),forced(5)
    integer :: flags(4),low(4),high(4),status
    logical :: valid,agreed
    if(present(last_id)) then
      allocate(watermark);watermark=last_id
    end if
    ok=.false.;accel=0
    if(present(acceleration)) accel=acceleration
    times=[time,interval,accel]
    valid=all(ieee_is_finite(times)).and.time>=0.and.interval>0.and.iterations>=1
    valid=valid.and.all(abs(accel)<=1e6_dp)
    valid=valid.and.stage>=0.and.stage<=iterations.and.size(boundary_change)==size(exterior)
    if(stage==0) then
      valid=valid.and.work%next_stage==0
    else
      valid=valid.and.stage==work%next_stage.and.iterations==work%iterations.and.target_level==work%depth
      valid=valid.and.all(accel==work%acceleration)
      valid=valid.and.(present(physical_source).eqv.work%exports_physical)
      valid=valid.and.time==work%time.and.interval==work%interval.and.allocated(work%old_boundary)
      if(allocated(work%old_boundary)) valid=valid.and.size(work%old_boundary)==size(exterior)
    end if
    call agree();if(.not.agreed) return
    call MPI_Allreduce(times,low_time,5,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(times,high_time,5,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low_time/=high_time)) return
    flags=[stage,iterations,target_level,merge(1,0,present(physical_source))]
    call MPI_Allreduce(flags,low,4,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,4,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    half_interval=.5_dp*interval
    if(half_interval<=0) return
    if(stage>0.and.stage<iterations) then
      ! PeleC initializes new spray sources from old, then particleMK skips
      ! every intermediate SDC iteration. No particle or gas sampling here.
      source=work%old_source;remainder=work%old_remainder
      if(present(physical_source)) physical_source=work%old_physical_source
      wall_change=work%old_wall;boundary_change=work%old_boundary
      work%next_stage=stage+1;ok=.true.;return
    end if
    if(present(physical_source)) allocate(physical_delta)
    candidate=work;drift=0
    if(stage==0) then
      trial=populations;drift=interval
      call apply_force()
      if(.not.agreed) return
    else
      trial=work%midpoint
    end if
    call advance_mpi_eb_population_stage(comm,species,transport,components,options,tree,owners,geometry, &
      state,temperature,width,target_level,grow,maximum_parcels,trial,half_interval,drift,faces, &
      change,leftover,wall,exterior,valid,nodal_phi,centroids,target_volume_fraction, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,maximum_received,physical_delta,watermark)
    if(.not.valid) return
    if(stage==iterations) then
      call apply_force()
      if(.not.agreed) return
    end if
    call to_rate(change)
    if(.not.agreed) return
    call to_rate(leftover)
    if(.not.agreed) return
    if(present(physical_source)) then
      call to_rate(physical_delta)
      if(.not.agreed) return
    end if
    if(stage==0) then
      candidate%acceleration=accel;candidate%force_change=forced
      candidate%exports_physical=present(physical_source)
      if(present(physical_source)) candidate%old_physical_source=physical_delta
      candidate%midpoint=trial;candidate%old_wall=wall;candidate%old_boundary=exterior
      candidate%old_source=change;candidate%old_remainder=leftover
      candidate%time=time;candidate%interval=interval;candidate%depth=target_level;candidate%iterations=iterations
    else
      valid=all(ieee_is_finite(candidate%old_wall)).and.all(ieee_is_finite(candidate%old_boundary))
      call agree();if(.not.agreed) return
      valid=.not.any(wall>huge(1._dp)-max(candidate%old_wall,0._dp)).and. &
        .not.any(wall< -huge(1._dp)-min(candidate%old_wall,0._dp))
      valid=valid.and..not.any(exterior>huge(1._dp)-max(candidate%old_boundary,0._dp)).and. &
        .not.any(exterior< -huge(1._dp)-min(candidate%old_boundary,0._dp))
      call agree();if(.not.agreed) return
      valid=.not.any(forced>huge(1._dp)-max(candidate%force_change,0._dp)).and. &
        .not.any(forced< -huge(1._dp)-min(candidate%force_change,0._dp))
      call agree();if(.not.agreed) return
      candidate%force_change=candidate%force_change+forced
      wall=wall+candidate%old_wall;exterior=exterior+candidate%old_boundary
    end if
    candidate%next_stage=stage+1
    if(stage==iterations) then
      candidate%next_stage=0
      populations=trial
    end if
    work=candidate
    if(present(physical_source)) call move_alloc(physical_delta%fields,physical_source%fields)
    call move_alloc(change%fields,source%fields);call move_alloc(leftover%fields,remainder%fields)
    if(present(last_id)) last_id=watermark
    wall_change=wall;boundary_change=exterior;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
    ! Accelerate every copy, but count external work only for real particles.
    ! trial and candidate stay private until the entire stage succeeds.
    subroutine apply_force()
      real(dp) :: local(5),dv(3),velocity(3),mass
      integer :: population,p
      valid=interval<=1e6_dp;call agree();if(.not.agreed) return
      dv=half_interval*accel;local=0;forced=0
      do population=1,3
        valid=allocated(trial(population)%parcels)
        if(.not.valid) exit
        do p=1,size(trial(population)%parcels)
          associate(parcel=>trial(population)%parcels(p))
            if(parcel%mass<=0) cycle
            velocity=parcel%velocity;mass=parcel%mass*parcel%multiplicity
            parcel%velocity=velocity+dv
            valid=valid.and.valid_spray_parcel(parcel)
            if(population==1) then
              local(2:4)=local(2:4)+mass*dv
              local(5)=local(5)+mass*dot_product(velocity+.5_dp*dv,dv)
            end if
          end associate
        end do
        if(.not.valid) exit
      end do
      valid=valid.and.all(ieee_is_finite(local))
      call agree();if(.not.agreed) return
      call MPI_Allreduce(local,forced,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
      agreed=status==MPI_SUCCESS.and.all(ieee_is_finite(forced))
    end subroutine apply_force
    subroutine to_rate(field)
      type(amr_tree_state_3d), intent(inout) :: field
      integer :: node
      valid=.true.
      if(half_interval<1) then
        do node=1,size(owners)
          if(.not.allocated(field%fields(node)%values)) cycle
          if(any(abs(field%fields(node)%values)>huge(1._dp)*half_interval)) valid=.false.
        end do
      end if
      call agree();if(.not.agreed) return
      do node=1,size(owners)
        if(allocated(field%fields(node)%values)) field%fields(node)%values=field%fields(node)%values/half_interval
      end do
    end subroutine to_rate
  end subroutine advance_mpi_eb_population_sdc_stage

  ! One level stage, populations ordered physical / virtual / ghost.
  ! All populations sample the same frozen stage-time gas. Only physical
  ! wall/outflow ledgers belong to the physical system inventory.
  subroutine advance_mpi_eb_population_stage(comm,species,transport,components,options,tree,owners,geometry, &
      state,temperature,width,target_level,grow,maximum_parcels,populations,kick_dt,drift_dt,faces, &
      source,remainder,wall_change,boundary_change,ok,nodal_phi,centroids,target_volume_fraction, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,maximum_received,physical_source,last_id)
    integer(int64), optional, intent(inout) :: last_id
    integer(int64), allocatable :: watermark
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: owners(:),width,target_level,grow,maximum_parcels
    type(spray_level_population), intent(inout) :: populations(3)
    type(amr_tree_state_3d), intent(in) :: state,temperature,centroids,target_volume_fraction
    real(dp), intent(in) :: kick_dt,drift_dt,nodal_phi(0:,0:,0:)
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5),boundary_change(:)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:),inflow_state(:,:),inflow_mach(6)
    real(dp), optional, intent(in) :: inflow_temperature(6),inflow_pressure(6)
    integer, optional, intent(in) :: maximum_received
    type(amr_tree_state_3d), optional, intent(inout) :: physical_source
    type(amr_tree_state_3d), allocatable :: physical_delta
    type(amr_tree_state_3d) :: exported
    type(spray_options) :: population_options
    type(spray_level_population) :: candidate(3)
    type(amr_tree_state_3d) :: halo,halo_t,volume,total,leftover,part,part_leftover
    real(dp) :: lower(3),spacing(3),wall(5),physical_wall(5)
    real(dp) :: exterior(reactive_nvar(size(species))),physical_exterior(reactive_nvar(size(species)))
    character(len=32) :: particle_faces(6)
    integer :: population,status
    logical :: valid,agreed,periodic(3)
    if(present(last_id)) then
      allocate(watermark);watermark=last_id
    end if
    ok=.false.
    valid=size(boundary_change)==size(exterior)
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    do population=1,3
      valid=valid.and.allocated(populations(population)%parcels).and. &
        allocated(populations(population)%component_mass).and.allocated(populations(population)%nodes)
    end do
    if(valid.and.target_level==0) valid=size(populations(3)%parcels)==0
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,halo,halo_t,comm,valid, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_mach,inflow_temperature,inflow_pressure,volume,width,target_volume_fraction, &
      faces,inflow_state,restrict_snapshot=.false.)
    if(.not.valid) return
    lower=[geometry(1)%x_lower,geometry(1)%y_lower,geometry(1)%z_lower]
    spacing=[geometry(1)%dx,geometry(1)%dy,geometry(1)%dz]
    periodic=[faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']
    particle_faces=faces
    where(particle_faces=='inflow') particle_faces='outflow'
    where(particle_faces=='isothermal_wall') particle_faces='no_slip_wall'
    if(present(physical_source)) allocate(physical_delta)
    candidate=populations
    do population=1,3
      population_options=options
      if(population/=1.or.drift_dt<=0) then
        population_options%khrt_children=.false.
        population_options%khrt_reference_diameter_cubed=0
      end if
      ! Root has no ghosts, but use retained routing for its allocated empty set.
      call advance_mpi_amr_population_halo(comm,species,transport,components,population_options,tree,owners,halo,halo_t, &
        width,target_level,grow,maximum_parcels,candidate(population)%nodes,candidate(population)%parcels, &
        candidate(population)%component_mass,lower,spacing,periodic,kick_dt,drift_dt,particle_faces,part,part_leftover, &
        wall,exterior,valid,volume,nodal_phi,centroids,geometry,maximum_received,population<3.or.target_level==0, &
        containing_cell=.true.,physical_source=physical_delta,last_id=watermark)
      if(.not.valid) return
      if(population==1) then
        if(present(physical_source)) then
          exported=physical_delta;deallocate(physical_delta)
        end if
        total=part;leftover=part_leftover;physical_wall=wall;physical_exterior=exterior
      else
        call add_field(total,part,valid)
        if(.not.valid) return
        call add_field(leftover,part_leftover,valid)
        if(.not.valid) return
      end if
    end do
    if(present(physical_source)) call move_alloc(exported%fields,physical_source%fields)
    populations=candidate
    call move_alloc(total%fields,source%fields);call move_alloc(leftover%fields,remainder%fields)
    if(present(last_id)) last_id=watermark
    wall_change=physical_wall;boundary_change=physical_exterior;ok=.true.
  contains
    subroutine add_field(destination,addition,good)
      type(amr_tree_state_3d), intent(inout) :: destination
      type(amr_tree_state_3d), intent(in) :: addition
      logical, intent(out) :: good
      integer :: node
      logical :: safe
      safe=.true.
      do node=1,size(owners)
        if(.not.allocated(destination%fields(node)%values)) cycle
        associate(a=>destination%fields(node)%values,b=>addition%fields(node)%values)
          safe=safe.and..not.any(a>huge(1._dp)-max(b,0._dp)).and..not.any(a< -huge(1._dp)-min(b,0._dp))
        end associate
      end do
      call MPI_Allreduce(safe,good,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) good=.false.
      if(.not.good) return
      do node=1,size(owners)
        if(allocated(destination%fields(node)%values)) &
          destination%fields(node)%values=destination%fields(node)%values+addition%fields(node)%values
      end do
    end subroutine add_field
  end subroutine advance_mpi_eb_population_stage

  ! Construct complete stage-time gas halos, then advance a ghost population.
  ! Ancestor interpolation never averages down the caller's stage snapshots.
  subroutine advance_mpi_eb_ghost_stage(comm,species,transport,components,options,tree,owners,geometry,state,temperature, &
      width,target_level,grow,maximum_ghosts,nodes,ghosts,component_mass,kick_dt,drift_dt,faces, &
      source,remainder,wall_change,boundary_change,ok,nodal_phi,centroids,target_volume_fraction, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,maximum_received)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: owners(:),width,target_level,grow,maximum_ghosts
    integer, allocatable, intent(inout) :: nodes(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature,centroids,target_volume_fraction
    type(spray_parcel), allocatable, intent(inout) :: ghosts(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: kick_dt,drift_dt,nodal_phi(0:,0:,0:)
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5),boundary_change(:)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:),inflow_state(:,:),inflow_mach(6)
    real(dp), optional, intent(in) :: inflow_temperature(6),inflow_pressure(6)
    integer, optional, intent(in) :: maximum_received
    type(amr_tree_state_3d) :: halo,halo_t,volume
    real(dp) :: lower(3),spacing(3)
    logical :: valid,periodic(3)
    ok=.false.
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,halo,halo_t,comm,valid, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_mach,inflow_temperature,inflow_pressure,volume,width,target_volume_fraction, &
      faces,inflow_state,restrict_snapshot=.false.)
    if(.not.valid) return
    lower=[geometry(1)%x_lower,geometry(1)%y_lower,geometry(1)%z_lower]
    spacing=[geometry(1)%dx,geometry(1)%dy,geometry(1)%dz]
    periodic=[faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']
    call advance_mpi_amr_ghost_halo(comm,species,transport,components,options,tree,owners,halo,halo_t, &
      width,target_level,grow,maximum_ghosts,nodes,ghosts,component_mass,lower,spacing,periodic,kick_dt,drift_dt,faces, &
      source,remainder,wall_change,boundary_change,ok,volume,nodal_phi,centroids,geometry,maximum_received)
  end subroutine advance_mpi_eb_ghost_stage

  ! One ghost-population stage: kick, drift, same-level reaction sum, then reassign.
  ! Ghost halo remainders must not be added to coarse physical-particle sources.
  ! All outputs remain private until the final population exchange succeeds.
  subroutine advance_mpi_amr_ghost_halo(comm,species,transport,components,options,tree,owners,gas_halo,temperature_halo, &
      width,target_level,grow,maximum_ghosts,nodes,ghosts,component_mass,lower,spacing,periodic,kick_dt,drift_dt,faces, &
      source,remainder,wall_change,boundary_change,ok,fluid_fraction,nodal_phi,centroids,eb_geometry,maximum_received)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),width,target_level,grow,maximum_ghosts
    integer, allocatable, intent(inout) :: nodes(:)
    type(amr_tree_state_3d), intent(in) :: gas_halo,temperature_halo
    type(spray_parcel), allocatable, intent(inout) :: ghosts(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),kick_dt,drift_dt
    logical, intent(in) :: periodic(3)
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5),boundary_change(:)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: fluid_fraction,centroids
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(eb_geometry_3d), optional, intent(in) :: eb_geometry(:)
    integer, optional, intent(in) :: maximum_received
    call advance_mpi_amr_population_halo(comm,species,transport,components,options,tree,owners, &
      gas_halo,temperature_halo,width,target_level,grow,maximum_ghosts,nodes,ghosts,component_mass, &
      lower,spacing,periodic,kick_dt,drift_dt,faces,source,remainder,wall_change,boundary_change,ok, &
      fluid_fraction,nodal_phi,centroids,eb_geometry,maximum_received,.false.)
  end subroutine advance_mpi_amr_ghost_halo

  ! Physical stage: reject loss outside the level grow region after kick/drift.
  subroutine advance_mpi_amr_level_halo(comm,species,transport,components,options,tree,owners,gas_halo,temperature_halo, &
      width,target_level,grow,maximum_parcels,nodes,parcels,component_mass,lower,spacing,periodic,kick_dt,drift_dt,faces, &
      source,remainder,wall_change,boundary_change,ok,fluid_fraction,nodal_phi,centroids,eb_geometry,maximum_received)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),width,target_level,grow,maximum_parcels
    integer, allocatable, intent(inout) :: nodes(:)
    type(amr_tree_state_3d), intent(in) :: gas_halo,temperature_halo
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),kick_dt,drift_dt
    logical, intent(in) :: periodic(3)
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5),boundary_change(:)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: fluid_fraction,centroids
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(eb_geometry_3d), optional, intent(in) :: eb_geometry(:)
    integer, optional, intent(in) :: maximum_received
    call advance_mpi_amr_population_halo(comm,species,transport,components,options,tree,owners, &
      gas_halo,temperature_halo,width,target_level,grow,maximum_parcels,nodes,parcels,component_mass, &
      lower,spacing,periodic,kick_dt,drift_dt,faces,source,remainder,wall_change,boundary_change,ok, &
      fluid_fraction,nodal_phi,centroids,eb_geometry,maximum_received,.true.)
  end subroutine advance_mpi_amr_level_halo

  subroutine advance_mpi_amr_population_halo(comm,species,transport,components,options,tree,owners,gas_halo,temperature_halo, &
      width,target_level,grow,maximum_ghosts,nodes,ghosts,component_mass,lower,spacing,periodic,kick_dt,drift_dt,faces, &
      source,remainder,wall_change,boundary_change,ok,fluid_fraction,nodal_phi,centroids,eb_geometry, &
      maximum_received,keep_physical,containing_cell,physical_source,last_id)
    integer(int64), optional, intent(inout) :: last_id
    integer(int64), allocatable :: watermark
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),width,target_level,grow,maximum_ghosts
    integer, allocatable, intent(inout) :: nodes(:)
    type(amr_tree_state_3d), intent(in) :: gas_halo,temperature_halo
    type(spray_parcel), allocatable, intent(inout) :: ghosts(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),kick_dt,drift_dt
    logical, intent(in) :: periodic(3)
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5),boundary_change(:)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: fluid_fraction,centroids
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(eb_geometry_3d), optional, intent(in) :: eb_geometry(:)
    integer, optional, intent(in) :: maximum_received
    type(amr_tree_state_3d), optional, intent(inout) :: physical_source
    type(amr_tree_state_3d), allocatable :: physical_delta
    logical, intent(in) :: keep_physical
    logical, optional, intent(in) :: containing_cell
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: masses(:,:)
    integer, allocatable :: locations(:)
    type(amr_tree_state_3d) :: raw,routed,leftover
    real(dp) :: wall(5),exterior(reactive_nvar(size(species)))
    integer :: status
    logical :: valid,agreed
    if(present(last_id)) then
      allocate(watermark);watermark=last_id
    end if
    ok=.false.
    valid=allocated(ghosts).and.allocated(component_mass).and.allocated(nodes)
    valid=valid.and.size(boundary_change)==size(exterior)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    if(present(physical_source)) allocate(physical_delta)
    candidate=ghosts;masses=component_mass;locations=nodes
    call kick_mpi_amr_mixture_halo(comm,species,transport,components,options,tree,owners, &
      gas_halo,temperature_halo,width,locations,candidate,masses,lower,spacing,periodic,kick_dt,raw,valid, &
      fluid_fraction,nodal_phi,centroids,drift_dt,faces,wall,exterior,target_level,routed,leftover,.false., &
      eb_geometry,maximum_received,containing_cell,physical_delta,watermark)
    if(.not.valid) return
    if(keep_physical) then
      call redistribute_mpi_amr_level_particles(comm,tree,owners,candidate,masses,locations,lower,spacing, &
        periodic,target_level,grow,maximum_ghosts,valid,nodal_phi)
    else
    call redistribute_mpi_amr_ghost_particles(comm,tree,owners,candidate,masses,locations,lower,spacing, &
      periodic,target_level,grow,maximum_ghosts,valid,nodal_phi)
    end if
    if(.not.valid) return
    call move_alloc(candidate,ghosts);call move_alloc(masses,component_mass);call move_alloc(locations,nodes)
    call move_alloc(routed%fields,source%fields);call move_alloc(leftover%fields,remainder%fields)
    if(present(physical_source)) physical_source=physical_delta
    if(present(last_id)) last_id=watermark
    wall_change=wall;boundary_change=exterior;ok=.true.
  end subroutine advance_mpi_amr_population_halo

  ! A single drag/thermal/evaporation kick on a caller-owned level population.
  ! Gas halos are frozen at the supplied stage time. Optional drift follows the
  ! kick. Optional KHRT children need last_id; source nodes describe the input
  ! parents. The caller must rebuild particle ownership after child creation.
  ! Species/model tables and options are identical across ranks by caller contract.
  ! Gas source and particle/component updates publish together after MPI success.
  subroutine kick_mpi_amr_mixture_halo(comm,species,transport,components,options,tree,owners, &
      gas_halo,temperature_halo,width,nodes,parcels,component_mass,lower,spacing,periodic,dt,source_halo,ok, &
      fluid_fraction,nodal_phi,centroids,drift_interval,boundary_faces,wall_change,boundary_change, &
      target_level,level_source,source_remainder,route_to_coarse,eb_geometry,maximum_received, &
      containing_cell,physical_source,last_id,wall_films,film_mass,film_settings,model_rejected)
    integer(int64), optional, intent(inout) :: last_id
    integer(int64) :: watermark
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    logical, optional, intent(out) :: model_rejected
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:),pressures(:),impact_pressure(:)
    type(spray_eb_wall_options) :: impact_settings
    integer :: film_flag,film_low,film_high,parent,lo,hi,mid
    logical :: rejected
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    real(dp) :: mw(size(components)),impact_exterior(5+size(components))
    type(spray_options), intent(in) :: options
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),width,nodes(:)
    type(amr_tree_state_3d), intent(in) :: gas_halo,temperature_halo
    type(amr_tree_state_3d), intent(inout) :: source_halo
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),dt
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: fluid_fraction,centroids
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:),drift_interval
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(inout) :: wall_change(5),boundary_change(:)
    real(dp) :: drift_dt,drift_low,drift_high,wall(5),global_wall(5)
    real(dp) :: exterior(reactive_nvar(size(species))),global_exterior(reactive_nvar(size(species)))
    real(dp) :: exterior5(5),upper(3),empty_phi(0:1,0:1,0:1),rho,cp,energy
    real(dp), allocatable :: energies(:)
    integer :: flags,low_flags,high_flags,hits,c
    character(len=32) :: faces(6),reference_faces(6)

    type(spray_parcel), allocatable :: candidate(:)
    integer, optional, intent(in) :: target_level,maximum_received
    type(amr_tree_state_3d), optional, intent(inout) :: level_source,source_remainder
    logical, optional, intent(in) :: route_to_coarse
    type(eb_geometry_3d), optional, intent(in) :: eb_geometry(:)
    type(amr_tree_state_3d), optional, intent(inout) :: physical_source
    type(amr_tree_state_3d) :: physical_delta
    integer :: export_flag,export_low,export_high
    logical, optional, intent(in) :: containing_cell
    type(amr_tree_state_3d) :: source,routed,remainder
    integer :: transfer_flags,transfer_low,transfer_high,node,rank,n(3)
    real(dp), allocatable :: masses(:,:),positions(:,:),sampled(:,:),sampled_t(:,:),change(:,:)
    real(dp) :: primitive(reactive_nvar(size(species))),tg,sound,low_dt,high_dt
    real(dp), allocatable :: gas_properties(:,:),reynolds(:),child_energy(:)
    logical, allocatable :: eligible(:)
    type(spray_parcel) :: preview(1)
    real(dp) :: preview_wall(5),preview_exterior(5)
    logical :: contact(1)
    type(spray_liquid) :: effective
    type(spray_film) :: film
    real(dp) :: old_mass,fractions(size(components)),rates(size(components)),hv(size(components))
    integer :: p,nv,nc,status
    logical :: valid,agreed
    ok=.false.;nv=reactive_nvar(size(species));nc=size(components)
    if(present(model_rejected)) model_rejected=.false.
    film_flag=merge(1,0,present(wall_films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))
    call MPI_Allreduce(film_flag,film_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,film_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.film_low/=film_high) return
    if(film_flag/=0.and.film_flag/=7) return
    if(present(wall_films)) then
      valid=allocated(wall_films).and.allocated(film_mass).and.present(last_id).and. &
        present(nodal_phi).and.present(drift_interval)
      call agree();if(.not.agreed) return
      valid=valid_wall_options(film_settings,[character(len=16):: &
        'slip_wall','slip_wall','slip_wall','slip_wall','slip_wall','slip_wall'])
      valid=valid.and.all(film_settings%temperature==film_settings%temperature(1))
      call agree();if(.not.agreed) return
      patches=wall_films;patch_mass=film_mass
      impact_settings%temperature=film_settings%temperature(1)
      impact_settings%contact_angle=film_settings%contact_angle
      impact_settings%maximum_films=film_settings%maximum_films
    end if
    valid=allocated(parcels).and.allocated(component_mass).and.valid_spray_options(options)
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    valid=valid.and.nc>=1.and.nc<=16.and.size(transport)==size(species)
    valid=valid.and.ieee_is_finite(dt)
    if(valid) valid=dt>=0.and.dt<=1e6_dp
    call agree();if(.not.agreed) return
    do c=1,nc
      valid=valid.and.valid_spray_liquid(components(c))
      valid=valid.and.components(c)%vapor_index>=1.and.components(c)%vapor_index<=size(species)
    end do
    valid=valid.and.size(nodes)==size(parcels).and.all(shape(component_mass)==[nc,size(parcels)])
    call agree();if(.not.agreed) return
    do p=1,size(parcels)
      valid=valid_component_inventory(parcels(p),component_mass(:,p))
      if(.not.valid) exit
    end do
    call agree();if(.not.agreed) return
    call MPI_Allreduce(dt,low_dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,high_dt,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low_dt/=high_dt) return
    export_flag=merge(1,0,options%khrt_children)+2*merge(1,0,present(last_id))
    call MPI_Allreduce(export_flag,export_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(export_flag,export_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.export_low/=export_high) return
    watermark=0
    if(present(last_id)) watermark=last_id
    flags=merge(1,0,present(drift_interval))+2*merge(1,0,present(boundary_faces))+ &
      4*merge(1,0,present(wall_change))+8*merge(1,0,present(boundary_change))
    call MPI_Allreduce(flags,low_flags,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high_flags,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low_flags/=high_flags) return
    if(flags/=0.and.flags/=15) return
    if(present(drift_interval)) then
      drift_dt=drift_interval;faces=boundary_faces
      valid=ieee_is_finite(drift_dt).and.size(boundary_change)==nv
      if(valid) valid=drift_dt>=0.and.drift_dt<=1e6_dp
      valid=valid.and.valid_spray_boundary_faces(faces)
      valid=valid.and.all(([faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']).eqv.periodic)
      call agree();if(.not.agreed) return
      call MPI_Allreduce(drift_dt,drift_low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(drift_dt,drift_high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.drift_low/=drift_high) return
      reference_faces=faces
      call MPI_Bcast(reference_faces,32*6,MPI_CHARACTER,0,comm,status)
      if(status/=MPI_SUCCESS) return
      valid=all(faces==reference_faces)
      call agree();if(.not.agreed) return
    end if
    transfer_flags=merge(1,0,present(target_level))+2*merge(1,0,present(level_source))+ &
      4*merge(1,0,present(source_remainder))
    call MPI_Allreduce(transfer_flags,transfer_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(transfer_flags,transfer_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.transfer_low/=transfer_high) return
    if(transfer_flags/=0.and.transfer_flags/=7) return
    valid=transfer_flags==7.or..not.(present(route_to_coarse).or.present(eb_geometry).or.present(maximum_received))
    if(present(target_level)) valid=valid.and.(present(eb_geometry).eqv.present(fluid_fraction))
    call agree();if(.not.agreed) return
    export_flag=merge(1,0,present(physical_source))
    call MPI_Allreduce(export_flag,export_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(export_flag,export_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.export_low/=export_high) return
    allocate(positions(3,size(parcels)),change(nv,size(parcels)))
    do p=1,size(parcels)
      positions(:,p)=parcels(p)%position
    end do
    call sample_mpi_amr_particle_halo(comm,tree,owners,gas_halo,nv,width,nodes,positions,lower,spacing, &
      periodic,sampled,valid,fluid_fraction,nodal_phi,centroids)
    if(.not.valid) return
    call sample_mpi_amr_particle_halo(comm,tree,owners,temperature_halo,1,width,nodes,positions,lower,spacing, &
      periodic,sampled_t,valid,fluid_fraction,nodal_phi,centroids)
    if(.not.valid) return
    if(present(target_level)) then
      valid=target_level>=0.and.target_level<=maxval(tree%nodes%level)
      do p=1,size(nodes)
        valid=valid.and.tree%nodes(nodes(p))%level==target_level
      end do
      call agree();if(.not.agreed) return
    end if
    candidate=parcels;masses=component_mass;change=0
    allocate(gas_properties(4,size(parcels)),reynolds(size(parcels)))
    gas_properties=0;reynolds=0
    allocate(pressures(size(parcels)),source=1._dp)
    do p=1,size(parcels)
      if(parcels(p)%mass<=0) cycle
      call reactive_conserved_to_primitive(species,sampled(:,p),sampled_t(1,p),primitive,tg,sound,valid)
      if(.not.valid) exit
      pressures(p)=primitive(5)
      call advance_mixture_frozen_gas(species,transport,primitive,tg,candidate(p),components,masses(:,p), &
        options,dt,change(:,p),valid)
      if(.not.valid) exit
      if(candidate(p)%mass>0.and.options%breakup_model/='none') then
        old_mass=candidate(p)%mass;fractions=masses(:,p)/old_mass
        call liquid_mixture_properties(components,fractions,candidate(p)%temperature,rho,cp,energy,valid)
        if(.not.valid) exit
        effective=components(1);effective%density=rho;effective%density_model='constant'
        effective%density_coefficients=0;effective%density_temperature_range=0
        effective%heat_capacity=cp;effective%heat_capacity_model='constant'
        effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
        effective%reference_temperature=candidate(p)%temperature;effective%reference_energy=energy
        effective%boiling_model='reject';effective%critical_temperature=0;effective%normal_boiling_temperature=0
        call build_mixture_film(species,transport,primitive,tg,candidate(p),components,fractions,options, &
          film,rates,hv,valid)
        if(.not.valid) exit
        if(options%khrt_children) then
          gas_properties(:,p)=primitive(:4)
          reynolds(p)=film%density*norm2(primitive(2:4)-candidate(p)%velocity)* &
            parcel_diameter(candidate(p),effective)/film%viscosity
        else
          call advance_spray_breakup(candidate(p),effective,options,primitive(1),primitive(2:4), &
            film%density,film%viscosity,dt,valid)
          if(.not.valid) exit
          masses(:,p)=masses(:,p)*(candidate(p)%mass/old_mass)
        end if
      end if
      candidate(p)%position=parcels(p)%position
    end do
    call agree();if(.not.agreed) return
    if(options%khrt_children) then
      allocate(eligible(size(candidate)));eligible=.true.
      ! A reflecting or EB contact takes precedence over release, as in the
      ! common-step coupling. Preview the same segmented flight as the drift.
      if(present(drift_interval)) then
        upper=lower+spacing*real(tree%nodes(1)%shape,dp);empty_phi=1
        do p=1,size(candidate)
          if(candidate(p)%mass<=0) cycle
          preview(1)=candidate(p);contact=.false.;preview_wall=0;preview_exterior=0;hits=0
          if(present(nodal_phi)) then
            call drift_spray_level_set(nodal_phi,lower,upper,drift_dt,32,preview,preview_wall,hits,valid, &
              faces,[0._dp],preview_exterior,contact)
          else
            call drift_spray_level_set(empty_phi,lower,upper,drift_dt,32,preview,preview_wall,hits,valid, &
              faces,[0._dp],preview_exterior,contact)
          end if
          if(.not.valid) exit
          eligible(p)=.not.contact(1)
        end do
        call agree();if(.not.agreed) return
      end if
      call advance_mpi_khrt_shedding(comm,components,options,candidate,masses,gas_properties(1,:), &
        gas_properties(2:4,:),reynolds,dt,watermark,child_energy,valid,eligible)
      if(.not.valid) return
      valid=all(ieee_is_finite(child_energy))
      if(valid) valid=.not.any(change(5,:)>huge(1._dp)-max(child_energy,0._dp)).and. &
        .not.any(change(5,:)< -huge(1._dp)-min(child_energy,0._dp))
      call agree();if(.not.agreed) return
      change(5,:)=change(5,:)+child_energy
    end if
    if(present(wall_films)) then
      ! Shedding children inherit the sampled pressure of their pre-drift parent,
      ! matching the common-step impact preparation. Sources retain parent sites.
      allocate(impact_pressure(size(candidate)))
      do p=1,size(candidate)
        parent=p
        if(p>size(parcels)) then
          parent=0
          lo=1;hi=size(parcels)
          do while(lo<=hi)
            mid=lo+(hi-lo)/2
            if(parcels(mid)%id==candidate(p)%parent_id) then
              parent=mid;exit
            else if(parcels(mid)%id<candidate(p)%parent_id) then
              lo=mid+1
            else
              hi=mid-1
            end if
          end do
          if(parent==0) then
            valid=.false.;exit
          end if
        end if
        impact_pressure(p)=pressures(parent)
      end do
      call agree();if(.not.agreed) return
      do c=1,nc
        mw(c)=species(components(c)%vapor_index)%molecular_weight
      end do
      upper=lower+spacing*real(tree%nodes(1)%shape,dp)
      global_wall=0;impact_exterior=0;hits=0
      call advance_mpi_level_set_impacts(comm,components,options,impact_settings,nodal_phi,mw,impact_pressure, &
        lower,upper,drift_dt,candidate,masses,patches,patch_mass,watermark,global_wall,hits,valid,rejected, &
        impact_exterior,faces)
      if(rejected.and.present(model_rejected)) model_rejected=.true.
      if(.not.valid) return
      ! Impact API reports inventory received by the wall, already MPI summed.
      global_wall=-global_wall;global_exterior=0;global_exterior(:5)=impact_exterior(:5)
      do c=1,nc
        global_exterior(5+components(c)%vapor_index)=global_exterior(5+components(c)%vapor_index)+impact_exterior(5+c)
      end do
    else if(present(drift_interval)) then
      allocate(energies(size(candidate)));energies=0
      do p=1,size(candidate)
        if(candidate(p)%mass<=0) cycle
        call liquid_mixture_properties(components,masses(:,p)/candidate(p)%mass, &
          candidate(p)%temperature,rho,cp,energy,valid)
        if(.not.valid) exit
        energies(p)=energy
      end do
      call agree();if(.not.agreed) return
      upper=lower+spacing*real(tree%nodes(1)%shape,dp)
      wall=0;exterior=0;exterior5=0;hits=0;empty_phi=1
      if(present(nodal_phi)) then
        call drift_spray_level_set(nodal_phi,lower,upper,drift_dt,32,candidate,wall,hits,valid, &
          faces,energies,exterior5)
      else
        call drift_spray_level_set(empty_phi,lower,upper,drift_dt,32,candidate,wall,hits,valid, &
          faces,energies,exterior5)
      end if
      call agree();if(.not.agreed) return
      exterior(:5)=exterior5
      do p=1,size(candidate)
        if(candidate(p)%mass>0) cycle
        do c=1,nc
          exterior(5+components(c)%vapor_index)=exterior(5+components(c)%vapor_index)- &
            masses(c,p)*candidate(p)%multiplicity
        end do
        masses(:,p)=0
      end do
      call MPI_Allreduce(wall,global_wall,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(exterior,global_exterior,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
      if(status/=MPI_SUCCESS) return
      if(.not.all(ieee_is_finite(global_wall)).or..not.all(ieee_is_finite(global_exterior))) return
    end if
    call deposit_mpi_amr_particle_halo(comm,tree,owners,gas_halo,nv,width,nodes,positions,change, &
      lower,spacing,periodic,source,valid,fluid_fraction,nodal_phi,centroids,containing_cell)
    if(.not.valid) return
    if(present(target_level)) then
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
      allocate(routed%fields(size(tree%nodes)))
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(routed%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
      end do
      call sum_mpi_amr_particle_source_halo(comm,tree,owners,source,nv,width,target_level,periodic, &
        routed,remainder,valid,maximum_received,route_to_coarse,eb_geometry,fluid_fraction)
      if(.not.valid) return
    end if
    if(present(physical_source)) then
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
      allocate(physical_delta%fields(size(owners)))
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(physical_delta%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
      end do
      call deposit_mpi_amr_tree_spray_field(comm,tree,owners,physical_delta,nv,positions,change,lower,spacing, &
        periodic,valid,maximum_received=maximum_received,eb_geometry=eb_geometry,nodal_phi=nodal_phi,containing_cell=.true.)
      if(.not.valid) return
    end if
    if(present(drift_interval)) then
      wall_change=global_wall;boundary_change=global_exterior
    end if
    if(present(target_level)) then
      call move_alloc(routed%fields,level_source%fields)
      call move_alloc(remainder%fields,source_remainder%fields)
    end if
    if(present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    if(present(last_id)) last_id=watermark
    call move_alloc(candidate,parcels);call move_alloc(masses,component_mass)
    if(present(physical_source)) call move_alloc(physical_delta%fields,physical_source%fields)
    call move_alloc(source%fields,source_halo%fields);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
  end subroutine kick_mpi_amr_mixture_halo

  ! Regular or nodal EB AMR: common particle substeps and whole-interval rollback.
  ! Species, transport, liquids and options must be identical on every rank.
  ! Optional gas_change sums accepted deposition increments per fluid volume;
  ! rejected trials and wall/outflow inventories do not enter this source.
  subroutine advance_mpi_amr_mixture_coupling(comm,species,transport,components,options,state,temperature, &
      tree,owners,parcels,component_mass,lower,spacing,interval,substeps,ok,last_id,boundary_faces,boundary_change, &
      wall_films,film_mass,film_settings,eb_geometry,nodal_phi,model_rejected,eb_outflow_change,gas_change, &
      motion_mode,maximum_level)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(amr_tree_state_3d), optional, intent(inout) :: gas_change
    type(amr_tree_state_3d) :: integrated_change
    type(amr_tree_state_3d), allocatable :: trial_change
    integer :: export_flag,export_low,export_high
    integer, optional, intent(in) :: maximum_level
    character(len=*), optional, intent(in) :: motion_mode
    character(len=16) :: motion
    integer :: motion_code,motion_low,motion_high
    real(dp) :: kick_fraction,kick_dt,motion_dt
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),interval
    integer, intent(out) :: substeps
    logical, intent(out) :: ok
    integer(int64), optional, intent(inout) :: last_id
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(inout) :: boundary_change(:)
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    type(eb_geometry_3d), optional, intent(in) :: eb_geometry(:)
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    real(dp), optional, intent(inout) :: eb_outflow_change(:)
    real(dp), allocatable :: impact_exterior(:),liquid_energies(:)
    character(len=16) :: outflow_faces(6)
    real(dp) :: exterior(reactive_nvar(size(species))),finish(3),delta(3),flight_fraction
    integer :: axis
    logical, optional, intent(out) :: model_rejected
    type(spray_eb_wall_options) :: eb_walls
    logical :: rejected
    integer :: eb_hits
    type(spray_wall_film), allocatable :: patches(:),trial_patches(:)
    real(dp), allocatable :: patch_mass(:,:),trial_patch_mass(:,:),impact_previous(:,:),impact_pressure(:)
    real(dp) :: mw(size(components)),impact_wall(5),wall_ledger(5)
    logical, allocatable :: hit(:),eligible(:)
    type(spray_parcel) :: preview(1)
    logical :: preview_contact(1)
    real(dp) :: preview_change(5),preview_exterior(5)
    integer :: preview_hits
    logical :: crossing
    real(dp) :: crossing_fraction,crossing_position(3),crossing_normal(3)
    integer :: crossing_cell(3)
    integer :: film_flag
    real(dp) :: ledger(reactive_nvar(size(species))),trial_ledger(reactive_nvar(size(species)))
    real(dp) :: global_ledger(reactive_nvar(size(species))),exchange(5)
    character(len=32) :: faces(6),reference_faces(6)
    integer(int64) :: accepted_id,trial_id
    real(dp), allocatable :: reynolds(:),child_energy(:)
    integer :: child_flag,low_flag,high_flag
    type(spray_parcel), allocatable :: drops(:),trial(:),survivors(:)
    type(spray_liquid) :: effective
    type(spray_film) :: film
    real(dp) :: old_mass
    type(amr_tree_state_3d) :: gas,temp,next_gas,next_temp,field
    real(dp), allocatable :: masses(:,:),next_mass(:,:),kept_mass(:,:),positions(:,:),samples(:,:)
    real(dp), allocatable :: primitive(:,:),tg(:),change(:,:)
    real(dp) :: elapsed,dt,local_dt,min_interval,max_interval,sound,invtau,heat,rho,cp,e,tdot,speed
    real(dp) :: y(size(components)),rates(size(components)),hv(size(components)),length(3)
    integer :: p,n,np,nvar,nc,ierr,attempt,count_steps,kept,min_nc,max_nc,node,shape3(3)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: minimum_spacing
    logical :: periodic(3)
    logical :: valid,agreed,accepted
    ok=.false.;substeps=0;faces='periodic';ledger=0;wall_ledger=0
    exterior=0
    if(present(model_rejected)) model_rejected=.false.
    export_flag=merge(1,0,present(gas_change))
    call MPI_Allreduce(export_flag,export_low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(export_flag,export_high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.export_low/=export_high) return
    motion='coupled';kick_fraction=1
    if(present(motion_mode)) motion=motion_mode
    select case(motion)
    case('coupled')
      motion_code=0
    case('kick_drift')
      motion_code=1;kick_fraction=.5_dp
    case('kick')
      motion_code=2;kick_fraction=.5_dp
    case default
      motion_code=-1
    end select
    call MPI_Allreduce(motion_code,motion_low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(motion_code,motion_high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.motion_low<0.or.motion_low/=motion_high) return
    if(present(eb_geometry)) faces="slip_wall"
    if (present(boundary_faces)) faces=boundary_faces
    valid=allocated(parcels).and.allocated(component_mass).and.valid_spray_options(options)
    valid=valid.and.(present(eb_geometry).eqv.present(nodal_phi))
    if(present(eb_geometry).and.present(boundary_faces)) then
      valid=valid.and.all(faces=='periodic'.or.faces=='outflow'.or.faces=='slip_wall'.or.faces=='no_slip_wall')
      valid=valid.and.present(eb_outflow_change)
    end if
    valid=valid.and.ieee_is_finite(interval)
    valid=valid.and.interval>=0.and.interval<=1e6_dp
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    valid=valid.and.size(components)>=1.and.size(components)<=16.and.valid_spray_boundary_faces(faces)
    if (present(boundary_faces)) then
      do n=1,6
        valid=valid.and.len_trim(boundary_faces(n))<=len(faces(n))
      end do
    end if
    if (present(boundary_change)) valid=valid.and.size(boundary_change)==size(ledger)
    if(present(eb_outflow_change)) then
      valid=valid.and.present(eb_geometry)
      valid=valid.and.size(eb_outflow_change)==size(ledger)
    end if
    valid=valid.and.(present(wall_films).eqv.present(film_mass))
    valid=valid.and.(present(wall_films).eqv.present(film_settings))
    if (present(wall_films)) valid=valid.and.allocated(wall_films).and.present(last_id).and. &
      (present(boundary_faces).or.present(eb_geometry))
    if (present(film_mass)) valid=valid.and.allocated(film_mass)
    if (present(film_settings)) then
      if(present(eb_geometry)) then
        outflow_faces='slip_wall'
        valid=valid.and.valid_wall_options(film_settings,outflow_faces)
      else
        valid=valid.and.valid_wall_options(film_settings,faces)
      end if
      if (present(eb_geometry)) then
        valid=valid.and.all(film_settings%temperature==film_settings%temperature(1))
        eb_walls%temperature=film_settings%temperature(1)
        eb_walls%contact_angle=film_settings%contact_angle
        eb_walls%maximum_films=film_settings%maximum_films
      end if
    end if
    call agree()
    if (.not.agreed) return
    film_flag=merge(1,0,present(wall_films))+2*merge(1,0,present(eb_geometry))+4*merge(1,0,present(nodal_phi))
    film_flag=film_flag+8*merge(1,0,present(eb_outflow_change))
    call MPI_Allreduce(film_flag,low_flag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,high_flag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_flag/=high_flag) return
    if(present(eb_outflow_change)) allocate(impact_exterior(5+size(components)))
    reference_faces=faces
    call MPI_Bcast(reference_faces,32*6,MPI_CHARACTER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(faces==reference_faces)
    call agree()
    if (.not.agreed) return
    periodic=[faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']
    child_flag=merge(1,0,options%khrt_children)
    call MPI_Allreduce(child_flag,low_flag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(child_flag,high_flag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_flag/=high_flag) return
    accepted_id=0
    if (present(last_id)) accepted_id=last_id
    np=size(parcels);nc=size(components);nvar=reactive_nvar(size(species))
    valid=size(component_mass,1)==nc.and.size(component_mass,2)==np.and.np<=125000
    valid=valid.and.size(transport)==size(species)
    call MPI_Allreduce(interval,min_interval,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(interval,max_interval,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nc,min_nc,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nc,max_nc,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=valid.and.min_interval==max_interval.and.min_nc==max_nc
    call agree()
    if (.not.agreed) return
    do p=1,np
      valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p))
      if (p>1) valid=valid.and.parcels(p)%id>parcels(p-1)%id
    end do
    do n=1,nc
      valid=valid.and.valid_spray_liquid(components(n)).and.components(n)%vapor_index<=size(species)
    end do
    call agree()
    if (.not.agreed) return
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    if (present(wall_films)) then
      patches=wall_films;patch_mass=film_mass
      call migrate_mpi_amr_wall_films(comm,tree,owners,patches,patch_mass,lower,spacing, &
        film_settings%maximum_films,valid,nodal_phi,maximum_level)
      if (.not.valid) return
      do n=1,nc
        mw(n)=species(components(n)%vapor_index)%molecular_weight
      end do
    end if
    allocate(positions(3,np),change(nvar,np),primitive(nvar,np),tg(np))
    do p=1,np
      positions(:,p)=drops(p)%position
    end do
    change=0
    ! This collective also validates geometry and the original gas on every rank.
    call deposit_mpi_amr_tree_spray_gas(comm,species,tree,owners,gas,temp,positions,change, &
      lower,spacing,periodic,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi,maximum_level=maximum_level)
    if (.not.valid) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree()
    if (.not.agreed) return
    if(present(gas_change)) then
      integrated_change=gas
      do node=1,size(owners)
        if(.not.allocated(integrated_change%fields(node)%values)) cycle
        integrated_change%fields(node)%values=0
      end do
      allocate(trial_change)
    end if
    length=spacing*real(tree%nodes(1)%shape,dp)
    node=size(domain,2)
    if(present(maximum_level)) node=maximum_level+1
    minimum_spacing=minval(length/real(domain(:,node),dp))
    elapsed=0;count_steps=0
    do while(elapsed<interval)
      valid=count_steps<options%maximum_substeps
      call agree()
      if (.not.agreed) return
      np=size(drops)
      deallocate(positions,change,primitive,tg)
      if (allocated(reynolds)) deallocate(reynolds)
      allocate(positions(3,np),change(nvar,np),primitive(nvar,np),tg(np),reynolds(np))
      primitive=0;tg=0
      ! Sample the temperature guess alongside conserved variables in one exchange.
      if (allocated(field%fields)) deallocate(field%fields)
      allocate(field%fields(size(tree%nodes)))
      do node=1,size(tree%nodes)
        if (.not.allocated(gas%fields(node)%values)) cycle
        shape3=tree%nodes(node)%shape
        allocate(field%fields(node)%values(nvar+1,shape3(1),shape3(2),shape3(3)))
        field%fields(node)%values(:nvar,:,:,:)=gas%fields(node)%values
        field%fields(node)%values(nvar+1,:,:,:)=temp%fields(node)%values(1,:,:,:)
      end do
      do p=1,np
        positions(:,p)=drops(p)%position
      end do
      call sample_mpi_amr_tree_spray_field(comm,tree,owners,field,nvar+1,positions, &
        lower,spacing,periodic,samples,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi,maximum_level=maximum_level)
      if (.not.valid) return
      local_dt=interval-elapsed
      do p=1,np
        if (drops(p)%mass<=0) cycle
        call reactive_conserved_to_primitive(species,samples(:nvar,p),samples(nvar+1,p),primitive(:,p),tg(p),sound,valid)
        if (.not.valid) exit
        y=masses(:,p)/drops(p)%mass
        call mixture_transfer_rates(species,transport,primitive(:,p),tg(p),drops(p),components,y,options, &
          invtau,heat,rates,hv,valid)
        if (.not.valid) exit
        call liquid_mixture_properties(components,y,drops(p)%temperature,rho,cp,e,valid)
        if (.not.valid) exit
        if (invtau>0) local_dt=min(local_dt,options%max_fractional_change/invtau/kick_fraction)
        tdot=heat*(tg(p)-drops(p)%temperature)
        do n=1,nc
          if (rates(n)>0.and.y(n)>mixture_trace_fraction) &
            local_dt=min(local_dt,options%max_fractional_change*masses(n,p)/rates(n)/kick_fraction)
          tdot=tdot-rates(n)*(hv(n)-liquid_energy(components(n),drops(p)%temperature))
        end do
        tdot=tdot/(drops(p)%mass*cp)
        if (abs(tdot)>0) local_dt=min(local_dt, &
          min(2.0_dp,options%max_fractional_change*drops(p)%temperature)/abs(tdot)/kick_fraction)
        speed=max(maxval(abs(drops(p)%velocity)),maxval(abs(primitive(2:4,p))))
        if (speed>0.and.motion/='kick') local_dt=min(local_dt,options%particle_cfl*minimum_spacing/speed)
      end do
      call agree()
      if (.not.agreed) return
      call MPI_Allreduce(local_dt,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      accepted=.false.
      do attempt=1,30
        if (dt<=0.or.elapsed+dt<=elapsed) exit
        kick_dt=kick_fraction*dt;motion_dt=dt
        if(motion=='kick') motion_dt=0
        trial=drops;next_mass=masses;change=0;valid=.true.;trial_id=accepted_id;reynolds=0;trial_ledger=0
        do p=1,np
          if (drops(p)%mass<=0) cycle
          call advance_mixture_frozen_gas(species,transport,primitive(:,p),tg(p),trial(p),components, &
            next_mass(:,p),options,kick_dt,change(:,p),valid)
          if (.not.valid) exit
          if (trial(p)%mass>0.and.options%breakup_model/='none') then
            old_mass=trial(p)%mass;y=next_mass(:,p)/old_mass
            call liquid_mixture_properties(components,y,trial(p)%temperature,rho,cp,e,valid)
            if (.not.valid) exit
            effective=components(1);effective%density=rho;effective%density_model='constant'
            effective%density_coefficients=0;effective%density_temperature_range=0
            effective%heat_capacity=cp;effective%heat_capacity_model='constant'
            effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
            effective%reference_temperature=trial(p)%temperature;effective%reference_energy=e
            effective%boiling_model='reject';effective%critical_temperature=0
            effective%normal_boiling_temperature=0
            call build_mixture_film(species,transport,primitive(:,p),tg(p),trial(p),components,y,options, &
              film,rates,hv,valid)
            if (.not.valid) exit
            if (options%khrt_children) then
              reynolds(p)=film%density*norm2(primitive(2:4,p)-trial(p)%velocity)* &
                parcel_diameter(trial(p),effective)/film%viscosity
            else
              call advance_spray_breakup(trial(p),effective,options,primitive(1,p),primitive(2:4,p), &
                film%density,film%viscosity,kick_dt,valid)
              if (.not.valid) exit
              next_mass(:,p)=next_mass(:,p)*(trial(p)%mass/old_mass)
            end if
          end if
          if(motion/='coupled') then
            trial(p)%position=drops(p)%position
            if(motion=='kick_drift'.and..not.present(eb_geometry)) &
              trial(p)%position=trial(p)%position+motion_dt*trial(p)%velocity
          end if
          if(present(eb_geometry)) then
            trial(p)%position=drops(p)%position
          else if(.not.present(boundary_faces)) then
            trial(p)%position=lower+modulo(trial(p)%position-lower,length)
          end if
        end do
        call agree()
        if (agreed.and.options%khrt_children) then
          if(allocated(eligible)) deallocate(eligible)
          allocate(eligible(np));eligible=.true.
          if(present(eb_geometry).and.motion_dt>0) then
            ! A parent crossing the EB takes the wall path, with no KH release
            ! on that substep. Unhit parents and their children then drift together.
            do p=1,np
              if(trial(p)%mass<=0) cycle
              if(present(boundary_faces)) then
                if(any(faces=='periodic')) then
                  ! Use the same segmented flight as the eventual drift, including
                  ! EB or reflecting-wall contacts after a periodic crossing.
                  preview(1)=trial(p);preview_contact=.false.
                  preview_change=0;preview_exterior=0;preview_hits=0
                  call drift_spray_level_set(nodal_phi,lower,lower+length,motion_dt,32,preview, &
                    preview_change,preview_hits,valid,faces,[0._dp],preview_exterior,preview_contact)
                  if(.not.valid) exit
                  eligible(p)=.not.preview_contact(1)
                  cycle
                end if
              end if
              crossing=.false.;crossing_fraction=0;crossing_position=0;crossing_normal=0;crossing_cell=0
              delta=motion_dt*trial(p)%velocity;finish=trial(p)%position+delta
              if(present(eb_outflow_change)) then
                flight_fraction=1
                do axis=1,3
                  if(delta(axis)>0.and.finish(axis)>=lower(axis)+length(axis)) then
                    flight_fraction=min(flight_fraction,(lower(axis)+length(axis)-trial(p)%position(axis))/delta(axis))
                  else if(delta(axis)<0.and.finish(axis)<=lower(axis)) then
                    flight_fraction=min(flight_fraction,(lower(axis)-trial(p)%position(axis))/delta(axis))
                  end if
                end do
                finish=max(lower,min(lower+length,trial(p)%position+flight_fraction*delta))
              end if
              call intersect_spray_level_set(nodal_phi,lower,lower+length,trial(p)%position, &
                finish,crossing,crossing_fraction, &
                crossing_position,crossing_normal,crossing_cell,valid)
              if(.not.valid) exit
              eligible(p)=.not.crossing
              if(present(boundary_faces)) then
                ! Domain-wall impact also takes precedence over KH release.
                do axis=1,3
                  if(delta(axis)<0.and.trial(p)%position(axis)+delta(axis)<=lower(axis)) then
                    if(faces(2*axis-1)=='slip_wall'.or.faces(2*axis-1)=='no_slip_wall') eligible(p)=.false.
                  else if(delta(axis)>0.and.trial(p)%position(axis)+delta(axis)>=lower(axis)+length(axis)) then
                    if(faces(2*axis)=='slip_wall'.or.faces(2*axis)=='no_slip_wall') eligible(p)=.false.
                  end if
                end do
              end if
            end do
            call agree()
          end if
          if(agreed) call advance_mpi_khrt_shedding(comm,components,options,trial,next_mass,primitive(1,:), &
            primitive(2:4,:),reynolds,kick_dt,trial_id,child_energy,valid,eligible)
          if (valid) valid=size(trial)<=125000
          call agree()
          if (agreed) change(5,:)=change(5,:)+child_energy
        end if
        if(allocated(impact_exterior)) impact_exterior=0
        if(agreed.and.present(eb_geometry).and..not.present(wall_films).and.motion_dt>0) then
          exchange=0;eb_hits=0
          if(allocated(impact_exterior)) then
            if(allocated(liquid_energies)) deallocate(liquid_energies)
            allocate(liquid_energies(size(trial)));liquid_energies=0;outflow_faces='outflow'
            if(present(boundary_faces)) outflow_faces=faces
            do p=1,size(trial)
              if(trial(p)%mass<=0) cycle
              call liquid_mixture_properties(components,next_mass(:,p)/trial(p)%mass,trial(p)%temperature,rho,cp,e,valid)
              if(.not.valid) exit
              liquid_energies(p)=e
            end do
            if(valid) call drift_spray_level_set(nodal_phi,lower,lower+length,motion_dt,32,trial,exchange,eb_hits,valid, &
              outflow_faces,liquid_energies,impact_exterior(:5))
            if(valid) then
              do p=1,size(trial)
                if(trial(p)%mass>0) cycle
                impact_exterior(6:)=impact_exterior(6:)-next_mass(:,p)*trial(p)%multiplicity
                next_mass(:,p)=0
              end do
            end if
            call agree()
            if(agreed) then
              call MPI_Allreduce(MPI_IN_PLACE,impact_exterior,size(impact_exterior),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
              valid=ierr==MPI_SUCCESS
            end if
          else
            call drift_spray_level_set(nodal_phi,lower,lower+length,motion_dt,32,trial,exchange,eb_hits,valid)
          end if
          if(valid) trial_ledger(:5)=exchange
          call agree()
        end if
        impact_wall=0
        if (agreed.and.present(wall_films)) then
          trial_patches=patches;trial_patch_mass=patch_mass
          call prepare_impacts()
          call agree()
          if (agreed.and.motion_dt>0) then
            if (present(eb_geometry)) then
              eb_hits=0
              call advance_mpi_level_set_impacts(comm,components,options,eb_walls,nodal_phi,mw,impact_pressure, &
                lower,lower+length,motion_dt,trial,next_mass,trial_patches,trial_patch_mass,trial_id,impact_wall, &
                eb_hits,valid,rejected,impact_exterior,boundary_faces)
              if (rejected) then
                if (present(model_rejected)) model_rejected=.true.
                return
              end if
            else
              call impact_mpi_planar_walls(comm,components,options,film_settings,faces,mw,impact_pressure,impact_previous, &
                lower,lower+length,trial,next_mass,trial_patches,trial_patch_mass,trial_id,hit,impact_wall,valid)
            end if
            if (valid) valid=size(trial)<=125000
            call agree()
          end if
        end if
        if (agreed) then
          if (present(boundary_faces).and.motion_dt>0) call apply_boundaries()
          call agree()
        end if
        if (agreed) then
          next_gas=gas;next_temp=temp
          call deposit_mpi_amr_tree_spray_gas(comm,species,tree,owners,next_gas,next_temp,positions,change, &
            lower,spacing,periodic,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi, &
            gas_change=trial_change,maximum_level=maximum_level)
          if (valid) then
            accepted=.true.;exit
          end if
        end if
        dt=.5_dp*dt
      end do
      if (.not.accepted) return
      if(present(gas_change)) then
        do node=1,size(owners)
          if(.not.allocated(integrated_change%fields(node)%values)) cycle
          associate(a=>integrated_change%fields(node)%values,b=>trial_change%fields(node)%values)
            if(any(a>huge(1._dp)-max(b,0._dp)).or.any(a< -huge(1._dp)-min(b,0._dp))) then
              valid=.false.;exit
            end if
            a=a+b
          end associate
        end do
        call agree()
        if(.not.agreed) return
      end if
      drops=trial;masses=next_mass;gas=next_gas;temp=next_temp;accepted_id=trial_id
      ledger=ledger+trial_ledger
      if(allocated(impact_exterior)) then
        exterior(:5)=exterior(:5)+impact_exterior(:5)
        do n=1,nc
          exterior(5+components(n)%vapor_index)=exterior(5+components(n)%vapor_index)+impact_exterior(5+n)
        end do
      end if
      if (present(wall_films)) then
        patches=trial_patches;patch_mass=trial_patch_mass;wall_ledger=wall_ledger+impact_wall
      end if
      if(present(wall_films).or.allocated(impact_exterior)) then
        kept=count(drops%mass>0);survivors=pack(drops,drops%mass>0)
        allocate(kept_mass(nc,kept));kept=0
        do p=1,size(drops)
          if (drops(p)%mass<=0) cycle
          kept=kept+1;kept_mass(:,kept)=masses(:,p)
        end do
        call move_alloc(survivors,drops);call move_alloc(kept_mass,masses)
      end if
      elapsed=elapsed+dt;count_steps=count_steps+1
    end do
    np=size(drops)
    kept=count(drops%mass>0);allocate(survivors(kept),kept_mass(nc,kept));kept=0
    do p=1,np
      if (drops(p)%mass<=0) cycle
      kept=kept+1;survivors(kept)=drops(p);kept_mass(:,kept)=masses(:,p)
    end do
    if(present(eb_geometry)) then
      call migrate_eb_tree_spray_parcels(comm,tree,owners,nodal_phi,survivors,lower,spacing, &
        options%maximum_parcels,valid,kept_mass,maximum_level=maximum_level)
    else
      call migrate_amr_tree_spray_parcels(comm,tree,owners,survivors,lower,spacing,periodic, &
        options%maximum_parcels,valid,kept_mass,maximum_level=maximum_level)
    end if
    if (.not.valid) return
    if (present(wall_films)) then
      call migrate_mpi_amr_wall_films(comm,tree,owners,patches,patch_mass,lower,spacing, &
        film_settings%maximum_films,valid,nodal_phi,maximum_level)
      if (.not.valid) return
    end if
    call MPI_Allreduce(ledger,global_ledger,size(ledger),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    global_ledger(:5)=global_ledger(:5)-wall_ledger
    if (.not.all(ieee_is_finite(global_ledger))) return
    if(present(gas_change)) call move_alloc(integrated_change%fields,gas_change%fields)
    state=gas;temperature=temp
    call move_alloc(survivors,parcels);call move_alloc(kept_mass,component_mass)
    if (present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    if (present(last_id)) last_id=accepted_id
    if (present(boundary_change)) boundary_change=global_ledger
    if(present(eb_outflow_change)) eb_outflow_change=exterior
    substeps=count_steps;ok=.true.
  contains
    subroutine prepare_impacts()
      integer :: q,parent,lo,hi,mid
      if (allocated(impact_previous)) deallocate(impact_previous,impact_pressure)
      allocate(impact_previous(3,size(trial)),impact_pressure(size(trial)))
      valid=.true.
      do q=1,size(trial)
        parent=q
        if (q>np) then
          lo=1;hi=np;parent=0
          do while(lo<=hi)
            mid=lo+(hi-lo)/2
            if (drops(mid)%id==trial(q)%parent_id) then
              parent=mid;exit
            else if (drops(mid)%id<trial(q)%parent_id) then
              lo=mid+1
            else
              hi=mid-1
            end if
          end do
          if (parent==0) then
            valid=.false.;return
          end if
        end if
        impact_previous(:,q)=positions(:,parent);impact_pressure(q)=1
        if (drops(parent)%mass>0) impact_pressure(q)=primitive(5,parent)
      end do
    end subroutine
    subroutine apply_boundaries()
      integer :: q,c
      valid=.true.
      do q=1,size(trial)
        if (trial(q)%mass<=0) cycle
        y=next_mass(:,q)/trial(q)%mass
        call liquid_mixture_properties(components,y,trial(q)%temperature,rho,cp,e,valid)
        if (.not.valid) return
        effective=components(1);effective%density=rho;effective%density_model='constant'
        effective%density_coefficients=0;effective%density_temperature_range=0
        effective%heat_capacity=cp;effective%heat_capacity_model='constant'
        effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
        effective%reference_temperature=trial(q)%temperature;effective%reference_energy=e
        effective%boiling_model='reject';effective%critical_temperature=0
        effective%normal_boiling_temperature=0
        call apply_parcel_boundaries(trial(q),effective,lower,lower+length,faces,exchange,valid)
        if (.not.valid) return
        trial_ledger(:5)=trial_ledger(:5)+exchange
        do c=1,nc
          trial_ledger(5+components(c)%vapor_index)=trial_ledger(5+components(c)%vapor_index)+exchange(1)*y(c)
        end do
        if (trial(q)%mass<=0) next_mass(:,q)=0
      end do
      valid=all(ieee_is_finite(trial_ledger))
    end subroutine
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
end module
