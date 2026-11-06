module spray_mpi_amr_step_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use mixture_thermo_mod, only: mixture_pressure
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory,liquid_mixture_properties, &
    advance_mixture_frozen_gas,mixture_transfer_rates,mixture_trace_fraction,build_mixture_film, &
    apply_mixture_breakup_viscosity
  use spray_breakup_mod, only: advance_spray_breakup,advance_tab_no_birth_khrt
  use spray_coupling_3d_mod, only: valid_spray_boundary_faces,apply_parcel_boundaries
  use spray_mpi_mod, only: advance_mpi_khrt_shedding,split_mpi_reference_etab_parcels,impact_mpi_planar_walls, &
    compute_mpi_spray_timestep
  use spray_amr_tree_partition_mod, only: locate_amr_tree_spray_cells
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,valid_wall_options
  use spray_mpi_amr_film_mod, only: migrate_mpi_amr_wall_films
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_remap_3d_mod, only: covered
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
  public :: advance_mpi_amr_population_sdc_stage,advance_mpi_amr_recursive_population_stage
  public :: compute_mpi_amr_spray_timestep
  public :: convert_mpi_amr_tree_rotating_energy,validate_mpi_spray_internal_energy_halo
  type, public :: spray_sdc_population_stage
    integer :: next_stage=0,iterations=0,depth=-1
    real(dp) :: native_flow_cfl=-1,time=0,interval=0
    type(spray_level_population) :: midpoint(3)
    type(amr_tree_state_3d) :: old_source,old_remainder,old_physical_source
    logical :: exports_physical=.false.,impacts=.false.
    type(spray_wall_film), allocatable :: midpoint_films(:)
    real(dp), allocatable :: midpoint_film_mass(:,:)
    real(dp) :: old_wall(5)=0,force_change(5)=0,acceleration(3)=0
    real(dp) :: frame_source_change=0,frame_boundary_change=0,frame_values(4)=0
    logical :: has_frame=.false.,has_internal_energy=.false.
    integer :: frame_axis=0
    real(dp), allocatable :: old_boundary(:)
  end type spray_sdc_population_stage
contains
  ! Validate the independent scalar on the same owner-local valid/halo
  ! layout before optional branches enter collective sampling. No projection
  ! of total energy is stored in the caller or used as a substitute for UEINT.
  subroutine validate_mpi_spray_internal_energy_halo(comm,owners,gas,held,ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(in) :: gas
    type(amr_tree_state_3d), optional, intent(in) :: held
    logical, intent(out) :: ok
    integer :: flag,lower,upper,status,rank,p,shape_gas(4),bounds_gas(4)
    logical :: valid
    ok=.false.;flag=merge(1,0,present(held))
    call MPI_Allreduce(flag,lower,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flag,upper,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.lower/=upper) return
    if(.not.present(held)) then
      ok=.true.;return
    endif
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    valid=allocated(gas%fields).and.allocated(held%fields)
    if(valid) valid=size(gas%fields)==size(owners).and.size(held%fields)==size(owners).and. &
      lbound(gas%fields,1)==1.and.lbound(held%fields,1)==1
    if(valid) then
      do p=1,size(owners)
        if(allocated(gas%fields(p)%values).neqv.allocated(held%fields(p)%values)) then
          valid=.false.;exit
        endif
        if(owners(p)/=rank) then
          if(allocated(held%fields(p)%values)) valid=.false.
          cycle
        endif
        if(.not.allocated(gas%fields(p)%values)) cycle
        shape_gas=shape(gas%fields(p)%values);bounds_gas=lbound(gas%fields(p)%values)
        valid=all(shape(held%fields(p)%values)==[1,shape_gas(2:4)]).and. &
          all(lbound(held%fields(p)%values)==[1,bounds_gas(2:4)])
        if(.not.valid) exit
        valid=all(ieee_is_finite(held%fields(p)%values))
        if(.not.valid) exit
      enddo
    endif
    call MPI_Allreduce(valid,ok,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS) ok=.false.
  end subroutine validate_mpi_spray_internal_energy_halo

  subroutine compute_mpi_amr_spray_timestep(comm,tree,parcels,lower,spacing,periodic,options, &
      maximum_injection_speed,subcycling,dt,ok)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    type(spray_parcel), intent(in) :: parcels(:)
    real(dp), intent(in) :: lower(3),spacing(3),maximum_injection_speed
    logical, intent(in) :: periodic(3),subcycling
    type(spray_options), intent(in) :: options
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),locations(:,:),levels(:)
    type(spray_parcel), allocatable :: active(:),selected(:)
    real(dp), allocatable :: positions(:,:)
    real(dp) :: h(3),bound,result,factor
    integer :: flags(5),minimum(5),maximum(5),level,p,status
    logical :: valid,agreed
    dt=0;ok=.false.
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    flags=[size(ratios),merge(1,0,subcycling),merge(1,0,periodic)]
    call MPI_Allreduce(flags,minimum,5,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,maximum,5,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(minimum/=maximum)) return
    allocate(selected(0))
    call compute_mpi_spray_timestep(comm,selected,spacing,options%particle_cfl,options%fixed_parts, &
      maximum_injection_speed,bound,valid)
    if(.not.valid) return
    if(options%fixed_parts) then
      dt=bound;ok=.true.;return
    endif
    active=pack(parcels,parcels%id>0.and.parcels%mass>0)
    allocate(positions(3,size(active)),levels(size(active)))
    do p=1,size(active)
      positions(:,p)=active(p)%position
    enddo
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,locations,valid)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    do p=1,size(active)
      levels(p)=tree%nodes(locations(1,p))%level+1
    enddo
    result=huge(1._dp);factor=1
    do level=1,size(ratios)
      h=spacing*real(domain(:,1),dp)/real(domain(:,level),dp)
      selected=pack(active,levels==level)
      call compute_mpi_spray_timestep(comm,selected,h,options%particle_cfl,.false., &
        maximum_injection_speed,bound,valid)
      if(.not.valid) return
      ! PeleC computeNewDt converts each level's bound to a root interval
      ! using the cumulative n_cycle product before taking the minimum.
      if(subcycling) factor=factor*real(ratios(level),dp)
      if(bound>huge(1._dp)/factor) then
        bound=huge(1._dp)
      else
        bound=bound*factor
      endif
      result=min(result,bound)
    enddo
    dt=result;ok=.true.
  end subroutine compute_mpi_amr_spray_timestep

  subroutine advance_mpi_amr_recursive_population_stage(comm,species,transport,components,options,tree,owners, &
      state,temperature,width,target_level,grow,maximum_parcels,physical,virtuals,ghosts,work, &
      time,interval,stage,iterations,faces, &
      source,remainder,wall_change,boundary_change,ok,lower,spacing,maximum_received,physical_source,acceleration,last_id, &
      wall_films,film_mass,film_settings,model_rejected,frame,common_step,flow_cfl, &
      internal_energy_halo)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_halo

    real(dp), optional, intent(in) :: flow_cfl
    logical, optional, intent(in) :: common_step
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    logical, optional, intent(out) :: model_rejected
    integer(int64), optional, intent(inout) :: last_id
    integer(int64), allocatable :: watermark
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),width,target_level,grow,maximum_parcels
    type(spray_level_population), intent(inout) :: physical(:),virtuals(:),ghosts(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: time,interval,lower(3),spacing(3)
    integer, intent(in) :: stage,iterations
    type(spray_sdc_population_stage), intent(inout) :: work(:)
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5),boundary_change(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_received
    real(dp), optional, intent(in) :: acceleration(3)
    type(amr_tree_state_3d), optional, intent(inout) :: physical_source
    type(amr_tree_state_3d), allocatable :: physical_delta
    type(spray_sdc_population_stage) :: candidate
    type(spray_level_population) :: group(3),new_virtual,next_ghost
    type(amr_tree_state_3d) :: change,leftover
    real(dp) :: wall(5),exterior(reactive_nvar(size(species))),tolerance,finish
    integer :: levels,level,p,status,flags(5),low(5),high(5)
    real(dp) :: times(2),low_time(2),high_time(2)
    logical :: valid,agreed,periodic(3),common
    if(present(last_id)) then
      allocate(watermark);watermark=last_id
    end if
    if(present(model_rejected)) model_rejected=.false.
    common=.false.
    if(present(common_step)) common=common_step
    ok=.false.;valid=tree%is_valid()
    call agree();if(.not.agreed) return
    call validate_mpi_spray_internal_energy_halo(comm,owners,state,internal_energy_halo,valid)
    if(.not.valid) return
    levels=tree%level_count();level=target_level+1
    valid=all([size(physical),size(virtuals),size(ghosts),size(work)]==levels)
    valid=valid.and.level>=1.and.level<=levels.and.all(ieee_is_finite([time,interval]))
    valid=valid.and.time>=0.and.interval>0.and.size(boundary_change)==size(exterior)
    call agree();if(.not.agreed) return
    flags=[levels,target_level,stage,iterations,merge(1,0,common)];times=[time,interval]
    call MPI_Allreduce(flags,low,5,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,5,MPI_INTEGER,MPI_MAX,comm,status)
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
        if(common) then
          ! Common intervals enter all old stages before any new stage.
          valid=valid.and.work(level-1)%next_stage==1.and.work(level-1)%iterations==iterations
          valid=valid.and.work(level-1)%time==time.and.work(level-1)%interval==interval
        else
          valid=valid.and.work(level-1)%iterations>0.and.work(level-1)%next_stage==0
          valid=valid.and.time>=work(level-1)%time-tolerance
          valid=valid.and.interval<=work(level-1)%time+work(level-1)%interval-time+tolerance
        end if
      end if
    end if
    call agree();if(.not.agreed) return
    periodic=[faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']
    if(present(physical_source)) allocate(physical_delta)
    candidate=work(level);group(1)=physical(level);group(2)=virtuals(level);group(3)=ghosts(level)
    if(stage==0) then
      call refresh_mpi_amr_level_copies(comm,tree,owners,physical,lower,spacing,periodic,target_level,grow, &
        maximum_parcels,new_virtual,next_ghost,valid)
      if(.not.valid) return
      group(2)=new_virtual
    end if
    call advance_mpi_amr_population_sdc_stage(comm,species,transport,components,options,tree,owners, &
      state,temperature,width,target_level,grow,maximum_parcels,group,candidate,time,interval,stage,iterations,faces, &
      change,leftover,wall,exterior,valid,lower,spacing,maximum_received,physical_delta,acceleration,watermark, &
      wall_films,film_mass,film_settings,model_rejected,frame,flow_cfl=flow_cfl, &
        internal_energy_halo=internal_energy_halo)
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
  end subroutine advance_mpi_amr_recursive_population_stage

  subroutine advance_mpi_amr_population_stage(comm,species,transport,components,options,tree,owners, &
      state,temperature,width,target_level,grow,maximum_parcels,populations,kick_dt,drift_dt,faces, &
      source,remainder,wall_change,boundary_change,ok,lower,spacing,maximum_received,physical_source,last_id, &
      wall_films,film_mass,film_settings,model_rejected,frame,acceleration,force_change,tab_interval,flow_cfl, &
      frame_source_change,frame_boundary_change, &
      internal_energy_halo)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_halo

    real(dp), optional, intent(inout) :: frame_source_change,frame_boundary_change
    real(dp) :: physical_conversion,physical_boundary_conversion,part_boundary_conversion
    integer :: ledger_flag,ledger_low,ledger_high
    real(dp), optional, intent(in) :: flow_cfl
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(in) :: tab_interval
    real(dp), optional, intent(out) :: force_change(5)
    integer(int64), optional, intent(inout) :: last_id
    integer(int64), allocatable :: watermark
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    logical, optional, intent(out) :: model_rejected
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    integer :: film_flag,film_low,film_high
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),width,target_level,grow,maximum_parcels
    type(spray_level_population), intent(inout) :: populations(3)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: kick_dt,drift_dt,lower(3),spacing(3)
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5),boundary_change(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_received
    type(amr_tree_state_3d), optional, intent(inout) :: physical_source
    type(amr_tree_state_3d), allocatable :: physical_delta
    type(amr_tree_state_3d) :: exported
    type(spray_options) :: population_options
    type(spray_level_population) :: candidate(3)
    type(amr_tree_state_3d) :: halo,halo_t,total,leftover,part,part_leftover
    real(dp) :: wall(5),physical_wall(5)
    real(dp) :: exterior(reactive_nvar(size(species))),physical_exterior(reactive_nvar(size(species)))
    character(len=32) :: particle_faces(6)
    integer :: population,status
    logical :: valid,agreed,periodic(3)
    real(dp) :: part_force(5),physical_force(5)
    if(present(force_change)) force_change=0
    ok=.false.
    physical_force=0;physical_conversion=0;physical_boundary_conversion=0
    ledger_flag=merge(1,0,present(frame_source_change))+2*merge(1,0,present(frame_boundary_change))
    call MPI_Allreduce(ledger_flag,ledger_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(ledger_flag,ledger_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.ledger_low/=ledger_high) return
    if(present(last_id)) then
      allocate(watermark);watermark=last_id
    end if
    ok=.false.
    if(present(model_rejected)) model_rejected=.false.
    film_flag=merge(1,0,present(wall_films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))
    call MPI_Allreduce(film_flag,film_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,film_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.film_low/=film_high) return
    if(film_flag/=0.and.film_flag/=7) return
    if(present(wall_films)) then
      valid=allocated(wall_films).and.allocated(film_mass).and.present(last_id)
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS.or..not.agreed) return
      patches=wall_films;patch_mass=film_mass
    end if
    valid=size(boundary_change)==size(exterior)
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    do population=1,3
      valid=valid.and.allocated(populations(population)%parcels).and. &
        allocated(populations(population)%component_mass).and.allocated(populations(population)%nodes)
    end do
    if(valid.and.target_level==0) valid=size(populations(3)%parcels)==0
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call validate_mpi_spray_internal_energy_halo(comm,owners,state,internal_energy_halo,valid)
    if(.not.valid) return
    halo=state;halo_t=temperature
    if(present(frame)) then
      call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,halo,lower,spacing,frame,.false.,valid,halo_width=width)
      if(.not.valid) return
    end if
    periodic=[faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']
    particle_faces=faces
    where(particle_faces=='inflow') particle_faces='outflow'
    where(particle_faces=='isothermal_wall') particle_faces='no_slip_wall'
    if(present(physical_source).or.present(frame)) allocate(physical_delta)
    candidate=populations
    do population=1,3
      population_options=options
      ! Root has no ghosts, but use retained routing for its allocated empty set.
      ! Pinned SprayParticles.cpp gates splash by isActive and do_move.
      if(population==1.and.drift_dt>0.and.present(wall_films)) then
      call advance_mpi_amr_population_halo(comm,species,transport,components,population_options,tree,owners,halo,halo_t, &
        width,target_level,grow,maximum_parcels,candidate(population)%nodes,candidate(population)%parcels, &
        candidate(population)%component_mass,lower,spacing,periodic,kick_dt,drift_dt,particle_faces,part,part_leftover, &
        wall,exterior,valid,maximum_received=maximum_received,keep_physical=population<3.or.target_level==0, &
        containing_cell=.true.,physical_source=physical_delta,last_id=watermark, &
        wall_films=patches,film_mass=patch_mass,film_settings=film_settings,model_rejected=model_rejected,frame=frame, &
        acceleration=acceleration,force_change=part_force,tab_interval=tab_interval, &
        etab_births=population==1.and.drift_dt>0,khrt_births=population==1.and.drift_dt>0,flow_cfl=flow_cfl, &
        frame_boundary_change=part_boundary_conversion,internal_energy_halo=internal_energy_halo)
      else
      call advance_mpi_amr_population_halo(comm,species,transport,components,population_options,tree,owners,halo,halo_t, &
        width,target_level,grow,maximum_parcels,candidate(population)%nodes,candidate(population)%parcels, &
        candidate(population)%component_mass,lower,spacing,periodic,kick_dt,drift_dt,particle_faces,part,part_leftover, &
        wall,exterior,valid,maximum_received=maximum_received,keep_physical=population<3.or.target_level==0, &
        containing_cell=.true.,physical_source=physical_delta,last_id=watermark,frame=frame, &
        acceleration=acceleration,force_change=part_force,tab_interval=tab_interval, &
        etab_births=population==1.and.drift_dt>0,khrt_births=population==1.and.drift_dt>0,flow_cfl=flow_cfl, &
        frame_boundary_change=part_boundary_conversion,internal_energy_halo=internal_energy_halo)
      end if
      if(.not.valid) return
      if(present(frame)) then
        call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,part,lower,spacing,frame,.true.,valid)
        if(.not.valid) return
        call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,part_leftover,lower,spacing, &
          frame,.true.,valid,halo_width=width)
        if(.not.valid) return
        if(allocated(physical_delta).and.population==1) then
          call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,physical_delta,lower,spacing, &
            frame,.true.,valid,energy_change=physical_conversion)
          if(.not.valid) return
        end if
      end if
      if(population==1) then
        physical_force=part_force;physical_boundary_conversion=part_boundary_conversion
        if(present(physical_source)) exported=physical_delta
        if(allocated(physical_delta)) deallocate(physical_delta)
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
    if(present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    if(present(last_id)) last_id=watermark
    if(present(force_change)) force_change=physical_force
    if(present(frame_source_change)) frame_source_change=physical_conversion
    if(present(frame_boundary_change)) frame_boundary_change=physical_boundary_conversion
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
  end subroutine advance_mpi_amr_population_stage

  subroutine advance_mpi_amr_population_sdc_stage(comm,species,transport,components,options,tree,owners, &
      state,temperature,width,target_level,grow,maximum_parcels,populations,work,time,interval,stage,iterations,faces, &
      source,remainder,wall_change,boundary_change,ok,lower,spacing,maximum_received,physical_source,acceleration,last_id, &
      wall_films,film_mass,film_settings,model_rejected,frame,flow_cfl, &
      internal_energy_halo)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_halo

    real(dp), optional, intent(in) :: flow_cfl
    type(rotating_frame_3d), optional, intent(in) :: frame
    integer(int64), optional, intent(inout) :: last_id
    integer(int64), allocatable :: watermark
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    logical, optional, intent(out) :: model_rejected
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    integer :: film_flag,film_low,film_high
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),width,target_level,grow,maximum_parcels
    type(spray_level_population), intent(inout) :: populations(3)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: time,interval,lower(3),spacing(3)
    integer, intent(in) :: stage,iterations
    type(spray_sdc_population_stage), intent(inout) :: work
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5),boundary_change(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_received
    real(dp), optional, intent(in) :: acceleration(3)
    type(amr_tree_state_3d), optional, intent(inout) :: physical_source
    type(amr_tree_state_3d), allocatable :: physical_delta
    type(spray_sdc_population_stage) :: candidate
    type(spray_level_population) :: trial(3)
    type(amr_tree_state_3d) :: change,leftover
    real(dp) :: half_interval,drift,wall(5),exterior(reactive_nvar(size(species)))
    real(dp) :: times(6),low_time(6),high_time(6),native_cfl,accel(3),forced(5),coupled_accel(3),stage_force(5)
    integer :: flags(5),low(5),high(5),status
    logical :: valid,agreed
    real(dp) :: frame_source,frame_boundary
    if(present(last_id)) then
      allocate(watermark);watermark=last_id
    end if
    ok=.false.;accel=0;forced=0;coupled_accel=0
    call validate_mpi_spray_internal_energy_halo(comm,owners,state,internal_energy_halo,valid)
    if(.not.valid) return
    if(present(model_rejected)) model_rejected=.false.
    film_flag=merge(1,0,present(wall_films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))
    call MPI_Allreduce(film_flag,film_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,film_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.film_low/=film_high) return
    if(film_flag/=0.and.film_flag/=7) return
    valid=.true.
    if(present(wall_films)) valid=allocated(wall_films).and.allocated(film_mass).and.present(last_id)
    if(stage>0) then
      valid=valid.and.(present(wall_films).eqv.work%impacts)
      if(present(wall_films)) valid=valid.and.allocated(work%midpoint_films).and.allocated(work%midpoint_film_mass)
    end if
    call agree();if(.not.agreed) return
    if(present(acceleration)) accel=acceleration
    if(options%mixture_time_scheme=='euler') coupled_accel=accel
    native_cfl=-1
    if(present(flow_cfl)) native_cfl=flow_cfl
    times=[time,interval,accel,native_cfl]
    valid=all(ieee_is_finite(times)).and.time>=0.and.interval>0.and.iterations>=1
    valid=valid.and.all(abs(accel)<=1e6_dp).and.native_cfl>=-1
    valid=valid.and.stage>=0.and.stage<=iterations.and.size(boundary_change)==size(exterior)
    if(stage==0) then
      valid=valid.and.work%next_stage==0
    else
      valid=valid.and.stage==work%next_stage.and.iterations==work%iterations.and.target_level==work%depth
      valid=valid.and.all(accel==work%acceleration).and.native_cfl==work%native_flow_cfl
      valid=valid.and.(present(physical_source).eqv.work%exports_physical)
      valid=valid.and.(present(frame).eqv.work%has_frame)
      valid=valid.and.(present(internal_energy_halo).eqv.work%has_internal_energy)
      if(present(frame)) then
        valid=valid.and.frame%axis==work%frame_axis
        valid=valid.and.all([frame%angular_speed,frame%axis_point]==work%frame_values)
      endif
      valid=valid.and.time==work%time.and.interval==work%interval.and.allocated(work%old_boundary)
      if(allocated(work%old_boundary)) valid=valid.and.size(work%old_boundary)==size(exterior)
    end if
    call agree();if(.not.agreed) return
    call MPI_Allreduce(times,low_time,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(times,high_time,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low_time/=high_time)) return
    flags=[stage,iterations,target_level,merge(1,0,present(physical_source)),merge(1,0,present(frame))]
    call MPI_Allreduce(flags,low,5,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,5,MPI_INTEGER,MPI_MAX,comm,status)
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
      if(options%mixture_time_scheme/='euler'.and..not.options%fixed_parts) then
        call apply_force()
        if(.not.agreed) return
      endif
    else
      trial=work%midpoint
    end if
    if(present(wall_films)) then
      if(stage==0) then
        patches=wall_films;patch_mass=film_mass
      else
        patches=work%midpoint_films;patch_mass=work%midpoint_film_mass
      end if
    call advance_mpi_amr_population_stage(comm,species,transport,components,options,tree,owners, &
      state,temperature,width,target_level,grow,maximum_parcels,trial,half_interval,drift,faces, &
      change,leftover,wall,exterior,valid,lower,spacing,maximum_received,physical_delta,watermark, &
      patches,patch_mass,film_settings,model_rejected,frame, &
        acceleration=coupled_accel,force_change=stage_force,tab_interval=interval,flow_cfl=flow_cfl, &
        frame_source_change=frame_source,frame_boundary_change=frame_boundary, &
        internal_energy_halo=internal_energy_halo)
    else
    call advance_mpi_amr_population_stage(comm,species,transport,components,options,tree,owners, &
      state,temperature,width,target_level,grow,maximum_parcels,trial,half_interval,drift,faces, &
      change,leftover,wall,exterior,valid,lower,spacing,maximum_received,physical_delta,watermark,frame=frame, &
        acceleration=coupled_accel,force_change=stage_force,tab_interval=interval,flow_cfl=flow_cfl, &
        frame_source_change=frame_source,frame_boundary_change=frame_boundary, &
        internal_energy_halo=internal_energy_halo)
    end if
    if(.not.valid) return
    if(options%mixture_time_scheme=='euler'.or.options%fixed_parts) forced=stage_force
    if(stage==iterations.and.options%mixture_time_scheme/='euler'.and..not.options%fixed_parts) then
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
      candidate%acceleration=accel;candidate%native_flow_cfl=native_cfl;candidate%force_change=forced
      candidate%frame_source_change=frame_source;candidate%frame_boundary_change=frame_boundary
      candidate%has_internal_energy=present(internal_energy_halo)
      candidate%has_frame=present(frame);candidate%frame_axis=0;candidate%frame_values=0
      if(present(frame)) then
        candidate%frame_axis=frame%axis;candidate%frame_values=[frame%angular_speed,frame%axis_point]
      endif
      candidate%impacts=present(wall_films)
      if(present(wall_films)) then
        candidate%midpoint_films=patches;candidate%midpoint_film_mass=patch_mass
      end if
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
      candidate%frame_source_change=candidate%frame_source_change+frame_source
      candidate%frame_boundary_change=candidate%frame_boundary_change+frame_boundary
      valid=ieee_is_finite(candidate%frame_source_change).and.ieee_is_finite(candidate%frame_boundary_change)
      call agree();if(.not.agreed) return
      candidate%force_change=candidate%force_change+forced
      wall=wall+candidate%old_wall;exterior=exterior+candidate%old_boundary
    end if
    candidate%next_stage=stage+1
    if(stage==iterations) then
      candidate%next_stage=0
      populations=trial
      if(present(wall_films)) then
        call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
      end if
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
  end subroutine advance_mpi_amr_population_sdc_stage

  ! Convert an owner-local gas field at its cell centers. The same operation
  ! applies to a deposited increment, whose density component may be negative.
  ! Publish the candidate only after every rank has accepted every cell.
  subroutine convert_mpi_amr_tree_rotating_energy(comm,tree,owners,state,lower,spacing,frame,to_rotating,ok,geometry, &
      halo_width,maximum_depth,energy_change)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(inout) :: state
    real(dp), intent(in) :: lower(3),spacing(3)
    type(rotating_frame_3d), intent(in) :: frame
    logical, intent(in) :: to_rotating
    logical, intent(out) :: ok
    type(eb_geometry_3d), optional, intent(in) :: geometry(:)
    integer, optional, intent(in) :: halo_width
    integer, optional, intent(in) :: maximum_depth
    ! Integrated change of the storage convention on exposed valid cells.
    ! This is coordinate bookkeeping, not work by a particle body force.
    real(dp), optional, intent(inout) :: energy_change
    type(amr_tree_state_3d) :: candidate
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: frame_values(4),low_values(4),high_values(4),length(3),position(3),potential,term
    integer :: rank,status,node,level,i,j,k,low_axis,high_axis,low_flags(5),high_flags(5),flags(5),width,n(3),max_depth
    real(dp) :: bounds(6),low_bounds(6),high_bounds(6)
    real(dp) :: local_change,global_change,volume,transfer,amount
    logical :: valid,agreed
    ok=.false.
    width=0
    if(present(halo_width)) width=halo_width
    max_depth=tree%level_count()-1
    if(present(maximum_depth)) max_depth=maximum_depth
    flags=[merge(1,0,present(geometry)),merge(1,0,to_rotating),width,max_depth,merge(1,0,present(energy_change))]
    call MPI_Allreduce(flags,low_flags,5,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high_flags,5,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low_flags/=high_flags)) return
    valid=frame%is_valid().and.tree%is_valid().and.allocated(state%fields)
    valid=valid.and.all(ieee_is_finite([lower,spacing])).and.all(spacing>0).and.width>=0.and.width<=16
    valid=valid.and.max_depth>=0.and.max_depth<tree%level_count()
    if(valid) valid=size(owners)==size(tree%nodes).and.size(state%fields)==size(owners)
    if(valid.and.present(geometry)) valid=size(geometry)==size(owners)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(frame%axis,low_axis,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(frame%axis,high_axis,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low_axis/=high_axis) return
    frame_values=[frame%angular_speed,frame%axis_point]
    call MPI_Allreduce(frame_values,low_values,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(frame_values,high_values,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low_values/=high_values)) return
    bounds=[lower,spacing]
    call MPI_Allreduce(bounds,low_bounds,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,high_bounds,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low_bounds/=high_bounds)) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    candidate=state
    length=spacing*real(tree%nodes(1)%shape,dp)
    valid=.true.;local_change=0
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      if(tree%nodes(node)%level>max_depth) cycle
      if(.not.allocated(candidate%fields(node)%values)) then
        valid=.false.;exit
      end if
      n=tree%nodes(node)%shape
      if(any([size(candidate%fields(node)%values,2),size(candidate%fields(node)%values,3), &
          size(candidate%fields(node)%values,4)]/=n+2*width)) then
        valid=.false.;exit
      end if
      if(any([lbound(candidate%fields(node)%values,2),lbound(candidate%fields(node)%values,3), &
          lbound(candidate%fields(node)%values,4)]/=1-width)) then
        valid=.false.;exit
      end if
      if(size(candidate%fields(node)%values,1)<5) then
        valid=.false.;exit
      end if
      level=tree%nodes(node)%level+1
      if(present(geometry).and.present(energy_change)) then
        if(.not.geometry(node)%is_valid()) then
          valid=.false.;exit
        endif
        if(any([geometry(node)%nx,geometry(node)%ny,geometry(node)%nz]/=n)) then
          valid=.false.;exit
        endif
      endif
      do k=1-width,n(3)+width
      do j=1-width,n(2)+width
      do i=1-width,n(1)+width
        position=lower+length*(real(origin(:,node)+int([i,j,k],int64),dp)-.5_dp)/real(domain(:,level),dp)
        if(present(geometry)) position=[geometry(node)%x_lower+(real(i,dp)-.5_dp)*geometry(node)%dx, &
          geometry(node)%y_lower+(real(j,dp)-.5_dp)*geometry(node)%dy, &
          geometry(node)%z_lower+(real(k,dp)-.5_dp)*geometry(node)%dz]
        call rotating_potential_3d(frame,position,potential,valid)
        if(.not.valid) exit
        term=candidate%fields(node)%values(1,i,j,k)
        if(.not.ieee_is_finite(term)) then
          valid=.false.;exit
        endif
        if(potential>1) then
          if(abs(term)>huge(1._dp)/potential) then
            valid=.false.;exit
          endif
        endif
        term=term*potential
        if(to_rotating) term=-term
        if(.not.ieee_is_finite(term)) then
          valid=.false.;exit
        end if
        transfer=term
        amount=candidate%fields(node)%values(5,i,j,k)
        if(.not.ieee_is_finite(amount)) then
          valid=.false.;exit
        endif
        if(term>0) then
          if(amount>huge(1._dp)-term) then
            valid=.false.;exit
          endif
        else
          if(amount< -huge(1._dp)-term) then
            valid=.false.;exit
          endif
        endif
        term=amount+term
        if(.not.ieee_is_finite(term)) then
          valid=.false.;exit
        end if
        candidate%fields(node)%values(5,i,j,k)=term
        if(present(energy_change)) then
          if(all([i,j,k]>=1).and.all([i,j,k]<=n)) then
            if(tree%nodes(node)%level<max_depth) then
              if(covered(tree,node,[i,j,k])) cycle
            endif
            volume=product(length/real(domain(:,level),dp))
            if(present(geometry)) volume=geometry(node)%volume_fraction(i,j,k)* &
              geometry(node)%dx*geometry(node)%dy*geometry(node)%dz
            if(.not.ieee_is_finite(volume).or.volume<0) then
              valid=.false.;exit
            endif
            if(volume>1) then
              if(abs(transfer)>huge(1._dp)/volume) then
                valid=.false.;exit
              endif
            endif
            amount=volume*transfer
            if(amount>0) then
              if(local_change>huge(1._dp)-amount) then
                valid=.false.;exit
              endif
            else
              if(local_change< -huge(1._dp)-amount) then
                valid=.false.;exit
              endif
            endif
            local_change=local_change+amount
          endif
        endif
      end do
      if(.not.valid) exit
      end do
      if(.not.valid) exit
      end do
      if(.not.valid) exit
    end do
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    if(present(energy_change)) then
      call MPI_Allreduce(local_change,global_change,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
      if(status/=MPI_SUCCESS.or..not.ieee_is_finite(global_change)) return
    endif
    call move_alloc(candidate%fields,state%fields)
    if(present(energy_change)) energy_change=global_change
    ok=.true.
  end subroutine convert_mpi_amr_tree_rotating_energy

  ! Advance a caller-owned hierarchy one recursive SDC callback at a time.
  ! The complete hierarchy must stay private until its outer gas step succeeds.
  subroutine advance_mpi_eb_recursive_population_stage(comm,species,transport,components,options,tree,owners,geometry, &
      state,temperature,width,target_level,grow,maximum_parcels,physical,virtuals,ghosts,work, &
      time,interval,stage,iterations,faces, &
      source,remainder,wall_change,boundary_change,ok,nodal_phi,centroids,target_volume_fraction, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,maximum_received,physical_source,acceleration,last_id, &
      wall_films,film_mass,film_settings,model_rejected,frame,common_step,flow_cfl, &
      internal_energy,ancestor_internal_energy_start,ancestor_internal_energy_end)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_internal_energy_start,ancestor_internal_energy_end

    real(dp), optional, intent(in) :: flow_cfl
    logical, optional, intent(in) :: common_step
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    logical, optional, intent(out) :: model_rejected
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
    integer :: levels,level,p,status,flags(5),low(5),high(5)
    real(dp) :: times(2),low_time(2),high_time(2)
    logical :: valid,agreed,periodic(3),common
    if(present(last_id)) then
      allocate(watermark);watermark=last_id
    end if
    if(present(model_rejected)) model_rejected=.false.
    common=.false.
    if(present(common_step)) common=common_step
    call validate_mpi_spray_internal_energy_halo(comm,owners,state,internal_energy,valid)
    if(.not.valid) then
      ok=.false.;return
    endif
    ok=.false.;valid=tree%is_valid()
    call agree();if(.not.agreed) return
    levels=tree%level_count();level=target_level+1
    valid=all([size(physical),size(virtuals),size(ghosts),size(work)]==levels)
    valid=valid.and.level>=1.and.level<=levels.and.all(ieee_is_finite([time,interval]))
    valid=valid.and.time>=0.and.interval>0.and.size(boundary_change)==size(exterior)
    call agree();if(.not.agreed) return
    flags=[levels,target_level,stage,iterations,merge(1,0,common)];times=[time,interval]
    call MPI_Allreduce(flags,low,5,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,5,MPI_INTEGER,MPI_MAX,comm,status)
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
        if(common) then
          ! Common intervals enter all old stages before any new stage.
          valid=valid.and.work(level-1)%next_stage==1.and.work(level-1)%iterations==iterations
          valid=valid.and.work(level-1)%time==time.and.work(level-1)%interval==interval
        else
          valid=valid.and.work(level-1)%iterations>0.and.work(level-1)%next_stage==0
          valid=valid.and.time>=work(level-1)%time-tolerance
          valid=valid.and.interval<=work(level-1)%time+work(level-1)%interval-time+tolerance
        end if
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
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,maximum_received,physical_delta,acceleration,watermark, &
      wall_films,film_mass,film_settings,model_rejected,frame,flow_cfl=flow_cfl, &
        internal_energy=internal_energy,ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end)
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
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,maximum_received,physical_source,acceleration,last_id, &
      wall_films,film_mass,film_settings,model_rejected,frame,flow_cfl, &
      internal_energy,ancestor_internal_energy_start,ancestor_internal_energy_end)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_internal_energy_start,ancestor_internal_energy_end

    real(dp), optional, intent(in) :: flow_cfl
    type(rotating_frame_3d), optional, intent(in) :: frame
    integer(int64), optional, intent(inout) :: last_id
    integer(int64), allocatable :: watermark
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    logical, optional, intent(out) :: model_rejected
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    integer :: film_flag,film_low,film_high
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
    real(dp) :: times(6),low_time(6),high_time(6),native_cfl,accel(3),forced(5),coupled_accel(3),stage_force(5)
    integer :: flags(5),low(5),high(5),status
    logical :: valid,agreed
    real(dp) :: frame_source,frame_boundary
    if(present(last_id)) then
      allocate(watermark);watermark=last_id
    end if
    call validate_mpi_spray_internal_energy_halo(comm,owners,state,internal_energy,valid)
    if(.not.valid) then
      ok=.false.;return
    endif
    ok=.false.;accel=0;forced=0;coupled_accel=0
    if(present(model_rejected)) model_rejected=.false.
    film_flag=merge(1,0,present(wall_films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))
    call MPI_Allreduce(film_flag,film_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,film_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.film_low/=film_high) return
    if(film_flag/=0.and.film_flag/=7) return
    valid=.true.
    if(present(wall_films)) valid=allocated(wall_films).and.allocated(film_mass).and.present(last_id)
    if(stage>0) then
      valid=valid.and.(present(wall_films).eqv.work%impacts)
      if(present(wall_films)) valid=valid.and.allocated(work%midpoint_films).and.allocated(work%midpoint_film_mass)
    end if
    call agree();if(.not.agreed) return
    if(present(acceleration)) accel=acceleration
    if(options%mixture_time_scheme=='euler') coupled_accel=accel
    native_cfl=-1
    if(present(flow_cfl)) native_cfl=flow_cfl
    times=[time,interval,accel,native_cfl]
    valid=all(ieee_is_finite(times)).and.time>=0.and.interval>0.and.iterations>=1
    valid=valid.and.all(abs(accel)<=1e6_dp).and.native_cfl>=-1
    valid=valid.and.stage>=0.and.stage<=iterations.and.size(boundary_change)==size(exterior)
    if(stage==0) then
      valid=valid.and.work%next_stage==0
    else
      valid=valid.and.stage==work%next_stage.and.iterations==work%iterations.and.target_level==work%depth
      valid=valid.and.all(accel==work%acceleration).and.native_cfl==work%native_flow_cfl
      valid=valid.and.(present(physical_source).eqv.work%exports_physical)
      valid=valid.and.(present(frame).eqv.work%has_frame)
      valid=valid.and.(present(internal_energy).eqv.work%has_internal_energy)
      if(present(frame)) then
        valid=valid.and.frame%axis==work%frame_axis
        valid=valid.and.all([frame%angular_speed,frame%axis_point]==work%frame_values)
      endif
      valid=valid.and.time==work%time.and.interval==work%interval.and.allocated(work%old_boundary)
      if(allocated(work%old_boundary)) valid=valid.and.size(work%old_boundary)==size(exterior)
    end if
    call agree();if(.not.agreed) return
    call MPI_Allreduce(times,low_time,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(times,high_time,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low_time/=high_time)) return
    flags=[stage,iterations,target_level,merge(1,0,present(physical_source)),merge(1,0,present(frame))]
    call MPI_Allreduce(flags,low,5,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,5,MPI_INTEGER,MPI_MAX,comm,status)
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
      if(options%mixture_time_scheme/='euler'.and..not.options%fixed_parts) then
        call apply_force()
        if(.not.agreed) return
      endif
    else
      trial=work%midpoint
    end if
    if(present(wall_films)) then
      if(stage==0) then
        patches=wall_films;patch_mass=film_mass
      else
        patches=work%midpoint_films;patch_mass=work%midpoint_film_mass
      end if
    call advance_mpi_eb_population_stage(comm,species,transport,components,options,tree,owners,geometry, &
      state,temperature,width,target_level,grow,maximum_parcels,trial,half_interval,drift,faces, &
      change,leftover,wall,exterior,valid,nodal_phi,centroids,target_volume_fraction, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,maximum_received,physical_delta,watermark, &
      patches,patch_mass,film_settings,model_rejected,frame, &
        acceleration=coupled_accel,force_change=stage_force,tab_interval=interval,flow_cfl=flow_cfl, &
        frame_source_change=frame_source,frame_boundary_change=frame_boundary, &
        internal_energy=internal_energy,ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end)
    else
    call advance_mpi_eb_population_stage(comm,species,transport,components,options,tree,owners,geometry, &
      state,temperature,width,target_level,grow,maximum_parcels,trial,half_interval,drift,faces, &
      change,leftover,wall,exterior,valid,nodal_phi,centroids,target_volume_fraction, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,maximum_received,physical_delta,watermark,frame=frame, &
        acceleration=coupled_accel,force_change=stage_force,tab_interval=interval,flow_cfl=flow_cfl, &
        frame_source_change=frame_source,frame_boundary_change=frame_boundary, &
        internal_energy=internal_energy,ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end)
    end if
    if(.not.valid) return
    if(options%mixture_time_scheme=='euler'.or.options%fixed_parts) forced=stage_force
    if(stage==iterations.and.options%mixture_time_scheme/='euler'.and..not.options%fixed_parts) then
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
      candidate%acceleration=accel;candidate%native_flow_cfl=native_cfl;candidate%force_change=forced
      candidate%frame_source_change=frame_source;candidate%frame_boundary_change=frame_boundary
      candidate%has_internal_energy=present(internal_energy)
      candidate%has_frame=present(frame);candidate%frame_axis=0;candidate%frame_values=0
      if(present(frame)) then
        candidate%frame_axis=frame%axis;candidate%frame_values=[frame%angular_speed,frame%axis_point]
      endif
      candidate%impacts=present(wall_films)
      if(present(wall_films)) then
        candidate%midpoint_films=patches;candidate%midpoint_film_mass=patch_mass
      end if
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
      candidate%frame_source_change=candidate%frame_source_change+frame_source
      candidate%frame_boundary_change=candidate%frame_boundary_change+frame_boundary
      valid=ieee_is_finite(candidate%frame_source_change).and.ieee_is_finite(candidate%frame_boundary_change)
      call agree();if(.not.agreed) return
      candidate%force_change=candidate%force_change+forced
      wall=wall+candidate%old_wall;exterior=exterior+candidate%old_boundary
    end if
    candidate%next_stage=stage+1
    if(stage==iterations) then
      candidate%next_stage=0
      populations=trial
      if(present(wall_films)) then
        call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
      end if
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
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,maximum_received,physical_source,last_id, &
      wall_films,film_mass,film_settings,model_rejected,frame,acceleration,force_change,tab_interval,flow_cfl, &
      frame_source_change,frame_boundary_change, &
      internal_energy,ancestor_internal_energy_start,ancestor_internal_energy_end)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_internal_energy_start,ancestor_internal_energy_end

    real(dp), optional, intent(inout) :: frame_source_change,frame_boundary_change
    real(dp) :: physical_conversion,physical_boundary_conversion,part_boundary_conversion
    integer :: ledger_flag,ledger_low,ledger_high
    real(dp), optional, intent(in) :: flow_cfl
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(in) :: tab_interval
    real(dp), optional, intent(out) :: force_change(5)
    integer(int64), optional, intent(inout) :: last_id
    integer(int64), allocatable :: watermark
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    logical, optional, intent(out) :: model_rejected
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    integer :: film_flag,film_low,film_high
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
    type(amr_tree_state_3d), allocatable :: held_halo
    integer :: held_allocation_status
    integer :: population,status
    logical :: valid,agreed,periodic(3)
    real(dp) :: part_force(5),physical_force(5)
    if(present(force_change)) force_change=0
    call validate_mpi_spray_internal_energy_halo(comm,owners,state,internal_energy,valid)
    if(.not.valid) then
      ok=.false.;return
    endif
    ok=.false.
    physical_force=0;physical_conversion=0;physical_boundary_conversion=0
    ledger_flag=merge(1,0,present(frame_source_change))+2*merge(1,0,present(frame_boundary_change))
    call MPI_Allreduce(ledger_flag,ledger_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(ledger_flag,ledger_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.ledger_low/=ledger_high) return
    if(present(last_id)) then
      allocate(watermark);watermark=last_id
    end if
    ok=.false.
    if(present(model_rejected)) model_rejected=.false.
    film_flag=merge(1,0,present(wall_films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))
    call MPI_Allreduce(film_flag,film_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,film_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.film_low/=film_high) return
    if(film_flag/=0.and.film_flag/=7) return
    if(present(wall_films)) then
      valid=allocated(wall_films).and.allocated(film_mass).and.present(last_id)
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS.or..not.agreed) return
      patches=wall_films;patch_mass=film_mass
    end if
    valid=size(boundary_change)==size(exterior)
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    do population=1,3
      valid=valid.and.allocated(populations(population)%parcels).and. &
        allocated(populations(population)%component_mass).and.allocated(populations(population)%nodes)
    end do
    if(valid.and.target_level==0) valid=size(populations(3)%parcels)==0
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    if(present(internal_energy)) then
      allocate(held_halo,stat=held_allocation_status)
      call MPI_Allreduce(held_allocation_status==0,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS.or..not.agreed) return
    endif
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,halo,halo_t,comm,valid, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_mach,inflow_temperature,inflow_pressure,volume,width,target_volume_fraction, &
      faces,inflow_state,restrict_snapshot=.false.,frame=frame, &
        internal_energy=internal_energy,ghost_internal_energy=held_halo, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end)
    if(.not.valid) return
    lower=[geometry(1)%x_lower,geometry(1)%y_lower,geometry(1)%z_lower]
    spacing=[geometry(1)%dx,geometry(1)%dy,geometry(1)%dz]
    if(present(frame)) then
      call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,halo,lower,spacing,frame,.false.,valid,geometry,width)
      if(.not.valid) return
    end if
    periodic=[faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']
    particle_faces=faces
    where(particle_faces=='inflow') particle_faces='outflow'
    where(particle_faces=='isothermal_wall') particle_faces='no_slip_wall'
    if(present(physical_source).or.present(frame)) allocate(physical_delta)
    candidate=populations
    do population=1,3
      population_options=options
      ! Root has no ghosts, but use retained routing for its allocated empty set.
      ! Pinned SprayParticles.cpp gates splash by isActive and do_move.
      if(population==1.and.drift_dt>0.and.present(wall_films)) then
      call advance_mpi_amr_population_halo(comm,species,transport,components,population_options,tree,owners,halo,halo_t, &
        width,target_level,grow,maximum_parcels,candidate(population)%nodes,candidate(population)%parcels, &
        candidate(population)%component_mass,lower,spacing,periodic,kick_dt,drift_dt,particle_faces,part,part_leftover, &
        wall,exterior,valid,volume,nodal_phi,centroids,geometry,maximum_received,population<3.or.target_level==0, &
        containing_cell=.true.,physical_source=physical_delta,last_id=watermark, &
        wall_films=patches,film_mass=patch_mass,film_settings=film_settings,model_rejected=model_rejected,frame=frame, &
        acceleration=acceleration,force_change=part_force,tab_interval=tab_interval, &
        etab_births=population==1.and.drift_dt>0,khrt_births=population==1.and.drift_dt>0,flow_cfl=flow_cfl, &
        frame_boundary_change=part_boundary_conversion, &
        internal_energy_halo=held_halo)
      else
      call advance_mpi_amr_population_halo(comm,species,transport,components,population_options,tree,owners,halo,halo_t, &
        width,target_level,grow,maximum_parcels,candidate(population)%nodes,candidate(population)%parcels, &
        candidate(population)%component_mass,lower,spacing,periodic,kick_dt,drift_dt,particle_faces,part,part_leftover, &
        wall,exterior,valid,volume,nodal_phi,centroids,geometry,maximum_received,population<3.or.target_level==0, &
        containing_cell=.true.,physical_source=physical_delta,last_id=watermark,frame=frame, &
        acceleration=acceleration,force_change=part_force,tab_interval=tab_interval, &
        etab_births=population==1.and.drift_dt>0,khrt_births=population==1.and.drift_dt>0,flow_cfl=flow_cfl, &
        frame_boundary_change=part_boundary_conversion, &
        internal_energy_halo=held_halo)
      end if
      if(.not.valid) return
      if(present(frame)) then
        call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,part,lower,spacing,frame,.true.,valid,geometry)
        if(.not.valid) return
        call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,part_leftover,lower,spacing, &
          frame,.true.,valid,geometry,width)
        if(.not.valid) return
        if(allocated(physical_delta).and.population==1) then
          call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,physical_delta,lower,spacing, &
            frame,.true.,valid,geometry,energy_change=physical_conversion)
          if(.not.valid) return
        end if
      end if
      if(population==1) then
        physical_force=part_force;physical_boundary_conversion=part_boundary_conversion
        if(present(physical_source)) exported=physical_delta
        if(allocated(physical_delta)) deallocate(physical_delta)
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
    if(present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    if(present(last_id)) last_id=watermark
    if(present(force_change)) force_change=physical_force
    if(present(frame_source_change)) frame_source_change=physical_conversion
    if(present(frame_boundary_change)) frame_boundary_change=physical_boundary_conversion
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
      source,remainder,wall_change,boundary_change,ok,fluid_fraction,nodal_phi,centroids,eb_geometry,maximum_received,flow_cfl, &
      internal_energy_halo)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_halo

    real(dp), optional, intent(in) :: flow_cfl
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
      fluid_fraction,nodal_phi,centroids,eb_geometry,maximum_received,.false.,flow_cfl=flow_cfl, &
      internal_energy_halo=internal_energy_halo)
  end subroutine advance_mpi_amr_ghost_halo

  ! Physical stage: reject loss outside the level grow region after kick/drift.
  subroutine advance_mpi_amr_level_halo(comm,species,transport,components,options,tree,owners,gas_halo,temperature_halo, &
      width,target_level,grow,maximum_parcels,nodes,parcels,component_mass,lower,spacing,periodic,kick_dt,drift_dt,faces, &
      source,remainder,wall_change,boundary_change,ok,fluid_fraction,nodal_phi,centroids,eb_geometry,maximum_received,flow_cfl, &
      internal_energy_halo)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_halo

    real(dp), optional, intent(in) :: flow_cfl
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
      fluid_fraction,nodal_phi,centroids,eb_geometry,maximum_received,.true.,flow_cfl=flow_cfl, &
      internal_energy_halo=internal_energy_halo)
  end subroutine advance_mpi_amr_level_halo

  subroutine advance_mpi_amr_population_halo(comm,species,transport,components,options,tree,owners,gas_halo,temperature_halo, &
      width,target_level,grow,maximum_ghosts,nodes,ghosts,component_mass,lower,spacing,periodic,kick_dt,drift_dt,faces, &
      source,remainder,wall_change,boundary_change,ok,fluid_fraction,nodal_phi,centroids,eb_geometry, &
      maximum_received,keep_physical,containing_cell,physical_source,last_id,wall_films,film_mass,film_settings, &
      model_rejected,frame,acceleration,force_change,tab_interval,etab_births,khrt_births,flow_cfl,frame_boundary_change, &
      internal_energy_halo)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_halo

    real(dp), optional, intent(inout) :: frame_boundary_change
    real(dp) :: pending_boundary_conversion,total_boundary_conversion
    real(dp), optional, intent(in) :: flow_cfl
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(in) :: tab_interval
    logical, optional, intent(in) :: etab_births,khrt_births
    real(dp), optional, intent(out) :: force_change(5)
    integer(int64), optional, intent(inout) :: last_id
    integer(int64), allocatable :: watermark
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    logical, optional, intent(out) :: model_rejected
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    integer :: film_flag,film_low,film_high
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
    type(amr_tree_state_3d) :: raw,routed,leftover,total_source,total_remainder,total_physical
    real(dp), allocatable :: history(:)
    integer :: cycle,cycles,node
    real(dp) :: intervals(4),low_intervals(4),high_intervals(4)
    real(dp) :: slice_kick,slice_drift,slice_tab,total_wall(5),total_exterior(reactive_nvar(size(species)))
    real(dp) :: total_force(5)
    real(dp) :: wall(5),exterior(reactive_nvar(size(species)))
    integer :: status
    logical :: valid,agreed
    real(dp) :: pending_force(5)
    if(present(force_change)) force_change=0
    if(present(last_id)) then
      allocate(watermark);watermark=last_id
    end if
    ok=.false.
    if(present(model_rejected)) model_rejected=.false.
    film_flag=merge(1,0,present(wall_films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))
    call MPI_Allreduce(film_flag,film_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,film_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.film_low/=film_high) return
    if(film_flag/=0.and.film_flag/=7) return
    if(present(wall_films)) then
      valid=allocated(wall_films).and.allocated(film_mass).and.present(last_id)
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS.or..not.agreed) return
      patches=wall_films;patch_mass=film_mass
    end if
    valid=allocated(ghosts).and.allocated(component_mass).and.allocated(nodes)
    valid=valid.and.size(boundary_change)==size(exterior)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    if(present(physical_source)) allocate(physical_delta)
    candidate=ghosts;masses=component_mass;locations=nodes
    ! Explicit native flow subcycles use a frozen stage gas field.
    ! Keep all particle, film and source updates private until every slice succeeds.
    intervals=[kick_dt,drift_dt,kick_dt,-1._dp]
    if(present(flow_cfl)) intervals(4)=flow_cfl
    if(present(tab_interval)) intervals(3)=tab_interval
    valid=all(ieee_is_finite(intervals))
    if(valid) valid=all(intervals(:3)>=0).and.all(intervals(:3)<=1e6_dp).and.intervals(4)>=-1
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(intervals,low_intervals,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(intervals,high_intervals,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low_intervals/=high_intervals)) return
    ! PeleC omits spray_cfl_lev at all three population calls, so the
    ! PelePhysics default (-1) disables subdivision. Honor a supplied CFL
    ! exactly as moveKickDrift does, and keep moveKick at one iteration.
    cycles=1
    valid=options%maximum_substeps>=1
    if(drift_dt>0.and.intervals(4)>.5_dp) then
      valid=valid.and.intervals(4)<=.5_dp*real(min(options%maximum_substeps,huge(cycles)-1),dp)
      if(valid) cycles=ceiling(intervals(4)/.5_dp)
    endif
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    slice_kick=kick_dt/real(cycles,dp);slice_drift=drift_dt/real(cycles,dp)
    slice_tab=intervals(3)/real(cycles,dp)
    total_wall=0;total_exterior=0;total_force=0;total_boundary_conversion=0
    if(options%etab_parcel_factor>=0) allocate(history(size(candidate)),source=0._dp)
    do cycle=1,cycles
      if(allocated(history)) then
        call advance_slice(history)
      else
        call advance_slice()
      end if
      if(.not.valid) return
      if(cycle==1) then
        total_source=routed;total_remainder=leftover
        if(present(physical_source)) total_physical=physical_delta
      else
        do node=1,size(owners)
          if(allocated(routed%fields(node)%values)) &
            total_source%fields(node)%values=total_source%fields(node)%values+routed%fields(node)%values
          if(allocated(leftover%fields(node)%values)) &
            total_remainder%fields(node)%values=total_remainder%fields(node)%values+leftover%fields(node)%values
          if(present(physical_source)) then
            if(allocated(physical_delta%fields(node)%values)) &
              total_physical%fields(node)%values=total_physical%fields(node)%values+physical_delta%fields(node)%values
          end if
        end do
      end if
      total_wall=total_wall+wall;total_exterior=total_exterior+exterior;total_force=total_force+pending_force
      total_boundary_conversion=total_boundary_conversion+pending_boundary_conversion
    end do
    valid=all(ieee_is_finite(total_wall)).and.all(ieee_is_finite(total_exterior)).and. &
      all(ieee_is_finite(total_force)).and.ieee_is_finite(total_boundary_conversion)
    do node=1,size(owners)
      if(allocated(total_source%fields(node)%values)) &
        valid=valid.and.all(ieee_is_finite(total_source%fields(node)%values))
      if(allocated(total_remainder%fields(node)%values)) &
        valid=valid.and.all(ieee_is_finite(total_remainder%fields(node)%values))
      if(present(physical_source)) then
        if(allocated(total_physical%fields(node)%values)) &
          valid=valid.and.all(ieee_is_finite(total_physical%fields(node)%values))
      endif
    enddo
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    if(present(force_change)) force_change=total_force
    if(present(frame_boundary_change)) frame_boundary_change=total_boundary_conversion
    call move_alloc(candidate,ghosts);call move_alloc(masses,component_mass);call move_alloc(locations,nodes)
    call move_alloc(total_source%fields,source%fields);call move_alloc(total_remainder%fields,remainder%fields)
    if(present(physical_source)) physical_source=total_physical
    if(present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    if(present(last_id)) last_id=watermark
    wall_change=total_wall;boundary_change=total_exterior;ok=.true.
  contains
    subroutine advance_slice(tangent)
      real(dp), allocatable, optional, intent(inout) :: tangent(:)
    if(present(wall_films)) then
    call kick_mpi_amr_mixture_halo(comm,species,transport,components,options,tree,owners, &
      gas_halo,temperature_halo,width,locations,candidate,masses,lower,spacing,periodic,slice_kick,raw,valid, &
      fluid_fraction,nodal_phi,centroids,slice_drift,faces,wall,exterior,target_level,routed,leftover,.false., &
      eb_geometry,maximum_received,containing_cell,physical_delta,watermark, &
      patches,patch_mass,film_settings,model_rejected,frame, &
        acceleration=acceleration,force_change=pending_force,tab_interval=slice_tab, &
        etab_births=etab_births,khrt_births=khrt_births,tab_tangent_history=tangent, &
        last_flow_cycle=cycle==cycles,khrt_interval=kick_dt,frame_boundary_change=pending_boundary_conversion, &
        internal_energy_halo=internal_energy_halo)
    else
    call kick_mpi_amr_mixture_halo(comm,species,transport,components,options,tree,owners, &
      gas_halo,temperature_halo,width,locations,candidate,masses,lower,spacing,periodic,slice_kick,raw,valid, &
      fluid_fraction,nodal_phi,centroids,slice_drift,faces,wall,exterior,target_level,routed,leftover,.false., &
      eb_geometry,maximum_received,containing_cell,physical_delta,watermark,frame=frame, &
        acceleration=acceleration,force_change=pending_force,tab_interval=slice_tab, &
        etab_births=etab_births,khrt_births=khrt_births,tab_tangent_history=tangent, &
        last_flow_cycle=cycle==cycles,khrt_interval=kick_dt,frame_boundary_change=pending_boundary_conversion, &
        internal_energy_halo=internal_energy_halo)
    end if
    if(.not.valid) return
    if(keep_physical) then
      call redistribute_mpi_amr_level_particles(comm,tree,owners,candidate,masses,locations,lower,spacing, &
        periodic,target_level,grow,maximum_ghosts,valid,nodal_phi,tangent)
    else
    call redistribute_mpi_amr_ghost_particles(comm,tree,owners,candidate,masses,locations,lower,spacing, &
      periodic,target_level,grow,maximum_ghosts,valid,nodal_phi,tangent)
    end if
    if(.not.valid) return
    end subroutine advance_slice
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
      containing_cell,physical_source,last_id,wall_films,film_mass,film_settings,model_rejected, &
      frame,acceleration,force_change,tab_interval,etab_births,khrt_births, &
      tab_tangent_history,last_flow_cycle,khrt_interval,frame_boundary_change, &
      internal_energy_halo)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_halo
    real(dp), optional, intent(inout) :: frame_boundary_change
    real(dp) :: boundary_conversion,global_boundary_conversion
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(in) :: tab_interval
    logical, optional, intent(in) :: etab_births,khrt_births,last_flow_cycle
    real(dp), allocatable, optional, intent(inout) :: tab_tangent_history(:)
    real(dp), optional, intent(in) :: khrt_interval
    real(dp), allocatable :: held_cell
    real(dp), allocatable :: next_tangent_history(:)
    real(dp) :: kh_dt,kh_low,kh_high
    logical :: finish_flow,history_sorted
    integer(int64) :: history_max_id
    real(dp), optional, intent(out) :: force_change(5)
    integer(int64), optional, intent(inout) :: last_id
    integer(int64) :: watermark
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    logical, optional, intent(out) :: model_rejected
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:),pressures(:),impact_pressure(:),impact_previous(:,:)
    type(spray_eb_wall_options) :: impact_settings
    integer :: film_flag,film_low,film_high,parent,lo,hi,mid,pass
    integer(int64) :: impact_parent_id
    logical :: rejected
    logical, allocatable :: planar_hits(:)
    real(dp) :: planar_change(5),potential
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    real(dp) :: mw(size(components)),impact_exterior(5+size(components))
    type(spray_options), intent(in) :: options
    type(spray_options) :: parcel_breakup_options
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
    real(dp) :: drift_dt,drift_low,drift_high,wall(5),global_wall(5),frame_data(4),frame_min(4),frame_max(4)
    real(dp) :: exterior(reactive_nvar(size(species))),global_exterior(reactive_nvar(size(species)))
    real(dp) :: exterior5(5),upper(3),empty_phi(0:1,0:1,0:1),rho,cp,energy
    real(dp), allocatable :: energies(:)
    integer :: flags,low_flags,high_flags,hits,c,frame_flag,frame_low,frame_high,axis_min,axis_max
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
    type(amr_tree_state_3d) :: source,routed,remainder,primitive_halo,primitive_temperature_halo
    integer :: transfer_flags,transfer_low,transfer_high,node,rank,n(3)
    real(dp), allocatable :: masses(:,:),positions(:,:),sampled(:,:),sampled_t(:,:),change(:,:)
    real(dp) :: primitive(reactive_nvar(size(species))),tg,sound,low_dt,high_dt
    real(dp), allocatable :: gas_properties(:,:),reynolds(:),child_momentum(:,:),child_energy(:)
    real(dp), allocatable :: etab_tangent(:),etab_positions(:,:),etab_feedback(:,:)
    real(dp) :: etab_delta
    logical, allocatable :: eligible(:)
    type(spray_parcel) :: preview(1)
    real(dp) :: preview_wall(5),preview_exterior(5)
    logical :: contact(1)
    type(spray_liquid) :: effective
    type(spray_film) :: film
    real(dp) :: old_mass,fractions(size(components)),rates(size(components)),hv(size(components))
    real(dp) :: tab_dt,tab_low,tab_high,breakup_dt
    integer :: p,nv,nc,status,ix,iy,iz
    logical :: valid,agreed
    real(dp) :: accel(3),accel_low(3),accel_high(3),local_force(5),global_force(5),parcel_force(5)
    integer :: force_flag,force_low,force_high
    logical :: allow_etab_births,allow_khrt_births
    ok=.false.;nv=reactive_nvar(size(species));nc=size(components)
    if(present(force_change)) force_change=0
    accel=0;local_force=0;global_force=0;boundary_conversion=0;global_boundary_conversion=0
    if(present(acceleration)) accel=acceleration
    allow_khrt_births=.true.
    if(present(khrt_births)) allow_khrt_births=khrt_births
    allow_etab_births=.true.
    if(present(etab_births)) allow_etab_births=etab_births
    finish_flow=.true.
    if(present(last_flow_cycle)) finish_flow=last_flow_cycle
    force_flag=merge(1,0,present(acceleration))+2*merge(1,0,present(force_change))+ &
      4*merge(1,0,present(tab_interval))+8*merge(1,0,present(etab_births))+16*merge(1,0,allow_etab_births)+ &
      32*merge(1,0,present(khrt_births))+64*merge(1,0,allow_khrt_births)+ &
      128*merge(1,0,present(tab_tangent_history))+256*merge(1,0,present(last_flow_cycle))+ &
      512*merge(1,0,finish_flow)+1024*merge(1,0,present(khrt_interval))
    call MPI_Allreduce(force_flag,force_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(force_flag,force_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.force_low/=force_high) return
    tab_dt=dt
    if(present(tab_interval)) tab_dt=tab_interval
    valid=ieee_is_finite(tab_dt)
    if(valid) valid=tab_dt>=0.and.tab_dt<=1e6_dp
    call agree();if(.not.agreed) return
    call MPI_Allreduce(tab_dt,tab_low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(tab_dt,tab_high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.tab_low/=tab_high) return
    kh_dt=dt
    if(present(khrt_interval)) kh_dt=khrt_interval
    valid=ieee_is_finite(kh_dt)
    if(valid) valid=kh_dt>=0.and.kh_dt<=1e6_dp
    call agree();if(.not.agreed) return
    call MPI_Allreduce(kh_dt,kh_low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(kh_dt,kh_high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.kh_low/=kh_high) return
    ! TABBreakup.H evolves the supplied flow subinterval directly. KHRT's
    ! caller applies dtmod explicitly, so retain the kick interval for KHRT.
    breakup_dt=kh_dt
    if(options%breakup_model=='etab') breakup_dt=tab_dt
    valid=all(ieee_is_finite(accel)).and.all(abs(accel)<=1e6_dp)
    if(any(accel/=0)) valid=valid.and.present(force_change).and. &
      (options%mixture_time_scheme=='euler'.or.options%fixed_parts)
    call agree();if(.not.agreed) return
    call MPI_Allreduce(accel,accel_low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(accel,accel_high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(accel_low/=accel_high)) return
    frame_flag=merge(1,0,present(frame))
    call MPI_Allreduce(frame_flag,frame_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(frame_flag,frame_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.frame_low/=frame_high) return
    if(present(frame)) then
      valid=frame%is_valid();call agree();if(.not.agreed) return
      call MPI_Allreduce(frame%axis,axis_min,1,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame%axis,axis_max,1,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.axis_min/=axis_max) return
      frame_data=[frame%angular_speed,frame%axis_point]
      call MPI_Allreduce(frame_data,frame_min,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame_data,frame_max,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(frame_min/=frame_max)) return
    end if
    if(present(model_rejected)) model_rejected=.false.
    film_flag=merge(1,0,present(wall_films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))
    call MPI_Allreduce(film_flag,film_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,film_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.film_low/=film_high) return
    if(film_flag/=0.and.film_flag/=7) return
    if(present(wall_films)) then
      valid=allocated(wall_films).and.allocated(film_mass).and.present(last_id).and. &
        present(drift_interval).and.(present(nodal_phi).or.present(boundary_faces))
      call agree();if(.not.agreed) return
      if(present(nodal_phi)) then
        valid=valid_wall_options(film_settings,[character(len=16):: &
          'slip_wall','slip_wall','slip_wall','slip_wall','slip_wall','slip_wall'])
        valid=valid.and.all(film_settings%temperature==film_settings%temperature(1))
      else
        valid=valid_wall_options(film_settings,boundary_faces)
      endif
      call agree();if(.not.agreed) return
      patches=wall_films;patch_mass=film_mass
      impact_settings%temperature=film_settings%temperature(1)
      impact_settings%contact_angle=film_settings%contact_angle
      impact_settings%maximum_films=film_settings%maximum_films
    end if
    valid=allocated(parcels).and.allocated(component_mass).and.valid_spray_options(options)
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    valid=valid.and.(options%etab_parcel_factor<0.or.present(last_id))
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
      if(options%fixed_parts) drift_dt=0
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
    export_flag=merge(1,0,present(physical_source))+2*merge(1,0,present(frame_boundary_change))
    call MPI_Allreduce(export_flag,export_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(export_flag,export_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.export_low/=export_high) return
    if(present(tab_tangent_history)) then
      valid=allocated(tab_tangent_history).and.options%etab_parcel_factor>=0
      call agree();if(.not.agreed) return
      valid=size(tab_tangent_history)==size(parcels)
      if(valid) valid=all(ieee_is_finite(tab_tangent_history)).and.all(tab_tangent_history>=0)
      call agree();if(.not.agreed) return
    endif
    history_max_id=maxval(parcels%id)
    history_sorted=.true.
    do p=2,size(parcels)
      if(parcels(p)%id<=parcels(p-1)%id) history_sorted=.false.
    enddo
    allocate(positions(3,size(parcels)),change(nv,size(parcels)))
    do p=1,size(parcels)
      positions(:,p)=parcels(p)%position
    end do
    call validate_mpi_spray_internal_energy_halo(comm,owners,gas_halo,internal_energy_halo,valid)
    if(.not.valid) return
    if(present(internal_energy_halo)) then
      allocate(held_cell,stat=status)
      call MPI_Allreduce(status==0,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS.or..not.agreed) return
    endif
    call sample_mpi_amr_particle_halo(comm,tree,owners,gas_halo,nv,width,nodes,positions,lower,spacing, &
      periodic,sampled,valid,fluid_fraction,nodal_phi,centroids)
    if(.not.valid) return
    call sample_mpi_amr_particle_halo(comm,tree,owners,temperature_halo,1,width,nodes,positions,lower,spacing, &
      periodic,sampled_t,valid,fluid_fraction,nodal_phi,centroids)
    if(.not.valid) return
    ! The first interpolation validates owner-local halo layouts and geometry.
    ! Pinned InterpolateGasPhase weights cell primitive rho/u/Y/T, rather
    ! than recovering temperature from weighted conserved kinetic energy.
    primitive_halo=gas_halo;primitive_temperature_halo=temperature_halo
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    valid=.true.
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      if(.not.allocated(primitive_halo%fields(node)%values)) cycle
      do iz=lbound(primitive_halo%fields(node)%values,4),ubound(primitive_halo%fields(node)%values,4)
      do iy=lbound(primitive_halo%fields(node)%values,3),ubound(primitive_halo%fields(node)%values,3)
      do ix=lbound(primitive_halo%fields(node)%values,2),ubound(primitive_halo%fields(node)%values,2)
        if(present(fluid_fraction)) then
          if(fluid_fraction%fields(node)%values(1,ix,iy,iz)<=0) then
            primitive_halo%fields(node)%values(:,ix,iy,iz)=0
            primitive_temperature_halo%fields(node)%values(1,ix,iy,iz)=0
            cycle
          endif
        endif
        if(present(internal_energy_halo)) held_cell=internal_energy_halo%fields(node)%values(1,ix,iy,iz)
        call reactive_conserved_to_primitive(species,gas_halo%fields(node)%values(:,ix,iy,iz), &
          temperature_halo%fields(node)%values(1,ix,iy,iz),primitive,tg,sound,valid,internal_energy_density=held_cell)
        if(.not.valid) exit
        primitive_halo%fields(node)%values(:,ix,iy,iz)=primitive
        primitive_temperature_halo%fields(node)%values(1,ix,iy,iz)=tg
      enddo
      if(.not.valid) exit
      enddo
      if(.not.valid) exit
      enddo
      if(.not.valid) exit
    enddo
    call agree();if(.not.agreed) return
    call sample_mpi_amr_particle_halo(comm,tree,owners,primitive_halo,nv,width,nodes,positions,lower,spacing, &
      periodic,sampled,valid,fluid_fraction,nodal_phi,centroids)
    if(.not.valid) return
    call sample_mpi_amr_particle_halo(comm,tree,owners,primitive_temperature_halo,1,width,nodes,positions,lower,spacing, &
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
    if(options%etab_parcel_factor>=0) then
      allocate(etab_tangent(size(parcels)));etab_tangent=0
      if(present(tab_tangent_history)) etab_tangent=tab_tangent_history
    endif
    allocate(pressures(size(parcels)),source=1._dp)
    do p=1,size(parcels)
      if(parcels(p)%mass<=0) cycle
      primitive=sampled(:,p);tg=sampled_t(1,p)
      primitive(5)=mixture_pressure(species,primitive(6:),primitive(1),tg,valid)
      if(.not.valid) exit
      pressures(p)=primitive(5)
      if(options%mixture_time_scheme=='euler') then
        ! SDC supplies a half-kick. Drag.H estimates nsub from flow_dt,
        ! then advances each parcel for dtmod*flow_dt with dtmod=1/2.
        call advance_mixture_frozen_gas(species,transport,primitive,tg,candidate(p),components,masses(:,p), &
          options,2*dt,change(:,p),valid,time_modifier=.5_dp,acceleration=accel,force_change=parcel_force, &
          reynolds_number=reynolds(p))
      else
        call advance_mixture_frozen_gas(species,transport,primitive,tg,candidate(p),components,masses(:,p), &
          options,dt,change(:,p),valid,acceleration=accel,force_change=parcel_force,reynolds_number=reynolds(p))
      endif
      if(.not.valid) exit
      local_force=local_force+parcel_force
      valid=all(ieee_is_finite(local_force))
      if(.not.valid) exit
      if(candidate(p)%mass>0.and.options%breakup_model/='none'.and. &
          (options%breakup_model=='etab'.or.finish_flow)) then
        old_mass=candidate(p)%mass;fractions=masses(:,p)/old_mass
        call liquid_mixture_properties(components,fractions,candidate(p)%temperature,rho,cp,energy,valid,primitive(5), &
          species(components%vapor_index)%molecular_weight)
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
        parcel_breakup_options=options
        call apply_mixture_breakup_viscosity(components,fractions,candidate(p)%temperature,parcel_breakup_options, &
          valid,primitive(5),species(components%vapor_index)%molecular_weight)
        if(.not.valid) exit
        if(options%khrt_children) then
          gas_properties(:,p)=primitive(:4)
        else
          if(options%etab_parcel_factor>=0) then
            call advance_spray_breakup(candidate(p),effective,parcel_breakup_options,primitive(1),primitive(2:4), &
              film%density,film%viscosity,breakup_dt,valid,etab_delta,reynolds_number=reynolds(p))
            if(valid) etab_tangent(p)=etab_tangent(p)+etab_delta
            if(valid.and.finish_flow.and..not.allow_etab_births) &
              call advance_tab_no_birth_khrt(candidate(p),effective,parcel_breakup_options,primitive(1),primitive(2:4), &
                reynolds(p),kh_dt,valid)
          else if(finish_flow.or.options%breakup_model=='etab') then
            call advance_spray_breakup(candidate(p),effective,parcel_breakup_options,primitive(1),primitive(2:4), &
              film%density,film%viscosity,breakup_dt,valid,reynolds_number=reynolds(p))
          endif
          if(.not.valid) exit
          masses(:,p)=masses(:,p)*(candidate(p)%mass/old_mass)
        end if
      end if
      candidate(p)%position=parcels(p)%position
    end do
    call agree();if(.not.agreed) return
    if(options%etab_parcel_factor>=0.and.dt>0.and.allow_etab_births.and.finish_flow) then
      call split_mpi_reference_etab_parcels(comm,candidate,masses,options%maximum_multiplicity, &
        min(125000,options%maximum_parcels),options%etab_parcel_factor,etab_tangent,0._dp, &
        watermark,etab_positions,etab_feedback,valid)
      if(.not.valid) return
      change(2:5,:)=change(2:5,:)+etab_feedback
    endif
    if(options%khrt_children.and.finish_flow) then
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
        gas_properties(2:4,:),reynolds,kh_dt,watermark,child_energy,valid,eligible,gas_momentum=child_momentum, &
        pressure=pressures,gas_molecular_weights=species(components%vapor_index)%molecular_weight, &
        child_drift_dt=0._dp,make_new_drops=allow_khrt_births)
      if(.not.valid) return
      valid=all(ieee_is_finite(child_energy))
      if(valid) valid=.not.any(change(5,:)>huge(1._dp)-max(child_energy,0._dp)).and. &
        .not.any(change(5,:)< -huge(1._dp)-min(child_energy,0._dp))
      call agree();if(.not.agreed) return
      change(2:4,:)=change(2:4,:)+child_momentum
      change(5,:)=change(5,:)+child_energy
    end if
    global_wall=0;global_exterior=0
    if(present(wall_films).and..not.options%fixed_parts) then
      ! Shedding children inherit the sampled pressure of their pre-drift parent,
      ! matching the common-step impact preparation. Sources retain parent sites.
      allocate(impact_pressure(size(candidate)))
      do p=1,size(candidate)
        parent=0;impact_parent_id=candidate(p)%id
        do pass=1,2
          lo=1;hi=size(parcels)
          do while(lo<=hi)
            mid=lo+(hi-lo)/2
            if(parcels(mid)%id==impact_parent_id) then
              parent=mid;exit
            else if(parcels(mid)%id<impact_parent_id) then
              lo=mid+1
            else
              hi=mid-1
            end if
          end do
          if(parent>0) exit
          impact_parent_id=candidate(p)%parent_id
        enddo
        if(parent==0) then
          valid=.false.;exit
        endif
        impact_pressure(p)=pressures(parent)
      end do
      call agree();if(.not.agreed) return
      do c=1,nc
        mw(c)=species(components(c)%vapor_index)%molecular_weight
      end do
      upper=lower+spacing*real(tree%nodes(1)%shape,dp)
      global_wall=0;impact_exterior=0;hits=0
      rejected=.false.
      if(present(nodal_phi)) then
      call advance_mpi_level_set_impacts(comm,components,options,impact_settings,nodal_phi,mw,impact_pressure, &
        lower,upper,drift_dt,candidate,masses,patches,patch_mass,watermark,global_wall,hits,valid,rejected, &
        impact_exterior,faces,frame,energy_conversion=global_boundary_conversion)
      else
        allocate(impact_previous(3,size(candidate)))
        do p=1,size(candidate)
          impact_previous(:,p)=candidate(p)%position
          if(candidate(p)%mass>0) candidate(p)%position=candidate(p)%position+drift_dt*candidate(p)%velocity
        enddo
        call impact_mpi_planar_walls(comm,components,options,film_settings,faces,mw,impact_pressure,impact_previous, &
          lower,upper,candidate,masses,patches,patch_mass,watermark,planar_hits,global_wall,valid)
        if(.not.valid) return
        ! The impact kernel resolves enabled walls. Apply periodic/outflow and
        ! remaining reflecting faces to the same completed flight.
        wall=0;exterior=0
        do p=1,size(candidate)
          if(candidate(p)%mass<=0) cycle
          fractions=masses(:,p)/candidate(p)%mass
          call liquid_mixture_properties(components,fractions,candidate(p)%temperature,rho,cp,energy,valid)
          if(.not.valid) exit
          effective=components(1);effective%density=rho;effective%density_model='constant'
          effective%density_coefficients=0;effective%density_temperature_range=0
          effective%heat_capacity=cp;effective%heat_capacity_model='constant'
          effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
          effective%reference_temperature=candidate(p)%temperature;effective%reference_energy=energy
          effective%boiling_model='reject';effective%critical_temperature=0;effective%normal_boiling_temperature=0
          call apply_parcel_boundaries(candidate(p),effective,lower,upper,faces,planar_change,valid)
          if(.not.valid) exit
          if(planar_change(1)<0) then
            if(present(frame)) then
              call rotating_potential_3d(frame,candidate(p)%position,potential,valid)
              if(.not.valid) exit
              planar_change(5)=planar_change(5)-planar_change(1)*potential
              boundary_conversion=boundary_conversion-planar_change(1)*potential
            endif
            exterior(:5)=exterior(:5)+planar_change
            do c=1,nc
              exterior(5+components(c)%vapor_index)=exterior(5+components(c)%vapor_index)+planar_change(1)*fractions(c)
            enddo
            masses(:,p)=0
          else
            wall=wall+planar_change
          endif
        enddo
        call agree();if(.not.agreed) return
        call MPI_Allreduce(wall,planar_change,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
        if(status/=MPI_SUCCESS) return
        ! Convert extra gas-side reflecting inventory back to wall-received form.
        global_wall=global_wall-planar_change
        call MPI_Allreduce(exterior,global_exterior,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
        if(status/=MPI_SUCCESS) return
      endif
      if(rejected.and.present(model_rejected)) model_rejected=.true.
      if(.not.valid) return
      ! Impact API reports inventory received by the wall, already MPI summed.
      global_wall=-global_wall
      if(present(nodal_phi)) then
        global_exterior=0;global_exterior(:5)=impact_exterior(:5)
        do c=1,nc
          global_exterior(5+components(c)%vapor_index)=global_exterior(5+components(c)%vapor_index)+impact_exterior(5+c)
        end do
      endif
    else if(present(drift_interval).and..not.options%fixed_parts) then
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
          faces,energies,exterior5,frame=frame,energy_conversion=boundary_conversion)
      else
        call drift_spray_level_set(empty_phi,lower,upper,drift_dt,32,candidate,wall,hits,valid, &
          faces,energies,exterior5,frame=frame,energy_conversion=boundary_conversion)
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
    call MPI_Allreduce(boundary_conversion,potential,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if(status/=MPI_SUCCESS) return
    global_boundary_conversion=global_boundary_conversion+potential
    valid=ieee_is_finite(global_boundary_conversion);call agree();if(.not.agreed) return
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
    call MPI_Allreduce(local_force,global_force,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if(status/=MPI_SUCCESS) return
    valid=all(ieee_is_finite(global_force));call agree();if(.not.agreed) return
    if(present(tab_tangent_history)) then
      allocate(next_tangent_history(size(candidate)),source=0._dp)
      do p=1,size(candidate)
        if(p<=size(parcels)) then
          if(candidate(p)%id==parcels(p)%id) then
            next_tangent_history(p)=etab_tangent(p);cycle
          endif
        endif
        ! Replacement children start a new tangent history. Retained original
        ! parcels keep their accumulated value after parent removal/reordering.
        if(candidate(p)%id>history_max_id) cycle
        if(history_sorted) then
          lo=1;hi=size(parcels)
          do while(lo<=hi)
            mid=lo+(hi-lo)/2
            if(candidate(p)%id==parcels(mid)%id) then
              next_tangent_history(p)=etab_tangent(mid);exit
            else if(candidate(p)%id>parcels(mid)%id) then
              lo=mid+1
            else
              hi=mid-1
            endif
          enddo
        else
          do parent=1,size(parcels)
            if(candidate(p)%id/=parcels(parent)%id) cycle
            next_tangent_history(p)=etab_tangent(parent);exit
          enddo
        endif
      enddo
      valid=all(ieee_is_finite(next_tangent_history))
      call agree();if(.not.agreed) return
    endif
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
    if(present(force_change)) force_change=global_force
    if(present(frame_boundary_change)) frame_boundary_change=global_boundary_conversion
    if(present(last_id)) last_id=watermark
    if(present(tab_tangent_history)) call move_alloc(next_tangent_history,tab_tangent_history)
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
      motion_mode,maximum_level,containing_cell,frame,acceleration,force_change,native_flow_cfl)
    real(dp), optional, intent(in) :: native_flow_cfl
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(out) :: force_change(5)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(amr_tree_state_3d), optional, intent(inout) :: gas_change
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d) :: integrated_change
    type(amr_tree_state_3d), allocatable :: trial_change,etab_trial_change
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
    logical, optional, intent(in) :: containing_cell
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
    real(dp), allocatable :: reynolds(:),child_momentum(:,:),child_energy(:)
    real(dp), allocatable :: etab_tangent(:),etab_positions(:,:),etab_feedback(:,:),etab_change(:,:)
    real(dp) :: etab_delta,birth_drift,breakup_dt,flow_cfl,global_flow_cfl,flow_slice
    real(dp), allocatable :: committed_tangent(:)
    integer(int64), allocatable :: tangent_ids(:)
    integer :: flow_cycles,tangent_index
    logical :: native_flow,last_flow_cycle
    integer :: child_flag,low_flag,high_flag
    type(spray_parcel), allocatable :: drops(:),trial(:),survivors(:)
    type(spray_options) :: breakup_options
    type(spray_liquid) :: effective
    type(spray_film) :: film
    real(dp) :: old_mass
    type(amr_tree_state_3d) :: gas,temp,next_gas,next_temp,field,stage_gas,stage_temp
    real(dp), allocatable :: masses(:,:),next_mass(:,:),kept_mass(:,:),positions(:,:),samples(:,:)
    real(dp), allocatable :: primitive(:,:),tg(:),change(:,:)
    real(dp) :: elapsed,dt,local_dt,min_interval,max_interval,sound,invtau,heat,rho,cp,e,tdot,speed,potential
    real(dp) :: accel(3),accel_low(3),accel_high(3),local_force(5),global_force(5),trial_force(5),parcel_force(5)
    integer :: force_flag,force_low,force_high
    real(dp) :: surface_energy
    real(dp) :: y(size(components)),rates(size(components)),hv(size(components)),length(3)
    integer :: p,n,np,nvar,nc,ierr,attempt,count_steps,kept,min_nc,max_nc,node,shape3(3),frame_flag,frame_low,frame_high
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: minimum_spacing,sample_temperature,sample_sound
    real(dp) :: sample_primitive(5+size(species))
    integer :: ix,iy,iz
    logical :: periodic(3)
    logical :: valid,agreed,accepted
    ok=.false.;substeps=0;faces='periodic';ledger=0;wall_ledger=0
    if(present(force_change)) force_change=0
    accel=0;local_force=0;global_force=0
    if(present(acceleration)) accel=acceleration
    force_flag=merge(1,0,present(acceleration))+2*merge(1,0,present(force_change))
    call MPI_Allreduce(force_flag,force_low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(force_flag,force_high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.force_low/=force_high) return
    valid=all(ieee_is_finite(accel)).and.all(abs(accel)<=1e6_dp)
    if(any(accel/=0)) valid=valid.and.present(force_change).and. &
      (options%mixture_time_scheme=='euler'.or.options%fixed_parts)
    call agree();if(.not.agreed) return
    call MPI_Allreduce(accel,accel_low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(accel,accel_high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(accel_low/=accel_high)) return
    frame_flag=merge(1,0,present(frame))
    call MPI_Allreduce(frame_flag,frame_low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(frame_flag,frame_high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.frame_low/=frame_high) return
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
    breakup_options=options
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
    valid=valid.and.(.not.breakup_options%khrt_children.or.present(last_id))
    valid=valid.and.(breakup_options%etab_parcel_factor<0.or.present(last_id))
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
    child_flag=merge(1,0,breakup_options%khrt_children)
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
    if(present(frame)) then
      call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,gas,lower,spacing,frame,.false.,valid,eb_geometry)
      if(.not.valid) return
    end if
    if(motion/='coupled') then
      stage_gas=gas;stage_temp=temp
    endif
    if (present(wall_films)) then
      patches=wall_films;patch_mass=film_mass
      call migrate_mpi_amr_wall_films(comm,tree,owners,patches,patch_mass,lower,spacing, &
        film_settings%maximum_films,valid,nodal_phi,maximum_level,periodic=periodic)
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
      lower,spacing,periodic,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi,maximum_level=maximum_level, &
      containing_cell=containing_cell)
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
    native_flow=motion/='coupled';flow_cycles=1;flow_cfl=-1
    if(present(native_flow_cfl)) flow_cfl=native_flow_cfl
    valid=ieee_is_finite(flow_cfl)
    if(valid) valid=flow_cfl>=-1
    call agree();if(.not.agreed) return
    call MPI_Allreduce(flow_cfl,global_flow_cfl,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=global_flow_cfl==flow_cfl
    call agree();if(.not.agreed) return
    if(native_flow.and.motion=='kick_drift'.and.flow_cfl>.5_dp) then
      valid=flow_cfl<=.5_dp*real(min(options%maximum_substeps,huge(flow_cycles)-1),dp)
      call agree();if(.not.agreed) return
      flow_cycles=ceiling(flow_cfl/.5_dp)
    endif
    flow_slice=interval/real(flow_cycles,dp)
    allocate(tangent_ids(size(drops)),committed_tangent(size(drops)))
    tangent_ids=drops%id;committed_tangent=0
    elapsed=0;count_steps=0
    do while((native_flow.and.count_steps<flow_cycles).or.(.not.native_flow.and.elapsed<interval))
      valid=count_steps<options%maximum_substeps
      call agree()
      if (.not.agreed) return
      np=size(drops)
      deallocate(positions,change,primitive,tg)
      if (allocated(reynolds)) deallocate(reynolds)
      allocate(positions(3,np),change(nvar,np),primitive(nvar,np),tg(np),reynolds(np))
      primitive=0;tg=0
      ! Native MK/MKD samples the supplied gas stage, frozen across particle
      ! slices. InterpolateGasPhase weights cell primitives and temperature.
      if(.not.native_flow.or.count_steps==0) then
      if (allocated(field%fields)) deallocate(field%fields)
      allocate(field%fields(size(tree%nodes)))
      do node=1,size(tree%nodes)
        if (.not.allocated(gas%fields(node)%values)) cycle
        shape3=tree%nodes(node)%shape
        allocate(field%fields(node)%values(nvar+1,shape3(1),shape3(2),shape3(3)))
        field%fields(node)%values(:nvar,:,:,:)=gas%fields(node)%values
        field%fields(node)%values(nvar+1,:,:,:)=temp%fields(node)%values(1,:,:,:)
        if(native_flow) then
          field%fields(node)%values(:nvar,:,:,:)=stage_gas%fields(node)%values
          field%fields(node)%values(nvar+1,:,:,:)=stage_temp%fields(node)%values(1,:,:,:)
          do iz=1,shape3(3)
          do iy=1,shape3(2)
          do ix=1,shape3(1)
            if(present(eb_geometry)) then
              if(eb_geometry(node)%volume_fraction(ix,iy,iz)<=0) then
                field%fields(node)%values(:,ix,iy,iz)=0;cycle
              endif
            endif
            call reactive_conserved_to_primitive(species,stage_gas%fields(node)%values(:,ix,iy,iz), &
              stage_temp%fields(node)%values(1,ix,iy,iz),sample_primitive,sample_temperature,sample_sound,valid)
            if(.not.valid) exit
            field%fields(node)%values(:nvar,ix,iy,iz)=sample_primitive
            field%fields(node)%values(nvar+1,ix,iy,iz)=sample_temperature
          enddo
          if(.not.valid) exit
          enddo
          if(.not.valid) exit
          enddo
          if(.not.valid) exit
        endif
      end do
      call agree();if(.not.agreed) return
      endif
      do p=1,np
        positions(:,p)=drops(p)%position
      end do
      call sample_mpi_amr_tree_spray_field(comm,tree,owners,field,nvar+1,positions, &
        lower,spacing,periodic,samples,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi,maximum_level=maximum_level)
      if (.not.valid) return
      local_dt=interval-elapsed
      do p=1,np
        if (drops(p)%mass<=0) cycle
        if(native_flow) then
          primitive(:,p)=samples(:nvar,p);tg(p)=samples(nvar+1,p)
          primitive(5,p)=mixture_pressure(species,primitive(6:,p),primitive(1,p),tg(p),valid)
        else
          call reactive_conserved_to_primitive(species,samples(:nvar,p),samples(nvar+1,p),primitive(:,p),tg(p),sound,valid)
        endif
        if (.not.valid) exit
        y=masses(:,p)/drops(p)%mass
        call mixture_transfer_rates(species,transport,primitive(:,p),tg(p),drops(p),components,y,options, &
          invtau,heat,rates,hv,valid)
        if (.not.valid) exit
        call liquid_mixture_properties(components,y,drops(p)%temperature,rho,cp,e,valid,primitive(5,p), &
          species(components%vapor_index)%molecular_weight)
        if (.not.valid) exit
        if(.not.native_flow) then
        if (invtau>0) local_dt=min(local_dt,options%max_fractional_change/invtau/kick_fraction)
        tdot=heat*(tg(p)-drops(p)%temperature)
        if(options%abramzon_sirignano.and.options%evaporation) tdot=heat*max(0._dp,tg(p)-drops(p)%temperature)
        do n=1,nc
          if (rates(n)>0.and.y(n)>mixture_trace_fraction) &
            local_dt=min(local_dt,options%max_fractional_change*masses(n,p)/rates(n)/kick_fraction)
          if(rates(n)==0) cycle
          surface_energy=liquid_energy(components(n),drops(p)%temperature)
          if(options%mixture_time_scheme=='euler') then
            call liquid_evaporation_energy(components(n),drops(p)%temperature,surface_energy,valid, &
              primitive(5,p),species(components(n)%vapor_index)%molecular_weight)
            if(.not.valid) exit
          endif
          tdot=tdot-rates(n)*(hv(n)-surface_energy)
        end do
        if(.not.valid) exit
        tdot=tdot/(drops(p)%mass*cp)
        if (abs(tdot)>0) local_dt=min(local_dt, &
          min(2.0_dp,options%max_fractional_change*drops(p)%temperature)/abs(tdot)/kick_fraction)
        speed=max(maxval(abs(drops(p)%velocity)),maxval(abs(primitive(2:4,p))))
        if (speed>0.and.motion/='kick'.and..not.options%fixed_parts) &
          local_dt=min(local_dt,options%particle_cfl*minimum_spacing/speed)
        endif
      end do
      call agree()
      if (.not.agreed) return
      if(native_flow) local_dt=flow_slice
      call MPI_Allreduce(local_dt,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      accepted=.false.
      do attempt=1,30
        if (dt<=0.or.elapsed+dt<=elapsed) exit
        kick_dt=kick_fraction*dt;motion_dt=dt
        last_flow_cycle=.not.native_flow.or.count_steps==flow_cycles-1
        breakup_dt=kick_dt
        if(native_flow) breakup_dt=kick_fraction*interval
        if(motion=='kick'.or.options%fixed_parts) motion_dt=0
        ! EB flight follows release and moves parents and children together.
        ! Free-space parents have already moved; construct children from loc0.
        birth_drift=motion_dt
        if(present(eb_geometry)) birth_drift=0
        trial=drops;next_mass=masses;change=0;valid=.true.;trial_id=accepted_id;reynolds=0;trial_ledger=0;trial_force=0
        if(allocated(etab_tangent)) deallocate(etab_tangent)
        if(breakup_options%etab_parcel_factor>=0) then
          allocate(etab_tangent(np));etab_tangent=0
          if(native_flow) then
            do p=1,np
              tangent_index=tangent_slot(drops(p)%id)
              if(tangent_index>0) etab_tangent(p)=committed_tangent(tangent_index)
            enddo
          endif
        endif
        do p=1,np
          if (drops(p)%mass<=0) cycle
          if(options%mixture_time_scheme=='euler') then
            call advance_mixture_frozen_gas(species,transport,primitive(:,p),tg(p),trial(p),components, &
              next_mass(:,p),options,dt,change(:,p),valid,time_modifier=kick_fraction, &
              acceleration=accel,force_change=parcel_force,reynolds_number=reynolds(p))
          else
            call advance_mixture_frozen_gas(species,transport,primitive(:,p),tg(p),trial(p),components, &
              next_mass(:,p),options,kick_dt,change(:,p),valid,acceleration=accel,force_change=parcel_force, &
              reynolds_number=reynolds(p))
          end if
          if (.not.valid) exit
          trial_force=trial_force+parcel_force
          valid=all(ieee_is_finite(trial_force))
          if(.not.valid) exit
          if (trial(p)%mass>0.and.options%breakup_model/='none') then
            old_mass=trial(p)%mass;y=next_mass(:,p)/old_mass
            call liquid_mixture_properties(components,y,trial(p)%temperature,rho,cp,e,valid,primitive(5,p), &
              species(components%vapor_index)%molecular_weight)
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
            breakup_options=options
            call apply_mixture_breakup_viscosity(components,y,trial(p)%temperature,breakup_options,valid, &
              primitive(5,p),species(components%vapor_index)%molecular_weight)
            if(.not.valid) exit
            if (.not.breakup_options%khrt_children) then
              if(breakup_options%etab_parcel_factor>=0) then
                call advance_spray_breakup(trial(p),effective,breakup_options,primitive(1,p),primitive(2:4,p), &
                  film%density,film%viscosity,merge(dt,kick_dt,native_flow),valid,etab_delta,reynolds_number=reynolds(p))
                if(valid) etab_tangent(p)=etab_tangent(p)+etab_delta
                if(valid.and.motion=='kick'.and.last_flow_cycle) &
                  call advance_tab_no_birth_khrt(trial(p),effective,breakup_options,primitive(1,p),primitive(2:4,p), &
                    reynolds(p),breakup_dt,valid)
              else
                if(last_flow_cycle.or.breakup_options%breakup_model=='etab') &
                  call advance_spray_breakup(trial(p),effective,breakup_options,primitive(1,p), &
                  primitive(2:4,p),film%density,film%viscosity, &
                    merge(dt,breakup_dt,native_flow.and.breakup_options%breakup_model=='etab'), &
                    valid,reynolds_number=reynolds(p))
              endif
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
        if(agreed.and.breakup_options%etab_parcel_factor>=0.and.motion/='kick'.and.last_flow_cycle) then
          call split_mpi_reference_etab_parcels(comm,trial,next_mass,options%maximum_multiplicity, &
            min(125000,options%maximum_parcels),breakup_options%etab_parcel_factor,etab_tangent,birth_drift, &
            trial_id,etab_positions,etab_feedback,valid,birth_positions=positions)
          if(valid) then
            if(allocated(etab_change)) deallocate(etab_change)
            allocate(etab_change(nvar,np));etab_change=0;etab_change(2:5,:)=etab_feedback
          endif
          call agree()
        endif
        if (agreed.and.breakup_options%khrt_children.and.last_flow_cycle) then
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
            primitive(2:4,:),reynolds,breakup_dt,trial_id,child_energy,valid,eligible,gas_momentum=child_momentum, &
            pressure=primitive(5,:),gas_molecular_weights=species(components%vapor_index)%molecular_weight, &
            birth_positions=positions,child_drift_dt=birth_drift,make_new_drops=motion/='kick')
          if (valid) valid=size(trial)<=125000
          call agree()
          if (agreed) then
            change(2:4,:)=change(2:4,:)+child_momentum
            change(5,:)=change(5,:)+child_energy
          endif
        end if
        if(agreed.and..not.present(eb_geometry).and..not.present(boundary_faces)) then
          ! Newly released children can cross the periodic edge independently
          ! of their parent. Normalize their completed flight before publication.
          do p=1,size(trial)
            trial(p)%position=lower+modulo(trial(p)%position-lower,length)
          enddo
        endif
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
              outflow_faces,liquid_energies,impact_exterior(:5),frame=frame)
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
                eb_hits,valid,rejected,impact_exterior,boundary_faces,frame)
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
            gas_change=trial_change,maximum_level=maximum_level,containing_cell=containing_cell)
          if(valid.and.breakup_options%etab_parcel_factor>=0.and.motion/='kick') then
            if(present(gas_change)) then
              if(allocated(etab_trial_change)) deallocate(etab_trial_change)
              allocate(etab_trial_change)
              call deposit_mpi_amr_tree_spray_gas(comm,species,tree,owners,next_gas,next_temp, &
                etab_positions,etab_change,lower,spacing,periodic,valid,eb_geometry=eb_geometry, &
                nodal_phi=nodal_phi,gas_change=etab_trial_change,maximum_level=maximum_level, &
                containing_cell=containing_cell)
              if(valid) then
                do node=1,size(owners)
                  if(.not.allocated(trial_change%fields(node)%values)) cycle
                  trial_change%fields(node)%values=trial_change%fields(node)%values+ &
                    etab_trial_change%fields(node)%values
                enddo
              endif
            else
              call deposit_mpi_amr_tree_spray_gas(comm,species,tree,owners,next_gas,next_temp, &
                etab_positions,etab_change,lower,spacing,periodic,valid,eb_geometry=eb_geometry, &
                nodal_phi=nodal_phi,maximum_level=maximum_level,containing_cell=containing_cell)
            endif
          endif
          if (valid) then
            accepted=.true.;exit
          end if
        end if
        if(native_flow) return
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
      if(native_flow.and.breakup_options%etab_parcel_factor>=0) then
        do p=1,np
          tangent_index=tangent_slot(drops(p)%id)
          if(tangent_index>0) committed_tangent(tangent_index)=etab_tangent(p)
        enddo
      endif
      drops=trial;masses=next_mass;gas=next_gas;temp=next_temp;accepted_id=trial_id
      ledger=ledger+trial_ledger
      local_force=local_force+trial_force
      valid=all(ieee_is_finite(local_force));call agree();if(.not.agreed) return
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
        options%maximum_parcels,valid,kept_mass,maximum_level=maximum_level,periodic=periodic)
    else
      call migrate_amr_tree_spray_parcels(comm,tree,owners,survivors,lower,spacing,periodic, &
        options%maximum_parcels,valid,kept_mass,maximum_level=maximum_level)
    end if
    if (.not.valid) return
    if (present(wall_films)) then
      call migrate_mpi_amr_wall_films(comm,tree,owners,patches,patch_mass,lower,spacing, &
        film_settings%maximum_films,valid,nodal_phi,maximum_level,periodic=periodic)
      if (.not.valid) return
    end if
    call MPI_Allreduce(ledger,global_ledger,size(ledger),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    global_ledger(:5)=global_ledger(:5)-wall_ledger
    if (.not.all(ieee_is_finite(global_ledger))) return
    call MPI_Allreduce(local_force,global_force,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=all(ieee_is_finite(global_force));call agree();if(.not.agreed) return
    if(present(frame)) then
      call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,gas,lower,spacing,frame,.true.,valid,eb_geometry)
      if(.not.valid) return
      if(present(gas_change)) then
        call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,integrated_change,lower,spacing, &
          frame,.true.,valid,eb_geometry)
        if(.not.valid) return
      end if
    end if
    if(present(gas_change)) call move_alloc(integrated_change%fields,gas_change%fields)
    state=gas;temperature=temp
    call move_alloc(survivors,parcels);call move_alloc(kept_mass,component_mass)
    if (present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    if (present(last_id)) last_id=accepted_id
    if (present(boundary_change)) boundary_change=global_ledger
    if(present(eb_outflow_change)) eb_outflow_change=exterior
    if(present(force_change)) force_change=global_force
    substeps=count_steps;ok=.true.
  contains
    integer function tangent_slot(id) result(slot)
      integer(int64), intent(in) :: id
      integer :: lo,hi,mid
      slot=0;lo=1;hi=size(tangent_ids)
      do while(lo<=hi)
        mid=lo+(hi-lo)/2
        if(tangent_ids(mid)==id) then
          slot=mid;return
        else if(tangent_ids(mid)<id) then
          lo=mid+1
        else
          hi=mid-1
        endif
      enddo
    end function tangent_slot
    subroutine prepare_impacts()
      integer :: q,parent,lo,hi,mid,pass
      integer(int64) :: target_id
      if (allocated(impact_previous)) deallocate(impact_previous,impact_pressure)
      allocate(impact_previous(3,size(trial)),impact_pressure(size(trial)))
      valid=.true.
      do q=1,size(trial)
        parent=0;target_id=trial(q)%id
        do pass=1,2
          lo=1;hi=np;parent=0
          do while(lo<=hi)
            mid=lo+(hi-lo)/2
            if (drops(mid)%id==target_id) then
              parent=mid;exit
            else if (drops(mid)%id<target_id) then
              lo=mid+1
            else
              hi=mid-1
            end if
          end do
          if(parent>0) exit
          target_id=trial(q)%parent_id
        enddo
        if(parent==0) then
          valid=.false.;return
        endif
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
        if(present(frame).and.exchange(1)<0) then
          call rotating_potential_3d(frame,trial(q)%position,potential,valid)
          if(.not.valid) return
          exchange(5)=exchange(5)-exchange(1)*potential
        end if
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
