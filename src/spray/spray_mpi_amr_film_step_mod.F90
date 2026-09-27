module spray_mpi_amr_film_step_mod
  use mpi_f08
  use spray_coupling_3d_mod, only: valid_spray_boundary_faces
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use spray_parcel_mod, only: spray_liquid,spray_options,valid_spray_options,valid_spray_liquid
  use spray_mixture_mod, only: mixture_trace_fraction
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,valid_wall_options, &
    valid_wall_component_inventory,wall_liquid_conductivity_at_temperature,mixture_wall_film_rates, &
    advance_mixture_wall_film_exchange,advance_mixture_wall_film_frozen
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use spray_mpi_amr_film_mod, only: migrate_mpi_amr_wall_films,spray_level_film_population,refresh_mpi_eb_film_copies
  use spray_mpi_amr_field_mod, only: sample_mpi_amr_tree_spray_field,deposit_mpi_amr_tree_spray_gas
  use spray_mpi_amr_field_mod, only: sample_mpi_amr_particle_halo,deposit_mpi_amr_particle_halo, &
    sum_mpi_amr_particle_source_halo,deposit_mpi_amr_tree_spray_field
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  implicit none
  private
  public :: advance_mpi_amr_film_coupling,kick_mpi_eb_film_halo,advance_mpi_eb_film_sdc_stage
  public :: advance_mpi_eb_recursive_film_stage,finish_mpi_eb_film_level,append_mpi_eb_film_births
  type, public :: spray_sdc_film_stage
    integer :: next_stage=0,iterations=0,depth=-1
    real(dp) :: time=0,interval=0,old_wall(5)=0
    logical :: exports_physical=.false.,births_attached=.false.
    type(spray_level_film_population) :: midpoint(3)
    type(amr_tree_state_3d) :: old_source,old_remainder,old_physical_source
  end type spray_sdc_film_stage
contains
  ! Deposits created after the old particle exchange receive only the final
  ! film half kick. Keep existing midpoint films, rates and copies unchanged.
  subroutine append_mpi_eb_film_births(comm,tree,owners,lower,spacing,grow,maximum_films,phi,films,mass,work,ok)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),grow,maximum_films
    real(dp), intent(in) :: lower(3),spacing(3),phi(0:,0:,0:)
    type(spray_wall_film), intent(in) :: films(:)
    real(dp), intent(in) :: mass(:,:)
    type(spray_sdc_film_stage), intent(inout) :: work
    logical, intent(out) :: ok
    type(spray_level_film_population) :: candidate
    integer :: nc,nold,nnew,status
    logical :: valid,agreed
    ok=.false.;nc=size(mass,1);nnew=size(films)
    valid=work%next_stage==1.and.work%iterations>=1.and..not.work%births_attached
    valid=valid.and.allocated(work%midpoint(1)%films).and.allocated(work%midpoint(1)%component_mass)
    valid=valid.and.nc>=1.and.nc<=16.and.size(mass,2)==nnew.and.nnew<=1000000
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    nold=size(work%midpoint(1)%films)
    valid=all(shape(work%midpoint(1)%component_mass)==[nc,nold]).and.nold<=1000000-nnew
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    allocate(candidate%films(nold+nnew),candidate%component_mass(nc,nold+nnew))
    candidate%films(:nold)=work%midpoint(1)%films;candidate%films(nold+1:)=films
    candidate%component_mass(:,:nold)=work%midpoint(1)%component_mass
    candidate%component_mass(:,nold+1:)=mass
    call migrate_mpi_amr_wall_films(comm,tree,owners,candidate%films,candidate%component_mass,lower,spacing, &
      maximum_films,valid,phi,nodes=candidate%nodes,probes=candidate%probes,logical_level=work%depth,grow=grow)
    if(.not.valid) return
    work%midpoint(1)=candidate;work%births_attached=.true.;ok=.true.
  end subroutine append_mpi_eb_film_births

  ! Invoke after descendant gas synchronization. Stationary films need no drift
  ! redistribution; virtuals expire each interval, incoming ghosts at the last cycle.
  subroutine finish_mpi_eb_film_level(comm,depth,time,interval,iteration,ncycle,virtuals,ghosts,work,ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: depth,iteration,ncycle
    real(dp), intent(in) :: time,interval
    type(spray_level_film_population), intent(inout) :: virtuals(:),ghosts(:)
    type(spray_sdc_film_stage), intent(in) :: work(:)
    logical, intent(out) :: ok
    integer :: level,p,status,flags(4),lowflags(4),highflags(4)
    real(dp) :: times(2),lowtime(2),hightime(2),finish,tolerance
    logical :: valid,agreed
    ok=.false.;level=depth+1;times=[time,interval]
    valid=level>=1.and.level<=size(work).and.size(virtuals)==size(work).and.size(ghosts)==size(work)
    valid=valid.and.iteration>=1.and.iteration<=ncycle.and.all(ieee_is_finite(times))
    valid=valid.and.time>=0.and.interval>0
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    flags=[depth,iteration,ncycle,size(work)]
    call MPI_Allreduce(flags,lowflags,4,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,highflags,4,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(lowflags/=highflags)) return
    call MPI_Allreduce(times,lowtime,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(times,hightime,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(lowtime/=hightime)) return
    valid=time<=huge(1._dp)-interval
    valid=valid.and.work(level)%time==time.and.work(level)%interval==interval
    if(valid) then
      finish=time+interval;tolerance=32*epsilon(1._dp)*max(finish,interval,tiny(1._dp))
      do p=level,size(work)
        valid=valid.and.work(p)%next_stage==0.and.work(p)%iterations>0
        valid=valid.and.abs(work(p)%time+work(p)%interval-finish)<=tolerance
      end do
      valid=valid.and.allocated(virtuals(level)%films).and.allocated(virtuals(level)%component_mass).and. &
        allocated(virtuals(level)%nodes).and.allocated(virtuals(level)%probes)
      valid=valid.and.allocated(ghosts(level)%films).and.allocated(ghosts(level)%component_mass).and. &
        allocated(ghosts(level)%nodes).and.allocated(ghosts(level)%probes)
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call clear(virtuals(level))
    if(iteration==ncycle) call clear(ghosts(level))
    ok=.true.
  contains
    subroutine clear(group)
      type(spray_level_film_population), intent(inout) :: group
      group%films=group%films(:0);group%component_mass=group%component_mass(:,:0)
      group%nodes=group%nodes(:0);group%probes=group%probes(:,:0)
    end subroutine clear
  end subroutine finish_mpi_eb_film_level

  ! The full film hierarchy remains private to the outer gas transaction.
  subroutine advance_mpi_eb_recursive_film_stage(comm,species,transport,components,options,wall_options,tree,owners,geometry, &
      state,temperature,width,target_level,grow,maximum_films,physical,virtuals,ghosts,work, &
      time,interval,stage,iterations,faces, &
      source,remainder,wall_change,ok,nodal_phi,centroids,target_volume_fraction, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,physical_source)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_wall_options), intent(in) :: wall_options
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: owners(:),width,target_level,grow,maximum_films
    type(spray_level_film_population), intent(inout) :: physical(:),virtuals(:),ghosts(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature,centroids,target_volume_fraction
    real(dp), intent(in) :: time,interval,nodal_phi(0:,0:,0:)
    integer, intent(in) :: stage,iterations
    type(spray_sdc_film_stage), intent(inout) :: work(:)
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:),inflow_state(:,:),inflow_mach(6)
    real(dp), optional, intent(in) :: inflow_temperature(6),inflow_pressure(6)
    type(amr_tree_state_3d), optional, intent(inout) :: physical_source
    type(amr_tree_state_3d), allocatable :: physical_delta
    type(spray_sdc_film_stage) :: candidate
    type(spray_level_film_population) :: group(3),new_virtual,next_ghost
    type(amr_tree_state_3d) :: change,leftover,halo,halo_t,volume
    real(dp) :: wall(5),lower(3),spacing(3),tolerance,finish
    integer :: levels,level,p,status,flags(4),low(4),high(4)
    real(dp) :: times(2),low_time(2),high_time(2)
    logical :: valid,agreed,periodic(3)
    ok=.false.;valid=tree%is_valid().and.size(geometry)==size(owners)
    if(valid) valid=size(owners)==size(tree%nodes)
    call agree();if(.not.agreed) return
    levels=tree%level_count();level=target_level+1
    valid=all([size(physical),size(virtuals),size(ghosts),size(work)]==levels)
    valid=valid.and.level>=1.and.level<=levels.and.all(ieee_is_finite([time,interval]))
    valid=valid.and.time>=0.and.interval>0.and..not.any(faces=='periodic')
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
      call refresh_mpi_eb_film_copies(comm,tree,owners,physical,lower,spacing,target_level,grow, &
        maximum_films,nodal_phi,new_virtual,next_ghost,valid)
      if(.not.valid) return
      group(2)=new_virtual
    end if
    if(stage==0.or.stage==iterations) then
      call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,halo,halo_t,comm,valid, &
        ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
        inflow_mach,inflow_temperature,inflow_pressure,volume,width,target_volume_fraction, &
        faces,inflow_state,restrict_snapshot=.false.)
      if(.not.valid) return
    end if
    call advance_mpi_eb_film_sdc_stage(comm,species,transport,components,options,wall_options,tree,owners,geometry, &
      halo,halo_t,width,target_level,group,candidate,lower,spacing,periodic,time,interval,stage,iterations, &
      volume,nodal_phi,centroids,change,leftover,wall,valid,physical_delta)
    if(.not.valid) return
    physical(level)=group(1);virtuals(level)=group(2);ghosts(level)=group(3);work(level)=candidate
    if(stage==0.and.level<levels) ghosts(level+1)=next_ghost
    if(present(physical_source)) call move_alloc(physical_delta%fields,physical_source%fields)
    call move_alloc(change%fields,source%fields);call move_alloc(leftover%fields,remainder%fields)
    wall_change=wall;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
  end subroutine advance_mpi_eb_recursive_film_stage

  ! Populations are physical / virtual / ghost. Their halos are supplied at
  ! the current gas stage time by the recursive owner. Sources returned are rates.
  subroutine advance_mpi_eb_film_sdc_stage(comm,species,transport,components,options,wall_options,tree,owners,geometry, &
      gas_halo,temperature_halo,width,target_level,populations,work,lower,spacing,periodic,time,interval,stage,iterations, &
      fluid_fraction,nodal_phi,centroids,source,remainder,wall_change,ok,physical_source)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_wall_options), intent(in) :: wall_options
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: owners(:),width,target_level,stage,iterations
    type(amr_tree_state_3d), intent(in) :: gas_halo,temperature_halo,fluid_fraction,centroids
    real(dp), intent(in) :: lower(3),spacing(3),time,interval,nodal_phi(0:,0:,0:)
    logical, intent(in) :: periodic(3)
    type(spray_level_film_population), intent(inout) :: populations(3)
    type(spray_sdc_film_stage), intent(inout) :: work
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: physical_source
    type(spray_level_film_population) :: trial(3)
    type(spray_sdc_film_stage) :: candidate
    type(amr_tree_state_3d) :: total,leftover,part,part_remainder,exported
    type(amr_tree_state_3d), allocatable :: physical_delta
    real(dp) :: half_interval,wall(5),physical_wall(5),times(2),lowtime(2),hightime(2)
    integer :: flags(5),lowflags(5),highflags(5),population,status
    logical :: valid,agreed
    ok=.false.;times=[time,interval]
    valid=all(ieee_is_finite(times)).and.time>=0.and.interval>0.and.interval<=1e6_dp
    valid=valid.and.iterations>=1.and.stage>=0.and.stage<=iterations
    if(stage==0) then
      valid=valid.and.work%next_stage==0
    else
      valid=valid.and.work%next_stage==stage.and.work%iterations==iterations.and.work%depth==target_level
      valid=valid.and.work%time==time.and.work%interval==interval
      valid=valid.and.(work%exports_physical.eqv.present(physical_source))
    end if
    call agree();if(.not.agreed) return
    flags=[width,target_level,stage,iterations,merge(1,0,present(physical_source))]
    call MPI_Allreduce(flags,lowflags,5,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,highflags,5,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(lowflags/=highflags)) return
    call MPI_Allreduce(times,lowtime,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(times,hightime,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(lowtime/=hightime)) return
    half_interval=.5_dp*interval
    if(half_interval<=0) return
    if(stage>0.and.stage<iterations) then
      source=work%old_source;remainder=work%old_remainder;wall_change=work%old_wall
      if(present(physical_source)) physical_source=work%old_physical_source
      work%next_stage=stage+1;ok=.true.;return
    end if
    candidate=work
    if(stage==0) then
      trial=populations
    else
      trial=work%midpoint
    end if
    valid=.true.
    do population=1,3
      valid=valid.and.allocated(trial(population)%films).and.allocated(trial(population)%component_mass).and. &
        allocated(trial(population)%probes).and.allocated(trial(population)%nodes)
    end do
    call agree();if(.not.agreed) return
    if(present(physical_source)) allocate(physical_delta)
    do population=1,3
      call kick_mpi_eb_film_halo(comm,species,transport,components,options,wall_options,tree,owners,geometry, &
        gas_halo,temperature_halo,width,target_level,trial(population)%nodes,trial(population)%probes, &
        trial(population)%films,trial(population)%component_mass,lower,spacing,periodic,half_interval, &
        fluid_fraction,nodal_phi,centroids,part,part_remainder,wall,valid,physical_delta)
      if(.not.valid) return
      call remove_depleted(trial(population))
      if(population==1) then
        total=part;leftover=part_remainder;physical_wall=wall
        if(present(physical_source)) then
          exported=physical_delta;deallocate(physical_delta)
        end if
      else
        call add_field(total,part);if(.not.agreed) return
        call add_field(leftover,part_remainder);if(.not.agreed) return
      end if
    end do
    call to_rate(total);if(.not.agreed) return
    call to_rate(leftover);if(.not.agreed) return
    if(present(physical_source)) then
      call to_rate(exported);if(.not.agreed) return
    end if
    if(stage==0) then
      candidate%births_attached=.false.
      candidate%time=time;candidate%interval=interval;candidate%depth=target_level;candidate%iterations=iterations
      candidate%midpoint=trial;candidate%old_source=total;candidate%old_remainder=leftover;candidate%old_wall=physical_wall
      candidate%exports_physical=present(physical_source)
      if(present(physical_source)) candidate%old_physical_source=exported
    else
      valid=all(ieee_is_finite(candidate%old_wall))
      if(valid) valid=.not.any(physical_wall>huge(1._dp)-max(candidate%old_wall,0._dp)).and. &
        .not.any(physical_wall< -huge(1._dp)-min(candidate%old_wall,0._dp))
      call agree();if(.not.agreed) return
      physical_wall=physical_wall+candidate%old_wall
    end if
    candidate%next_stage=stage+1
    if(stage==iterations) then
      candidate%next_stage=0;populations=trial
    end if
    work=candidate;wall_change=physical_wall
    call move_alloc(total%fields,source%fields);call move_alloc(leftover%fields,remainder%fields)
    if(present(physical_source)) call move_alloc(exported%fields,physical_source%fields)
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
    subroutine remove_depleted(group)
      type(spray_level_film_population), intent(inout) :: group
      type(spray_level_film_population) :: kept
      logical :: live(size(group%films))
      integer :: c
      live=group%films%mass>0
      kept%films=pack(group%films,live);kept%nodes=pack(group%nodes,live)
      allocate(kept%component_mass(size(components),count(live)),kept%probes(3,count(live)))
      do c=1,size(components)
        kept%component_mass(c,:)=pack(group%component_mass(c,:),live)
      end do
      do c=1,3
        kept%probes(c,:)=pack(group%probes(c,:),live)
      end do
      group=kept
    end subroutine remove_depleted
    subroutine add_field(destination,addition)
      type(amr_tree_state_3d), intent(inout) :: destination
      type(amr_tree_state_3d), intent(in) :: addition
      integer :: node
      valid=.true.
      do node=1,size(owners)
        if(.not.allocated(destination%fields(node)%values)) cycle
        associate(a=>destination%fields(node)%values,b=>addition%fields(node)%values)
          valid=valid.and.all(ieee_is_finite(a)).and.all(ieee_is_finite(b))
          if(valid) valid=.not.any(a>huge(1._dp)-max(b,0._dp)).and..not.any(a< -huge(1._dp)-min(b,0._dp))
        end associate
      end do
      call agree();if(.not.agreed) return
      do node=1,size(owners)
        if(allocated(destination%fields(node)%values)) &
          destination%fields(node)%values=destination%fields(node)%values+addition%fields(node)%values
      end do
    end subroutine add_field
    subroutine to_rate(field)
      type(amr_tree_state_3d), intent(inout) :: field
      integer :: node
      valid=.true.
      do node=1,size(owners)
        if(.not.allocated(field%fields(node)%values)) cycle
        valid=valid.and.all(ieee_is_finite(field%fields(node)%values))
        if(half_interval<1) valid=valid.and.all(abs(field%fields(node)%values)<=huge(1._dp)*half_interval)
      end do
      call agree();if(.not.agreed) return
      do node=1,size(owners)
        if(allocated(field%fields(node)%values)) field%fields(node)%values=field%fields(node)%values/half_interval
      end do
    end subroutine to_rate
  end subroutine advance_mpi_eb_film_sdc_stage

  ! One frozen-gas film half step. Caller supplies owner-level nodes and probes
  ! displaced into the fluid, and handles copy selection and stage publication.
  ! Zero-mass outputs must be removed before a subsequent film call.
  subroutine kick_mpi_eb_film_halo(comm,species,transport,components,options,wall_options,tree,owners,geometry, &
      gas_halo,temperature_halo,width,target_level,nodes,probes,films,film_mass,lower,spacing,periodic,dt, &
      fluid_fraction,nodal_phi,centroids,source,remainder,wall_change,ok,physical_source)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_wall_options), intent(in) :: wall_options
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: owners(:),width,target_level,nodes(:)
    type(amr_tree_state_3d), intent(in) :: gas_halo,temperature_halo,fluid_fraction,centroids
    real(dp), intent(in) :: probes(:,:),lower(3),spacing(3),dt,nodal_phi(0:,0:,0:)
    logical, intent(in) :: periodic(3)
    type(spray_wall_film), intent(inout) :: films(:)
    real(dp), intent(inout) :: film_mass(:,:)
    type(amr_tree_state_3d), intent(inout) :: source,remainder
    real(dp), intent(inout) :: wall_change(5)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: physical_source
    type(spray_wall_film), allocatable :: candidate(:)
    real(dp), allocatable :: masses(:,:),sampled(:,:),sampled_t(:,:),change(:,:)
    type(amr_tree_state_3d) :: raw,routed,leftover,physical_delta
    real(dp) :: primitive(reactive_nvar(size(species))),tg,sound,conductivity,wall_energy
    real(dp) :: ledger,global_ledger,controls(3),low(3),high(3)
    integer :: flags(4),lowflags(4),highflags(4),nv,nc,p,c,node,n(3),rank,status
    logical :: valid,agreed
    character(len=9) :: wall_faces(6)
    ok=.false.;nv=reactive_nvar(size(species));nc=size(components);wall_faces='slip_wall'
    valid=tree%is_valid().and.valid_spray_options(options).and.valid_wall_options(wall_options,wall_faces)
    valid=valid.and.all(wall_options%temperature==wall_options%temperature(1)).and.ieee_is_finite(dt)
    valid=valid.and.dt>=0.and.dt<=1e6_dp.and.nc>=1.and.nc<=min(16,size(species))
    valid=valid.and.size(transport)==size(species).and.size(films)<=wall_options%maximum_films
    valid=valid.and.all(shape(film_mass)==[nc,size(films)]).and.all(shape(probes)==[3,size(films)])
    valid=valid.and.size(nodes)==size(films).and.target_level>=0
    call agree();if(.not.agreed) return
    valid=all(nodes>=1).and.all(nodes<=size(tree%nodes)).and.target_level<tree%level_count()
    call agree();if(.not.agreed) return
    do p=1,size(films)
      valid=valid.and.valid_wall_component_inventory(films(p),film_mass(:,p))
      valid=valid.and.tree%nodes(nodes(p))%level==target_level
    end do
    do c=1,nc
      valid=valid.and.valid_spray_liquid(components(c)).and.components(c)%vapor_index<=size(species)
    end do
    call agree();if(.not.agreed) return
    flags=[width,target_level,nc,merge(1,0,present(physical_source))]
    controls=[dt,wall_options%temperature(1),merge(1._dp,0._dp,wall_options%evolve)]
    call MPI_Allreduce(flags,lowflags,4,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,highflags,4,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(lowflags/=highflags)) return
    call MPI_Allreduce(controls,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    call sample_mpi_amr_particle_halo(comm,tree,owners,gas_halo,nv,width,nodes,probes,lower,spacing, &
      periodic,sampled,valid,fluid_fraction,nodal_phi,centroids)
    if(.not.valid) return
    call sample_mpi_amr_particle_halo(comm,tree,owners,temperature_halo,1,width,nodes,probes,lower,spacing, &
      periodic,sampled_t,valid,fluid_fraction,nodal_phi,centroids)
    if(.not.valid) return
    candidate=films;masses=film_mass;allocate(change(nv,size(films)));change=0;ledger=0
    if(wall_options%evolve) then
      do p=1,size(films)
        call reactive_conserved_to_primitive(species,sampled(:,p),sampled_t(1,p),primitive,tg,sound,valid)
        if(.not.valid) exit
        conductivity=wall_liquid_conductivity_at_temperature(wall_options,films(p)%material_id,films(p)%temperature)
        call advance_mixture_wall_film_frozen(candidate(p),components,masses(:,p),species,transport,primitive,tg, &
          options,conductivity,wall_options%temperature(1),dt,wall_options%minimum_height,change(:,p),wall_energy,valid)
        if(.not.valid) exit
        if(candidate(p)%mass>0) valid= &
          wall_liquid_conductivity_at_temperature(wall_options,candidate(p)%material_id,candidate(p)%temperature)>0
        if(.not.valid) exit
        valid=ieee_is_finite(wall_energy)
        if(valid) valid=ledger<=huge(1._dp)+min(wall_energy,0._dp).and. &
          ledger>=-huge(1._dp)+max(wall_energy,0._dp)
        if(.not.valid) exit
        ledger=ledger-wall_energy
      end do
    end if
    call agree();if(.not.agreed) return
    call deposit_mpi_amr_particle_halo(comm,tree,owners,gas_halo,nv,width,nodes,probes,change,lower,spacing, &
      periodic,raw,valid,fluid_fraction,nodal_phi,centroids,containing_cell=.true.)
    if(.not.valid) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    allocate(routed%fields(size(owners)))
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(routed%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
    end do
    call sum_mpi_amr_particle_source_halo(comm,tree,owners,raw,nv,width,target_level,periodic, &
      routed,leftover,valid,route_to_coarse=.false.,eb_geometry=geometry,halo_volume_fraction=fluid_fraction)
    if(.not.valid) return
    if(present(physical_source)) then
      physical_delta=routed
      do node=1,size(owners)
        if(owners(node)==rank) physical_delta%fields(node)%values=0
      end do
      call deposit_mpi_amr_tree_spray_field(comm,tree,owners,physical_delta,nv,probes,change,lower,spacing, &
        periodic,valid,eb_geometry=geometry,nodal_phi=nodal_phi,containing_cell=.true.)
      if(.not.valid) return
    end if
    call MPI_Allreduce(ledger,global_ledger,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if(status/=MPI_SUCCESS.or..not.ieee_is_finite(global_ledger)) return
    films=candidate;film_mass=masses;wall_change=0;wall_change(5)=global_ledger
    call move_alloc(routed%fields,source%fields);call move_alloc(leftover%fields,remainder%fields)
    if(present(physical_source)) call move_alloc(physical_delta%fields,physical_source%fields)
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
  end subroutine kick_mpi_eb_film_halo
  ! Common film substeps, owner-local gas feedback and whole-interval publication.
  ! All model tables/options are identical across ranks by caller contract.
  subroutine advance_mpi_amr_film_coupling(comm,species,transport,components,options,wall_options,state,temperature, &
      tree,owners,films,film_mass,lower,spacing,interval,faces,substeps,boundary_change,ok,eb_geometry,nodal_phi, &
      gas_change,maximum_level)
    integer, optional, intent(in) :: maximum_level
    type(amr_tree_state_3d), optional, intent(inout) :: gas_change
    type(amr_tree_state_3d) :: integrated_change
    type(amr_tree_state_3d), allocatable :: trial_change
    integer :: export_flag,export_low,export_high
    type(MPI_Comm), intent(in) :: comm
    type(eb_geometry_3d), optional, intent(in) :: eb_geometry(:)
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_wall_options), intent(in) :: wall_options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(amr_patch_tree_3d), intent(in) :: tree
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    real(dp), allocatable, intent(inout) :: film_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),interval
    integer, intent(in) :: owners(:)
    character(len=*), intent(in) :: faces(6)
    integer, intent(out) :: substeps
    real(dp), intent(inout) :: boundary_change(:)
    logical, intent(out) :: ok
    type(spray_wall_film), allocatable :: patches(:),trial(:),kept(:)
    type(amr_tree_state_3d) :: gas,temp,next_gas,next_temp,field
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    logical :: periodic(3)
    character(len=32) :: reference_faces(6)
    real(dp) :: upper(3),minimum_h(3)
    real(dp), allocatable :: masses(:,:),next_mass(:,:),kept_mass(:,:),positions(:,:),samples(:,:),change(:,:)
    real(dp), allocatable :: rates(:,:),hv(:,:),qg(:),qw(:)
    real(dp) :: primitive(reactive_nvar(size(species))),tg,sound,conductivity,tdot,elapsed,dt,local_dt,wall_energy
    real(dp) :: ledger,trial_ledger,global_ledger,controls(3),low(3),high(3)
    integer :: np,nc,nvar,p,c,axis,face,attempt,count_steps,ierr,node,shape3(3)
    logical :: valid,agreed,accepted,eb
    integer :: flag,lowflag,highflag
    real(dp) :: probe_offset,wall_temperature
    ok=.false.;substeps=0
    export_flag=merge(1,0,present(gas_change))
    call MPI_Allreduce(export_flag,export_low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(export_flag,export_high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.export_low/=export_high) return
    flag=merge(1,0,present(eb_geometry))+2*merge(1,0,present(nodal_phi))
    call MPI_Allreduce(flag,lowflag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(flag,highflag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.lowflag/=highflag) return
    if (flag/=0.and.flag/=3) return
    eb=flag==3
    valid=allocated(films).and.allocated(film_mass).and.valid_spray_options(options)
    reference_faces=faces
    if (eb) reference_faces='slip_wall'
    valid=valid.and.valid_wall_options(wall_options,reference_faces).and.ieee_is_finite(interval)
    if (eb) valid=valid.and.all(wall_options%temperature==wall_options%temperature(1)).and. &
      valid_spray_boundary_faces(faces)
    valid=valid.and.interval>=0.and.interval<=1e6_dp.and.size(transport)==size(species)
    call agree()
    if (.not.agreed) return
    np=size(films);nc=size(components);nvar=reactive_nvar(size(species))
    valid=nc>=1.and.nc<=min(16,size(species)).and.np<=min(125000,wall_options%maximum_films)
    valid=valid.and.all(shape(film_mass)==[nc,np]).and.size(boundary_change)==nvar
    call agree()
    if (.not.agreed) return
    controls=[interval,real(nc,dp),merge(1.0_dp,0.0_dp,wall_options%evolve)]
    call MPI_Allreduce(controls,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.any(low/=high)) return
    do p=1,np
      valid=valid.and.valid_wall_component_inventory(films(p),film_mass(:,p))
    end do
    do c=1,nc
      valid=valid.and.valid_spray_liquid(components(c)).and.components(c)%vapor_index<=size(species)
    end do
    call agree()
    if (.not.agreed) return
    reference_faces=faces
    call MPI_Bcast(reference_faces,32*6,MPI_CHARACTER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(faces==reference_faces).and.all(len_trim(faces)<=32)
    call agree();if (.not.agreed) return
    periodic=faces([1,3,5])=='periodic'
    gas=state;temp=temperature;patches=films;masses=film_mass
    ! Validate the tree and closed film positions before any stencil lookup.
    call migrate_mpi_amr_wall_films(comm,tree,owners,patches,masses,lower,spacing,wall_options%maximum_films,valid, &
      nodal_phi,maximum_level)
    if (.not.valid) return
    np=size(patches)
    valid=np<=125000
    call agree();if (.not.agreed) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree();if (.not.agreed) return
    upper=lower+spacing*real(tree%nodes(1)%shape,dp)
    node=size(ratios)
    if(present(maximum_level)) node=maximum_level+1
    minimum_h=spacing*real(domain(:,1),dp)/real(domain(:,node),dp)

    allocate(positions(3,np),change(nvar,np),rates(nc,np),hv(nc,np),qg(np),qw(np))
    if (eb) then
      ! Migration already checks surface position/normal and rank agreement.
      probe_offset=max(1e-8_dp*minval((upper-lower)/real(shape(nodal_phi)-1,dp)), &
        128*epsilon(1._dp)*maxval(abs([lower,upper])))
      do p=1,np
        positions(:,p)=patches(p)%position+probe_offset*patches(p)%normal
      end do
    else
    do p=1,np
      axis=maxloc(abs(patches(p)%normal),dim=1)
      face=2*axis-merge(1,0,patches(p)%normal(axis)>0)
      valid=wall_options%enabled(face).and.abs(patches(p)%normal(axis))==1
      if (mod(face,2)==1) then
        valid=valid.and.patches(p)%position(axis)==lower(axis)
      else
        valid=valid.and.patches(p)%position(axis)==upper(axis)
      end if
      if (.not.valid) exit
      ! Probe adjacent fluid cells without moving the exact stationary wall patch.
      positions(:,p)=max(lower+.5_dp*minimum_h,min(upper-.5_dp*minimum_h,patches(p)%position))
    end do
    end if
    call agree();if (.not.agreed) return
    change=0
    call deposit_mpi_amr_tree_spray_gas(comm,species,tree,owners,gas,temp,positions,change, &
      lower,spacing,periodic,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi,maximum_level=maximum_level)
    if (.not.valid) return
    if(present(gas_change)) then
      integrated_change=gas
      do node=1,size(owners)
        if(allocated(integrated_change%fields(node)%values)) integrated_change%fields(node)%values=0
      end do
      allocate(trial_change)
    end if
    elapsed=0;ledger=0;count_steps=0
    do while(wall_options%evolve.and.elapsed<interval)
      valid=count_steps<options%maximum_substeps
      call agree()
      if (.not.agreed) return
      if (allocated(field%fields)) deallocate(field%fields)
      allocate(field%fields(size(tree%nodes)))
      do node=1,size(tree%nodes)
        if (.not.allocated(gas%fields(node)%values)) cycle
        shape3=tree%nodes(node)%shape
        allocate(field%fields(node)%values(nvar+1,shape3(1),shape3(2),shape3(3)))
        field%fields(node)%values(:nvar,:,:,:)=gas%fields(node)%values
        field%fields(node)%values(nvar+1,:,:,:)=temp%fields(node)%values(1,:,:,:)
      end do
      call sample_mpi_amr_tree_spray_field(comm,tree,owners,field,nvar+1,positions, &
        lower,spacing,periodic,samples,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi,maximum_level=maximum_level)
      if (.not.valid) return
      local_dt=interval-elapsed;rates=0;hv=0;qg=0;qw=0
      do p=1,np
        if (patches(p)%mass<=0) cycle
        if (eb) then
          wall_temperature=wall_options%temperature(1)
        else
          axis=maxloc(abs(patches(p)%normal),dim=1)
          face=2*axis-merge(1,0,patches(p)%normal(axis)>0)
          valid=wall_options%enabled(face)
          if (.not.valid) exit
          wall_temperature=wall_options%temperature(face)
        end if
        call reactive_conserved_to_primitive(species,samples(:nvar,p),samples(nvar+1,p),primitive,tg,sound,valid)
        if (.not.valid) exit
        conductivity=wall_liquid_conductivity_at_temperature(wall_options,patches(p)%material_id,patches(p)%temperature)
        call mixture_wall_film_rates(patches(p),components,masses(:,p),species,transport,primitive,tg,options, &
          conductivity,wall_temperature,rates(:,p),hv(:,p),qg(p),qw(p),tdot,valid)
        if (.not.valid) exit
        do c=1,nc
          if (rates(c,p)>0.and.masses(c,p)>mixture_trace_fraction*patches(p)%mass) &
            local_dt=min(local_dt,options%max_fractional_change*masses(c,p)/rates(c,p))
        end do
        if (abs(tdot)>0) local_dt=min(local_dt,min(2.0_dp,options%max_fractional_change*patches(p)%temperature)/abs(tdot))
      end do
      call agree()
      if (.not.agreed) return
      call MPI_Allreduce(local_dt,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      accepted=.false.
      do attempt=1,30
        if (dt<=0.or.elapsed+dt<=elapsed) exit
        trial=patches;next_mass=masses;change=0;trial_ledger=0;valid=.true.
        do p=1,np
          if (patches(p)%mass<=0) cycle
          call advance_mixture_wall_film_exchange(trial(p),components,next_mass(:,p),rates(:,p)*dt,hv(:,p), &
            qg(p)*dt,qw(p)*dt,wall_options%minimum_height,change(:,p),wall_energy,valid,deplete_traces=.true.)
          if (valid.and.trial(p)%mass>0) valid= &
            wall_liquid_conductivity_at_temperature(wall_options,trial(p)%material_id,trial(p)%temperature)>0
          if (.not.valid) exit
          trial_ledger=trial_ledger-wall_energy
        end do
        call agree()
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
      gas=next_gas;temp=next_temp;patches=trial;masses=next_mass;ledger=ledger+trial_ledger
      elapsed=elapsed+dt;count_steps=count_steps+1
    end do
    kept=pack(patches,patches%mass>0)
    allocate(kept_mass(nc,size(kept)))
    do c=1,nc
      kept_mass(c,:)=pack(masses(c,:),patches%mass>0)
    end do
    call migrate_mpi_amr_wall_films(comm,tree,owners,kept,kept_mass,lower,spacing,wall_options%maximum_films,valid, &
      nodal_phi,maximum_level)
    if (.not.valid) return
    call MPI_Allreduce(ledger,global_ledger,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (.not.ieee_is_finite(global_ledger)) return
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(kept,films);call move_alloc(kept_mass,film_mass)
    if(present(gas_change)) call move_alloc(integrated_change%fields,gas_change%fields)
    boundary_change=0;boundary_change(5)=global_ledger;substeps=count_steps;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
end module
