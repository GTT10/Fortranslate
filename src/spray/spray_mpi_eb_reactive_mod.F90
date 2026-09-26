module spray_mpi_eb_reactive_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options,valid_spray_options,valid_spray_parcel
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_mpi_amr_step_mod, only: advance_mpi_amr_mixture_coupling, &
    advance_mpi_eb_recursive_population_stage,spray_sdc_population_stage
  use spray_mpi_amr_tree_mod, only: spray_level_population,initialize_mpi_amr_sdc_populations, &
    collect_mpi_amr_level_particles,redistribute_mpi_amr_subtree_particles
  use mpi_eb_tree_halo_3d_mod, only: fill_sdc_mpi_eb_particle_predictor_halo_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_full_3d,advance_mpi_eb_tree_sdc_3d, &
    advance_mpi_eb_tree_sdc_subcycle_3d
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options
  use spray_mpi_amr_film_step_mod, only: advance_mpi_amr_film_coupling
  use spray_mpi_mod, only: split_mpi_spray_parcels
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d, &
    integrate_mpi_eb_tree_3d,average_down_mpi_eb_tree_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  implicit none
  private
  public :: advance_mpi_eb_spray_reactive,advance_mpi_eb_spray_sdc,advance_mpi_eb_spray_sdc_subcycle
contains
  ! Recursive frozen-gas drag/heat/evaporation core. Breakup, films and particle
  ! forcing require their own recursive lifecycle and are not enabled here.
  subroutine advance_mpi_eb_spray_sdc_subcycle(comm,tree,owners,geometry,extended_geometry,phi, &
      species,reactions,transport,liquids,options,state,temperature,parcels,component_mass,lower,spacing, &
      dt,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces,boundary_change,wall_change, &
      chemical_change,source_change,minimum_theta,particle_steps,last_id,ok,reaction_source,difmag, &
      redistribution,target_volume_fraction,gas_particle_change,level_steps,gas_acceleration, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),iterations
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),spacing(3),dt,force(3),rtol,atol
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    logical, intent(in) :: chemistry,transport_flags(4)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),minimum_theta
    integer, intent(out) :: particle_steps
    integer(int64), intent(inout) :: last_id
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), optional, intent(in) :: difmag,target_volume_fraction,gas_acceleration(3)
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(inout) :: gas_particle_change(:)
    integer(int64), optional, intent(inout) :: level_steps(:)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: wall_temperature(6)
    type(spray_level_population), allocatable :: physical(:),virtuals(:),ghosts(:)
    type(spray_sdc_population_stage), allocatable :: work(:)
    type(amr_tree_state_3d) :: gas,temp,metrics,centroids,particle_first,particle_last,source_correction
    type(amr_tree_state_3d), allocatable :: ir
    type(spray_parcel), allocatable :: final(:)
    real(dp), allocatable :: masses(:,:)
    integer(int64), allocatable :: counts(:)
    integer(int64) :: watermark,low_id,high_id,maximum_id,global_id,total_steps
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: external(size(source_change)),exchange(size(source_change)),particle_wall(size(source_change))
    real(dp) :: exterior(size(source_change)),theta,h(3),origin(3)
    integer :: nv,levels,rank,status,node,n(3),i,j,k,p,ios
    logical :: valid,agreed
    ok=.false.;particle_steps=0;nv=reactive_nvar(size(species))
    valid=allocated(parcels).and.allocated(component_mass).and.tree%is_valid().and.valid_spray_options(options)
    valid=valid.and.options%breakup_model=='none'.and.all(boundary_faces=='outflow'.or. &
      boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or. &
      boundary_faces=='inflow'.or.boundary_faces=='isothermal_wall')
    valid=valid.and.all([size(boundary_change),size(wall_change),size(chemical_change),size(source_change)]==nv)
    if(present(gas_particle_change)) valid=valid.and.size(gas_particle_change)==nv
    call agree();if(.not.agreed) return
    levels=tree%level_count()
    if(present(level_steps)) valid=size(level_steps)==levels
    call agree();if(.not.agreed) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(last_id,low_id,1,MPI_INTEGER8,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(last_id,high_id,1,MPI_INTEGER8,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low_id/=high_id) return
    maximum_id=max(0_int64,maxval(parcels%id))
    call MPI_Allreduce(maximum_id,global_id,1,MPI_INTEGER8,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.last_id<global_id) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,agreed)
    if(.not.agreed) return
    valid=all(lower==[geometry(1)%x_lower,geometry(1)%y_lower,geometry(1)%z_lower]).and. &
      all(spacing==[geometry(1)%dx,geometry(1)%dy,geometry(1)%dz])
    call agree();if(.not.agreed) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,agreed)
    if(.not.agreed) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,agreed)
    if(.not.agreed) return
    valid=size(extended_geometry)==size(owners)
    call agree();if(.not.agreed) return
    do node=1,size(owners)
      valid=extended_geometry(node)%is_valid()
      if(valid) valid=all([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]== &
        tree%nodes(node)%shape+12)
      if(.not.valid) exit
    end do
    call agree();if(.not.agreed) return
    allocate(work(levels),counts(levels),metrics%fields(size(owners)),centroids%fields(size(owners)),stat=ios)
    valid=ios==0;call agree();if(.not.agreed) return
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape;h=[geometry(node)%dx,geometry(node)%dy,geometry(node)%dz]
      origin=[geometry(node)%x_lower,geometry(node)%y_lower,geometry(node)%z_lower]
      allocate(metrics%fields(node)%values(1,-2:n(1)+3,-2:n(2)+3,-2:n(3)+3), &
        centroids%fields(node)%values(3,-2:n(1)+3,-2:n(2)+3,-2:n(3)+3),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      metrics%fields(node)%values(1,:,:,:)=extended_geometry(node)%volume_fraction(4:n(1)+9,4:n(2)+9,4:n(3)+9)
      do k=-2,n(3)+3
      do j=-2,n(2)+3
      do i=-2,n(1)+3
        centroids%fields(node)%values(:,i,j,k)=origin+(real([i,j,k],dp)-.5_dp+ &
          [extended_geometry(node)%cell_centroid_x(i+6,j+6,k+6), &
           extended_geometry(node)%cell_centroid_y(i+6,j+6,k+6), &
           extended_geometry(node)%cell_centroid_z(i+6,j+6,k+6)])*h
      end do
      end do
      end do
    end do
    call agree();if(.not.agreed) return
    call initialize_mpi_amr_sdc_populations(comm,tree,owners,parcels,component_mass,lower,spacing, &
      [.false.,.false.,.false.],2,options%maximum_parcels,physical,virtuals,ghosts,agreed,phi)
    if(.not.agreed) return
    gas=state;temp=temperature;particle_first=state
    do node=1,size(owners)
      if(owners(node)==rank) particle_first%fields(node)%values=0
    end do
    source_correction=particle_first
    particle_last=particle_first;particle_wall=0;exterior=0;watermark=last_id
    if(present(reaction_source)) then
      allocate(ir);ir=reaction_source
    end if
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,geometry,gas,temp, &
      dt,iterations,chemistry,transport_flags,force,rtol,atol,extended_geometry,boundary_faces, &
      boundary,wall,chemical,external,theta,comm,agreed,ir,difmag,level_steps=counts,redistribution=redistribution, &
      target_volume_fraction=target_volume_fraction,halo_source=particle_source,auxiliary_change=exchange, &
      gas_acceleration=gas_acceleration,predictor_builder=particle_predictor,inflow_state=inflow_state, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      wall_temperature=wall_temperature,post_level=particle_post_level)
    if(.not.agreed) return
    call synchronize_particle_sources()
    if(.not.agreed) return
    total_steps=0
    do p=1,levels
      valid=work(p)%next_stage==0.and.abs(work(p)%time+work(p)%interval-dt)<=32*epsilon(dt)*dt
      valid=valid.and.counts(p)<=int(huge(particle_steps)/2,int64)-total_steps
      call agree();if(.not.agreed) return
      total_steps=total_steps+counts(p)
    end do
    call collect_mpi_amr_level_particles(comm,tree,owners,physical,lower,spacing,[.false.,.false.,.false.], &
      options%maximum_parcels,final,masses,agreed,phi)
    if(.not.agreed) return
    valid=safe_sum(boundary,exterior).and.safe_sum(wall,particle_wall)
    call agree();if(.not.agreed) return
    boundary=boundary+exterior;wall=wall+particle_wall
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(final,parcels);call move_alloc(masses,component_mass)
    if(present(reaction_source)) reaction_source=ir
    if(present(gas_particle_change)) gas_particle_change=exchange
    if(present(level_steps)) level_steps=counts
    boundary_change=boundary;wall_change=wall;chemical_change=chemical;source_change=external
    minimum_theta=theta;last_id=watermark;particle_steps=int(2*total_steps);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
    logical function safe_sum(a,b)
      real(dp), intent(in) :: a(:),b(:)
      safe_sum=all(ieee_is_finite(a)).and.all(ieee_is_finite(b))
      if(safe_sum) safe_sum=.not.any(a>huge(1._dp)-max(b,0._dp)).and..not.any(a< -huge(1._dp)-min(b,0._dp))
    end function safe_sum
    subroutine add_source_rate(destination,rate,factor,good)
      real(dp), intent(inout) :: destination(:,:,:,:)
      real(dp), intent(in) :: rate(:,:,:,:),factor
      logical, intent(out) :: good
      real(dp) :: term(size(rate,1))
      integer :: x,y,z
      good=ieee_is_finite(factor).and.all(ieee_is_finite(rate))
      if(.not.good) return
      if(abs(factor)>1) then
        good=all(abs(rate)<=huge(1._dp)/abs(factor))
        if(.not.good) return
      end if
      do z=1,size(rate,4)
      do y=1,size(rate,3)
      do x=1,size(rate,2)
        term=factor*rate(:,x,y,z)
        good=safe_sum(destination(:,x,y,z),term)
        if(.not.good) return
        destination(:,x,y,z)=destination(:,x,y,z)+term
      end do
      end do
      end do
    end subroutine add_source_rate
    ! Conservative source synchronization: replace provisional copy-supported
    ! leaf exchange by the actual physical half-kick exchange at its receiver.
    ! This is a conservation extension to the pinned same-level source transfer.
    subroutine synchronize_particle_sources()
      real(dp) :: correction(nv),primitive(reactive_nprim(size(species))),recovered,sound
      integer :: p,x,y,z,d(3)
      call integrate_mpi_eb_tree_3d(tree,owners,geometry,source_correction,correction,comm,agreed)
      if(.not.agreed) return
      valid=safe_sum(exchange,correction)
      call agree();if(.not.agreed) return
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(covered(tree,p,[x,y,z]).or.geometry(p)%volume_fraction(x,y,z)<=0) cycle
          if(all(source_correction%fields(p)%values(:,x,y,z)==0)) cycle
          valid=safe_sum(gas%fields(p)%values(:,x,y,z),source_correction%fields(p)%values(:,x,y,z))
          if(.not.valid) exit
          gas%fields(p)%values(:,x,y,z)=gas%fields(p)%values(:,x,y,z)+source_correction%fields(p)%values(:,x,y,z)
          call reactive_conserved_to_primitive(species,gas%fields(p)%values(:,x,y,z),temp%fields(p)%values(1,x,y,z), &
            primitive,recovered,sound,valid)
          if(.not.valid) exit
          temp%fields(p)%values(1,x,y,z)=recovered
        end do
        if(.not.valid) exit
        end do
        if(.not.valid) exit
        end do
        if(.not.valid) exit
      end do
      call agree();if(.not.agreed) return
      call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,gas,temp,comm,agreed)
      if(.not.agreed) return
      exchange=exchange+correction
    end subroutine synchronize_particle_sources
    subroutine particle_post_level(depth,time,interval,iteration,ncycle,post_comm,good)
      integer, intent(in) :: depth,iteration,ncycle
      real(dp), intent(in) :: time,interval
      type(MPI_Comm), intent(in) :: post_comm
      logical, intent(out) :: good
      logical :: synchronized
      integer :: p,status
      synchronized=.true.
      do p=depth+1,levels
        synchronized=synchronized.and.work(p)%next_stage==0.and. &
          abs(work(p)%time+work(p)%interval-time-interval)<=32*epsilon(dt)*max(dt,abs(time+interval))
      end do
      call MPI_Allreduce(synchronized,good,1,MPI_LOGICAL,MPI_LAND,post_comm,status)
      if(status/=MPI_SUCCESS) good=.false.
      if(.not.good) return
      if(depth==0.or.(iteration<ncycle.and.depth<levels-1)) then
        call redistribute_mpi_amr_subtree_particles(post_comm,tree,owners,physical,lower,spacing, &
          [.false.,.false.,.false.],depth,iteration,options%maximum_parcels,good,phi)
        if(.not.good) return
      end if
      virtuals(depth+1)%parcels=virtuals(depth+1)%parcels(:0)
      virtuals(depth+1)%component_mass=virtuals(depth+1)%component_mass(:,:0)
      virtuals(depth+1)%nodes=virtuals(depth+1)%nodes(:0)
      if(iteration==ncycle) then
        ghosts(depth+1)%parcels=ghosts(depth+1)%parcels(:0)
        ghosts(depth+1)%component_mass=ghosts(depth+1)%component_mass(:,:0)
        ghosts(depth+1)%nodes=ghosts(depth+1)%nodes(:0)
      end if
    end subroutine particle_post_level
    subroutine particle_source(stage_gas,stage_temp,interval,stage,depth,time,rate,source_comm,good, &
        ancestor_start,ancestor_end,temperature_start,temperature_end,fractions)
      type(amr_tree_state_3d), intent(in) :: stage_gas,stage_temp,ancestor_start,ancestor_end,temperature_start,temperature_end
      real(dp), intent(in) :: interval,time,fractions(:)
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: source_comm
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: remainder,physical_rate
      real(dp) :: wall_delta(5),outside(nv),wall_full(nv)
      integer :: p
      call advance_mpi_eb_recursive_population_stage(source_comm,species,transport,liquids,options,tree,owners,geometry, &
        stage_gas,stage_temp,3,depth,2,options%maximum_parcels,physical,virtuals,ghosts,work,time,interval,stage,iterations, &
        boundary_faces,rate,remainder,wall_delta,outside,good,phi,centroids,metrics, &
        ancestor_start,ancestor_end,temperature_start,temperature_end,fractions, &
        inflow_state,inflow_mach,inflow_temperature,inflow_pressure,physical_source=physical_rate)
      if(.not.good) return
      do p=1,size(owners)
        if(owners(p)/=rank.or.tree%nodes(p)%level/=depth) cycle
        if(stage==0) particle_first%fields(p)%values=rate%fields(p)%values
        particle_last%fields(p)%values=rate%fields(p)%values
      end do
      if(stage/=iterations) return
      ! All fields are rates. Count each accepted physical half kick once;
      ! remove only this logical level's provisional valid-cell source.
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        associate(delta=>source_correction%fields(p)%values,actual=>physical_rate%fields(p)%values, &
          old_actual=>work(depth+1)%old_physical_source%fields(p)%values)
          call add_source_rate(delta,old_actual,.5_dp*interval,valid)
          if(valid) call add_source_rate(delta,actual,.5_dp*interval,valid)
          if(valid.and.tree%nodes(p)%level==depth) then
            call add_source_rate(delta,work(depth+1)%old_source%fields(p)%values,-.5_dp*interval,valid)
            if(valid) call add_source_rate(delta,rate%fields(p)%values,-.5_dp*interval,valid)
          end if
        end associate
        if(.not.valid) exit
      end do
      call agree();good=agreed
      if(.not.good) return
      wall_full=0;wall_full(:5)=wall_delta
      valid=safe_sum(particle_wall,wall_full).and.safe_sum(exterior,outside)
      call agree();good=agreed
      if(.not.good) return
      particle_wall=particle_wall+wall_full;exterior=exterior+outside
    end subroutine particle_source
    subroutine particle_predictor(rate,interval,iteration,depth,time,halo,source_comm,good, &
        source_start,source_end,fractions)
      type(amr_tree_state_3d), intent(in) :: rate,source_start,source_end
      real(dp), intent(in) :: interval,time,fractions(:)
      integer, intent(in) :: iteration,depth
      type(amr_tree_state_3d), intent(inout) :: halo
      type(MPI_Comm), intent(in) :: source_comm
      logical, intent(out) :: good
      valid=work(depth+1)%time==time.and.work(depth+1)%interval==interval.and.work(depth+1)%next_stage==iteration
      call agree();good=agreed
      if(.not.good) return
      call fill_sdc_mpi_eb_particle_predictor_halo_3d(tree,owners,geometry,rate,work(depth+1)%old_source, &
        work(depth+1)%old_remainder,nv,3,halo,source_comm,good,boundary_faces,metrics, &
        source_start,source_end,particle_first,particle_last,fractions)
    end subroutine particle_predictor
  end subroutine advance_mpi_eb_spray_sdc_subcycle

  ! S/2 - F/2 - EB reaction/transport/hydro - F/2 - S/2; atomic publication.
  ! Supplying last_id with film state enables EB impacts in both particle halves.
  ! Replicated tables/options and a consistent geometry/nodal field are caller contracts.
  ! Ledgers use change-in-system signs, including the particle wall impulse.
  subroutine advance_mpi_eb_spray_reactive(comm,tree,owners,geometry,phi,species,reactions,transport, &
      liquids,options,state,temperature,parcels,component_mass,lower,spacing,dt,solver, &
      chemistry,transport_enabled,viscosity,conduction,diffusion,barodiffusion,rtol,atol, &
      particle_steps,minimum_theta,boundary_change,wall_change,chemical_change,ok,integrator,subcycling, &
      films,film_mass,film_settings,film_steps,last_id,model_rejected,split_created,particle_outflow, &
      acceleration,source_change,inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry, &
      boundary_faces,inflow_state,wall_temperature)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),spacing(3),dt,rtol,atol
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    character(len=*), intent(in) :: solver
    logical, intent(in) :: chemistry,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    integer, intent(out) :: particle_steps
    real(dp), intent(inout) :: minimum_theta,boundary_change(:),wall_change(:),chemical_change(:)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    real(dp), optional, intent(in) :: difmag
    real(dp) :: coefficient,low_coefficient,high_coefficient
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),wall_temperature(6)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    logical, optional, intent(in) :: subcycling
    type(spray_wall_film), allocatable, optional, intent(inout) :: films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    integer, optional, intent(out) :: film_steps
    integer(int64), optional, intent(inout) :: last_id
    logical, optional, intent(out) :: model_rejected
    integer(int64), optional, intent(out) :: split_created
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(inout) :: source_change(:)
    real(dp) :: accel(3),low_accel(3),high_accel(3),force_ledger(reactive_nvar(size(species)))
    logical, optional, intent(in) :: particle_outflow
    character(len=32), allocatable :: particle_faces(:)
    logical :: outflow
    real(dp), allocatable :: half_exterior(:)
    real(dp) :: exterior(reactive_nvar(size(species)))
    integer(int64) :: before_split
    real(dp) :: low_multiplicity,high_multiplicity
    integer(int64) :: watermark
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    real(dp) :: first_film(reactive_nvar(size(species))),second_film(reactive_nvar(size(species)))
    integer :: flags,lowflags,highflags,first_film_steps,second_film_steps
    character(len=32) :: film_faces(6)
    type(amr_tree_state_3d) :: gas,temp
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    real(dp) :: chemical(reactive_nvar(size(species))),first(reactive_nvar(size(species)))
    real(dp) :: second(reactive_nvar(size(species))),theta
    integer :: first_steps,second_steps,ierr,nv
    logical :: valid,agreed
    ok=.false.;particle_steps=0;nv=size(boundary)
    if(present(film_steps)) film_steps=0
    if(present(split_created)) split_created=0
    if(present(model_rejected)) model_rejected=.false.
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    valid=ieee_is_finite(coefficient).and.coefficient>=0
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(coefficient,low_coefficient,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(coefficient,high_coefficient,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.low_coefficient/=high_coefficient) return
    accel=0;force_ledger=0
    if(present(acceleration)) accel=acceleration
    valid=all(ieee_is_finite(accel))
    valid=valid.and.all(abs(accel)<=1e6_dp)
    if(any(accel/=0)) valid=valid.and.present(source_change)
    if(present(source_change)) valid=valid.and.size(source_change)==nv
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(accel,low_accel,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(accel,high_accel,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(low_accel/=high_accel)) return
    outflow=.false.;exterior=0
    if(present(particle_outflow)) outflow=particle_outflow
    valid=.true.
    if(present(boundary_faces)) then
      valid=all(boundary_faces=='inflow'.or.boundary_faces=='outflow'.or.boundary_faces=='periodic'.or. &
        boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or.boundary_faces=='isothermal_wall')
      if(any(boundary_faces=='periodic'.or.boundary_faces=='slip_wall'.or. &
        boundary_faces=='no_slip_wall'.or.boundary_faces=='isothermal_wall')) then
        particle_faces=boundary_faces
        where(particle_faces=='inflow') particle_faces='outflow'
        where(particle_faces=='isothermal_wall') particle_faces='no_slip_wall'
      end if
      valid=valid.and.outflow
    else
      valid=.not.present(inflow_state)
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    watermark=0
    if(present(last_id)) watermark=last_id
    flags=merge(1,0,present(films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))+ &
      8*merge(1,0,present(last_id))+16*merge(1,0,outflow)
    call MPI_Allreduce(flags,lowflags,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,highflags,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.lowflags/=highflags) return
    if(outflow) then
      if(flags/=16.and.flags/=24.and.flags/=31) return
      allocate(half_exterior(nv));half_exterior=0
      flags=flags-16
    end if
    if(flags/=0.and.flags/=7.and.flags/=8.and.flags/=15) return
    valid=valid_spray_options(options)
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    valid=valid.and.(options%maximum_multiplicity==0.or.present(last_id))
    if(present(films)) valid=valid.and.allocated(films).and.allocated(film_mass)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(options%maximum_multiplicity,low_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(options%maximum_multiplicity,high_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.low_multiplicity/=high_multiplicity) return
    if(present(films)) then
      patches=films;patch_mass=film_mass
    end if
    first_film=0;second_film=0;first_film_steps=0;second_film_steps=0;film_faces='slip_wall'
    if(present(boundary_faces)) then
      where(boundary_faces=='periodic') film_faces='periodic'
    end if
    valid=allocated(parcels).and.allocated(component_mass).and.ieee_is_finite(dt)
    if(valid) valid=dt>0.and.dt<=1e6_dp
    valid=valid.and.size(boundary_change)==nv.and.size(wall_change)==nv.and.size(chemical_change)==nv
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    boundary=0;wall=0;chemical=0;first=0;second=0;theta=1
    if(any(accel/=0)) then
      call force_kick()
      if(.not.valid) return
    end if
    call particle_half(first_steps,first)
    if(.not.valid) return
    if(present(films)) then
      call advance_mpi_amr_film_coupling(comm,species,transport,liquids,options,film_settings,gas,temp, &
        tree,owners,patches,patch_mass,lower,spacing,.5_dp*dt,film_faces,first_film_steps,first_film,valid,geometry,phi)
      if(.not.valid) return
    end if
    call advance_mpi_eb_tree_full_3d(species,reactions,transport,tree,owners,geometry,solver,dt, &
      chemistry,rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion, &
      gas,temp,boundary,wall,chemical,theta,comm,valid,integrator,subcycling, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      difmag=coefficient,extended_geometry=extended_geometry,boundary_faces=boundary_faces, &
      inflow_state=inflow_state,wall_temperature=wall_temperature)
    if(.not.valid) return
    if(present(films)) then
      call advance_mpi_amr_film_coupling(comm,species,transport,liquids,options,film_settings,gas,temp, &
        tree,owners,patches,patch_mass,lower,spacing,.5_dp*dt,film_faces,second_film_steps,second_film,valid,geometry,phi)
      if(.not.valid) return
    end if
    call particle_half(second_steps,second)
    if(.not.valid) return
    before_split=watermark
    if(options%maximum_multiplicity>0) then
      call split_mpi_spray_parcels(comm,drops,options%maximum_multiplicity,options%maximum_parcels,watermark,valid,masses)
      if(.not.valid) return
    end if
    if(any(accel/=0)) then
      call force_kick()
      if(.not.valid) return
    end if
    state=gas;temperature=temp
    call move_alloc(drops,parcels)
    call move_alloc(masses,component_mass)
    if(present(films)) then
      call move_alloc(patches,films);call move_alloc(patch_mass,film_mass)
    end if
    if(present(last_id)) last_id=watermark
    if(present(split_created)) split_created=watermark-before_split
    if(present(film_steps)) film_steps=first_film_steps+second_film_steps
    boundary_change=boundary+exterior;wall_change=wall+first+second+first_film+second_film;chemical_change=chemical
    if(present(source_change)) source_change=force_ledger
    minimum_theta=theta;particle_steps=first_steps+second_steps;ok=.true.
  contains
    subroutine force_kick()
      real(dp) :: local(nv),global(nv)
      call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,gas,nv,comm,valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temp,1,comm,valid)
      if(.not.valid) return
      call local_force_kick(local)
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      valid=ierr==MPI_SUCCESS.and.agreed
      if(.not.valid) return
      call MPI_Allreduce(local,global,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      valid=ierr==MPI_SUCCESS
      if(valid) force_ledger=force_ledger+global
    end subroutine
    subroutine local_force_kick(local)
      real(dp), intent(out) :: local(nv)
      real(dp) :: dv(3),mass,momentum(3),energy,volume,velocity(3),prim(nv),checked_t,sound
      integer :: rank,status,node,i,j,k,p,n(3)
      valid=.false.;local=0;dv=.5_dp*dt*accel
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
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
          volume=geometry(node)%volume_fraction(i,j,k)*geometry(node)%dx*geometry(node)%dy*geometry(node)%dz
          local(2:4)=local(2:4)+volume*momentum;local(5)=local(5)+volume*energy
        end do
        end do
        end do
      end do
      do p=1,size(drops)
        if(drops(p)%mass<=0) cycle
        mass=drops(p)%mass*drops(p)%multiplicity;velocity=drops(p)%velocity
        drops(p)%velocity=velocity+dv
        if(.not.valid_spray_parcel(drops(p))) return
        local(2:4)=local(2:4)+mass*dv
        local(5)=local(5)+mass*dot_product(velocity+.5_dp*dv,dv)
      end do
      valid=all(ieee_is_finite(local))
    end subroutine
    subroutine particle_half(steps,exchange)
      integer, intent(out) :: steps
      real(dp), intent(inout) :: exchange(:)
      logical :: rejected
      if(present(last_id).and.present(films)) then
        call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
          drops,masses,lower,spacing,.5_dp*dt,steps,valid,last_id=watermark,boundary_change=exchange, &
          wall_films=patches,film_mass=patch_mass,film_settings=film_settings, &
          eb_geometry=geometry,nodal_phi=phi,model_rejected=rejected,eb_outflow_change=half_exterior,boundary_faces=particle_faces)
        if(present(model_rejected)) model_rejected=rejected
      else if(present(last_id)) then
        call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
          drops,masses,lower,spacing,.5_dp*dt,steps,valid,last_id=watermark,boundary_change=exchange, &
          eb_geometry=geometry,nodal_phi=phi,eb_outflow_change=half_exterior,boundary_faces=particle_faces)
      else
        call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
          drops,masses,lower,spacing,.5_dp*dt,steps,valid,boundary_change=exchange, &
          eb_geometry=geometry,nodal_phi=phi,eb_outflow_change=half_exterior,boundary_faces=particle_faces)
      end if
      if(valid.and.outflow) exterior=exterior+half_exterior
    end subroutine
  end subroutine
  ! Common-step gas/particle SDC. New kicks are reevaluated from one saved
  ! kick-drift endpoint, never cumulatively advanced across SDC iterations.
  ! Existing stationary films evolve thermally; particle contacts remain specular.
  subroutine advance_mpi_eb_spray_sdc(comm,tree,owners,geometry,extended_geometry,phi,species,reactions,transport, &
      liquids,options,state,temperature,parcels,component_mass,lower,spacing,dt,iterations,chemistry,transport_flags, &
      force,rtol,atol,boundary_faces,boundary_change,wall_change,chemical_change,source_change,minimum_theta, &
      particle_steps,last_id,ok,reaction_source,difmag,redistribution,target_volume_fraction,gas_particle_change, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature,split_created,particle_acceleration, &
      films,film_mass,film_settings,film_steps,film_impacts,model_rejected,gas_acceleration)
    type(spray_wall_film), allocatable, optional, intent(inout) :: films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    logical, optional, intent(in) :: film_impacts
    logical, optional, intent(out) :: model_rejected
    logical :: impacts
    integer, optional, intent(out) :: film_steps
    type(spray_wall_film), allocatable :: midpoint_films(:),final_films(:)
    real(dp), allocatable :: midpoint_film_mass(:,:),final_film_mass(:,:)
    integer :: film_flag,film_low,film_high,old_film_steps,new_film_steps
    real(dp), optional, intent(in) :: gas_acceleration(3)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),iterations
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),spacing(3),dt,force(3),rtol,atol
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    logical, intent(in) :: chemistry,transport_flags(4)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),minimum_theta
    integer, intent(out) :: particle_steps
    integer(int64), intent(inout) :: last_id
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), optional, intent(in) :: difmag,target_volume_fraction
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(inout) :: gas_particle_change(:)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: wall_temperature(6)
    real(dp), optional, intent(in) :: particle_acceleration(3)
    integer(int64), optional, intent(out) :: split_created
    type(amr_tree_state_3d) :: gas,temp,old_particle_rate
    type(amr_tree_state_3d), allocatable :: ir
    type(spray_parcel), allocatable :: midpoint(:),final(:)
    real(dp), allocatable :: midpoint_mass(:,:),final_mass(:,:)
    integer(int64) :: midpoint_id,final_id,before_split
    real(dp) :: accel(3),low_accel(3),high_accel(3)
    real(dp) :: old_force(size(source_change)),new_force(size(source_change))
    real(dp) :: low_multiplicity,high_multiplicity
    integer :: nv,ierr,old_steps,new_steps
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: external(size(source_change)),exchange(size(source_change)),theta
    real(dp) :: old_wall(size(source_change)),new_wall(size(source_change))
    real(dp) :: old_outflow(size(source_change)),new_outflow(size(source_change))
    character(len=32) :: particle_faces(6)
    logical :: valid,ready
    ok=.false.;particle_steps=0;nv=reactive_nvar(size(species))
    if(present(split_created)) split_created=0
    if(present(model_rejected)) model_rejected=.false.
    impacts=.false.
    if(present(film_impacts)) impacts=film_impacts
    if(present(film_steps)) film_steps=0
    film_flag=merge(1,0,present(films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))+8*merge(1,0,impacts)
    call MPI_Allreduce(film_flag,film_low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,film_high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.film_low/=film_high) return
    if(film_flag/=0.and.film_flag/=7.and.film_flag/=15) return
    if(present(split_created)) split_created=0
    accel=0
    if(present(particle_acceleration)) accel=particle_acceleration
    valid=allocated(parcels).and.allocated(component_mass).and.valid_spray_options(options)
    if(present(films)) valid=valid.and.allocated(films).and.allocated(film_mass)
    valid=valid.and.all([size(boundary_change),size(wall_change),size(chemical_change),size(source_change)]==nv)
    valid=valid.and.all(ieee_is_finite(accel)).and.all(abs(accel)<=1e6_dp)
    if(present(gas_particle_change)) valid=valid.and.size(gas_particle_change)==nv
    call MPI_Allreduce(valid,ready,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.ready) return
    call MPI_Allreduce(options%maximum_multiplicity,low_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(options%maximum_multiplicity,high_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.low_multiplicity/=high_multiplicity) return
    call MPI_Allreduce(accel,low_accel,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(accel,high_accel,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(low_accel/=high_accel)) return
    particle_faces=boundary_faces
    where(particle_faces=='inflow') particle_faces='outflow'
    where(particle_faces=='isothermal_wall') particle_faces='no_slip_wall'
    gas=state;temp=temperature;old_steps=0;new_steps=0
    old_film_steps=0;new_film_steps=0
    old_force=0;new_force=0
    old_wall=0;new_wall=0;old_outflow=0;new_outflow=0
    if(present(reaction_source)) then
      allocate(ir);ir=reaction_source
    end if
    call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,gas,temp,extended_geometry, &
      dt,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces,boundary,wall,chemical,external,theta, &
      comm,ready,ir,difmag,redistribution=redistribution,target_volume_fraction=target_volume_fraction, &
      auxiliary_source=particle_source,auxiliary_change=exchange,inflow_state=inflow_state, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      wall_temperature=wall_temperature,gas_acceleration=gas_acceleration)
    if(.not.ready) return
    before_split=final_id
    if(options%maximum_multiplicity>0) then
      call split_mpi_spray_parcels(comm,final,options%maximum_multiplicity,options%maximum_parcels,final_id,ready,final_mass)
      if(.not.ready) return
    end if
    valid=old_steps<=huge(1)-new_steps.and.old_film_steps<=huge(1)-new_film_steps
    call MPI_Allreduce(valid,ready,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.ready) return
    boundary=boundary+old_outflow+new_outflow;wall=wall+old_wall+new_wall
    external=external+old_force+new_force
    valid=all(ieee_is_finite(boundary)).and.all(ieee_is_finite(wall)).and.all(ieee_is_finite(external))
    call MPI_Allreduce(valid,ready,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.ready) return
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(final,parcels);call move_alloc(final_mass,component_mass)
    if(present(films)) then
      call move_alloc(final_films,films);call move_alloc(final_film_mass,film_mass)
    end if
    if(present(film_steps)) film_steps=old_film_steps+new_film_steps
    last_id=final_id
    if(present(split_created)) split_created=final_id-before_split
    if(present(reaction_source)) reaction_source=ir
    boundary_change=boundary;wall_change=wall;chemical_change=chemical;source_change=external
    if(present(gas_particle_change)) gas_particle_change=exchange
    minimum_theta=theta;particle_steps=old_steps+new_steps;ok=.true.
  contains
    ! External particle half kicks bracket the coupled interval. The second
    ! kick is evaluated only at the final iteration, as in pinned particleMK.
    subroutine particle_force(population,interval,ledger,force_comm,good)
      type(spray_parcel), intent(inout) :: population(:)
      real(dp), intent(in) :: interval
      real(dp), intent(out) :: ledger(:)
      type(MPI_Comm), intent(in) :: force_comm
      logical, intent(out) :: good
      real(dp) :: local(size(ledger)),dv(3),velocity(3),mass
      integer :: p,status
      logical :: valid,agreed
      local=0;ledger=0;good=.false.;valid=.true.;dv=interval*accel
      do p=1,size(population)
        if(population(p)%mass<=0) cycle
        velocity=population(p)%velocity;mass=population(p)%mass*population(p)%multiplicity
        population(p)%velocity=velocity+dv
        valid=valid.and.valid_spray_parcel(population(p))
        local(2:4)=local(2:4)+mass*dv
        local(5)=local(5)+mass*dot_product(velocity+.5_dp*dv,dv)
      end do
      valid=valid.and.all(ieee_is_finite(local))
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,force_comm,status)
      if(status/=MPI_SUCCESS.or..not.agreed) return
      call MPI_Allreduce(local,ledger,size(ledger),MPI_DOUBLE_PRECISION,MPI_SUM,force_comm,status)
      good=status==MPI_SUCCESS.and.all(ieee_is_finite(ledger))
    end subroutine particle_force
    subroutine particle_source(stage_gas,stage_temperature,interval,stage,depth,rate,source_comm,good)
      type(amr_tree_state_3d), intent(in) :: stage_gas,stage_temperature
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: source_comm
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: work_gas,work_temperature,change,film_change
      type(spray_parcel), allocatable :: trial(:)
      type(spray_wall_film), allocatable :: patches(:)
      real(dp), allocatable :: patch_mass(:,:)
      real(dp) :: film_wall(nv)
      integer :: film_count
      real(dp), allocatable :: masses(:,:)
      real(dp) :: wall_delta(nv),exterior(nv),half_interval,force_delta(nv)
      integer(int64) :: watermark
      integer :: steps,p,status
      character(len=16) :: motion
      logical :: accepted,all_ok,rejected
      good=.false.
      if(depth/=-1) return
      if(stage>0.and.stage<iterations) then
        rate=old_particle_rate;good=.true.;return
      end if
      half_interval=.5_dp*interval
      if(half_interval<=0) return
      work_gas=stage_gas;work_temperature=stage_temperature
      if(stage==0) then
        trial=parcels;masses=component_mass;watermark=last_id;motion='kick_drift'
      else
        trial=midpoint;masses=midpoint_mass;watermark=midpoint_id;motion='kick'
      end if
      film_count=0;film_wall=0
      if(present(films)) then
        if(stage==0) then
          patches=films;patch_mass=film_mass
        else
          patches=midpoint_films;patch_mass=midpoint_film_mass
        end if
      end if
      force_delta=0
      if(stage==0) then
        call particle_force(trial,half_interval,force_delta,source_comm,accepted)
        if(.not.accepted) return
      end if
      if(impacts) then
      call advance_mpi_amr_mixture_coupling(source_comm,species,transport,liquids,options,work_gas,work_temperature, &
        tree,owners,trial,masses,lower,spacing,interval,steps,accepted,last_id=watermark, &
        boundary_change=wall_delta,eb_geometry=geometry,nodal_phi=phi,eb_outflow_change=exterior, &
        boundary_faces=particle_faces,gas_change=change,motion_mode=motion, &
        wall_films=patches,film_mass=patch_mass,film_settings=film_settings,model_rejected=rejected)
        if(present(model_rejected)) model_rejected=rejected
      else
      call advance_mpi_amr_mixture_coupling(source_comm,species,transport,liquids,options,work_gas,work_temperature, &
        tree,owners,trial,masses,lower,spacing,interval,steps,accepted,last_id=watermark, &
        boundary_change=wall_delta,eb_geometry=geometry,nodal_phi=phi,eb_outflow_change=exterior, &
        boundary_faces=particle_faces,gas_change=change,motion_mode=motion)
      end if
      if(.not.accepted) return
      if(stage>0) then
        call particle_force(trial,half_interval,force_delta,source_comm,accepted)
        if(.not.accepted) return
      end if
      if(present(films)) then
        call advance_mpi_amr_film_coupling(source_comm,species,transport,liquids,options,film_settings, &
          work_gas,work_temperature,tree,owners,patches,patch_mass,lower,spacing,half_interval, &
          particle_faces,film_count,film_wall,accepted,geometry,phi,film_change)
        if(.not.accepted) return
        accepted=.true.
        do p=1,size(owners)
          if(.not.allocated(change%fields(p)%values)) cycle
          associate(a=>change%fields(p)%values,b=>film_change%fields(p)%values)
            if(any(a>huge(1._dp)-max(b,0._dp)).or.any(a< -huge(1._dp)-min(b,0._dp))) then
              accepted=.false.;exit
            end if
            a=a+b
          end associate
        end do
        call MPI_Allreduce(accepted,all_ok,1,MPI_LOGICAL,MPI_LAND,source_comm,status)
        if(status/=MPI_SUCCESS.or..not.all_ok) return
        wall_delta=wall_delta+film_wall
      end if
      accepted=.true.
      do p=1,size(owners)
        if(.not.allocated(change%fields(p)%values)) cycle
        if(half_interval<1) then
          if(any(abs(change%fields(p)%values)>huge(1._dp)*half_interval)) accepted=.false.
        end if
      end do
      call MPI_Allreduce(accepted,all_ok,1,MPI_LOGICAL,MPI_LAND,source_comm,status)
      if(status/=MPI_SUCCESS.or..not.all_ok) return
      do p=1,size(owners)
        if(.not.allocated(change%fields(p)%values)) cycle
        change%fields(p)%values=change%fields(p)%values/half_interval
      end do
      rate=change
      if(stage==0) then
        old_particle_rate=change
        call move_alloc(trial,midpoint);call move_alloc(masses,midpoint_mass)
        if(present(films)) then
          call move_alloc(patches,midpoint_films);call move_alloc(patch_mass,midpoint_film_mass)
        end if
        old_film_steps=film_count
        midpoint_id=watermark;old_steps=steps;old_wall=wall_delta;old_outflow=exterior;old_force=force_delta
      else
        call move_alloc(trial,final);call move_alloc(masses,final_mass)
        if(present(films)) then
          call move_alloc(patches,final_films);call move_alloc(patch_mass,final_film_mass)
        end if
        new_film_steps=film_count
        final_id=watermark;new_steps=steps;new_wall=wall_delta;new_outflow=exterior;new_force=force_delta
      end if
      good=.true.
    end subroutine particle_source
  end subroutine advance_mpi_eb_spray_sdc
end module
