module mpi_eb_tree_full_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use reactive_3d_mod, only: advance_reactive_chemistry_3d
  use reactive_eb_3d_driver_mod, only: react_sdc_eb_state_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,extract_eb_geometry_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d,average_down_mpi_eb_tree_3d, &
    validate_mpi_eb_tree_geometry_3d
  use mpi_eb_tree_hydro_3d_mod, only: advance_mpi_eb_tree_hydro_3d,compute_sdc_mpi_eb_tree_hydro_3d, &
    compute_mol_mpi_eb_tree_hydro_3d
  use mpi_eb_tree_transport_3d_mod, only: advance_mpi_eb_tree_transport_3d,compute_sdc_mpi_eb_tree_transport_3d
  use mpi_eb_tree_subcycle_hydro_3d_mod, only: advance_mpi_eb_tree_subcycle_hydro_3d
  use mpi_eb_tree_subcycle_transport_3d_mod, only: advance_mpi_eb_tree_subcycle_transport_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,valid_eb_integrated_flux_3d
  use eb_reactive_wall_flux_3d_mod, only: reactive_eb_slip_wall_source_3d
  use mpi_eb_tree_redistribution_3d_mod, only: redistribute_sdc_mpi_eb_tree_flux_rate_3d, &
    redistribute_sdc_mpi_eb_tree_level_rate_3d,redistribute_mpi_eb_tree_state_3d, &
    redistribute_sdc_mpi_eb_tree_state_rate_3d,redistribute_sdc_mpi_eb_tree_state_level_rate_3d
  use mpi_eb_tree_flux_3d_mod, only: synchronize_mpi_eb_tree_flux_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  implicit none
  private
  public :: react_mpi_eb_tree_3d,advance_mpi_eb_tree_full_3d,react_sdc_mpi_eb_tree_3d
  public :: advance_mpi_eb_tree_sdc_3d,advance_mpi_eb_tree_sdc_subcycle_3d
  ! Auxiliary sources receive the old state at stage 0 and provisional new
  ! states at stages 1..iterations. Rates are per fluid volume per time.
  ! Callback-owned physical state must remain private until the driver succeeds.
  public :: mpi_eb_sdc_diffusion_source
  abstract interface
    ! Rates per fluid volume/time; flux integrated over interval per open face area.
    subroutine mpi_eb_sdc_diffusion_source(state,temperature,interval,depth,rate,flux,theta,comm,ok, &
        ancestor_start,ancestor_end,temperature_start,temperature_end,fractions)
      import :: amr_tree_state_3d,eb_integrated_flux_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: state,temperature
      real(dp), intent(in) :: interval
      integer, intent(in) :: depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: theta
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
      type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
      real(dp), optional, intent(in) :: fractions(:)
    end subroutine

    ! Called after descendants and gas synchronization, before publication.
    ! Callback state must remain private until the entire outer step succeeds.
    subroutine mpi_eb_sdc_post_level(depth,time,interval,iteration,ncycle,comm,ok)
      import :: dp,MPI_Comm
      integer, intent(in) :: depth,iteration,ncycle
      real(dp), intent(in) :: time,interval
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    ! Build the complete three-cell predictor halo for iteration 1..iterations.
    ! Valid fluid cells must equal rate; callback state remains caller-owned.
    subroutine mpi_eb_sdc_predictor_builder(rate,interval,iteration,depth,halo,comm,ok)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: rate
      real(dp), intent(in) :: interval
      integer, intent(in) :: iteration,depth
      type(amr_tree_state_3d), intent(inout) :: halo
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    ! Recursive predictor donors use the midpoint of the level interval.
    subroutine mpi_eb_sdc_recursive_predictor_builder(rate,interval,iteration,depth,time,halo,comm,ok, &
        source_start,source_end,source_fraction)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: rate,source_start,source_end
      real(dp), intent(in) :: interval,time,source_fraction(:)
      integer, intent(in) :: iteration,depth
      type(amr_tree_state_3d), intent(inout) :: halo
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    subroutine mpi_eb_sdc_auxiliary_source(state,temperature,interval,stage,depth,rate,comm,ok)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: state,temperature
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    ! Time is the level interval start relative to the outer step, for all stages.
    subroutine mpi_eb_sdc_recursive_source(state,temperature,interval,stage,depth,time,rate,comm,ok)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: state,temperature
      real(dp), intent(in) :: interval,time
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    subroutine mpi_eb_sdc_recursive_halo_source(state,temperature,interval,stage,depth,time,rate,comm,ok, &
        ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: state,temperature
      type(amr_tree_state_3d), intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
      real(dp), intent(in) :: interval,time,ancestor_fraction(:)
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
  end interface
  public :: mpi_eb_sdc_post_level
  public :: mpi_eb_sdc_predictor_builder,mpi_eb_sdc_recursive_predictor_builder
  public :: mpi_eb_sdc_auxiliary_source,mpi_eb_sdc_recursive_source,mpi_eb_sdc_recursive_halo_source
contains
  subroutine agree(condition,comm,ok)
    logical, intent(in) :: condition
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: status
    call MPI_Allreduce(condition,ok,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if (status/=MPI_SUCCESS) ok=.false.
  end subroutine

  subroutine react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old_state,old_temperature, &
      nonreacting_rates,dt,rtol,atol,state,temperature,reaction_source,comm,ok,target_depth, &
      evaluate_reaction)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: old_state,old_temperature,nonreacting_rates
    real(dp), intent(in) :: dt,rtol,atol
    type(amr_tree_state_3d), intent(inout) :: state,temperature,reaction_source
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d) :: candidate,temps,sources
    real(dp) :: controls(3),low(3),high(3),stamp
    integer :: nv,depth,counts(4),min_counts(4),max_counts(4),rank,status,node,n(3)
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species));depth=-1;controls=[dt,rtol,atol]
    valid=all(ieee_is_finite(controls)).and.min(dt,rtol,atol)>0.and.size(reactions)>0
    if(present(target_depth)) then
      depth=target_depth
      valid=valid.and.depth>=0
    end if
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(controls,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    counts=[size(species),size(reactions),depth,merge(1,0,present(evaluate_reaction))]
    call MPI_Allreduce(counts,min_counts,4,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(counts,max_counts,4,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(min_counts/=max_counts)) return
    call validate_mpi_amr_tree_local_3d(tree,owners,old_state,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,old_temperature,1,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,nonreacting_rates,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,reaction_source,nv-4,comm,ready)
    if(.not.ready) return
    call agree(depth<tree%level_count().and.size(geometry)==size(tree%nodes),comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    candidate=old_state;temps=old_temperature;sources=reaction_source
    if(depth>=0) then
      call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ready)
      if(.not.ready) return
      candidate=state;temps=temperature
    end if
    valid=.true.
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
      n=tree%nodes(node)%shape
      if(.not.geometry(node)%is_valid()) then
        valid=.false.;exit
      end if
      if(any([geometry(node)%nx,geometry(node)%ny,geometry(node)%nz]/=n)) then
        valid=.false.;exit
      end if
      stamp=tree_work_start()
      call react_sdc_eb_state_3d(species,reactions,old_state%fields(node)%values, &
        old_temperature%fields(node)%values(1,:,:,:),geometry(node),nonreacting_rates%fields(node)%values, &
        dt,rtol,atol,candidate%fields(node)%values,temps%fields(node)%values(1,:,:,:), &
        sources%fields(node)%values,valid, &
          evaluate_reaction=evaluate_reaction)
      call tree_work_stop(node,stamp)
      if(.not.valid) exit
    end do
    call agree(valid,comm,ready)
    if(.not.ready) return
    state=candidate;temperature=temps;reaction_source=sources;ok=.true.
  end subroutine react_sdc_mpi_eb_tree_3d

  subroutine react_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,dt,rtol,atol,state,temperature, &
      chemical_change,comm,ok,integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(in) :: dt,rtol,atol
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(inout) :: chemical_change(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    type(amr_tree_state_3d) :: candidate,ct
    real(dp) :: work_stamp
    real(dp) :: before(size(chemical_change)),after(size(chemical_change)),parameters(3),low(3),high(3)
    logical, allocatable :: active(:,:,:)
    integer :: nv,rank,status,node,child,n(3),lo(3),hi(3),ios
    character(len=32) :: method,reference
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species));parameters=[dt,rtol,atol]
    valid=all(ieee_is_finite(parameters))
    if (valid) valid=all(parameters>0._dp)
    valid=valid.and.size(chemical_change)==nv.and.size(reactions)>0
    method='@default'
    if (present(integrator)) then
      valid=valid.and.len_trim(integrator)<=len(method)
      method=integrator
    end if
    call agree(valid,comm,ready)
    if (.not.ready) return
    call MPI_Allreduce(parameters,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(parameters,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.any(low/=high)) return
    reference=method
    call MPI_Bcast(reference,len(reference),MPI_CHARACTER,0,comm,status)
    call agree(status==MPI_SUCCESS.and.method==reference,comm,ready)
    if (.not.ready) return
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,state,before,comm,ready)
    if (.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    candidate=state;ct=temperature
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(active(n(1),n(2),n(3)),stat=ios)
      valid=ios==0
      if (.not.valid) exit
      active=geometry(node)%volume_fraction>0._dp
      do child=node+1,size(owners)
        if (tree%nodes(child)%parent/=node) cycle
        associate(p=>tree%nodes(child)%patch)
        lo=[p%coarse_i_lower,p%coarse_j_lower,p%coarse_k_lower]
        hi=[p%coarse_i_upper,p%coarse_j_upper,p%coarse_k_upper]
        end associate
        active(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      work_stamp=tree_work_start()
      if (any(active)) call advance_reactive_chemistry_3d(species,reactions,candidate%fields(node)%values, &
        ct%fields(node)%values(1,:,:,:),n(1),n(2),n(3),dt,rtol,atol,valid,active,integrator)
      call tree_work_stop(node,work_stamp)
      deallocate(active)
      if (.not.valid) exit
    end do
    call agree(valid,comm,ready)
    if (.not.ready) return
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,ready)
    if (.not.ready) return
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,candidate,after,comm,ready)
    if (.not.ready) return
    call agree(all(ieee_is_finite(after-before)),comm,ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(ct%fields,temperature%fields)
    chemical_change=after-before;ok=.true.
  end subroutine

  ! Transactional R/2 T/2 H T/2 R/2, common or recursive H/T, owner-local fields.
  subroutine advance_mpi_eb_tree_full_3d(species,reactions,transport,tree,owners,geometry,solver,dt, &
      chemistry_enabled,rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion, &
      state,temperature,boundary_change,wall_change,chemical_change,minimum_theta,comm,ok,integrator,subcycling, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: chemistry_enabled,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),minimum_theta
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    real(dp), optional, intent(in) :: difmag
    real(dp) :: coefficient,coefficient_low,coefficient_high
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    logical, optional, intent(in) :: subcycling
    logical :: recursive_steps
    type(amr_tree_state_3d) :: candidate,ct
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: exchange(size(boundary_change)),transport_sum(size(boundary_change)),theta,stage_theta
    integer :: flags(3),low(3),high(3),status,nv,half
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    valid=ieee_is_finite(dt).and.ieee_is_finite(coefficient).and.coefficient>=0
    if (valid) valid=dt>0._dp
    if(present(wall_temperature)) then
      valid=valid.and.present(boundary_faces)
      valid=valid.and.all(ieee_is_finite(wall_temperature)).and.all(wall_temperature>0)
    end if
    if(present(boundary_faces)) then
      if(any(boundary_faces=='isothermal_wall')) then
        valid=valid.and.present(wall_temperature).and.transport_enabled.and.conduction
      end if
    end if
    call agree(valid.and.size(boundary_change)==nv.and.size(wall_change)==nv.and.size(chemical_change)==nv,comm,ready)
    if (.not.ready) return
    call MPI_Allreduce(coefficient,coefficient_low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(coefficient,coefficient_high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.coefficient_low/=coefficient_high) return
    recursive_steps=.false.
    if (present(subcycling)) recursive_steps=subcycling
    flags=merge(1,0,[chemistry_enabled,transport_enabled,recursive_steps])
    call MPI_Allreduce(flags,low,3,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,3,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.any(low/=high)) return
    candidate=state;ct=temperature
    chemical=0._dp;transport_sum=0._dp;theta=1._dp
    do half=1,2
      if (half==1.and.chemistry_enabled) then
        call reaction_half(valid)
        if (.not.valid) return
      end if
      if (transport_enabled) then
        if (recursive_steps) then
          call advance_mpi_eb_tree_subcycle_transport_3d(species,transport,tree,owners,geometry,0.5_dp*dt, &
            viscosity,conduction,diffusion,barodiffusion,candidate,ct,exchange,stage_theta,comm,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,wall_temperature=wall_temperature)
        else
          call advance_mpi_eb_tree_transport_3d(species,transport,tree,owners,geometry,0.5_dp*dt, &
            viscosity,conduction,diffusion,barodiffusion,candidate,ct,exchange,stage_theta,comm,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,wall_temperature=wall_temperature)
        end if
        if (.not.valid) return
        transport_sum=transport_sum+exchange;theta=min(theta,stage_theta)
      end if
      if (half==1) then
        if (recursive_steps) then
          call advance_mpi_eb_tree_subcycle_hydro_3d(species,tree,owners,geometry,solver,dt,candidate,ct, &
            boundary,wall,comm,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state, &
            difmag=coefficient,extended_geometry=extended_geometry)
        else
          call advance_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,solver,dt,candidate,ct,boundary,wall,comm,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state, &
            difmag=coefficient,extended_geometry=extended_geometry)
        end if
        if (.not.valid) return
      else if (chemistry_enabled) then
        call reaction_half(valid)
        if (.not.valid) return
      end if
    end do
    boundary=boundary+transport_sum
    call agree(all(ieee_is_finite(boundary)).and.all(ieee_is_finite(chemical)),comm,valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(ct%fields,temperature%fields)
    boundary_change=boundary;wall_change=wall;chemical_change=chemical;minimum_theta=theta;ok=.true.
  contains
    subroutine reaction_half(good)
      logical, intent(out) :: good
      call react_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,0.5_dp*dt,rtol,atol, &
        candidate,ct,exchange,comm,good,integrator)
      if (good) chemical=chemical+exchange
    end subroutine
  end subroutine
  ! Common-step or selected-level SDC; stages and publication are collective.
  ! Selected levels export fluxes without initial restriction or final reflux.
  subroutine advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry, &
      state,temperature,extended_geometry,dt,iterations,chemistry,transport_flags,force,rtol,atol, &
      boundary_faces,boundary_change,wall_change,chemical_change,source_change,minimum_theta,comm,ok, &
      reaction_source,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
      target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start,fractions_end, &
      source_start,source_end,integrated_flux,old_predictor,new_predictor,redistribution,redistribution_transfer, &
      target_volume_fraction,auxiliary_source,auxiliary_change,gas_acceleration,predictor_builder,evaluate_diffusion, &
      evaluate_reaction,embedded_wall_temperature,embedded_noslip,mol_mode)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(mpi_eb_sdc_diffusion_source), optional :: evaluate_diffusion
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),iterations
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), intent(in) :: dt,force(3),rtol,atol
    logical, intent(in) :: chemistry,transport_flags(4)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),minimum_theta
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: fractions_start(:),fractions_end(:)
    type(amr_tree_state_3d), optional, intent(in) :: source_start,source_end
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: integrated_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: old_predictor,new_predictor
    type(amr_tree_state_3d), optional, intent(inout) :: redistribution_transfer
    type(amr_tree_state_3d) :: transfer0,transfer1,hydro_transfer,final_transfer
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(in) :: target_volume_fraction
    procedure(mpi_eb_sdc_auxiliary_source), optional :: auxiliary_source
    real(dp), optional, intent(inout) :: auxiliary_change(:)
    procedure(mpi_eb_sdc_predictor_builder), optional :: predictor_builder
    type(amr_tree_state_3d), allocatable :: predictor_halo
    type(amr_tree_state_3d) :: a0,a1
    real(dp) :: auxiliary(size(source_change))
    real(dp), optional, intent(in) :: gas_acceleration(3)
    real(dp) :: accel(3)
    real(dp) :: redist_target
    character(len=16) :: redistribution_method
    type(amr_tree_state_3d) :: redistributed,before_sync,sync_rate
    type(amr_tree_state_3d) :: op,np
    real(dp), allocatable :: midpoint(:),fraction_low(:),fraction_high(:)
    integer :: selected_depth,selection(6),minselection(6),maxselection(6),levels
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: wall_temperature(6),embedded_wall_temperature
    type(eb_geometry_3d), allocatable :: heat_geometry(:),mol_geometry(:)
    logical, optional, intent(in) :: embedded_noslip,mol_mode
    logical :: noslip,wall_active,mol
    real(dp), allocatable :: heat_stage,heat_frequency,wall_impulse(:)
    real(dp) :: impulse_old(3),impulse_new(3)
    real(dp) :: heat_old,heat_new,mol_wall_old(size(wall_change)),mol_wall_new(size(wall_change))
    type(amr_tree_state_3d) :: old,ot,candidate,ct,e0,e1,d0,d1,rates,hydro,ir,ledger
    type(eb_integrated_flux_3d), allocatable :: f0(:),f1(:),hf(:)
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: external(size(source_change)),v0(size(source_change)),v1(size(source_change)),theta,stage_theta
    real(dp) :: controls(15),low(15),high(15),recovered,speed
    real(dp), allocatable :: primitive(:),wall_rate(:,:,:,:)
    integer :: counts(11),mincounts(11),maxcounts(11),nv,rank,status,node,n(3),it,i,j,k,c,ios
    logical :: valid,ready,diffusive,has_diffusion
    ok=.false.;nv=reactive_nvar(size(species));theta=1
    has_diffusion=present(evaluate_diffusion)
    accel=0
    if(present(gas_acceleration)) accel=gas_acceleration
    redistribution_method='none'
    if(present(redistribution)) redistribution_method=redistribution
    redist_target=.5_dp
    if(present(target_volume_fraction)) redist_target=target_volume_fraction
    call agree(redistribution_method=='none'.or.redistribution_method=='flux_redist'.or. &
      redistribution_method=='state_redist',comm,ready)
    if(.not.ready) return
    controls(:11)=[dt,force,rtol,atol,merge(1._dp,0._dp,redistribution_method=='flux_redist')+ &
      merge(2._dp,0._dp,redistribution_method=='state_redist'),redist_target,accel]
    mol=.false.
    if(present(mol_mode)) mol=mol_mode
    controls(15)=merge(1._dp,0._dp,mol)
    mol_wall_old=0;mol_wall_new=0
    noslip=.false.
    if(present(embedded_noslip)) noslip=embedded_noslip
    wall_active=noslip.or.present(embedded_wall_temperature)
    controls(14)=merge(1._dp,0._dp,noslip)
    controls(12:13)=0
    if(present(embedded_wall_temperature)) controls(12:13)=[1._dp,embedded_wall_temperature]
    heat_old=0;heat_new=0;impulse_old=0;impulse_new=0
    valid=all(ieee_is_finite(controls)).and.min(dt,rtol,atol)>0.and.iterations>=1
    valid=valid.and.redist_target>0.and.redist_target<=1.and.all(abs(accel)<=1e6_dp)
    valid=valid.and.all([size(boundary_change),size(wall_change),size(chemical_change),size(source_change)]==nv)
    valid=valid.and..not.(transport_flags(4).and..not.transport_flags(3))
    valid=valid.and..not.(present(reaction_source).and..not.chemistry)
    valid=valid.and.(present(auxiliary_source).eqv.present(auxiliary_change))
    if(present(auxiliary_change)) valid=valid.and.size(auxiliary_change)==nv
    valid=valid.and.(chemistry.or..not.present(evaluate_reaction))
    if(present(embedded_wall_temperature)) valid=valid.and.transport_flags(2).and.embedded_wall_temperature>0
    if(noslip) valid=valid.and.transport_flags(1)
    ! Particle predictor callbacks still require their own MOL staging.
    if(mol) valid=valid.and..not.present(predictor_builder)
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(controls,low,15,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,15,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    counts=[iterations,merge(1,0,chemistry),merge(1,0,transport_flags),merge(1,0,present(reaction_source)), &
      merge(1,0,present(auxiliary_source)),merge(1,0,present(predictor_builder)),merge(1,0,present(evaluate_diffusion)), &
      merge(1,0,present(evaluate_reaction))]
    call MPI_Allreduce(counts,mincounts,11,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(counts,maxcounts,11,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(mincounts/=maxcounts)) return
    selected_depth=-1
    if(present(target_depth)) selected_depth=target_depth
    selection=[selected_depth,sum(merge(1,0,[present(ancestor_start),present(ancestor_end), &
      present(temperature_start),present(temperature_end),present(fractions_start),present(fractions_end)])), &
      sum(merge(1,0,[present(source_start),present(source_end)])),merge(1,0,present(integrated_flux)), &
      2*merge(1,0,present(old_predictor))+merge(1,0,present(new_predictor)), &
      merge(1,0,present(redistribution_transfer))]
    valid=selected_depth>=-1
    if(present(target_depth)) valid=selected_depth>=0.and.present(integrated_flux)
    valid=valid.and.(selection(2)==0.or.selection(2)==6).and.(selection(3)==0.or.selection(3)==2)
    if(selected_depth>=0.and.redistribution_method/='none') then
      valid=valid.and.present(redistribution_transfer)
    else
      valid=valid.and..not.present(redistribution_transfer)
    end if
    if(selection(2)>0) valid=valid.and.selected_depth>=0
    if(selection(3)>0) valid=valid.and.selection(2)==6
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(selection,minselection,6,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(selection,maxselection,6,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(minselection/=maxselection)) return
    old=state;ot=temperature
    if(selected_depth<0) then
      call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,old,ot,comm,ready)
      if(.not.ready) return
    else
      call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,old,nv,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,ot,1,comm,ready)
      if(.not.ready) return
      if(selected_depth>=tree%level_count()) return
    end if
    if(mol) then
      call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,ready,padding=6)
      if(.not.ready) return
      allocate(mol_geometry(size(owners)),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      do node=1,size(owners)
        n=tree%nodes(node)%shape
        call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],n+9,mol_geometry(node),valid)
        if(.not.valid) exit
      end do
      call agree(valid,comm,ready)
      if(.not.ready) return
    end if
    if(wall_active) then
      call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,ready,padding=6)
      if(.not.ready) return
      if(noslip) then
        allocate(wall_impulse(3),stat=ios)
        call agree(ios==0,comm,ready)
        if(.not.ready) return
      end if
      allocate(heat_geometry(size(owners)),heat_stage,heat_frequency,stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      do node=1,size(owners)
        n=tree%nodes(node)%shape
        call extract_eb_geometry_3d(extended_geometry(node),[5,5,5],n+8,heat_geometry(node),valid)
        if(.not.valid) exit
      end do
      call agree(valid,comm,ready)
      if(.not.ready) return
    end if
    if(selection(2)==6) then
      levels=tree%level_count()-1
      call agree(size(fractions_start)==levels.and.size(fractions_end)==levels,comm,ready)
      if(.not.ready) return
      call agree(all(ieee_is_finite(fractions_start)).and.all(ieee_is_finite(fractions_end)),comm,ready)
      if(.not.ready) return
      call agree(all(fractions_start>=0).and.all(fractions_end<=1).and.all(fractions_end>=fractions_start),comm,ready)
      if(.not.ready) return
      allocate(fraction_low(2*levels),fraction_high(2*levels),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      call MPI_Allreduce([fractions_start,fractions_end],fraction_low,2*levels,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce([fractions_start,fractions_end],fraction_high,2*levels,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(fraction_low/=fraction_high)) return
      if(selection(3)==2) midpoint=.5_dp*(fractions_start+fractions_end)
    end if
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    allocate(ir%fields(size(owners)),primitive(reactive_nprim(size(species))),stat=ios)
    call agree(ios==0,comm,ready)
    if(.not.ready) return
    if(present(reaction_source)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,reaction_source,nv-4,comm,ready)
      if(.not.ready) return
      ir=reaction_source
    else
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)),source=0._dp,stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        end if
      end do
      call agree(valid,comm,ready)
      if(.not.ready) return
    end if
    do node=1,size(owners)
      if(owners(node)/=rank.or..not.active(node)) cycle
      do c=1,nv-4
        where(geometry(node)%volume_fraction<=0) ir%fields(node)%values(c,:,:,:)=0
      end do
      if(.not.all(ieee_is_finite(ir%fields(node)%values))) valid=.false.
    end do
    call agree(valid,comm,ready)
    if(.not.ready) return
    candidate=old;ct=ot;diffusive=any(transport_flags(:3))
    call sources(old,ot,e0,d0,a0,0,f0,transfer0,ready,fractions_start)
    if(.not.ready) return
    e1=e0;d1=d0;heat_new=heat_old;impulse_new=impulse_old;mol_wall_new=mol_wall_old
    if(present(auxiliary_source)) a1=a0
    call predictor_field(e0,d0,a0,op,ready)
    if(.not.ready) return
    if(mol) then
      hydro=d0;hf=f0
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        hydro%fields(node)%values=0
        if(.not.active(node)) cycle
        hf(node)%x=0;hf(node)%y=0;hf(node)%z=0
      end do
      if(present(redistribution_transfer)) hydro_transfer=hydro
    end if
    do it=1,merge(1,iterations,mol.and..not.chemistry)
      ! Later reacting MOL iterations re-evaluate the previous reacted state;
      ! only the first iteration forms the Euler predictor.
      if(.not.mol.or.it==1) then
      call combined_rates(.true.)
      if(present(predictor_builder)) then
        if(allocated(predictor_halo)) deallocate(predictor_halo)
        allocate(predictor_halo,stat=ios)
        call agree(ios==0,comm,ready)
        if(.not.ready) return
        call predictor_builder(rates,dt,it,selected_depth,predictor_halo,comm,valid)
        call agree(valid,comm,ready)
        if(.not.ready) return
      end if
      if(.not.mol) then
      call compute_sdc_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,old,ot,rates,extended_geometry, &
        dt,hydro,hf,comm,ready,difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start, &
        source_start,source_end,midpoint,predictor_halo)
      if(.not.ready) return
      if(redistribution_method=='state_redist') then
        if(selected_depth>=0) then
          call redistribute_sdc_mpi_eb_tree_state_level_rate_3d(tree,owners,geometry,old,hydro,nv,selected_depth,dt, &
            redistributed,hydro_transfer,comm,ready,redist_target,boundary_faces)
        else
          call redistribute_sdc_mpi_eb_tree_state_rate_3d(tree,owners,geometry,old,hydro,nv,dt,redistributed, &
            comm,ready,redist_target,boundary_faces)
        end if
        if(.not.ready) return
        hydro=redistributed
      else if(redistribution_method=='flux_redist') then
        if(selected_depth>=0) then
          call redistribute_sdc_mpi_eb_tree_level_rate_3d(tree,owners,geometry,hydro,nv,selected_depth,dt, &
            redistributed,hydro_transfer,comm,ready,boundary_faces)
        else
          call redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree,owners,geometry,hydro,nv,redistributed, &
            comm,ready,boundary_faces)
        end if
        if(.not.ready) return
        hydro=redistributed
      end if
      end if
      do node=1,size(owners)
        if(owners(node)/=rank.or..not.active(node)) cycle
        candidate%fields(node)%values=old%fields(node)%values+dt*(rates%fields(node)%values+hydro%fields(node)%values)
      end do
      call recover(ready)
      if(.not.ready) return
      end if
      call sources(candidate,ct,e1,d1,a1,it,f1,transfer1,ready,fractions_end)
      if(.not.ready) return
      call combined_rates(.false.)
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        rates%fields(node)%values=rates%fields(node)%values+hydro%fields(node)%values
      end do
      if(chemistry) then
        call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old,ot,rates,dt,rtol,atol, &
          candidate,ct,ir,comm,ready,target_depth, &
          evaluate_reaction=evaluate_reaction)
        if(.not.ready) return
      else
        do node=1,size(owners)
          if(owners(node)/=rank.or..not.active(node)) cycle
          candidate%fields(node)%values=old%fields(node)%values+dt*rates%fields(node)%values
        end do
        call recover(ready)
        if(.not.ready) return
      end if
    end do
    ledger=old
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      ledger%fields(node)%values=0
      if(.not.active(node)) cycle
      ledger%fields(node)%values=candidate%fields(node)%values-(old%fields(node)%values+dt*rates%fields(node)%values)
    end do
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,ledger,chemical,comm,ready)
    if(.not.ready) return
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,e0,v0,comm,ready)
    if(.not.ready) return
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,e1,v1,comm,ready)
    if(.not.ready) return
    external=.5_dp*dt*(v0+v1)
    if(present(auxiliary_source)) then
      call integrate_mpi_eb_tree_3d(tree,owners,geometry,a0,v0,comm,ready)
      if(.not.ready) return
      call integrate_mpi_eb_tree_3d(tree,owners,geometry,a1,v1,comm,ready)
      if(.not.ready) return
      auxiliary=.5_dp*dt*(v0+v1)
      call agree(all(ieee_is_finite(auxiliary)),comm,ready)
      if(.not.ready) return
    end if
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      if(.not.active(node)) cycle
      n=tree%nodes(node)%shape
      allocate(wall_rate(nv,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      call reactive_eb_slip_wall_source_3d(species,old%fields(node)%values,ot%fields(node)%values(1,:,:,:), &
        geometry(node),wall_rate,valid)
      if(.not.valid) exit
      ledger%fields(node)%values=dt*wall_rate
      deallocate(wall_rate)
      if(diffusive.or.has_diffusion.or.mol) then
        hf(node)%x=hf(node)%x+.5_dp*(f0(node)%x+f1(node)%x)
        hf(node)%y=hf(node)%y+.5_dp*(f0(node)%y+f1(node)%y)
        hf(node)%z=hf(node)%z+.5_dp*(f0(node)%z+f1(node)%z)
      end if
    end do
    call agree(valid,comm,ready)
    if(.not.ready) return
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,ledger,wall,comm,ready)
    if(.not.ready) return
    if(mol) wall=.5_dp*(mol_wall_old+mol_wall_new)
    wall(5)=wall(5)+.5_dp*(heat_old+heat_new)
    wall(2:4)=wall(2:4)+.5_dp*(impulse_old+impulse_new)
    call agree(all(ieee_is_finite(wall)),comm,ready)
    if(.not.ready) return
    if(selected_depth<0) then
      if(redistribution_method/='none') before_sync=candidate
      call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,hf,comm,ready, &
        boundary_faces=boundary_faces)
      if(.not.ready) return
      if(redistribution_method/='none') then
        sync_rate=hydro
        do node=1,size(owners)
          if(owners(node)/=rank) cycle
          sync_rate%fields(node)%values=(candidate%fields(node)%values-before_sync%fields(node)%values)/dt
        end do
        if(redistribution_method=='state_redist') then
          call redistribute_mpi_eb_tree_state_3d(tree,owners,geometry,sync_rate,nv,redistributed, &
            comm,ready,redist_target,boundary_faces)
        else
          call redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree,owners,geometry,sync_rate,nv,redistributed, &
            comm,ready,boundary_faces)
        end if
        if(.not.ready) return
        do node=1,size(owners)
          if(owners(node)/=rank) cycle
          candidate%fields(node)%values=before_sync%fields(node)%values+dt*redistributed%fields(node)%values
        end do
        call recover(ready)
        if(.not.ready) return
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,ready)
        if(.not.ready) return
      end if
    end if
    boundary=0
    if(rank==owners(1).and.selected_depth<=0) then
      n=tree%nodes(1)%shape
      associate(g=>geometry(1),f=>hf(1))
      do c=1,nv
        boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
          g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
          g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
          g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
          g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
          g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
      end do
      end associate
    end if
    call MPI_Bcast(boundary,nv,MPI_DOUBLE_PRECISION,owners(1),comm,status)
    call agree(status==MPI_SUCCESS.and.all(ieee_is_finite(boundary)).and.all(ieee_is_finite(external)),comm,ready)
    if(.not.ready) return
    call predictor_field(e1,d1,a1,np,ready)
    if(.not.ready) return
    if(present(redistribution_transfer)) then
      final_transfer=hydro_transfer
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        final_transfer%fields(node)%values=hydro_transfer%fields(node)%values+ &
          .5_dp*(transfer0%fields(node)%values+transfer1%fields(node)%values)
        valid=valid.and.all(ieee_is_finite(final_transfer%fields(node)%values))
      end do
      call agree(valid,comm,ready)
      if(.not.ready) return
    end if
    if(present(reaction_source)) then
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        do c=1,nv-4
          where(geometry(node)%volume_fraction<=0) &
            ir%fields(node)%values(c,:,:,:)=reaction_source%fields(node)%values(c,:,:,:)
        end do
      end do
      reaction_source=ir
    end if
    if(present(auxiliary_change)) auxiliary_change=auxiliary
    if(present(redistribution_transfer)) redistribution_transfer=final_transfer
    if(present(integrated_flux)) integrated_flux=hf
    if(present(old_predictor)) old_predictor=op
    if(present(new_predictor)) new_predictor=np
    call move_alloc(candidate%fields,state%fields);call move_alloc(ct%fields,temperature%fields)
    boundary_change=boundary;wall_change=wall;chemical_change=chemical;source_change=external
    minimum_theta=theta;ok=.true.
  contains
    logical function active(p)
      integer, intent(in) :: p
      active=selected_depth<0.or.tree%nodes(p)%level==selected_depth
    end function
    subroutine predictor_field(ext,diff,aux,output,good)
      type(amr_tree_state_3d), intent(in) :: ext,diff,aux
      type(amr_tree_state_3d), intent(out) :: output
      logical, intent(out) :: good
      logical :: finite
      integer :: p
      output=ext;finite=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        output%fields(p)%values=ext%fields(p)%values+diff%fields(p)%values
        if(present(auxiliary_source)) output%fields(p)%values=output%fields(p)%values+aux%fields(p)%values
        output%fields(p)%values(6:,:,:,:)=output%fields(p)%values(6:,:,:,:)+ir%fields(p)%values(:nv-5,:,:,:)
        output%fields(p)%values(5,:,:,:)=output%fields(p)%values(5,:,:,:)+ir%fields(p)%values(nv-4,:,:,:)
        if(.not.all(ieee_is_finite(output%fields(p)%values))) finite=.false.
      end do
      call agree(finite,comm,good)
    end subroutine
    subroutine combined_rates(with_reaction)
      logical, intent(in) :: with_reaction
      integer :: p
      rates=e0
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        rates%fields(p)%values=.5_dp*(e0%fields(p)%values+e1%fields(p)%values+d0%fields(p)%values+d1%fields(p)%values)
        if(present(auxiliary_source)) rates%fields(p)%values=rates%fields(p)%values+ &
          .5_dp*(a0%fields(p)%values+a1%fields(p)%values)
        if(.not.with_reaction) cycle
        rates%fields(p)%values(6:,:,:,:)=rates%fields(p)%values(6:,:,:,:)+ir%fields(p)%values(:nv-5,:,:,:)
        rates%fields(p)%values(5,:,:,:)=rates%fields(p)%values(5,:,:,:)+ir%fields(p)%values(nv-4,:,:,:)
      end do
    end subroutine
    subroutine sources(gas,temp,external_rate,transport_rate,auxiliary_rate,stage,flux,exchange,good,fractions)
      type(amr_tree_state_3d), intent(in) :: gas,temp
      type(amr_tree_state_3d), intent(out) :: external_rate,transport_rate,auxiliary_rate,exchange
      integer, intent(in) :: stage
      type(eb_integrated_flux_3d), allocatable, intent(out) :: flux(:)
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fractions(:)
      integer :: p,d(3),x,y,z,alloc_status
      logical :: local_ok,all_ok
      type(amr_tree_state_3d) :: extra
      type(eb_integrated_flux_3d), allocatable :: extra_flux(:)
      real(dp) :: extra_theta,global_theta
      local_ok=.true.;good=.false.
      allocate(external_rate%fields(size(owners)),stat=alloc_status)
      call agree(alloc_status==0,comm,all_ok)
      if(.not.all_ok) return
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        d=tree%nodes(p)%shape
        allocate(external_rate%fields(p)%values(nv,d(1),d(2),d(3)),source=0._dp,stat=alloc_status)
        if(alloc_status/=0) then
          local_ok=.false.;exit
        end if
        if(.not.active(p)) cycle
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          if(.not.all(ieee_is_finite(gas%fields(p)%values(:,x,y,z)))) then
            local_ok=.false.;cycle
          end if
          if(gas%fields(p)%values(1,x,y,z)<=0) then
            local_ok=.false.;cycle
          end if
          ! force is force per volume; accel is acceleration evaluated at this stage.
          external_rate%fields(p)%values(2:4,x,y,z)=force+gas%fields(p)%values(1,x,y,z)*accel
          external_rate%fields(p)%values(5,x,y,z)=dot_product( &
            external_rate%fields(p)%values(2:4,x,y,z),gas%fields(p)%values(2:4,x,y,z))/ &
            gas%fields(p)%values(1,x,y,z)
        end do
        end do
        end do
      end do
      call agree(local_ok,comm,all_ok)
      if(.not.all_ok) return
      transport_rate=external_rate
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        transport_rate%fields(p)%values=0
      end do
      if(present(auxiliary_source)) then
        auxiliary_rate=transport_rate
        call auxiliary_source(gas,temp,dt,stage,selected_depth,auxiliary_rate,comm,local_ok)
        call agree(local_ok,comm,all_ok)
        if(.not.all_ok) return
        call validate_mpi_amr_tree_local_3d(tree,owners,auxiliary_rate,nv,comm,all_ok)
        if(.not.all_ok) return
        local_ok=.true.
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          if(.not.active(p)) then
            auxiliary_rate%fields(p)%values=0;cycle
          end if
          d=tree%nodes(p)%shape
          do z=1,d(3)
          do y=1,d(2)
          do x=1,d(1)
            if(geometry(p)%volume_fraction(x,y,z)<=0) auxiliary_rate%fields(p)%values(:,x,y,z)=0
          end do
          end do
          end do
          local_ok=local_ok.and.all(ieee_is_finite(auxiliary_rate%fields(p)%values))
        end do
        call agree(local_ok,comm,all_ok)
        if(.not.all_ok) return
      end if
      if(selected_depth>=0.and.redistribution_method/='none') exchange=transport_rate
      if(diffusive) then
        call compute_sdc_mpi_eb_tree_transport_3d(species,transport,tree,owners,geometry,dt, &
          transport_flags(1),transport_flags(2),transport_flags(3),transport_flags(4),gas,temp, &
          transport_rate,flux,stage_theta,comm,all_ok,inflow_mach,inflow_temperature,inflow_pressure, &
          boundary_faces,inflow_state,wall_temperature, &
          target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions, &
          embedded_wall_temperature=embedded_wall_temperature,extended_geometry=heat_geometry, &
          embedded_wall_energy=heat_stage,embedded_frequency=heat_frequency, &
          embedded_noslip=embedded_noslip,embedded_wall_impulse=wall_impulse)
        if(wall_active.and.all_ok) then
          if(stage==0) then
            heat_old=heat_stage
            if(noslip) impulse_old=wall_impulse
          else
            heat_new=heat_stage
            if(noslip) impulse_new=wall_impulse
          end if
        end if
        if(.not.all_ok) return
      else
        stage_theta=1
      end if
      if(has_diffusion) then
        extra=transport_rate;extra_theta=1
        call evaluate_diffusion(gas,temp,dt,selected_depth,extra,extra_flux,extra_theta,comm,local_ok, &
          ancestor_start,ancestor_end,temperature_start,temperature_end,fractions)
        call agree(local_ok,comm,all_ok)
        if(.not.all_ok) return
        call validate_mpi_amr_tree_local_3d(tree,owners,extra,nv,comm,all_ok)
        if(.not.all_ok) return
        local_ok=allocated(extra_flux).and.ieee_is_finite(extra_theta)
        if(local_ok) local_ok=size(extra_flux)==size(owners).and.extra_theta>=0.and.extra_theta<=1
        call agree(local_ok,comm,all_ok)
        if(.not.all_ok) return
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          local_ok=valid_eb_integrated_flux_3d(extra_flux(p),geometry(p),nv)
          local_ok=local_ok.and.all(ieee_is_finite(extra%fields(p)%values))
          if(.not.local_ok) exit
          transport_rate%fields(p)%values=transport_rate%fields(p)%values+extra%fields(p)%values
          if(diffusive) then
            flux(p)%x=flux(p)%x+extra_flux(p)%x
            flux(p)%y=flux(p)%y+extra_flux(p)%y
            flux(p)%z=flux(p)%z+extra_flux(p)%z
          end if
        end do
        call agree(local_ok,comm,all_ok)
        if(.not.all_ok) return
        call MPI_Allreduce(extra_theta,global_theta,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
        if(status/=MPI_SUCCESS) return
        stage_theta=min(stage_theta,global_theta)
        if(.not.diffusive) call move_alloc(extra_flux,flux)
      end if
      if(mol) then
        call compute_mol_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,mol_geometry,gas,temp,dt,2, &
          extra,extra_flux,comm,all_ok,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions)
        if(.not.all_ok) return
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          transport_rate%fields(p)%values=transport_rate%fields(p)%values+extra%fields(p)%values
          if(diffusive.or.has_diffusion) then
            flux(p)%x=flux(p)%x+extra_flux(p)%x
            flux(p)%y=flux(p)%y+extra_flux(p)%y
            flux(p)%z=flux(p)%z+extra_flux(p)%z
          end if
        end do
        if(.not.(diffusive.or.has_diffusion)) call move_alloc(extra_flux,flux)
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          extra%fields(p)%values=0
          if(.not.active(p)) cycle
          call reactive_eb_slip_wall_source_3d(species,gas%fields(p)%values,temp%fields(p)%values(1,:,:,:), &
            geometry(p),extra%fields(p)%values,local_ok)
          if(.not.local_ok) exit
        end do
        call agree(local_ok,comm,all_ok)
        if(.not.all_ok) return
        if(stage==0) then
          call integrate_mpi_eb_tree_3d(tree,owners,geometry,extra,mol_wall_old,comm,all_ok)
          mol_wall_old=dt*mol_wall_old
        else
          call integrate_mpi_eb_tree_3d(tree,owners,geometry,extra,mol_wall_new,comm,all_ok)
          mol_wall_new=dt*mol_wall_new
        end if
        if(.not.all_ok) return
      end if
      if(diffusive.or.has_diffusion.or.mol) then
        if(redistribution_method=='state_redist') then
          if(selected_depth>=0) then
            call redistribute_sdc_mpi_eb_tree_state_level_rate_3d(tree,owners,geometry,gas,transport_rate,nv,selected_depth,dt, &
              redistributed,exchange,comm,all_ok,redist_target,boundary_faces)
          else
            call redistribute_sdc_mpi_eb_tree_state_rate_3d(tree,owners,geometry,gas,transport_rate,nv,dt,redistributed, &
              comm,all_ok,redist_target,boundary_faces)
          end if
          if(.not.all_ok) return
          transport_rate=redistributed
        else if(redistribution_method=='flux_redist') then
          if(selected_depth>=0) then
            call redistribute_sdc_mpi_eb_tree_level_rate_3d(tree,owners,geometry,transport_rate,nv,selected_depth,dt, &
              redistributed,exchange,comm,all_ok,boundary_faces)
          else
            call redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree,owners,geometry,transport_rate,nv,redistributed, &
              comm,all_ok,boundary_faces)
          end if
          if(.not.all_ok) return
          transport_rate=redistributed
        end if
        theta=min(theta,stage_theta)
      end if
      good=.true.
    end subroutine
    subroutine recover(good)
      logical, intent(out) :: good
      logical :: local_ok,cell_ok
      local_ok=.true.
      do node=1,size(owners)
        if(owners(node)/=rank.or..not.active(node)) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          call reactive_conserved_to_primitive(species,candidate%fields(node)%values(:,i,j,k), &
            ct%fields(node)%values(1,i,j,k),primitive,recovered,speed,cell_ok)
          if(.not.cell_ok) then
            local_ok=.false.;cycle
          end if
          ct%fields(node)%values(1,i,j,k)=recovered
        end do
        end do
        end do
      end do
      call agree(local_ok,comm,good)
    end subroutine
  end subroutine advance_mpi_eb_tree_sdc_3d
  subroutine advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,geometry,state,temperature, &
      dt,iterations,chemistry, &
      transport_flags,force,rtol,atol,extended_geometry,boundary_faces, &
      boundary_change,wall_change,chemical_change,source_change,minimum_theta,comm,ok,reaction_source, &
      difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature,level_steps,redistribution, &
      target_volume_fraction,auxiliary_source,auxiliary_change,route_auxiliary,gas_acceleration,halo_source, &
      predictor_builder,post_level,evaluate_diffusion, &
      evaluate_reaction,embedded_wall_temperature,embedded_noslip,mol_mode)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(mpi_eb_sdc_diffusion_source), optional :: evaluate_diffusion
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: dt,force(3),rtol,atol
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry,transport_flags(4)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(MPI_Comm), intent(in) :: comm
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),minimum_theta
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: wall_temperature(6),gas_acceleration(3),embedded_wall_temperature
    logical, optional, intent(in) :: embedded_noslip,mol_mode
    integer(int64), optional, intent(inout) :: level_steps(:)
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(in) :: target_volume_fraction
    procedure(mpi_eb_sdc_recursive_source), optional :: auxiliary_source
    procedure(mpi_eb_sdc_recursive_halo_source), optional :: halo_source
    procedure(mpi_eb_sdc_recursive_predictor_builder), optional :: predictor_builder
    procedure(mpi_eb_sdc_post_level), optional :: post_level
    logical :: has_auxiliary
    real(dp), optional, intent(inout) :: auxiliary_change(:)
    real(dp) :: auxiliary(size(source_change)),source_time
    logical, optional, intent(in) :: route_auxiliary
    logical :: routing
    type(amr_tree_state_3d) :: off_old,off_new
    type(amr_tree_state_3d), allocatable :: auxiliary_pending(:)
    real(dp) :: redist_target
    character(len=16) :: method
    integer :: method_code,method_low,method_high
    type(amr_tree_state_3d), allocatable :: pending(:)
    type(amr_tree_state_3d) :: candidate,ct,starts,ends,ts,te
    type(amr_tree_state_3d), allocatable :: ir
    type(amr_tree_state_3d) :: source0,source1
    type(eb_integrated_flux_3d), allocatable :: total(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: clock_start(:),clock_dt(:)
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: external(size(source_change)),theta
    integer :: nv,node,n(3),levels,c,rank,status,ios
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    has_auxiliary=present(auxiliary_source).or.present(halo_source)
    routing=.false.
    if(present(route_auxiliary)) routing=route_auxiliary
    valid=ieee_is_finite(dt).and..not.(present(auxiliary_source).and.present(halo_source))
    if(valid) valid=dt>0
    valid=valid.and.all([size(boundary_change),size(wall_change),size(chemical_change),size(source_change)]==nv)
    if(present(reaction_source).and..not.chemistry) valid=.false.
    valid=valid.and.(.not.routing.or.has_auxiliary)
    valid=valid.and.(has_auxiliary.eqv.present(auxiliary_change))
    if(present(auxiliary_change)) valid=valid.and.size(auxiliary_change)==nv
    valid=valid.and.(chemistry.or..not.present(evaluate_reaction))
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(merge(1,0,has_auxiliary)+merge(2,0,present(halo_source)),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,has_auxiliary)+merge(2,0,present(halo_source)),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call MPI_Allreduce(merge(1,0,present(predictor_builder)),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,present(predictor_builder)),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call MPI_Allreduce(merge(1,0,present(post_level)),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,present(post_level)),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call MPI_Allreduce(merge(1,0,routing),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,routing),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call MPI_Allreduce(merge(1,0,present(evaluate_reaction)),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,present(evaluate_reaction)),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    method='none'
    if(present(redistribution)) method=redistribution
    call agree(method=='none'.or.method=='flux_redist'.or.method=='state_redist',comm,ready)
    if(.not.ready) return
    method_code=merge(1,0,method=='flux_redist')+merge(2,0,method=='state_redist')
    call MPI_Allreduce(method_code,method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(method_code,method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    redist_target=.5_dp
    if(present(target_volume_fraction)) redist_target=target_volume_fraction
    candidate=state;ct=temperature
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,ready)
    if(.not.ready) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,comm,ready)
    if(.not.ready) return
    levels=size(ratios)
    if(present(level_steps)) then
      call agree(size(level_steps)==levels,comm,ready)
    else
      call agree(.true.,comm,ready)
    end if
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    starts=candidate;ends=candidate;ts=ct;te=ct
    allocate(clock_start(levels),clock_dt(levels),counts(levels),source0%fields(size(tree%nodes)),stat=ios)
    call agree(ios==0,comm,ready)
    if(.not.ready) return
    counts=0;clock_start=0;clock_dt=dt
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(source0%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp,stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
    end do
    call agree(valid,comm,ready)
    if(.not.ready) return
    source1=source0
    if(routing) then
      allocate(auxiliary_pending(0:levels-1))
      do c=0,levels-1
        auxiliary_pending(c)=source0
      end do
    end if
    if(method/='none') then
      allocate(pending(0:levels-1),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      do c=0,levels-1
        pending(c)=source0
      end do
    end if
    if(chemistry) then
      allocate(ir)
      if(present(reaction_source)) then
        ir=reaction_source
      else
        allocate(ir%fields(size(tree%nodes)))
        do node=1,size(tree%nodes)
          if(owners(node)/=rank) cycle
          n=tree%nodes(node)%shape
          allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)),source=0._dp,stat=ios)
          if(ios/=0) then
            valid=.false.;exit
          end if
        end do
      end if
    end if
    call agree(valid,comm,ready)
    if(.not.ready) return
    auxiliary=0;source_time=0
    wall=0;chemical=0;external=0;theta=1
    call advance_level(0,0._dp,dt,1,1,total,valid)
    if(.not.valid) return
    boundary=0
    if(rank==owners(1)) then
    n=tree%nodes(1)%shape
    associate(g=>geometry(1),f=>total(1))
    do c=1,nv
      boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
        g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
        g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
        g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
        g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
        g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
    end do
    end associate
    end if
    call MPI_Bcast(boundary,nv,MPI_DOUBLE_PRECISION,owners(1),comm,status)
    call agree(status==MPI_SUCCESS.and.all(ieee_is_finite(boundary)).and.all(ieee_is_finite(wall)).and. &
      all(ieee_is_finite(chemical)).and.all(ieee_is_finite(external)),comm,ready)
    if(.not.ready) return
    call move_alloc(candidate%fields,state%fields);call move_alloc(ct%fields,temperature%fields)
    boundary_change=boundary;wall_change=wall;chemical_change=chemical;source_change=external
    if(present(auxiliary_change)) auxiliary_change=auxiliary
    minimum_theta=theta
    if(present(reaction_source)) reaction_source=ir
    if(present(level_steps)) level_steps=counts
    ok=.true.
  contains
    subroutine predictor_adapter(rate,interval,iteration,depth,halo,comm,ok)
      type(amr_tree_state_3d), intent(in) :: rate
      real(dp), intent(in) :: interval
      integer, intent(in) :: iteration,depth
      type(amr_tree_state_3d), intent(inout) :: halo
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
      real(dp) :: first(levels-1),last(levels-1)
      first=0;last=0
      first(:depth)=max(0._dp,min(1._dp,(source_time-clock_start(:depth))/clock_dt(:depth)))
      last(:depth)=max(first(:depth),min(1._dp,(source_time+interval-clock_start(:depth))/clock_dt(:depth)))
      call predictor_builder(rate,interval,iteration,depth,source_time,halo,comm,ok, &
        source0,source1,.5_dp*(first+last))
    end subroutine predictor_adapter
    subroutine source_adapter(state,temperature,interval,stage,depth,rate,comm,ok)
      type(amr_tree_state_3d), intent(in) :: state,temperature
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
      type(amr_tree_state_3d) :: off
      integer :: p,x,y,z,d(3)
      logical :: valid,agreed
      real(dp) :: fraction(levels-1),sample_time
      if(present(halo_source)) then
        fraction=0;sample_time=source_time
        if(stage>0) sample_time=sample_time+interval
        fraction(:depth)=max(0._dp,min(1._dp,(sample_time-clock_start(:depth))/clock_dt(:depth)))
        call halo_source(state,temperature,interval,stage,depth,source_time,rate,comm,ok, &
          starts,ends,ts,te,fraction)
      else
        call auxiliary_source(state,temperature,interval,stage,depth,source_time,rate,comm,ok)
      end if
      if(.not.routing) return
      call agree(ok,comm,agreed)
      ok=.false.
      if(.not.agreed) return
      call validate_mpi_amr_tree_local_3d(tree,owners,rate,nv,comm,agreed)
      if(.not.agreed) return
      off=source0;valid=.true.
      do p=1,size(tree%nodes)
        if(allocated(off%fields(p)%values)) off%fields(p)%values=0
      end do
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        if(tree%nodes(p)%level==depth) cycle
        d=tree%nodes(p)%shape
        off%fields(p)%values=0
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0.or.covered(tree,p,[x,y,z])) cycle
          off%fields(p)%values(:,x,y,z)=rate%fields(p)%values(:,x,y,z)
        end do
        end do
        end do
        valid=valid.and.all(ieee_is_finite(off%fields(p)%values))
      end do
      call agree(valid,comm,agreed)
      if(.not.agreed) return
      if(stage==0) then
        off_old=off
      else
        off_new=off
      end if
      ok=.true.

    end subroutine source_adapter
    recursive subroutine advance_level(level,time,interval,iteration,ncycle,flux,good)
      integer, intent(in) :: level,iteration,ncycle
      real(dp), intent(in) :: time,interval
      type(eb_integrated_flux_3d), allocatable, intent(out) :: flux(:)
      logical, intent(out) :: good
      type(eb_integrated_flux_3d), allocatable :: fine_flux(:)
      type(amr_tree_state_3d) :: op,np,correction,redistributed,before_sync,before_t
      type(amr_tree_state_3d) :: off_transfer
      real(dp) :: off_inventory(nv)
      real(dp), allocatable :: step_aux(:)
      type(amr_tree_state_3d), allocatable :: transfer
      real(dp) :: fstart(levels-1),fend(levels-1),b(nv),w(nv),ch(nv),ext(nv),limiter,elapsed,sub_dt
      integer :: p,d(3),step,allocation_status,target,bucket
      logical :: accepted,local_ok
      good=.false.
      if(counts(level+1)==huge(1_int64)) return
      clock_start(level+1)=time;clock_dt(level+1)=interval
      fstart=0;fend=0
      fstart(:level)=max(0._dp,min(1._dp,(time-clock_start(:level))/clock_dt(:level)))
      fend(:level)=max(fstart(:level),min(1._dp,(time+interval-clock_start(:level))/clock_dt(:level)))
      do p=1,size(tree%nodes)
        if(owners(p)/=rank.or.tree%nodes(p)%level/=level) cycle
        starts%fields(p)%values=candidate%fields(p)%values
        ts%fields(p)%values=ct%fields(p)%values
      end do
      if(method/='none') then
        allocate(transfer,stat=allocation_status)
        call agree(allocation_status==0,comm,accepted)
        if(.not.accepted) return
      end if
      if(has_auxiliary) allocate(step_aux(nv))
      source_time=time
      if(present(predictor_builder)) then
      if(has_auxiliary) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,method,transfer, &
        redist_target,source_adapter,step_aux,gas_acceleration=gas_acceleration,predictor_builder=predictor_adapter, &
        evaluate_diffusion=evaluate_diffusion, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode)
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,method,transfer, &
        redist_target,gas_acceleration=gas_acceleration,predictor_builder=predictor_adapter, &
        evaluate_diffusion=evaluate_diffusion, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode)
      end if
      else
      if(has_auxiliary) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,method,transfer, &
        redist_target,source_adapter,step_aux,gas_acceleration=gas_acceleration,evaluate_diffusion=evaluate_diffusion, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode)
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,method,transfer, &
        redist_target,gas_acceleration=gas_acceleration,evaluate_diffusion=evaluate_diffusion, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode)
      end if
      end if
      if(.not.accepted) return
      if(method/='none') then
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          target=tree%nodes(p)%level;bucket=min(level,target)
          pending(bucket)%fields(p)%values=pending(bucket)%fields(p)%values+transfer%fields(p)%values
        end do
      end if
      if(routing) then
        off_transfer=off_old
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          target=tree%nodes(p)%level;bucket=min(level,target)
          off_transfer%fields(p)%values=.5_dp*interval*(off_old%fields(p)%values+off_new%fields(p)%values)
          auxiliary_pending(bucket)%fields(p)%values=auxiliary_pending(bucket)%fields(p)%values+ &
            off_transfer%fields(p)%values
        end do
        call integrate_mpi_eb_tree_3d(tree,owners,geometry,off_transfer,off_inventory,comm,accepted)
        if(.not.accepted) return
        step_aux=step_aux+off_inventory
      end if
      if(has_auxiliary) then
        auxiliary=auxiliary+step_aux
        call agree(all(ieee_is_finite(auxiliary)),comm,accepted)
        if(.not.accepted) return
      end if
      counts(level+1)=counts(level+1)+1
      wall=wall+w;chemical=chemical+ch;external=external+ext;theta=min(theta,limiter)
      local_ok=.true.
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        if(tree%nodes(p)%level==level) then
          ends%fields(p)%values=candidate%fields(p)%values
          te%fields(p)%values=ct%fields(p)%values
          source0%fields(p)%values=op%fields(p)%values;source1%fields(p)%values=np%fields(p)%values
        else
          d=tree%nodes(p)%shape
          allocate(flux(p)%x(nv,0:d(1),d(2),d(3)),flux(p)%y(nv,d(1),0:d(2),d(3)), &
            flux(p)%z(nv,d(1),d(2),0:d(3)),stat=allocation_status)
          if(allocation_status/=0) then
            local_ok=.false.;exit
          end if
          flux(p)%x=0;flux(p)%y=0;flux(p)%z=0
        end if
      end do
      call agree(local_ok,comm,accepted)
      if(.not.accepted) return
      if(level+1<levels) then
        elapsed=0
        do step=1,ratios(level+2)
          sub_dt=interval/real(ratios(level+2),dp)
          if(step==ratios(level+2)) sub_dt=interval-elapsed
          if(sub_dt<=0.or.time+elapsed+sub_dt<=time+elapsed) return
          call advance_level(level+1,time+elapsed,sub_dt,step,ratios(level+2),fine_flux,accepted)
          if(.not.accepted) return
          do p=1,size(tree%nodes)
            if(owners(p)/=rank.or.tree%nodes(p)%level/=level+1) cycle
            flux(p)%x=flux(p)%x+fine_flux(p)%x
            flux(p)%y=flux(p)%y+fine_flux(p)%y
            flux(p)%z=flux(p)%z+fine_flux(p)%z
          end do
          elapsed=elapsed+sub_dt
        end do
      end if
      if(method/='none') then
        before_sync=candidate;before_t=ct
      end if
      call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,flux,comm,accepted, &
        minimum_level=level,maximum_level=min(level+1,levels-1),boundary_faces=boundary_faces)
      if(.not.accepted) return
      if(method/='none') then
        correction=source0
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          correction%fields(p)%values=0
          target=tree%nodes(p)%level
          if(target<level.or.target>level+1) cycle
          correction%fields(p)%values=candidate%fields(p)%values-before_sync%fields(p)%values
        end do
        if(method=='state_redist') then
          call redistribute_mpi_eb_tree_state_3d(tree,owners,geometry,correction,nv,redistributed, &
            comm,accepted,redist_target,boundary_faces)
        else
          call redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree,owners,geometry,correction,nv,redistributed, &
            comm,accepted,boundary_faces)
        end if
        if(.not.accepted) return
        candidate=before_sync;ct=before_t
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          target=tree%nodes(p)%level
          if(target<level) then
            pending(target)%fields(p)%values=pending(target)%fields(p)%values+redistributed%fields(p)%values
          else
            candidate%fields(p)%values=candidate%fields(p)%values+redistributed%fields(p)%values+ &
              pending(level)%fields(p)%values
          end if
          pending(level)%fields(p)%values=0
        end do
        call recover_transfer(accepted)
        if(.not.accepted) return
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,accepted)
        if(.not.accepted) return
      end if
      if(routing) then
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          candidate%fields(p)%values=candidate%fields(p)%values+auxiliary_pending(level)%fields(p)%values
          auxiliary_pending(level)%fields(p)%values=0
        end do
        call recover_transfer(accepted)
        if(.not.accepted) return
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,accepted)
        if(.not.accepted) return
      end if
      if(present(post_level)) then
        call post_level(level,time,interval,iteration,ncycle,comm,local_ok)
        call agree(local_ok,comm,accepted)
        if(.not.accepted) return
      end if
      good=.true.
    end subroutine advance_level
    subroutine recover_transfer(good)
      logical, intent(out) :: good
      real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed
      integer :: p,x,y,z,d(3)
      logical :: accepted,local_ok
      local_ok=.true.
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          call reactive_conserved_to_primitive(species,candidate%fields(p)%values(:,x,y,z), &
            ct%fields(p)%values(1,x,y,z),primitive,recovered,speed,accepted)
          if(.not.accepted) then
            local_ok=.false.;cycle
          end if
          ct%fields(p)%values(1,x,y,z)=recovered
        end do
        end do
        end do
      end do
      call agree(local_ok,comm,good)
    end subroutine recover_transfer
  end subroutine advance_mpi_eb_tree_sdc_subcycle_3d
end module mpi_eb_tree_full_3d_mod
