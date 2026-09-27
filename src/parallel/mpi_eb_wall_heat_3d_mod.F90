module mpi_eb_wall_heat_3d_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use state_indices_mod, only: iet
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_two_level_3d_mod, only: matching_eb_metrics_3d
  use eb_wall_heat_source_3d_mod, only: compute_eb_wall_transport_source_3d
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d,integrate_mpi_eb_tree_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d
  implicit none
  private
  public :: compute_mpi_eb_wall_heat_3d,compute_mpi_eb_wall_transport_3d
contains
  ! Two same-resolution donor layers cover every quadratic stencil. Geometry
  ! is replicated; temperature/state are exchanged from owners at each call.
  subroutine compute_mpi_eb_wall_transport_3d(species,transport,tree,owners,geometry,extended_geometry, &
      state,temperature,conduction,noslip,wall_temperature,rate,wall_exchange,frequency,comm,ok, &
      boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,target_depth, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    logical, intent(in) :: conduction,noslip
    real(dp), intent(in) :: wall_temperature
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(inout) :: wall_exchange(:),frequency
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(amr_tree_state_3d) :: metrics,ghost,gt,vf,candidate
    type(eb_geometry_3d) :: center
    real(dp), allocatable :: source(:,:,:,:)
    real(dp) :: controls(3),low(3),high(3),local_frequency,global_frequency,node_frequency
    real(dp) :: inventory(reactive_nvar(size(species))),node_exchange(reactive_nvar(size(species)))
    integer :: node,n(3),rank,status,ios,nv,depth,lo_depth,hi_depth
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species));depth=-1
    if(present(target_depth)) depth=target_depth
    valid=ieee_is_finite(wall_temperature).and.size(transport)==size(species).and.size(wall_exchange)==nv
    if(conduction) valid=valid.and.wall_temperature>0
    controls=[merge(1._dp,0._dp,conduction),merge(1._dp,0._dp,noslip),wall_temperature]
    valid=valid.and.depth>=-1
    if(present(target_depth)) valid=valid.and.depth>=0
    call agree(valid,ready)
    if(.not.ready) return
    call MPI_Allreduce(controls,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    call MPI_Allreduce(depth,lo_depth,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(depth,hi_depth,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.lo_depth/=hi_depth) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,ready,padding=2)
    if(.not.ready) return
    call agree(depth<tree%level_count(),ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    allocate(metrics%fields(size(owners)),candidate%fields(size(owners)),stat=ios)
    call agree(ios==0,ready)
    if(.not.ready) return
    do node=1,size(owners)
      n=tree%nodes(node)%shape
      call extract_eb_geometry_3d(extended_geometry(node),[3,3,3],n+2,center,valid)
      if(valid) valid=matching_eb_metrics_3d(center,geometry(node))
      if(.not.valid) exit
      if(owners(node)/=rank) cycle
      allocate(metrics%fields(node)%values(1,-1:n(1)+2,-1:n(2)+2,-1:n(3)+2), &
        candidate%fields(node)%values(nv,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      metrics%fields(node)%values(1,:,:,:)=extended_geometry(node)%volume_fraction
      candidate%fields(node)%values=0
    end do
    call agree(valid,ready)
    if(.not.ready) return
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,ghost,gt,comm,ready, &
      halo_width=2,target_volume_fraction=metrics,ghost_volume_fraction=vf, &
      restrict_snapshot=.false.,boundary_faces=boundary_faces, &
      inflow_state=inflow_state,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
      temperature_end=temperature_end,ancestor_fraction=ancestor_fraction)
    if(.not.ready) return
    local_frequency=0
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
      n=tree%nodes(node)%shape
      allocate(source(nv,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      call compute_eb_wall_transport_source_3d(species,transport,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:), &
        extended_geometry(node),[3,3,3],n+2,conduction,noslip,wall_temperature,source,node_exchange,node_frequency,valid)
      if(.not.valid) exit
      candidate%fields(node)%values=source
      local_frequency=max(local_frequency,node_frequency)
      deallocate(source)
    end do
    call agree(valid,ready)
    if(.not.ready) return
    call MPI_Allreduce(local_frequency,global_frequency,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    ! Composite integration excludes coarse cells covered by children.
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,candidate,inventory,comm,ready)
    if(.not.ready) return
    call move_alloc(candidate%fields,rate%fields)
    wall_exchange=inventory;frequency=global_frequency;ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine
  ! Preserve the heat-only API and its scalar ledger.
  subroutine compute_mpi_eb_wall_heat_3d(species,transport,tree,owners,geometry,extended_geometry, &
      state,temperature,wall_temperature,rate,wall_power,frequency,comm,ok, &
      boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,target_depth, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: wall_temperature
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(inout) :: wall_power,frequency
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    real(dp) :: exchange(reactive_nvar(size(species)))
    call compute_mpi_eb_wall_transport_3d(species,transport,tree,owners,geometry,extended_geometry,state,temperature, &
      .true.,.false.,wall_temperature,rate,exchange,frequency,comm,ok,boundary_faces,inflow_state,inflow_mach, &
      inflow_temperature,inflow_pressure,target_depth,ancestor_start,ancestor_end,temperature_start, &
      temperature_end,ancestor_fraction)
    if(ok) wall_power=exchange(iet)
  end subroutine
end module
