module mpi_eb_tree_full_3d_mod
  use mpi_f08
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: advance_reactive_chemistry_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d,average_down_mpi_eb_tree_3d
  use mpi_eb_tree_hydro_3d_mod, only: advance_mpi_eb_tree_hydro_3d
  use mpi_eb_tree_transport_3d_mod, only: advance_mpi_eb_tree_transport_3d
  use mpi_eb_tree_subcycle_hydro_3d_mod, only: advance_mpi_eb_tree_subcycle_hydro_3d
  use mpi_eb_tree_subcycle_transport_3d_mod, only: advance_mpi_eb_tree_subcycle_transport_3d
  implicit none
  private
  public :: react_mpi_eb_tree_3d,advance_mpi_eb_tree_full_3d
contains
  subroutine agree(condition,comm,ok)
    logical, intent(in) :: condition
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: status
    call MPI_Allreduce(condition,ok,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if (status/=MPI_SUCCESS) ok=.false.
  end subroutine

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
      state,temperature,boundary_change,wall_change,chemical_change,minimum_theta,comm,ok,integrator,subcycling)
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
    valid=ieee_is_finite(dt)
    if (valid) valid=dt>0._dp
    call agree(valid.and.size(boundary_change)==nv.and.size(wall_change)==nv.and.size(chemical_change)==nv,comm,ready)
    if (.not.ready) return
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
            viscosity,conduction,diffusion,barodiffusion,candidate,ct,exchange,stage_theta,comm,valid)
        else
          call advance_mpi_eb_tree_transport_3d(species,transport,tree,owners,geometry,0.5_dp*dt, &
            viscosity,conduction,diffusion,barodiffusion,candidate,ct,exchange,stage_theta,comm,valid)
        end if
        if (.not.valid) return
        transport_sum=transport_sum+exchange;theta=min(theta,stage_theta)
      end if
      if (half==1) then
        if (recursive_steps) then
          call advance_mpi_eb_tree_subcycle_hydro_3d(species,tree,owners,geometry,solver,dt,candidate,ct, &
            boundary,wall,comm,valid)
        else
          call advance_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,solver,dt,candidate,ct,boundary,wall,comm,valid)
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
end module mpi_eb_tree_full_3d_mod
