module mpi_amr_patch_tree_chemistry_3d_mod
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: advance_reactive_chemistry_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use mpi_amr_patch_tree_hydro_3d_mod, only: advance_mpi_amr_tree_hydro_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  implicit none
  private
  public :: advance_mpi_amr_tree_chemistry_3d,advance_mpi_amr_tree_strang_3d
contains
  subroutine advance_mpi_amr_tree_chemistry_3d(tree,owners,species,reactions,state,temperature,interval,rtol,atol, &
      comm,ok,chemistry_integrator)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: interval,rtol,atol
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: chemistry_integrator
    type(amr_tree_state_3d) :: candidate,temps
    logical, allocatable :: active(:,:,:)
    real(dp) :: controls(3),low(3),high(3)
    character(len=64) :: method,reference_method
    integer :: rank,status,ios,counts(2),min_counts(2),max_counts(2)
    integer :: node,child,n(3),lo(3),hi(3)
    logical :: valid
    ok=.false.
    controls=[interval,rtol,atol]
    valid=all(ieee_is_finite(controls))
    if (valid) valid=interval>=0.and.min(rtol,atol)>0
    method='(default)'
    if (present(chemistry_integrator)) then
      if (len_trim(chemistry_integrator)>64) valid=.false.
      method=chemistry_integrator
    end if
    call agree(comm,valid)
    if (.not.valid) return
    call MPI_Allreduce(controls,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    counts=[size(species),size(reactions)]
    call MPI_Allreduce(counts,min_counts,2,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(counts,max_counts,2,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(min_counts/=max_counts)) return
    reference_method=method
    call MPI_Bcast(reference_method,64,MPI_CHARACTER,0,comm,status)
    if (status/=MPI_SUCCESS) return
    valid=method==reference_method;call agree(comm,valid)
    if (.not.valid) return
    call clone_inputs(tree,owners,size(species),state,temperature,candidate,temps,comm,rank,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(active(n(1),n(2),n(3)),stat=ios)
      if (ios/=0) then
        valid=.false.;exit
      end if
      active=.true.
      ! Only exposed cells own chemistry. Covered parents are supplied by
      ! restriction after every leaf and uncovered ancestor has succeeded.
      do child=node+1,size(tree%nodes)
        if (tree%nodes(child)%parent/=node) cycle
        lo=[tree%nodes(child)%patch%coarse_i_lower,tree%nodes(child)%patch%coarse_j_lower, &
          tree%nodes(child)%patch%coarse_k_lower]
        hi=[tree%nodes(child)%patch%coarse_i_upper,tree%nodes(child)%patch%coarse_j_upper, &
          tree%nodes(child)%patch%coarse_k_upper]
        active(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      call advance_reactive_chemistry_3d(species,reactions,candidate%fields(node)%values, &
        temps%fields(node)%values(1,:,:,:),n(1),n(2),n(3),interval,rtol,atol,valid,active,chemistry_integrator)
      if (.not.valid) exit
      deallocate(active)
    end do
    call agree(comm,valid)
    if (.not.valid) return
    call average_down_mpi_amr_tree_3d(tree,owners,candidate,reactive_nvar(size(species)),comm,valid)
    if (.not.valid) return
    call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,candidate,temps,comm,valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    ok=.true.
  end subroutine advance_mpi_amr_tree_chemistry_3d

  subroutine advance_mpi_amr_tree_strang_3d(tree,owners,species,reactions,state,temperature,root_spacing,dt,faces, &
      reconstruction,limiter,solver,rtol,atol,boundary_change,maximum_correction,comm,ok, &
      chemistry_integrator,inflow_state,inflow_temperature)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),dt,rtol,atol
    character(len=*), intent(in) :: faces(6),reconstruction,limiter,solver
    real(dp), intent(inout) :: boundary_change(:),maximum_correction
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: chemistry_integrator
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:)
    type(amr_tree_state_3d) :: candidate,temps
    real(dp) :: budget(size(boundary_change)),largest
    integer :: rank
    logical :: valid
    ok=.false.
    valid=ieee_is_finite(dt).and.dt>0
    call agree(comm,valid)
    if (.not.valid) return
    call clone_inputs(tree,owners,size(species),state,temperature,candidate,temps,comm,rank,valid)
    if (.not.valid) return
    call advance_mpi_amr_tree_chemistry_3d(tree,owners,species,reactions,candidate,temps,.5_dp*dt,rtol,atol, &
      comm,valid,chemistry_integrator)
    if (.not.valid) return
    budget=0;largest=0
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,candidate,temps,root_spacing,dt,faces, &
      reconstruction,limiter,solver,budget,largest,comm,valid,inflow_state,inflow_temperature)
    if (.not.valid) return
    call advance_mpi_amr_tree_chemistry_3d(tree,owners,species,reactions,candidate,temps,.5_dp*dt,rtol,atol, &
      comm,valid,chemistry_integrator)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    boundary_change=budget;maximum_correction=largest;ok=.true.
  end subroutine advance_mpi_amr_tree_strang_3d
  subroutine clone_inputs(tree,owners,nspecies,state,temperature,candidate,temps,comm,rank,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nspecies
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(amr_tree_state_3d), intent(out) :: candidate,temps
    type(MPI_Comm), intent(in) :: comm
    integer, intent(out) :: rank
    logical, intent(out) :: ok
    integer :: node,status,ios
    call validate_mpi_amr_tree_local_3d(tree,owners,state,reactive_nvar(nspecies),comm,ok)
    if (.not.ok) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ok)
    if (.not.ok) return
    ok=.false.
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    ok=.true.
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(state%fields(node)%values))) ok=.false.
      if (.not.all(ieee_is_finite(temperature%fields(node)%values))) ok=.false.
      if (any(temperature%fields(node)%values<=0)) ok=.false.
    end do
    call agree(comm,ok)
    if (.not.ok) return
    allocate(candidate%fields(size(tree%nodes)),temps%fields(size(tree%nodes)),stat=ios)
    if (ios==0) then
      do node=1,size(tree%nodes)
        if (owners(node)/=rank) cycle
        allocate(candidate%fields(node)%values,source=state%fields(node)%values,stat=ios)
        if (ios/=0) exit
        allocate(temps%fields(node)%values,source=temperature%fields(node)%values,stat=ios)
        if (ios/=0) exit
      end do
    end if
    ok=ios==0;call agree(comm,ok)
  end subroutine
  subroutine agree(comm,condition)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(inout) :: condition
    logical :: result
    integer :: ierr
    call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    condition=.false.
    if (ierr==MPI_SUCCESS) condition=result
  end subroutine
end module mpi_amr_patch_tree_chemistry_3d_mod
