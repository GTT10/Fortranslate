module mpi_amr_patch_tree_reflux_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,face_marks,get_face
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use mpi_amr_patch_tree_flux_peers_3d_mod, only: synchronize_mpi_amr_tree_flux_peers_3d,valid_local_tree_flux_level_3d
  use mpi_amr_patch_tree_flux_restriction_3d_mod, only: restrict_mpi_amr_tree_flux_level_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_reflux_budget_3d_mod, only: balance_mpi_amr_reflux_species_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  public :: reflux_mpi_amr_tree_3d
contains
  subroutine reflux_mpi_amr_tree_3d(tree,owners,state,nvar,flux,root_spacing,faces,boundary_correction, &
      maximum_correction,comm,ok,minimum_depth,synchronized_flux,first_species,last_species, &
      auxiliary_state,auxiliary_flux,auxiliary_synchronized_flux)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(amr_tree_state_3d), intent(inout) :: state
    type(amr_tree_flux_3d), intent(in) :: flux(:)
    real(dp), intent(in) :: root_spacing(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(inout) :: boundary_correction(:),maximum_correction
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_depth
    integer, optional, intent(in) :: first_species,last_species
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: synchronized_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: auxiliary_state
    type(amr_tree_flux_3d), optional, intent(in) :: auxiliary_flux(:)
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: auxiliary_synchronized_flux(:)
    ! The reported maximum belongs to the caller's conserved components.
    ! Appending held internal energy for synchronization must not change it.
    call reflux_mpi_amr_tree_core(tree,owners,state,nvar,flux,root_spacing,faces,boundary_correction, &
      maximum_correction,comm,ok,minimum_depth,synchronized_flux,first_species,last_species, &
      auxiliary_state,auxiliary_flux,auxiliary_synchronized_flux,nvar)
  end subroutine reflux_mpi_amr_tree_3d

  recursive subroutine reflux_mpi_amr_tree_core(tree,owners,state,nvar,flux,root_spacing,faces,boundary_correction, &
      maximum_correction,comm,ok,minimum_depth,synchronized_flux,first_species,last_species, &
      auxiliary_state,auxiliary_flux,auxiliary_synchronized_flux,diagnostic_components)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar,diagnostic_components
    type(amr_tree_state_3d), intent(inout) :: state
    type(amr_tree_flux_3d), intent(in) :: flux(:)
    real(dp), intent(in) :: root_spacing(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(inout) :: boundary_correction(:),maximum_correction
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_depth
    integer, optional, intent(in) :: first_species,last_species
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: synchronized_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: auxiliary_state
    type(amr_tree_flux_3d), optional, intent(in) :: auxiliary_flux(:)
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: auxiliary_synchronized_flux(:)
    type(amr_tree_state_3d) :: candidate
    type(amr_tree_state_3d) :: coupled,coupled_auxiliary
    type(amr_tree_flux_3d), allocatable :: coupled_flux(:),coupled_effective(:),gas_effective(:),scalar_effective(:)
    type(amr_tree_flux_3d), allocatable :: work(:)
    type(face_marks), allocatable :: marks(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: spacing(:,:),budget(:),global_budget(:),correction(:),delta_low(:),delta_high(:)
    logical, allocatable :: covered(:,:,:)
    real(dp) :: largest,global_largest,area,low_spacing(3),high_spacing(3)
    integer :: rank,status,ios,node,child,level,n(3),lo(3),hi(3),i,j,k,axis,idx(3),upper_idx(3),floor_depth,low,high
    logical :: valid,ready
    integer :: bounds(2),bounds_low(2),bounds_high(2)
    integer :: coupled_flags,coupled_low,coupled_high
    real(dp), allocatable :: coupled_budget(:)
    ok=.false.
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nvar,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    floor_depth=0
    if (present(minimum_depth)) floor_depth=minimum_depth
    call MPI_Allreduce(floor_depth,low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(floor_depth,high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (low/=high.or.floor_depth<0.or.floor_depth>=tree%level_count()) return
    coupled_flags=0
    if(present(auxiliary_state)) coupled_flags=coupled_flags+1
    if(present(auxiliary_flux)) coupled_flags=coupled_flags+2
    if(present(auxiliary_synchronized_flux)) coupled_flags=coupled_flags+4
    call MPI_Allreduce(coupled_flags,coupled_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(coupled_flags,coupled_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(coupled_low/=coupled_high) return
    if(coupled_flags/=0.and.coupled_flags/=3.and.coupled_flags/=7) return
    bounds=0
    if(present(first_species).neqv.present(last_species)) bounds=-1
    if(present(first_species).and.present(last_species)) bounds=[first_species,last_species]
    call MPI_Allreduce(bounds,bounds_low,2,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,bounds_high,2,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(bounds_low/=bounds_high).or.any(bounds<0)) return
    if(any(bounds/=0)) then
      if(bounds(1)<2.or.bounds(2)<bounds(1).or.bounds(2)>nvar) return
    endif
    valid=size(flux)==size(tree%nodes).and.size(boundary_correction)==nvar.and.valid_reactive_boundaries_3d(faces)
    if (.not.all(ieee_is_finite(root_spacing)).or.any(root_spacing<=0)) valid=.false.
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(state%fields(node)%values))) valid=.false.
      if (any(tree%nodes(node)%shape==huge(1))) valid=.false.
    end do
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(root_spacing,low_spacing,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(root_spacing,high_spacing,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low_spacing/=high_spacing)) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready)
    if (.not.ready) return
    if(present(auxiliary_state)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,auxiliary_state,1,comm,ready)
      if(.not.ready) return
      valid=nvar<huge(1).and.size(auxiliary_flux)==size(flux)
      call agree(valid,ready)
      if(.not.ready) return
      allocate(coupled%fields(size(flux)),coupled_auxiliary%fields(size(flux)),coupled_flux(size(flux)), &
        gas_effective(size(flux)),scalar_effective(size(flux)),coupled_budget(nvar+1),stat=ios)
      call agree(ios==0,ready)
      if(.not.ready) return
      valid=.true.
      do node=1,size(flux)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(coupled%fields(node)%values(nvar+1,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        endif
        coupled%fields(node)%values(:nvar,:,:,:)=state%fields(node)%values
        coupled%fields(node)%values(nvar+1,:,:,:)=auxiliary_state%fields(node)%values(1,:,:,:)
        if(tree%nodes(node)%level<floor_depth) cycle
        if(.not.allocated(flux(node)%x).or..not.allocated(flux(node)%y).or..not.allocated(flux(node)%z)) then
          valid=.false.;exit
        endif
        if(.not.allocated(auxiliary_flux(node)%x).or..not.allocated(auxiliary_flux(node)%y).or. &
            .not.allocated(auxiliary_flux(node)%z)) then
          valid=.false.;exit
        endif
        if(any(shape(flux(node)%x)/=[nvar,n(1)+1,n(2),n(3)]).or. &
          any(shape(flux(node)%y)/=[nvar,n(1),n(2)+1,n(3)]).or. &
          any(shape(flux(node)%z)/=[nvar,n(1),n(2),n(3)+1]).or. &
          any(shape(auxiliary_flux(node)%x)/=[1,n(1)+1,n(2),n(3)]).or. &
          any(shape(auxiliary_flux(node)%y)/=[1,n(1),n(2)+1,n(3)]).or. &
          any(shape(auxiliary_flux(node)%z)/=[1,n(1),n(2),n(3)+1])) then
          valid=.false.;exit
        endif
        allocate(coupled_flux(node)%x(nvar+1,n(1)+1,n(2),n(3)), &
          coupled_flux(node)%y(nvar+1,n(1),n(2)+1,n(3)),coupled_flux(node)%z(nvar+1,n(1),n(2),n(3)+1),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        endif
        coupled_flux(node)%x(:nvar,:,:,:)=flux(node)%x;coupled_flux(node)%x(nvar+1,:,:,:)=auxiliary_flux(node)%x(1,:,:,:)
        coupled_flux(node)%y(:nvar,:,:,:)=flux(node)%y;coupled_flux(node)%y(nvar+1,:,:,:)=auxiliary_flux(node)%y(1,:,:,:)
        coupled_flux(node)%z(:nvar,:,:,:)=flux(node)%z;coupled_flux(node)%z(nvar+1,:,:,:)=auxiliary_flux(node)%z(1,:,:,:)
      enddo
      call agree(valid,ready)
      if(.not.ready) return
      coupled_budget=0;global_largest=0
      call reflux_mpi_amr_tree_core(tree,owners,coupled,nvar+1,coupled_flux,root_spacing,faces,coupled_budget, &
        global_largest,comm,valid,minimum_depth=floor_depth,synchronized_flux=coupled_effective, &
        first_species=first_species,last_species=last_species,diagnostic_components=diagnostic_components)
      if(.not.valid) return
      ! Allocate all output copies before publishing either conserved field.
      candidate=state;coupled_auxiliary=auxiliary_state
      do node=1,size(flux)
        if(owners(node)/=rank) cycle
        candidate%fields(node)%values=coupled%fields(node)%values(:nvar,:,:,:)
        coupled_auxiliary%fields(node)%values(1,:,:,:)=coupled%fields(node)%values(nvar+1,:,:,:)
        if(tree%nodes(node)%level<floor_depth) cycle
        gas_effective(node)%x=coupled_effective(node)%x(:nvar,:,:,:)
        gas_effective(node)%y=coupled_effective(node)%y(:nvar,:,:,:)
        gas_effective(node)%z=coupled_effective(node)%z(:nvar,:,:,:)
        scalar_effective(node)%x=coupled_effective(node)%x(nvar+1:nvar+1,:,:,:)
        scalar_effective(node)%y=coupled_effective(node)%y(nvar+1:nvar+1,:,:,:)
        scalar_effective(node)%z=coupled_effective(node)%z(nvar+1:nvar+1,:,:,:)
      enddo
      call move_alloc(candidate%fields,state%fields)
      call move_alloc(coupled_auxiliary%fields,auxiliary_state%fields)
      if(present(synchronized_flux)) call move_alloc(gas_effective,synchronized_flux)
      if(present(auxiliary_synchronized_flux)) call move_alloc(scalar_effective,auxiliary_synchronized_flux)
      boundary_correction=coupled_budget(:nvar);maximum_correction=global_largest;ok=.true.;return
    endif
    allocate(spacing(3,size(ratios)),budget(nvar),global_budget(nvar),correction(nvar),delta_low(nvar),delta_high(nvar), &
      work(size(flux)),marks(size(flux)),candidate%fields(size(flux)),stat=ios)
    call agree(ios==0,ready)
    if (.not.ready) return
    spacing(:,1)=root_spacing
    do level=2,size(ratios)
      spacing(:,level)=spacing(:,level-1)/real(ratios(level),dp)
      if (any(spacing(:,level)<=0)) return
    end do
    do node=1,size(flux)
      if (owners(node)/=rank) cycle
      allocate(candidate%fields(node)%values,source=state%fields(node)%values,stat=ios)
      if (ios/=0) exit
      if (tree%nodes(node)%level<floor_depth) cycle
      n=tree%nodes(node)%shape
      allocate(marks(node)%x(n(1)+1,n(2),n(3)),marks(node)%y(n(1),n(2)+1,n(3)), &
        marks(node)%z(n(1),n(2),n(3)+1),stat=ios)
      if (ios/=0) exit
      marks(node)%x=.false.;marks(node)%y=.false.;marks(node)%z=.false.
    end do
    call agree(ios==0,ready)
    if (.not.ready) return
    valid=.true.
    do level=floor_depth,tree%level_count()-1
      if (.not.valid_local_tree_flux_level_3d(tree,owners,rank,nvar,level,flux,marks)) valid=.false.
    end do
    call agree(valid,ready)
    if (.not.ready) return
    do node=1,size(flux)
      if (owners(node)/=rank.or.tree%nodes(node)%level<floor_depth) cycle
      allocate(work(node)%x,source=flux(node)%x,stat=ios)
      if (ios/=0) exit
      allocate(work(node)%y,source=flux(node)%y,stat=ios)
      if (ios/=0) exit
      allocate(work(node)%z,source=flux(node)%z,stat=ios)
      if (ios/=0) exit
    end do
    call agree(ios==0,ready)
    if (.not.ready) return
    budget=0;largest=0
    do level=size(ratios),floor_depth+1,-1
      call synchronize_mpi_amr_tree_flux_peers_3d(tree,owners,nvar,level-1,faces,work,marks,comm,valid)
      if (.not.valid) return
      valid=.true.
      do node=1,size(flux)
        if (owners(node)/=rank.or.tree%nodes(node)%level+1/=level) cycle
        n=tree%nodes(node)%shape
        allocate(covered(n(1),n(2),n(3)),stat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
        covered=.false.
        do child=node+1,size(flux)
          if (tree%nodes(child)%parent/=node) cycle
          lo=[tree%nodes(child)%patch%coarse_i_lower,tree%nodes(child)%patch%coarse_j_lower, &
            tree%nodes(child)%patch%coarse_k_lower]
          hi=[tree%nodes(child)%patch%coarse_i_upper,tree%nodes(child)%patch%coarse_j_upper, &
            tree%nodes(child)%patch%coarse_k_upper]
          covered(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
        end do
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              if (covered(i,j,k)) cycle
              idx=[i,j,k];correction=0
              do axis=1,3
                upper_idx=idx;upper_idx(axis)=upper_idx(axis)+1
                delta_low=get_face(work(node),axis,idx)-get_face(flux(node),axis,idx)
                delta_high=get_face(work(node),axis,upper_idx)-get_face(flux(node),axis,upper_idx)
                correction=correction+(delta_low-delta_high)/spacing(axis,level)
              end do
              candidate%fields(node)%values(:,i,j,k)=candidate%fields(node)%values(:,i,j,k)+correction
              largest=max(largest,maxval(abs(correction(:diagnostic_components))))
            end do
          end do
        end do
        deallocate(covered)
        if (.not.all(ieee_is_finite(candidate%fields(node)%values))) valid=.false.
      end do
      call agree(valid,ready)
      if (.not.ready) return
      if (level>floor_depth+1) then
        call restrict_mpi_amr_tree_flux_level_3d(tree,owners,nvar,level-1,work,marks,comm,valid)
        if (.not.valid) return
      end if
    end do
    if(bounds(1)>0) then
      call balance_mpi_amr_reflux_species_3d(tree,owners,state,candidate,nvar,work,root_spacing, &
        faces,bounds(1),bounds(2),floor_depth,comm,valid)
      if(.not.valid) return
      ! Preserve the directly accumulated flux correction for unchanged
      ! components; subtracting large total energies loses low-order bits.
      ! Only the repaired species can acquire an additional correction here.
      do node=1,size(flux)
        if(owners(node)/=rank.or.tree%nodes(node)%level<floor_depth) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              largest=max(largest,maxval(abs(candidate%fields(node)%values(bounds(1):bounds(2),i,j,k)- &
                state%fields(node)%values(bounds(1):bounds(2),i,j,k))))
            enddo
          enddo
        enddo
      enddo
    endif
    if (floor_depth==0.and.rank==owners(1)) then
      n=tree%nodes(1)%shape
      do axis=1,3
        if (faces(2*axis-1)=='periodic') cycle
        area=product(root_spacing)/root_spacing(axis)
        hi=n;hi(axis)=1
        do k=1,hi(3)
          do j=1,hi(2)
            do i=1,hi(1)
              idx=[i,j,k];upper_idx=idx;upper_idx(axis)=n(axis)+1
              budget=budget+area*(get_face(work(1),axis,idx)-get_face(flux(1),axis,idx)- &
                get_face(work(1),axis,upper_idx)+get_face(flux(1),axis,upper_idx))
            end do
          end do
        end do
      end do
    end if
    call agree(all(ieee_is_finite(budget)).and.ieee_is_finite(largest),ready)
    if (.not.ready) return
    call MPI_Allreduce(budget,global_budget,nvar,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(largest,global_largest,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    call average_down_mpi_amr_tree_3d(tree,owners,candidate,nvar,comm,valid,minimum_depth=floor_depth)
    if (.not.valid) return
    if (present(synchronized_flux)) call move_alloc(work,synchronized_flux)
    call move_alloc(candidate%fields,state%fields)
    boundary_correction=global_budget;maximum_correction=global_largest;ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine
end module
