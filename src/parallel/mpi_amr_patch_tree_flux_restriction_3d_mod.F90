module mpi_amr_patch_tree_flux_restriction_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,face_marks,get_face,set_face,marked,set_mark,matches
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_metadata_3d
  use mpi_amr_patch_tree_flux_peers_3d_mod, only: valid_local_tree_flux_level_3d
  implicit none
  private
  public :: restrict_mpi_amr_tree_flux_level_3d
contains
  ! Restrict one child depth into its parents. Peer synchronization is separate.
  subroutine restrict_mpi_amr_tree_flux_level_3d(tree,owners,nvar,depth,flux,refined,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar,depth
    type(amr_tree_flux_3d), intent(inout) :: flux(:)
    type(face_marks), intent(inout) :: refined(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_flux_3d), allocatable :: work(:)
    type(face_marks), allocatable :: marks(:)
    real(dp), allocatable :: buffer(:,:,:,:)
    real(dp) :: average(nvar)
    integer :: rank,status,ios,node,parent,low_depth,high_depth,axis,b,c,r,lo(3),hi(3),n(3),p(3),f(3)
    integer :: i,j,k,s,t,dimension,count
    integer(int64) :: elements
    logical :: valid,ready
    ok=.false.
    call validate_mpi_amr_tree_metadata_3d(tree,owners,nvar,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(depth,low_depth,1,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(depth,high_depth,1,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (low_depth/=high_depth.or.depth<1.or.depth>=tree%level_count()) return
    valid=valid_local_tree_flux_level_3d(tree,owners,rank,nvar,depth,flux,refined).and. &
      valid_local_tree_flux_level_3d(tree,owners,rank,nvar,depth-1,flux,refined)
    call agree(valid,ready)
    if (.not.ready) return
    allocate(work(size(flux)),marks(size(flux)),stat=ios)
    call agree(ios==0,ready)
    if (.not.ready) return
    do node=1,size(flux)
      if (tree%nodes(node)%level/=depth-1.or.owners(node)/=rank) cycle
      allocate(work(node)%x,source=flux(node)%x,stat=ios)
      if (ios/=0) exit
      allocate(work(node)%y,source=flux(node)%y,stat=ios)
      if (ios/=0) exit
      allocate(work(node)%z,source=flux(node)%z,stat=ios)
      if (ios/=0) exit
      allocate(marks(node)%x,source=refined(node)%x,stat=ios)
      if (ios/=0) exit
      allocate(marks(node)%y,source=refined(node)%y,stat=ios)
      if (ios/=0) exit
      allocate(marks(node)%z,source=refined(node)%z,stat=ios)
      if (ios/=0) exit
    end do
    call agree(ios==0,ready)
    if (.not.ready) return
    do node=2,size(flux)
      if (tree%nodes(node)%level/=depth) cycle
      parent=tree%nodes(node)%parent
      associate(box=>tree%nodes(node)%patch)
        lo=[box%coarse_i_lower,box%coarse_j_lower,box%coarse_k_lower]
        hi=[box%coarse_i_upper,box%coarse_j_upper,box%coarse_k_upper]
        r=box%refinement_ratio
      end associate
      do axis=1,3
        b=modulo(axis,3)+1;c=modulo(axis+1,3)+1
        n=hi-lo+1;n(axis)=n(axis)+1
        elements=int(nvar,int64)
        do dimension=1,3
          if (elements>int(huge(1),int64)/int(n(dimension),int64)) exit
          elements=elements*int(n(dimension),int64)
        end do
        if (dimension<=3) return ! Count geometry agrees on all ranks.
        count=int(elements);valid=.true.
        if (rank==owners(node).or.rank==owners(parent)) then
          allocate(buffer(nvar,n(1),n(2),n(3)),stat=ios)
          valid=ios==0
        end if
        call agree(valid,ready)
        if (.not.ready) return
        if (rank==owners(node)) then
          do k=1,n(3)
            do j=1,n(2)
              do i=1,n(1)
                p=[i,j,k];average=0
                do t=0,r-1
                  do s=0,r-1
                    f=(p-1)*r+1;f(b)=f(b)+s;f(c)=f(c)+t
                    average=average+get_face(flux(node),axis,f)/real(r,dp)/real(r,dp)
                  end do
                end do
                buffer(:,i,j,k)=average
              end do
            end do
          end do
          valid=all(ieee_is_finite(buffer))
        end if
        call agree(valid,ready)
        if (.not.ready) return
        status=MPI_SUCCESS
        if (owners(node)/=owners(parent)) then
          if (rank==owners(node)) call MPI_Send(buffer,count,MPI_DOUBLE_PRECISION,owners(parent),28406,comm,status)
          if (rank==owners(parent)) call MPI_Recv(buffer,count,MPI_DOUBLE_PRECISION,owners(node),28406, &
            comm,MPI_STATUS_IGNORE,status)
        end if
        valid=status==MPI_SUCCESS
        if (rank==owners(parent).and.valid) then
          do k=1,n(3)
            do j=1,n(2)
              do i=1,n(1)
                p=lo+[i,j,k]-1
                if (marked(marks(parent),axis,p)) then
                  if (.not.matches(buffer(:,i,j,k),get_face(work(parent),axis,p))) valid=.false.
                else
                  call set_face(work(parent),axis,p,buffer(:,i,j,k))
                  call set_mark(marks(parent),axis,p)
                end if
              end do
            end do
          end do
        end if
        call agree(valid,ready)
        if (.not.ready) return
        if (allocated(buffer)) deallocate(buffer)
      end do
    end do
    do node=1,size(flux)
      if (tree%nodes(node)%level/=depth-1.or.owners(node)/=rank) cycle
      call move_alloc(work(node)%x,flux(node)%x)
      call move_alloc(work(node)%y,flux(node)%y)
      call move_alloc(work(node)%z,flux(node)%z)
      call move_alloc(marks(node)%x,refined(node)%x)
      call move_alloc(marks(node)%y,refined(node)%y)
      call move_alloc(marks(node)%z,refined(node)%z)
    end do
    ok=.true.
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
