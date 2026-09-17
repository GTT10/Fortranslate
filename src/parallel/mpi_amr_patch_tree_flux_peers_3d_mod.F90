module mpi_amr_patch_tree_flux_peers_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,face_marks
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_metadata_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  public :: synchronize_mpi_amr_tree_flux_peers_3d
contains
  subroutine synchronize_mpi_amr_tree_flux_peers_3d(tree,owners,nvar,depth,faces,flux,refined,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar,depth
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_flux_3d), intent(inout) :: flux(:)
    type(face_marks), intent(inout) :: refined(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_flux_3d), allocatable :: work(:)
    type(face_marks), allocatable :: marks(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: a(:,:),b(:,:)
    integer :: rank,status,ios,node,left,right,axis,side,other,n(3),dims(3),hi(3),p(3),q(3),i,j,k,t
    integer :: low,high,periodic(3),minimum(3),maximum(3),count
    integer(int64) :: normal,opposite,lower(3),upper(3),elements
    logical :: valid,ready
    ok=.false.
    call validate_mpi_amr_tree_metadata_3d(tree,owners,nvar,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(depth,low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(depth,high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (low/=high.or.depth<0.or.depth>=tree%level_count()) return
    valid=valid_reactive_boundaries_3d(faces).and.size(flux)==size(tree%nodes).and.size(refined)==size(tree%nodes)
    valid=valid.and.nvar<huge(1)
    call agree(valid,ready)
    if (.not.ready) return
    periodic=merge(1,0,faces(1:6:2)=='periodic')
    call MPI_Allreduce(periodic,minimum,3,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(periodic,maximum,3,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(minimum/=maximum)) return
    valid=.true.
    do node=1,size(flux)
      if (tree%nodes(node)%level/=depth) cycle
      if (owners(node)/=rank) then
        if (allocated(flux(node)%x).or.allocated(flux(node)%y).or.allocated(flux(node)%z)) valid=.false.
        if (allocated(refined(node)%x).or.allocated(refined(node)%y).or.allocated(refined(node)%z)) valid=.false.
        cycle
      end if
      n=tree%nodes(node)%shape
      if (any(n==huge(1))) then
        valid=.false.;cycle
      end if
      if (.not.allocated(flux(node)%x).or..not.allocated(flux(node)%y).or..not.allocated(flux(node)%z).or. &
          .not.allocated(refined(node)%x).or..not.allocated(refined(node)%y).or..not.allocated(refined(node)%z)) then
        valid=.false.;cycle
      end if
      if (any(shape(flux(node)%x)/=[nvar,n(1)+1,n(2),n(3)]).or. &
          any(shape(flux(node)%y)/=[nvar,n(1),n(2)+1,n(3)]).or. &
          any(shape(flux(node)%z)/=[nvar,n(1),n(2),n(3)+1])) valid=.false.
      if (any(shape(refined(node)%x)/=[n(1)+1,n(2),n(3)]).or. &
          any(shape(refined(node)%y)/=[n(1),n(2)+1,n(3)]).or. &
          any(shape(refined(node)%z)/=[n(1),n(2),n(3)+1])) valid=.false.
      if (.not.all(ieee_is_finite(flux(node)%x)).or..not.all(ieee_is_finite(flux(node)%y)).or. &
          .not.all(ieee_is_finite(flux(node)%z))) valid=.false.
    end do
    call agree(valid,ready)
    if (.not.ready) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready)
    if (.not.ready) return
    ! Candidates contain owned active fields only; inactive caller fields are untouched.
    allocate(work(size(flux)),marks(size(flux)),stat=ios)
    call agree(ios==0,ready)
    if (.not.ready) return
    valid=.true.
    do node=1,size(flux)
      if (tree%nodes(node)%level/=depth.or.owners(node)/=rank) cycle
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
    do left=1,size(flux)
      if (tree%nodes(left)%level/=depth) cycle
      do axis=1,3
        do side=1,2
          normal=origin(axis,left)+int((side-1)*tree%nodes(left)%shape(axis),int64)
          if (periodic(axis)==1) normal=modulo(normal,domain(axis,depth+1))
          do right=left,size(flux)
            if (tree%nodes(right)%level/=depth) cycle
            do other=1,2
              if (left==right.and.other<=side) cycle
              opposite=origin(axis,right)+int((other-1)*tree%nodes(right)%shape(axis),int64)
              if (periodic(axis)==1) opposite=modulo(opposite,domain(axis,depth+1))
              if (normal/=opposite) cycle
              lower=max(origin(:,left),origin(:,right))
              upper=min(origin(:,left)+int(tree%nodes(left)%shape,int64), &
                origin(:,right)+int(tree%nodes(right)%shape,int64))-1
              lower(axis)=0;upper(axis)=0
              if (any(upper<lower)) cycle
              dims=int(upper-lower)+1
              elements=int(nvar+1,int64)
              do t=1,3
                if (elements>int(huge(1),int64)/int(dims(t),int64)) exit
                elements=elements*int(dims(t),int64)
              end do
              if (t<=3) return ! Identical geometry/count on every rank.
              count=int(elements);valid=.true.
              if (rank==owners(left).or.rank==owners(right)) then
                allocate(a(nvar+1,product(dims)),b(nvar+1,product(dims)),stat=ios)
                valid=ios==0
              end if
              call agree(valid,ready)
              if (.not.ready) return
              hi=dims
              t=0
              do k=1,hi(3)
                do j=1,hi(2)
                  do i=1,hi(1)
                    t=t+1
                    p=face_index(left,side,[i,j,k])
                    q=face_index(right,other,[i,j,k])
                    if (rank==owners(left)) call read_face(work(left),marks(left),axis,p,a(:,t))
                    if (rank==owners(right)) call read_face(work(right),marks(right),axis,q,b(:,t))
                  end do
                end do
              end do
              status=MPI_SUCCESS
              if (owners(left)/=owners(right)) then
                if (rank==owners(right)) call MPI_Send(b,count,MPI_DOUBLE_PRECISION,owners(left),28404,comm,status)
                if (rank==owners(left)) call MPI_Recv(b,count,MPI_DOUBLE_PRECISION,owners(right),28404, &
                  comm,MPI_STATUS_IGNORE,status)
              end if
              valid=status==MPI_SUCCESS
              if (rank==owners(left).and.valid) then
                do t=1,size(a,2)
                  if (a(nvar+1,t)==b(nvar+1,t)) then
                    if (maxval(abs(a(1:nvar,t)-b(1:nvar,t)))>64*epsilon(1.0_dp)* &
                        max(1.0_dp,maxval(abs(a(1:nvar,t))),maxval(abs(b(1:nvar,t))))) valid=.false.
                  else if (b(nvar+1,t)==1.0_dp) then
                    a(:,t)=b(:,t)
                  end if
                end do
              end if
              call agree(valid,ready)
              if (.not.ready) return
              status=MPI_SUCCESS
              if (owners(left)/=owners(right)) then
                if (rank==owners(left)) call MPI_Send(a,count,MPI_DOUBLE_PRECISION,owners(right),28405,comm,status)
                if (rank==owners(right)) call MPI_Recv(a,count,MPI_DOUBLE_PRECISION,owners(left),28405, &
                  comm,MPI_STATUS_IGNORE,status)
              end if
              call agree(status==MPI_SUCCESS,ready)
              if (.not.ready) return
              t=0
              do k=1,hi(3)
                do j=1,hi(2)
                  do i=1,hi(1)
                    t=t+1
                    p=face_index(left,side,[i,j,k])
                    q=face_index(right,other,[i,j,k])
                    if (rank==owners(left)) call write_face(work(left),marks(left),axis,p,a(:,t))
                    if (rank==owners(right)) call write_face(work(right),marks(right),axis,q,a(:,t))
                  end do
                end do
              end do
              if (allocated(a)) deallocate(a,b)
            end do
          end do
        end do
      end do
    end do
    do node=1,size(flux)
      if (tree%nodes(node)%level/=depth.or.owners(node)/=rank) cycle
      call move_alloc(work(node)%x,flux(node)%x)
      call move_alloc(work(node)%y,flux(node)%y)
      call move_alloc(work(node)%z,flux(node)%z)
      call move_alloc(marks(node)%x,refined(node)%x)
      call move_alloc(marks(node)%y,refined(node)%y)
      call move_alloc(marks(node)%z,refined(node)%z)
    end do
    ok=.true.
  contains
    function face_index(peer,boundary,offset) result(index)
      integer, intent(in) :: peer,boundary,offset(3)
      integer :: index(3),direction
      do direction=1,3
        if (direction==axis) then
          index(direction)=1+(boundary-1)*tree%nodes(peer)%shape(direction)
        else
          index(direction)=int(lower(direction)-origin(direction,peer))+offset(direction)
        end if
      end do
    end function
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine
  subroutine read_face(field,mark,axis,p,value)
    type(amr_tree_flux_3d), intent(in) :: field
    type(face_marks), intent(in) :: mark
    integer, intent(in) :: axis,p(3)
    real(dp), intent(out) :: value(:)
    integer :: n
    n=size(value)-1
    select case(axis)
    case(1)
      value(:n)=field%x(:,p(1),p(2),p(3));value(n+1)=merge(1.0_dp,0.0_dp,mark%x(p(1),p(2),p(3)))
    case(2)
      value(:n)=field%y(:,p(1),p(2),p(3));value(n+1)=merge(1.0_dp,0.0_dp,mark%y(p(1),p(2),p(3)))
    case(3)
      value(:n)=field%z(:,p(1),p(2),p(3));value(n+1)=merge(1.0_dp,0.0_dp,mark%z(p(1),p(2),p(3)))
    end select
  end subroutine
  subroutine write_face(field,mark,axis,p,value)
    type(amr_tree_flux_3d), intent(inout) :: field
    type(face_marks), intent(inout) :: mark
    integer, intent(in) :: axis,p(3)
    real(dp), intent(in) :: value(:)
    integer :: n
    n=size(value)-1
    select case(axis)
    case(1)
      field%x(:,p(1),p(2),p(3))=value(:n);mark%x(p(1),p(2),p(3))=value(n+1)==1
    case(2)
      field%y(:,p(1),p(2),p(3))=value(:n);mark%y(p(1),p(2),p(3))=value(n+1)==1
    case(3)
      field%z(:,p(1),p(2),p(3))=value(:n);mark%z(p(1),p(2),p(3))=value(n+1)==1
    end select
  end subroutine
end module
