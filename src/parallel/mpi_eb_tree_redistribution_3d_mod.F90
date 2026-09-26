module mpi_eb_tree_redistribution_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,periodic_eb_faces_match
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: fetch_mpi_amr_tree_cells_3d,validate_mpi_amr_tree_local_3d
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d
  implicit none
  private
  public :: redistribute_sdc_mpi_eb_tree_flux_rate_3d,redistribute_sdc_mpi_eb_tree_level_rate_3d
contains
  ! Two sparse exchanges: neighbor rates for each owned cut leaf, followed by
  ! donor excess for each owned recipient. No gas/rate tree is gathered.
  subroutine redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree,owners,geometry,rhs,nvar,rate,comm,ok,boundary_faces)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: rhs
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    character(len=32) :: faces(6)
    type(amr_tree_state_3d) :: candidate,active,excess,delta
    integer(int64), allocatable :: origin(:,:),domain(:,:),scale(:)
    integer, allocatable :: ratios(:),requests(:,:),reverse(:,:),destination(:,:)
    integer(int64) :: cell_lo(3),cell_hi(3),nf,nr
    real(dp), allocatable :: values(:,:),incoming(:,:),average(:)
    integer :: p,child,n(3),lo(3),hi(3),i,j,k,stage,rank,status,ios,q
    integer :: controls(6),low(6),high(6)
    real(dp) :: kappa,volume,neighbors,total
    logical :: valid,ready
    ok=.false.;faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    call agree(nvar>0.and.valid_reactive_boundaries_3d(faces),ready)
    if(.not.ready) return
    controls=merge(1,0,faces=='periodic')+merge(2,0,faces=='inflow')+merge(3,0,faces=='slip_wall')+ &
      merge(4,0,faces=='no_slip_wall')+merge(5,0,faces=='isothermal_wall')
    call MPI_Allreduce(controls,low,6,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,6,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    call validate_mpi_amr_tree_local_3d(tree,owners,rhs,nvar,comm,ready)
    if(.not.ready) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nvar,comm,ready)
    if(.not.ready) return
    call agree(periodic_eb_faces_match(geometry(1),faces),ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready)
    if(.not.ready) return
    candidate=rhs
    allocate(active%fields(size(tree%nodes)),excess%fields(size(tree%nodes)), &
      delta%fields(size(tree%nodes)),scale(size(ratios)),average(nvar),stat=ios)
    call agree(ios==0,ready)
    if(.not.ready) return
    scale=domain(1,size(ratios))/domain(1,:)
    valid=.true.
    do p=1,size(tree%nodes)
      n=tree%nodes(p)%shape
      allocate(active%fields(p)%values(1,n(1),n(2),n(3)),source=1._dp,stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      if(owners(p)==rank) then
        allocate(excess%fields(p)%values(nvar,n(1),n(2),n(3)), &
          delta%fields(p)%values(nvar,n(1),n(2),n(3)),source=0._dp,stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        end if
      end if
      do child=p+1,size(tree%nodes)
        if(tree%nodes(child)%parent/=p) cycle
        associate(b=>tree%nodes(child)%patch)
          lo=[b%coarse_i_lower,b%coarse_j_lower,b%coarse_k_lower]
          hi=[b%coarse_i_upper,b%coarse_j_upper,b%coarse_k_upper]
        end associate
        active%fields(p)%values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=0
      end do
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(geometry(p)%volume_fraction(i,j,k)<=0) then
          active%fields(p)%values(1,i,j,k)=0
          if(owners(p)==rank) candidate%fields(p)%values(:,i,j,k)=0
        else if(owners(p)==rank) then
          if(.not.all(ieee_is_finite(rhs%fields(p)%values(:,i,j,k)))) valid=.false.
        end if
      end do
      end do
      end do
    end do
    call agree(valid,ready)
    if(.not.ready) return
    do stage=0,2
      nf=0;nr=0
      do p=1,size(tree%nodes)
        n=tree%nodes(p)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(active%fields(p)%values(1,i,j,k)==0) cycle
          kappa=geometry(p)%volume_fraction(i,j,k)
          if(kappa>=1) cycle
          cell_lo=(origin(:,p)+[i-1,j-1,k-1])*scale(tree%nodes(p)%level+1)
          cell_hi=cell_lo+scale(tree%nodes(p)%level+1)
          volume=kappa*geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
          neighbors=0
          if(stage==2.and.owners(p)==rank) average=volume*rhs%fields(p)%values(:,i,j,k)
          call visit_neighbors()
          total=volume+neighbors
          if(.not.ieee_is_finite(total).or.volume<=0.or.neighbors<=0) valid=.false.
          if(stage==2.and.owners(p)==rank.and.valid) then
            average=average/total
            candidate%fields(p)%values(:,i,j,k)=kappa*rhs%fields(p)%values(:,i,j,k)+(1-kappa)*average
            excess%fields(p)%values(:,i,j,k)=volume*(1-kappa)*(rhs%fields(p)%values(:,i,j,k)-average)/neighbors
          end if
        end do
        end do
        end do
      end do
      call agree(valid.and.max(nf,nr)<=int(huge(1)/max(4,nvar),int64),ready)
      if(.not.ready) return
      if(stage==0) then
        allocate(requests(4,int(nf)),reverse(4,int(nr)),destination(4,int(nr)),stat=ios)
        call agree(ios==0,ready)
        if(.not.ready) return
      else if(stage==1) then
        call fetch_mpi_amr_tree_cells_3d(tree,owners,rhs,nvar,requests,values,comm,ready)
        if(.not.ready) return
      end if
    end do
    call fetch_mpi_amr_tree_cells_3d(tree,owners,excess,nvar,reverse,incoming,comm,ready)
    if(.not.ready) return
    do q=1,size(destination,2)
      p=destination(1,q);i=destination(2,q);j=destination(3,q);k=destination(4,q)
      delta%fields(p)%values(:,i,j,k)=delta%fields(p)%values(:,i,j,k)+incoming(:,q)
    end do
    do p=1,size(tree%nodes)
      if(owners(p)/=rank) cycle
      candidate%fields(p)%values=candidate%fields(p)%values+delta%fields(p)%values
      valid=valid.and.all(ieee_is_finite(candidate%fields(p)%values))
    end do
    call agree(valid,ready)
    if(.not.ready) return
    call move_alloc(candidate%fields,rate%fields)
    ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine agree
    subroutine visit_neighbors()
      integer :: axis,side,q,direction,lower(3),upper(3),xyz(3),a,b,c,face_cell
      integer(int64) :: face,s,qlo(3),qhi(3),first(3),last(3)
      real(dp) :: neighbor_volume
      do axis=1,3
      do side=0,1
        if(.not.open_face(geometry(p),[i,j,k],axis,side)) cycle
        face=cell_lo(axis)
        if(side==1) face=cell_hi(axis)
        if(faces(2*axis-1)=='periodic') then
          if(side==0.and.face==0) face=domain(axis,size(ratios))
          if(side==1.and.face==domain(axis,size(ratios))) face=0
        end if
        do q=1,size(tree%nodes)
          s=scale(tree%nodes(q)%level+1)
          ! A leaf can touch this face only on an integer cell boundary.
          if(modulo(face,s)/=0) cycle
          qlo=origin(:,q)*s;qhi=(origin(:,q)+tree%nodes(q)%shape)*s
          if(side==0.and.(face<=qlo(axis).or.face>qhi(axis))) cycle
          if(side==1.and.(face<qlo(axis).or.face>=qhi(axis))) cycle
          first=0;last=0
          do direction=1,3
            if(direction==axis) cycle
            first(direction)=max(cell_lo(direction),qlo(direction))
            last(direction)=min(cell_hi(direction),qhi(direction))
          end do
          if(any([(last(direction)<=first(direction).and.direction/=axis,direction=1,3)])) cycle
          do direction=1,3
            if(direction==axis) cycle
            lower(direction)=int(first(direction)/s-origin(direction,q))+1
            upper(direction)=int((last(direction)-1)/s-origin(direction,q))+1
          end do
          face_cell=int(face/s-origin(axis,q))+side
          lower(axis)=face_cell;upper(axis)=face_cell
          do c=lower(3),upper(3)
          do b=lower(2),upper(2)
          do a=lower(1),upper(1)
            xyz=[a,b,c]
            if(active%fields(q)%values(1,a,b,c)==0) cycle
            if(.not.open_face(geometry(q),xyz,axis,1-side)) cycle
            neighbor_volume=geometry(q)%volume_fraction(a,b,c)* &
              geometry(q)%dx*geometry(q)%dy*geometry(q)%dz
            neighbors=neighbors+neighbor_volume
            if(owners(p)==rank) then
              nf=nf+1
              if(stage==1) requests(:,nf)=[q,a,b,c]
              if(stage==2) average=average+neighbor_volume*values(:,nf)
            end if
            if(owners(q)==rank) then
              nr=nr+1
              if(stage==1) then
                reverse(:,nr)=[p,i,j,k]
                destination(:,nr)=[q,a,b,c]
              end if
            end if
          end do
          end do
          end do
        end do
      end do
      end do
    end subroutine visit_neighbors
    logical function open_face(g,cell,axis,side) result(opened)
      type(eb_geometry_3d), intent(in) :: g
      integer, intent(in) :: cell(3),axis,side
      integer :: f(3)
      f=cell;f(axis)=f(axis)-1+side
      select case(axis)
      case(1)
        opened=g%x_face_fraction(f(1),f(2),f(3))>0
      case(2)
        opened=g%y_face_fraction(f(1),f(2),f(3))>0
      case(3)
        opened=g%z_face_fraction(f(1),f(2),f(3))>0
      end select
    end function open_face
  end subroutine redistribute_sdc_mpi_eb_tree_flux_rate_3d
  subroutine redistribute_sdc_mpi_eb_tree_level_rate_3d(tree,owners,geometry,rhs,nvar,depth,dt, &
      rate,transfer,comm,ok,boundary_faces)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar,depth
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: rhs
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate,transfer
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: selected,result,exchange
    integer :: p,rank,status,lower,upper
    real(dp) :: dt_low,dt_high
    logical :: valid,ready
    ok=.false.
    valid=ieee_is_finite(dt)
    if(valid) valid=dt>0
    call agree(valid,ready)
    if(.not.ready) return
    call MPI_Allreduce(depth,lower,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(depth,upper,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.lower/=upper) return
    call MPI_Allreduce(dt,dt_low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,dt_high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.dt_low/=dt_high) return
    call validate_mpi_amr_tree_local_3d(tree,owners,rhs,nvar,comm,ready)
    if(.not.ready) return
    call agree(depth>=0.and.depth<tree%level_count(),ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    selected=rhs
    do p=1,size(tree%nodes)
      if(owners(p)/=rank) cycle
      if(tree%nodes(p)%level/=depth) selected%fields(p)%values=0
    end do
    call redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree,owners,geometry,selected,nvar,result,comm,ready,boundary_faces)
    if(.not.ready) return
    exchange=result;valid=.true.
    do p=1,size(tree%nodes)
      if(owners(p)/=rank) cycle
      if(tree%nodes(p)%level==depth) then
        exchange%fields(p)%values=0
      else
        if(dt>1) then
          if(any(abs(exchange%fields(p)%values)>huge(1._dp)/dt)) valid=.false.
        end if
      end if
    end do
    call agree(valid,ready)
    if(.not.ready) return
    do p=1,size(tree%nodes)
      if(owners(p)/=rank.or.tree%nodes(p)%level==depth) cycle
      exchange%fields(p)%values=dt*exchange%fields(p)%values
      result%fields(p)%values=0
    end do
    call move_alloc(result%fields,rate%fields)
    call move_alloc(exchange%fields,transfer%fields)
    ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine agree
  end subroutine redistribute_sdc_mpi_eb_tree_level_rate_3d
end module mpi_eb_tree_redistribution_3d_mod
