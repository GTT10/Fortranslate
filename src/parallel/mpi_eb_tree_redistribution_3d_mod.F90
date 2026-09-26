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
  public :: redistribute_sdc_mpi_eb_tree_state_level_rate_3d
  public :: redistribute_mpi_eb_tree_state_3d,redistribute_sdc_mpi_eb_tree_state_rate_3d
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
  subroutine redistribute_sdc_mpi_eb_tree_state_level_rate_3d(tree,owners,geometry,state,rhs,nvar,depth,dt, &
      rate,transfer,comm,ok,target_volume_fraction,boundary_faces)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar,depth
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,rhs
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate,transfer
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: target_volume_fraction
    type(amr_tree_state_3d) :: selected,result,exchange,correction
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
    call redistribute_mpi_eb_tree_state_3d(tree,owners,geometry,selected,nvar,result,comm,ready, &
      target_volume_fraction,boundary_faces)
    if(.not.ready) return
    call redistribute_mpi_eb_tree_state_3d(tree,owners,geometry,state,nvar,correction,comm,ready, &
      target_volume_fraction,boundary_faces,depth)
    if(.not.ready) return
    exchange=result;valid=.true.
    do p=1,size(tree%nodes)
      if(owners(p)/=rank) cycle
      if(tree%nodes(p)%level==depth) then
        exchange%fields(p)%values=0
        result%fields(p)%values=result%fields(p)%values+correction%fields(p)%values/dt
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
      exchange%fields(p)%values=dt*exchange%fields(p)%values+correction%fields(p)%values
      result%fields(p)%values=0
    end do
    do p=1,size(tree%nodes)
      if(owners(p)/=rank) cycle
      valid=valid.and.all(ieee_is_finite(result%fields(p)%values))
      valid=valid.and.all(ieee_is_finite(exchange%fields(p)%values))
    end do
    call agree(valid,ready)
    if(.not.ready) return
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
  end subroutine redistribute_sdc_mpi_eb_tree_state_level_rate_3d
  ! Replicated geometry determines identical weighted neighborhoods. State and
  ! neighborhood means stay owner-local; two sparse fetches supply remote values.
  ! correction_depth returns signed displacements from neighborhoods at one
  ! collectively selected depth, rather than the full redistributed state.
  subroutine redistribute_mpi_eb_tree_state_3d(tree,owners,geometry,input,nv,result,comm,ok, &
      target_volume_fraction,boundary_faces,correction_depth)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nv
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(MPI_Comm), intent(in) :: comm
    type(amr_tree_state_3d), intent(in) :: input
    type(amr_tree_state_3d), intent(inout) :: result
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: target_volume_fraction
    character(len=*), optional, intent(in) :: boundary_faces(6)
    integer, optional, intent(in) :: correction_depth
    integer :: depth,depth_low,depth_high
    type :: neighborhood
      integer, allocatable :: member(:)
    end type
    type(neighborhood), allocatable :: groups(:)
    type(amr_tree_state_3d) :: ids,candidate,averages
    integer(int64), allocatable :: origin(:,:),domain(:,:),scale(:)
    integer, allocatable :: ratios(:),cells(:,:),receiver(:),memberships(:),queue(:),seen(:),neighbors(:),directions(:)
    real(dp), allocatable :: volume(:),target(:),beta(:),self(:),values(:,:),incoming(:,:),mean(:)
    integer(int64) :: maximum,cell_lo(3),cell_hi(3)
    character(len=32) :: faces(6)
    integer :: p,n(3),i,j,k,child,lo(3),hi(3),nc,root,other,m,count,head,tail,current,allocation_status
    real(dp) :: fraction,total,weight,normal(3),score,best
    logical :: valid,fallback,ready
    integer :: rank,status,q,request_count,controls(6),low(6),high(6)
    integer, allocatable :: requests(:,:),lookup(:),reverse_lookup(:)
    logical, allocatable :: needed(:)
    real(dp) :: fraction_low,fraction_high
    ok=.false.;fraction=.5_dp;faces='outflow'
    if(present(target_volume_fraction)) fraction=target_volume_fraction
    if(present(boundary_faces)) faces=boundary_faces
    valid=ieee_is_finite(fraction)
    if(valid) valid=fraction>0.and.fraction<=1
    call agree(valid.and.nv>0.and.valid_reactive_boundaries_3d(faces),ready)
    if(.not.ready) return
    call MPI_Allreduce(fraction,fraction_low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(fraction,fraction_high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.fraction_low/=fraction_high) return
    controls=merge(1,0,faces=='periodic')+merge(2,0,faces=='inflow')+merge(3,0,faces=='slip_wall')+ &
      merge(4,0,faces=='no_slip_wall')+merge(5,0,faces=='isothermal_wall')
    call MPI_Allreduce(controls,low,6,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,6,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    call validate_mpi_amr_tree_local_3d(tree,owners,input,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if(.not.ready) return
    call agree(periodic_eb_faces_match(geometry(1),faces),ready)
    if(.not.ready) return
    depth=-1
    if(present(correction_depth)) depth=correction_depth
    valid=depth>=-1.and.depth<tree%level_count()
    if(present(correction_depth)) valid=valid.and.depth>=0
    call agree(valid,ready)
    if(.not.ready) return
    call MPI_Allreduce(depth,depth_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(depth,depth_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.depth_low/=depth_high) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready)
    if(.not.ready) return
    maximum=0
    do p=1,size(tree%nodes)
      n=tree%nodes(p)%shape
      if(int(n(1),int64)*n(2)>int(huge(1),int64)/n(3)) return
      maximum=maximum+int(n(1),int64)*n(2)*n(3)
      if(maximum>huge(1)) return
    end do
    allocate(ids%fields(size(tree%nodes)),scale(size(ratios)),cells(4,int(maximum)),stat=allocation_status)
    call agree(allocation_status==0,ready)
    if(.not.ready) return
    scale=domain(1,size(ratios))/domain(1,:)
    candidate=input;averages=input;nc=0;valid=.true.
    do p=1,size(tree%nodes)
      n=tree%nodes(p)%shape
      allocate(ids%fields(p)%values(1,n(1),n(2),n(3)),source=0._dp,stat=allocation_status)
      call agree(allocation_status==0,ready)
      if(.not.ready) return
      if(owners(p)==rank) then
        averages%fields(p)%values=0
        if(depth>=0) candidate%fields(p)%values=0
      end if
      do child=p+1,size(tree%nodes)
        if(tree%nodes(child)%parent/=p) cycle
        associate(b=>tree%nodes(child)%patch)
          lo=[b%coarse_i_lower,b%coarse_j_lower,b%coarse_k_lower]
          hi=[b%coarse_i_upper,b%coarse_j_upper,b%coarse_k_upper]
        end associate
        ids%fields(p)%values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=-1
      end do
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(geometry(p)%volume_fraction(i,j,k)<=0) then
          if(owners(p)==rank) candidate%fields(p)%values(:,i,j,k)=0
          ids%fields(p)%values(1,i,j,k)=0
          cycle
        end if
        if(owners(p)==rank) then
          if(.not.all(ieee_is_finite(input%fields(p)%values(:,i,j,k)))) valid=.false.
        end if
        if(ids%fields(p)%values(1,i,j,k)<0) cycle
        nc=nc+1;cells(:,nc)=[p,i,j,k];ids%fields(p)%values(1,i,j,k)=real(nc,dp)
      end do
      end do
      end do
    end do
    call agree(valid,ready)
    if(.not.ready) return
    if(nc>huge(1)/max(6,nv)) return
    allocate(groups(nc),receiver(nc),memberships(nc),queue(nc),seen(nc),neighbors(6*nc),directions(6*nc), &
      volume(nc),target(nc),beta(nc),self(nc),mean(nv),lookup(nc),reverse_lookup(nc),needed(nc),stat=allocation_status)
    call agree(allocation_status==0,ready)
    if(.not.ready) return
    do root=1,nc
      p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
      associate(g=>geometry(p))
        target(root)=fraction*g%dx*g%dy*g%dz
        volume(root)=g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz
      end associate
      if(min(volume(root),target(root))<=0) return
    end do
    receiver=0;memberships=1;self=1;beta=0;fallback=.false.
    do root=1,nc
      if(volume(root)>=target(root)) cycle
      current=root
      call adjacent(count)
      p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
      associate(g=>geometry(p))
        normal=[g%boundary_normal_x(i,j,k),g%boundary_normal_y(i,j,k),g%boundary_normal_z(i,j,k)]
      end associate
      best=-1
      do m=1,count
        other=neighbors(m)
        if(volume(other)<max(target(root),target(other))) cycle
        score=normal(abs(directions(m)))*real(sign(1,directions(m)),dp)
        if(score<=512*epsilon(1._dp).or.score<=best) cycle
        best=score;receiver(root)=other
      end do
      if(receiver(root)==0) then
        fallback=.true.;exit
      end if
      other=receiver(root)
      beta(root)=(target(root)-volume(root))/volume(other)
      memberships(other)=memberships(other)+1
    end do
    if(.not.fallback) then
      do root=1,nc
        other=receiver(root)
        if(other==0) cycle
        self(other)=self(other)-beta(root)/real(memberships(other),dp)
      end do
    else
      memberships=1;seen=0
      do root=1,nc
        if(volume(root)>=target(root)) cycle
        head=1;tail=1;queue(1)=root;seen(root)=root;total=volume(root)
        do while(total<target(root).and.head<=tail)
          current=queue(head)
          call adjacent(count)
          do m=1,count
            if(total>=target(root)) exit
            other=neighbors(m)
            if(seen(other)==root) cycle
            tail=tail+1;queue(tail)=other;seen(other)=root;total=total+volume(other)
          end do
          head=head+1
        end do
        if(total<target(root)) return
        allocate(groups(root)%member(tail),stat=allocation_status)
        call agree(allocation_status==0,ready)
        if(.not.ready) return
        groups(root)%member=queue(:tail)
        do m=2,tail
          other=queue(m);memberships(other)=memberships(other)+1
        end do
      end do
    end if
    ! Fetch only members required for neighborhood averages owned by this rank.
    needed=.false.
    do root=1,nc
      if(owners(cells(1,root))/=rank) cycle
      if(depth>=0.and.tree%nodes(cells(1,root))%level/=depth) cycle
      needed(root)=.true.
      if(fallback) then
        if(allocated(groups(root)%member)) needed(groups(root)%member)=.true.
      else
        other=receiver(root)
        if(other>0) needed(other)=.true.
      end if
    end do
    call make_requests(lookup,ready)
    if(.not.ready) return
    call fetch_mpi_amr_tree_cells_3d(tree,owners,input,nv,requests,values,comm,ready)
    if(.not.ready) return
    valid=.true.
    do root=1,nc
      if(owners(cells(1,root))/=rank) cycle
      if(depth>=0.and.tree%nodes(cells(1,root))%level/=depth) cycle
      if(fallback) then
        count=1
        if(allocated(groups(root)%member)) count=size(groups(root)%member)
        total=0;mean=0
        do m=1,count
          other=root
          if(allocated(groups(root)%member)) other=groups(root)%member(m)
          weight=volume(other)/real(memberships(other),dp)
          total=total+weight;mean=mean+weight*values(:,lookup(other))
        end do
      else
        total=self(root)*volume(root);mean=total*values(:,lookup(root))
        other=receiver(root)
        if(other>0) then
          weight=beta(root)*volume(other)/real(memberships(other),dp)
          total=total+weight;mean=mean+weight*values(:,lookup(other))
        end if
        if(self(root)<0) valid=.false.
      end if
      if(total<=0) then
        valid=.false.;cycle
      end if
      p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
      averages%fields(p)%values(:,i,j,k)=mean/total
      valid=valid.and.all(ieee_is_finite(averages%fields(p)%values(:,i,j,k)))
    end do
    call agree(valid,ready)
    if(.not.ready) return
    ! Fetch neighborhood means whose scatter reaches an owned recipient.
    needed=.false.
    do root=1,nc
      if(depth>=0.and.tree%nodes(cells(1,root))%level/=depth) cycle
      needed(root)=owners(cells(1,root))==rank
      if(fallback) then
        if(allocated(groups(root)%member)) then
          do m=1,size(groups(root)%member)
            other=groups(root)%member(m)
            if(owners(cells(1,other))==rank) needed(root)=.true.
          end do
        end if
      else
        other=receiver(root)
        if(other>0) then
          if(owners(cells(1,other))==rank) needed(root)=.true.
        end if
      end if
    end do
    call make_requests(reverse_lookup,ready)
    if(.not.ready) return
    call fetch_mpi_amr_tree_cells_3d(tree,owners,averages,nv,requests,incoming,comm,ready)
    if(.not.ready) return
    do root=1,nc
      p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
      if(owners(p)==rank) candidate%fields(p)%values(:,i,j,k)=0
    end do
    do root=1,nc
      q=reverse_lookup(root)
      if(q==0) cycle
      if(fallback) then
        count=1
        if(allocated(groups(root)%member)) count=size(groups(root)%member)
        do m=1,count
          other=root
          if(allocated(groups(root)%member)) other=groups(root)%member(m)
          p=cells(1,other);i=cells(2,other);j=cells(3,other);k=cells(4,other)
          if(owners(p)/=rank) cycle
          mean=incoming(:,q)
          if(depth>=0) mean=mean-input%fields(p)%values(:,i,j,k)
          candidate%fields(p)%values(:,i,j,k)=candidate%fields(p)%values(:,i,j,k)+ &
            mean/real(memberships(other),dp)
        end do
      else
        p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
        if(owners(p)==rank) then
          mean=incoming(:,q)
          if(depth>=0) mean=mean-input%fields(p)%values(:,i,j,k)
          candidate%fields(p)%values(:,i,j,k)=candidate%fields(p)%values(:,i,j,k)+ &
            self(root)*real(memberships(root),dp)*mean
        end if
        other=receiver(root)
        if(other>0) then
          p=cells(1,other);i=cells(2,other);j=cells(3,other);k=cells(4,other)
          if(owners(p)==rank) then
            mean=incoming(:,q)
            if(depth>=0) mean=mean-input%fields(p)%values(:,i,j,k)
            candidate%fields(p)%values(:,i,j,k)=candidate%fields(p)%values(:,i,j,k)+beta(root)*mean
          end if
        end if
      end if
    end do
    if(.not.fallback) then
      do root=1,nc
        p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
        if(owners(p)==rank) candidate%fields(p)%values(:,i,j,k)= &
          candidate%fields(p)%values(:,i,j,k)/real(memberships(root),dp)
      end do
    end if
    valid=.true.
    do p=1,size(tree%nodes)
      if(owners(p)/=rank) cycle
      valid=valid.and.all(ieee_is_finite(candidate%fields(p)%values))
    end do
    call agree(valid,ready)
    if(.not.ready) return
    call move_alloc(candidate%fields,result%fields);ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine agree
    subroutine make_requests(map,success)
      integer, intent(out) :: map(:)
      logical, intent(out) :: success
      integer :: t
      request_count=count_needed()
      if(allocated(requests)) deallocate(requests)
      allocate(requests(4,request_count),stat=allocation_status)
      call agree(allocation_status==0,success)
      if(.not.success) return
      map=0;request_count=0
      do t=1,nc
        if(.not.needed(t)) cycle
        request_count=request_count+1;map(t)=request_count
        requests(:,request_count)=cells(:,t)
      end do
    end subroutine make_requests
    integer function count_needed() result(number)
      integer :: t
      number=0
      do t=1,nc
        if(needed(t)) number=number+1
      end do
    end function count_needed
    subroutine adjacent(count)
      integer, intent(out) :: count
      integer :: axis,side,q,direction,lower(3),upper(3),xyz(3),a,b,c,face_cell
      integer(int64) :: face,s,qlo(3),qhi(3),first(3),last(3)
      integer :: neighbor
      p=cells(1,current);i=cells(2,current);j=cells(3,current);k=cells(4,current)
      cell_lo=(origin(:,p)+[i-1,j-1,k-1])*scale(tree%nodes(p)%level+1)
      cell_hi=cell_lo+scale(tree%nodes(p)%level+1)
      count=0
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
            if(ids%fields(q)%values(1,a,b,c)<=0) cycle
            if(.not.open_face(geometry(q),xyz,axis,1-side)) cycle
            neighbor=int(ids%fields(q)%values(1,a,b,c))
            count=count+1;neighbors(count)=neighbor;directions(count)=axis*(2*side-1)
          end do
          end do
          end do
        end do
      end do
      end do
    end subroutine adjacent
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
  end subroutine redistribute_mpi_eb_tree_state_3d
  ! Affine SDC convention, matching S(U+dt*R)-U over dt in the serial kernel.
  subroutine redistribute_sdc_mpi_eb_tree_state_rate_3d(tree,owners,geometry,state,rhs,nv,dt, &
      rate,comm,ok,target_volume_fraction,boundary_faces)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nv
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,rhs
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: target_volume_fraction
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: provisional,candidate
    integer :: p,n(3),i,j,k,rank,status
    real(dp) :: low,high
    logical :: valid,ready
    ok=.false.;valid=ieee_is_finite(dt)
    if(valid) valid=dt>0
    call agree(valid,ready)
    if(.not.ready) return
    call MPI_Allreduce(dt,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low/=high) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,rhs,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    provisional=rhs;valid=.true.
    do p=1,size(tree%nodes)
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        provisional%fields(p)%values(:,i,j,k)=0
        if(geometry(p)%volume_fraction(i,j,k)<=0) cycle
        if(.not.all(ieee_is_finite(rhs%fields(p)%values(:,i,j,k))).or. &
           .not.all(ieee_is_finite(state%fields(p)%values(:,i,j,k)))) then
          valid=.false.;cycle
        end if
        provisional%fields(p)%values(:,i,j,k)=state%fields(p)%values(:,i,j,k)+dt*rhs%fields(p)%values(:,i,j,k)
      end do
      end do
      end do
    end do
    call agree(valid,ready)
    if(.not.ready) return
    call redistribute_mpi_eb_tree_state_3d(tree,owners,geometry,provisional,nv,candidate,comm,ready, &
      target_volume_fraction,boundary_faces)
    if(.not.ready) return
    do p=1,size(tree%nodes)
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(geometry(p)%volume_fraction(i,j,k)<=0) cycle
        candidate%fields(p)%values(:,i,j,k)=(candidate%fields(p)%values(:,i,j,k)-state%fields(p)%values(:,i,j,k))/dt
      end do
      end do
      end do
      valid=valid.and.all(ieee_is_finite(candidate%fields(p)%values))
    end do
    call agree(valid,ready)
    if(.not.ready) return
    call move_alloc(candidate%fields,rate%fields);ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine agree
  end subroutine redistribute_sdc_mpi_eb_tree_state_rate_3d
end module mpi_eb_tree_redistribution_3d_mod
