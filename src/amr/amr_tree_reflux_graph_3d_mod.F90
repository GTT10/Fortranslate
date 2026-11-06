module amr_tree_reflux_graph_3d_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  type :: cell_lookup
    integer, allocatable :: id(:,:,:)
  end type
  public :: build_amr_reflux_graph_3d
contains
  ! Metadata only: gather no fluid field and make no assumptions about owners.
  ! Integer face intersections retain mixed-level and cross-parent adjacency.
  ! The returned graph contains only the connected fine region touching one
  ! face of target=[node,i,j,k]. Covered cells never become donors.
  subroutine build_amr_reflux_graph_3d(tree,root_spacing,faces,target,axis,side,requests,volumes, &
      edges,edge_side,seed,cell_lower,cell_upper,lattice_spacing,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    real(dp), intent(in) :: root_spacing(3)
    character(len=*), intent(in) :: faces(6)
    integer, intent(in) :: target(4),axis,side
    integer, allocatable, intent(out) :: requests(:,:),edges(:,:),edge_side(:)
    real(dp), allocatable, intent(out) :: volumes(:),seed(:)
    integer(int64), allocatable, intent(out) :: cell_lower(:,:),cell_upper(:,:)
    real(dp), intent(out) :: lattice_spacing(3)
    logical, intent(out) :: ok
    type(cell_lookup), allocatable :: lookup(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:),lower(:,:),upper(:,:)
    integer, allocatable :: ratios(:),all_requests(:,:),all_edges(:,:),all_sides(:), &
      degree(:),offset(:),cursor(:),adjacent(:),queue(:),index(:)
    real(dp), allocatable :: all_seed(:)
    integer(int64) :: total,scale,extent(3),target_lo(3),target_hi(3)
    integer :: node,depth,nc,ne,n(3),i,j,k,c,a,s,e,u,v,ios,head,tail,slot,count
    logical :: ready
    ok=.false.;lattice_spacing=0
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ready)
    if(.not.ready) return
    if(target(1)<1.or.target(1)>size(tree%nodes).or.axis<1.or.axis>3.or.abs(side)/=1) return
    node=target(1);depth=tree%nodes(node)%level
    if(any(target(2:4)<1).or.any(target(2:4)>tree%nodes(node)%shape)) return
    if(covered(tree,node,target(2:4))) return
    if(size(faces)/=6.or.any(root_spacing<=0).or..not.all(ieee_is_finite(root_spacing))) return
    if(.not.valid_reactive_boundaries_3d(faces)) return
    extent=domain(:,size(ratios));lattice_spacing=root_spacing*real(domain(:,1),dp)/real(extent,dp)
    total=0
    do node=1,size(tree%nodes)
      if(tree%nodes(node)%level<=depth) cycle
      n=tree%nodes(node)%shape
      do a=1,3
        if(int(n(a),int64)>int(huge(1),int64)) return
      enddo
      ! Include covered cells in the allocation bound; their ids remain zero.
      scale=int(n(1),int64)*int(n(2),int64)
      if(scale>int(huge(1),int64)/int(n(3),int64)) return
      total=total+scale*int(n(3),int64)
      if(total> (int(huge(1),int64)-1_int64)/12_int64) return
    enddo
    if(total==0) return
    count=int(total)
    allocate(lookup(size(tree%nodes)),all_requests(4,count),lower(3,count),upper(3,count), &
      all_seed(count),all_edges(2,6*count),all_sides(6*count),stat=ios)
    if(ios/=0) return
    nc=0;ne=0;all_seed=0
    do node=1,size(tree%nodes)
      if(tree%nodes(node)%level<=depth) cycle
      n=tree%nodes(node)%shape
      allocate(lookup(node)%id(n(1),n(2),n(3)),stat=ios)
      if(ios/=0) return
      lookup(node)%id=0
      scale=extent(1)/domain(1,tree%nodes(node)%level+1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            if(covered(tree,node,[i,j,k])) cycle
            nc=nc+1;lookup(node)%id(i,j,k)=nc;all_requests(:,nc)=[node,i,j,k]
            lower(:,nc)=(origin(:,node)+int([i,j,k]-1,int64))*scale
            upper(:,nc)=lower(:,nc)+scale
          enddo
        enddo
      enddo
    enddo
    node=target(1);scale=extent(1)/domain(1,depth+1)
    target_lo=(origin(:,node)+int(target(2:4)-1,int64))*scale;target_hi=target_lo+scale
    call neighbors(target_lo,target_hi,axis,side,0)
    if(.not.any(all_seed(:nc)>0)) return
    do c=1,nc
      do a=1,3
        do s=-1,1,2
          call neighbors(lower(:,c),upper(:,c),a,s,c)
        enddo
      enddo
    enddo
    allocate(degree(nc),offset(nc+1),cursor(nc),adjacent(2*ne),queue(nc),index(nc),stat=ios)
    if(ios/=0) return
    degree=0
    do e=1,ne
      degree(all_edges(1,e))=degree(all_edges(1,e))+1
      degree(all_edges(2,e))=degree(all_edges(2,e))+1
    enddo
    offset(1)=1
    do c=1,nc
      offset(c+1)=offset(c)+degree(c)
    enddo
    cursor=offset(:nc)
    do e=1,ne
      u=all_edges(1,e);v=all_edges(2,e)
      slot=cursor(u);adjacent(slot)=v;cursor(u)=slot+1
      slot=cursor(v);adjacent(slot)=u;cursor(v)=slot+1
    enddo
    index=0;queue(1)=findloc(all_seed(:nc)>0,.true.,dim=1);index(queue(1))=1;head=1;tail=1
    do while(head<=tail)
      u=queue(head);head=head+1
      do slot=offset(u),offset(u+1)-1
        v=adjacent(slot)
        if(index(v)/=0) cycle
        tail=tail+1;queue(tail)=v;index(v)=tail
      enddo
    enddo
    ! A face whose seeds belong to disconnected regions needs separate returns.
    if(any(all_seed(:nc)>0.and.index==0)) return
    count=0
    do e=1,ne
      if(index(all_edges(1,e))>0) count=count+1
    enddo
    allocate(requests(4,tail),volumes(tail),seed(tail),cell_lower(3,tail),cell_upper(3,tail), &
      edges(2,count),edge_side(count),stat=ios)
    if(ios/=0) return
    do c=1,tail
      u=queue(c);requests(:,c)=all_requests(:,u);seed(c)=all_seed(u)
      cell_lower(:,c)=lower(:,u);cell_upper(:,c)=upper(:,u)
      volumes(c)=product(real(upper(:,u)-lower(:,u),dp)*lattice_spacing)
    enddo
    count=0
    do e=1,ne
      if(index(all_edges(1,e))==0) cycle
      count=count+1;edges(:,count)=index(all_edges(:,e));edge_side(count)=all_sides(e)
    enddo
    ok=.true.
  contains
    subroutine neighbors(lo,hi,normal,direction,source)
      integer(int64), intent(in) :: lo(3),hi(3)
      integer, intent(in) :: normal,direction,source
      integer(int64) :: plane,w,base(3),end(3),low_cell(3),high_cell(3),overlap(3)
      integer :: peer,lev,ii,jj,kk,dest,b,t
      real(dp) :: area
      plane=lo(normal)
      if(direction>0) plane=hi(normal)
      if(plane==0.or.plane==extent(normal)) then
        if(faces(2*normal-1)/='periodic') return
        if(plane==0) then
          plane=extent(normal)
        else
          plane=0
        endif
      endif
      b=modulo(normal,3)+1;t=modulo(normal+1,3)+1
      do peer=1,size(tree%nodes)
        lev=tree%nodes(peer)%level+1
        if(lev<=depth+1) cycle
        w=extent(1)/domain(1,lev);base=origin(:,peer)*w;end=base+int(tree%nodes(peer)%shape,int64)*w
        if(plane<base(normal).or.plane>end(normal)) cycle
        if(modulo(plane-base(normal),w)/=0) cycle
        if(hi(b)<=base(b).or.lo(b)>=end(b).or.hi(t)<=base(t).or.lo(t)>=end(t)) cycle
        low_cell=(max(lo,base)-base)/w+1;high_cell=(min(hi,end)-1-base)/w+1
        low_cell(normal)=(plane-base(normal))/w+1
        if(direction<0) low_cell(normal)=low_cell(normal)-1
        high_cell(normal)=low_cell(normal)
        if(low_cell(normal)<1_int64.or.low_cell(normal)>int(tree%nodes(peer)%shape(normal),int64)) cycle
        do kk=int(low_cell(3)),int(high_cell(3))
          do jj=int(low_cell(2)),int(high_cell(2))
            do ii=int(low_cell(1)),int(high_cell(1))
              dest=lookup(peer)%id(ii,jj,kk)
              if(dest==0) cycle
              if(source==0) then
                overlap=min(hi,upper(:,dest))-max(lo,lower(:,dest))
                area=real(overlap(b),dp)*lattice_spacing(b)*real(overlap(t),dp)*lattice_spacing(t)
                all_seed(dest)=all_seed(dest)+area
              else if(source<dest) then
                ne=ne+1;all_edges(:,ne)=[source,dest];all_sides(ne)=direction*normal
              endif
            enddo
          enddo
        enddo
      enddo
    end subroutine neighbors
  end subroutine build_amr_reflux_graph_3d
end module amr_tree_reflux_graph_3d_mod
