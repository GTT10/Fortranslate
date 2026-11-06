module conservative_reflux_budget_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: spread_reflux_return
contains
  ! Return one integrated conserved vector to a connected fine-cell stencil.
  ! Positive edge integral travels edges(1,e) -> edges(2,e). The boundary
  ! supplies integrated_return at seed cells; internal integrals redistribute
  ! it to donors with sufficient density/species inventory. Volumes may differ,
  ! so a composite AMR or cut-cell graph uses the same budget calculation.
  ! No species is clipped and no component is independently rescaled.
  ! Caller fields and edge integrals are published together only on success.
  subroutine spread_reflux_return(state,volumes,edges,seed,integrated_return,first_species,last_species, &
      corrected_state,interior_integrals,ok)
    real(dp), intent(in) :: state(:,:),volumes(:),seed(:),integrated_return(:)
    integer, intent(in) :: edges(:,:),first_species,last_species
    real(dp), intent(inout) :: corrected_state(:,:),interior_integrals(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: capacity(:),weights(:),boundary(:),subtree(:),candidate(:,:),flux(:,:)
    integer, allocatable :: degree(:),offset(:),cursor(:),adjacent(:),edge_id(:),queue(:),parent(:),parent_edge(:)
    real(dp) :: available,required,limit,total,amount,delta
    integer :: nc,nv,ne,ios,i,c,e,a,b,root,anchor,head,tail,slot,v
    logical :: good
    ok=.false.;nv=size(state,1);nc=size(state,2);ne=size(edges,2)
    if(nv<6.or.nc<1.or.nc==huge(1)) return
    if(size(edges,1)/=2.or.ne>(huge(1)-1)/2) return
    if(first_species<2.or.last_species<first_species.or.last_species>nv) return
    if(size(volumes)/=nc.or.size(seed)/=nc.or.size(integrated_return)/=nv) return
    if(any(shape(corrected_state)/=[nv,nc]).or.any(shape(interior_integrals)/=[nv,ne])) return
    if(.not.all(ieee_is_finite(state)).or..not.all(ieee_is_finite(volumes))) return
    if(.not.all(ieee_is_finite(seed)).or..not.all(ieee_is_finite(integrated_return))) return
    if(any(volumes<=0).or.any(seed<0).or.any(state(1,:)<=0)) return
    if(any(state(first_species:last_species,:)<0)) return
    if(.not.any(seed>0)) return
    do e=1,ne
      if(any(edges(:,e)<1).or.any(edges(:,e)>nc)) return
      if(edges(1,e)==edges(2,e)) return
    enddo
    allocate(capacity(nc),weights(nc),boundary(nc),subtree(nc),candidate(nv,nc),flux(nv,ne), &
      degree(nc),offset(nc+1),cursor(nc),adjacent(2*ne),edge_id(2*ne),queue(nc),parent(nc), &
      parent_edge(nc),stat=ios)
    if(ios/=0) return
    ! Scale before summation to avoid overflow for physical face-area weights.
    boundary=seed/maxval(seed);total=sum(boundary)
    boundary=boundary/total
    anchor=maxloc(boundary,dim=1)
    boundary(anchor)=0;boundary(anchor)=1-sum(boundary)
    degree=0
    do e=1,ne
      degree(edges(1,e))=degree(edges(1,e))+1
      degree(edges(2,e))=degree(edges(2,e))+1
    enddo
    offset(1)=1
    do i=1,nc
      offset(i+1)=offset(i)+degree(i)
    enddo
    cursor=offset(:nc)
    do e=1,ne
      a=edges(1,e);b=edges(2,e)
      slot=cursor(a);adjacent(slot)=b;edge_id(slot)=e;cursor(a)=slot+1
      slot=cursor(b);adjacent(slot)=a;edge_id(slot)=e;cursor(b)=slot+1
    enddo
    ! Deterministic spanning tree. Disconnected donors cannot fund this face.
    root=findloc(seed>0,.true.,dim=1);queue(1)=root;parent=0;parent_edge=0;parent(root)=root
    head=1;tail=1
    do while(head<=tail)
      a=queue(head);head=head+1
      do slot=offset(a),offset(a+1)-1
        b=adjacent(slot)
        if(parent(b)/=0) cycle
        tail=tail+1;queue(tail)=b;parent(b)=a;parent_edge(b)=edge_id(slot)
      enddo
    enddo
    if(tail/=nc) return
    capacity=1
    do i=1,nc
      do c=1,nv
        if(c/=1.and.(c<first_species.or.c>last_species)) cycle
        if(integrated_return(c)>=0) cycle
        required=-integrated_return(c)
        call multiply(state(c,i),volumes(i),available,good)
        if(.not.good) return
        if(available>required) then
          limit=1
        else
          limit=available/required
          ! Density must remain strictly positive. Species may reach zero.
          if(c==1) limit=limit*(1-16*epsilon(1._dp))
        endif
        capacity(i)=min(capacity(i),limit)
      enddo
    enddo
    total=sum(capacity)
    if(total<1) return
    if(any(integrated_return<0)) then
      weights=capacity/total
    else
      weights=volumes/maxval(volumes);weights=weights/sum(weights)
    endif
    anchor=maxloc(weights,dim=1)
    weights(anchor)=0;weights(anchor)=1-sum(weights)
    if(any(weights<0).or.any(weights>capacity)) return
    candidate=state
    do i=1,nc
      do c=1,nv
        call multiply(weights(i),integrated_return(c),amount,good)
        if(.not.good) return
        if(volumes(i)<1) then
          if(abs(amount)>huge(1._dp)*volumes(i)) return
        endif
        delta=amount/volumes(i)
        call add(state(c,i),delta,candidate(c,i),good)
        if(.not.good) return
      enddo
    enddo
    if(any(candidate(1,:)<=0).or.any(candidate(first_species:last_species,:)<0)) return
    ! Each subtree's desired increment minus its external boundary increment
    ! is carried through its single parent edge. All conserved components use
    ! this same coefficient, retaining density/species and energy relations.
    subtree=weights-boundary;flux=0
    do slot=nc,2,-1
      v=queue(slot);a=parent(v);e=parent_edge(v)
      amount=subtree(v)
      if(edges(2,e)/=v) amount=-amount
      do c=1,nv
        call multiply(amount,integrated_return(c),flux(c,e),good)
        if(.not.good) return
      enddo
      subtree(a)=subtree(a)+subtree(v)
    enddo
    corrected_state=candidate;interior_integrals=flux;ok=.true.
  end subroutine spread_reflux_return

  pure subroutine multiply(a,b,value,ok)
    real(dp), intent(in) :: a,b
    real(dp), intent(out) :: value
    logical, intent(out) :: ok
    ok=.false.;value=0
    if(.not.ieee_is_finite(a).or..not.ieee_is_finite(b)) return
    if(abs(a)>1) then
      if(abs(b)>huge(1._dp)/abs(a)) return
    else if(abs(b)>1) then
      if(abs(a)>huge(1._dp)/abs(b)) return
    endif
    value=a*b;ok=ieee_is_finite(value)
  end subroutine multiply
  pure subroutine add(a,b,value,ok)
    real(dp), intent(in) :: a,b
    real(dp), intent(out) :: value
    logical, intent(out) :: ok
    ok=.false.;value=0
    if(.not.ieee_is_finite(a).or..not.ieee_is_finite(b)) return
    if(b>0) then
      if(a>huge(1._dp)-b) return
    else
      if(a< -huge(1._dp)-b) return
    endif
    value=a+b;ok=ieee_is_finite(value)
  end subroutine add
end module conservative_reflux_budget_mod
