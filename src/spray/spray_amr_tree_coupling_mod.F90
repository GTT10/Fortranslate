module spray_amr_tree_coupling_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d, &
    average_down_amr_patch_tree_3d
  use spray_amr_tree_partition_mod, only: locate_amr_tree_spray_cells
  implicit none
  private
  public :: amr_tree_spray_stencil,sample_amr_tree_spray_field,deposit_amr_tree_spray_field
  public :: deposit_amr_tree_spray_gas
contains
  ! Keep conserved fields and EOS temperatures private until every patch accepts
  ! the feedback, including parents refreshed by composite deposition.
  subroutine deposit_amr_tree_spray_gas(species,tree,state,temperature,positions,contributions, &
      lower,spacing,periodic,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: positions(:,:),contributions(:,:),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate,recovered
    integer :: node,n(3)
    logical :: valid
    ok=.false.
    if (.not.state%is_valid(tree).or..not.temperature%is_valid(tree)) return
    if (size(state%fields(1)%values,1)/=reactive_nvar(size(species))) return
    if (size(temperature%fields(1)%values,1)/=1) return
    recovered=temperature
    do node=1,size(tree%nodes)
      n=tree%nodes(node)%shape
      call recover_reactive_temperatures_3d(species,state%fields(node)%values, &
        temperature%fields(node)%values(1,:,:,:),n(1),n(2),n(3),recovered%fields(node)%values(1,:,:,:),valid)
      if (.not.valid) return
    end do
    candidate=state
    call deposit_amr_tree_spray_field(tree,candidate,positions,contributions,lower,spacing,periodic,valid)
    if (.not.valid) return
    if (size(positions,2)==0) then
      ok=.true.;return
    end if
    do node=1,size(tree%nodes)
      n=tree%nodes(node)%shape
      call recover_reactive_temperatures_3d(species,candidate%fields(node)%values, &
        temperature%fields(node)%values(1,:,:,:),n(1),n(2),n(3),recovered%fields(node)%values(1,:,:,:),valid)
      if (.not.valid) return
    end do
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(recovered%fields,temperature%fields)
    ok=.true.
  end subroutine

  ! A rectangular cloud one containing-cell wide overlaps only exposed cells.
  ! On a uniform mesh this is the usual trilinear cloud-in-cell stencil.
  subroutine amr_tree_spray_stencil(tree,position,lower,spacing,periodic,cells,weights,volumes,ok,maximum_cells)
    type(amr_patch_tree_3d), intent(in) :: tree
    real(dp), intent(in) :: position(3),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    integer, allocatable, intent(inout) :: cells(:,:)
    real(dp), allocatable, intent(inout) :: weights(:),volumes(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_cells
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),location(:,:),candidate(:,:)
    real(dp), allocatable :: w(:),v(:)
    real(dp) :: length(3),center(3),width(3),clo(3),chi(3),lo(3),hi(3),d(3), &
      cell_lo(3),cell_hi(3),overlap(3),weight,total,expected,volume
    integer :: shifts(3),sx,sy,sz,pass,node,level,child,first(3),last(3),i,j,k,n,capacity,limit,ios,idx(3)
    logical :: valid,covered
    ok=.false.;limit=1000000
    if (present(maximum_cells)) limit=maximum_cells
    if (limit<1.or.limit>1000000) return
    call locate_amr_tree_spray_cells(tree,reshape(position,[3,1]),lower,spacing,periodic,location,valid)
    if (.not.valid) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if (.not.valid) return
    length=spacing*real(tree%nodes(1)%shape,dp);center=(position-lower)/length
    where(periodic) center=modulo(position-lower,length)/length
    level=tree%nodes(location(1,1))%level+1;width=1._dp/real(domain(:,level),dp)
    expected=product((min(1._dp,center+.5_dp*width)-max(0._dp,center-.5_dp*width))/width,mask=.not.periodic)
    shifts=merge(1,0,periodic);capacity=0
    do pass=1,2
      n=0;total=0._dp
      do sz=-shifts(3),shifts(3)
      do sy=-shifts(2),shifts(2)
      do sx=-shifts(1),shifts(1)
        clo=max(0._dp,center+real([sx,sy,sz],dp)-.5_dp*width)
        chi=min(1._dp,center+real([sx,sy,sz],dp)+.5_dp*width)
        if (any(chi<=clo)) cycle
        do node=1,size(tree%nodes)
          level=tree%nodes(node)%level+1;d=real(domain(:,level),dp)
          lo=max(clo,real(origin(:,node),dp)/d)
          hi=min(chi,real(origin(:,node)+int(tree%nodes(node)%shape,int64),dp)/d)
          if (any(hi<=lo)) cycle
          first=int(max(1_int64,min(int(tree%nodes(node)%shape,int64),floor(lo*d,kind=int64)-origin(:,node)+1)))
          last=int(max(1_int64,min(int(tree%nodes(node)%shape,int64),ceiling(hi*d,kind=int64)-origin(:,node))))
          volume=product(length/d)
          do k=first(3),last(3)
          do j=first(2),last(2)
          do i=first(1),last(1)
            idx=[i,j,k];covered=.false.
            do child=node+1,size(tree%nodes)
              if (tree%nodes(child)%parent/=node) cycle
              associate(p=>tree%nodes(child)%patch)
              if (all(idx>=[p%coarse_i_lower,p%coarse_j_lower,p%coarse_k_lower]).and. &
                  all(idx<=[p%coarse_i_upper,p%coarse_j_upper,p%coarse_k_upper])) covered=.true.
              end associate
              if (covered) exit
            end do
            if (covered) cycle
            cell_lo=real(origin(:,node)+int(idx-1,int64),dp)/d
            cell_hi=real(origin(:,node)+int(idx,int64),dp)/d
            overlap=max(0._dp,min(cell_hi,chi)-max(cell_lo,clo))/width
            weight=product(overlap)
            if (weight<=0._dp) cycle
            if (n==limit) return
            n=n+1;total=total+weight
            if (pass==2) then
              if (n>capacity) return
              candidate(:,n)=[node,idx];w(n)=weight;v(n)=volume
            end if
          end do
          end do
          end do
        end do
      end do
      end do
      end do
      if (n<1.or.total<=0._dp) return
      if (abs(total-expected)>1.e-10_dp*expected) return
      if (pass==1) then
        capacity=n
        allocate(candidate(4,n),w(n),v(n),stat=ios)
        if (ios/=0) return
      else
        if (n/=capacity) return
      end if
    end do
    w=w/total
    call move_alloc(candidate,cells);call move_alloc(w,weights);call move_alloc(v,volumes)
    ok=.true.
  end subroutine

  subroutine sample_amr_tree_spray_field(tree,state,positions,lower,spacing,periodic,samples,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: state
    real(dp), intent(in) :: positions(:,:),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    real(dp), allocatable, intent(inout) :: samples(:,:)
    logical, intent(out) :: ok
    integer, allocatable :: cells(:,:),lookup(:,:)
    real(dp), allocatable :: weights(:),volumes(:),candidate(:,:),term(:)
    integer :: p,q,c,node,i,j,k,nv,ios
    logical :: valid
    ok=.false.
    if (.not.state%is_valid(tree).or.size(positions,1)/=3) return
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,lookup,valid)
    if (.not.valid) return
    nv=size(state%fields(1)%values,1)
    allocate(candidate(nv,size(positions,2)),term(nv),stat=ios)
    if (ios/=0) return
    candidate=0._dp
    do p=1,size(positions,2)
      call amr_tree_spray_stencil(tree,positions(:,p),lower,spacing,periodic,cells,weights,volumes,valid)
      if (.not.valid) return
      do q=1,size(weights)
        node=cells(1,q);i=cells(2,q);j=cells(3,q);k=cells(4,q)
        term=weights(q)*state%fields(node)%values(:,i,j,k)
        do c=1,nv
          if (.not.can_add(candidate(c,p),term(c))) return
        end do
        candidate(:,p)=candidate(:,p)+term
      end do
    end do
    call move_alloc(candidate,samples);ok=.true.
  end subroutine

  ! Contributions are integrated quantities; state holds volume densities.
  ! Refresh hidden parents only after every exposed-cell update is accepted.
  subroutine deposit_amr_tree_spray_field(tree,state,positions,contributions,lower,spacing,periodic,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(inout) :: state
    real(dp), intent(in) :: positions(:,:),contributions(:,:),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    integer, allocatable :: cells(:,:),lookup(:,:)
    real(dp), allocatable :: weights(:),volumes(:),term(:)
    integer :: p,q,c,node,i,j,k,nv
    logical :: valid
    ok=.false.
    if (.not.state%is_valid(tree).or.size(positions,1)/=3) return
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,lookup,valid)
    if (.not.valid) return
    nv=size(state%fields(1)%values,1)
    if (any(shape(contributions)/=[nv,size(positions,2)])) return
    if (.not.all(ieee_is_finite(contributions))) return
    if (size(positions,2)==0) then
      ok=.true.;return
    end if
    candidate=state;allocate(term(nv))
    do p=1,size(positions,2)
      call amr_tree_spray_stencil(tree,positions(:,p),lower,spacing,periodic,cells,weights,volumes,valid)
      if (.not.valid) return
      do q=1,size(weights)
        node=cells(1,q);i=cells(2,q);j=cells(3,q);k=cells(4,q)
        term=weights(q)*contributions(:,p)
        if (any(abs(term)>huge(1._dp)*min(1._dp,volumes(q)))) return
        term=term/volumes(q)
        do c=1,nv
          if (.not.can_add(candidate%fields(node)%values(c,i,j,k),term(c))) return
        end do
        candidate%fields(node)%values(:,i,j,k)=candidate%fields(node)%values(:,i,j,k)+term
      end do
    end do
    call average_down_amr_patch_tree_3d(tree,candidate,valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields);ok=.true.
  end subroutine

  logical function can_add(a,b) result(valid)
    real(dp), intent(in) :: a,b
    valid=.false.
    if (b>0._dp) then
      if (a>huge(a)-b) return
    else if (b<0._dp) then
      if (a< -huge(a)-b) return
    end if
    valid=.true.
  end function
end module
