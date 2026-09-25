module spray_eb_tree_stencil_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,periodic_eb_faces_match
  use spray_amr_tree_coupling_mod, only: amr_tree_spray_stencil
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  implicit none
  private
  public :: eb_tree_spray_stencil,compatible_eb_nodal_lattice
contains
  pure logical function compatible_eb_nodal_lattice(dims,finest) result(ok)
    integer, intent(in) :: dims(3)
    integer(int64), intent(in) :: finest(3)
    integer(int64) :: factors(3)
    ok=.false.
    if(any(dims<1).or.any(dims>1000000).or.any(finest<1)) return
    if(product(int(dims,int64))>1000000_int64) return
    if(any(mod(int(dims,int64),finest)/=0)) return
    factors=int(dims,int64)/finest
    ok=all(factors>=1).and.all(factors==factors(1))
  end function
  ! Composite cloud weights, fluid-volume weighted, with visibility to fluid
  ! centroids. Geometry/phi must originate from the same finest nodal field.
  subroutine eb_tree_spray_stencil(tree,geometry,phi,position,lower,spacing,cells,weights,volumes,ok,periodic)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(in) :: phi(0:,0:,0:),position(3),lower(3),spacing(3)
    integer, allocatable, intent(inout) :: cells(:,:)
    real(dp), allocatable, intent(inout) :: weights(:),volumes(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: periodic(3)
    logical :: wrap(3)
    character(len=16) :: faces(6)
    real(dp) :: start(3),extent(3)
    integer :: axis
    integer, allocatable :: candidate(:,:),packed(:,:),ratios(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    real(dp), allocatable :: w(:),v(:),pw(:),pv(:)
    real(dp) :: upper(3),centroid(3),normal(3),point(3),fraction,f,total,h(3),lo(3),actual(3),tolerance(3)
    integer :: node,q,i,j,k,cell(3),n,ios,level
    logical :: valid,hit
    ok=.false.;wrap=.false.
    if(present(periodic)) wrap=periodic
    if(.not.tree%is_valid()) return
    if(size(geometry)/=size(tree%nodes)) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if(.not.valid) return
    if(.not.compatible_eb_nodal_lattice(shape(phi)-1,domain(:,size(domain,2)))) return
    call amr_tree_spray_stencil(tree,position,lower,spacing,wrap,candidate,w,v,valid)
    if(.not.valid) return
    upper=lower+spacing*real(tree%nodes(1)%shape,dp)
    extent=upper-lower;start=position
    faces='outflow'
    do axis=1,3
      if(.not.wrap(axis)) cycle
      faces(2*axis-1:2*axis)='periodic'
      start(axis)=lower(axis)+modulo(position(axis)-lower(axis),extent(axis))
    end do
    if(.not.periodic_eb_faces_match(geometry(1),faces)) return
    if(wrap(1)) then
      if(any(phi(0,:,:)/=phi(ubound(phi,1),:,:))) return
    end if
    if(wrap(2)) then
      if(any(phi(:,0,:)/=phi(:,ubound(phi,2),:))) return
    end if
    if(wrap(3)) then
      if(any(phi(:,:,0)/=phi(:,:,ubound(phi,3)))) return
    end if
    do node=1,size(geometry)
      if(.not.geometry(node)%is_valid()) return
      if(any([geometry(node)%nx,geometry(node)%ny,geometry(node)%nz]/=tree%nodes(node)%shape)) return
      level=tree%nodes(node)%level+1
      h=(upper-lower)/real(domain(:,level),dp)
      lo=lower+real(origin(:,node),dp)*h
      actual=[geometry(node)%dx,geometry(node)%dy,geometry(node)%dz]
      if(any(abs(actual-h)>64*epsilon(1._dp)*h)) return
      actual=[geometry(node)%x_lower,geometry(node)%y_lower,geometry(node)%z_lower]
      tolerance=64*epsilon(1._dp)*max(abs(lo),h)
      if(any(abs(actual-lo)>tolerance)) return
    end do
    hit=.false.;fraction=0;point=0;normal=0;cell=0
    call intersect_spray_level_set(phi,lower,upper,start,start,hit,fraction,point,normal,cell,valid)
    if(.not.valid) return
    do q=1,size(w)
      node=candidate(1,q);i=candidate(2,q);j=candidate(3,q);k=candidate(4,q)
      associate(g=>geometry(node))
        f=g%volume_fraction(i,j,k)
        w(q)=w(q)*f;v(q)=v(q)*f
        if(f<=0) cycle
        centroid=[g%x_lower+(real(i,dp)-.5_dp+g%cell_centroid_x(i,j,k))*g%dx, &
          g%y_lower+(real(j,dp)-.5_dp+g%cell_centroid_y(i,j,k))*g%dy, &
          g%z_lower+(real(k,dp)-.5_dp+g%cell_centroid_z(i,j,k))*g%dz]
      end associate
      call visible_centroid(centroid,hit,valid)
      if(.not.valid) return
      if(hit) w(q)=0
    end do
    total=sum(w);n=count(w>0)
    if(n<1.or.total<=0) return
    allocate(packed(4,n),pw(n),pv(n),stat=ios)
    if(ios/=0) return
    n=0
    do q=1,size(w)
      if(w(q)<=0) cycle
      n=n+1;packed(:,n)=candidate(:,q);pw(n)=w(q)/total;pv(n)=v(q)
    end do
    call move_alloc(packed,cells);call move_alloc(pw,weights);call move_alloc(pv,volumes)
    ok=.true.
  contains
    subroutine visible_centroid(target,blocked,accepted)
      real(dp), intent(in) :: target(3)
      logical, intent(out) :: blocked,accepted
      real(dp) :: a(3),b(3),delta(3),finish(3),alpha,crossing(3),limit
      integer :: segment,direction
      a=start;b=target
      do direction=1,3
        if(wrap(direction)) b(direction)=b(direction)- &
          anint((b(direction)-a(direction))/extent(direction))*extent(direction)
      end do
      accepted=.false.;blocked=.false.
      ! A nearest periodic image crosses each root face at most once.
      do segment=1,4
        delta=b-a;crossing=2;alpha=1
        do direction=1,3
          if(.not.wrap(direction)) cycle
          if(b(direction)>upper(direction).and.delta(direction)>0) then
            crossing(direction)=(upper(direction)-a(direction))/delta(direction)
          else if(b(direction)<lower(direction).and.delta(direction)<0) then
            crossing(direction)=(lower(direction)-a(direction))/delta(direction)
          end if
          alpha=min(alpha,crossing(direction))
        end do
        finish=max(lower,min(upper,a+alpha*delta))
        call intersect_spray_level_set(phi,lower,upper,a,finish,blocked,fraction,point,normal,cell,accepted)
        if(.not.accepted.or.blocked.or.alpha==1) return
        a=finish
        do direction=1,3
          if(crossing(direction)/=alpha) cycle
          limit=sign(extent(direction),delta(direction))
          a(direction)=a(direction)-limit;b(direction)=b(direction)-limit
        end do
      end do
      accepted=.false.
    end subroutine visible_centroid
  end subroutine
end module
