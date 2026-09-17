module spray_eb_tree_stencil_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_amr_tree_coupling_mod, only: amr_tree_spray_stencil
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  implicit none
  private
  public :: eb_tree_spray_stencil
contains
  ! Composite cloud weights, fluid-volume weighted, with visibility to fluid
  ! centroids. Geometry/phi must originate from the same finest nodal field.
  subroutine eb_tree_spray_stencil(tree,geometry,phi,position,lower,spacing,cells,weights,volumes,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(in) :: phi(0:,0:,0:),position(3),lower(3),spacing(3)
    integer, allocatable, intent(inout) :: cells(:,:)
    real(dp), allocatable, intent(inout) :: weights(:),volumes(:)
    logical, intent(out) :: ok
    integer, allocatable :: candidate(:,:),packed(:,:),ratios(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    real(dp), allocatable :: w(:),v(:),pw(:),pv(:)
    real(dp) :: upper(3),centroid(3),normal(3),point(3),fraction,f,total,h(3),lo(3),actual(3),tolerance(3)
    integer :: node,q,i,j,k,cell(3),n,ios,level
    logical :: valid,hit
    ok=.false.
    if(.not.tree%is_valid()) return
    if(size(geometry)/=size(tree%nodes)) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if(.not.valid) return
    if(any(int(shape(phi)-1,int64)/=domain(:,size(domain,2)))) return
    call amr_tree_spray_stencil(tree,position,lower,spacing,[.false.,.false.,.false.],candidate,w,v,valid)
    if(.not.valid) return
    upper=lower+spacing*real(tree%nodes(1)%shape,dp)
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
    call intersect_spray_level_set(phi,lower,upper,position,position,hit,fraction,point,normal,cell,valid)
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
      call intersect_spray_level_set(phi,lower,upper,position,centroid,hit,fraction,point,normal,cell,valid)
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
  end subroutine
end module
