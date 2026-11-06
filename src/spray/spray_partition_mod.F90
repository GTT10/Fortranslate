! Deterministic Cartesian ownership and send-buffer ordering for particle migration.
module spray_partition_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use spray_parcel_mod, only: spray_parcel,valid_spray_parcel
  use spray_mixture_mod, only: valid_component_inventory
  implicit none
  private
  public :: spray_cartesian_owner,partition_spray_parcels,spray_cartesian_block
contains
  subroutine spray_cartesian_block(cells,parts,rank,first,counts,ok)
    integer, intent(in) :: cells(3),parts(3),rank
    integer, intent(out) :: first(3),counts(3)
    logical, intent(out) :: ok
    integer :: coordinate(3),a,base,extra
    first=0;counts=0;ok=.false.
    if (spray_cartesian_owner([0.0_dp,0.0_dp,0.0_dp],[0.0_dp,0.0_dp,0.0_dp], &
      [1.0_dp,1.0_dp,1.0_dp],cells,parts,[.false.,.false.,.false.])<0) return
    if (rank<0.or.rank>=product(parts)) return
    coordinate=[modulo(rank,parts(1)),modulo(rank/parts(1),parts(2)),rank/(parts(1)*parts(2))]
    do a=1,3
      base=cells(a)/parts(a);extra=modulo(cells(a),parts(a))
      counts(a)=base+merge(1,0,coordinate(a)<extra)
      first(a)=coordinate(a)*base+min(coordinate(a),extra)
    end do
    ok=.true.
  end subroutine

  integer function spray_cartesian_owner(position,lower,spacing,cells,parts,periodic) result(rank)
    real(dp), intent(in) :: position(3),lower(3),spacing(3)
    integer, intent(in) :: cells(3),parts(3)
    logical, intent(in) :: periodic(3)
    real(dp) :: coordinate,length
    integer :: axis,cell,base,extra,cut,owner(3)
    rank=-1
    if (.not.all(ieee_is_finite([position,lower,spacing]))) return
    if (maxval(abs(position))>1e9_dp.or.maxval(abs(lower))>1e9_dp) return
    if (any(spacing<1e-12_dp).or.any(spacing>1e9_dp)) return
    if (any(cells<1).or.any(cells>100000000).or.any(parts<1).or.any(parts>cells)) return
    if (any(parts>1000000)) return
    if (product(int(parts,int64))>1000000_int64) return
    do axis=1,3
      length=spacing(axis)*cells(axis);coordinate=position(axis)-lower(axis)
      if (periodic(axis)) then
        coordinate=modulo(coordinate,length)
      else
        if (coordinate<0.or.coordinate>=length) return
      end if
      cell=min(cells(axis)-1,int(coordinate/spacing(axis)))
      base=cells(axis)/parts(axis);extra=modulo(cells(axis),parts(axis));cut=(base+1)*extra
      if (cell<cut) then
        owner(axis)=cell/(base+1)
      else
        owner(axis)=extra+(cell-cut)/base
      end if
    end do
    rank=owner(1)+parts(1)*(owner(2)+parts(2)*owner(3))
  end function

  subroutine partition_spray_parcels(parcels,lower,spacing,cells,parts,periodic,counts,offsets,ok, &
      component_mass,tab_tangent_history)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(in) :: lower(3),spacing(3)
    integer, intent(in) :: cells(3),parts(3)
    logical, intent(in) :: periodic(3)
    integer, allocatable, intent(inout) :: counts(:),offsets(:)
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    type(spray_parcel), allocatable :: ordered(:)
    real(dp), allocatable, optional, intent(inout) :: tab_tangent_history(:)
    real(dp), allocatable :: history(:)
    real(dp), allocatable :: masses(:,:)
    integer, allocatable :: owner(:),bins(:),starts(:),next(:)
    integer :: p,r,n,destination
    ok=.false.
    if (.not.allocated(parcels)) return
    if (spray_cartesian_owner(lower,lower,spacing,cells,parts,periodic)<0) return
    n=product(parts)
    if (present(component_mass)) then
      if (.not.allocated(component_mass)) return
      if (size(component_mass,2)/=size(parcels)) return
      if (size(component_mass,1)<1.or.size(component_mass,1)>16) return
    end if
    if(present(tab_tangent_history)) then
      if(.not.allocated(tab_tangent_history)) return
      if(size(tab_tangent_history)/=size(parcels)) return
      if(.not.all(ieee_is_finite(tab_tangent_history)).or.any(tab_tangent_history<0)) return
    endif
    allocate(owner(size(parcels)),bins(n),starts(n),next(n),ordered(size(parcels)));bins=0
    do p=1,size(parcels)
      if (.not.valid_spray_parcel(parcels(p))) return
      if (p>1) then
        if (parcels(p)%id<=parcels(p-1)%id) return
      end if
      if (present(component_mass)) then
        if (.not.valid_component_inventory(parcels(p),component_mass(:,p))) return
      end if
      r=spray_cartesian_owner(parcels(p)%position,lower,spacing,cells,parts,periodic)
      if (r<0) return
      owner(p)=r+1;bins(r+1)=bins(r+1)+1
    end do
    starts(1)=0
    do r=2,n
      starts(r)=starts(r-1)+bins(r-1)
    end do
    next=starts
    if (present(component_mass)) allocate(masses(size(component_mass,1),size(parcels)))
    if(present(tab_tangent_history)) allocate(history(size(parcels)))
    do p=1,size(parcels)
      r=owner(p);destination=next(r)+1;next(r)=destination
      ordered(destination)=parcels(p)
      if(present(tab_tangent_history)) history(destination)=tab_tangent_history(p)
      if (present(component_mass)) masses(:,destination)=component_mass(:,p)
    end do
    call move_alloc(ordered,parcels);call move_alloc(bins,counts);call move_alloc(starts,offsets)
    if (present(component_mass)) call move_alloc(masses,component_mass)
    if(present(tab_tangent_history)) call move_alloc(history,tab_tangent_history)
    ok=.true.
  end subroutine
end module
