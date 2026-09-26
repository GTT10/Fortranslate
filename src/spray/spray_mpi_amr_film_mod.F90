module spray_mpi_amr_film_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use spray_wall_impact_mod, only: spray_wall_film,valid_wall_component_inventory
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_integer_layout_3d
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  use spray_eb_tree_stencil_mod, only: compatible_eb_nodal_lattice
  implicit none
  private
  public :: migrate_mpi_amr_wall_films
contains
  ! Stationary films use the finest containing cell, with closed physical bounds.
  ! Optional maximum_level caps ownership without changing the wall geometry.
  ! Preserve exact wall coordinates and every film/component field. Films have
  ! no global ID; receive order is source rank followed by original local order.
  subroutine migrate_mpi_amr_wall_films(comm,tree,node_owners,films,component_mass,lower,spacing,maximum_films,ok, &
      nodal_phi,maximum_level)
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: node_owners(:),maximum_films
    integer, optional, intent(in) :: maximum_level
    logical, intent(out) :: ok
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: node,nodes,descriptor(13),reference(13),low_nodes,high_nodes,allocation_status
    type(spray_wall_film), allocatable :: candidate(:)
    real(dp), allocatable :: sent(:,:),received(:,:),masses(:,:)
    integer, allocatable :: counts(:),offsets(:),incoming(:),starts(:),cursor(:),owners(:),material(:),got_material(:)
    real(dp) :: bounds(6),lowbounds(6),highbounds(6),upper(3)
    integer :: config(4),lowconfig(4),highconfig(4),top_level,nranks,ierr,nc,np,p,r,q,nreceived,width
    logical :: valid,agreed
    integer :: dims(3),reference_dims(3)
    real(dp), allocatable :: reference_phi(:,:,:)
    real(dp) :: probe_offset
    ok=.false.;nc=0
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=allocated(films).and.allocated(component_mass)
    if (allocated(component_mass)) nc=size(component_mass,1)
    call agree()
    if (.not.agreed) return
    np=size(films)
    valid=nc>=1.and.nc<=16.and.size(component_mass,2)==np.and.np<=1000000
    valid=valid.and.tree%is_valid()
    if (valid) valid=size(node_owners)==size(tree%nodes)
    if (valid) valid=all(node_owners>=0).and.all(node_owners<nranks)
    valid=valid.and.maximum_films>=0.and.maximum_films<=1000000
    valid=valid.and.all(ieee_is_finite([lower,spacing]))
    if (valid) valid=all(abs(lower)<=1.e9_dp).and.all(spacing>=1.e-12_dp).and.all(spacing<=1.e9_dp)
    call agree()
    if (.not.agreed) return
    top_level=maxval(tree%nodes%level)
    if(present(maximum_level)) top_level=maximum_level
    valid=top_level>=0.and.top_level<=maxval(tree%nodes%level)
    call agree();if(.not.agreed) return
    config=[nc,maximum_films,merge(1,0,present(nodal_phi)),top_level];bounds=[lower,spacing]
    call MPI_Allreduce(config,lowconfig,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,highconfig,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,lowbounds,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,highbounds,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(lowconfig/=highconfig).or.any(lowbounds/=highbounds)) return
    nodes=size(tree%nodes)
    call MPI_Allreduce(nodes,low_nodes,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nodes,high_nodes,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_nodes/=high_nodes) return
    do node=1,nodes
      associate(n=>tree%nodes(node),patch=>tree%nodes(node)%patch)
        descriptor=[node_owners(node),n%parent,n%level,n%shape,patch%refinement_ratio, &
          patch%coarse_i_lower,patch%coarse_i_upper,patch%coarse_j_lower,patch%coarse_j_upper, &
          patch%coarse_k_lower,patch%coarse_k_upper]
      end associate
      reference=descriptor
      call MPI_Bcast(reference,13,MPI_INTEGER,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      valid=valid.and.all(reference==descriptor)
    end do
    call agree();if (.not.agreed) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if (valid) valid=all(domain<=2_int64**52)
    call agree();if (.not.agreed) return
    if (present(nodal_phi)) then
      dims=shape(nodal_phi)-1
      valid=compatible_eb_nodal_lattice(dims,domain(:,size(ratios)))
      valid=valid.and.all(ieee_is_finite(nodal_phi))
      if (valid) valid=maxval(abs(nodal_phi))<=1e100_dp
      call agree();if (.not.agreed) return
      reference_dims=dims
      call MPI_Bcast(reference_dims,3,MPI_INTEGER,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      valid=all(dims==reference_dims)
      call agree();if (.not.agreed) return
      allocate(reference_phi(0:dims(1),0:dims(2),0:dims(3)),stat=allocation_status)
      valid=allocation_status==0;call agree();if (.not.agreed) return
      reference_phi=nodal_phi
      call MPI_Bcast(reference_phi,size(reference_phi),MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      valid=all(reference_phi==nodal_phi)
      call agree();if (.not.agreed) return
      deallocate(reference_phi)
    end if
    upper=lower+spacing*real(tree%nodes(1)%shape,dp)
    valid=all(upper>lower)
    call agree();if (.not.agreed) return
    if (present(nodal_phi)) then
      probe_offset=max(1e-8_dp*minval((upper-lower)/real(dims,dp)), &
        128*epsilon(1._dp)*maxval(abs([lower,upper])))
    end if
    allocate(counts(nranks),offsets(nranks),incoming(nranks),starts(nranks),cursor(nranks),owners(np),stat=allocation_status)
    valid=allocation_status==0;call agree();if (.not.agreed) return
    counts=0
    do p=1,np
      valid=valid_wall_component_inventory(films(p),component_mass(:,p))
      if (.not.valid) exit
      r=film_owner(films(p));valid=r>=0
      if (.not.valid) exit
      owners(p)=r+1;counts(r+1)=counts(r+1)+1
    end do
    call agree()
    if (.not.agreed) return
    call MPI_Alltoall(counts,1,MPI_INTEGER,incoming,1,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=sum(int(incoming,int64))<=maximum_films
    call agree()
    if (.not.agreed) return
    offsets(1)=0;starts(1)=0
    do r=2,nranks
      offsets(r)=offsets(r-1)+counts(r-1);starts(r)=starts(r-1)+incoming(r-1)
    end do
    cursor=offsets;nreceived=sum(incoming);width=10+nc
    allocate(sent(width,np),received(width,nreceived),material(np),got_material(nreceived),stat=allocation_status)
    valid=allocation_status==0;call agree();if (.not.agreed) return
    do p=1,np
      r=owners(p);cursor(r)=cursor(r)+1;q=cursor(r)
      sent(:,q)=[films(p)%position,films(p)%normal,films(p)%mass,films(p)%temperature, &
        films(p)%diameter,films(p)%height,component_mass(:,p)]
      material(q)=films(p)%material_id
    end do
    call MPI_Alltoallv(sent,width*counts,width*offsets,MPI_DOUBLE_PRECISION,received,width*incoming, &
      width*starts,MPI_DOUBLE_PRECISION,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Alltoallv(material,counts,offsets,MPI_INTEGER,got_material,incoming,starts,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    allocate(candidate(nreceived),masses(nc,nreceived),stat=allocation_status)
    valid=allocation_status==0;call agree();if (.not.agreed) return
    call MPI_Comm_rank(comm,r,ierr)
    if (ierr/=MPI_SUCCESS) return
    do p=1,nreceived
      candidate(p)=spray_wall_film(received(1:3,p),received(4:6,p),received(7,p),received(8,p), &
        received(9,p),received(10,p),got_material(p))
      masses(:,p)=received(11:,p)
      valid=valid_wall_component_inventory(candidate(p),masses(:,p)).and.film_owner(candidate(p))==r
      if (.not.valid) exit
    end do
    call agree()
    if (.not.agreed) return
    call move_alloc(candidate,films);call move_alloc(masses,component_mass);ok=.true.
  contains
    integer function film_owner(film) result(owner_rank)
      type(spray_wall_film), intent(in) :: film
      real(dp) :: probe(3),point(3),normal(3),fraction
      integer :: cell(3)
      logical :: hit,good
      owner_rank=-1
      if (present(nodal_phi)) then
        probe=film%position+probe_offset*film%normal
        point=0;normal=0;fraction=0;cell=0;hit=.false.
        call intersect_spray_level_set(nodal_phi,lower,upper,probe,film%position-probe_offset*film%normal, &
          hit,fraction,point,normal,cell,good)
        if (.not.good.or..not.hit) return
        if (norm2(point-film%position)>.01_dp*probe_offset) return
        if (dot_product(normal,film%normal)<1-1e-10_dp) return
        owner_rank=owner(probe)
      else
        owner_rank=owner(film%position)
      end if
    end function
    integer function owner(position) result(owner_rank)
      real(dp), intent(in) :: position(3)
      real(dp) :: fraction(3)
      integer(int64) :: cell(3)
      integer :: peer,level,best
      owner_rank=-1
      if (.not.all(ieee_is_finite(position))) return
      if (any(position<lower).or.any(position>upper).or.any(abs(position)>1.e9_dp)) return
      fraction=(position-lower)/(spacing*real(tree%nodes(1)%shape,dp));best=-1
      do peer=1,size(tree%nodes)
        level=tree%nodes(peer)%level+1
        if (level<=best.or.level>top_level+1) cycle
        ! Closed physical bounds: an exact upper wall belongs to the last cell.
        cell=min(domain(:,level)-1,floor(fraction*real(domain(:,level),dp),kind=int64))
        if (any(cell<origin(:,peer)).or.any(cell-origin(:,peer)>=int(tree%nodes(peer)%shape,int64))) cycle
        best=level;owner_rank=node_owners(peer)
      end do
    end function
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
end module
