module spray_mpi_amr_film_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use spray_wall_impact_mod, only: spray_wall_film,valid_wall_component_inventory
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_integer_layout_3d
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  use spray_eb_tree_stencil_mod, only: compatible_eb_nodal_lattice
  use spray_amr_tree_partition_mod, only: locate_amr_tree_ghost_cells
  implicit none
  private
  public :: migrate_mpi_amr_wall_films,initialize_mpi_eb_film_levels,refresh_mpi_eb_film_copies
  type, public :: spray_level_film_population
    type(spray_wall_film), allocatable :: films(:)
    real(dp), allocatable :: component_mass(:,:),probes(:,:)
    integer, allocatable :: nodes(:)
  end type spray_level_film_population
contains
  ! Rebuild a level's virtual films and the next level's incoming ghost films.
  ! Selection follows logical physical levels; neither output changes the input.
  subroutine refresh_mpi_eb_film_copies(comm,tree,owners,physical,lower,spacing,target_level,grow,maximum_films, &
      nodal_phi,virtuals,next_ghosts,ok,periodic)
    logical, optional, intent(in) :: periodic(3)
    logical :: wrap(3)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),target_level,grow,maximum_films
    type(spray_level_film_population), intent(in) :: physical(:)
    real(dp), intent(in) :: lower(3),spacing(3),nodal_phi(0:,0:,0:)
    type(spray_level_film_population), intent(inout) :: virtuals,next_ghosts
    logical, intent(out) :: ok
    type(spray_level_film_population) :: coarse_copy,fine_copy
    integer, allocatable :: locations(:,:)
    integer :: level,levels,nc,p,n,cursor,rank,nranks,status,flags(5),low(5),high(5)
    integer(int64) :: count_all
    real(dp) :: probe_offset,upper(3)
    logical :: valid,agreed
    wrap=.false.
    if(present(periodic)) wrap=periodic
    ok=.false.;valid=tree%is_valid()
    call agree();if(.not.agreed) return
    levels=tree%level_count()
    valid=size(physical)==levels.and.target_level>=0.and.target_level<levels.and.grow>=0.and.grow<=1000000
    valid=valid.and.maximum_films>=0.and.maximum_films<=1000000.and.size(owners)==size(tree%nodes)
    call agree();if(.not.agreed) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,status)
    if(status/=MPI_SUCCESS) return
    flags=[levels,target_level,grow,maximum_films,size(owners)]
    call MPI_Allreduce(flags,low,5,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,5,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    nc=0;count_all=0
    do level=1,levels
      associate(group=>physical(level))
        valid=allocated(group%films).and.allocated(group%component_mass).and. &
          allocated(group%nodes).and.allocated(group%probes)
        if(.not.valid) exit
        if(level==1) nc=size(group%component_mass,1)
        n=size(group%films);count_all=count_all+int(n,int64)
        valid=nc>=1.and.nc<=16.and.all(shape(group%component_mass)==[nc,n]).and. &
          all(shape(group%probes)==[3,n]).and.size(group%nodes)==n
        if(.not.valid) exit
        valid=all(group%nodes>=1).and.all(group%nodes<=size(owners)).and.all(ieee_is_finite(group%probes))
        if(.not.valid) exit
        do p=1,n
          valid=valid_wall_component_inventory(group%films(p),group%component_mass(:,p)).and. &
            tree%nodes(group%nodes(p))%level==level-1.and.owners(group%nodes(p))==rank
          if(.not.valid) exit
        end do
      end associate
      if(.not.valid) exit
    end do
    valid=valid.and.count_all<=1000000
    call agree();if(.not.agreed) return
    n=0
    do level=target_level+2,levels
      n=n+size(physical(level)%films)
    end do
    allocate(coarse_copy%films(n),coarse_copy%component_mass(nc,n));cursor=0
    do level=target_level+2,levels
      n=size(physical(level)%films)
      coarse_copy%films(cursor+1:cursor+n)=physical(level)%films
      coarse_copy%component_mass(:,cursor+1:cursor+n)=physical(level)%component_mass
      cursor=cursor+n
    end do
    ! This collective validates replicated topology/geometry before ghost routing.
    call migrate_mpi_amr_wall_films(comm,tree,owners,coarse_copy%films,coarse_copy%component_mass,lower,spacing, &
      maximum_films,valid,nodal_phi,maximum_level=target_level,nodes=coarse_copy%nodes,probes=coarse_copy%probes,periodic=wrap)
    if(.not.valid) return
    upper=lower+spacing*real(tree%nodes(1)%shape,dp)
    probe_offset=max(1e-8_dp*minval((upper-lower)/real(shape(nodal_phi)-1,dp)), &
      128*epsilon(1._dp)*maxval(abs([lower,upper])))
    valid=.true.
    do level=1,levels
      do p=1,size(physical(level)%films)
        valid=valid.and.all(physical(level)%probes(:,p)== &
          physical(level)%films(p)%position+probe_offset*physical(level)%films(p)%normal)
      end do
    end do
    call agree();if(.not.agreed) return
    if(target_level+1<levels) then
      call locate_amr_tree_ghost_cells(tree,physical(target_level+1)%probes,lower,spacing, &
        wrap,target_level+1,grow,locations,valid)
      call agree();if(.not.agreed) return
      call route_copies(physical(target_level+1),locations(1,:),fine_copy)
      if(.not.agreed) return
    else
      allocate(fine_copy%films(0),fine_copy%component_mass(nc,0),fine_copy%nodes(0),fine_copy%probes(3,0))
    end if
    virtuals=coarse_copy;next_ghosts=fine_copy;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
    ! Input hierarchy and replicated metadata have already been validated.
    subroutine route_copies(input,destinations,output)
      type(spray_level_film_population), intent(in) :: input
      integer, intent(in) :: destinations(:)
      type(spray_level_film_population), intent(out) :: output
      integer, allocatable :: counts(:),offsets(:),incoming(:),starts(:),next(:),sent_meta(:,:),received_meta(:,:)
      real(dp), allocatable :: sent(:,:),received(:,:)
      integer :: p,r,q,width,received_count
      allocate(counts(nranks),offsets(nranks),incoming(nranks),starts(nranks),next(nranks));counts=0
      do p=1,size(destinations)
        if(destinations(p)==0) cycle
        r=owners(destinations(p))+1;counts(r)=counts(r)+1
      end do
      call MPI_Alltoall(counts,1,MPI_INTEGER,incoming,1,MPI_INTEGER,comm,status)
      if(status/=MPI_SUCCESS) then
        agreed=.false.;return
      end if
      valid=sum(int(incoming,int64))<=maximum_films
      call agree();if(.not.agreed) return
      offsets(1)=0;starts(1)=0
      do r=2,nranks
        offsets(r)=offsets(r-1)+counts(r-1);starts(r)=starts(r-1)+incoming(r-1)
      end do
      next=offsets;received_count=sum(incoming);width=13+nc
      allocate(sent(width,sum(counts)),received(width,received_count), &
        sent_meta(2,sum(counts)),received_meta(2,received_count))
      do p=1,size(destinations)
        if(destinations(p)==0) cycle
        r=owners(destinations(p))+1;next(r)=next(r)+1;q=next(r)
        associate(f=>input%films(p))
          sent(:,q)=[f%position,f%normal,f%mass,f%temperature,f%diameter,f%height, &
            input%probes(:,p),input%component_mass(:,p)]
          sent_meta(:,q)=[f%material_id,destinations(p)]
        end associate
      end do
      call MPI_Alltoallv(sent,width*counts,width*offsets,MPI_DOUBLE_PRECISION,received,width*incoming, &
        width*starts,MPI_DOUBLE_PRECISION,comm,status)
      if(status/=MPI_SUCCESS) then
        agreed=.false.;return
      end if
      call MPI_Alltoallv(sent_meta,2*counts,2*offsets,MPI_INTEGER,received_meta,2*incoming,2*starts,MPI_INTEGER,comm,status)
      if(status/=MPI_SUCCESS) then
        agreed=.false.;return
      end if
      allocate(output%films(received_count),output%component_mass(nc,received_count), &
        output%probes(3,received_count),output%nodes(received_count))
      valid=.true.
      do p=1,received_count
        output%films(p)=spray_wall_film(received(1:3,p),received(4:6,p),received(7,p),received(8,p), &
          received(9,p),received(10,p),received_meta(1,p))
        output%probes(:,p)=received(11:13,p);output%component_mass(:,p)=received(14:,p)
        output%nodes(p)=received_meta(2,p)
        valid=valid_wall_component_inventory(output%films(p),output%component_mass(:,p)).and. &
          owners(output%nodes(p))==rank.and.tree%nodes(output%nodes(p))%level==target_level+1
        if(.not.valid) exit
      end do
      call agree()
    end subroutine route_copies
  end subroutine refresh_mpi_eb_film_copies

  ! Build physical level populations without modifying the flat input films.
  subroutine initialize_mpi_eb_film_levels(comm,tree,owners,films,component_mass,lower,spacing,maximum_films, &
      nodal_phi,physical,ok,periodic)
    logical, optional, intent(in) :: periodic(3)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),maximum_films
    type(spray_wall_film), allocatable, intent(in) :: films(:)
    real(dp), allocatable, intent(in) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),nodal_phi(0:,0:,0:)
    type(spray_level_film_population), allocatable, intent(inout) :: physical(:)
    logical, intent(out) :: ok
    type(spray_wall_film), allocatable :: candidate(:)
    real(dp), allocatable :: mass(:,:),probes(:,:)
    integer, allocatable :: nodes(:)
    logical, allocatable :: selected(:)
    type(spray_level_film_population), allocatable :: levels(:)
    integer :: level,p,c,status
    logical :: valid,agreed
    ok=.false.;valid=allocated(films).and.allocated(component_mass)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    candidate=films;mass=component_mass
    call migrate_mpi_amr_wall_films(comm,tree,owners,candidate,mass,lower,spacing,maximum_films,valid, &
      nodal_phi,nodes=nodes,probes=probes,periodic=periodic)
    if(.not.valid) return
    allocate(levels(tree%level_count()),selected(size(candidate)))
    do level=1,size(levels)
      do p=1,size(candidate)
        selected(p)=tree%nodes(nodes(p))%level==level-1
      end do
      levels(level)%films=pack(candidate,selected);levels(level)%nodes=pack(nodes,selected)
      allocate(levels(level)%component_mass(size(mass,1),count(selected)),levels(level)%probes(3,count(selected)))
      do c=1,size(mass,1)
        levels(level)%component_mass(c,:)=pack(mass(c,:),selected)
      end do
      do c=1,3
        levels(level)%probes(c,:)=pack(probes(c,:),selected)
      end do
    end do
    call move_alloc(levels,physical);ok=.true.
  end subroutine initialize_mpi_eb_film_levels

  ! Stationary films use the finest containing cell, with closed physical bounds.
  ! Optional maximum_level caps ownership without changing the wall geometry.
  ! Paired logical_level/grow instead selects that level's valid or grown box;
  ! films outside all such boxes reject, and physical bounds remain enforced.
  ! Preserve exact wall coordinates and every film/component field. Films have
  ! no global ID; receive order is source rank followed by original local order.
  subroutine migrate_mpi_amr_wall_films(comm,tree,node_owners,films,component_mass,lower,spacing,maximum_films,ok, &
      nodal_phi,maximum_level,nodes,probes,logical_level,grow,periodic)
    logical, optional, intent(in) :: periodic(3)
    logical :: wrap(3)
    integer, allocatable, optional, intent(inout) :: nodes(:)
    real(dp), allocatable, optional, intent(inout) :: probes(:,:)
    integer, allocatable :: next_nodes(:)
    real(dp), allocatable :: next_probes(:,:)
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: node_owners(:),maximum_films
    integer, optional, intent(in) :: maximum_level,logical_level,grow
    logical, intent(out) :: ok
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: node,node_count,descriptor(13),reference(13),low_nodes,high_nodes,allocation_status
    type(spray_wall_film), allocatable :: candidate(:)
    real(dp), allocatable :: sent(:,:),received(:,:),masses(:,:)
    integer, allocatable :: counts(:),offsets(:),incoming(:),starts(:),cursor(:),owners(:),material(:),got_material(:)
    real(dp) :: bounds(6),lowbounds(6),highbounds(6),upper(3)
    integer :: config(10),lowconfig(10),highconfig(10),top_level,selected_level,extent,nranks,ierr,nc,np,p,r,q,nreceived,width
    logical :: valid,agreed
    integer :: dims(3),reference_dims(3)
    real(dp), allocatable :: reference_phi(:,:,:)
    real(dp) :: probe_offset
    ok=.false.;nc=0;wrap=.false.
    if(present(periodic)) wrap=periodic
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=allocated(films).and.allocated(component_mass).and.(present(nodes).eqv.present(probes))
    valid=valid.and.(present(logical_level).eqv.present(grow))
    valid=valid.and..not.(present(logical_level).and.present(maximum_level))
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
    selected_level=-1;extent=0
    if(present(logical_level)) then
      selected_level=logical_level;extent=grow
    end if
    valid=top_level>=0.and.top_level<=maxval(tree%nodes%level)
    if(present(logical_level)) valid=valid.and.selected_level>=0.and.selected_level<=top_level.and.extent>=0.and.extent<=1000000
    call agree();if(.not.agreed) return
    config=[nc,maximum_films,merge(1,0,present(nodal_phi)),top_level, &
      merge(1,0,present(nodes)),selected_level,extent,merge(1,0,wrap)]
    bounds=[lower,spacing]
    call MPI_Allreduce(config,lowconfig,10,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,highconfig,10,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,lowbounds,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,highbounds,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(lowconfig/=highconfig).or.any(lowbounds/=highbounds)) return
    node_count=size(tree%nodes)
    call MPI_Allreduce(node_count,low_nodes,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(node_count,high_nodes,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_nodes/=high_nodes) return
    do node=1,node_count
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
      if(wrap(1)) valid=valid.and.all(nodal_phi(0,:,:)==nodal_phi(dims(1),:,:))
      if(wrap(2)) valid=valid.and.all(nodal_phi(:,0,:)==nodal_phi(:,dims(2),:))
      if(wrap(3)) valid=valid.and.all(nodal_phi(:,:,0)==nodal_phi(:,:,dims(3)))
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
    if(present(nodes)) then
      allocate(next_nodes(nreceived),next_probes(3,nreceived),stat=allocation_status)
      valid=allocation_status==0;call agree();if(.not.agreed) return
    end if
    do p=1,nreceived
      candidate(p)=spray_wall_film(received(1:3,p),received(4:6,p),received(7,p),received(8,p), &
        received(9,p),received(10,p),got_material(p))
      masses(:,p)=received(11:,p)
      valid=valid_wall_component_inventory(candidate(p),masses(:,p))
      if(present(nodes)) then
        valid=valid.and.film_owner(candidate(p),next_nodes(p),next_probes(:,p))==r
      else
        valid=valid.and.film_owner(candidate(p))==r
      end if
      if (.not.valid) exit
    end do
    call agree()
    if (.not.agreed) return
    if(present(nodes)) then
      call move_alloc(next_nodes,nodes);call move_alloc(next_probes,probes)
    end if
    call move_alloc(candidate,films);call move_alloc(masses,component_mass);ok=.true.
  contains
    integer function film_owner(film,found_node,fluid_probe) result(owner_rank)
      integer, optional, intent(out) :: found_node
      real(dp), optional, intent(out) :: fluid_probe(3)
      type(spray_wall_film), intent(in) :: film
      real(dp) :: probe(3),point(3),normal(3),fraction
      integer :: cell(3)
      logical :: hit,good
      owner_rank=-1
      if(present(found_node)) found_node=0
      if(present(fluid_probe)) fluid_probe=0
      if(any(film%position<lower).or.any(film%position>upper)) return
      if (present(nodal_phi)) then
        probe=film%position+probe_offset*film%normal
        point=0;normal=0;fraction=0;cell=0;hit=.false.
        call intersect_film_probe(nodal_phi,lower,upper,probe,film%position-probe_offset*film%normal, &
          wrap,hit,fraction,point,normal,cell,good)
        if (.not.good.or..not.hit) return
        if (norm2(point-film%position)>.01_dp*probe_offset) return
        if (dot_product(normal,film%normal)<1-1e-10_dp) return
        owner_rank=owner(probe,found_node)
        if(present(fluid_probe)) fluid_probe=probe
      else
        owner_rank=owner(film%position,found_node)
        if(present(fluid_probe)) fluid_probe=film%position
      end if
    end function
    integer function owner(position,found_node) result(owner_rank)
      integer, optional, intent(out) :: found_node
      real(dp), intent(in) :: position(3)
      real(dp) :: fraction(3),offset(3),length(3)
      integer(int64) :: cell(3)
      integer :: peer,level,best,axis
      integer, allocatable :: location(:,:)
      logical :: located
      owner_rank=-1
      if(present(found_node)) found_node=0
      if (.not.all(ieee_is_finite(position))) return
      if (any((position<lower.or.position>upper).and..not.wrap).or.any(abs(position)>1.e9_dp)) return
      if(selected_level==0) then
        owner_rank=node_owners(1)
        if(present(found_node)) found_node=1
        return
      else if(selected_level>0) then
        ! Keep the generating logical level until synchronization, including
        ! deposits in its grown region. A miss rejects rather than losing film.
        call locate_amr_tree_ghost_cells(tree,reshape(position,[3,1]),lower,spacing, &
          wrap,selected_level,extent,location,located)
        if(.not.located) return
        peer=location(1,1)
        if(peer==0) return
        owner_rank=node_owners(peer)
        if(present(found_node)) found_node=peer
        return
      end if
      length=spacing*real(tree%nodes(1)%shape,dp);offset=position-lower
      do axis=1,3
        if(wrap(axis)) offset(axis)=modulo(offset(axis),length(axis))
      end do
      fraction=offset/length;best=-1
      do peer=1,size(tree%nodes)
        level=tree%nodes(peer)%level+1
        if (level<=best.or.level>top_level+1) cycle
        ! Closed physical bounds: an exact upper wall belongs to the last cell.
        cell=min(domain(:,level)-1,floor(fraction*real(domain(:,level),dp),kind=int64))
        if (any(cell<origin(:,peer)).or.any(cell-origin(:,peer)>=int(tree%nodes(peer)%shape,int64))) cycle
        best=level;owner_rank=node_owners(peer)
        if(present(found_node)) found_node=peer
      end do
    end function
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
  ! Split the short normal probe at periodic faces. Each piece uses the same
  ! tetrahedral intersection as interior films; wall coordinates stay unwrapped.
  subroutine intersect_film_probe(phi,lower,upper,start,finish,wrap,hit,fraction,point,normal,cell,ok)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),upper(3),start(3),finish(3)
    logical, intent(in) :: wrap(3)
    logical, intent(out) :: hit,ok
    real(dp), intent(out) :: fraction,point(3),normal(3)
    integer, intent(out) :: cell(3)
    real(dp) :: a(3),b(3),end(3),delta(3),length(3),crossing(3),alpha,shift,offset(3),local_fraction
    integer :: axis,segment
    ok=.false.;hit=.false.;fraction=0;point=0;normal=0;cell=0
    if(.not.any(wrap)) then
      call intersect_spray_level_set(phi,lower,upper,start,finish,hit,fraction,point,normal,cell,ok)
      return
    end if
    length=upper-lower
    if(any(abs(finish-start)>length)) return
    if(any((start<lower.or.start>upper.or.finish<lower.or.finish>upper).and..not.wrap)) return
    a=start
    do axis=1,3
      if(wrap(axis)) a(axis)=lower(axis)+modulo(start(axis)-lower(axis),length(axis))
    end do
    offset=start-a;b=finish-offset
    ! A probe shorter than one domain length crosses at most three faces.
    do segment=1,4
      delta=b-a;crossing=2;alpha=1
      do axis=1,3
        if(.not.wrap(axis)) cycle
        if(b(axis)>upper(axis).and.delta(axis)>0) crossing(axis)=(upper(axis)-a(axis))/delta(axis)
        if(b(axis)<lower(axis).and.delta(axis)<0) crossing(axis)=(lower(axis)-a(axis))/delta(axis)
        alpha=min(alpha,crossing(axis))
      end do
      end=max(lower,min(upper,a+alpha*delta))
      call intersect_spray_level_set(phi,lower,upper,a,end,hit,local_fraction,point,normal,cell,ok)
      if(.not.ok) return
      if(hit) then
        point=point+offset
        axis=maxloc(abs(finish-start),dim=1)
        if(finish(axis)/=start(axis)) fraction=(point(axis)-start(axis))/(finish(axis)-start(axis))
        return
      end if
      if(alpha==1) then
        fraction=1;point=finish
        return
      end if
      a=end
      do axis=1,3
        if(crossing(axis)/=alpha) cycle
        shift=sign(length(axis),delta(axis))
        a(axis)=a(axis)-shift;b(axis)=b(axis)-shift;offset(axis)=offset(axis)+shift
      end do
    end do
    ok=.false.
  end subroutine intersect_film_probe
end module
