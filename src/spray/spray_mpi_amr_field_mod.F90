module spray_mpi_amr_field_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d,fetch_mpi_amr_tree_cells_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use spray_amr_tree_partition_mod, only: locate_amr_tree_spray_cells
  use spray_amr_tree_coupling_mod, only: amr_tree_spray_stencil
  use eb_geometry_3d_mod, only: eb_geometry_3d,periodic_eb_faces_match
  use spray_eb_tree_stencil_mod, only: eb_tree_spray_stencil,compatible_eb_nodal_lattice
  implicit none
  private
  type :: particle_stencil
    integer, allocatable :: cells(:,:)
    real(dp), allocatable :: weights(:),volumes(:)
  end type
  public :: sample_mpi_amr_tree_spray_field,deposit_mpi_amr_tree_spray_field
  public :: deposit_mpi_amr_tree_spray_gas
contains
  ! Borrow the identical ordered species table on every rank. U/T publish only
  ! after owner-routed deposition, parent restriction and collective EOS success.
  ! Optional gas_change is the deposited increment per fluid volume, accumulated
  ! independently of the gas background and published with U/T.
  subroutine deposit_mpi_amr_tree_spray_gas(comm,species,tree,owners,state,temperature,positions,contributions, &
      lower,spacing,periodic,ok,maximum_requests,maximum_received,eb_geometry,nodal_phi,gas_change)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: positions(:,:),contributions(:,:),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_requests,maximum_received
    type(eb_geometry_3d), optional, intent(in) :: eb_geometry(:)
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(amr_tree_state_3d) :: candidate,recovered
    type(amr_tree_state_3d), optional, intent(inout) :: gas_change
    type(amr_tree_state_3d), allocatable :: change
    integer :: ierr,flag,low,high
    logical :: valid,has_particles,any_particles
    ok=.false.;recovered=temperature
    flag=merge(1,0,present(eb_geometry))+2*merge(1,0,present(nodal_phi))
    call MPI_Allreduce(flag,low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(flag,high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.low/=high.or.(flag/=0.and.flag/=3)) return
    if(present(eb_geometry)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,state,reactive_nvar(size(species)),comm,valid)
      if(.not.valid) return
      call validate_eb_sampling_geometry(comm,tree,eb_geometry,nodal_phi,valid,periodic)
      if(.not.valid) return
    end if
    call recover_gas(state,recovered,valid)
    if (.not.valid) return
    candidate=state
    if(present(gas_change)) allocate(change)
    call deposit_mpi_amr_tree_spray_field(comm,tree,owners,candidate,reactive_nvar(size(species)), &
      positions,contributions,lower,spacing,periodic,valid,maximum_requests,maximum_received,eb_geometry,nodal_phi,change)
    if (.not.valid) return
    has_particles=size(positions,2)>0
    call MPI_Allreduce(has_particles,any_particles,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (.not.any_particles) then
      if(present(gas_change)) call move_alloc(change%fields,gas_change%fields)
      ok=.true.;return
    end if
    recovered=temperature
    call recover_gas(candidate,recovered,valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(recovered%fields,temperature%fields)
    if(present(gas_change)) call move_alloc(change%fields,gas_change%fields)
    ok=.true.
  contains
    subroutine recover_gas(gas,temp,accepted)
      type(amr_tree_state_3d), intent(in) :: gas
      type(amr_tree_state_3d), intent(inout) :: temp
      logical, intent(out) :: accepted
      if(present(eb_geometry)) then
        call recover_eb_particle_gas(comm,tree,owners,eb_geometry,species,gas,temp,accepted)
      else
        call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,gas,temp,comm,accepted)
      end if
    end subroutine
  end subroutine

  ! Each rank supplies its own particle positions, possibly empty. Only requested
  ! exposed-cell values are exchanged; the owner-local field is never gathered.
  subroutine sample_mpi_amr_tree_spray_field(comm,tree,owners,state,nvar,positions,lower,spacing,periodic, &
      samples,ok,maximum_requests,eb_geometry,nodal_phi,maximum_level)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(amr_tree_state_3d), intent(in) :: state
    real(dp), intent(in) :: positions(:,:),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    real(dp), allocatable, intent(inout) :: samples(:,:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_requests,maximum_level
    type(eb_geometry_3d), optional, intent(in) :: eb_geometry(:)
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(particle_stencil), allocatable :: stencils(:)
    integer, allocatable :: requests(:,:)
    real(dp), allocatable :: values(:,:),candidate(:,:),term(:)
    integer :: ierr,p,q,c,total,ios
    logical :: valid,agreed
    ok=.false.
    call prepare_requests(comm,tree,owners,state,nvar,positions,lower,spacing,periodic, &
      stencils,requests,valid,maximum_requests,eb_geometry,nodal_phi,maximum_level)
    if (.not.valid) return
    allocate(candidate(nvar,size(positions,2)),term(nvar),stat=ios)
    valid=ios==0;call agree()
    if (.not.agreed) return
    call fetch_mpi_amr_tree_cells_3d(tree,owners,state,nvar,requests,values,comm,valid)
    if (.not.valid) return
    total=0;candidate=0._dp
    do p=1,size(stencils)
      do q=1,size(stencils(p)%weights)
        total=total+1;term=stencils(p)%weights(q)*values(:,total)
        do c=1,nvar
          if (term(c)>0._dp) then
            if (candidate(c,p)>huge(1._dp)-term(c)) valid=.false.
          else if (term(c)<0._dp) then
            if (candidate(c,p)< -huge(1._dp)-term(c)) valid=.false.
          end if
        end do
        if (.not.valid) exit
        candidate(:,p)=candidate(:,p)+term
      end do
      if (.not.valid) exit
    end do
    call agree()
    if (.not.agreed) return
    call move_alloc(candidate,samples);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
  ! Contributions are integrated quantities. Destination bins preserve source
  ! particle/stencil order; receives are applied in increasing source-rank order.
  ! An explicit maximum_level deposits raw level-view increments without parent
  ! restriction; use this for temporary source fields, not physical gas updates.
  subroutine deposit_mpi_amr_tree_spray_field(comm,tree,owners,state,nvar,positions,contributions, &
      lower,spacing,periodic,ok,maximum_requests,maximum_received,eb_geometry,nodal_phi,gas_change,maximum_level)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(amr_tree_state_3d), intent(inout) :: state
    real(dp), intent(in) :: positions(:,:),contributions(:,:),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_requests,maximum_received,maximum_level
    type(eb_geometry_3d), optional, intent(in) :: eb_geometry(:)
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(particle_stencil), allocatable :: stencils(:)
    type(amr_tree_state_3d) :: candidate,change
    type(amr_tree_state_3d), optional, intent(inout) :: gas_change
    integer :: export_flag,export_low,export_high
    integer, allocatable :: requests(:,:),send_cells(:,:),recv_cells(:,:)
    integer, allocatable :: counts(:),received(:),offsets(:),displacements(:),cursor(:)
    real(dp), allocatable :: send_values(:,:),recv_values(:,:),term(:)
    integer :: rank,nranks,ierr,ios,capacity,total,incoming,p,q,c,r,slot,node,idx(3)
    logical :: valid,agreed,has_particles,any_particles
    ok=.false.;capacity=1000000
    if (present(maximum_received)) capacity=maximum_received
    export_flag=merge(1,0,present(gas_change))+2*merge(1,0,present(maximum_level))
    call MPI_Allreduce(export_flag,export_low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(export_flag,export_high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.export_low/=export_high) return
    call prepare_requests(comm,tree,owners,state,nvar,positions,lower,spacing,periodic, &
      stencils,requests,valid,maximum_requests,eb_geometry,nodal_phi,maximum_level)
    if (.not.valid) return
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(shape(contributions)==[nvar,size(positions,2)]).and.capacity>=0.and.capacity<=1000000
    valid=valid.and.all(ieee_is_finite(contributions))
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(state%fields(node)%values))) valid=.false.
    end do
    call agree()
    if (.not.agreed) return
    if(present(gas_change)) then
      change=state
      do node=1,size(owners)
        if(owners(node)==rank) change%fields(node)%values=0
      end do
    end if
    has_particles=size(positions,2)>0
    call MPI_Allreduce(has_particles,any_particles,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (.not.any_particles) then
      if(present(gas_change)) call move_alloc(change%fields,gas_change%fields)
      ok=.true.;return
    end if
    total=size(requests,2)
    allocate(counts(nranks),received(nranks),offsets(nranks),displacements(nranks),cursor(nranks), &
      send_cells(4,total),send_values(nvar,total),term(nvar),stat=ios)
    valid=ios==0;call agree()
    if (.not.agreed) return
    counts=0
    do q=1,total
      r=owners(requests(1,q))+1;counts(r)=counts(r)+1
    end do
    call MPI_Alltoall(counts,1,MPI_INTEGER,received,1,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=sum(int(received,int64))<=int(min(capacity,huge(1)/max(4,nvar)),int64)
    call agree()
    if (.not.agreed) return
    incoming=sum(received);offsets(1)=0;displacements(1)=0
    do r=2,nranks
      offsets(r)=offsets(r-1)+counts(r-1)
      displacements(r)=displacements(r-1)+received(r-1)
    end do
    cursor=offsets
    do p=1,size(stencils)
      do q=1,size(stencils(p)%weights)
        node=stencils(p)%cells(1,q);r=owners(node)+1
        slot=cursor(r)+1;cursor(r)=slot
        send_cells(:,slot)=stencils(p)%cells(:,q)
        term=stencils(p)%weights(q)*contributions(:,p)
        if (any(abs(term)>huge(1._dp)*min(1._dp,stencils(p)%volumes(q)))) then
          valid=.false.;exit
        end if
        send_values(:,slot)=term/stencils(p)%volumes(q)
      end do
      if (.not.valid) exit
    end do
    call agree()
    if (.not.agreed) return
    allocate(recv_cells(4,incoming),recv_values(nvar,incoming),stat=ios)
    valid=ios==0;call agree()
    if (.not.agreed) return
    call MPI_Alltoallv(send_cells,4*counts,4*offsets,MPI_INTEGER, &
      recv_cells,4*received,4*displacements,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Alltoallv(send_values,nvar*counts,nvar*offsets,MPI_DOUBLE_PRECISION, &
      recv_values,nvar*received,nvar*displacements,MPI_DOUBLE_PRECISION,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    candidate=state
    do q=1,incoming
      node=recv_cells(1,q);idx=recv_cells(2:4,q);term=recv_values(:,q)
      do c=1,nvar
        associate(value=>candidate%fields(node)%values(c,idx(1),idx(2),idx(3)))
        if (term(c)>0._dp) then
          if (value>huge(1._dp)-term(c)) valid=.false.
        else if (term(c)<0._dp) then
          if (value< -huge(1._dp)-term(c)) valid=.false.
        end if
        end associate
        if(present(gas_change)) then
          associate(value=>change%fields(node)%values(c,idx(1),idx(2),idx(3)))
          if(term(c)>0) then
            if(value>huge(1._dp)-term(c)) valid=.false.
          else if(term(c)<0) then
            if(value< -huge(1._dp)-term(c)) valid=.false.
          end if
          end associate
        end if
      end do
      if (.not.valid) exit
      if(present(gas_change)) change%fields(node)%values(:,idx(1),idx(2),idx(3))= &
        change%fields(node)%values(:,idx(1),idx(2),idx(3))+term
      candidate%fields(node)%values(:,idx(1),idx(2),idx(3))= &
        candidate%fields(node)%values(:,idx(1),idx(2),idx(3))+term
    end do
    call agree()
    if (.not.agreed) return
    if(.not.present(maximum_level)) then
    if(present(eb_geometry)) then
      call average_eb_particle_field(comm,tree,owners,eb_geometry,candidate,nvar,valid)
    else
      call average_down_mpi_amr_tree_3d(tree,owners,candidate,nvar,comm,valid)
    end if
    if (.not.valid) return
    end if
    if(present(gas_change)) then
      if(.not.present(maximum_level)) then
      if(present(eb_geometry)) then
        call average_eb_particle_field(comm,tree,owners,eb_geometry,change,nvar,valid)
      else
        call average_down_mpi_amr_tree_3d(tree,owners,change,nvar,comm,valid)
      end if
      if(.not.valid) return
      end if
      call move_alloc(change%fields,gas_change%fields)
    end if
    call move_alloc(candidate%fields,state%fields);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine prepare_requests(comm,tree,owners,state,nvar,positions,lower,spacing,periodic, &
      stencils,requests,ok,maximum_requests,eb_geometry,nodal_phi,maximum_level)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(amr_tree_state_3d), intent(in) :: state
    real(dp), intent(in) :: positions(:,:),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    type(particle_stencil), allocatable, intent(out) :: stencils(:)
    integer, allocatable, intent(out) :: requests(:,:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_requests,maximum_level
    type(eb_geometry_3d), optional, intent(in) :: eb_geometry(:)
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    integer, allocatable :: lookup(:,:)
    real(dp) :: bounds(6),reference_bounds(6)
    integer :: flags(6),reference_flags(6),ierr,limit,p,total,n,ios,top_level
    logical :: valid,agreed
    ok=.false.;limit=1000000
    if (present(maximum_requests)) limit=maximum_requests
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nvar,comm,valid)
    if (.not.valid) return
    top_level=tree%level_count()-1
    if(present(maximum_level)) top_level=maximum_level
    valid=top_level>=0.and.top_level<tree%level_count()
    valid=valid.and.size(positions,1)==3.and.limit>=0.and.limit<=1000000
    valid=valid.and.all(ieee_is_finite([lower,spacing]))
    valid=valid.and.(present(eb_geometry).eqv.present(nodal_phi))
    if (valid) valid=size(positions,2)<=limit
    call agree()
    if (.not.agreed) return
    bounds=[lower,spacing];reference_bounds=bounds
    call MPI_Bcast(reference_bounds,6,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    flags=[merge(1,0,periodic),merge(1,0,present(eb_geometry)),merge(1,0,present(nodal_phi)),top_level]
    reference_flags=flags
    call MPI_Bcast(reference_flags,6,MPI_INTEGER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(bounds==reference_bounds).and.all(flags==reference_flags)
    call agree()
    if (.not.agreed) return
    if(present(eb_geometry)) then
      call validate_eb_sampling_geometry(comm,tree,eb_geometry,nodal_phi,valid,periodic)
      if(.not.valid) return
    end if
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,lookup,valid,maximum_level)
    call agree()
    if (.not.agreed) return
    allocate(stencils(size(positions,2)),stat=ios)
    valid=ios==0;call agree()
    if (.not.agreed) return
    total=0
    do p=1,size(stencils)
      if (total==limit) then
        valid=.false.;exit
      end if
      if(present(eb_geometry)) then
        call eb_tree_spray_stencil(tree,eb_geometry,nodal_phi,positions(:,p),lower,spacing, &
          stencils(p)%cells,stencils(p)%weights,stencils(p)%volumes,valid,periodic=periodic,maximum_level=maximum_level)
        if(valid) valid=size(stencils(p)%weights)<=limit-total
      else
        call amr_tree_spray_stencil(tree,positions(:,p),lower,spacing,periodic, &
          stencils(p)%cells,stencils(p)%weights,stencils(p)%volumes,valid,maximum_cells=limit-total,maximum_level=maximum_level)
      end if
      if (.not.valid) exit
      total=total+size(stencils(p)%weights)
    end do
    valid=valid.and.total<=huge(1)/max(4,nvar)
    call agree()
    if (.not.agreed) return
    allocate(requests(4,total),stat=ios)
    valid=ios==0;call agree()
    if (.not.agreed) return
    total=0
    do p=1,size(stencils)
      n=size(stencils(p)%weights)
      requests(:,total+1:total+n)=stencils(p)%cells;total=total+n
    end do
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
  ! Compare the replicated geometry actually used by the EB sampler, including
  ! ranks with no queries. Gas field values remain owner-local.
  subroutine validate_eb_sampling_geometry(comm,tree,geometry,phi,ok,periodic)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(in) :: phi(0:,0:,0:)
    logical, intent(out) :: ok
    logical, intent(in) :: periodic(3)
    character(len=16) :: faces(6)
    integer :: axis
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: buffer(:),local(:)
    real(dp) :: bounds(9),reference(9)
    integer :: node,n,ios,ierr,dims(3)
    logical :: valid,agreed
    ok=.false.
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if(valid) valid=size(geometry)==size(tree%nodes)
    call agree()
    if(.not.agreed) return
    valid=compatible_eb_nodal_lattice(shape(phi)-1,domain(:,size(domain,2)))
    valid=valid.and.all(shape(phi)<=1000001).and.all(shape(phi)>1)
    call agree()
    if(.not.agreed) return
    valid=product(int(shape(phi)-1,int64))<=1000000_int64.and.all(ieee_is_finite(phi))
    if(valid) valid=maxval(abs(phi))<=1e100_dp
    do node=1,size(geometry)
      if(.not.valid) exit
      valid=geometry(node)%is_valid()
      if(valid) valid=all([geometry(node)%nx,geometry(node)%ny,geometry(node)%nz]==tree%nodes(node)%shape)
    end do
    call agree()
    if(.not.agreed) return
    faces='outflow'
    do axis=1,3
      if(periodic(axis)) faces(2*axis-1:2*axis)='periodic'
    end do
    valid=periodic_eb_faces_match(geometry(1),faces)
    if(periodic(1)) valid=valid.and.all(phi(0,:,:)==phi(ubound(phi,1),:,:))
    if(periodic(2)) valid=valid.and.all(phi(:,0,:)==phi(:,ubound(phi,2),:))
    if(periodic(3)) valid=valid.and.all(phi(:,:,0)==phi(:,:,ubound(phi,3)))
    call agree()
    if(.not.agreed) return
    dims=shape(phi)
    call MPI_Bcast(dims,3,MPI_INTEGER,0,comm,ierr)
    valid=ierr==MPI_SUCCESS.and.all(dims==shape(phi));call agree()
    if(.not.agreed) return
    allocate(buffer(size(phi)),stat=ios)
    valid=ios==0;call agree()
    if(.not.agreed) return
    buffer=reshape(phi,[size(phi)])
    call MPI_Bcast(buffer,size(buffer),MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=all(buffer==reshape(phi,[size(phi)]))
    call agree()
    if(.not.agreed) return
    deallocate(buffer)
    do node=1,size(geometry)
      associate(g=>geometry(node))
        bounds=[g%x_lower,g%x_upper,g%y_lower,g%y_upper,g%z_lower,g%z_upper,g%dx,g%dy,g%dz]
        reference=bounds
        call MPI_Bcast(reference,9,MPI_DOUBLE_PRECISION,0,comm,ierr)
        if(ierr/=MPI_SUCCESS) return
        valid=all(reference==bounds);call agree()
        if(.not.agreed) return
        n=size(g%volume_fraction)
        allocate(buffer(4*n),local(4*n),stat=ios)
        valid=ios==0;call agree()
        if(.not.agreed) return
        local(:n)=reshape(g%volume_fraction,[n])
        local(n+1:2*n)=reshape(g%cell_centroid_x,[n])
        local(2*n+1:3*n)=reshape(g%cell_centroid_y,[n])
        local(3*n+1:)=reshape(g%cell_centroid_z,[n])
        buffer=local
        call MPI_Bcast(buffer,size(buffer),MPI_DOUBLE_PRECISION,0,comm,ierr)
        if(ierr/=MPI_SUCCESS) return
        valid=all(buffer==local);call agree()
        if(.not.agreed) return
        deallocate(buffer,local)
      end associate
    end do
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine average_eb_particle_field(comm,tree,owners,geometry,state,nvar,ok)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: weighted,candidate
    integer :: rank,ierr,node,n(3),c,i,j,k,ios,child
    real(dp) :: fraction,expected
    logical :: valid,agreed,covered
    ok=.false.
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS) return
    allocate(weighted%fields(size(tree%nodes)),stat=ios)
    valid=ios==0
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(weighted%fields(node)%values(nvar+1,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      do c=1,nvar
        weighted%fields(node)%values(c,:,:,:)=state%fields(node)%values(c,:,:,:)*geometry(node)%volume_fraction
      end do
      weighted%fields(node)%values(nvar+1,:,:,:)=geometry(node)%volume_fraction
    end do
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    call average_down_mpi_amr_tree_3d(tree,owners,weighted,nvar+1,comm,valid)
    if(.not.valid) return
    candidate=state
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        covered=.false.
        do child=node+1,size(tree%nodes)
          if(tree%nodes(child)%parent/=node) cycle
          associate(box=>tree%nodes(child)%patch)
            if(i>=box%coarse_i_lower.and.i<=box%coarse_i_upper.and. &
               j>=box%coarse_j_lower.and.j<=box%coarse_j_upper.and. &
               k>=box%coarse_k_lower.and.k<=box%coarse_k_upper) covered=.true.
          end associate
        end do
        if(.not.covered) cycle
        fraction=weighted%fields(node)%values(nvar+1,i,j,k)
        expected=geometry(node)%volume_fraction(i,j,k)
        if(abs(fraction-expected)>128*epsilon(1._dp)*max(fraction,expected)) valid=.false.
        if(expected<=0.or.fraction<=0) cycle
        if(any(abs(weighted%fields(node)%values(:nvar,i,j,k))>huge(1._dp)*min(1._dp,fraction))) then
          valid=.false.;cycle
        end if
        candidate%fields(node)%values(:,i,j,k)=weighted%fields(node)%values(:nvar,i,j,k)/fraction
      end do
      end do
      end do
    end do
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    call move_alloc(candidate%fields,state%fields);ok=.true.
  end subroutine

  subroutine recover_eb_particle_gas(comm,tree,owners,geometry,species,state,temperature,ok)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state
    type(amr_tree_state_3d), intent(inout) :: temperature
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    real(dp) :: primitive(reactive_nvar(size(species))),t,sound
    integer :: rank,ierr,node,i,j,k,n(3)
    logical :: valid,agreed,cell_ok
    ok=.false.
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,valid)
    if(.not.valid) return
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS) return
    candidate=temperature;valid=.true.
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
        call reactive_conserved_to_primitive(species,state%fields(node)%values(:,i,j,k), &
          temperature%fields(node)%values(1,i,j,k),primitive,t,sound,cell_ok)
        if(.not.cell_ok) then
          valid=.false.
        else
          candidate%fields(node)%values(1,i,j,k)=t
        end if
      end do
      end do
      end do
    end do
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    call move_alloc(candidate%fields,temperature%fields);ok=.true.
  end subroutine

end module
