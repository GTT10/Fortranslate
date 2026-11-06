module mpi_eb_tree_redistribution_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_reactive_redistribution_3d_mod, only: state_redistribution_graph_slopes_3d,limit_state_recipient_slopes, &
    state_redistribution_boundary_nodes_3d,state_redist_register_3d,compute_state_redist_level_stage_3d, &
    initialize_state_redist_tree_register_3d,restrict_state_redist_patch_source_3d,srd_coarse_boundary, &
    assemble_state_redist_patch_rereflux_3d,state_redist_patch_spill_requests_3d, &
    state_redist_patch_rereflux_budget_3d, &
    state_redist_rereflux_budget_valid_3d,state_redist_box_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,periodic_eb_faces_match,extract_eb_geometry_3d
  use eb_reactive_two_level_3d_mod, only: matching_eb_metrics_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: fetch_mpi_amr_tree_cells_3d,validate_mpi_amr_tree_local_3d, &
    validate_mpi_amr_tree_metadata_3d,add_mpi_amr_tree_cells_3d
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d
  implicit none
  private
  public :: compute_state_redist_mpi_level_stage_3d,assemble_state_redist_mpi_interface_source_3d
  public :: assemble_state_redist_mpi_interface_rereflux_3d
  public :: redistribute_sdc_mpi_eb_tree_flux_rate_3d,redistribute_sdc_mpi_eb_tree_level_rate_3d
  public :: redistribute_sdc_mpi_eb_tree_state_level_rate_3d
  public :: redistribute_mpi_eb_tree_state_3d,redistribute_sdc_mpi_eb_tree_state_rate_3d
contains
  subroutine assemble_state_redist_mpi_interface_rereflux_3d(tree,owners,geometry,registers,nv,depth, &
      increment,comm,ok,boundary_faces,spill,spill_active)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nv,depth
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(state_redist_register_3d), intent(in) :: registers(:)
    type(amr_tree_state_3d), intent(inout) :: increment
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d), optional, intent(inout) :: spill
    logical, optional, intent(out) :: spill_active
    type(amr_tree_state_3d) :: source,pending,pending_spill
    type(state_redist_box_3d) :: box
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: h
    logical :: periodic(3)
    real(dp), allocatable :: patch_increment(:,:,:,:),halo_spill(:,:,:,:),values(:,:),all_values(:,:)
    integer, allocatable :: requests(:,:),all_requests(:,:)
    real(dp) :: local_budget(3,nv),global_budget(3,nv),patch_budget(3,nv),spill_budget(2,nv)
    integer :: rank,status,p,n(3),capacity,count,m
    logical :: good,valid,ready,local_active,global_active
    ok=.false.
    if(present(spill_active)) spill_active=.false.
    call assemble_state_redist_mpi_interface_source_3d(tree,owners,geometry,registers,nv,depth, &
      source,comm,ready,boundary_faces)
    if(.not.ready) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
    call MPI_Allreduce(good,ready,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.ready) return
    periodic=.false.
    if(present(boundary_faces)) periodic=[boundary_faces(1)=='periodic',boundary_faces(3)=='periodic', &
      boundary_faces(5)=='periodic']
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    allocate(pending%fields(size(tree%nodes)));valid=.true.;local_budget=0._dp
    capacity=0;count=0
    if(present(spill).and.depth>0) then
      allocate(pending_spill%fields(size(tree%nodes)))
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        n=tree%nodes(p)%shape
        allocate(pending_spill%fields(p)%values(nv,n(1),n(2),n(3)),source=0._dp)
        if(tree%nodes(p)%level==depth) capacity=capacity+product(n+2)-product(n)
      enddo
      allocate(all_requests(4,capacity),all_values(nv,capacity))
    endif
    do p=1,size(tree%nodes)
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape
      allocate(pending%fields(p)%values(nv,n(1),n(2),n(3)),source=0._dp)
      if(tree%nodes(p)%level/=depth) cycle
      h=(geometry(p)%nx-n(1))/2
      box%valid_lower=h+1;box%valid_upper=h+n
      box%domain_lower=1_int64+int(h,int64)-origin(:,p)
      box%domain_upper=domain(:,depth+1)+int(h,int64)-origin(:,p)
      box%periodic=periodic
      if(present(spill).and.depth>0) then
        call assemble_state_redist_patch_rereflux_3d(geometry(p),registers(p),n,source%fields(p)%values, &
          patch_increment,good,halo_spill,box)
      else
        call assemble_state_redist_patch_rereflux_3d(geometry(p),registers(p),n,source%fields(p)%values, &
          patch_increment,good,box=box)
      endif
      if(.not.good) then
        valid=.false.;exit
      endif
      call state_redist_patch_rereflux_budget_3d(geometry(p),n,source%fields(p)%values, &
        patch_increment,patch_budget,good)
      if(.not.good) then
        valid=.false.;exit
      endif
      local_budget=local_budget+patch_budget
      if(present(spill).and.depth>0) then
        call state_redist_patch_spill_requests_3d(tree,p,geometry,halo_spill,requests,values,spill_budget,good,boundary_faces)
        if(.not.good) then
          valid=.false.;exit
        endif
        m=size(requests,2)
        all_requests(:,count+1:count+m)=requests
        all_values(:,count+1:count+m)=values
        count=count+m
        local_budget(2,:)=local_budget(2,:)+spill_budget(1,:)
        local_budget(3,:)=local_budget(3,:)+spill_budget(2,:)
      endif
      call move_alloc(patch_increment,pending%fields(p)%values)
    enddo
    call MPI_Allreduce(valid,ready,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.ready) return
    if(present(spill).and.depth>0) then
      call add_mpi_amr_tree_cells_3d(tree,owners,pending_spill,nv,all_requests(:,:count), &
        all_values(:,:count),comm,ready,maximum_depth=depth-1)
      if(.not.ready) return
    endif
    call MPI_Allreduce(local_budget,global_budget,size(local_budget),MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if(status/=MPI_SUCCESS) return
    if(.not.state_redist_rereflux_budget_valid_3d(global_budget)) return
    if(present(spill)) then
      if(depth>0) then
        local_active=count>0
        call MPI_Allreduce(local_active,global_active,1,MPI_LOGICAL,MPI_LOR,comm,status)
        if(status/=MPI_SUCCESS) return
        call move_alloc(pending_spill%fields,spill%fields)
        if(present(spill_active)) spill_active=global_active
      endif
    endif
    call move_alloc(pending%fields,increment%fields);ok=.true.
  end subroutine assemble_state_redist_mpi_interface_rereflux_3d

  ! Owner-local interface ledger restriction, sparse deposit, then exact
  ! coarse peer/periodic halo refresh. No rank gathers a complete register.
  subroutine assemble_state_redist_mpi_interface_source_3d(tree,owners,geometry,registers,nv,depth, &
      source_halos,comm,ok,boundary_faces)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nv,depth
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(state_redist_register_3d), intent(in) :: registers(:)
    type(amr_tree_state_3d), intent(inout) :: source_halos
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type records
      integer, allocatable :: target(:,:)
      real(dp), allocatable :: value(:,:)
    end type
    type(records), allocatable :: patches(:)
    type(amr_tree_state_3d) :: sources,pending
    type(state_redist_register_3d) :: expected
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),requests(:,:),images(:,:)
    real(dp), allocatable :: contributions(:,:),values(:,:)
    integer :: rank,status,p,q,n(3),dims(3),extra(3),width,i,j,k,c,idx(3),count,offset,query
    integer :: controls(7),low(7),high(7)
    integer(int64) :: total,cell(3)
    character(len=32) :: faces(6)
    logical :: good,ready,valid,inside,periodic(3)
    ok=.false.;faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    call validate_mpi_amr_tree_metadata_3d(tree,owners,nv,comm,ready)
    if(.not.ready) return
    valid=valid_reactive_boundaries_3d(faces).and.depth>=0.and.depth+1<tree%level_count()
    valid=valid.and.size(geometry)==size(tree%nodes).and.size(registers)==size(tree%nodes)
    call agree(valid,ready)
    if(.not.ready) return
    controls(1)=depth
    controls(2:7)=merge(1,0,faces=='periodic')
    call MPI_Allreduce(controls,low,7,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,7,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    extra=[geometry(1)%nx,geometry(1)%ny,geometry(1)%nz]-tree%nodes(1)%shape
    valid=all(extra>=2).and.all(modulo(extra,2)==0).and.all(extra==extra(1))
    call agree(valid,ready)
    if(.not.ready) return
    width=extra(1)/2
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready,padding=width)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
    call agree(good,ready)
    if(.not.ready) return
    count=size(tree%nodes);periodic=faces(1:6:2)=='periodic'
    allocate(sources%fields(count),pending%fields(count),patches(count));valid=.true.;total=0
    do p=1,count
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape
      ! The sparse cell API checks all owned valid fields; inactive levels
      ! carry zeros and are not transmitted unless explicitly requested.
      allocate(sources%fields(p)%values(nv,n(1),n(2),n(3)),source=0._dp)
      if(tree%nodes(p)%level/=depth.and.tree%nodes(p)%level/=depth+1) cycle
      dims=n+2*width
      if(.not.registers(p)%is_valid(nv,dims)) then
        valid=.false.;exit
      endif
      call initialize_state_redist_tree_register_3d(tree,p,nv,width,expected,good,faces,registers(p)%scale)
      if(.not.good) then
        valid=.false.;exit
      endif
      if(allocated(expected%coarse_flags).neqv.allocated(registers(p)%coarse_flags)) then
        valid=.false.;exit
      endif
      if(allocated(expected%fine_ghost).neqv.allocated(registers(p)%fine_ghost)) then
        valid=.false.;exit
      endif
      if(.not.allocated(registers(p)%valid_cells)) then
        valid=.false.;exit
      endif
      if(any(expected%valid_cells.neqv.registers(p)%valid_cells)) then
        valid=.false.;exit
      endif
      if(allocated(expected%coarse_flags)) then
        if(any(expected%coarse_flags/=registers(p)%coarse_flags)) then
          valid=.false.;exit
        endif
      endif
      if(allocated(expected%fine_ghost)) then
        if(any(expected%fine_ghost.neqv.registers(p)%fine_ghost)) then
          valid=.false.;exit
        endif
      endif
      if(tree%nodes(p)%level==depth) then
        allocate(pending%fields(p)%values(nv,dims(1),dims(2),dims(3)),source=0._dp)
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          idx=[i,j,k]+width
          if(registers(p)%coarse_flags(idx(1),idx(2),idx(3))/=srd_coarse_boundary) cycle
          sources%fields(p)%values(:,i,j,k)=registers(p)%coarse_delta(:,idx(1),idx(2),idx(3))
        enddo
        enddo
        enddo
      else
        call restrict_state_redist_patch_source_3d(tree,p,geometry(p),registers(p),nv, &
          patches(p)%target,patches(p)%value,good,faces)
        if(.not.good) then
          valid=.false.;exit
        endif
        total=total+int(size(patches(p)%target,2),int64)
      endif
    enddo
    valid=valid.and.total<=int(huge(1)/max(4,nv),int64)
    call agree(valid,ready)
    if(.not.ready) return
    allocate(requests(4,int(total)),contributions(nv,int(total)));offset=0
    do p=1,count
      if(.not.allocated(patches(p)%target)) cycle
      query=size(patches(p)%target,2)
      requests(:,offset+1:offset+query)=patches(p)%target
      contributions(:,offset+1:offset+query)=patches(p)%value
      offset=offset+query
    enddo
    ! Check remote targets from topology on the sender. No remote mask or
    ! coarse payload must be replicated to establish boundary ownership.
    do q=1,count
      if(.not.any(requests(1,:)==q)) cycle
      call initialize_state_redist_tree_register_3d(tree,q,nv,1,expected,good,faces)
      if(.not.good) then
        valid=.false.;exit
      endif
      do query=1,size(requests,2)
        if(requests(1,query)/=q) cycle
        idx=requests(2:4,query)+1
        if(expected%coarse_flags(idx(1),idx(2),idx(3))/=srd_coarse_boundary) valid=.false.
      enddo
    enddo
    call agree(valid,ready)
    if(.not.ready) return
    call add_mpi_amr_tree_cells_3d(tree,owners,sources,nv,requests,contributions,comm,ready)
    if(.not.ready) return
    total=0
    do p=1,count
      if(owners(p)/=rank.or.tree%nodes(p)%level/=depth) cycle
      n=tree%nodes(p)%shape
      do k=width,width+n(3)+1
      do j=width,width+n(2)+1
      do i=width,width+n(1)+1
        call locate(p,[i,j,k],q,idx,inside)
        if(inside) total=total+1
      enddo
      enddo
      enddo
    enddo
    call agree(total<=int(huge(1)/max(4,nv),int64),ready)
    if(.not.ready) return
    allocate(images(4,int(total)));query=0
    do p=1,count
      if(owners(p)/=rank.or.tree%nodes(p)%level/=depth) cycle
      n=tree%nodes(p)%shape
      do k=width,width+n(3)+1
      do j=width,width+n(2)+1
      do i=width,width+n(1)+1
        call locate(p,[i,j,k],q,idx,inside)
        if(.not.inside) cycle
        query=query+1;images(:,query)=[q,idx]
      enddo
      enddo
      enddo
    enddo
    call fetch_mpi_amr_tree_cells_3d(tree,owners,sources,nv,images,values,comm,ready)
    if(.not.ready) return
    query=0
    do p=1,count
      if(owners(p)/=rank.or.tree%nodes(p)%level/=depth) cycle
      n=tree%nodes(p)%shape
      do k=width,width+n(3)+1
      do j=width,width+n(2)+1
      do i=width,width+n(1)+1
        call locate(p,[i,j,k],q,idx,inside)
        if(.not.inside) cycle
        query=query+1;pending%fields(p)%values(:,i,j,k)=values(:,query)
      enddo
      enddo
      enddo
    enddo
    call move_alloc(pending%fields,source_halos%fields);ok=.true.
  contains
    subroutine locate(p,image,q,index,inside)
      integer, intent(in) :: p,image(3)
      integer, intent(out) :: q,index(3)
      logical, intent(out) :: inside
      integer :: c,node
      cell=origin(:,p)+int(image-width-1,int64);inside=.false.;q=0;index=0
      do c=1,3
        if(periodic(c)) then
          cell(c)=modulo(cell(c),domain(c,depth+1))
        else if(cell(c)<0.or.cell(c)>=domain(c,depth+1)) then
          return
        endif
      enddo
      do node=1,count
        if(tree%nodes(node)%level/=depth) cycle
        if(any(cell<origin(:,node)).or.any(cell>=origin(:,node)+int(tree%nodes(node)%shape,int64))) cycle
        q=node;index=int(cell-origin(:,node))+1;inside=.true.;return
      enddo
    end subroutine
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine assemble_state_redist_mpi_interface_source_3d

  ! Prepared state/RHS halos are owner-local. The shared patch backend computes
  ! private stage outputs; every rank agrees before either output is published.
  subroutine compute_state_redist_mpi_level_stage_3d(tree,owners,geometry,old_halos,rhs_halos,nv,depth,dt, &
      rate,registers,comm,ok,boundary_faces,target_volume_fraction,max_order,register_scale,update_scales)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nv,depth
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: old_halos,rhs_halos
    type(amr_tree_state_3d), optional, intent(in) :: update_scales
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate
    type(state_redist_register_3d), allocatable, intent(inout) :: registers(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: target_volume_fraction,register_scale
    integer, optional, intent(in) :: max_order
    type(amr_tree_state_3d) :: pending_rate
    type(state_redist_register_3d), allocatable :: pending_registers(:)
    character(len=32) :: faces(6)
    real(dp) :: target,scale,values(3),low_values(3),high_values(3)
    integer :: order,rank,status,p,width,extra(3),controls(13),low(13),high(13)
    logical :: valid,ready
    ok=.false.;faces='outflow';target=.5_dp;scale=1._dp;order=0
    if(present(boundary_faces)) faces=boundary_faces
    if(present(target_volume_fraction)) target=target_volume_fraction
    if(present(register_scale)) scale=register_scale
    if(present(max_order)) order=max_order
    valid=ieee_is_finite(dt).and.ieee_is_finite(target).and.ieee_is_finite(scale)
    valid=valid.and.valid_reactive_boundaries_3d(faces)
    call agree(valid,ready)
    if(.not.ready) return
    controls(:6)=merge(1,0,faces=='periodic')+merge(2,0,faces=='inflow')+merge(3,0,faces=='slip_wall')+ &
      merge(4,0,faces=='no_slip_wall')+merge(5,0,faces=='isothermal_wall')
    controls(7:12)=[depth,order,merge(1,0,present(boundary_faces)),merge(1,0,present(max_order)), &
      merge(1,0,present(target_volume_fraction)),merge(1,0,present(register_scale))]
    controls(13)=merge(1,0,present(update_scales))
    call MPI_Allreduce(controls,low,13,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,13,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    values=[dt,target,scale]
    call MPI_Allreduce(values,low_values,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(values,high_values,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low_values/=high_values)) return
    call validate_mpi_amr_tree_local_3d(tree,owners,rate,nv,comm,ready)
    if(.not.ready) return
    valid=size(geometry)==size(tree%nodes)
    call agree(valid,ready)
    if(.not.ready) return
    extra=[geometry(1)%nx,geometry(1)%ny,geometry(1)%nz]-tree%nodes(1)%shape
    valid=all(extra>=10).and.all(modulo(extra,2)==0).and.all(extra==extra(1))
    call agree(valid,ready)
    if(.not.ready) return
    width=extra(1)/2
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready,padding=width)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    valid=allocated(old_halos%fields).and.allocated(rhs_halos%fields)
    if(valid) valid=size(old_halos%fields)==size(tree%nodes).and.size(rhs_halos%fields)==size(tree%nodes)
    if(valid) then
      do p=1,size(tree%nodes)
        if(tree%nodes(p)%level/=depth.or.owners(p)==rank) cycle
        if(allocated(old_halos%fields(p)%values).or.allocated(rhs_halos%fields(p)%values)) valid=.false.
      enddo
    endif
    if(present(update_scales)) then
      valid=valid.and.allocated(update_scales%fields)
      if(valid) valid=size(update_scales%fields)==size(tree%nodes)
      if(valid) then
        do p=1,size(tree%nodes)
          if(tree%nodes(p)%level/=depth.or.owners(p)==rank) cycle
          if(allocated(update_scales%fields(p)%values)) valid=.false.
        enddo
      endif
    endif
    call agree(valid,ready)
    if(.not.ready) return
    pending_rate=rate
    if(allocated(registers)) pending_registers=registers
    call compute_state_redist_level_stage_3d(tree,geometry,old_halos,rhs_halos,nv,depth,dt,pending_rate, &
      pending_registers,valid,faces,target,order,scale,owners==rank,update_scales)
    call agree(valid,ready)
    if(.not.ready) return
    rate=pending_rate
    call move_alloc(pending_registers,registers)
    ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine agree
  end subroutine compute_state_redist_mpi_level_stage_3d
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
      rate,transfer,comm,ok,target_volume_fraction,boundary_faces,max_order, &
      boundary_halos,boundary_geometry)
    type(amr_tree_state_3d), optional, intent(in) :: boundary_halos
    type(eb_geometry_3d), optional, intent(in) :: boundary_geometry(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar,depth
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,rhs
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate,transfer
    type(MPI_Comm), intent(in) :: comm
    integer, optional, intent(in) :: max_order
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: target_volume_fraction
    type(amr_tree_state_3d) :: selected,result,exchange,correction,provisional
    integer :: p,rank,status,lower,upper,order,order_low,order_high,i,j,k,n(3)
    real(dp) :: dt_low,dt_high
    logical :: valid,ready
    ok=.false.;order=0
    if(present(max_order)) order=max_order
    call agree(order==0.or.order==2.or.order==4,ready)
    if(.not.ready) return
    call MPI_Allreduce(order,order_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(order,order_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.order_low/=order_high) return
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
    if(order>0) then
      call validate_mpi_amr_tree_local_3d(tree,owners,state,nvar,comm,ready)
      if(.not.ready) return
      call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nvar,comm,ready)
      if(.not.ready) return
      provisional=state;valid=.true.
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        n=tree%nodes(p)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          provisional%fields(p)%values(:,i,j,k)=0
          if(geometry(p)%volume_fraction(i,j,k)<=0) then
            selected%fields(p)%values(:,i,j,k)=0;cycle
          endif
          if(any(.not.ieee_is_finite(state%fields(p)%values(:,i,j,k))).or. &
            any(.not.ieee_is_finite(rhs%fields(p)%values(:,i,j,k)))) then
            valid=.false.;cycle
          endif
          provisional%fields(p)%values(:,i,j,k)=state%fields(p)%values(:,i,j,k)+dt*rhs%fields(p)%values(:,i,j,k)
        enddo
        enddo
        enddo
      enddo
      call agree(valid,ready)
      if(.not.ready) return
      call redistribute_mpi_eb_tree_state_3d(tree,owners,geometry,provisional,nvar,correction,comm,ready, &
        target_volume_fraction,boundary_faces,correction_depth=depth,max_order=order, &
        boundary_halos=boundary_halos,boundary_geometry=boundary_geometry)
      if(.not.ready) return
      result=selected;exchange=correction;valid=.true.
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        if(tree%nodes(p)%level==depth) then
          result%fields(p)%values=result%fields(p)%values+correction%fields(p)%values/dt
          exchange%fields(p)%values=0
        endif
        valid=valid.and.all(ieee_is_finite(result%fields(p)%values)).and. &
          all(ieee_is_finite(exchange%fields(p)%values))
      enddo
      call agree(valid,ready)
      if(.not.ready) return
      call move_alloc(result%fields,rate%fields);call move_alloc(exchange%fields,transfer%fields)
      ok=.true.;return
    endif
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
      target_volume_fraction,boundary_faces,correction_depth,max_order, &
      boundary_halos,boundary_geometry)
    type(amr_tree_state_3d), optional, intent(in) :: boundary_halos
    type(eb_geometry_3d), optional, intent(in) :: boundary_geometry(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nv
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(MPI_Comm), intent(in) :: comm
    type(amr_tree_state_3d), intent(in) :: input
    type(amr_tree_state_3d), intent(inout) :: result
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: target_volume_fraction
    character(len=*), optional, intent(in) :: boundary_faces(6)
    integer, optional, intent(in) :: max_order
    integer, optional, intent(in) :: correction_depth
    integer :: depth,depth_low,depth_high
    type :: neighborhood
      integer, allocatable :: member(:)
    end type
    type(neighborhood), allocatable :: groups(:)
    type(amr_tree_state_3d) :: ids,candidate,averages
    type(eb_geometry_3d) :: checked_geometry
    integer(int64), allocatable :: origin(:,:),domain(:,:),scale(:)
    integer, allocatable :: ratios(:),cells(:,:),receiver(:),memberships(:),queue(:),seen(:),neighbors(:),directions(:)
    real(dp), allocatable :: volume(:),target(:),beta(:),self(:),values(:,:),incoming(:,:),mean(:)
    integer(int64) :: maximum,cell_lo(3),cell_hi(3)
    character(len=32) :: faces(6)
    integer :: p,n(3),i,j,k,child,lo(3),hi(3),nc,root,other,m,count,head,tail,current,allocation_status
    real(dp) :: fraction,total,weight,normal(3),score,best
    logical, allocatable :: stencil_connected(:)
    logical :: valid,fallback,ready
    integer :: rank,status,q,request_count,controls(6),low(6),high(6)
    integer, allocatable :: requests(:,:),lookup(:),reverse_lookup(:)
    logical, allocatable :: needed(:)
    real(dp) :: fraction_low,fraction_high
    real(dp), allocatable :: center(:,:),lattice(:,:),centroid(:,:),width(:,:),fraction_cell(:), &
      stencil_centroid(:,:),stencil_center(:,:),stencil_width(:,:),positions(:,:),lower_bound(:),upper_bound(:), &
      graph_values(:,:),root_slopes(:,:),local_bound(:)
    real(dp), allocatable :: stencil_value(:,:),stencil_fraction(:),stencil_fit(:,:),ghost_value(:,:)
    real(dp) :: ghost_point(3,125),ghost_fit(3,125),ghost_center(3,125),ghost_fraction(125),root_centroid(3)
    logical :: ghost_near(125),extdir(6)
    integer :: ghost_count,global_cell(3),global_dims(3),extdir_bounds(2,3),pool_size,halo_mode,halo_low,halo_high
    integer :: selected_order,order_low,order_high,packet_nv
    ok=.false.;fraction=.5_dp;faces='outflow'
    selected_order=0
    if(present(max_order)) selected_order=max_order
    call agree(selected_order==0.or.selected_order==2.or.selected_order==4,ready)
    if(.not.ready) return
    call MPI_Allreduce(selected_order,order_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(selected_order,order_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.order_low/=order_high) return
    packet_nv=nv
    if(selected_order>0) then
      call agree(nv>0.and.nv<=huge(1)/4,ready)
      if(.not.ready) return
      packet_nv=4*nv
    endif
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
    ! Presence and validation are collective, including idle ranks.
    halo_mode=0
    if(present(boundary_halos).neqv.present(boundary_geometry)) halo_mode=-1
    if(present(boundary_halos).and.present(boundary_geometry)) halo_mode=1
    call MPI_Allreduce(halo_mode,halo_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(halo_mode,halo_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.halo_low<0.or.halo_low/=halo_high) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    valid=.true.
    if(present(boundary_halos)) then
      valid=allocated(boundary_halos%fields)
      if(valid) valid=size(boundary_halos%fields)==size(tree%nodes).and.size(boundary_geometry)==size(tree%nodes)
      if(valid) then
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          if(.not.allocated(boundary_halos%fields(p)%values)) then
            valid=.false.;cycle
          endif
          n=tree%nodes(p)%shape
          if(any(shape(boundary_halos%fields(p)%values)/=[nv,n+4])) then
            valid=.false.;cycle
          endif
          if(any(.not.ieee_is_finite(boundary_halos%fields(p)%values))) valid=.false.
          if(any([boundary_geometry(p)%nx,boundary_geometry(p)%ny,boundary_geometry(p)%nz]/=n+4)) then
            valid=.false.;cycle
          endif
          call extract_eb_geometry_3d(boundary_geometry(p),[3,3,3],n+2,checked_geometry,ready)
          if(.not.ready) then
            valid=.false.;cycle
          endif
          if(.not.matching_eb_metrics_3d(geometry(p),checked_geometry)) valid=.false.
        enddo
      endif
    endif
    call agree(valid,ready)
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
      allocation_status=0
      if(owners(p)==rank.and.selected_order>0) then
        deallocate(averages%fields(p)%values)
        allocate(averages%fields(p)%values(packet_nv,n(1),n(2),n(3)),stat=allocation_status)
      endif
      call agree(allocation_status==0,ready)
      if(.not.ready) return
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
      if(selected_order==0.and.depth>=0.and.tree%nodes(cells(1,root))%level/=depth) cycle
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
      if(selected_order==0.and.depth>=0.and.tree%nodes(cells(1,root))%level/=depth) cycle
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
      averages%fields(p)%values(1:nv,i,j,k)=mean/total
      valid=valid.and.all(ieee_is_finite(averages%fields(p)%values(:,i,j,k)))
    end do
    call agree(valid,ready)
    if(.not.ready) return
    if(selected_order>0) then
      call reconstruct_owned_graph_slopes(ready)
      if(.not.ready) return
    endif
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
    call fetch_mpi_amr_tree_cells_3d(tree,owners,averages,packet_nv,requests,incoming,comm,ready)
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
          call reconstructed_mean(root,other,q,mean)
          if(depth>=0) mean=mean-input%fields(p)%values(:,i,j,k)
          if(selected_order==0) then
            candidate%fields(p)%values(:,i,j,k)=candidate%fields(p)%values(:,i,j,k)+mean/real(memberships(other),dp)
          else
            candidate%fields(p)%values(:,i,j,k)=candidate%fields(p)%values(:,i,j,k)+ &
              (1._dp/real(memberships(other),dp))*mean
          endif
        end do
      else
        p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
        if(owners(p)==rank) then
          call reconstructed_mean(root,root,q,mean)
          if(depth>=0) mean=mean-input%fields(p)%values(:,i,j,k)
          candidate%fields(p)%values(:,i,j,k)=candidate%fields(p)%values(:,i,j,k)+ &
            self(root)*real(merge(memberships(root),1,selected_order==0),dp)*mean
        end if
        other=receiver(root)
        if(other>0) then
          p=cells(1,other);i=cells(2,other);j=cells(3,other);k=cells(4,other)
          if(owners(p)==rank) then
            call reconstructed_mean(root,other,q,mean)
            if(depth>=0) mean=mean-input%fields(p)%values(:,i,j,k)
            candidate%fields(p)%values(:,i,j,k)=candidate%fields(p)%values(:,i,j,k)+ &
              (beta(root)/real(merge(1,memberships(other),selected_order==0),dp))*mean
          end if
        end if
      end if
    end do
    if(.not.fallback.and.selected_order==0) then
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
    function displacement(other,root) result(delta)
      integer, intent(in) :: other,root
      real(dp) :: delta(3),length
      integer :: a
      delta=center(:,other)-center(:,root)
      do a=1,3
        if(faces(2*a-1)/='periodic') cycle
        length=real(domain(a,size(ratios)),dp)
        delta(a)=modulo(delta(a)+.5_dp*length,length)-.5_dp*length
      enddo
    end function displacement
    subroutine reconstructed_mean(root,other,q,mean)
      integer, intent(in) :: root,other,q
      real(dp), intent(out) :: mean(:)
      real(dp) :: offset(3)
      integer :: c
      mean=incoming(1:nv,q)
      if(selected_order==0) return
      offset=displacement(other,root)-(centroid(:,root)-center(:,root))
      do c=1,nv
        mean(c)=mean(c)+dot_product(offset,incoming(nv+3*c-2:nv+3*c,q))
      enddo
    end subroutine reconstructed_mean
    subroutine build_graph_stencil(r,last)
      integer, intent(in) :: r
      integer, intent(out) :: last
      integer :: h,t,next,member_count
      real(dp) :: offset(3),raw_offset(3)
      seen=0
      h=1;last=1;queue(1)=r;seen(r)=1
        do while(h<=last)
          current=queue(h);call adjacent(member_count)
          do t=1,member_count
            next=neighbors(t)
            if(seen(next)==1) cycle
            raw_offset=displacement(next,r)-(center(:,next)-lattice(:,next))+(center(:,r)-lattice(:,r))
            if(any(abs(raw_offset/width(:,r))>2+512*epsilon(1._dp))) cycle
            last=last+1;queue(last)=next;seen(next)=1
          enddo
          h=h+1
        enddo
        stencil_connected(:last)=.true.
        ! Include noncovered leaf values for the grown AMReX stencil even
        ! when an aperture path does not connect them to the root.
        do next=1,nc
          if(seen(next)==1) cycle
          raw_offset=displacement(next,r)-(center(:,next)-lattice(:,next))+(center(:,r)-lattice(:,r))
          if(any(abs(raw_offset/width(:,r))>2+512*epsilon(1._dp))) cycle
          last=last+1;queue(last)=next;stencil_connected(last)=.false.
        enddo
        do t=1,last
          next=queue(t);offset=displacement(next,r)
          stencil_centroid(:,t)=(offset+centroid(:,next)-center(:,next))/width(:,r)
          stencil_center(:,t)=(offset-(center(:,next)-lattice(:,next))+(center(:,r)-lattice(:,r)))/width(:,r)
          stencil_width(:,t)=width(:,next)/width(:,r)
        enddo
    end subroutine build_graph_stencil
    subroutine reconstruct_owned_graph_slopes(success)
      logical, intent(out) :: success
      integer :: r,t,last,next,nmember
      real(dp) :: w,sum_weight,shift(3)
      logical :: reconstructed,condition
      success=.false.
      call agree(nc<=huge(1)-125,condition)
      if(.not.condition) return
      pool_size=nc+125
      allocate(stencil_value(nv,pool_size),stencil_fraction(pool_size),stencil_fit(3,pool_size), &
        ghost_value(nv,125),stat=allocation_status)
      call agree(allocation_status==0,condition)
      if(.not.condition) return
      allocate(center(3,nc),lattice(3,nc),centroid(3,nc),width(3,nc),fraction_cell(nc), &
        stencil_centroid(3,pool_size),stencil_center(3,pool_size),stencil_width(3,pool_size), &
        positions(3,nc),stencil_connected(pool_size), &
        lower_bound(nv),upper_bound(nv),local_bound(nv),root_slopes(3,nv),stat=allocation_status)
      call agree(allocation_status==0,condition)
      if(.not.condition) return
      do r=1,nc
        p=cells(1,r);i=cells(2,r);j=cells(3,r);k=cells(4,r)
        w=real(scale(tree%nodes(p)%level+1),dp)
        width(:,r)=w;lattice(:,r)=(real(origin(:,p),dp)+real([i,j,k],dp)-.5_dp)*w
        associate(g=>geometry(p))
          center(:,r)=lattice(:,r)+w*[g%cell_centroid_x(i,j,k),g%cell_centroid_y(i,j,k),g%cell_centroid_z(i,j,k)]
          fraction_cell(r)=g%volume_fraction(i,j,k)
        end associate
      enddo
      do r=1,nc
        shift=0;sum_weight=0
        if(fallback) then
          nmember=1
          if(allocated(groups(r)%member)) nmember=size(groups(r)%member)
          do t=1,nmember
            next=r
            if(allocated(groups(r)%member)) next=groups(r)%member(t)
            w=volume(next)/real(memberships(next),dp)
            shift=shift+w*displacement(next,r);sum_weight=sum_weight+w
          enddo
        else
          sum_weight=self(r)*volume(r);next=receiver(r)
          if(next>0) then
            w=beta(r)*volume(next)/real(memberships(next),dp)
            shift=w*displacement(next,r);sum_weight=sum_weight+w
          endif
        endif
        if(sum_weight<=0) return
        centroid(:,r)=center(:,r)+shift/sum_weight
      enddo
      local_bound=huge(1._dp)
      do r=1,nc
        p=cells(1,r);i=cells(2,r);j=cells(3,r);k=cells(4,r)
        if(owners(p)/=rank) cycle
        local_bound=min(local_bound,input%fields(p)%values(:,i,j,k))
      enddo
      call MPI_Allreduce(local_bound,lower_bound,nv,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      local_bound=-huge(1._dp)
      do r=1,nc
        p=cells(1,r);i=cells(2,r);j=cells(3,r);k=cells(4,r)
        if(owners(p)/=rank) cycle
        local_bound=max(local_bound,input%fields(p)%values(:,i,j,k))
      enddo
      call MPI_Allreduce(local_bound,upper_bound,nv,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS) return
      ! Fetch only the weighted means in stencils rooted on this owner.
      needed=.false.
      do r=1,nc
        if(owners(cells(1,r))/=rank) cycle
        if(depth>=0.and.tree%nodes(cells(1,r))%level/=depth) cycle
        call build_graph_stencil(r,last)
        needed(queue(:last))=.true.
      enddo
      call make_requests(lookup,condition)
      if(.not.condition) return
      call fetch_mpi_amr_tree_cells_3d(tree,owners,averages,packet_nv,requests,graph_values,comm,condition)
      if(.not.condition) return
      condition=.true.
      do r=1,nc
        if(owners(cells(1,r))/=rank) cycle
        if(depth>=0.and.tree%nodes(cells(1,r))%level/=depth) cycle
        call build_graph_stencil(r,last)
        stencil_value(:,:last)=graph_values(1:nv,lookup(queue(:last)))
        stencil_fraction(:last)=fraction_cell(queue(:last));stencil_fit(:,:last)=stencil_centroid(:,:last)
        extdir=.false.;extdir_bounds=0
        if(present(boundary_halos)) then
          p=cells(1,r);global_cell=int(origin(:,p))+cells(2:4,r)
          global_dims=int(domain(:,tree%nodes(p)%level+1))
          extdir=faces=='inflow';extdir_bounds(1,:)=1-global_cell;extdir_bounds(2,:)=global_dims-global_cell
          call state_redistribution_boundary_nodes_3d(boundary_geometry(p),boundary_halos%fields(p)%values, &
            cells(2:4,r)+2,global_cell,global_dims,faces,ghost_value,ghost_point,ghost_fit,ghost_center, &
            ghost_fraction,ghost_near,ghost_count,reconstructed)
          if(.not.reconstructed) then
            condition=.false.;cycle
          endif
          root_centroid=(center(:,r)-lattice(:,r))/width(:,r)
          do t=1,ghost_count
            last=last+1;stencil_value(:,last)=ghost_value(:,t)
            stencil_fraction(last)=ghost_fraction(t);stencil_connected(last)=ghost_near(t)
            stencil_centroid(:,last)=ghost_point(:,t)-root_centroid
            stencil_fit(:,last)=ghost_fit(:,t)-root_centroid
            stencil_center(:,last)=ghost_center(:,t);stencil_width(:,last)=1
          enddo
        endif
        call state_redistribution_graph_slopes_3d(stencil_value(:,:last),stencil_centroid(:,:last), &
          stencil_center(:,:last),stencil_width(:,:last),stencil_fraction(:last),selected_order,root_slopes,reconstructed, &
          near_connected=stencil_connected(:last),extdir_faces=extdir,extdir_bounds=extdir_bounds, &
          fit_centroid=stencil_fit(:,:last))
        if(.not.reconstructed) then
          condition=.false.;cycle
        endif
        p=cells(1,r);i=cells(2,r);j=cells(3,r);k=cells(4,r)
        if(fallback) then
          nmember=1
          if(allocated(groups(r)%member)) nmember=size(groups(r)%member)
          do t=1,nmember
            next=r
            if(allocated(groups(r)%member)) next=groups(r)%member(t)
            positions(:,t)=displacement(next,r)/width(:,r)
          enddo
        else
          nmember=1;positions(:,1)=0;next=receiver(r)
          if(next>0) then
            nmember=2;positions(:,2)=displacement(next,r)/width(:,r)
          endif
        endif
        shift=(centroid(:,r)-center(:,r))/width(:,r)
        call limit_state_recipient_slopes(averages%fields(p)%values(1:nv,i,j,k),shift,root_slopes, &
          positions(:,:nmember),lower_bound,upper_bound)
        root_slopes=root_slopes/spread(width(:,r),2,nv)
        averages%fields(p)%values(nv+1:packet_nv,i,j,k)=reshape(root_slopes,[3*nv])
      enddo
      call agree(condition,success)
    end subroutine reconstruct_owned_graph_slopes
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
      rate,comm,ok,target_volume_fraction,boundary_faces,max_order, &
      boundary_halos,boundary_geometry)
    type(amr_tree_state_3d), optional, intent(in) :: boundary_halos
    type(eb_geometry_3d), optional, intent(in) :: boundary_geometry(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nv
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,rhs
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    integer, optional, intent(in) :: max_order
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
      target_volume_fraction,boundary_faces,max_order=max_order, &
      boundary_halos=boundary_halos,boundary_geometry=boundary_geometry)
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
