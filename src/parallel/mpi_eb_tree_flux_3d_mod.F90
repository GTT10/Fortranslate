module mpi_eb_tree_flux_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use eb_held_funding_queue_3d_mod, only: pack_eb_held_funding_faces_3d,physical_eb_held_funding_debt_3d, &
    pack_eb_held_funding_increment_3d
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,periodic_eb_faces_match
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,valid_eb_integrated_flux_3d, &
    reconcile_reactive_eb_peer_flux_3d,restrict_eb_integrated_flux_3d,reflux_reactive_eb_state_3d, &
    compute_reactive_eb_reflux_increment_3d,reconcile_reactive_eb_interface_budget_3d, &
    supply_reactive_eb_interface_donors_3d,supply_reactive_eb_endpoint_donor_3d, &
    supply_reactive_eb_held_eos_donor_3d,supply_reactive_eb_gas_eos_donor_3d, &
    reactive_eb_funding_admission, &
    reactive_eb_exchange_interval_3d,add_eb_rereflux_to_fine_3d
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use eb_reactive_redistribution_3d_mod, only: state_redist_register_3d
  use mpi_eb_tree_redistribution_3d_mod, only: assemble_state_redist_mpi_interface_rereflux_3d
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d,average_down_mpi_eb_tree_3d, &
    average_down_mpi_eb_tree_field_3d,integrate_mpi_eb_tree_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d,fetch_mpi_amr_tree_cells_3d
  implicit none
  private
  public :: synchronize_mpi_eb_tree_flux_3d,assemble_mpi_eb_tree_reflux_3d
  public :: assemble_mpi_eb_held_funding_debt_3d
  public :: repair_mpi_eb_tree_interface_budget_3d,synchronize_mpi_eb_tree_register_flux_3d
  public :: repair_mpi_eb_tree_held_interface_budget_3d
  public :: synchronize_mpi_eb_tree_scalar_flux_3d,synchronize_mpi_eb_tree_scalar_register_flux_3d
contains
  subroutine agree(condition,comm,ok)
    logical, intent(in) :: condition
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: status
    call MPI_Allreduce(condition,ok,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if (status/=MPI_SUCCESS) ok=.false.
  end subroutine

  subroutine bounds(tree,node,lo,hi)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: node
    integer, intent(out) :: lo(3),hi(3)
    associate(p=>tree%nodes(node)%patch)
    lo=[p%coarse_i_lower,p%coarse_j_lower,p%coarse_k_lower]
    hi=[p%coarse_i_upper,p%coarse_j_upper,p%coarse_k_upper]
    end associate
  end subroutine

  ! Close exactly one coarse/fine source interval. The accepted flux core
  ! injects fine-covered corrections into descendants and jointly recovers
  ! EOS/average-down. Only then consume the two ledgers for this interface.
  subroutine synchronize_mpi_eb_tree_register_flux_3d(species,tree,owners,geometry,extended_geometry, &
      state,temperature,flux,registers,depth,comm,ok,boundary_faces,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),depth
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    type(state_redist_register_3d), intent(inout) :: registers(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d) :: correction,candidate,temps,spill
    type(eb_integrated_flux_3d), allocatable :: working(:)
    integer :: rank,status,p,nv,first_depth
    logical :: ready,has_spill
    ok=.false.;nv=reactive_nvar(size(species))
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call assemble_state_redist_mpi_interface_rereflux_3d(tree,owners,extended_geometry,registers,nv, &
      depth,correction,comm,ready,boundary_faces,spill,has_spill)
    if(.not.ready) return
    first_depth=depth
    if(has_spill) then
      first_depth=depth-1
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        if(tree%nodes(p)%level==first_depth) correction%fields(p)%values=spill%fields(p)%values
      enddo
    endif
    call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,state,temperature,flux,comm,ready, &
      depth,depth+1,boundary_faces=boundary_faces,frame=frame,rereflux_delta=correction, &
      rereflux_first_level=first_depth)
    if(.not.ready) then
      ! Keep the accepted state and live ledgers private during deficit repair.
      candidate=state;temps=temperature;working=flux
      call repair_mpi_eb_tree_interface_budget_3d(species,tree,owners,geometry,candidate,temps,working, &
        comm,ready,depth,boundary_faces,frame,rereflux_delta=correction,rereflux_first_level=first_depth)
      if(.not.ready) return
      call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,temps,working,comm,ready, &
        depth,depth+1,boundary_faces=boundary_faces,frame=frame,rereflux_delta=correction, &
        rereflux_first_level=first_depth)
      if(.not.ready) return
      call move_alloc(candidate%fields,state%fields);call move_alloc(temps%fields,temperature%fields)
      flux=working
    end if
    do p=1,size(tree%nodes)
      if(owners(p)/=rank) cycle
      if(tree%nodes(p)%level==depth) then
        registers(p)%coarse_delta=0
      else if(tree%nodes(p)%level==depth+1) then
        registers(p)%fine_delta=0
      endif
    enddo
    ok=.true.
  end subroutine synchronize_mpi_eb_tree_register_flux_3d

  ! Match serial peer ordering, deepest-first replacement, and exposed-cell reflux.
  ! State, temperature and flux are owner-local; geometry is replicated.
  ! All three outputs commit only after collective success. The peers-only mode
  ! corrects provisional endpoints without restriction or descendant averaging.
  subroutine synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,state,temperature,flux,comm,ok, &
      minimum_level,maximum_level,peers_only,boundary_faces,frame,rereflux_delta,rereflux_first_level)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_level,maximum_level
    logical, optional, intent(in) :: peers_only
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    integer, optional, intent(in) :: rereflux_first_level
    call synchronize_mpi_eb_tree_flux_core_3d(species,tree,owners,geometry,state,temperature,flux,comm,ok, &
      minimum_level,maximum_level,peers_only,boundary_faces,frame,rereflux_delta=rereflux_delta, &
      rereflux_first_level=rereflux_first_level)
  end subroutine synchronize_mpi_eb_tree_flux_3d

  ! Assemble owner-local signed reflux increments and canonical replacement
  ! faces without advancing, recovering or averaging any gas state. A debt
  ! remains signed for the composite correction; these are not accepted states.
  ! Both outputs publish only on global success; inactive levels have zero
  ! increments and retain their original faces, including idle-rank ownership.
  ! Assemble packed gas/actual-H debt without applying it to physical state.
  ! All local candidates stay private until every owner and idle rank agrees.
  subroutine assemble_mpi_eb_held_funding_debt_3d(species,tree,owners,geometry,state,temperature,stored, &
      flux,scalar_flux,clock,debt,comm,ok,boundary_faces,frame,rereflux_delta,scalar_rereflux_delta)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),clock
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature,stored
    type(eb_integrated_flux_3d), intent(in) :: flux(:),scalar_flux(:)
    type(amr_tree_state_3d), intent(inout) :: debt
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta,scalar_rereflux_delta
    type(amr_tree_state_3d), allocatable :: joint_rereflux
    type(amr_tree_state_3d) :: packed,raw,physical
    type(eb_integrated_flux_3d), allocatable :: packed_faces(:),canonical(:)
    logical, allocatable :: local(:)
    integer :: nv,p,n(3),rank,status,allocation_status,low,high,controls(2),control_low(2),control_high(2)
    logical :: good,local_ok
    ok=.false.;nv=reactive_nvar(size(species))
    call MPI_Allreduce(clock,low,1,MPI_INTEGER,MPI_MIN,comm,status)
    call agree(status==MPI_SUCCESS,comm,good);if(.not.good) return
    call MPI_Allreduce(clock,high,1,MPI_INTEGER,MPI_MAX,comm,status)
    call agree(status==MPI_SUCCESS,comm,good);if(.not.good) return
    if(low/=high) return
    controls=[merge(1,0,present(rereflux_delta)),merge(1,0,present(scalar_rereflux_delta))]
    call MPI_Allreduce(controls,control_low,2,MPI_INTEGER,MPI_MIN,comm,status)
    call agree(status==MPI_SUCCESS,comm,good);if(.not.good) return
    call MPI_Allreduce(controls,control_high,2,MPI_INTEGER,MPI_MAX,comm,status)
    call agree(status==MPI_SUCCESS,comm,good);if(.not.good) return
    if(any(control_low/=control_high).or.control_low(1)/=control_low(2)) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,good)
    if(.not.good) return
    call agree(clock>=0.and.clock<tree%level_count(),comm,good)
    if(.not.good) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,good)
    if(.not.good) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,good)
    if(.not.good) return
    call validate_mpi_amr_tree_local_3d(tree,owners,stored,1,comm,good)
    if(.not.good) return
    call MPI_Comm_rank(comm,rank,status)
    call agree(status==MPI_SUCCESS,comm,good);if(.not.good) return
    allocate(local(size(tree%nodes)),packed%fields(size(tree%nodes)),stat=allocation_status)
    call agree(allocation_status==0,comm,good);if(.not.good) return
    local=owners==rank
    if(present(rereflux_delta)) then
      allocate(joint_rereflux,stat=allocation_status)
      call agree(allocation_status==0,comm,good);if(.not.good) return
      call pack_eb_held_funding_increment_3d(tree,geometry,nv,local,rereflux_delta, &
        scalar_rereflux_delta,joint_rereflux,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
    endif
    call pack_eb_held_funding_faces_3d(tree,geometry,nv,local,flux,scalar_flux,packed_faces,local_ok)
    call agree(local_ok,comm,good);if(.not.good) return
    local_ok=.true.
    do p=1,size(tree%nodes)
      if(.not.local(p)) cycle
      n=tree%nodes(p)%shape
      allocate(packed%fields(p)%values(nv+1,n(1),n(2),n(3)),stat=allocation_status)
      if(allocation_status/=0) then
        local_ok=.false.;exit
      endif
      packed%fields(p)%values(:nv,:,:,:)=state%fields(p)%values
      packed%fields(p)%values(nv+1:nv+1,:,:,:)=stored%fields(p)%values
    enddo
    call agree(local_ok,comm,good);if(.not.good) return
    call assemble_mpi_eb_tree_reflux_3d(species,tree,owners,geometry,packed,temperature,packed_faces,raw,canonical, &
      comm,good,0,clock,boundary_faces,frame,rereflux_delta=joint_rereflux,passive_components=1)
    if(.not.good) return
    call physical_eb_held_funding_debt_3d(tree,geometry,nv,local,clock,raw,physical,local_ok)
    call agree(local_ok,comm,good);if(.not.good) return
    call move_alloc(physical%fields,debt%fields);ok=.true.
  end subroutine assemble_mpi_eb_held_funding_debt_3d

  subroutine assemble_mpi_eb_tree_reflux_3d(species,tree,owners,geometry,state,temperature,flux,increment,replacement, &
      comm,ok,minimum_level,maximum_level,boundary_faces,frame,rereflux_delta,rereflux_first_level, &
      passive_components)
    integer, optional, intent(in) :: passive_components
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(eb_integrated_flux_3d), intent(in) :: flux(:)
    type(amr_tree_state_3d), intent(inout) :: increment
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: replacement(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_level,maximum_level
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    type(amr_tree_state_3d) :: private_state,private_temperature,private_increment
    type(eb_integrated_flux_3d), allocatable :: private_flux(:)
    integer, optional, intent(in) :: rereflux_first_level
    private_state=state;private_temperature=temperature;private_flux=flux
    call synchronize_mpi_eb_tree_flux_core_3d(species,tree,owners,geometry,private_state,private_temperature, &
      private_flux,comm,ok,minimum_level,maximum_level,boundary_faces=boundary_faces,frame=frame, &
      reflux_increment=private_increment,rereflux_delta=rereflux_delta,rereflux_first_level=rereflux_first_level, &
        passive_components=passive_components)
    if(.not.ok) return
    call move_alloc(private_increment%fields,increment%fields)
    call move_alloc(private_flux,replacement)
  end subroutine assemble_mpi_eb_tree_reflux_3d

  ! Repair an adjacent-level coarse/fine deficit before the accepted reflux.
  ! Owner-local state/flux publication is atomic. Only fine boundary states
  ! and faces are changed here; coarse states receive the matching replacement
  ! through the subsequent strict synchronization. No species are clipped.
  ! Connected interior donors fund empty interface intervals. Insufficient
  ! inventories or unsupported paths reject the entire transaction.
  subroutine repair_mpi_eb_tree_interface_budget_3d(species,tree,owners,geometry,state,temperature,flux,comm,ok, &
      minimum_level,boundary_faces,frame,rereflux_delta,rereflux_first_level)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),minimum_level
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    integer, optional, intent(in) :: rereflux_first_level
    integer :: register_first_level,controls(3),low(3),high(3),status
    ok=.false.
    register_first_level=minimum_level
    if(present(rereflux_first_level)) register_first_level=rereflux_first_level
    controls=[minimum_level,register_first_level, &
      merge(1,0,present(rereflux_delta))+2*merge(1,0,present(rereflux_first_level))]
    call MPI_Allreduce(controls,low,3,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,3,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    if(register_first_level<0.or.register_first_level>minimum_level) return
    call repair_mpi_eb_tree_interface_components_3d(species,tree,owners,geometry,state,temperature,flux,comm,ok, &
      minimum_level,boundary_faces,frame,rereflux_delta,0,register_first_level)
  end subroutine repair_mpi_eb_tree_interface_budget_3d

  ! Repair paired gas/UEINT inventories with the same donor packets and face
  ! coefficients. The scalar stays independent of UEDEN during every EOS query.
  ! Only a successful complete repair publishes either family.
  subroutine repair_mpi_eb_tree_held_interface_budget_3d(species,tree,owners,geometry,state,temperature, &
      internal_energy,flux,internal_energy_flux,comm,ok,minimum_level,boundary_faces,frame, &
      rereflux_delta,internal_energy_rereflux_delta, &
      rereflux_first_level,reactor_energy_bias,native_projection,check_native_endpoints)
    procedure(reactive_eb_funding_admission), optional :: native_projection
    logical, optional, intent(in) :: check_native_endpoints
    logical :: native_endpoints
    type(amr_tree_state_3d), optional, intent(in) :: reactor_energy_bias
    integer, optional, intent(in) :: rereflux_first_level
    integer :: register_first_level
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),minimum_level
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature,internal_energy
    type(eb_integrated_flux_3d), intent(inout) :: flux(:),internal_energy_flux(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta,internal_energy_rereflux_delta
    type(amr_tree_state_3d) :: packed,temps,gas_work,held_work
    type(amr_tree_state_3d), allocatable :: packed_delta
    type(eb_integrated_flux_3d), allocatable :: packets(:),gas_faces(:),held_faces(:)
    integer :: nv,p,rank,status,ios,n(3),flags,flag_min,flag_max
    integer :: first_min,first_max
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    native_endpoints=.false.
    if(present(check_native_endpoints)) native_endpoints=check_native_endpoints
    flags=merge(1,0,present(rereflux_delta))+2*merge(1,0,present(internal_energy_rereflux_delta))+ &
      4*merge(1,0,present(rereflux_first_level))+8*merge(1,0,present(reactor_energy_bias))+ &
      16*merge(1,0,present(native_projection))+32*merge(1,0,native_endpoints)
    call MPI_Allreduce(flags,flag_min,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,flag_max,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.flag_min/=flag_max) return
    if(native_endpoints.and..not.present(native_projection)) return
    if(iand(flags,3)/=0.and.iand(flags,3)/=3) return
    register_first_level=minimum_level
    if(present(rereflux_first_level)) register_first_level=rereflux_first_level
    call MPI_Allreduce(register_first_level,first_min,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(register_first_level,first_max,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.first_min/=first_max) return
    if(register_first_level<0.or.register_first_level>minimum_level) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,ready)
    if(.not.ready) return
    if(present(reactor_energy_bias)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,reactor_energy_bias,1,comm,ready)
      if(.not.ready) return
    endif
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call pack_fields(state,internal_energy,packed,ready)
    if(.not.ready) return
    if(present(rereflux_delta)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,rereflux_delta,nv,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy_rereflux_delta,1,comm,ready)
      if(.not.ready) return
      allocate(packed_delta,stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      call pack_fields(rereflux_delta,internal_energy_rereflux_delta,packed_delta,ready)
      if(.not.ready) return
    endif
    valid=size(flux)==size(owners).and.size(internal_energy_flux)==size(owners)
    call agree(valid,comm,ready)
    if(.not.ready) return
    allocate(packets(size(owners)),stat=ios)
    call agree(ios==0,comm,ready)
    if(.not.ready) return
    valid=.true.
    do p=1,size(owners)
      if(owners(p)/=rank) then
        valid=.not.(allocated(flux(p)%x).or.allocated(flux(p)%y).or.allocated(flux(p)%z).or. &
          allocated(internal_energy_flux(p)%x).or.allocated(internal_energy_flux(p)%y).or. &
          allocated(internal_energy_flux(p)%z))
        if(.not.valid) exit
        cycle
      endif
      valid=valid_eb_integrated_flux_3d(flux(p),geometry(p),nv).and. &
        valid_eb_integrated_flux_3d(internal_energy_flux(p),geometry(p),1)
      if(.not.valid) exit
      n=tree%nodes(p)%shape
      allocate(packets(p)%x(nv+1,0:n(1),n(2),n(3)),packets(p)%y(nv+1,n(1),0:n(2),n(3)), &
        packets(p)%z(nv+1,n(1),n(2),0:n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      endif
      packets(p)%x(:nv,:,:,:)=flux(p)%x;packets(p)%x(nv+1:nv+1,:,:,:)=internal_energy_flux(p)%x
      packets(p)%y(:nv,:,:,:)=flux(p)%y;packets(p)%y(nv+1:nv+1,:,:,:)=internal_energy_flux(p)%y
      packets(p)%z(:nv,:,:,:)=flux(p)%z;packets(p)%z(nv+1:nv+1,:,:,:)=internal_energy_flux(p)%z
    enddo
    call agree(valid,comm,ready)
    if(.not.ready) return
    temps=temperature
    call repair_mpi_eb_tree_interface_components_3d(species,tree,owners,geometry,packed,temps,packets,comm,ready, &
      minimum_level,boundary_faces,frame,packed_delta,1,register_first_level,reactor_energy_bias,native_projection, &
      check_native_endpoints)
    if(.not.ready) return
    gas_work=state;held_work=internal_energy;gas_faces=flux;held_faces=internal_energy_flux
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      gas_work%fields(p)%values=packed%fields(p)%values(:nv,:,:,:)
      held_work%fields(p)%values=packed%fields(p)%values(nv+1:nv+1,:,:,:)
      gas_faces(p)%x=packets(p)%x(:nv,:,:,:);held_faces(p)%x=packets(p)%x(nv+1:nv+1,:,:,:)
      gas_faces(p)%y=packets(p)%y(:nv,:,:,:);held_faces(p)%y=packets(p)%y(nv+1:nv+1,:,:,:)
      gas_faces(p)%z=packets(p)%z(:nv,:,:,:);held_faces(p)%z=packets(p)%z(nv+1:nv+1,:,:,:)
    enddo
    call move_alloc(gas_work%fields,state%fields);call move_alloc(held_work%fields,internal_energy%fields)
    call move_alloc(temps%fields,temperature%fields)
    flux=gas_faces;internal_energy_flux=held_faces;ok=.true.
  contains
    subroutine pack_fields(gas,held,result,good)
      type(amr_tree_state_3d), intent(in) :: gas,held
      type(amr_tree_state_3d), intent(out) :: result
      logical, intent(out) :: good
      integer :: p,n(3),allocation_status
      logical :: local_ok
      allocate(result%fields(size(owners)),stat=allocation_status)
      call agree(allocation_status==0,comm,good)
      if(.not.good) return
      local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        local_ok=all(ieee_is_finite(gas%fields(p)%values)).and.all(ieee_is_finite(held%fields(p)%values))
        if(.not.local_ok) exit
        n=tree%nodes(p)%shape
        allocate(result%fields(p)%values(nv+1,n(1),n(2),n(3)),stat=allocation_status)
        if(allocation_status/=0) then
          local_ok=.false.;exit
        endif
        result%fields(p)%values(:nv,:,:,:)=gas%fields(p)%values
        result%fields(p)%values(nv+1:nv+1,:,:,:)=held%fields(p)%values
      enddo
      call agree(local_ok,comm,good)
    end subroutine pack_fields
  end subroutine repair_mpi_eb_tree_held_interface_budget_3d

  subroutine repair_mpi_eb_tree_interface_components_3d(species,tree,owners,geometry,state,temperature,flux,comm,ok, &
      minimum_level,boundary_faces,frame,rereflux_delta, &
      passive_components,register_first_level,reactor_energy_bias,native_projection,check_native_endpoints)
    procedure(reactive_eb_funding_admission), optional :: native_projection
    logical, optional, intent(in) :: check_native_endpoints
    logical :: native_endpoints
    type(amr_tree_state_3d), optional, intent(in) :: reactor_energy_bias
    integer, intent(in) :: passive_components,register_first_level
    integer :: active_components
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),minimum_level
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    type(amr_tree_state_3d) :: candidate,temps,increment,routed_register,endpoint_budget,bias_view
    type(eb_integrated_flux_3d), allocatable :: effective(:),working(:)
    real(dp), allocatable :: fine_state(:,:),volumes(:),transfers(:,:),weights(:),areas(:),old(:),base(:),endpoint(:)
    real(dp), allocatable :: corrected(:,:),dc(:),df(:,:),primitive(:)
    real(dp), allocatable :: donor_state(:,:),donor_volume(:),funded_fine(:,:),funded_donor(:,:),exchange(:,:)
    real(dp), allocatable :: packet_values(:,:),low_packet(:),high_packet(:)
    integer, allocatable :: packet_keys(:,:)
    integer :: packet_count,packet_capacity,addition
    integer, allocatable :: used_donors(:,:),next_used_donors(:,:)
    integer, allocatable :: route_node(:),route_cell(:,:),route_axis(:),route_step(:)
    integer :: route_length,register_flag,register_low,register_high
    integer, allocatable :: positions(:,:),ratios(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    logical, allocatable :: changed(:)
    integer :: nv,rank,status,ios,child,parent,other,axis,side,ratio,lo(3),hi(3),n(3),first(3),last(3)
    integer :: i,j,k,x,y,z,q,nf,capacity,face_index(3),coarse_cell(3),fine_cell(3),fine_face(3),fl(3),fh(3)
    integer :: lower(3),upper(3),coarse_owner,fine_owner,p,donor_cell(3),internal_face(3)
    integer :: order(3),search_axis,direction,distance,path,donor_axis,donor_step,r,neighbor(3),cursor(3),mirror(3)
    integer :: cursor_node,neighbor_node,donor_node,selected,ancestor_level,ancestor_start
    real(dp) :: coarse_volume,coarse_area,orientation,theta,recovered,speed,maximum_correction,theta_floor
    real(dp) :: coarse_guess
    real(dp) :: coarse_lower,coarse_limit,required_correction,composite_budget(reactive_nvar(size(species))+passive_components)
    logical :: ready,valid,unsafe,covered,local_ok,found,connected_route,periodic(3)
    logical :: route_improved,route_complete,route_can_reroute,immediate_funding
    ok=.false.;immediate_funding=present(native_projection)
    native_endpoints=.false.
    if(present(check_native_endpoints)) native_endpoints=check_native_endpoints
    if(present(reactor_energy_bias).and.passive_components/=1) return
    if(passive_components<0.or.passive_components>1) return
    active_components=reactive_nvar(size(species));nv=active_components+passive_components
    register_flag=merge(1,0,present(rereflux_delta))+2*merge(1,0,present(reactor_energy_bias))+ &
      4*merge(1,0,immediate_funding)+8*merge(1,0,native_endpoints)
    call MPI_Allreduce(register_flag,register_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(register_flag,register_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.register_low/=register_high) return
    if(native_endpoints.and..not.immediate_funding) return
    if(immediate_funding.and..not.present(rereflux_delta)) return
    if(present(rereflux_delta)) then
      call route_mpi_eb_tree_rereflux_3d(tree,owners,rereflux_delta,nv,register_first_level,minimum_level+1, &
        routed_register,comm,ready)
      if(.not.ready) return
    end if
    call assemble_mpi_eb_tree_reflux_3d(species,tree,owners,geometry,state,temperature,flux,increment,effective, &
      comm,ready,minimum_level,minimum_level+1,boundary_faces,frame, &
        passive_components=passive_components)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    if(immediate_funding) then
      ! Native source admission spends current gas/UEINT stock, without a
      ! future coarse/fine reflux credit before the donor clock is advanced.
      do p=1,size(owners)
        if(owners(p)==rank) increment%fields(p)%values=0
      enddo
      effective=flux
    endif
    if(present(reactor_energy_bias)) then
      ! Recovery also runs in owner-only loops. Replicate this read-only
      ! predicate field once in a rank-uniform order, never from a nested EOS
      ! call; mutable state/flux ownership remains unchanged.
      allocate(bias_view%fields(size(owners)),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      local_ok=.true.
      do p=1,size(owners)
        n=tree%nodes(p)%shape
        allocate(bias_view%fields(p)%values(1,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) then
          local_ok=.false.;cycle
        endif
        if(owners(p)==rank) then
          bias_view%fields(p)%values=reactor_energy_bias%fields(p)%values
          local_ok=local_ok.and.all(ieee_is_finite(bias_view%fields(p)%values))
        endif
      enddo
      call agree(local_ok,comm,ready)
      if(.not.ready) return
      do p=1,size(owners)
        n=tree%nodes(p)%shape
        call MPI_Bcast(bias_view%fields(p)%values,product(n),MPI_DOUBLE_PRECISION,owners(p),comm,status)
        if(status/=MPI_SUCCESS) return
      enddo
    endif
    block
      ! Include ordinary reflux as well as redistribution-register debt. A
      ! negative composite inventory cannot be funded by conservative packets.
      endpoint_budget=state;local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        n=tree%nodes(p)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          call add_mpi_eb_delta_cell_3d(endpoint_budget%fields(p)%values(:,i,j,k), &
            increment%fields(p)%values(:,i,j,k),valid)
          if(valid.and.present(rereflux_delta)) call add_mpi_eb_delta_cell_3d( &
            endpoint_budget%fields(p)%values(:,i,j,k),routed_register%fields(p)%values(:,i,j,k),valid)
          local_ok=local_ok.and.valid
        end do
        end do
        end do
      end do
      call agree(local_ok,comm,ready)
      if(.not.ready) return
      call integrate_mpi_eb_tree_3d(tree,owners,geometry,endpoint_budget,composite_budget,comm,ready)
      if(.not.ready) return
      if(composite_budget(1)<=0.or.any(composite_budget(6:active_components)<0)) return
    end block
    periodic=.false.
    if(present(boundary_faces)) periodic=[boundary_faces(1)=='periodic',boundary_faces(3)=='periodic', &
      boundary_faces(5)=='periodic']
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,comm,ready)
    if(.not.ready) return
    candidate=state;temps=temperature;working=flux
    allocate(old(nv),base(nv),endpoint(nv),low_packet(nv),high_packet(nv), &
      primitive(reactive_nprim(size(species))),changed(size(owners)),stat=ios)
    call agree(ios==0,comm,ready)
    if(.not.ready) return
    changed=.false.
    ! Apply the pending peer-face divergence to the private fine baseline before
    ! spending donor inventories. Publishing canonical fine faces prevents the
    ! subsequent strict synchronization from applying that divergence twice.
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      if(tree%nodes(p)%level/=minimum_level+1) cycle
      candidate%fields(p)%values=candidate%fields(p)%values+increment%fields(p)%values
      changed(p)=any(increment%fields(p)%values/=0)
      working(p)=effective(p)
    enddo
    if(present(rereflux_delta)) then
      ! Donor inventories include the pending source register once. Descendant
      ! receivers are routed before applying the same conserved flux packets.
      local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or.tree%nodes(p)%level<=minimum_level) cycle
        n=tree%nodes(p)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          call add_mpi_eb_delta_cell_3d(candidate%fields(p)%values(:,i,j,k), &
            routed_register%fields(p)%values(:,i,j,k),valid)
          local_ok=local_ok.and.valid
        end do
        end do
        end do
        changed(p)=changed(p).or.any(routed_register%fields(p)%values/=0)
      end do
      call agree(local_ok,comm,ready)
      if(.not.ready) return
    end if
    if(.not.immediate_funding) then
    do child=1,size(owners)
      if(tree%nodes(child)%level/=minimum_level+1) cycle
      parent=tree%nodes(child)%parent
      ratio=tree%nodes(child)%patch%refinement_ratio
      capacity=ratio*ratio;n=tree%nodes(parent)%shape
      call bounds(tree,child,lo,hi)
      allocate(fine_state(nv,capacity),volumes(capacity),transfers(nv,capacity),weights(capacity),areas(capacity), &
        positions(3,capacity), &
        donor_state(nv,capacity),donor_volume(capacity),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      coarse_owner=owners(parent);fine_owner=owners(child)
      do axis=1,3
      do side=0,1
        orientation=real(2*side-1,dp)
        first=lo;last=hi
        first(axis)=merge(lo(axis)-1,hi(axis)+1,side==0);last(axis)=first(axis)
        do k=first(3),last(3)
        do j=first(2),last(2)
        do i=first(1),last(1)
          coarse_cell=[i,j,k]
          if(any(coarse_cell<1).or.any(coarse_cell>n)) cycle
          if(geometry(parent)%volume_fraction(i,j,k)<=0) cycle
          covered=.false.
          do other=1,size(owners)
            if(tree%nodes(other)%parent/=parent) cycle
            call bounds(tree,other,lower,upper)
            if(all(coarse_cell>=lower).and.all(coarse_cell<=upper)) covered=.true.
          enddo
          if(covered) cycle
          face_index=coarse_cell;face_index(axis)=face_index(axis)-side
          coarse_volume=geometry(parent)%volume_fraction(i,j,k)* &
            geometry(parent)%dx*geometry(parent)%dy*geometry(parent)%dz
          coarse_area=open_area(parent,axis,face_index)
          if(coarse_area<=0) cycle
          unsafe=.false.
          if(rank==coarse_owner) then
            call cell_increment(parent,coarse_cell,[0,0,0],0,endpoint)
            endpoint=candidate%fields(parent)%values(:,i,j,k)+endpoint
            unsafe=any(endpoint(6:active_components)<0).or.endpoint(1)<=0
          endif
          call MPI_Bcast(unsafe,1,MPI_LOGICAL,coarse_owner,comm,status)
          if(status/=MPI_SUCCESS) return
          if(.not.unsafe) cycle
          if(rank==coarse_owner) then
            call cell_increment(parent,coarse_cell,face_index,axis,base)
            base=candidate%fields(parent)%values(:,i,j,k)+base
            old=orientation*coarse_area*face_value(working(parent),axis,face_index)
            if(present(rereflux_delta)) then
              ! Encode the coarse register as an offset of this face's old
              ! extensive transfer. The provisional base remains physical;
              ! base + fine - old still equals the complete signed endpoint.
              old=old-coarse_volume*routed_register%fields(parent)%values(:,i,j,k)
            end if
          endif
          call MPI_Bcast(base,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,status)
          if(status/=MPI_SUCCESS) return
          call MPI_Bcast(old,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,status)
          if(status/=MPI_SUCCESS) return
          fl=(coarse_cell-lo)*ratio+1;fh=fl+ratio-1
          fl(axis)=merge(1,tree%nodes(child)%shape(axis),side==0);fh(axis)=fl(axis)
          nf=0
          do z=fl(3),fh(3)
          do y=fl(2),fh(2)
          do x=fl(1),fh(1)
            fine_cell=[x,y,z];fine_face=fine_cell;fine_face(axis)=fine_face(axis)-(1-side)
            if(geometry(child)%volume_fraction(x,y,z)<=0) cycle
            if(open_area(child,axis,fine_face)<=0) cycle
            nf=nf+1;positions(:,nf)=fine_cell
            volumes(nf)=geometry(child)%volume_fraction(x,y,z)*geometry(child)%dx*geometry(child)%dy*geometry(child)%dz
            areas(nf)=open_area(child,axis,fine_face);weights(nf)=areas(nf)/coarse_area
            if(rank==fine_owner) then
              fine_state(:,nf)=candidate%fields(child)%values(:,x,y,z)
              transfers(:,nf)=orientation*areas(nf)*face_value(effective(child),axis,fine_face)
            endif
          enddo
          enddo
          enddo
          valid=nf>0
          call agree(valid,comm,ready)
          if(.not.ready) return
          call MPI_Bcast(fine_state(:,:nf),nv*nf,MPI_DOUBLE_PRECISION,fine_owner,comm,status)
          if(status/=MPI_SUCCESS) return
          call MPI_Bcast(transfers(:,:nf),nv*nf,MPI_DOUBLE_PRECISION,fine_owner,comm,status)
          if(status/=MPI_SUCCESS) return
          theta_floor=0
          do q=1,nf
            maximum_correction=1;low_packet=0;high_packet=0
            if(side==0) low_packet=transfers(:,q)-weights(q)*old
            if(side==1) high_packet=weights(q)*old-transfers(:,q)
            call propagate_cell_packet(child,positions(:,q),axis,low_packet,high_packet,.true.,valid,maximum_correction)
            call agree(valid,comm,ready)
            if(.not.ready) return
            theta_floor=max(theta_floor,1-maximum_correction)
          enddo
          theta=1
          call reconcile_reactive_eb_interface_budget_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
            transfers(:,:nf),weights(:nf),corrected,dc,df,theta,valid,theta_floor, &
        passive_components=passive_components)
          call agree(valid,comm,ready)
          if(.not.ready) then
            call fund_exposed_coarse_endpoint(valid)
            if(valid) then
              if(rank==coarse_owner) then
                call cell_increment(parent,coarse_cell,face_index,axis,base)
                base=candidate%fields(parent)%values(:,i,j,k)+base
              endif
              call MPI_Bcast(base,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,status)
              if(status/=MPI_SUCCESS) return
              call reconcile_reactive_eb_interface_budget_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
                transfers(:,:nf),weights(:nf),corrected,dc,df,theta,valid,theta_floor, &
        passive_components=passive_components)
              call agree(valid,comm,ready)
            endif
          endif
          if(.not.ready) then
            ! Each subface chooses its own connected path. A shared donor
            ! is debited immediately in the private state, so later searches
            ! see its remaining inventory. Any failure rolls back the group
            ! and every earlier group before public arrays are published.
            coarse_lower=0;coarse_limit=1
            call reactive_eb_exchange_interval_3d(reshape(base,[nv,1]), &
              reshape((sum(transfers(:,:nf),dim=2)-old)/coarse_volume,[nv,1]),coarse_lower,coarse_limit,valid, &
        passive_components=passive_components)
            if(.not.valid.or.coarse_lower>0) return
            required_correction=1-coarse_limit*(1-32*epsilon(1._dp))
            order=[axis,mod(axis,3)+1,mod(axis+1,3)+1]
            do selected=1,nf
              if(allocated(used_donors)) deallocate(used_donors)
              allocate(used_donors(4,0),stat=ios)
              call agree(ios==0,comm,ready)
              if(.not.ready) return
              funding: do
              donor_state(:,:nf)=fine_state(:,:nf);donor_volume(:nf)=volumes(:nf)
              call supply_reactive_eb_interface_donors_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
                transfers(:,:nf),weights(:nf),donor_state(:,:nf),donor_volume(:nf),funded_fine,funded_donor,exchange, &
                valid,target_subface=selected,force_exchange=covering_child(child,positions(:,selected))/=0, &
        passive_components=passive_components)
              call agree(valid,comm,ready)
              if(covering_child(child,positions(:,selected))/=0) then
                call recipient_allowance(selected,maximum_correction,valid)
                if(.not.valid) return
                ready=ready.and.maximum_correction>=required_correction
              endif
              if(ready) exit funding
              found=.false.;connected_route=.false.
              search: do search_axis=1,3
                donor_axis=order(search_axis)
                do direction=0,1
                  donor_step=2*direction-1
                  if(donor_axis==axis) donor_step=(2*direction-1)*(2*side-1)
                  do distance=1,int(min(domain(donor_axis,minimum_level+2)-1,int(huge(distance)-1,int64)))
                    cursor=positions(:,selected);cursor_node=child;valid=.true.
                    do path=1,distance
                      call path_neighbor(cursor_node,cursor,donor_axis,donor_step,neighbor_node,neighbor, &
                        internal_face,mirror,valid)
                      if(.not.valid) exit
                      cursor=neighbor;cursor_node=neighbor_node
                    enddo
                    if(.not.valid) cycle
                    donor_cell=cursor;donor_node=cursor_node
                    do r=1,nf
                      if(donor_node==child.and.all(donor_cell==positions(:,r))) valid=.false.
                    enddo
                    do r=1,size(used_donors,2)
                      if(all(used_donors(:,r)==[donor_node,donor_cell])) valid=.false.
                    enddo
                    if(.not.valid) cycle
                    x=donor_cell(1);y=donor_cell(2);z=donor_cell(3)
                    donor_volume(selected)=cell_volume(donor_node,donor_cell)
                    if(rank==owners(donor_node)) donor_state(:,selected)=candidate%fields(donor_node)%values(:,x,y,z)
                    call MPI_Bcast(donor_state(:,selected),nv,MPI_DOUBLE_PRECISION,owners(donor_node),comm,status)
                    if(status/=MPI_SUCCESS) return
                    call supply_reactive_eb_interface_donors_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
                      transfers(:,:nf),weights(:nf),donor_state(:,:nf),donor_volume(:nf),funded_fine,funded_donor,exchange, &
                      valid,target_subface=selected,force_exchange=covering_child(child,positions(:,selected))/=0, &
        passive_components=passive_components)
                    call agree(valid,comm,ready)
                    if(.not.ready) cycle
                    call build_straight_route(selected,distance,donor_axis,donor_step,valid)
                    if(.not.valid) return
                    call route_accepts_packet(selected,valid)
                    if(.not.valid.or..not.route_improved) cycle
                    if(route_complete) call minimize_recipient_packet(selected,valid)
                    if(.not.valid) cycle
                    connected_route=.true.;found=.true.;exit search
                  enddo
                enddo
              enddo search
              if(.not.found) then
                call search_connected_donor(selected,found,valid)
                if(.not.valid) return
                if(.not.found) call search_connected_donor(selected,found,valid,allow_partial=.true.)
                if(.not.valid.or..not.found) return
                connected_route=.true.
              endif
              allocate(next_used_donors(4,size(used_donors,2)+1),stat=ios)
              call agree(ios==0,comm,ready)
              if(.not.ready) return
              next_used_donors(:,:size(used_donors,2))=used_donors
              next_used_donors(:,size(next_used_donors,2))=[donor_node,donor_cell]
              call move_alloc(next_used_donors,used_donors)
              fine_state(:,:nf)=funded_fine;fine_cell=positions(:,selected)
              x=donor_cell(1);y=donor_cell(2);z=donor_cell(3)
              if(rank==owners(donor_node).and.covering_child(donor_node,donor_cell)==0) then
                candidate%fields(donor_node)%values(:,x,y,z)=funded_donor(:,selected)
                changed(donor_node)=.true.
              endif
              if(rank==fine_owner.and.covering_child(child,fine_cell)==0) then
                candidate%fields(child)%values(:,fine_cell(1),fine_cell(2),fine_cell(3))=funded_fine(:,selected)
                changed(child)=.true.
              endif
              cursor=fine_cell;cursor_node=child
              route_length=distance
              if(connected_route) then
                route_length=size(route_node)
              else
                if(allocated(route_node)) deallocate(route_node,route_cell,route_axis,route_step)
                allocate(route_node(route_length),route_cell(3,route_length),route_axis(route_length), &
                  route_step(route_length),stat=ios)
                call agree(ios==0,comm,ready)
                if(.not.ready) return
              endif
              do path=1,route_length
                if(connected_route) then
                  cursor_node=route_node(path);cursor=route_cell(:,path)
                  donor_axis=route_axis(path);donor_step=route_step(path)
                else
                  route_node(path)=cursor_node;route_cell(:,path)=cursor
                  route_axis(path)=donor_axis;route_step(path)=donor_step
                endif
                call path_neighbor(cursor_node,cursor,donor_axis,donor_step,neighbor_node,neighbor, &
                  internal_face,mirror,valid)
                if(.not.valid) return
                if(rank==owners(cursor_node)) then
                  call add_path_packet(cursor_node,donor_axis,internal_face,-real(donor_step,dp)*exchange(:,selected))
                endif
                if(neighbor_node/=cursor_node.or.any(mirror/=internal_face)) then
                  if(rank==owners(neighbor_node)) then
                    call add_path_packet(neighbor_node,donor_axis,mirror,-real(donor_step,dp)*exchange(:,selected))
                  endif
                endif
                cursor=neighbor;cursor_node=neighbor_node
              enddo
              call propagate_route_interior(exchange(:,selected),valid)
              if(.not.valid) return
              enddo funding
            enddo
            theta_floor=0
            do q=1,nf
              maximum_correction=1;low_packet=0;high_packet=0
              if(side==0) low_packet=transfers(:,q)-weights(q)*old
              if(side==1) high_packet=weights(q)*old-transfers(:,q)
              call propagate_cell_packet(child,positions(:,q),axis,low_packet,high_packet,.true.,valid,maximum_correction)
              call agree(valid,comm,ready)
              if(.not.ready) return
              theta_floor=max(theta_floor,1-maximum_correction)
            enddo
            call reconcile_reactive_eb_interface_budget_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
              transfers(:,:nf),weights(:nf),corrected,dc,df,theta,valid,theta_floor, &
        passive_components=passive_components)
            call agree(valid,comm,ready)
            if(.not.ready) return
          endif
          if(present(rereflux_delta)) then
            if(rank==coarse_owner) coarse_guess=temperature%fields(parent)%values(1,i,j,k)
            call MPI_Bcast(coarse_guess,1,MPI_DOUBLE_PRECISION,coarse_owner,comm,status)
            if(status/=MPI_SUCCESS) return
            call constrain_coarse_eos(valid)
            call agree(valid,comm,ready)
            if(.not.ready) return
          end if
          if(theta==1) cycle
          packet_capacity=0
          do q=1,nf
            addition=transfer_cell_count(child,positions(:,q))
            if(addition<0.or.addition>(huge(packet_capacity)-packet_capacity)/2) return
            packet_capacity=packet_capacity+2*addition
          enddo
          allocate(packet_keys(5,packet_capacity),packet_values(nv,packet_capacity),stat=ios)
          call agree(ios==0,comm,ready)
          if(.not.ready) return
          packet_count=0
          do q=1,nf
            fine_cell=positions(:,q);fine_face=fine_cell;fine_face(axis)=fine_face(axis)-(1-side)
            low_packet=0;high_packet=0
            if(side==0) low_packet=volumes(q)*df(:,q)
            if(side==1) high_packet=-volumes(q)*df(:,q)
            call propagate_cell_packet(child,fine_cell,axis,low_packet,high_packet,.true.,valid)
            call agree(valid,comm,ready)
            if(.not.ready) return
            if(rank==fine_owner) then
              call set_face(working(child),axis,fine_face,orientation*corrected(:,q)/areas(q))
              call set_face(effective(child),axis,fine_face,face_value(working(child),axis,fine_face))
            endif
          enddo
          deallocate(packet_keys,packet_values)
          if(rank==coarse_owner) call set_face(effective(parent),axis,face_index, &
            orientation*sum(corrected,dim=2)/coarse_area)
        enddo
        enddo
        enddo
      enddo
      enddo
      deallocate(fine_state,volumes,transfers,weights,areas,positions,donor_state,donor_volume)
    enddo
    endif
    ! StateRedist debt also reaches exposed cells away from fine interfaces
    ! on the closing level and older ancestors. The registered path scans all
    ! such endpoints; flux-only reconciliation keeps its established range.
    ! Candidate state remains pre-register and strict closure charges the
    ! original read-only debt once after these conservative gas/held packets.
    ancestor_start=minimum_level
    do ancestor_level=ancestor_start,register_first_level,-1
      do parent=1,size(owners)
        if(tree%nodes(parent)%level/=ancestor_level) cycle
        n=tree%nodes(parent)%shape
        coarse_owner=owners(parent)
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          coarse_cell=[i,j,k]
          if(geometry(parent)%volume_fraction(i,j,k)<=0) cycle
          if(covering_child(parent,coarse_cell)/=0) cycle
          unsafe=.false.
          if(rank==coarse_owner) then
            call cell_increment(parent,coarse_cell,[0,0,0],0,endpoint)
            endpoint=candidate%fields(parent)%values(:,i,j,k)+endpoint
            unsafe=endpoint(1)<=0.or.any(endpoint(6:active_components)<0)
            ! The new flux-only closing-level scan funds thermal EOS deficits;
            ! existing density/species reconciliation keeps its original range.
            if(.not.present(rereflux_delta).and.ancestor_level==minimum_level.and.unsafe) then
              unsafe=.false.
            else if(.not.unsafe) then
              call recover_repair_cell(species,geometry(parent),coarse_cell,endpoint, &
                temps%fields(parent)%values(1,i,j,k),primitive,recovered,speed,valid,frame, &
          repair_node=parent)
              unsafe=.not.valid
            endif
          endif
          call MPI_Bcast(unsafe,1,MPI_LOGICAL,coarse_owner,comm,status)
          if(status/=MPI_SUCCESS) return
          if(native_endpoints.and..not.unsafe) then
            ! The callback reconstructs the pending source once from the
            ! pre-register inventory. Every rank enters its collectives.
            if(rank==coarse_owner) endpoint=candidate%fields(parent)%values(:,i,j,k)
            call MPI_Bcast(endpoint,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,status)
            if(status/=MPI_SUCCESS) return
            call native_projection(parent,coarse_cell,endpoint,.true.,valid)
            call agree(valid,comm,ready)
            unsafe=.not.ready
          endif
          if(.not.unsafe) cycle
          coarse_volume=cell_volume(parent,coarse_cell)
          call fund_exposed_coarse_endpoint(valid)
          call agree(valid,comm,ready)
          if(.not.ready) return
        enddo
        enddo
        enddo
      enddo
    enddo
    ! Fine candidate inventories already contain routed register debt. Fund
    ! exposed physical leaves before EOS publication, retaining both pre/post
    ! states so strict closure can still charge the original debt once.
      do parent=1,size(tree%nodes)
        if(tree%nodes(parent)%level<=minimum_level) cycle
        n=tree%nodes(parent)%shape
        coarse_owner=owners(parent)
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          coarse_cell=[i,j,k]
          if(cell_volume(parent,coarse_cell)<=0.or.covering_child(parent,coarse_cell)/=0) cycle
          unsafe=.false.
          if(rank==coarse_owner) then
          ! Registered physical source debt can exhaust density or species.
          ! Raw zero-debt thermal funding still requires a physical endpoint.
          if(present(rereflux_delta)) then
            if(any(routed_register%fields(parent)%values(:,i,j,k)/=0)) &
              unsafe=candidate%fields(parent)%values(1,i,j,k)<=0.or. &
                any(candidate%fields(parent)%values(6:active_components,i,j,k)<0)
          endif
          if(candidate%fields(parent)%values(1,i,j,k)>0.and. &
            all(candidate%fields(parent)%values(6:active_components,i,j,k)>=0)) then
            call recover_repair_cell(species,geometry(parent),coarse_cell,candidate%fields(parent)%values(:,i,j,k), &
              temps%fields(parent)%values(1,i,j,k),primitive,recovered,speed,valid,frame, &
          repair_node=parent)
            unsafe=.not.valid
          endif
          endif
          call MPI_Bcast(unsafe,1,MPI_LOGICAL,coarse_owner,comm,status)
          if(status/=MPI_SUCCESS) return
          if(native_endpoints.and..not.unsafe) then
            ! The callback reconstructs the pending source once from the
            ! pre-register inventory. Every rank enters its collectives.
            valid=.true.
            if(rank==coarse_owner) then
              endpoint=candidate%fields(parent)%values(:,i,j,k)
              if(present(rereflux_delta)) call add_mpi_eb_delta_cell_3d(endpoint, &
                -routed_register%fields(parent)%values(:,i,j,k),valid)
            endif
            call agree(valid,comm,ready)
            if(.not.ready) return
            call MPI_Bcast(endpoint,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,status)
            if(status/=MPI_SUCCESS) return
            call native_projection(parent,coarse_cell,endpoint,.true.,valid)
            call agree(valid,comm,ready)
            unsafe=.not.ready
          endif
          if(.not.unsafe) cycle
          coarse_volume=cell_volume(parent,coarse_cell)
          call fund_exposed_fine_endpoint(valid)
          call agree(valid,comm,ready)
          if(.not.ready) return
        enddo
        enddo
        enddo
      enddo
    ! Verify every exposed coarse endpoint after all coupled face adjustments,
    ! and recover changed fine gas states before publishing anything.
    local_ok=.true.
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(geometry(p)%volume_fraction(i,j,k)<=0) cycle
        if(tree%nodes(p)%level>=register_first_level.and.tree%nodes(p)%level<=minimum_level) then
          covered=.false.
          do child=1,size(owners)
            if(tree%nodes(child)%parent/=p) cycle
            call bounds(tree,child,lo,hi)
            if(all([i,j,k]>=lo).and.all([i,j,k]<=hi)) covered=.true.
          enddo
          if(covered) cycle
          call cell_increment(p,[i,j,k],[0,0,0],0,endpoint)
          endpoint=candidate%fields(p)%values(:,i,j,k)+endpoint
          local_ok=local_ok.and.all(ieee_is_finite(endpoint)).and.endpoint(1)>0.and.all(endpoint(6:active_components)>=0)
          call recover_repair_cell(species,geometry(p),[i,j,k],endpoint, &
            temps%fields(p)%values(1,i,j,k),primitive,recovered,speed,valid,frame, &
          repair_node=p)
          local_ok=local_ok.and.valid
        endif
        if(.not.changed(p)) cycle
        if(any(candidate%fields(p)%values(6:active_components,i,j,k)<0)) then
          local_ok=.false.;cycle
        endif
        call recover_repair_cell(species,geometry(p),[i,j,k],candidate%fields(p)%values(:,i,j,k), &
          temperature%fields(p)%values(1,i,j,k),primitive,recovered,speed,valid,frame, &
          repair_node=p)
        if(.not.valid) then
          local_ok=.false.;cycle
        endif
        temps%fields(p)%values(1,i,j,k)=recovered
      enddo
      enddo
      enddo
    enddo
    call agree(local_ok,comm,ready)
    if(.not.ready) return
    if(present(rereflux_delta)) then
      ! Publish the repaired pre-register fine state. Strict closure will apply
      ! the original register exactly once; it remains unconsumed on failure.
      local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or.tree%nodes(p)%level<=minimum_level) cycle
        n=tree%nodes(p)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(all(routed_register%fields(p)%values(:,i,j,k)==0)) cycle
          call add_mpi_eb_delta_cell_3d(candidate%fields(p)%values(:,i,j,k), &
            -routed_register%fields(p)%values(:,i,j,k),valid)
          if(.not.valid) then
            local_ok=.false.;cycle
          end if
          if(geometry(p)%volume_fraction(i,j,k)<=0) cycle
          if(any(candidate%fields(p)%values(6:active_components,i,j,k)<0)) then
            local_ok=.false.;cycle
          end if
          call recover_repair_cell(species,geometry(p),[i,j,k],candidate%fields(p)%values(:,i,j,k), &
            temps%fields(p)%values(1,i,j,k),primitive,recovered,speed,valid,frame, &
          repair_node=p)
          local_ok=local_ok.and.valid
          if(valid) temps%fields(p)%values(1,i,j,k)=recovered
        end do
        end do
        end do
      end do
      call agree(local_ok,comm,ready)
      if(.not.ready) return
    end if
    call move_alloc(candidate%fields,state%fields);call move_alloc(temps%fields,temperature%fields)
    flux=working;ok=.true.
  contains
    subroutine recover_repair_cell(table,metrics,index,values,guess,primitive,recovered,speed,valid,rotation,repair_node)
      type(nasa7_species), intent(in) :: table(:)
      type(eb_geometry_3d), intent(in) :: metrics
      integer, intent(in) :: index(3),repair_node
      real(dp), intent(in) :: values(:),guess
      real(dp), intent(out) :: primitive(:),recovered,speed
      logical, intent(out) :: valid
      type(rotating_frame_3d), optional, intent(in) :: rotation
      real(dp) :: gas(active_components),delta(active_components),qp(reactive_nprim(size(table))),tg,sound
      if(passive_components==0) then
        call recover_reactive_eb_cell_3d(table,metrics,index,values,guess,primitive,recovered,speed,valid,rotation)
      else
        call recover_reactive_eb_cell_3d(table,metrics,index,values(:active_components),guess, &
          primitive,recovered,speed,valid,rotation,internal_energy_density=values(active_components+1))
      endif
      if(.not.valid.or..not.present(reactor_energy_bias)) return
      ! Independent UEINT remains the returned temperature authority. Gas and
      ! reactor energy admission are additional queries with private outputs.
      gas=values(:active_components)
      call recover_reactive_eb_cell_3d(table,metrics,index,gas,guess,qp,tg,sound,valid,rotation)
      if(.not.valid) return
      delta=0
      delta(5)=bias_view%fields(repair_node)%values(1,index(1),index(2),index(3))
      call add_mpi_eb_delta_cell_3d(gas,delta,valid)
      if(.not.valid) return
      call recover_reactive_eb_cell_3d(table,metrics,index,gas,guess,qp,tg,sound,valid,rotation)
    end subroutine recover_repair_cell
    subroutine supply_endpoint_packet(target,delta,target_volume,donor,donor_volume, &
        funded_target,funded_donor,packet,complete,valid,passive_components,allow_inventory_packet, &
        donor_node,donor_cell,donor_delta,target_guess,donor_guess,prefer_inventory)
      real(dp), intent(in) :: target(nv),delta(nv),target_volume,donor(nv),donor_volume,donor_delta(nv)
      real(dp), intent(in) :: target_guess,donor_guess
      real(dp), intent(inout) :: funded_target(nv),funded_donor(nv),packet(nv)
      integer, intent(in) :: passive_components,donor_node,donor_cell(3)
      logical, intent(in) :: allow_inventory_packet
      logical, optional, intent(in) :: prefer_inventory
      logical, intent(out) :: complete,valid
      real(dp) :: value(nv),qp(reactive_nprim(size(species))),tg,sound
      logical :: eos_ok,inventory_required
      if(present(reactor_energy_bias)) then
        call supply_reactive_eb_gas_eos_donor_3d(species,geometry(parent),coarse_cell,target_guess, &
          geometry(donor_node),donor_cell,donor_guess,target,delta,target_volume,donor,donor_delta,donor_volume, &
          funded_target,funded_donor,packet,complete,valid,frame,allow_raw_target=all(delta==0), &
          reactor_energy_bias=[bias_view%fields(parent)%values(1,coarse_cell(1),coarse_cell(2),coarse_cell(3)), &
            bias_view%fields(donor_node)%values(1,donor_cell(1),donor_cell(2),donor_cell(3))], &
          native_projection=native_projection,inventory_nodes=[parent,donor_node],check_native_target=native_endpoints)
        return
      endif
      inventory_required=present(native_projection)
      if(present(prefer_inventory)) inventory_required=inventory_required.or.prefer_inventory
      if(inventory_required.and.allow_inventory_packet) then
          ! Similar donor/receiver density makes a contrast packet too small
          ! for a material sink. Search the complete inventory interval with
          ! physical gas and actual held-energy queries before taking a loan.
          ! Zero offset adds no reactor-energy constraint in the held case.
          if(passive_components==0) then
            call supply_reactive_eb_gas_eos_donor_3d(species,geometry(parent),coarse_cell,target_guess, &
              geometry(donor_node),donor_cell,donor_guess,target,delta,target_volume,donor,donor_delta,donor_volume, &
              funded_target,funded_donor,packet,complete,valid,frame,allow_raw_target=all(delta==0), &
              native_projection=native_projection,inventory_nodes=[parent,donor_node],check_native_target=native_endpoints)
          else
            call supply_reactive_eb_gas_eos_donor_3d(species,geometry(parent),coarse_cell,target_guess, &
              geometry(donor_node),donor_cell,donor_guess,target,delta,target_volume,donor,donor_delta,donor_volume, &
              funded_target,funded_donor,packet,complete,valid,frame,allow_raw_target=all(delta==0), &
              reactor_energy_bias=[0._dp,0._dp], &
              native_projection=native_projection,inventory_nodes=[parent,donor_node],check_native_target=native_endpoints)
          endif
          return
      endif
      call supply_reactive_eb_endpoint_donor_3d(target,delta,target_volume,donor,donor_volume, &
        funded_target,funded_donor,packet,complete,valid,passive_components,allow_inventory_packet)
      if(valid.and..not.complete) return
      if(valid) then
        value=funded_target
        call add_mpi_eb_delta_cell_3d(value,delta,eos_ok)
        if(eos_ok) call recover_repair_cell(species,geometry(parent),coarse_cell,value, &
          target_guess,qp,tg,sound,eos_ok,frame, &
          repair_node=parent)
        if(eos_ok) return
      endif
      if(.not.allow_inventory_packet) then
        ! Flux-only gas/species deficits retain the established interface
        ! contrast/covered donor route. New caloric funding requires a
        ! physical density/composition endpoint before thermal admission.
        value=target
        call add_mpi_eb_delta_cell_3d(value,delta,eos_ok)
        valid=.false.;complete=.false.
        if(.not.eos_ok) return
        if(value(1)<=0.or.any(value(6:active_components)<0)) return
      endif
      if(passive_components==0) then
        call supply_reactive_eb_gas_eos_donor_3d(species,geometry(parent),coarse_cell,target_guess, &
          geometry(donor_node),donor_cell,donor_guess,target,delta,target_volume,donor,donor_delta,donor_volume, &
          funded_target,funded_donor,packet,complete,valid,frame,allow_raw_target=all(delta==0), &
              native_projection=native_projection,inventory_nodes=[parent,donor_node],check_native_target=native_endpoints)
      else
        call supply_reactive_eb_held_eos_donor_3d(species,geometry(parent),coarse_cell,target_guess, &
          geometry(donor_node),donor_cell,donor_guess,target,delta,target_volume,donor,donor_delta,donor_volume, &
          funded_target,funded_donor,packet,complete,valid,frame,allow_raw_target=all(delta==0))
      endif
    end subroutine supply_endpoint_packet
    subroutine fund_exposed_fine_endpoint(funded)
      logical, intent(out) :: funded
      real(dp) :: source(nv),delta(nv),donor(nv),donor_delta(nv),target_next(nv),target_post(nv)
      real(dp) :: donor_next(nv),donor_pre(nv),packet(nv),face_packet(nv),mirror_packet(nv)
      real(dp) :: source_temp,donor_temp,recovered_target,recovered_donor,sound,area,mirror_area,volume
      integer :: a,step,distance,h,limit,node,next_node,cell(3),next_cell(3),face(3),mirror(3),terminal,terminal_cell(3)
      logical :: connected,path_ok,packet_ok,complete,good,post_good,faces_ok,covered_terminal
      integer :: broadcast_status
      logical :: all_faces_ok
      funded=.false.
      search: do a=1,3
      do step=-1,1,2
        limit=int(min(domain(a,tree%nodes(parent)%level+1)-1,int(huge(distance)-1,int64)))
        do distance=1,limit
          node=parent;cell=coarse_cell;path_ok=.true.
          do h=1,distance
            call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
            if(.not.connected) then
              path_ok=.false.;exit
            endif
            node=next_node;cell=next_cell
          enddo
          if(.not.path_ok) exit
          terminal=node;terminal_cell=cell
          if(rank==coarse_owner) then
          source=candidate%fields(parent)%values(:,i,j,k)
          delta=0
          if(present(rereflux_delta)) delta=routed_register%fields(parent)%values(:,i,j,k)
          source_temp=temps%fields(parent)%values(1,i,j,k)
          endif
          call MPI_Bcast(source,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
          if(broadcast_status/=MPI_SUCCESS) return
          call MPI_Bcast(delta,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
          if(broadcast_status/=MPI_SUCCESS) return
          call MPI_Bcast(source_temp,1,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
          if(broadcast_status/=MPI_SUCCESS) return
          call add_mpi_eb_delta_cell_3d(source,-delta,good)
          if(.not.good) return
          if(rank==owners(terminal)) then
          donor=candidate%fields(terminal)%values(:,cell(1),cell(2),cell(3))
          donor_delta=0
          if(present(rereflux_delta)) donor_delta=routed_register%fields(terminal)%values(:,cell(1),cell(2),cell(3))
          donor_temp=temps%fields(terminal)%values(1,cell(1),cell(2),cell(3))
          endif
          call MPI_Bcast(donor,nv,MPI_DOUBLE_PRECISION,owners(terminal),comm,broadcast_status)
          if(broadcast_status/=MPI_SUCCESS) return
          call MPI_Bcast(donor_delta,nv,MPI_DOUBLE_PRECISION,owners(terminal),comm,broadcast_status)
          if(broadcast_status/=MPI_SUCCESS) return
          call MPI_Bcast(donor_temp,1,MPI_DOUBLE_PRECISION,owners(terminal),comm,broadcast_status)
          if(broadcast_status/=MPI_SUCCESS) return
          covered_terminal=covering_child(terminal,cell)/=0
          if(covered_terminal) then
            donor_pre=donor
            call covered_ancestor_inventory(terminal,cell,donor,good)
            if(.not.good) cycle
            donor_delta=donor
            call add_mpi_eb_delta_cell_3d(donor_delta,-donor_pre,good)
            if(.not.good) cycle
          endif
          volume=cell_volume(terminal,cell)
          call supply_endpoint_packet(source,delta,coarse_volume,donor,volume, &
            target_next,donor_next,packet,complete,packet_ok,passive_components,.true., &
            terminal,cell,donor_delta,source_temp,donor_temp)
          if(.not.packet_ok) cycle
          call recover_repair_cell(species,geometry(parent),coarse_cell,target_next, &
            source_temp,primitive,recovered_target,sound,good,frame, &
          repair_node=parent)
          if(.not.good.and.(complete.or.any(delta/=0))) cycle
          target_post=target_next
          call add_mpi_eb_delta_cell_3d(target_post,delta,good)
          if(.not.good) cycle
          call recover_repair_cell(species,geometry(parent),coarse_cell,target_post, &
            source_temp,primitive,recovered_target,sound,post_good,frame, &
          repair_node=parent)
          if(complete.and..not.post_good) cycle
          call restore_coarse_donor(terminal,cell,donor_next,donor_delta,donor_temp, &
            donor_pre,recovered_donor,good)
          if(.not.good) cycle
          ! Every route face carries the same extensive packet. Exposed
          ! transit cells receive zero net divergence, including peer seams.
          node=parent;cell=coarse_cell;faces_ok=.true.
          do h=1,distance
            call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
            area=open_area(node,a,face);mirror_area=open_area(next_node,a,mirror)
            if(abs(area-mirror_area)>64*epsilon(1._dp)*max(area,mirror_area)) then
              faces_ok=.false.;exit
            endif
            if(area<1) then
              if(any(abs(packet)>huge(1._dp)*area)) then
                faces_ok=.false.;exit
              endif
            endif
            if(mirror_area<1) then
              if(any(abs(packet)>huge(1._dp)*mirror_area)) then
                faces_ok=.false.;exit
              endif
            endif
            face_packet=-real(step,dp)*packet/area
            mirror_packet=-real(step,dp)*packet/mirror_area
            if(any(.not.ieee_is_finite(face_packet)).or.any(.not.ieee_is_finite(mirror_packet))) then
              faces_ok=.false.;exit
            endif
            if(rank==owners(node)) then
            if(any(abs(face_value(working(node),a,face))>huge(1._dp)-abs(face_packet)).or. &
              any(abs(face_value(effective(node),a,face))>huge(1._dp)-abs(face_packet))) faces_ok=.false.
            endif
            if(node/=next_node.or.any(face/=mirror)) then
              if(rank==owners(next_node)) then
              if(any(abs(face_value(working(next_node),a,mirror))>huge(1._dp)-abs(mirror_packet)).or. &
                any(abs(face_value(effective(next_node),a,mirror))>huge(1._dp)-abs(mirror_packet))) faces_ok=.false.
              endif
            endif
            node=next_node;cell=next_cell
          enddo
          call agree(faces_ok,comm,all_faces_ok)
          faces_ok=all_faces_ok
          if(.not.faces_ok) cycle
          ! Allocate covered transit/terminal packets to physical descendants.
          ! This checks every changed leaf's pre/post registered EOS and
          ! rolls back a refused route before trying another physical donor.
          call propagate_coarse_route(packet,a,step,distance,good)
          if(.not.good) cycle
          if(covered_terminal) donor_next=donor_pre
          if(rank==coarse_owner) then
          candidate%fields(parent)%values(:,i,j,k)=target_post
          if(post_good) temps%fields(parent)%values(1,i,j,k)=recovered_target
          changed(parent)=.true.
          endif
          if(rank==owners(terminal)) then
          candidate%fields(terminal)%values(:,terminal_cell(1),terminal_cell(2),terminal_cell(3))=donor_next
          changed(terminal)=.true.
          endif
          node=parent;cell=coarse_cell
          do h=1,distance
            call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
            if(rank==owners(node)) then
            call add_path_packet(node,a,face,-real(step,dp)*packet)
            endif
            if(node/=next_node.or.any(face/=mirror)) then
              if(rank==owners(next_node)) then
              call add_path_packet(next_node,a,mirror,-real(step,dp)*packet)
              endif
            endif
            node=next_node;cell=next_cell
          enddo
          if(complete) then
            funded=.true.;exit search
          endif
        enddo
      enddo
      enddo search
      if(.not.funded.and.present(rereflux_delta)) call fund_adjacent_coarse_endpoint(funded)
    end subroutine fund_exposed_fine_endpoint
    ! A fine leaf can borrow from the exposed physical coarse cell across
    ! its own CF boundary after same-level stock is exhausted. Both faces
    ! carry one extensive gas/held packet; different face areas do not imply
    ! different inventories. Virtual covered coarse cells remain excluded.
    subroutine fund_adjacent_coarse_endpoint(funded)
      logical, intent(out) :: funded
      real(dp) :: source(nv),delta(nv),donor(nv),donor_delta(nv),target_next(nv),target_post(nv)
      real(dp) :: donor_next(nv),donor_pre(nv),packet(nv),face_packet(nv),mirror_packet(nv)
      real(dp) :: source_temp,donor_temp,target_temp,donor_temp_next,sound,area,mirror_area,volume
      integer :: a,step,peer,donor_node,depth,r,cell(3),face(3),mirror(3)
      integer(int64) :: neighbor_global(3),coarse_global(3)
      logical :: packet_ok,complete,good,post_good,faces_ok
      integer :: broadcast_status
      logical :: all_faces_ok
      funded=.false.;depth=tree%nodes(parent)%level
      if(depth<1) return
      r=tree%nodes(parent)%patch%refinement_ratio
      if(r<2) return
      do a=1,3
      do step=-1,1,2
        if(step<0.and.coarse_cell(a)/=1) cycle
        if(step>0.and.coarse_cell(a)/=tree%nodes(parent)%shape(a)) cycle
        neighbor_global=origin(:,parent)+int(coarse_cell-1,int64)
        neighbor_global(a)=neighbor_global(a)+step
        if(periodic(a)) neighbor_global(a)=modulo(neighbor_global(a),domain(a,depth+1))
        if(any(neighbor_global<0).or.any(neighbor_global>=domain(:,depth+1))) cycle
        coarse_global=neighbor_global/int(r,int64);donor_node=0
        do peer=1,size(tree%nodes)
          if(tree%nodes(peer)%level/=depth-1) cycle
          if(any(coarse_global<origin(:,peer)).or. &
            any(coarse_global>=origin(:,peer)+int(tree%nodes(peer)%shape,int64))) cycle
          donor_node=peer;cell=int(coarse_global-origin(:,peer))+1;exit
        enddo
        if(donor_node==0) cycle
        if(covering_child(donor_node,cell)/=0) cycle
        volume=cell_volume(donor_node,cell)
        if(volume<=0) cycle
        face=coarse_cell;face(a)=face(a)+min(step,0)
        mirror=cell;mirror(a)=mirror(a)+min(-step,0)
        area=open_area(parent,a,face);mirror_area=open_area(donor_node,a,mirror)
        if(area<=0.or.mirror_area<=0) cycle
        if(area>mirror_area) then
          if(area-mirror_area>64*epsilon(1._dp)*area) cycle
        endif
        if(rank==owners(parent)) then
        source=candidate%fields(parent)%values(:,coarse_cell(1),coarse_cell(2),coarse_cell(3));delta=0
        if(present(rereflux_delta)) delta=routed_register%fields(parent)%values(:,coarse_cell(1),coarse_cell(2),coarse_cell(3))
        source_temp=temps%fields(parent)%values(1,coarse_cell(1),coarse_cell(2),coarse_cell(3))
        endif
        call MPI_Bcast(source,nv,MPI_DOUBLE_PRECISION,owners(parent),comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call MPI_Bcast(delta,nv,MPI_DOUBLE_PRECISION,owners(parent),comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call MPI_Bcast(source_temp,1,MPI_DOUBLE_PRECISION,owners(parent),comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call add_mpi_eb_delta_cell_3d(source,-delta,good)
        if(.not.good) return
        if(rank==owners(donor_node)) then
        donor=candidate%fields(donor_node)%values(:,cell(1),cell(2),cell(3))
        call cell_increment(donor_node,cell,[0,0,0],0,donor_delta)
        donor_temp=temps%fields(donor_node)%values(1,cell(1),cell(2),cell(3))
        endif
        call MPI_Bcast(donor,nv,MPI_DOUBLE_PRECISION,owners(donor_node),comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call MPI_Bcast(donor_delta,nv,MPI_DOUBLE_PRECISION,owners(donor_node),comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call MPI_Bcast(donor_temp,1,MPI_DOUBLE_PRECISION,owners(donor_node),comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        ! Donor kernels receive its complete post-debt inventory. They
        ! also retain the pre-debt endpoint by subtracting this readonly debt.
        call add_mpi_eb_delta_cell_3d(donor,donor_delta,good)
        if(.not.good) cycle
        if(donor(1)<=0.or.any(donor(6:active_components)<0)) cycle
        call supply_endpoint_packet(source,delta,coarse_volume,donor,volume,target_next,donor_next,packet, &
          complete,packet_ok,passive_components,.true.,donor_node,cell,donor_delta,source_temp,donor_temp, &
          prefer_inventory=.true.)
        if(.not.packet_ok) cycle
        call recover_repair_cell(species,geometry(parent),coarse_cell,target_next,source_temp, &
          primitive,target_temp,sound,good,frame,repair_node=parent)
        if(.not.good.and.(complete.or.any(delta/=0))) cycle
        target_post=target_next
        call add_mpi_eb_delta_cell_3d(target_post,delta,good)
        if(.not.good) cycle
        call recover_repair_cell(species,geometry(parent),coarse_cell,target_post,source_temp, &
          primitive,target_temp,sound,post_good,frame,repair_node=parent)
        if(complete.and..not.post_good) cycle
        call restore_coarse_donor(donor_node,cell,donor_next,donor_delta,donor_temp,donor_pre,donor_temp_next,good)
        if(.not.good) cycle
        if(area<1) then
          if(any(abs(packet)>huge(1._dp)*area)) cycle
        endif
        if(mirror_area<1) then
          if(any(abs(packet)>huge(1._dp)*mirror_area)) cycle
        endif
        face_packet=-real(step,dp)*packet/area;mirror_packet=-real(step,dp)*packet/mirror_area
        faces_ok=all(ieee_is_finite(face_packet)).and.all(ieee_is_finite(mirror_packet))
        if(rank==owners(parent)) then
        faces_ok=faces_ok.and.all(abs(face_value(working(parent),a,face))<=huge(1._dp)-abs(face_packet)).and. &
          all(abs(face_value(effective(parent),a,face))<=huge(1._dp)-abs(face_packet))
        endif
        if(rank==owners(donor_node)) then
        faces_ok=faces_ok.and.all(abs(face_value(working(donor_node),a,mirror))<=huge(1._dp)-abs(mirror_packet)).and. &
          all(abs(face_value(effective(donor_node),a,mirror))<=huge(1._dp)-abs(mirror_packet))
        endif
        call agree(faces_ok,comm,all_faces_ok)
        faces_ok=all_faces_ok
        if(.not.faces_ok) cycle
        if(rank==owners(parent)) then
        candidate%fields(parent)%values(:,coarse_cell(1),coarse_cell(2),coarse_cell(3))=target_post
        if(post_good) temps%fields(parent)%values(1,coarse_cell(1),coarse_cell(2),coarse_cell(3))=target_temp
        changed(parent)=.true.
        call add_path_packet(parent,a,face,-real(step,dp)*packet)
        endif
        if(rank==owners(donor_node)) then
        candidate%fields(donor_node)%values(:,cell(1),cell(2),cell(3))=donor_pre
        temps%fields(donor_node)%values(1,cell(1),cell(2),cell(3))=donor_temp_next
        changed(donor_node)=.true.
        call add_path_packet(donor_node,a,mirror,-real(step,dp)*packet)
        endif
        if(complete) then
          funded=.true.;return
        endif
      enddo
      enddo
    end subroutine fund_adjacent_coarse_endpoint

    subroutine fund_exposed_coarse_endpoint(funded)
      logical, intent(out) :: funded
      real(dp) :: correction(nv),donor_correction(nv),donor_budget(nv),target_next(nv),donor_next(nv),donor_pre(nv),packet(nv)
      real(dp) :: replacement(nv),replacement_effective(nv),other_volume,area,sound,face_packet(nv)
      real(dp) :: source(nv),donor(nv),source_temp,donor_temp,mirror_area,mirror_packet(nv)
      real(dp) :: mirror_replacement(nv),mirror_effective(nv)
      real(dp) :: recovered_target,recovered_donor
      integer :: a,step,neighbor(3),face(3),mirror(3),peer,peer_owner,broadcast_status
      logical :: complete,packet_complete,valid_target,valid_donor,valid_packet,connected,face_ok,all_face_ok
      logical :: distant_funded
      logical :: use_inventory
      funded=.false.;complete=.false.;use_inventory=.false.
      if(rank==coarse_owner) then
        if(present(rereflux_delta)) use_inventory=any(routed_register%fields(parent)%values(:,i,j,k)/=0)
      endif
      call MPI_Bcast(use_inventory,1,MPI_LOGICAL,coarse_owner,comm,broadcast_status)
      if(broadcast_status/=MPI_SUCCESS) return
      if(rank==coarse_owner.or.immediate_funding) then
        if(rank==coarse_owner) call cell_increment(parent,coarse_cell,[0,0,0],0,correction)
        if(immediate_funding) then
          call MPI_Bcast(correction,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
          if(broadcast_status/=MPI_SUCCESS) return
        endif
        axes: do a=1,3
        do step=-1,1,2
          neighbor=coarse_cell;neighbor(a)=neighbor(a)+step
          if(any(neighbor<1).or.any(neighbor>n)) cycle
          if(covering_child(parent,neighbor)/=0) cycle
          if(geometry(parent)%volume_fraction(neighbor(1),neighbor(2),neighbor(3))<=0) cycle
          face=coarse_cell;face(a)=face(a)+min(step,0)
          area=open_area(parent,a,face)
          if(area<=0) cycle
          if(rank==coarse_owner) then
            source=candidate%fields(parent)%values(:,i,j,k)
            donor=candidate%fields(parent)%values(:,neighbor(1),neighbor(2),neighbor(3))
            source_temp=temps%fields(parent)%values(1,i,j,k)
            donor_temp=temps%fields(parent)%values(1,neighbor(1),neighbor(2),neighbor(3))
            call cell_increment(parent,neighbor,[0,0,0],0,donor_correction)
          endif
          if(immediate_funding) then
            ! The native packet query contains collective reactor admission.
            ! Replicate this read-only proposal before every rank enters it;
            ! candidate fields and face replacements remain owner-local.
            call MPI_Bcast(source,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
            if(broadcast_status/=MPI_SUCCESS) return
            call MPI_Bcast(donor,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
            if(broadcast_status/=MPI_SUCCESS) return
            call MPI_Bcast(donor_correction,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
            if(broadcast_status/=MPI_SUCCESS) return
            call MPI_Bcast(source_temp,1,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
            if(broadcast_status/=MPI_SUCCESS) return
            call MPI_Bcast(donor_temp,1,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
            if(broadcast_status/=MPI_SUCCESS) return
          endif
          call pending_coarse_donor(donor,donor_correction,donor_budget,valid_donor)
          if(.not.valid_donor) cycle
          other_volume=cell_volume(parent,neighbor)
          target_next=0;donor_next=0;packet=0
          call supply_endpoint_packet(source,correction, &
            coarse_volume,donor_budget,other_volume, &
            target_next,donor_next,packet,packet_complete,valid_packet, &
        passive_components=passive_components,allow_inventory_packet=use_inventory, &
          donor_node=parent,donor_cell=neighbor,donor_delta=donor_correction, &
          target_guess=source_temp,donor_guess=donor_temp)
          if(.not.valid_packet) cycle
          call recover_repair_cell(species,geometry(parent),coarse_cell,target_next, &
            source_temp,primitive,recovered_target,sound,valid_target,frame, &
          repair_node=parent)
          if(.not.valid_target.and.(packet_complete.or.any(correction/=0))) cycle
          call restore_coarse_donor(parent,neighbor,donor_next,donor_correction, &
            donor_temp,donor_pre,recovered_donor,valid_donor)
          if(.not.valid_donor) cycle
          if(area<1) then
            if(any(abs(packet)>huge(1._dp)*area)) cycle
          endif
          face_packet=-real(step,dp)*packet/area
          face_ok=all(ieee_is_finite(face_packet))
          if(rank==coarse_owner) then
            face_ok=face_ok.and.all(abs(face_value(working(parent),a,face))<=huge(1._dp)-abs(face_packet)).and. &
              all(abs(face_value(effective(parent),a,face))<=huge(1._dp)-abs(face_packet))
            if(face_ok) then
              replacement=face_value(working(parent),a,face)+face_packet
              replacement_effective=face_value(effective(parent),a,face)+face_packet
              face_ok=all(ieee_is_finite(replacement)).and.all(ieee_is_finite(replacement_effective))
            endif
          endif
          if(immediate_funding) then
            call MPI_Bcast(face_ok,1,MPI_LOGICAL,coarse_owner,comm,broadcast_status)
            if(broadcast_status/=MPI_SUCCESS) return
          endif
          if(.not.face_ok) cycle
          if(rank==coarse_owner) then
            candidate%fields(parent)%values(:,i,j,k)=target_next
            candidate%fields(parent)%values(:,neighbor(1),neighbor(2),neighbor(3))=donor_pre
            if(valid_target) temps%fields(parent)%values(1,i,j,k)=recovered_target
            temps%fields(parent)%values(1,neighbor(1),neighbor(2),neighbor(3))=recovered_donor
            call set_face(working(parent),a,face,replacement)
            call set_face(effective(parent),a,face,replacement_effective)
            changed(parent)=.true.
          endif
          funded=.true.;complete=packet_complete
          if(complete) exit axes
        enddo
        enddo axes
      endif
      call MPI_Bcast(funded,1,MPI_LOGICAL,coarse_owner,comm,broadcast_status)
      if(broadcast_status/=MPI_SUCCESS) funded=.false.
      call MPI_Bcast(complete,1,MPI_LOGICAL,coarse_owner,comm,broadcast_status)
      if(broadcast_status/=MPI_SUCCESS.or.complete) return
      ! Candidate inventories live only on their owners. Broadcast both
      ! endpoints before evaluating a peer-patch packet on every rank.
      peers: do a=1,3
      do step=-1,1,2
        call path_neighbor(parent,coarse_cell,a,step,peer,neighbor,face,mirror,connected)
        if(.not.connected.or.peer==parent) cycle
        if(covering_child(peer,neighbor)/=0) cycle
        area=open_area(parent,a,face);mirror_area=open_area(peer,a,mirror)
        if(abs(area-mirror_area)>64*epsilon(1._dp)*max(area,mirror_area)) cycle
        peer_owner=owners(peer)
        if(rank==coarse_owner) then
          source=candidate%fields(parent)%values(:,i,j,k)
          source_temp=temps%fields(parent)%values(1,i,j,k)
          call cell_increment(parent,coarse_cell,[0,0,0],0,correction)
        endif
        call MPI_Bcast(source,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call MPI_Bcast(source_temp,1,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call MPI_Bcast(correction,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        if(rank==peer_owner) then
          donor=candidate%fields(peer)%values(:,neighbor(1),neighbor(2),neighbor(3))
          donor_temp=temps%fields(peer)%values(1,neighbor(1),neighbor(2),neighbor(3))
          call cell_increment(peer,neighbor,[0,0,0],0,donor_correction)
        endif
        call MPI_Bcast(donor,nv,MPI_DOUBLE_PRECISION,peer_owner,comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call MPI_Bcast(donor_temp,1,MPI_DOUBLE_PRECISION,peer_owner,comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call MPI_Bcast(donor_correction,nv,MPI_DOUBLE_PRECISION,peer_owner,comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call pending_coarse_donor(donor,donor_correction,donor_budget,valid_donor)
        if(.not.valid_donor) cycle
        other_volume=cell_volume(peer,neighbor)
        call supply_endpoint_packet(source,correction,coarse_volume,donor_budget,other_volume, &
          target_next,donor_next,packet,packet_complete,valid_packet, &
        passive_components=passive_components,allow_inventory_packet=use_inventory, &
          donor_node=peer,donor_cell=neighbor,donor_delta=donor_correction, &
          target_guess=source_temp,donor_guess=donor_temp)
        if(.not.valid_packet) cycle
        call recover_repair_cell(species,geometry(parent),coarse_cell,target_next, &
          source_temp,primitive,recovered_target,sound,valid_target,frame, &
          repair_node=parent)
        if(.not.valid_target.and.(packet_complete.or.any(correction/=0))) cycle
        call restore_coarse_donor(peer,neighbor,donor_next,donor_correction, &
          donor_temp,donor_pre,recovered_donor,valid_donor)
        if(.not.valid_donor) cycle
        if(area<1) then
          if(any(abs(packet)>huge(1._dp)*area)) cycle
        endif
        face_packet=-real(step,dp)*packet/area
        mirror_packet=-real(step,dp)*packet/mirror_area
        if(any(.not.ieee_is_finite(face_packet)).or.any(.not.ieee_is_finite(mirror_packet))) cycle
        face_ok=.true.
        if(rank==coarse_owner) then
          face_ok=all(abs(face_value(working(parent),a,face))<=huge(1._dp)-abs(face_packet)).and. &
            all(abs(face_value(effective(parent),a,face))<=huge(1._dp)-abs(face_packet))
          if(face_ok) then
            replacement=face_value(working(parent),a,face)+face_packet
            replacement_effective=face_value(effective(parent),a,face)+face_packet
            face_ok=all(ieee_is_finite(replacement)).and.all(ieee_is_finite(replacement_effective))
          endif
        endif
        if(rank==peer_owner) then
          face_ok=face_ok.and.all(abs(face_value(working(peer),a,mirror))<=huge(1._dp)-abs(mirror_packet)).and. &
            all(abs(face_value(effective(peer),a,mirror))<=huge(1._dp)-abs(mirror_packet))
          if(face_ok) then
            mirror_replacement=face_value(working(peer),a,mirror)+mirror_packet
            mirror_effective=face_value(effective(peer),a,mirror)+mirror_packet
            face_ok=all(ieee_is_finite(mirror_replacement)).and.all(ieee_is_finite(mirror_effective))
          endif
        endif
        call agree(face_ok,comm,all_face_ok)
        if(.not.all_face_ok) cycle
        if(rank==coarse_owner) then
          candidate%fields(parent)%values(:,i,j,k)=target_next
          if(valid_target) temps%fields(parent)%values(1,i,j,k)=recovered_target
          call set_face(working(parent),a,face,replacement)
          call set_face(effective(parent),a,face,replacement_effective)
          changed(parent)=.true.
        endif
        if(rank==peer_owner) then
          candidate%fields(peer)%values(:,neighbor(1),neighbor(2),neighbor(3))=donor_pre
          temps%fields(peer)%values(1,neighbor(1),neighbor(2),neighbor(3))=recovered_donor
          call set_face(working(peer),a,mirror,mirror_replacement)
          call set_face(effective(peer),a,mirror,mirror_effective)
          changed(peer)=.true.
        endif
        funded=.true.;complete=packet_complete
        if(complete) exit peers
      enddo
      enddo peers
      if(.not.complete) then
        call fund_distant_coarse_endpoint(distant_funded,use_inventory)
        funded=funded.or.distant_funded
      endif
    end subroutine fund_exposed_coarse_endpoint

    subroutine fund_distant_coarse_endpoint(funded,use_inventory)
      logical, intent(out) :: funded
      logical, intent(in) :: use_inventory
      real(dp) :: correction(nv),donor_correction(nv),donor_budget(nv),source(nv),donor(nv)
      real(dp) :: target_next(nv),donor_next(nv),donor_pre(nv),packet(nv)
      real(dp) :: source_temp,donor_temp,area,mirror_area,volume,face_packet(nv),mirror_packet(nv)
      real(dp) :: sound,recovered_target,recovered_donor
      integer :: a,step,distance,h,limit,node,next_node,cell(3),next_cell(3),face(3),mirror(3)
      integer :: terminal_node,terminal_cell(3),terminal_owner,broadcast_status
      logical :: connected,valid_path,all_path_ok,valid_packet,complete,valid_target,valid_donor,transit_valid
      funded=.false.
      search: do a=1,3
      do step=-1,1,2
      limit=int(min(domain(a,tree%nodes(parent)%level+1)-1,int(huge(distance)-1,int64)))
      do distance=1,limit
        node=parent;cell=coarse_cell;valid_path=.true.
        do h=1,distance
          call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
          if(.not.connected) then
            valid_path=.false.;exit
          endif
          node=next_node;cell=next_cell
        enddo
        call agree(valid_path,comm,all_path_ok)
        if(.not.all_path_ok) exit
        terminal_node=node;terminal_cell=cell;terminal_owner=owners(node)
        if(distance==1.and.covering_child(terminal_node,terminal_cell)==0) cycle
        ! Adjacent-interface covered cells retain the existing fine-level
        ! route. A farther covered donor is reached through this coarse path;
        ! the interface search does not visit an exposed recipient away from it.
        if(tree%nodes(parent)%level==minimum_level.and.terminal_node==parent.and. &
          covering_child(terminal_node,terminal_cell)/=0) then
          ! Native gas/H packet admission is opt-in; flux-only reconciliation
          ! retains its existing finest-level inventory rejection contract.
          if(.not.immediate_funding) cycle
        endif
        if(rank==coarse_owner) then
          source=candidate%fields(parent)%values(:,i,j,k)
          source_temp=temps%fields(parent)%values(1,i,j,k)
          call cell_increment(parent,coarse_cell,[0,0,0],0,correction)
        endif
        call MPI_Bcast(source,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call MPI_Bcast(source_temp,1,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call MPI_Bcast(correction,nv,MPI_DOUBLE_PRECISION,coarse_owner,comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        if(rank==terminal_owner) then
          donor=candidate%fields(node)%values(:,cell(1),cell(2),cell(3))
          donor_temp=temps%fields(node)%values(1,cell(1),cell(2),cell(3))
          call cell_increment(node,cell,[0,0,0],0,donor_correction)
        endif
        call MPI_Bcast(donor,nv,MPI_DOUBLE_PRECISION,terminal_owner,comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call MPI_Bcast(donor_temp,1,MPI_DOUBLE_PRECISION,terminal_owner,comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        call MPI_Bcast(donor_correction,nv,MPI_DOUBLE_PRECISION,terminal_owner,comm,broadcast_status)
        if(broadcast_status/=MPI_SUCCESS) return
        if((tree%nodes(parent)%level<minimum_level.or.present(native_projection)).and. &
          covering_child(terminal_node,terminal_cell)/=0) then
          ! The collective leaf traversal counts pending debt once; owner
          ! broadcasts above supply only the virtual publication anchor.
          call covered_ancestor_inventory(terminal_node,terminal_cell,donor_budget,valid_donor)
          if(valid_donor) then
            donor_correction=donor_budget
            call add_mpi_eb_delta_cell_3d(donor_correction,-donor,valid_donor)
          endif
        else
          call pending_coarse_donor(donor,donor_correction,donor_budget,valid_donor)
        endif
        if(.not.valid_donor) cycle
        volume=cell_volume(terminal_node,terminal_cell)
        call supply_endpoint_packet(source,correction,coarse_volume,donor_budget,volume, &
          target_next,donor_next,packet,complete,valid_packet, &
        passive_components=passive_components,allow_inventory_packet=use_inventory, &
          donor_node=terminal_node,donor_cell=terminal_cell,donor_delta=donor_correction, &
          target_guess=source_temp,donor_guess=donor_temp)
        if(.not.valid_packet) cycle
        call recover_repair_cell(species,geometry(parent),coarse_cell,target_next, &
          source_temp,primitive,recovered_target,sound,valid_target,frame, &
          repair_node=parent)
        if(.not.valid_target.and.(complete.or.any(correction/=0))) cycle
        call restore_coarse_donor(terminal_node,terminal_cell,donor_next,donor_correction, &
          donor_temp,donor_pre,recovered_donor,valid_donor)
        if(.not.valid_donor) cycle
        node=parent;cell=coarse_cell;valid_path=.true.
        do h=1,distance
          call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
          area=open_area(node,a,face);mirror_area=open_area(next_node,a,mirror)
          if(abs(area-mirror_area)>64*epsilon(1._dp)*max(area,mirror_area)) valid_path=.false.
          if(area<1) then
            if(any(abs(packet)>huge(1._dp)*area)) valid_path=.false.
          endif
          if(.not.valid_path) exit
          face_packet=-real(step,dp)*packet/area
          mirror_packet=-real(step,dp)*packet/mirror_area
          if(any(.not.ieee_is_finite(face_packet)).or.any(.not.ieee_is_finite(mirror_packet))) then
            valid_path=.false.;exit
          endif
          if(rank==owners(node)) then
            if(any(.not.ieee_is_finite(face_value(working(node),a,face))).or. &
              any(.not.ieee_is_finite(face_value(effective(node),a,face)))) valid_path=.false.
            if(any(abs(face_value(working(node),a,face))>huge(1._dp)-abs(face_packet)).or. &
              any(abs(face_value(effective(node),a,face))>huge(1._dp)-abs(face_packet))) valid_path=.false.
          endif
          if(node/=next_node.or.any(face/=mirror)) then
            if(rank==owners(next_node)) then
              if(any(.not.ieee_is_finite(face_value(working(next_node),a,mirror))).or. &
                any(.not.ieee_is_finite(face_value(effective(next_node),a,mirror)))) valid_path=.false.
              if(any(abs(face_value(working(next_node),a,mirror))>huge(1._dp)-abs(mirror_packet)).or. &
                any(abs(face_value(effective(next_node),a,mirror))>huge(1._dp)-abs(mirror_packet))) valid_path=.false.
            endif
          endif
          if(.not.valid_path) exit
          node=next_node;cell=next_cell
        enddo
        call agree(valid_path,comm,all_path_ok)
        if(.not.all_path_ok) cycle
        call propagate_coarse_route(packet,a,step,distance,transit_valid)
        if(.not.transit_valid) cycle
        if(rank==coarse_owner) then
          candidate%fields(parent)%values(:,i,j,k)=target_next
          if(valid_target) temps%fields(parent)%values(1,i,j,k)=recovered_target
        endif
        if(rank==terminal_owner) then
          candidate%fields(terminal_node)%values(:,cell(1),cell(2),cell(3))=donor_pre
          temps%fields(terminal_node)%values(1,cell(1),cell(2),cell(3))=recovered_donor
        endif
        node=parent;cell=coarse_cell
        do h=1,distance
          call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
          if(rank==owners(node)) then
            call add_path_packet(node,a,face,-real(step,dp)*packet)
            changed(node)=.true.
          endif
          if(node/=next_node.or.any(face/=mirror)) then
            if(rank==owners(next_node)) then
              call add_path_packet(next_node,a,mirror,-real(step,dp)*packet)
              changed(next_node)=.true.
            endif
          endif
          node=next_node;cell=next_cell
        enddo
        funded=.true.
        if(complete) exit search
      enddo
      enddo
      enddo search
    end subroutine fund_distant_coarse_endpoint
    ! Covered ancestor averages are not spendable inventories. Integrate
    ! physical leaves after their pending register/reflux debt, once; fine
    ! leaves already carry that debt privately in this transaction.
    recursive subroutine covered_ancestor_inventory(node,cell,budget,valid)
      integer, intent(in) :: node,cell(3)
      real(dp), intent(out) :: budget(nv)
      logical, intent(out) :: valid
      real(dp) :: child_budget(nv),term(nv),delta(nv),volume,other_volume,total_volume,weight
      real(dp) :: qp(reactive_nprim(size(species))),tg,sound
      integer :: covered_by,lower(3),upper(3),first(3),last(3),r,ii,jj,kk
      integer :: bcast_status
      budget=0;valid=.false.;volume=cell_volume(node,cell)
      if(volume<=0.or..not.ieee_is_finite(volume)) return
      covered_by=covering_child(node,cell)
      if(covered_by==0) then
        if(rank==owners(node)) then
          budget=candidate%fields(node)%values(:,cell(1),cell(2),cell(3))
          if(tree%nodes(node)%level<=minimum_level) then
            call cell_increment(node,cell,[0,0,0],0,delta)
            call add_mpi_eb_delta_cell_3d(budget,delta,valid)
          else
            valid=.true.
          endif
          if(valid) valid=all(ieee_is_finite(budget)).and.budget(1)>0.and.all(budget(6:active_components)>=0)
          if(valid) call recover_repair_cell(species,geometry(node),cell,budget, &
            temps%fields(node)%values(1,cell(1),cell(2),cell(3)),qp,tg,sound,valid,frame, &
          repair_node=node)
        endif
        call MPI_Bcast(valid,1,MPI_LOGICAL,owners(node),comm,bcast_status)
        if(bcast_status/=MPI_SUCCESS.or..not.valid) then
          valid=.false.;return
        endif
        call MPI_Bcast(budget,nv,MPI_DOUBLE_PRECISION,owners(node),comm,bcast_status)
        if(bcast_status/=MPI_SUCCESS) valid=.false.
        return
      endif
      call bounds(tree,covered_by,lower,upper)
      r=tree%nodes(covered_by)%patch%refinement_ratio
      first=(cell-lower)*r+1;last=first+r-1;total_volume=0
      do kk=first(3),last(3)
      do jj=first(2),last(2)
      do ii=first(1),last(1)
        other_volume=cell_volume(covered_by,[ii,jj,kk])
        if(other_volume<=0) cycle
        if(.not.ieee_is_finite(other_volume).or.total_volume>huge(1._dp)-other_volume) then
          valid=.false.;return
        endif
        call covered_ancestor_inventory(covered_by,[ii,jj,kk],child_budget,valid)
        if(.not.valid) return
        total_volume=total_volume+other_volume;weight=other_volume/volume
        if(weight>1) then
          if(any(abs(child_budget)>huge(1._dp)/weight)) then
            valid=.false.;return
          endif
        endif
        term=weight*child_budget
        call add_mpi_eb_delta_cell_3d(budget,term,valid)
        if(.not.valid) return
      enddo
      enddo
      enddo
      valid=abs(total_volume-volume)<=64*epsilon(1._dp)*volume
    end subroutine covered_ancestor_inventory
    subroutine qualify_ancestor_route(baseline,valid)
      type(amr_tree_state_3d), intent(in) :: baseline
      logical, intent(out) :: valid
      real(dp) :: value(nv),delta(nv),qp(reactive_nprim(size(species))),tg,sound
      integer :: node,d(3),ii,jj,kk
      logical :: good,local_ok
      local_ok=.true.
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(.not.changed(node)) cycle
        d=tree%nodes(node)%shape
        do kk=1,d(3)
        do jj=1,d(2)
        do ii=1,d(1)
          if(cell_volume(node,[ii,jj,kk])<=0.or.covering_child(node,[ii,jj,kk])/=0) cycle
          value=candidate%fields(node)%values(:,ii,jj,kk)
          if(all(value==baseline%fields(node)%values(:,ii,jj,kk))) cycle
          good=all(ieee_is_finite(value)).and.value(1)>0.and.all(value(6:active_components)>=0)
          if(good) call recover_repair_cell(species,geometry(node),[ii,jj,kk],value, &
            temps%fields(node)%values(1,ii,jj,kk),qp,tg,sound,good,frame, &
          repair_node=node)
          local_ok=local_ok.and.good
          if(tree%nodes(node)%level>minimum_level) then
            if(tree%nodes(parent)%level<=minimum_level) cycle
            if(.not.present(rereflux_delta)) cycle
            call add_mpi_eb_delta_cell_3d(value,-routed_register%fields(node)%values(:,ii,jj,kk),good)
          else
            call cell_increment(node,[ii,jj,kk],[0,0,0],0,delta)
            call add_mpi_eb_delta_cell_3d(value,delta,good)
          endif
          if(good) good=value(1)>0.and.all(value(6:active_components)>=0)
          if(good) call recover_repair_cell(species,geometry(node),[ii,jj,kk],value, &
            temps%fields(node)%values(1,ii,jj,kk),qp,tg,sound,good,frame, &
          repair_node=node)
          local_ok=local_ok.and.good
        enddo
        enddo
        enddo
      enddo
      call agree(local_ok,comm,valid)
    end subroutine qualify_ancestor_route
    subroutine pending_coarse_donor(raw,delta,budget,valid)
      real(dp), intent(in) :: raw(nv),delta(nv)
      real(dp), intent(out) :: budget(nv)
      logical, intent(out) :: valid
      budget=raw
      call add_mpi_eb_delta_cell_3d(budget,delta,valid)
      if(valid) valid=budget(1)>0.and.all(budget(6:active_components)>=0)
    end subroutine pending_coarse_donor
    subroutine restore_coarse_donor(node,cell,post,delta,guess,pre,recovered,valid)
      integer, intent(in) :: node,cell(3)
      real(dp), intent(in) :: post(nv),delta(nv),guess
      real(dp), intent(out) :: pre(nv),recovered
      logical, intent(out) :: valid
      real(dp) :: post_temp,sound
      call recover_repair_cell(species,geometry(node),cell,post,guess,primitive,post_temp,sound,valid,frame, &
          repair_node=node)
      if(.not.valid) return
      pre=post
      call add_mpi_eb_delta_cell_3d(pre,-delta,valid)
      if(.not.valid) return
      valid=pre(1)>0.and.all(pre(6:active_components)>=0)
      if(.not.valid) return
      call recover_repair_cell(species,geometry(node),cell,pre,guess,primitive,recovered,sound,valid,frame, &
          repair_node=node)
    end subroutine restore_coarse_donor
    subroutine propagate_coarse_route(packet,a,step,distance,valid)
      real(dp), intent(in) :: packet(nv)
      integer, intent(in) :: a,step,distance
      logical, intent(out) :: valid
      type(amr_tree_state_3d) :: saved_candidate
      type(eb_integrated_flux_3d), allocatable :: saved_working(:),saved_effective(:)
      logical, allocatable :: saved_changed(:)
      real(dp) :: low(nv),high(nv),maximum
      integer :: h,node,next_node,cell(3),next_cell(3),face(3),mirror(3),capacity,count,ios,p
      logical :: connected,ready,local_ok
      valid=.false.;capacity=0
      node=parent;cell=coarse_cell
      do h=1,distance
        call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
        if(covering_child(next_node,next_cell)/=0) then
          low=-real(step,dp)*packet;high=low
          if(h==distance) then
            low=0;high=0
            if(step>0) low=-packet
            if(step<0) high=packet
          endif
          count=transfer_cell_count(next_node,next_cell)
          if(count<0.or.count>(huge(capacity)-capacity)/6) return
          capacity=capacity+6*count
          maximum=1
          call propagate_cell_packet(next_node,next_cell,a,low,high,.true.,local_ok,maximum)
          call agree(local_ok.and.maximum>=1,comm,ready)
          if(.not.ready) return
        endif
        node=next_node;cell=next_cell
      enddo
      if(capacity==0) then
        valid=.true.;return
      endif
      allocate(packet_keys(5,capacity),packet_values(nv,capacity),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) then
        if(allocated(packet_keys)) deallocate(packet_keys)
        if(allocated(packet_values)) deallocate(packet_values)
        return
      endif
      packet_count=0
      saved_candidate=candidate;saved_working=working;saved_effective=effective;saved_changed=changed
      node=parent;cell=coarse_cell
      do h=1,distance
        call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
        if(covering_child(next_node,next_cell)/=0) then
          low=-real(step,dp)*packet;high=low
          if(h==distance) then
            low=0;high=0
            if(step>0) low=-packet
            if(step<0) high=packet
          endif
          call propagate_cell_packet(next_node,next_cell,a,low,high,.true.,local_ok)
          call agree(local_ok,comm,ready)
          if(.not.ready) exit
        endif
        node=next_node;cell=next_cell
      enddo
      local_ok=ready
      if(local_ok) then
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.changed(p)) cycle
          if(any(.not.ieee_is_finite(candidate%fields(p)%values)).or. &
            any(.not.ieee_is_finite(working(p)%x)).or.any(.not.ieee_is_finite(working(p)%y)).or. &
            any(.not.ieee_is_finite(working(p)%z)).or.any(.not.ieee_is_finite(effective(p)%x)).or. &
            any(.not.ieee_is_finite(effective(p)%y)).or.any(.not.ieee_is_finite(effective(p)%z))) local_ok=.false.
        enddo
      endif
      call agree(local_ok,comm,valid)
      if(valid.and.tree%nodes(parent)%level/=minimum_level) call qualify_ancestor_route(saved_candidate,valid)
      if(.not.valid) then
        candidate=saved_candidate;working=saved_working;effective=saved_effective;changed=saved_changed
      endif
      deallocate(packet_keys,packet_values)
    end subroutine propagate_coarse_route

    ! Tighten the existing conservative face blend when a nearly empty
    ! coarse endpoint cannot satisfy the unchanged mixture EOS contract.
    ! Donor/species bounds stay enforced by the same interval kernel.
    subroutine constrain_coarse_eos(good)
      logical, intent(out) :: good
      real(dp), allocatable :: best_faces(:,:),best_coarse(:),best_fine(:,:)
      real(dp) :: low,high,cap,best_theta
      logical :: accepted,found
      integer :: attempt
      call qualify_coarse_endpoint(good)
      if(good) return
      low=theta_floor;high=theta;found=.false.
      do attempt=1,64
        cap=low+.5_dp*(high-low)
        if(cap<=low.or.cap>=high) exit
        call reconcile_reactive_eb_interface_budget_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
          transfers(:,:nf),weights(:nf),corrected,dc,df,theta,accepted,theta_floor,cap, &
        passive_components=passive_components)
        if(.not.accepted) then
          low=cap
          cycle
        end if
        call qualify_coarse_endpoint(accepted)
        if(accepted) then
          found=.true.;best_theta=theta
          best_faces=corrected;best_coarse=dc;best_fine=df
          low=cap
        else
          high=cap
        end if
      end do
      good=found
      if(.not.good) return
      corrected=best_faces;dc=best_coarse;df=best_fine;theta=best_theta
    end subroutine constrain_coarse_eos
    subroutine qualify_coarse_endpoint(good)
      logical, intent(out) :: good
      type(eb_integrated_flux_3d) :: trial_flux
      real(dp), allocatable :: delta(:,:,:,:)
      logical, allocatable :: refined(:,:,:)
      real(dp) :: trial_value(nv),qp(reactive_nprim(size(species))),tg,sound
      integer :: child,bxlo(3),bxhi(3),d(3),ios
      good=.false.
      if(rank==coarse_owner) then
        d=tree%nodes(parent)%shape
        allocate(refined(d(1),d(2),d(3)),stat=ios)
        if(ios==0) then
          refined=.false.
          do child=1,size(owners)
            if(tree%nodes(child)%parent/=parent) cycle
            call bounds(tree,child,bxlo,bxhi)
            refined(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))=.true.
          end do
          trial_flux=effective(parent)
          call set_face(trial_flux,axis,face_index,orientation*sum(corrected,dim=2)/coarse_area)
          call compute_reactive_eb_reflux_increment_3d(geometry(parent),refined,working(parent),trial_flux,delta, &
            good,all_faces=.true.,redistribution_delta=routed_register%fields(parent)%values)
          if(good) then
            trial_value=candidate%fields(parent)%values(:,coarse_cell(1),coarse_cell(2),coarse_cell(3))+ &
              delta(:,coarse_cell(1),coarse_cell(2),coarse_cell(3))
            call recover_repair_cell(species,geometry(parent),coarse_cell,trial_value, &
              coarse_guess,qp,tg,sound,good,frame, &
          repair_node=parent)
          end if
        end if
      end if
      call MPI_Bcast(good,1,MPI_LOGICAL,coarse_owner,comm,status)
      if(status/=MPI_SUCCESS) good=.false.
    end subroutine qualify_coarse_endpoint
    ! Breadth-first fallback over open same-level faces. Geometry and route
    ! metadata are shared; only the donor endpoint vector is broadcast by MPI.
    recursive subroutine search_connected_donor(target,found,valid,allow_partial,excluded_edges,fixed_donor,rejected_route)
      integer, intent(in) :: target
      integer, optional, intent(in) :: excluded_edges(:,:),fixed_donor(4)
      logical, optional, intent(out) :: rejected_route
      logical, intent(out) :: found,valid
      logical, optional, intent(in) :: allow_partial
      integer, allocatable :: queue_node(:),queue_cell(:,:),previous(:),edge_axis(:),edge_step(:)
      integer, allocatable :: new_node(:),new_cell(:,:),new_previous(:),new_axis(:),new_step(:)
      integer(int64), allocatable :: node_offset(:),visited_bits(:)
      integer :: capacity,limit,next_capacity,word,bit,head,tail,entry,node,next_node,cell(3),next_cell(3),face(3),next_face(3)
      integer :: peer,a_index,a,step,direction,r,ios,path_entry,bcast_status
      integer(int64) :: total,product_cells,vertex
      integer :: excluded_peer,excluded_cell(3),excluded_face(3),excluded_mirror(3),blocked
      logical :: edge_ok,protected,ready,alternative_found
      found=.false.;valid=.false.;total=0
      if(present(rejected_route)) rejected_route=.false.
      if(present(excluded_edges)) then
        if(size(excluded_edges,1)/=6) return
      endif
      allocate(node_offset(size(tree%nodes)),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      node_offset=-1
      do peer=1,size(tree%nodes)
        if(tree%nodes(peer)%level/=minimum_level+1) cycle
        node_offset(peer)=total;product_cells=1
        do a=1,3
          if(product_cells>int(huge(capacity),int64)/tree%nodes(peer)%shape(a)) return
          product_cells=product_cells*tree%nodes(peer)%shape(a)
        enddo
        if(total>int(huge(capacity),int64)-product_cells) return
        total=total+product_cells
      enddo
      if(total<1) return
      limit=int(total);capacity=min(limit,64)
      ! One bit per possible vertex replaces the scan through every queued
      ! cell. Route arrays grow only when another frontier slot is needed.
      allocate(visited_bits(int((total+63_int64)/64_int64)),queue_node(capacity),queue_cell(3,capacity), &
        previous(capacity),edge_axis(capacity),edge_step(capacity),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      visited_bits=0
      head=1;tail=1;queue_node(1)=child;queue_cell(:,1)=positions(:,target);previous(1)=0
      edge_axis(1)=0;edge_step(1)=0;cell=positions(:,target)
      vertex=node_offset(child)+(int(cell(3)-1,int64)*int(tree%nodes(child)%shape(2),int64)+cell(2)-1)* &
        int(tree%nodes(child)%shape(1),int64)+cell(1)-1
      word=int(vertex/64_int64)+1;bit=int(modulo(vertex,64_int64))
      visited_bits(word)=ibset(visited_bits(word),bit)
      do while(head<=tail)
        node=queue_node(head);cell=queue_cell(:,head)
        do a_index=1,3
          a=order(a_index)
          do direction=0,1
            step=2*direction-1
            if(a==axis) step=step*(2*side-1)
            call path_neighbor(node,cell,a,step,next_node,next_cell,face,next_face,edge_ok)
            if(.not.edge_ok) cycle
            if(present(excluded_edges)) then
              protected=.false.
              do blocked=1,size(excluded_edges,2)
                if(a/=excluded_edges(5,blocked)) cycle
                if(node==excluded_edges(1,blocked).and.all(cell==excluded_edges(2:4,blocked)).and. &
                  step==excluded_edges(6,blocked)) then
                  protected=.true.;exit
                endif
                call path_neighbor(excluded_edges(1,blocked),excluded_edges(2:4,blocked), &
                  excluded_edges(5,blocked),excluded_edges(6,blocked),excluded_peer,excluded_cell, &
                  excluded_face,excluded_mirror,edge_ok)
                if(.not.edge_ok) return
                if(node==excluded_peer.and.all(cell==excluded_cell).and.step==-excluded_edges(6,blocked)) then
                  protected=.true.;exit
                endif
              enddo
              if(protected) cycle
            endif
            vertex=node_offset(next_node)+(int(next_cell(3)-1,int64)* &
              int(tree%nodes(next_node)%shape(2),int64)+next_cell(2)-1)* &
              int(tree%nodes(next_node)%shape(1),int64)+next_cell(1)-1
            word=int(vertex/64_int64)+1;bit=int(modulo(vertex,64_int64))
            if(btest(visited_bits(word),bit)) cycle
            if(tail==capacity) then
              if(capacity==limit) return
              next_capacity=limit
              if(capacity<=limit/2) next_capacity=2*capacity
              allocate(new_node(next_capacity),new_cell(3,next_capacity),new_previous(next_capacity), &
                new_axis(next_capacity),new_step(next_capacity),stat=ios)
              call agree(ios==0,comm,ready)
              if(.not.ready) return
              new_node(:tail)=queue_node;new_cell(:,:tail)=queue_cell;new_previous(:tail)=previous
              new_axis(:tail)=edge_axis;new_step(:tail)=edge_step
              call move_alloc(new_node,queue_node);call move_alloc(new_cell,queue_cell)
              call move_alloc(new_previous,previous);call move_alloc(new_axis,edge_axis)
              call move_alloc(new_step,edge_step);capacity=next_capacity
            endif
            visited_bits(word)=ibset(visited_bits(word),bit)
            tail=tail+1;queue_node(tail)=next_node;queue_cell(:,tail)=next_cell;previous(tail)=head
            edge_axis(tail)=a;edge_step(tail)=step
            protected=.false.
            do r=1,nf
              if(next_node==child.and.all(next_cell==positions(:,r))) protected=.true.
            enddo
            do r=1,size(used_donors,2)
              if(all(used_donors(:,r)==[next_node,next_cell])) protected=.true.
            enddo
            if(protected) cycle
            if(present(fixed_donor)) then
              if(any([next_node,next_cell]/=fixed_donor)) cycle
            endif
            donor_volume(target)=cell_volume(next_node,next_cell)
            if(rank==owners(next_node)) donor_state(:,target)= &
              candidate%fields(next_node)%values(:,next_cell(1),next_cell(2),next_cell(3))
            call MPI_Bcast(donor_state(:,target),nv,MPI_DOUBLE_PRECISION,owners(next_node),comm,bcast_status)
            if(bcast_status/=MPI_SUCCESS) return
            call supply_reactive_eb_interface_donors_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
              transfers(:,:nf),weights(:nf),donor_state(:,:nf),donor_volume(:nf),funded_fine,funded_donor,exchange, &
              edge_ok,target_subface=target,allow_partial=allow_partial, &
              force_exchange=covering_child(child,positions(:,target))/=0, &
        passive_components=passive_components)
            call agree(edge_ok,comm,ready)
            if(.not.ready) cycle
            donor_node=next_node;donor_cell=next_cell
            if(allocated(route_node)) deallocate(route_node,route_cell,route_axis,route_step)
            route_length=0;path_entry=tail
            do while(previous(path_entry)>0)
              route_length=route_length+1;path_entry=previous(path_entry)
            enddo
            allocate(route_node(route_length),route_cell(3,route_length),route_axis(route_length), &
              route_step(route_length),stat=ios)
            call agree(ios==0,comm,ready)
            if(.not.ready) return
            path_entry=tail
            do entry=route_length,1,-1
              route_axis(entry)=edge_axis(path_entry);route_step(entry)=edge_step(path_entry)
              path_entry=previous(path_entry)
              route_node(entry)=queue_node(path_entry);route_cell(:,entry)=queue_cell(:,path_entry)
            enddo
            call route_accepts_packet(target,edge_ok,allow_partial)
            if(.not.edge_ok.or..not.route_improved) then
              if(present(rejected_route)) rejected_route=route_can_reroute.and.any(exchange(:,target)/=0)
              if(.not.present(fixed_donor).and.route_can_reroute.and.any(exchange(:,target)/=0)) then
                call find_alternate_donor_route(target,alternative_found,valid,allow_partial)
                if(.not.valid) return
                if(alternative_found) then
                  found=.true.;return
                endif
              endif
              cycle
            endif
            if(route_complete) call minimize_recipient_packet(target,edge_ok)
            if(.not.edge_ok) cycle
            found=.true.;valid=.true.;return
          enddo
        enddo
        head=head+1
      enddo
      valid=.true.
    end subroutine search_connected_donor
    ! Explore canonical exclusion sets in breadth-first order. A rejected
    ! rerouted path contributes its own edges, so a feasible route may avoid
    ! multiple independently inadmissible corridors. Duplicate sets are
    ! visited once; gas/flux publication remains outside this metadata search.
    recursive subroutine find_alternate_donor_route(target,found,valid,allow_partial)
      integer, intent(in) :: target
      logical, intent(out) :: found,valid
      logical, optional, intent(in) :: allow_partial
      type :: retry_state
        integer, allocatable :: blocked(:,:),path(:,:)
      end type retry_state
      type(retry_state), allocatable :: frontier(:),grown(:)
      integer, allocatable :: saved_node(:),saved_cell(:,:),saved_axis(:),saved_step(:),blocked_edges(:,:)
      real(dp) :: saved_fine(nv),saved_donor(nv),saved_exchange(nv),missing(nv),allowance
      integer :: edge,count,ios,endpoint(4),key(6),temporary(6),head,tail,capacity,next_capacity,nblocked,entry,prior
      integer :: incident_node,incident_cell(3),incident_face(3),incident_mirror(3)
      logical :: ready,rejected,duplicate
      found=.false.;valid=.false.;count=route_length;endpoint=[donor_node,donor_cell]
      call recipient_allowance(target,allowance,valid,missing)
      if(.not.valid) return
      ! A covered endpoint receives a nonnegative geometric distribution of
      ! this same packet, independently of the remaining route. If an actual
      ! density/species deficit exists and no deficient component receives a
      ! positive increment, changing incident faces cannot improve its budget.
      ! Positive total mass alone does not replenish an absent species. Keep
      ! the prior conservative gate when no physical deficit is established;
      ! EOS-only rejection is not such a certificate. Held energy is signed
      ! and does not participate in this species-admission certificate.
      if(missing(1)>0.or.any(missing(6:active_components)>0)) then
        if(.not.(missing(1)>0.and.exchange(1,target)>0).and. &
          .not.any(missing(6:active_components)>0.and.exchange(6:active_components,target)>0)) return
      else
        if(exchange(1,target)<=0) return
      endif
      valid=.false.
      if(count>huge(count)/6) return
      capacity=16
      allocate(saved_node(count),saved_cell(3,count),saved_axis(count),saved_step(count),frontier(capacity), &
        stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) then
        valid=.false.;return
      endif
      saved_node=route_node;saved_cell=route_cell;saved_axis=route_axis;saved_step=route_step
      saved_fine=funded_fine(:,target);saved_donor=funded_donor(:,target);saved_exchange=exchange(:,target)
      allocate(frontier(1)%blocked(6,0),frontier(1)%path(6,count),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) then
        valid=.false.;return
      endif
      do edge=1,count
        frontier(1)%path(:,edge)=[saved_node(edge),saved_cell(:,edge),saved_axis(edge),saved_step(edge)]
      enddo
      head=1;tail=1
      do while(head<=tail)
        nblocked=size(frontier(head)%blocked,2)
        if(nblocked>=huge(nblocked)/6) then
          valid=.false.;return
        endif
        do edge=1,size(frontier(head)%path,2)
          ! Rerouting ordinary exposed interiors changes only the face ledger:
          ! their incoming/outgoing complete packets cancel without a state
          ! update. Only edges incident to a covered coarse cell can change
          ! descendant spreading and the admission of this rejected packet.
          ! Retain all such edges, including multiple covered corridors.
          call path_neighbor(frontier(head)%path(1,edge),frontier(head)%path(2:4,edge), &
            frontier(head)%path(5,edge),frontier(head)%path(6,edge),incident_node,incident_cell, &
            incident_face,incident_mirror,valid)
          if(.not.valid) return
          if(covering_child(frontier(head)%path(1,edge),frontier(head)%path(2:4,edge))==0.and. &
            covering_child(incident_node,incident_cell)==0) cycle
          call canonical_route_edge(frontier(head)%path(:,edge),key,valid)
          if(.not.valid) return
          duplicate=.false.
          do entry=1,nblocked
            if(all(key==frontier(head)%blocked(:,entry))) duplicate=.true.
          enddo
          if(duplicate) cycle
          allocate(blocked_edges(6,nblocked+1),stat=ios)
      call agree(ios==0,comm,ready)
          if(.not.ready) then
        valid=.false.;return
      endif
          blocked_edges(:,:nblocked)=frontier(head)%blocked;blocked_edges(:,nblocked+1)=key
          ! Insertion maintains one representation for every undirected set.
          do entry=nblocked+1,2,-1
            if(.not.route_edge_less(blocked_edges(:,entry),blocked_edges(:,entry-1))) exit
            temporary=blocked_edges(:,entry-1);blocked_edges(:,entry-1)=blocked_edges(:,entry)
            blocked_edges(:,entry)=temporary
          enddo
          duplicate=.false.
          do prior=1,tail
            if(size(frontier(prior)%blocked,2)/=nblocked+1) cycle
            if(all(frontier(prior)%blocked==blocked_edges)) then
              duplicate=.true.;exit
            endif
          enddo
          if(duplicate) then
            deallocate(blocked_edges);cycle
          endif
          if(tail==capacity) then
            if(capacity==huge(capacity)) then
              valid=.false.;return
            endif
            next_capacity=huge(capacity)
            if(capacity<=huge(capacity)/2) next_capacity=2*capacity
            allocate(grown(next_capacity),stat=ios)
      call agree(ios==0,comm,ready)
            if(.not.ready) then
        valid=.false.;return
      endif
            do prior=1,tail
              call move_alloc(frontier(prior)%blocked,grown(prior)%blocked)
              call move_alloc(frontier(prior)%path,grown(prior)%path)
            enddo
            call move_alloc(grown,frontier);capacity=next_capacity
          endif
          tail=tail+1;call move_alloc(blocked_edges,frontier(tail)%blocked)
          call search_connected_donor(target,found,valid,allow_partial,frontier(tail)%blocked,endpoint,rejected)
          if(.not.valid.or.found) return
          count=0
          if(rejected) count=route_length
          if(count>huge(count)/6) then
            valid=.false.;return
          endif
          allocate(frontier(tail)%path(6,count),stat=ios)
      call agree(ios==0,comm,ready)
          if(.not.ready) then
        valid=.false.;return
      endif
          do entry=1,count
            frontier(tail)%path(:,entry)=[route_node(entry),route_cell(:,entry),route_axis(entry),route_step(entry)]
          enddo
        enddo
        head=head+1
      enddo
      ! Restore the failed original candidate before resuming other donors.
      count=size(saved_node)
      if(allocated(route_node)) deallocate(route_node,route_cell,route_axis,route_step)
      allocate(route_node(count),route_cell(3,count),route_axis(count),route_step(count),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) then
        valid=.false.;return
      endif
      route_node=saved_node;route_cell=saved_cell;route_axis=saved_axis;route_step=saved_step;route_length=count
      donor_node=endpoint(1);donor_cell=endpoint(2:4)
      funded_fine(:,target)=saved_fine;funded_donor(:,target)=saved_donor;exchange(:,target)=saved_exchange
      valid=.true.
    end subroutine find_alternate_donor_route
    subroutine canonical_route_edge(edge,key,valid)
      integer, intent(in) :: edge(6)
      integer, intent(out) :: key(6)
      logical, intent(out) :: valid
      integer :: next_node,next_cell(3),face(3),mirror(3),reverse(6)
      call path_neighbor(edge(1),edge(2:4),edge(5),edge(6),next_node,next_cell,face,mirror,valid)
      if(.not.valid) return
      reverse=[next_node,next_cell,edge(5),-edge(6)];key=edge
      if(route_edge_less(reverse,key)) key=reverse
    end subroutine canonical_route_edge
    pure logical function route_edge_less(left,right) result(less)
      integer, intent(in) :: left(6),right(6)
      integer :: component
      less=.false.
      do component=1,6
        if(left(component)==right(component)) cycle
        less=left(component)<right(component);return
      enddo
    end function route_edge_less
    recursive logical function regular_route_endpoint(node,cell) result(regular)
      integer, intent(in) :: node,cell(3)
      integer :: descendant,lo(3),hi(3),first(3),last(3),r,ii,jj,kk,a,face(3)
      regular=.false.
      if(geometry(node)%volume_fraction(cell(1),cell(2),cell(3))/=1) return
      do a=1,3
        face=cell;face(a)=face(a)-1
        if(open_area(node,a,face)/=open_area_full(node,a)) return
        if(open_area(node,a,cell)/=open_area_full(node,a)) return
      enddo
      descendant=covering_child(node,cell)
      if(descendant/=0) then
        call bounds(tree,descendant,lo,hi)
        r=tree%nodes(descendant)%patch%refinement_ratio
        first=(cell-lo)*r+1;last=first+r-1
        do kk=first(3),last(3)
        do jj=first(2),last(2)
        do ii=first(1),last(1)
          if(.not.regular_route_endpoint(descendant,[ii,jj,kk])) return
        enddo
        enddo
        enddo
      endif
      regular=.true.
    end function regular_route_endpoint
    real(dp) function open_area_full(node,a) result(area)
      integer, intent(in) :: node,a
      select case(a)
      case(1);area=geometry(node)%dy*geometry(node)%dz
      case(2);area=geometry(node)%dx*geometry(node)%dz
      case(3);area=geometry(node)%dx*geometry(node)%dy
      end select
    end function open_area_full
    ! Fully regular endpoint spreading is uniform in every direction.
    ! A failed endpoint cannot be repaired by changing its incident edge.
    subroutine exclude_regular_endpoint_reroute(node,cell,target,regular)
      integer, intent(in) :: node,cell(3),target
      logical, intent(in) :: regular(2)
      integer :: endpoint,root_node,root_cell(3),ancestor,index(3),lo(3),hi(3),parent_node,r
      do endpoint=1,2
        root_node=child;root_cell=positions(:,target)
        if(endpoint==2) then
          root_node=donor_node;root_cell=donor_cell
        endif
        if(.not.regular(endpoint)) cycle
        ancestor=node;index=cell
        do while(ancestor/=root_node)
          parent_node=tree%nodes(ancestor)%parent
          if(parent_node==0) exit
          call bounds(tree,ancestor,lo,hi)
          r=tree%nodes(ancestor)%patch%refinement_ratio
          index=lo+(index-1)/r;ancestor=parent_node
        enddo
        if(ancestor==root_node.and.all(index==root_cell)) route_can_reroute=.false.
      enddo
    end subroutine exclude_regular_endpoint_reroute
    subroutine recipient_allowance(target,maximum,valid,deficit)
      integer, intent(in) :: target
      real(dp), intent(out) :: maximum
      logical, intent(out) :: valid
      real(dp), optional, intent(out) :: deficit(nv)
      real(dp) :: low(nv),high(nv)
      low=0;high=0;maximum=1
      if(present(deficit)) deficit=0
      if(side==0) low=transfers(:,target)-weights(target)*old
      if(side==1) high=weights(target)*old-transfers(:,target)
      call propagate_cell_packet(child,positions(:,target),axis,low,high,.true.,valid,maximum,deficit)
    end subroutine recipient_allowance
    subroutine build_straight_route(target,length,a,step,valid)
      integer, intent(in) :: target,length,a,step
      logical, intent(out) :: valid
      integer :: entry,node,cell(3),next_node,next_cell(3),face(3),next_face(3),ios
      logical :: ready
      valid=.false.
      if(allocated(route_node)) deallocate(route_node,route_cell,route_axis,route_step)
      allocate(route_node(length),route_cell(3,length),route_axis(length),route_step(length),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      route_length=length;node=child;cell=positions(:,target)
      do entry=1,length
        route_node(entry)=node;route_cell(:,entry)=cell;route_axis(entry)=a;route_step(entry)=step
        call path_neighbor(node,cell,a,step,next_node,next_cell,face,next_face,valid)
        if(.not.valid) return
        node=next_node;cell=next_cell
      enddo
      valid=.true.
    end subroutine build_straight_route
    ! Test the complete donor/recipient/interior change before selecting a
    ! route. Trial descendants and faces are restored on either outcome, so
    ! a rejected route does not consume inventories needed by later searches.
    recursive subroutine route_accepts_packet(target,valid,allow_partial)
      integer, intent(in) :: target
      logical, intent(out) :: valid
      logical, optional, intent(in) :: allow_partial
      type(amr_tree_state_3d) :: saved_candidate
      type(eb_integrated_flux_3d), allocatable :: saved_working(:),saved_effective(:)
      logical, allocatable :: saved_changed(:)
      real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed,value(nv),lower,upper,maximum,baseline(nv),direction(nv)
      real(dp) :: before_allowance,after_allowance,prior_allowance,before_deficit(nv),after_deficit(nv)
      logical :: allowance_ok,endpoint_regular(2)
      integer :: node,cell(3),d(3),ii,jj,kk,entry,prior,status
      logical :: good,local_ok,needs_trial,propagation_ok,interval_ok
      route_complete=.false.;route_can_reroute=.false.
      route_improved=covering_child(child,positions(:,target))==0
      before_allowance=0
      if(.not.route_improved) then
        call recipient_allowance(target,before_allowance,allowance_ok,before_deficit)
        if(.not.allowance_ok) then
          valid=.false.;return
        endif
      endif
      local_ok=.true.
      do entry=1,2
        node=child;cell=positions(:,target);value=funded_fine(:,target)
        if(entry==2) then
          node=donor_node;cell=donor_cell;value=funded_donor(:,target)
        endif
        if(rank/=owners(node)) cycle
        if(.not.all(ieee_is_finite(value)).or.value(1)<=0.or.any(value(6:active_components)<0)) then
          local_ok=.false.;cycle
        endif
        call recover_repair_cell(species,geometry(node),cell,value, &
          temperature%fields(node)%values(1,cell(1),cell(2),cell(3)),primitive,recovered,speed,good,frame, &
          repair_node=node)
        local_ok=local_ok.and.good
      enddo
      call agree(local_ok,comm,valid)
      if(.not.valid) then
        if(present(allow_partial)) then
          if(allow_partial) call find_partial_route_packet(target,1._dp,valid)
        endif
        return
      endif
      needs_trial=covering_child(donor_node,donor_cell)/=0.or.covering_child(child,positions(:,target))/=0
      do entry=2,route_length
        if(covering_child(route_node(entry),route_cell(:,entry))/=0) needs_trial=.true.
      enddo
      ! Only geometrically different endpoint spreading or a covered
      ! interior can change acceptance when the donor/path is rerouted.
      endpoint_regular(1)=regular_route_endpoint(child,positions(:,target))
      endpoint_regular(2)=regular_route_endpoint(donor_node,donor_cell)
      route_can_reroute=(covering_child(child,positions(:,target))/=0.and..not.endpoint_regular(1)).or. &
        (covering_child(donor_node,donor_cell)/=0.and..not.endpoint_regular(2))
      do entry=2,route_length
        if(covering_child(route_node(entry),route_cell(:,entry))/=0) route_can_reroute=.true.
      enddo
      ! Ordinary same-level routes leave every interior state unchanged.
      if(.not.needs_trial) return
      saved_candidate=candidate;saved_working=working;saved_effective=effective;saved_changed=changed
      cell=positions(:,target)
      if(rank==owners(child).and.covering_child(child,cell)==0) then
        candidate%fields(child)%values(:,cell(1),cell(2),cell(3))=funded_fine(:,target)
        changed(child)=.true.
      endif
      if(rank==owners(donor_node).and.covering_child(donor_node,donor_cell)==0) then
        candidate%fields(donor_node)%values(:,donor_cell(1),donor_cell(2),donor_cell(3))=funded_donor(:,target)
        changed(donor_node)=.true.
      endif
      call propagate_route_interior(exchange(:,target),good)
      local_ok=good;propagation_ok=good;maximum=1
      if(good.and.covering_child(child,positions(:,target))/=0) then
        call recipient_allowance(target,after_allowance,allowance_ok,after_deficit)
        ! A different species can remain the limiting constraint while this
        ! packet reduces another deficit. Accept that complementary progress
        ! only when no component's missing physical inventory increases.
        route_improved=allowance_ok.and.(after_allowance>before_allowance.or. &
          (all(after_deficit<=before_deficit).and.any(after_deficit<before_deficit)))
      endif
      ! A route can redistribute descendants of an earlier recipient.
      ! Preserve its funded correction budget before spending this packet.
      if(good) then
        do prior=1,target-1
          call recipient_allowance(prior,prior_allowance,good)
          if(.not.good.or.prior_allowance<required_correction) local_ok=.false.
        enddo
      endif
      if(propagation_ok) then
        do node=1,size(tree%nodes)
          if(rank/=owners(node)) cycle
          if(.not.changed(node)) cycle
          d=tree%nodes(node)%shape
          do kk=1,d(3)
          do jj=1,d(2)
          do ii=1,d(1)
            if(geometry(node)%volume_fraction(ii,jj,kk)<=0) cycle
            associate(trial_value=>candidate%fields(node)%values(:,ii,jj,kk))
              if(all(trial_value==saved_candidate%fields(node)%values(:,ii,jj,kk))) cycle
              baseline=saved_candidate%fields(node)%values(:,ii,jj,kk);direction=trial_value-baseline
              lower=0;upper=1
              call reactive_eb_exchange_interval_3d(reshape(baseline,[nv,1]),reshape(direction,[nv,1]), &
                lower,upper,interval_ok, &
        passive_components=passive_components)
              if(.not.interval_ok.or.lower>0) upper=0
              maximum=min(maximum,upper)
              if(upper<1) call exclude_regular_endpoint_reroute(node,[ii,jj,kk],target,endpoint_regular)
              if(.not.all(ieee_is_finite(trial_value)).or.trial_value(1)<=0.or.any(trial_value(6:active_components)<0)) then
                local_ok=.false.;cycle
              endif
              call recover_repair_cell(species,geometry(node),[ii,jj,kk],trial_value, &
                temperature%fields(node)%values(1,ii,jj,kk),primitive,recovered,speed,good,frame, &
          repair_node=node)
              local_ok=local_ok.and.good
              if(.not.good) call exclude_regular_endpoint_reroute(node,[ii,jj,kk],target,endpoint_regular)
            end associate
          enddo
          enddo
          enddo
        enddo
      endif
      call agree(local_ok,comm,good)
      if(good.and.covering_child(child,positions(:,target))/=0) then
        if(allowance_ok.and.after_allowance>=required_correction) then
          do prior=1,target
            call project_recipient_endpoint(prior,good)
            if(.not.good) exit
          enddo
          if(good) call accept_route_candidate(saved_candidate,good)
          local_ok=good;route_complete=good
        endif
      endif
      ! Propagation failures may leave a partially populated private ledger.
      if(allocated(packet_keys)) deallocate(packet_keys)
      if(allocated(packet_values)) deallocate(packet_values)
      candidate=saved_candidate;working=saved_working;effective=saved_effective;changed=saved_changed
      call agree(local_ok,comm,valid)
      call agree(route_can_reroute,comm,good)
      route_can_reroute=good
      ! A parent average may admit more withdrawal than its least stocked
      ! descendant. In the pooling search, retain the physical fraction and
      ! verify its represented states and EOS again before spending it.
      call MPI_Allreduce(MPI_IN_PLACE,maximum,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) then
        valid=.false.;return
      endif
      if(valid.or..not.propagation_ok.or.maximum<=0) return
      if(.not.present(allow_partial)) return
      if(.not.allow_partial) return
      call find_partial_route_packet(target,maximum,valid)
    end subroutine route_accepts_packet
    ! Retain donor inventory once the covered recipient can accept its final
    ! interface correction. Each trial uses the same complete conservative
    ! direction and validates physical descendants plus the corrected EOS.
    subroutine minimize_recipient_packet(target,valid)
      integer, intent(in) :: target
      logical, intent(out) :: valid
      real(dp) :: original_fine(nv),original_donor(nv),original_exchange(nv)
      real(dp) :: lower,upper,middle
      integer :: iteration
      logical :: accepted
      original_fine=funded_fine(:,target);original_donor=funded_donor(:,target);original_exchange=exchange(:,target)
      lower=0;upper=1
      do iteration=1,64
        middle=lower+(upper-lower)/2
        if(middle<=lower.or.middle>=upper) exit
        exchange(:,target)=middle*original_exchange
        funded_fine(:,target)=fine_state(:,target)+exchange(:,target)/volumes(target)
        funded_donor(:,target)=donor_state(:,target)-exchange(:,target)/donor_volume(target)
        call route_accepts_packet(target,accepted)
        if(accepted.and.route_complete) then
          upper=middle
        else
          lower=middle
        endif
      enddo
      ! Retain a represented margin for the later group reconciliation,
      ! whose component sums may round differently from this endpoint trial.
      upper=min(1._dp,upper+128*epsilon(1._dp))
      exchange(:,target)=upper*original_exchange
      funded_fine(:,target)=fine_state(:,target)+exchange(:,target)/volumes(target)
      funded_donor(:,target)=donor_state(:,target)-exchange(:,target)/donor_volume(target)
      call route_accepts_packet(target,valid)
      valid=valid.and.route_complete.and.route_improved
      if(valid) return
      ! Preserve the already accepted full packet if represented endpoint
      ! rounding prevents acceptance of the reduced packet.
      funded_fine(:,target)=original_fine;funded_donor(:,target)=original_donor;exchange(:,target)=original_exchange
      call route_accepts_packet(target,valid)
      valid=valid.and.route_complete.and.route_improved
    end subroutine minimize_recipient_packet
    subroutine project_recipient_endpoint(target,valid)
      integer, intent(in) :: target
      logical, intent(out) :: valid
      real(dp) :: low(nv),high(nv)
      integer :: count,ios
      logical :: ready
      valid=.false.;count=transfer_cell_count(child,positions(:,target))
      if(count<0.or.count>huge(count)/2) return
      allocate(packet_keys(5,2*count),packet_values(nv,2*count),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      packet_count=0;low=0;high=0
      if(side==0) low=required_correction*(transfers(:,target)-weights(target)*old)
      if(side==1) high=required_correction*(weights(target)*old-transfers(:,target))
      call propagate_cell_packet(child,positions(:,target),axis,low,high,.true.,valid)
      deallocate(packet_keys,packet_values)
    end subroutine project_recipient_endpoint
    subroutine accept_route_candidate(baseline,valid)
      type(amr_tree_state_3d), intent(in) :: baseline
      logical, intent(out) :: valid
      integer :: node,d(3),ii,jj,kk
      real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed
      logical :: good,local_ok
      local_ok=.true.
      do node=1,size(tree%nodes)
        if(rank/=owners(node)) cycle
        if(.not.changed(node)) cycle
        d=tree%nodes(node)%shape
        do kk=1,d(3)
        do jj=1,d(2)
        do ii=1,d(1)
          if(geometry(node)%volume_fraction(ii,jj,kk)<=0) cycle
          associate(value=>candidate%fields(node)%values(:,ii,jj,kk))
            if(all(value==baseline%fields(node)%values(:,ii,jj,kk))) cycle
            if(.not.all(ieee_is_finite(value)).or.value(1)<=0.or.any(value(6:active_components)<0)) then
              local_ok=.false.;cycle
            endif
            call recover_repair_cell(species,geometry(node),[ii,jj,kk],value, &
              temperature%fields(node)%values(1,ii,jj,kk),primitive,recovered,speed,good,frame, &
          repair_node=node)
            local_ok=local_ok.and.good
          end associate
        enddo
        enddo
        enddo
      enddo
      call agree(local_ok,comm,valid)
    end subroutine accept_route_candidate
    ! Keep all conserved components on the same exchange direction. The
    ! density/species bound supplies an upper trial; EOS recovery then brackets
    ! a represented admissible fraction without altering its tolerances.
    recursive subroutine find_partial_route_packet(target,maximum,valid)
      integer, intent(in) :: target
      real(dp), intent(in) :: maximum
      logical, intent(out) :: valid
      real(dp) :: original_fine(nv),original_donor(nv),original_exchange(nv)
      real(dp) :: lower,upper,middle
      integer :: iteration
      logical :: accepted
      valid=.false.
      original_fine=funded_fine(:,target);original_donor=funded_donor(:,target);original_exchange=exchange(:,target)
      upper=maximum*(1-32*epsilon(1._dp));lower=0
      exchange(:,target)=(upper)*original_exchange
      funded_fine(:,target)=fine_state(:,target)+exchange(:,target)/volumes(target)
      funded_donor(:,target)=donor_state(:,target)-exchange(:,target)/donor_volume(target)
      call route_accepts_packet(target,accepted)
      if(accepted.and.any(exchange(:,target)/=0)) then
        valid=.true.;return
      endif
      exchange(:,target)=(0._dp)*original_exchange
      funded_fine(:,target)=fine_state(:,target)+exchange(:,target)/volumes(target)
      funded_donor(:,target)=donor_state(:,target)-exchange(:,target)/donor_volume(target)
      call route_accepts_packet(target,accepted)
      if(accepted) then
        do iteration=1,64
          middle=lower+(upper-lower)/2
          if(middle<=lower.or.middle>=upper) exit
          exchange(:,target)=(middle)*original_exchange
          funded_fine(:,target)=fine_state(:,target)+exchange(:,target)/volumes(target)
          funded_donor(:,target)=donor_state(:,target)-exchange(:,target)/donor_volume(target)
          call route_accepts_packet(target,accepted)
          if(accepted) then
            lower=middle
          else
            upper=middle
          endif
        enddo
        if(lower>0) then
          exchange(:,target)=(lower)*original_exchange
          funded_fine(:,target)=fine_state(:,target)+exchange(:,target)/volumes(target)
          funded_donor(:,target)=donor_state(:,target)-exchange(:,target)/donor_volume(target)
          call route_accepts_packet(target,valid)
          valid=valid.and.any(exchange(:,target)/=0)
        endif
      endif
      if(.not.valid) then
        funded_fine(:,target)=original_fine;funded_donor(:,target)=original_donor;exchange(:,target)=original_exchange
      endif
    end subroutine find_partial_route_packet
    ! A turning route can change descendant inventories even when its
    ! parent cell has zero net divergence. Propagate each axis once, with both
    ! incident faces assembled, and retain the transaction's final EOS gate.
    subroutine propagate_route_interior(packet,valid)
      real(dp), intent(in) :: packet(nv)
      logical, intent(out) :: valid
      integer :: entry,node,cell(3),a,count,capacity,ios
      real(dp) :: low(nv),high(nv)
      logical :: ready
      valid=.false.;capacity=0
      do entry=2,route_length
        node=route_node(entry);cell=route_cell(:,entry)
        if(covering_child(node,cell)==0) cycle
        count=transfer_cell_count(node,cell)
        if(count<0.or.count>(huge(capacity)-capacity)/6) return
        capacity=capacity+6*count
      enddo
      if(covering_child(child,route_cell(:,1))/=0) then
        count=transfer_cell_count(child,route_cell(:,1))
        if(count<0.or.count>(huge(capacity)-capacity)/2) return
        capacity=capacity+2*count
      endif
      if(covering_child(donor_node,donor_cell)/=0) then
        count=transfer_cell_count(donor_node,donor_cell)
        if(count<0.or.count>(huge(capacity)-capacity)/2) return
        capacity=capacity+2*count
      endif
      if(capacity==0) then
        valid=.true.;return
      endif
      allocate(packet_keys(5,capacity),packet_values(nv,capacity),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      packet_count=0
      if(covering_child(child,route_cell(:,1))/=0) then
        a=route_axis(1);low=0;high=0
        if(route_step(1)>0) high=-packet
        if(route_step(1)<0) low=packet
        call propagate_cell_packet(child,route_cell(:,1),a,low,high,.true.,valid)
        if(.not.valid) return
      endif
      do entry=2,route_length
        node=route_node(entry);cell=route_cell(:,entry)
        if(covering_child(node,cell)==0) cycle
        do a=1,3
          low=0;high=0
          if(route_axis(entry-1)==a) then
            if(route_step(entry-1)>0) low=-packet
            if(route_step(entry-1)<0) high=packet
          endif
          if(route_axis(entry)==a) then
            if(route_step(entry)>0) high=-packet
            if(route_step(entry)<0) low=packet
          endif
          if(all(low==0).and.all(high==0)) cycle
          call propagate_cell_packet(node,cell,a,low,high,.true.,valid)
          if(.not.valid) return
        enddo
      enddo
      if(covering_child(donor_node,donor_cell)/=0) then
        ! The last edge enters the donor. Its incident packet withdraws the
        ! exchange from every physical descendant and updates their faces.
        a=route_axis(route_length);low=0;high=0
        if(route_step(route_length)>0) low=-packet
        if(route_step(route_length)<0) high=packet
        call propagate_cell_packet(donor_node,donor_cell,a,low,high,.true.,valid)
        if(.not.valid) return
      endif
      deallocate(packet_keys,packet_values)
      valid=.true.
    end subroutine propagate_route_interior
    integer function covering_child(node,cell) result(covered_by)
      integer, intent(in) :: node,cell(3)
      integer :: peer,lower(3),upper(3)
      covered_by=0
      do peer=1,size(tree%nodes)
        if(tree%nodes(peer)%parent/=node) cycle
        call bounds(tree,peer,lower,upper)
        if(all(cell>=lower).and.all(cell<=upper)) then
          covered_by=peer;return
        endif
      enddo
    end function covering_child
    recursive integer function transfer_cell_count(node,cell) result(count)
      integer, intent(in) :: node,cell(3)
      integer :: covered_by,lower(3),upper(3),first(3),last(3),r,ii,jj,kk,addition
      count=1;covered_by=covering_child(node,cell)
      if(covered_by==0) return
      call bounds(tree,covered_by,lower,upper)
      r=tree%nodes(covered_by)%patch%refinement_ratio
      first=(cell-lower)*r+1;last=first+r-1
      do kk=first(3),last(3)
      do jj=first(2),last(2)
      do ii=first(1),last(1)
        if(cell_volume(covered_by,[ii,jj,kk])<=0) cycle
        addition=transfer_cell_count(covered_by,[ii,jj,kk])
        if(addition<0.or.count>huge(count)-addition) then
          count=-1;return
        endif
        count=count+addition
      enddo
      enddo
      enddo
    end function transfer_cell_count
    real(dp) function cell_volume(node,cell) result(volume)
      integer, intent(in) :: node,cell(3)
      volume=geometry(node)%volume_fraction(cell(1),cell(2),cell(3))* &
        geometry(node)%dx*geometry(node)%dy*geometry(node)%dz
    end function cell_volume
    ! Record each shared internal face once. Adjacent refined cells supply the
    ! same packet to their common face; a disagreement aborts the transaction.
    subroutine record_refined_packet(node,a,face,packet,valid)
      integer, intent(in) :: node,a,face(3)
      real(dp), intent(in) :: packet(nv)
      logical, intent(out) :: valid
      integer :: entry,key(5)
      valid=.true.
      if(all(packet==0)) return
      if(open_area(node,a,face)<=0) then
        valid=.false.;return
      endif
      key=[node,a,face]
      do entry=1,packet_count
        if(any(packet_keys(:,entry)/=key)) cycle
        valid=all(packet_values(:,entry)==packet);return
      enddo
      if(packet_count==size(packet_keys,2)) then
        valid=.false.;return
      endif
      packet_count=packet_count+1;packet_keys(:,packet_count)=key;packet_values(:,packet_count)=packet
      if(rank==owners(node)) call add_path_packet(node,a,face,packet)
    end subroutine record_refined_packet
    recursive subroutine propagate_cell_packet(node,cell,a,low_packet,high_packet,root,valid,maximum_correction,deficit)
      integer, intent(in) :: node,cell(3),a
      real(dp), intent(in) :: low_packet(nv),high_packet(nv)
      logical, intent(in) :: root
      logical, intent(out) :: valid
      real(dp), optional, intent(inout) :: maximum_correction,deficit(nv)
      integer :: covered_by,lower(3),upper(3),first(3),last(3),index(3),face(3),r,ii,jj,kk,along
      real(dp) :: volume,total_volume,column_volume,used_volume,low_area,high_area,all_low_area,all_high_area
      real(dp) :: incoming(nv),outgoing(nv),previous(nv),next(nv),budget(nv),slope(nv),delta(nv),exchange_lower,exchange_upper
      integer :: bcast_status
      valid=.false.;volume=cell_volume(node,cell)
      if(volume<=0) return
      if(present(maximum_correction)) then
        if(rank==owners(node)) budget=candidate%fields(node)%values(:,cell(1),cell(2),cell(3))
        call MPI_Bcast(budget,nv,MPI_DOUBLE_PRECISION,owners(node),comm,bcast_status)
        if(bcast_status/=MPI_SUCCESS) return
        if(.not.all(ieee_is_finite(budget)).or.budget(1)<=0.or.any(budget(6:active_components)<0)) return
        slope=(low_packet-high_packet)/volume
        if(.not.all(ieee_is_finite(slope))) return
        exchange_lower=0;exchange_upper=1
        call reactive_eb_exchange_interval_3d(reshape(budget,[nv,1]),reshape(slope,[nv,1]), &
          exchange_lower,exchange_upper,valid, &
        passive_components=passive_components)
        if(.not.valid.or.exchange_lower>0) return
        maximum_correction=min(maximum_correction,exchange_upper)
        if(tree%nodes(parent)%level<minimum_level.and.tree%nodes(node)%level<=minimum_level) then
          if(rank==owners(node)) call cell_increment(node,cell,[0,0,0],0,delta)
          call MPI_Bcast(delta,nv,MPI_DOUBLE_PRECISION,owners(node),comm,bcast_status)
          if(bcast_status/=MPI_SUCCESS) return
          call add_mpi_eb_delta_cell_3d(budget,delta,valid)
          if(.not.valid) return
          exchange_lower=0;exchange_upper=1
          call reactive_eb_exchange_interval_3d(reshape(budget,[nv,1]),reshape(slope,[nv,1]), &
            exchange_lower,exchange_upper,valid,passive_components=passive_components)
          if(.not.valid.or.exchange_lower>0) return
          maximum_correction=min(maximum_correction,exchange_upper)
        endif
        if(.not.root) then
          face=cell;face(a)=face(a)-1
          if(open_area(node,a,face)<=0.and.any(low_packet/=0)) return
          if(open_area(node,a,cell)<=0.and.any(high_packet/=0)) return
        endif
      else
        if(rank==owners(node)) then
          candidate%fields(node)%values(:,cell(1),cell(2),cell(3))= &
            candidate%fields(node)%values(:,cell(1),cell(2),cell(3))+(low_packet-high_packet)/volume
          changed(node)=.true.
        endif
        if(.not.root) then
          face=cell;face(a)=face(a)-1
          call record_refined_packet(node,a,face,low_packet,valid)
          if(.not.valid) return
          call record_refined_packet(node,a,cell,high_packet,valid)
          if(.not.valid) return
        endif
      endif
      covered_by=covering_child(node,cell)
      valid=.true.
      if(covered_by==0) then
        ! Count physical leaves once; covered averages are not inventories.
        ! MPI broadcasts above keep this traversal identical on every rank.
        if(present(deficit)) then
          if(.not.present(maximum_correction)) then
            valid=.false.;return
          endif
          budget=budget+required_correction*slope
          if(.not.all(ieee_is_finite(budget))) then
            valid=.false.;return
          endif
          deficit(1)=deficit(1)+volume*max(0._dp,-budget(1))
          deficit(6:active_components)=deficit(6:active_components)+volume*max(0._dp,-budget(6:active_components))
          valid=all(ieee_is_finite(deficit))
        endif
        return
      endif
      call bounds(tree,covered_by,lower,upper)
      r=tree%nodes(covered_by)%patch%refinement_ratio
      first=(cell-lower)*r+1;last=first+r-1
      total_volume=0;all_low_area=0;all_high_area=0
      do kk=first(3),last(3)
      do jj=first(2),last(2)
      do ii=first(1),last(1)
        index=[ii,jj,kk];total_volume=total_volume+cell_volume(covered_by,index)
        if(index(a)/=first(a)) cycle
        face=index;face(a)=first(a)-1;all_low_area=all_low_area+open_area(covered_by,a,face)
        face(a)=last(a);all_high_area=all_high_area+open_area(covered_by,a,face)
      enddo
      enddo
      enddo
      valid=abs(total_volume-volume)<=64*epsilon(1._dp)*volume
      if(.not.valid) return
      if((all_low_area<=0.and.any(low_packet/=0)).or.(all_high_area<=0.and.any(high_packet/=0))) then
        valid=.false.;return
      endif
      do kk=first(3),last(3)
      do jj=first(2),last(2)
      do ii=first(1),last(1)
        index=[ii,jj,kk]
        if(index(a)/=first(a)) cycle
        face=index;face(a)=first(a)-1;low_area=open_area(covered_by,a,face)
        face(a)=last(a);high_area=open_area(covered_by,a,face)
        incoming=0;outgoing=0
        if(all_low_area>0) incoming=(low_area/all_low_area)*low_packet
        if(all_high_area>0) outgoing=(high_area/all_high_area)*high_packet
        column_volume=0
        do along=first(a),last(a)
          index(a)=along;column_volume=column_volume+cell_volume(covered_by,index)
        enddo
        if(column_volume<=0) then
          valid=all(incoming==0).and.all(outgoing==0)
          if(.not.valid) return
          cycle
        endif
        used_volume=0;previous=incoming
        do along=first(a),last(a)
          index(a)=along;used_volume=used_volume+cell_volume(covered_by,index)
          if(along==last(a)) then
            next=outgoing
          else
            next=(1-used_volume/column_volume)*incoming+(used_volume/column_volume)*outgoing
          endif
          if(cell_volume(covered_by,index)>0) then
            call propagate_cell_packet(covered_by,index,a,previous,next,.false.,valid,maximum_correction,deficit)
            if(.not.valid) return
          else if(any(previous/=next)) then
            valid=.false.;return
          endif
          previous=next
        enddo
      enddo
      enddo
      enddo
    end subroutine propagate_cell_packet
    ! Locate a straight path in the integer layout, including periodic seams.
    ! Coarse gaps and solid cells stop a route. Covered cells carry or donate
    ! through their physical descendants, with the same trial acceptance gate.
    ! Both copies of a patch-boundary face carry the same extensive packet.
    subroutine path_neighbor(node,cell,a,step,next_node,next_cell,face,next_face,valid)
      integer, intent(in) :: node,cell(3),a,step
      integer, intent(out) :: next_node,next_cell(3),face(3),next_face(3)
      logical, intent(out) :: valid
      integer(int64) :: global_cell(3)
      integer :: peer,depth
      valid=.false.;next_node=0;next_cell=0;face=cell;face(a)=face(a)+min(step,0)
      depth=tree%nodes(node)%level+1
      if(.not.active_path_cell(node,cell)) return
      global_cell=origin(:,node)+int(cell-1,int64);global_cell(a)=global_cell(a)+step
      if(periodic(a)) global_cell(a)=modulo(global_cell(a),domain(a,depth))
      if(any(global_cell<0).or.any(global_cell>=domain(:,depth))) return
      do peer=1,size(tree%nodes)
        if(tree%nodes(peer)%level/=depth-1) cycle
        if(any(global_cell<origin(:,peer)).or.any(global_cell>=origin(:,peer)+tree%nodes(peer)%shape)) cycle
        next_node=peer;next_cell=int(global_cell-origin(:,peer))+1;exit
      enddo
      if(next_node==0) return
      if(.not.active_path_cell(next_node,next_cell)) return
      next_face=next_cell;next_face(a)=next_face(a)+min(-step,0)
      valid=open_area(node,a,face)>0.and.open_area(next_node,a,next_face)>0
    end subroutine path_neighbor
    logical function active_path_cell(node,cell) result(active)
      integer, intent(in) :: node,cell(3)
      active=geometry(node)%volume_fraction(cell(1),cell(2),cell(3))>0
    end function active_path_cell
    subroutine add_path_packet(node,a,face,packet)
      integer, intent(in) :: node,a,face(3)
      real(dp), intent(in) :: packet(nv)
      real(dp) :: intensive(nv)
      intensive=packet/open_area(node,a,face)
      call set_face(working(node),a,face,face_value(working(node),a,face)+intensive)
      call set_face(effective(node),a,face,face_value(effective(node),a,face)+intensive)
    end subroutine add_path_packet
    real(dp) function open_area(p,a,index) result(value)
      integer, intent(in) :: p,a,index(3)
      associate(g=>geometry(p))
      select case(a)
      case(1);value=g%dy*g%dz*g%x_face_fraction(index(1),index(2),index(3))
      case(2);value=g%dx*g%dz*g%y_face_fraction(index(1),index(2),index(3))
      case(3);value=g%dx*g%dy*g%z_face_fraction(index(1),index(2),index(3))
      end select
      end associate
    end function open_area
    function face_value(f,a,index) result(value)
      type(eb_integrated_flux_3d), intent(in) :: f
      integer, intent(in) :: a,index(3)
      real(dp) :: value(nv)
      select case(a)
      case(1);value=f%x(:,index(1),index(2),index(3))
      case(2);value=f%y(:,index(1),index(2),index(3))
      case(3);value=f%z(:,index(1),index(2),index(3))
      end select
    end function face_value
    subroutine set_face(f,a,index,value)
      type(eb_integrated_flux_3d), intent(inout) :: f
      integer, intent(in) :: a,index(3)
      real(dp), intent(in) :: value(nv)
      select case(a)
      case(1);f%x(:,index(1),index(2),index(3))=value
      case(2);f%y(:,index(1),index(2),index(3))=value
      case(3);f%z(:,index(1),index(2),index(3))=value
      end select
    end subroutine set_face
    subroutine cell_increment(p,cell,excluded,excluded_axis,value)
      integer, intent(in) :: p,cell(3),excluded(3),excluded_axis
      real(dp), intent(out) :: value(nv)
      integer :: a,step,index(3)
      real(dp) :: volume,scale
      value=0
      if(present(rereflux_delta).and.excluded_axis==0) &
        value=routed_register%fields(p)%values(:,cell(1),cell(2),cell(3))
      if(immediate_funding) return
      volume=geometry(p)%volume_fraction(cell(1),cell(2),cell(3))*geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
      do a=1,3
      do step=-1,1,2
        index=cell;index(a)=index(a)+min(step,0)
        if(a==excluded_axis.and.all(index==excluded)) cycle
        scale=-real(step,dp)*open_area(p,a,index)/volume
        value=value+scale*(face_value(effective(p),a,index)-face_value(working(p),a,index))
      enddo
      enddo
    end subroutine cell_increment
  end subroutine repair_mpi_eb_tree_interface_components_3d

  subroutine synchronize_mpi_eb_tree_flux_core_3d(species,tree,owners,geometry,state,temperature,flux,comm,ok, &
      minimum_level,maximum_level,peers_only,boundary_faces,frame,reflux_increment,rereflux_delta,rereflux_first_level, &
      passive_components)
    integer, optional, intent(in) :: passive_components
    integer :: passives,passive_min,passive_max
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_level,maximum_level
    logical, optional, intent(in) :: peers_only
    character(len=*), optional, intent(in) :: boundary_faces(6)
    character(len=16) :: faces(6)
    integer :: codes(6),low_codes(6),high_codes(6),seam,side
    real(dp) :: period(3),length(3)
    type(amr_tree_state_3d), optional, intent(inout) :: reflux_increment
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    integer, optional, intent(in) :: rereflux_first_level
    type(amr_tree_state_3d) :: candidate,temps,candidate_increment,routed_delta
    type(eb_integrated_flux_3d), allocatable :: effective(:)
    type(eb_integrated_flux_3d) :: received,restricted
    type(eb_geometry_3d) :: section
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: offset(3),section_origin(3)
    logical, allocatable :: refined(:,:,:)
    integer :: nv,rank,status,node,peer,owner,receiver,axis,depth,parent,child,lo(3),hi(3),n(3),ratio,ios,i,j,k
    integer :: levels(4),low(4),high(4)
    integer :: supplied,minimum_supplied,maximum_supplied,axis_min,axis_max
    real(dp) :: frame_values(4),frame_min(4),frame_max(4)
    real(dp) :: primitive(reactive_nprim(size(species))),new_temperature,speed
    logical :: valid,ready
    ok=.false.;passives=0
    if(present(passive_components)) passives=passive_components
    call MPI_Allreduce(passives,passive_min,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(passives,passive_max,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.passive_min/=passive_max) return
    if(passives<0.or.passives>1) return
    if(passives>0.and..not.present(reflux_increment)) return
    nv=reactive_nvar(size(species))+passives
    supplied=merge(1,0,present(frame))+2*merge(1,0,present(reflux_increment))+4*merge(1,0,present(rereflux_delta))+ &
      8*merge(1,0,present(passive_components))
    call MPI_Allreduce(supplied,minimum_supplied,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(supplied,maximum_supplied,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.minimum_supplied/=maximum_supplied) return
    if(present(frame)) then
      call agree(frame%is_valid(),comm,ready)
      if(.not.ready) return
      call MPI_Allreduce(frame%axis,axis_min,1,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame%axis,axis_max,1,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS) return
      frame_values=[frame%angular_speed,frame%axis_point]
      call MPI_Allreduce(frame_values,frame_min,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame_values,frame_max,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.axis_min/=axis_max.or.any(frame_min/=frame_max)) return
    end if
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if (.not.ready) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    valid=all(faces=='outflow'.or.faces=='inflow'.or.faces=='periodic'.or. &
      faces=='slip_wall'.or.faces=='no_slip_wall'.or.faces=='isothermal_wall')
    valid=valid.and.periodic_eb_faces_match(geometry(1),faces)
    call agree(valid,comm,ready)
    if(.not.ready) return
    codes=merge(1,0,faces=='inflow')+merge(2,0,faces=='periodic')+ &
      merge(3,0,faces=='slip_wall')+merge(4,0,faces=='no_slip_wall')+merge(5,0,faces=='isothermal_wall')
    call MPI_Allreduce(codes,low_codes,6,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(codes,high_codes,6,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(low_codes/=high_codes)) return
    associate(g=>geometry(1))
    length=[g%x_upper-g%x_lower,g%y_upper-g%y_lower,g%z_upper-g%z_lower]
    end associate
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,ready)
    if (.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    levels=[0,tree%level_count()-1,0,0]
    if (present(peers_only)) levels(3)=merge(1,0,peers_only)
    if (present(minimum_level)) levels(1)=minimum_level
    if (present(maximum_level)) levels(2)=maximum_level
    levels(4)=levels(1)
    if(present(rereflux_first_level)) levels(4)=rereflux_first_level
    call MPI_Allreduce(levels,low,4,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(levels,high,4,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high).or.levels(1)<0.or.levels(2)<levels(1).or.levels(2)>=tree%level_count()) return
    if(levels(4)<0.or.levels(4)>levels(1)) return
    valid=size(flux)==size(owners)
    if (valid) then
      do node=1,size(owners)
        if (owners(node)==rank) then
          if (.not.valid_eb_integrated_flux_3d(flux(node),geometry(node),nv)) valid=.false.
        else
          if (allocated(flux(node)%x).or.allocated(flux(node)%y).or.allocated(flux(node)%z)) valid=.false.
        end if
      end do
    end if
    call agree(valid,comm,ready)
    if (.not.ready) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,comm,ready)
    if (.not.ready) return
    candidate=state;temps=temperature;effective=flux
    if(present(rereflux_delta)) then
      if(levels(3)/=0) return
      call route_mpi_eb_tree_rereflux_3d(tree,owners,rereflux_delta,nv,levels(4),levels(2),routed_delta,comm,ready)
      if(.not.ready) return
    endif
    if(present(reflux_increment).or.present(rereflux_delta)) then
      if(levels(3)/=0) return
      allocate(candidate_increment%fields(size(owners)),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(candidate_increment%fields(node)%values(nv,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        endif
        candidate_increment%fields(node)%values=0
      enddo
      call agree(valid,comm,ready)
      if(.not.ready) return
    endif
    do node=1,size(owners)
    do peer=node+1,size(owners)
      if (tree%nodes(node)%level<levels(1).or.tree%nodes(node)%level>levels(2)) cycle
      if (tree%nodes(peer)%level/=tree%nodes(node)%level) cycle
      owner=node;receiver=peer
      do axis=1,3
        if (origin(axis,node)==origin(axis,peer)) cycle
        if (origin(axis,node)>origin(axis,peer)) then
          owner=peer;receiver=node
        end if
        exit
      end do
      offset=origin(:,receiver)-origin(:,owner)
      if (any(offset < -int(tree%nodes(receiver)%shape,int64)).or. &
        any(offset>int(tree%nodes(owner)%shape,int64))) cycle
      call send_flux(effective(owner),owners(owner),owners(receiver),tree%nodes(owner)%shape,nv,received,comm,ready)
      if (.not.ready) return
      valid=.true.
      if (rank==owners(receiver)) then
        if (levels(3)==1) then
          call reconcile_reactive_eb_peer_flux_3d(species,geometry(owner),geometry(receiver),int(offset), &
            received,effective(receiver),candidate%fields(receiver)%values,temps%fields(receiver)%values(1,:,:,:),valid,frame=frame)
        else
          call reconcile_reactive_eb_peer_flux_3d(species,geometry(owner),geometry(receiver),int(offset), &
            received,effective(receiver),ok=valid,component_count=nv)
        end if
      end if
      call agree(valid,comm,ready)
      if (.not.ready) return
    end do
    end do
    do seam=1,3
      if(faces(2*seam-1)/='periodic') cycle
      period=0;period(seam)=length(seam)
      do owner=1,size(owners)
        depth=tree%nodes(owner)%level
        if(depth<levels(1).or.depth>levels(2)) cycle
        if(origin(seam,owner)+tree%nodes(owner)%shape(seam)/=domain(seam,depth+1)) cycle
        do receiver=1,size(owners)
          if(tree%nodes(receiver)%level/=depth.or.origin(seam,receiver)/=0) cycle
          offset=origin(:,receiver)-origin(:,owner)
          offset(seam)=offset(seam)+domain(seam,depth+1)
          if(.not.seam_overlap(offset,tree%nodes(owner)%shape,tree%nodes(receiver)%shape,seam)) cycle
          call send_flux(effective(owner),owners(owner),owners(receiver),tree%nodes(owner)%shape,nv,received,comm,ready)
          if(.not.ready) return
          valid=.true.
          if(rank==owners(receiver)) then
            if(levels(3)==1) then
              call reconcile_reactive_eb_peer_flux_3d(species,geometry(owner),geometry(receiver),int(offset), &
                received,effective(receiver),candidate%fields(receiver)%values, &
                temps%fields(receiver)%values(1,:,:,:),valid,periodic_length=period,frame=frame)
            else
              call reconcile_reactive_eb_peer_flux_3d(species,geometry(owner),geometry(receiver),int(offset), &
                received,effective(receiver),ok=valid,periodic_length=period,component_count=nv)
            end if
          end if
          call agree(valid,comm,ready)
          if(.not.ready) return
        end do
      end do
    end do
    ! Provisional level endpoints need peer corrections before child substeps,
    ! without restriction or averaging from descendants at another time.
    if (levels(3)==1) then
      call move_alloc(candidate%fields,state%fields)
      call move_alloc(temps%fields,temperature%fields)
      flux=effective;ok=.true.
      return
    end if
    do depth=levels(2),levels(1),-1
    do node=1,size(owners)
      if (tree%nodes(node)%level/=depth) cycle
      valid=.true.
      if (rank==owners(node)) then
        n=tree%nodes(node)%shape
        allocate(refined(n(1),n(2),n(3)),stat=ios)
        valid=ios==0
        if (valid) then
          refined=.false.
          do child=node+1,size(owners)
            if (tree%nodes(child)%parent/=node) cycle
            call bounds(tree,child,lo,hi)
            refined(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
          end do
          if(present(reflux_increment).or.present(rereflux_delta)) then
            call compute_reactive_eb_reflux_increment_3d(geometry(node),refined,flux(node),effective(node), &
              candidate_increment%fields(node)%values,valid,all_faces=.true.)
          else
          call reflux_reactive_eb_state_3d(species,geometry(node),refined,flux(node),effective(node), &
            candidate%fields(node)%values,temps%fields(node)%values(1,:,:,:),valid,all_faces=.true.,frame=frame)
          endif
          deallocate(refined)
        end if
      end if
      call agree(valid,comm,ready)
      if (.not.ready) return
      if (depth==levels(1)) cycle
      parent=tree%nodes(node)%parent;ratio=tree%nodes(node)%patch%refinement_ratio
      call bounds(tree,node,lo,hi)
      if (rank==owners(node)) call restrict_eb_integrated_flux_3d(geometry(node),ratio,effective(node),restricted,valid)
      call agree(valid,comm,ready)
      if (.not.ready) return
      call send_flux(restricted,owners(node),owners(parent),hi-lo+1,nv,received,comm,ready)
      if (.not.ready) return
      if (rank==owners(parent)) then
        effective(parent)%x(:,lo(1)-1:hi(1),lo(2):hi(2),lo(3):hi(3))=received%x
        effective(parent)%y(:,lo(1):hi(1),lo(2)-1:hi(2),lo(3):hi(3))=received%y
        effective(parent)%z(:,lo(1):hi(1),lo(2):hi(2),lo(3)-1:hi(3))=received%z
      end if
      call coarsen_eb_geometry_3d(geometry(node),ratio,section,valid)
      call agree(valid,comm,ready)
      if (.not.ready) return
      section_origin=origin(:,node)/int(ratio,int64)
      do peer=1,size(owners)
        if (peer==parent.or.tree%nodes(peer)%level/=depth-1) cycle
        offset=origin(:,peer)-section_origin
        if (any(offset < -int(tree%nodes(peer)%shape,int64)).or.any(offset>int(hi-lo+1,int64))) cycle
        call send_flux(restricted,owners(node),owners(peer),hi-lo+1,nv,received,comm,ready)
        if (.not.ready) return
        valid=.true.
        if (rank==owners(peer)) then
          call reconcile_reactive_eb_peer_flux_3d(species,section,geometry(peer),int(offset), &
            received,effective(peer),ok=valid,component_count=nv)
        end if
        call agree(valid,comm,ready)
        if (.not.ready) return
      end do
      do seam=1,3
        if(faces(2*seam-1)/='periodic') cycle
        period=0;period(seam)=length(seam)
        do side=0,1
          if(side==0.and.section_origin(seam)/=0) cycle
          if(side==1.and.section_origin(seam)+hi(seam)-lo(seam)+1/=domain(seam,depth)) cycle
          do peer=1,size(owners)
            if(tree%nodes(peer)%level/=depth-1) cycle
            if(side==0.and.origin(seam,peer)+tree%nodes(peer)%shape(seam)/=domain(seam,depth)) cycle
            if(side==1.and.origin(seam,peer)/=0) cycle
            offset=origin(:,peer)-section_origin
            offset(seam)=offset(seam)+int(2*side-1,int64)*domain(seam,depth)
            if(.not.seam_overlap(offset,hi-lo+1,tree%nodes(peer)%shape,seam)) cycle
            call send_flux(restricted,owners(node),owners(peer),hi-lo+1,nv,received,comm,ready)
            if(.not.ready) return
            valid=.true.
            if(rank==owners(peer)) call reconcile_reactive_eb_peer_flux_3d(species,section,geometry(peer),int(offset), &
              received,effective(peer),ok=valid,periodic_length=period,component_count=nv)
            call agree(valid,comm,ready)
            if(.not.ready) return
          end do
        end do
      end do
    end do
    end do
    if(present(rereflux_delta)) then
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          call add_mpi_eb_delta_cell_3d(candidate_increment%fields(node)%values(:,i,j,k), &
            routed_delta%fields(node)%values(:,i,j,k),ready)
          if(.not.ready) then
            valid=.false.;cycle
          endif
          if(present(reflux_increment)) cycle
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          if(all(candidate_increment%fields(node)%values(:,i,j,k)==0)) cycle
          call add_mpi_eb_delta_cell_3d(candidate%fields(node)%values(:,i,j,k), &
            candidate_increment%fields(node)%values(:,i,j,k),ready)
          if(.not.ready) then
            valid=.false.;cycle
          endif
          if(any(candidate%fields(node)%values(6:,i,j,k)<0)) then
            valid=.false.;cycle
          endif
          call recover_reactive_eb_cell_3d(species,geometry(node),[i,j,k],candidate%fields(node)%values(:,i,j,k), &
            temps%fields(node)%values(1,i,j,k),primitive,new_temperature,speed,ready,frame)
          if(.not.ready) then
            valid=.false.;cycle
          endif
          temps%fields(node)%values(1,i,j,k)=new_temperature
        enddo
        enddo
        enddo
      enddo
      call agree(valid,comm,ready)
      if(.not.ready) return
    endif
    if(present(reflux_increment)) then
      call move_alloc(candidate_increment%fields,reflux_increment%fields)
      flux=effective;ok=.true.
      return
    endif
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,temps,comm,ready,frame)
    if (.not.ready) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    flux=effective;ok=.true.
  end subroutine

  ! Owner-local one-component EB flux synchronization. Publish the field and
  ! effective fluxes only after every rank has accepted the full interval.
  subroutine synchronize_mpi_eb_tree_scalar_flux_3d(species,tree,owners,geometry,field,flux,comm,ok, &
      minimum_level,maximum_level,boundary_faces,rereflux_delta,rereflux_first_level)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: field
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_level,maximum_level
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    integer, optional, intent(in) :: rereflux_first_level
    type(amr_tree_state_3d) :: candidate,routed
    type(eb_integrated_flux_3d), allocatable :: effective(:)
    type(eb_integrated_flux_3d) :: received,restricted
    type(eb_geometry_3d) :: section
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: offset(3),section_origin(3)
    real(dp), allocatable :: delta(:,:,:,:)
    real(dp) :: period(3),length(3)
    logical, allocatable :: refined(:,:,:)
    character(len=16) :: faces(6)
    integer :: codes(6),low_codes(6),high_codes(6),levels(4),low(4),high(4)
    integer :: rank,status,node,peer,owner,receiver,axis,depth,parent,child,lo(3),hi(3),n(3),ratio,ios, &
      i,j,k,seam,side
    logical :: valid,ready
    ok=.false.
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,reactive_nvar(size(species)),comm,ready)
    if(.not.ready) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    valid=all(faces=='outflow'.or.faces=='inflow'.or.faces=='periodic'.or. &
      faces=='slip_wall'.or.faces=='no_slip_wall'.or.faces=='isothermal_wall')
    valid=valid.and.periodic_eb_faces_match(geometry(1),faces)
    call agree(valid,comm,ready)
    if(.not.ready) return
    codes=merge(1,0,faces=='inflow')+merge(2,0,faces=='periodic')+ &
      merge(3,0,faces=='slip_wall')+merge(4,0,faces=='no_slip_wall')+merge(5,0,faces=='isothermal_wall')
    call MPI_Allreduce(codes,low_codes,6,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(codes,high_codes,6,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low_codes/=high_codes)) return
    associate(g=>geometry(1))
    length=[g%x_upper-g%x_lower,g%y_upper-g%y_lower,g%z_upper-g%z_lower]
    end associate
    call validate_mpi_amr_tree_local_3d(tree,owners,field,1,comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    levels=[0,tree%level_count()-1,merge(1,0,present(rereflux_delta)),0]
    if(present(minimum_level)) levels(1)=minimum_level
    if(present(maximum_level)) levels(2)=maximum_level
    levels(4)=levels(1)
    if(present(rereflux_first_level)) levels(4)=rereflux_first_level
    call MPI_Allreduce(levels,low,4,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(levels,high,4,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    if(levels(1)<0.or.levels(2)<levels(1).or.levels(2)>=tree%level_count()) return
    if(levels(4)<0.or.levels(4)>levels(1)) return
    valid=size(flux)==size(owners)
    if(valid) then
      do node=1,size(owners)
        if(owners(node)==rank) then
          if(.not.valid_eb_integrated_flux_3d(flux(node),geometry(node),1)) valid=.false.
        else
          if(allocated(flux(node)%x).or.allocated(flux(node)%y).or.allocated(flux(node)%z)) valid=.false.
        endif
      enddo
    endif
    call agree(valid,comm,ready)
    if(.not.ready) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,comm,ready)
    if(.not.ready) return
    candidate=field;effective=flux
    if(present(rereflux_delta)) then
      call route_mpi_eb_tree_rereflux_3d(tree,owners,rereflux_delta,1,levels(4),levels(2),routed,comm,ready)
      if(.not.ready) return
    endif
    do node=1,size(owners)
    do peer=node+1,size(owners)
      if(tree%nodes(node)%level<levels(1).or.tree%nodes(node)%level>levels(2)) cycle
      if(tree%nodes(peer)%level/=tree%nodes(node)%level) cycle
      owner=node;receiver=peer
      do axis=1,3
        if(origin(axis,node)==origin(axis,peer)) cycle
        if(origin(axis,node)>origin(axis,peer)) then
          owner=peer;receiver=node
        endif
        exit
      enddo
      offset=origin(:,receiver)-origin(:,owner)
      if(any(offset < -int(tree%nodes(receiver)%shape,int64)).or. &
        any(offset>int(tree%nodes(owner)%shape,int64))) cycle
      call send_flux(effective(owner),owners(owner),owners(receiver),tree%nodes(owner)%shape,1,received,comm,ready)
      if(.not.ready) return
      valid=.true.
      if(rank==owners(receiver)) call reconcile_reactive_eb_peer_flux_3d(species,geometry(owner), &
        geometry(receiver),int(offset),received,effective(receiver),ok=valid,component_count=1)
      call agree(valid,comm,ready)
      if(.not.ready) return
    enddo
    enddo
    do seam=1,3
      if(faces(2*seam-1)/='periodic') cycle
      period=0;period(seam)=length(seam)
      do owner=1,size(owners)
        depth=tree%nodes(owner)%level
        if(depth<levels(1).or.depth>levels(2)) cycle
        if(origin(seam,owner)+tree%nodes(owner)%shape(seam)/=domain(seam,depth+1)) cycle
        do receiver=1,size(owners)
          if(tree%nodes(receiver)%level/=depth.or.origin(seam,receiver)/=0) cycle
          offset=origin(:,receiver)-origin(:,owner)
          offset(seam)=offset(seam)+domain(seam,depth+1)
          if(.not.seam_overlap(offset,tree%nodes(owner)%shape,tree%nodes(receiver)%shape,seam)) cycle
          call send_flux(effective(owner),owners(owner),owners(receiver),tree%nodes(owner)%shape,1,received,comm,ready)
          if(.not.ready) return
          valid=.true.
          if(rank==owners(receiver)) call reconcile_reactive_eb_peer_flux_3d(species,geometry(owner), &
            geometry(receiver),int(offset),received,effective(receiver),ok=valid, &
            periodic_length=period,component_count=1)
          call agree(valid,comm,ready)
          if(.not.ready) return
        enddo
      enddo
    enddo
    do depth=levels(2),levels(1),-1
    do node=1,size(owners)
      if(tree%nodes(node)%level/=depth) cycle
      valid=.true.
      if(rank==owners(node)) then
        n=tree%nodes(node)%shape
        allocate(refined(n(1),n(2),n(3)),stat=ios)
        valid=ios==0
        if(valid) then
          refined=.false.
          do child=node+1,size(owners)
            if(tree%nodes(child)%parent/=node) cycle
            call bounds(tree,child,lo,hi)
            refined(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
          enddo
          call compute_reactive_eb_reflux_increment_3d(geometry(node),refined,flux(node),effective(node), &
            delta,valid,all_faces=.true.)
          if(valid) valid=all(ieee_is_finite(delta))
          if(valid) then
            do k=1,n(3)
            do j=1,n(2)
            do i=1,n(1)
              if(refined(i,j,k).or.geometry(node)%volume_fraction(i,j,k)<=0) cycle
              call add_mpi_eb_delta_cell_3d(candidate%fields(node)%values(:,i,j,k),delta(:,i,j,k),ready)
              if(.not.ready) valid=.false.
            enddo
            enddo
            enddo
          endif
          deallocate(refined)
          if(allocated(delta)) deallocate(delta)
        endif
      endif
      call agree(valid,comm,ready)
      if(.not.ready) return
      if(depth==levels(1)) cycle
      parent=tree%nodes(node)%parent;ratio=tree%nodes(node)%patch%refinement_ratio
      call bounds(tree,node,lo,hi)
      if(rank==owners(node)) call restrict_eb_integrated_flux_3d(geometry(node),ratio,effective(node),restricted,valid)
      call agree(valid,comm,ready)
      if(.not.ready) return
      call send_flux(restricted,owners(node),owners(parent),hi-lo+1,1,received,comm,ready)
      if(.not.ready) return
      if(rank==owners(parent)) then
        effective(parent)%x(:,lo(1)-1:hi(1),lo(2):hi(2),lo(3):hi(3))=received%x
        effective(parent)%y(:,lo(1):hi(1),lo(2)-1:hi(2),lo(3):hi(3))=received%y
        effective(parent)%z(:,lo(1):hi(1),lo(2):hi(2),lo(3)-1:hi(3))=received%z
      endif
      call coarsen_eb_geometry_3d(geometry(node),ratio,section,valid)
      call agree(valid,comm,ready)
      if(.not.ready) return
      section_origin=origin(:,node)/int(ratio,int64)
      do peer=1,size(owners)
        if(peer==parent.or.tree%nodes(peer)%level/=depth-1) cycle
        offset=origin(:,peer)-section_origin
        if(any(offset < -int(tree%nodes(peer)%shape,int64)).or.any(offset>int(hi-lo+1,int64))) cycle
        call send_flux(restricted,owners(node),owners(peer),hi-lo+1,1,received,comm,ready)
        if(.not.ready) return
        valid=.true.
        if(rank==owners(peer)) call reconcile_reactive_eb_peer_flux_3d(species,section,geometry(peer), &
          int(offset),received,effective(peer),ok=valid,component_count=1)
        call agree(valid,comm,ready)
        if(.not.ready) return
      enddo
      do seam=1,3
        if(faces(2*seam-1)/='periodic') cycle
        period=0;period(seam)=length(seam)
        do side=0,1
          if(side==0.and.section_origin(seam)/=0) cycle
          if(side==1.and.section_origin(seam)+hi(seam)-lo(seam)+1/=domain(seam,depth)) cycle
          do peer=1,size(owners)
            if(tree%nodes(peer)%level/=depth-1) cycle
            if(side==0.and.origin(seam,peer)+tree%nodes(peer)%shape(seam)/=domain(seam,depth)) cycle
            if(side==1.and.origin(seam,peer)/=0) cycle
            offset=origin(:,peer)-section_origin
            offset(seam)=offset(seam)+int(2*side-1,int64)*domain(seam,depth)
            if(.not.seam_overlap(offset,hi-lo+1,tree%nodes(peer)%shape,seam)) cycle
            call send_flux(restricted,owners(node),owners(peer),hi-lo+1,1,received,comm,ready)
            if(.not.ready) return
            valid=.true.
            if(rank==owners(peer)) call reconcile_reactive_eb_peer_flux_3d(species,section,geometry(peer), &
              int(offset),received,effective(peer),ok=valid,periodic_length=period,component_count=1)
            call agree(valid,comm,ready)
            if(.not.ready) return
          enddo
        enddo
      enddo
    enddo
    enddo
    if(present(rereflux_delta)) then
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          call add_mpi_eb_delta_cell_3d(candidate%fields(node)%values(:,i,j,k), &
            routed%fields(node)%values(:,i,j,k),ready)
          if(.not.ready) valid=.false.
        enddo
        enddo
        enddo
      enddo
      call agree(valid,comm,ready)
      if(.not.ready) return
    endif
    call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,candidate,1,comm,ready)
    if(.not.ready) return
    call move_alloc(candidate%fields,field%fields)
    flux=effective;ok=.true.
  end subroutine synchronize_mpi_eb_tree_scalar_flux_3d

  subroutine synchronize_mpi_eb_tree_scalar_register_flux_3d(species,tree,owners,geometry,extended_geometry, &
      field,flux,registers,depth,comm,ok,boundary_faces)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),depth
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    type(amr_tree_state_3d), intent(inout) :: field
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    type(state_redist_register_3d), intent(inout) :: registers(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: correction,spill,candidate
    type(eb_integrated_flux_3d), allocatable :: working(:)
    integer :: rank,status,node,first_depth
    logical :: ready,has_spill
    ok=.false.
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call assemble_state_redist_mpi_interface_rereflux_3d(tree,owners,extended_geometry,registers,1, &
      depth,correction,comm,ready,boundary_faces,spill,has_spill)
    if(.not.ready) return
    first_depth=depth
    if(has_spill) then
      first_depth=depth-1
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(tree%nodes(node)%level==first_depth) correction%fields(node)%values=spill%fields(node)%values
      enddo
    endif
    candidate=field;working=flux
    call synchronize_mpi_eb_tree_scalar_flux_3d(species,tree,owners,geometry,candidate,working,comm,ready, &
      depth,depth+1,boundary_faces,correction,rereflux_first_level=first_depth)
    if(.not.ready) return
    field=candidate;flux=working
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      if(tree%nodes(node)%level==depth) then
        registers(node)%coarse_delta=0
      else if(tree%nodes(node)%level==depth+1) then
        registers(node)%fine_delta=0
      endif
    enddo
    ok=.true.
  end subroutine synchronize_mpi_eb_tree_scalar_register_flux_3d

  ! Sparse parent-cell requests carry only the coarse block needed by a
  ! child owner. Idle ranks participate without allocating gas/register fields.
  subroutine route_mpi_eb_tree_rereflux_3d(tree,owners,source,nv,first_level,last_level,increment,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nv,first_level,last_level
    type(amr_tree_state_3d), intent(in) :: source
    type(amr_tree_state_3d), intent(inout) :: increment
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    integer, allocatable :: requests(:,:)
    real(dp), allocatable :: values(:,:),block(:,:,:,:)
    logical, allocatable :: inside(:,:,:)
    integer :: rank,status,node,child,lo(3),hi(3),n(3),ratio,nq,i,j,k,q,ios
    logical :: valid,ready
    ok=.false.
    call validate_mpi_amr_tree_local_3d(tree,owners,source,nv,comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    valid=.true.
    do node=1,size(owners)
      if(rank/=owners(node)) cycle
      if(any(.not.ieee_is_finite(source%fields(node)%values))) valid=.false.
      if(tree%nodes(node)%level<first_level.or.tree%nodes(node)%level>last_level) then
        if(any(source%fields(node)%values/=0)) valid=.false.
      endif
    enddo
    call agree(valid,comm,ready)
    if(.not.ready) return
    candidate=source
    do child=2,size(owners)
      node=tree%nodes(child)%parent
      call bounds(tree,child,lo,hi);n=hi-lo+1;ratio=tree%nodes(child)%patch%refinement_ratio
      valid=product(int(n,int64))<=int(huge(1)/max(4,nv),int64)
      call agree(valid,comm,ready)
      if(.not.ready) return
      nq=0
      if(rank==owners(child)) nq=product(n)
      allocate(requests(4,nq),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      q=0
      if(rank==owners(child)) then
        do k=lo(3),hi(3)
        do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q=q+1;requests(:,q)=[node,i,j,k]
        enddo
        enddo
        enddo
      endif
      call fetch_mpi_amr_tree_cells_3d(tree,owners,candidate,nv,requests,values,comm,ready)
      if(.not.ready) return
      deallocate(requests)
      valid=.true.
      if(rank==owners(child)) then
        allocate(block(nv,n(1),n(2),n(3)),inside(n(1),n(2),n(3)),stat=ios)
        valid=ios==0
        if(valid) then
          block=reshape(values,[nv,n]);inside=.true.
          call add_eb_rereflux_to_fine_3d(block,inside,[ratio,ratio,ratio],[0,0,0],[0,0,0], &
            candidate%fields(child)%values,valid)
          deallocate(block,inside)
        endif
      endif
      call agree(valid,comm,ready)
      if(.not.ready) return
      if(rank==owners(node)) candidate%fields(node)%values(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=0
    enddo
    call move_alloc(candidate%fields,increment%fields);ok=.true.
  end subroutine route_mpi_eb_tree_rereflux_3d

  subroutine add_mpi_eb_delta_cell_3d(value,delta,ok)
    real(dp), intent(inout) :: value(:)
    real(dp), intent(in) :: delta(:)
    logical, intent(out) :: ok
    integer :: c
    ok=.false.
    if(size(value)/=size(delta)) return
    if(any(.not.ieee_is_finite(value)).or.any(.not.ieee_is_finite(delta))) return
    do c=1,size(value)
      if(delta(c)>0) then
        if(value(c)>huge(1._dp)-delta(c)) return
      else if(delta(c)<0) then
        if(value(c)<-huge(1._dp)-delta(c)) return
      endif
    enddo
    value=value+delta;ok=.true.
  end subroutine add_mpi_eb_delta_cell_3d

  pure logical function seam_overlap(displacement,owner_shape,receiver_shape,normal) result(overlap)
    integer(int64), intent(in) :: displacement(3)
    integer, intent(in) :: owner_shape(3),receiver_shape(3),normal
    integer :: tangent
    overlap=.true.
    do tangent=1,3
      if(tangent==normal) cycle
      if(displacement(tangent)>=owner_shape(tangent).or. &
        displacement(tangent)+receiver_shape(tangent)<=0) overlap=.false.
    end do
  end function seam_overlap

  ! A single temporary patch flux travels only from source owner to receiver.
  ! No other rank allocates its payload. Counts are checked before MPI calls.
  subroutine send_flux(source,sender,receiver,n,nv,received,comm,ok)
    type(eb_integrated_flux_3d), intent(in) :: source
    integer, intent(in) :: sender,receiver,n(3),nv
    type(eb_integrated_flux_3d), intent(inout) :: received
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: rank,status,counts(3),axis,d,shape4(4),ios
    integer(int64) :: elements
    logical :: valid,ready
    ok=.false.;valid=all(n>0).and.all(n<huge(1)).and.nv>0
    if (valid) then
      do axis=1,3
        shape4=[nv,n];shape4(axis+1)=shape4(axis+1)+1;elements=1_int64
        do d=1,4
          if (elements>int(huge(1),int64)/int(shape4(d),int64)) then
            valid=.false.;exit
          end if
          elements=elements*int(shape4(d),int64)
        end do
        if (.not.valid) exit
        counts(axis)=int(elements)
      end do
    end if
    call agree(valid,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    received=eb_integrated_flux_3d()
    if (rank==receiver) then
      allocate(received%x(nv,0:n(1),n(2),n(3)),received%y(nv,n(1),0:n(2),n(3)), &
        received%z(nv,n(1),n(2),0:n(3)),stat=ios)
      valid=ios==0
    end if
    call agree(valid,comm,ready)
    if (.not.ready) return
    if (sender==receiver) then
      if (rank==receiver) received=source
      ok=.true.;return
    end if
    do axis=1,3
      status=MPI_SUCCESS
      if (rank==sender) then
        select case(axis)
        case(1);call MPI_Send(source%x,counts(axis),MPI_DOUBLE_PRECISION,receiver,28510+axis,comm,status)
        case(2);call MPI_Send(source%y,counts(axis),MPI_DOUBLE_PRECISION,receiver,28510+axis,comm,status)
        case(3);call MPI_Send(source%z,counts(axis),MPI_DOUBLE_PRECISION,receiver,28510+axis,comm,status)
        end select
      else if (rank==receiver) then
        select case(axis)
        case(1)
          call MPI_Recv(received%x,counts(axis),MPI_DOUBLE_PRECISION,sender,28510+axis,comm,MPI_STATUS_IGNORE,status)
        case(2)
          call MPI_Recv(received%y,counts(axis),MPI_DOUBLE_PRECISION,sender,28510+axis,comm,MPI_STATUS_IGNORE,status)
        case(3)
          call MPI_Recv(received%z,counts(axis),MPI_DOUBLE_PRECISION,sender,28510+axis,comm,MPI_STATUS_IGNORE,status)
        end select
      end if
      call agree(status==MPI_SUCCESS,comm,ready)
      if (.not.ready) return
    end do
    ok=.true.
  end subroutine
end module mpi_eb_tree_flux_3d_mod
