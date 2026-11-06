module eb_held_funding_queue_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,valid_eb_integrated_flux_3d
  implicit none
  private
  public :: initialize_eb_held_funding_queue_3d,queue_eb_held_funding_receipt_3d
  public :: read_eb_held_funding_stock_3d,take_eb_held_funding_bucket_3d
  public :: take_eb_held_funding_target_3d
  public :: pack_eb_held_funding_faces_3d,physical_eb_held_funding_debt_3d
  public :: pack_eb_held_funding_increment_3d
  public :: initialize_eb_held_funding_context_3d,copy_eb_held_funding_context_3d
  public :: set_eb_held_funding_clock_faces_3d,close_eb_held_funding_clock_faces_3d
  public :: read_eb_held_funding_open_faces_3d,move_eb_held_funding_context_3d
  public :: valid_eb_held_funding_context_3d
  type, public :: eb_held_funding_family_3d
    integer :: clock=-1
    type(eb_integrated_flux_3d), allocatable :: gas(:),scalar(:)
  end type eb_held_funding_family_3d
  type, public :: eb_held_funding_context_3d
    type(amr_patch_tree_3d) :: topology
    type(eb_geometry_3d), allocatable :: geometry(:)
    type(amr_tree_state_3d), allocatable :: pending(:)
    type(eb_held_funding_family_3d), allocatable :: open(:)
  end type eb_held_funding_context_3d
  ! Each packet is a density impulse: gas conserved variables, then actual H.
  ! Buckets are native clock depths, not a second physical inventory. MPI
  ! callers supply their owner mask and agree before publishing a candidate.
contains
  ! Snapshot the entire pending receipt queue and every open ancestor
  ! family. A child provisional family lives separately until its accepted
  ! closure replaces it; accepted sibling families already belong to parent.
  subroutine initialize_eb_held_funding_context_3d(tree,geometry,nv,local,context,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv
    logical, intent(in) :: local(:)
    type(eb_held_funding_context_3d), intent(inout) :: context
    logical, intent(out) :: ok
    type(eb_held_funding_context_3d) :: candidate
    integer :: status
    logical :: good
    ok=.false.
    call initialize_eb_held_funding_queue_3d(tree,geometry,nv,local,candidate%pending,good)
    if(.not.good) return
    allocate(candidate%topology%nodes,source=tree%nodes,stat=status);if(status/=0) return
    allocate(candidate%geometry,source=geometry,stat=status);if(status/=0) return
    allocate(candidate%open(0:tree%level_count()-1),stat=status)
    if(status/=0) return
    call move_eb_held_funding_context_3d(candidate,context);ok=.true.
  end subroutine initialize_eb_held_funding_context_3d

  subroutine move_eb_held_funding_context_3d(source,destination)
    type(eb_held_funding_context_3d), intent(inout) :: source,destination
    call move_alloc(source%topology%nodes,destination%topology%nodes)
    call move_alloc(source%geometry,destination%geometry)
    call move_alloc(source%pending,destination%pending)
    call move_alloc(source%open,destination%open)
  end subroutine move_eb_held_funding_context_3d

  subroutine copy_eb_held_funding_context_3d(tree,geometry,nv,local,source,destination,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv
    logical, intent(in) :: local(:)
    type(eb_held_funding_context_3d), intent(in) :: source
    type(eb_held_funding_context_3d), intent(inout) :: destination
    logical, intent(out) :: ok
    type(eb_held_funding_context_3d) :: candidate
    integer :: b,status
    logical :: good
    ok=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    if(.not.valid_eb_held_funding_context_3d(tree,geometry,nv,local,source)) return
    call clone_queue(source%pending,candidate%pending,good);if(.not.good) return
    allocate(candidate%topology%nodes,source=tree%nodes,stat=status);if(status/=0) return
    allocate(candidate%geometry,source=geometry,stat=status);if(status/=0) return
    allocate(candidate%open(0:tree%level_count()-1),stat=status)
    if(status/=0) return
    do b=0,ubound(source%open,1)
      if(source%open(b)%clock<0) cycle
      call clone_faces(source%open(b)%gas,candidate%open(b)%gas,good);if(.not.good) return
      call clone_faces(source%open(b)%scalar,candidate%open(b)%scalar,good);if(.not.good) return
      candidate%open(b)%clock=b
    enddo
    call move_eb_held_funding_context_3d(candidate,destination);ok=.true.
  end subroutine copy_eb_held_funding_context_3d

  ! Set the provisional or final family of the deepest open native clock.
  ! It does not add that family to its parent or consume pending receipts.
  subroutine set_eb_held_funding_clock_faces_3d(tree,geometry,nv,local,clock,gas,scalar,context,result,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv,clock
    logical, intent(in) :: local(:)
    type(eb_integrated_flux_3d), allocatable, intent(in) :: gas(:),scalar(:)
    type(eb_held_funding_context_3d), intent(in) :: context
    type(eb_held_funding_context_3d), intent(inout) :: result
    logical, intent(out) :: ok
    type(eb_held_funding_context_3d) :: candidate
    logical :: good
    ok=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    if(clock<0.or.clock>=tree%level_count()) return
    if(.not.valid_face_pair(tree,geometry,nv,local,gas,scalar)) return
    call copy_eb_held_funding_context_3d(tree,geometry,nv,local,context,candidate,good)
    if(.not.good) return
    if(any(candidate%open(:clock-1)%clock<0).or.any(candidate%open(clock+1:)%clock>=0)) return
    call clone_faces(gas,candidate%open(clock)%gas,good);if(.not.good) return
    call clone_faces(scalar,candidate%open(clock)%scalar,good);if(.not.good) return
    candidate%open(clock)%clock=clock
    call move_eb_held_funding_context_3d(candidate,result);ok=.true.
  end subroutine set_eb_held_funding_clock_faces_3d

  ! The supplied closed child family replaces its provisional family. Only
  ! this accepted replacement is appended to the open parent. No subtraction
  ! of large nearly equal provisional families or double registration occurs.
  subroutine close_eb_held_funding_clock_faces_3d(tree,geometry,nv,local,clock,gas,scalar,context,result,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv,clock
    logical, intent(in) :: local(:)
    type(eb_integrated_flux_3d), allocatable, intent(in) :: gas(:),scalar(:)
    type(eb_held_funding_context_3d), intent(in) :: context
    type(eb_held_funding_context_3d), intent(inout) :: result
    logical, intent(out) :: ok
    type(eb_held_funding_context_3d) :: candidate
    integer :: p
    logical :: good
    ok=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    if(clock<0.or.clock>=tree%level_count()) return
    if(.not.valid_face_pair(tree,geometry,nv,local,gas,scalar)) return
    call copy_eb_held_funding_context_3d(tree,geometry,nv,local,context,candidate,good)
    if(.not.good) return
    if(candidate%open(clock)%clock/=clock.or.any(candidate%open(clock+1:)%clock>=0)) return
    ! Receipt payment and EOS acceptance belong to the physical transaction.
    ! Refuse to retire a clock while its own bucket remains unsettled.
    do p=1,size(tree%nodes)
      if(.not.local(p)) cycle
      if(any(candidate%pending(clock)%fields(p)%values/=0)) return
      if(clock==0) cycle
      call add_faces(candidate%open(clock-1)%gas(p),gas(p),good);if(.not.good) return
      call add_faces(candidate%open(clock-1)%scalar(p),scalar(p),good);if(.not.good) return
    enddo
    deallocate(candidate%open(clock)%gas,candidate%open(clock)%scalar)
    candidate%open(clock)%clock=-1
    call move_eb_held_funding_context_3d(candidate,result);ok=.true.
  end subroutine close_eb_held_funding_clock_faces_3d

  ! Read only ancestors of the upcoming native clock. Active current/deeper
  ! families are invalid here; the native operator owns its current face debt.
  subroutine read_eb_held_funding_open_faces_3d(tree,geometry,nv,local,clock,context,gas,scalar,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv,clock
    logical, intent(in) :: local(:)
    type(eb_held_funding_context_3d), intent(in) :: context
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: gas(:),scalar(:)
    logical, intent(out) :: ok
    type(eb_integrated_flux_3d), allocatable :: gas_candidate(:),scalar_candidate(:)
    integer :: b,p
    logical :: good
    ok=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    if(clock<0.or.clock>=tree%level_count()) return
    if(.not.valid_eb_held_funding_context_3d(tree,geometry,nv,local,context)) return
    if(any(context%open(:clock-1)%clock<0).or.any(context%open(clock:)%clock>=0)) return
    call zero_faces(tree,nv,local,gas_candidate,good);if(.not.good) return
    call zero_faces(tree,1,local,scalar_candidate,good);if(.not.good) return
    do b=0,clock-1
    do p=1,size(tree%nodes)
      if(.not.local(p)) cycle
      call add_faces(gas_candidate(p),context%open(b)%gas(p),good);if(.not.good) return
      call add_faces(scalar_candidate(p),context%open(b)%scalar(p),good);if(.not.good) return
    enddo
    enddo
    call move_alloc(gas_candidate,gas);call move_alloc(scalar_candidate,scalar);ok=.true.
  end subroutine read_eb_held_funding_open_faces_3d

  logical function valid_eb_held_funding_context_3d(tree,geometry,nv,local,context,clock) result(good)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv
    logical, intent(in) :: local(:)
    type(eb_held_funding_context_3d), intent(in) :: context
    integer, optional, intent(in) :: clock
    integer :: b,p
    logical :: retired
    good=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    if(.not.allocated(context%geometry)) return
    if(.not.valid_layout(context%topology,context%geometry,nv,local)) return
    if(size(context%topology%nodes)/=size(tree%nodes)) return
    do p=1,size(tree%nodes)
      associate(a=>tree%nodes(p),b=>context%topology%nodes(p))
      if(a%parent/=b%parent.or.a%level/=b%level.or.any(a%shape/=b%shape)) return
      if(any([a%patch%coarse_nx,a%patch%coarse_ny,a%patch%coarse_nz,a%patch%refinement_ratio, &
        a%patch%coarse_i_lower,a%patch%coarse_i_upper,a%patch%coarse_j_lower,a%patch%coarse_j_upper, &
        a%patch%coarse_k_lower,a%patch%coarse_k_upper]/=[b%patch%coarse_nx,b%patch%coarse_ny, &
        b%patch%coarse_nz,b%patch%refinement_ratio,b%patch%coarse_i_lower,b%patch%coarse_i_upper, &
        b%patch%coarse_j_lower,b%patch%coarse_j_upper,b%patch%coarse_k_lower,b%patch%coarse_k_upper])) return
      end associate
      if(.not.same_context_geometry(geometry(p),context%geometry(p))) return
    enddo
    if(.not.valid_queue(tree,geometry,nv,local,context%pending)) return
    if(.not.allocated(context%open)) return
    if(lbound(context%open,1)/=0.or.size(context%open)/=tree%level_count()) return
    retired=.false.
    do b=0,ubound(context%open,1)
      if(context%open(b)%clock==-1) then
        retired=.true.
        if(allocated(context%open(b)%gas).or.allocated(context%open(b)%scalar)) return
        cycle
      endif
      if(retired.or.context%open(b)%clock/=b) return
      if(.not.valid_face_pair(tree,geometry,nv,local,context%open(b)%gas,context%open(b)%scalar)) return
    enddo
    if(present(clock)) then
      if(clock<0.or.clock>=tree%level_count()) return
      if(any(context%open(:clock-1)%clock<0).or.any(context%open(clock:)%clock>=0)) return
    endif
    good=.true.
  end function valid_eb_held_funding_context_3d

  logical function same_context_geometry(a,b) result(good)
    type(eb_geometry_3d), intent(in) :: a,b
    good=.false.
    ! Both layouts have already passed is_valid. Replay binds actual input
    ! metrics exactly, rather than admitting a changed cached geometry.
    if(any([a%x_lower,a%x_upper,a%y_lower,a%y_upper,a%z_lower,a%z_upper,a%dx,a%dy,a%dz]/= &
      [b%x_lower,b%x_upper,b%y_lower,b%y_upper,b%z_lower,b%z_upper,b%dx,b%dy,b%dz])) return
    if(any(a%volume_fraction/=b%volume_fraction)) return
    if(any(a%cell_centroid_x/=b%cell_centroid_x)) return
    if(any(a%cell_centroid_y/=b%cell_centroid_y)) return
    if(any(a%cell_centroid_z/=b%cell_centroid_z)) return
    if(any(a%x_face_fraction/=b%x_face_fraction)) return
    if(any(a%y_face_fraction/=b%y_face_fraction)) return
    if(any(a%z_face_fraction/=b%z_face_fraction)) return
    if(any(a%x_face_centroid_y/=b%x_face_centroid_y)) return
    if(any(a%x_face_centroid_z/=b%x_face_centroid_z)) return
    if(any(a%y_face_centroid_x/=b%y_face_centroid_x)) return
    if(any(a%y_face_centroid_z/=b%y_face_centroid_z)) return
    if(any(a%z_face_centroid_x/=b%z_face_centroid_x)) return
    if(any(a%z_face_centroid_y/=b%z_face_centroid_y)) return
    if(any(a%boundary_area/=b%boundary_area)) return
    if(any(a%boundary_centroid_x/=b%boundary_centroid_x)) return
    if(any(a%boundary_centroid_y/=b%boundary_centroid_y)) return
    if(any(a%boundary_centroid_z/=b%boundary_centroid_z)) return
    if(any(a%boundary_normal_x/=b%boundary_normal_x)) return
    if(any(a%boundary_normal_y/=b%boundary_normal_y)) return
    if(any(a%boundary_normal_z/=b%boundary_normal_z)) return
    if(any(a%boundary_normal_integral_x/=b%boundary_normal_integral_x)) return
    if(any(a%boundary_normal_integral_y/=b%boundary_normal_integral_y)) return
    if(any(a%boundary_normal_integral_z/=b%boundary_normal_integral_z)) return
    if(any(a%cell_type/=b%cell_type)) return
    if(allocated(a%halo_volume_fraction).neqv.allocated(b%halo_volume_fraction)) return
    if(allocated(a%halo_volume_fraction)) then
      if(any(shape(a%halo_volume_fraction)/=shape(b%halo_volume_fraction)).or. &
        any(lbound(a%halo_volume_fraction)/=lbound(b%halo_volume_fraction))) return
      if(any(a%halo_volume_fraction/=b%halo_volume_fraction)) return
    endif
    good=.true.
  end function same_context_geometry

  logical function valid_face_pair(tree,geometry,nv,local,gas,scalar) result(good)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv
    logical, intent(in) :: local(:)
    type(eb_integrated_flux_3d), allocatable, intent(in) :: gas(:),scalar(:)
    integer :: p
    good=.false.
    if(.not.allocated(gas).or..not.allocated(scalar)) return
    if(size(gas)/=size(tree%nodes).or.size(scalar)/=size(tree%nodes)) return
    if(lbound(gas,1)/=1.or.lbound(scalar,1)/=1) return
    do p=1,size(tree%nodes)
      if(.not.local(p)) then
        if(allocated(gas(p)%x).or.allocated(gas(p)%y).or.allocated(gas(p)%z).or. &
          allocated(scalar(p)%x).or.allocated(scalar(p)%y).or.allocated(scalar(p)%z)) return
        cycle
      endif
      if(.not.valid_eb_integrated_flux_3d(gas(p),geometry(p),nv)) return
      if(.not.valid_eb_integrated_flux_3d(scalar(p),geometry(p),1)) return
    enddo
    good=.true.
  end function valid_face_pair

  subroutine clone_faces(source,destination,ok)
    type(eb_integrated_flux_3d), intent(in) :: source(:)
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: destination(:)
    logical, intent(out) :: ok
    type(eb_integrated_flux_3d), allocatable :: candidate(:)
    integer :: p,status
    ok=.false.
    allocate(candidate(size(source)),stat=status);if(status/=0) return
    do p=1,size(source)
      if(.not.allocated(source(p)%x)) cycle
      allocate(candidate(p)%x,source=source(p)%x,stat=status);if(status/=0) return
      allocate(candidate(p)%y,source=source(p)%y,stat=status);if(status/=0) return
      allocate(candidate(p)%z,source=source(p)%z,stat=status);if(status/=0) return
    enddo
    call move_alloc(candidate,destination);ok=.true.
  end subroutine clone_faces

  subroutine zero_faces(tree,components,local,faces,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: components
    logical, intent(in) :: local(:)
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: faces(:)
    logical, intent(out) :: ok
    type(eb_integrated_flux_3d), allocatable :: candidate(:)
    integer :: p,n(3),status
    ok=.false.
    allocate(candidate(size(tree%nodes)),stat=status);if(status/=0) return
    do p=1,size(tree%nodes)
      if(.not.local(p)) cycle
      n=tree%nodes(p)%shape
      allocate(candidate(p)%x(components,0:n(1),n(2),n(3)),candidate(p)%y(components,n(1),0:n(2),n(3)), &
        candidate(p)%z(components,n(1),n(2),0:n(3)),source=0._dp,stat=status)
      if(status/=0) return
    enddo
    call move_alloc(candidate,faces);ok=.true.
  end subroutine zero_faces

  subroutine add_faces(value,increment,ok)
    type(eb_integrated_flux_3d), intent(inout) :: value
    type(eb_integrated_flux_3d), intent(in) :: increment
    logical, intent(out) :: ok
    ! Inputs are validated full owner families. Private axis additions may
    ! fail partway; the enclosing context never publishes on that failure.
    ok=.false.
    if(any(increment%x/=0)) then
      call add_packet(value%x,increment%x,ok);if(.not.ok) return
    endif
    if(any(increment%y/=0)) then
      call add_packet(value%y,increment%y,ok);if(.not.ok) return
    endif
    if(any(increment%z/=0)) then
      call add_packet(value%z,increment%z,ok);if(.not.ok) return
    endif
    ok=.true.
  end subroutine add_faces

  subroutine pack_eb_held_funding_faces_3d(tree,geometry,nv,local,gas,scalar,packed,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv
    logical, intent(in) :: local(:)
    type(eb_integrated_flux_3d), intent(in) :: gas(:),scalar(:)
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: packed(:)
    logical, intent(out) :: ok
    type(eb_integrated_flux_3d), allocatable :: candidate(:)
    integer :: p,n(3),status
    ok=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    if(size(gas)/=size(tree%nodes).or.size(scalar)/=size(tree%nodes)) return
    allocate(candidate(size(tree%nodes)),stat=status)
    if(status/=0) return
    do p=1,size(tree%nodes)
      if(.not.local(p)) then
        if(allocated(gas(p)%x).or.allocated(gas(p)%y).or.allocated(gas(p)%z).or. &
          allocated(scalar(p)%x).or.allocated(scalar(p)%y).or.allocated(scalar(p)%z)) return
        cycle
      endif
      if(.not.valid_eb_integrated_flux_3d(gas(p),geometry(p),nv)) return
      if(.not.valid_eb_integrated_flux_3d(scalar(p),geometry(p),1)) return
      n=tree%nodes(p)%shape
      allocate(candidate(p)%x(nv+1,0:n(1),n(2),n(3)),candidate(p)%y(nv+1,n(1),0:n(2),n(3)), &
        candidate(p)%z(nv+1,n(1),n(2),0:n(3)),stat=status)
      if(status/=0) return
      candidate(p)%x(:nv,:,:,:)=gas(p)%x;candidate(p)%x(nv+1:nv+1,:,:,:)=scalar(p)%x
      candidate(p)%y(:nv,:,:,:)=gas(p)%y;candidate(p)%y(nv+1:nv+1,:,:,:)=scalar(p)%y
      candidate(p)%z(:nv,:,:,:)=gas(p)%z;candidate(p)%z(nv+1:nv+1,:,:,:)=scalar(p)%z
    enddo
    call move_alloc(candidate,packed);ok=.true.
  end subroutine pack_eb_held_funding_faces_3d

  ! Register re-reflux has virtual coarse entries needed by the existing
  ! covered-cell propagation. Pack them without the physical-stock mask;
  ! masking belongs after the joint conservative gas/H assembly.
  subroutine pack_eb_held_funding_increment_3d(tree,geometry,nv,local,gas,scalar,packed,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv
    logical, intent(in) :: local(:)
    type(amr_tree_state_3d), intent(in) :: gas,scalar
    type(amr_tree_state_3d), intent(inout) :: packed
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    integer :: p,n(3),status
    ok=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    if(.not.allocated(gas%fields).or..not.allocated(scalar%fields)) return
    if(size(gas%fields)/=size(tree%nodes).or.size(scalar%fields)/=size(tree%nodes)) return
    if(lbound(gas%fields,1)/=1.or.lbound(scalar%fields,1)/=1) return
    allocate(candidate%fields(size(tree%nodes)),stat=status)
    if(status/=0) return
    do p=1,size(tree%nodes)
      if(allocated(gas%fields(p)%values).neqv.local(p)) return
      if(allocated(scalar%fields(p)%values).neqv.local(p)) return
      if(.not.local(p)) cycle
      n=tree%nodes(p)%shape
      if(any(shape(gas%fields(p)%values)/=[nv,n]).or.any(shape(scalar%fields(p)%values)/=[1,n])) return
      if(any(lbound(gas%fields(p)%values)/=1).or.any(lbound(scalar%fields(p)%values)/=1)) return
      if(.not.all(ieee_is_finite(gas%fields(p)%values)).or. &
        .not.all(ieee_is_finite(scalar%fields(p)%values))) return
      allocate(candidate%fields(p)%values(nv+1,n(1),n(2),n(3)),stat=status)
      if(status/=0) return
      candidate%fields(p)%values(:nv,:,:,:)=gas%fields(p)%values
      candidate%fields(p)%values(nv+1:nv+1,:,:,:)=scalar%fields(p)%values
    enddo
    call move_alloc(candidate%fields,packed%fields);ok=.true.
  end subroutine pack_eb_held_funding_increment_3d

  subroutine physical_eb_held_funding_debt_3d(tree,geometry,nv,local,clock,source,debt,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv,clock
    logical, intent(in) :: local(:)
    type(amr_tree_state_3d), intent(in) :: source
    type(amr_tree_state_3d), intent(inout) :: debt
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    integer :: p,n(3),x,y,z
    logical :: good
    ok=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    if(clock<0.or.clock>=tree%level_count()) return
    ! Reflux assembly retains virtual coarse entries for propagation. They
    ! are not inventory available to a native physical donor transaction.
    if(.not.valid_packet(tree,geometry,nv,local,source,physical_only=.false.)) return
    call clone_packet(source,candidate,good);if(.not.good) return
    do p=1,size(tree%nodes)
      if(.not.local(p)) cycle
      if(tree%nodes(p)%level==clock) then
        candidate%fields(p)%values=0;cycle
      endif
      n=tree%nodes(p)%shape
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        if(geometry(p)%volume_fraction(x,y,z)>0.and..not.covered(tree,p,[x,y,z])) cycle
        candidate%fields(p)%values(:,x,y,z)=0
      enddo
      enddo
      enddo
    enddo
    call move_alloc(candidate%fields,debt%fields);ok=.true.
  end subroutine physical_eb_held_funding_debt_3d

  subroutine initialize_eb_held_funding_queue_3d(tree,geometry,nv,local,queue,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv
    logical, intent(in) :: local(:)
    type(amr_tree_state_3d), allocatable, intent(inout) :: queue(:)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), allocatable :: candidate(:)
    integer :: b,p,n(3),status
    ok=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    allocate(candidate(0:tree%level_count()-1),stat=status)
    if(status/=0) return
    do b=0,ubound(candidate,1)
      allocate(candidate(b)%fields(size(tree%nodes)),stat=status)
      if(status/=0) return
      do p=1,size(tree%nodes)
        if(.not.local(p)) cycle
        n=tree%nodes(p)%shape
        allocate(candidate(b)%fields(p)%values(nv+1,n(1),n(2),n(3)),source=0._dp,stat=status)
        if(status/=0) return
      enddo
    enddo
    call move_alloc(candidate,queue);ok=.true.
  end subroutine initialize_eb_held_funding_queue_3d

  subroutine queue_eb_held_funding_receipt_3d(tree,geometry,nv,local,clock,queue,receipt,result,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv,clock
    logical, intent(in) :: local(:)
    type(amr_tree_state_3d), allocatable, intent(in) :: queue(:)
    type(amr_tree_state_3d), intent(in) :: receipt
    type(amr_tree_state_3d), allocatable, intent(inout) :: result(:)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), allocatable :: candidate(:)
    integer :: p,b
    logical :: good
    ok=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    if(clock<0.or.clock>=tree%level_count()) return
    if(.not.valid_queue(tree,geometry,nv,local,queue)) return
    if(.not.valid_packet(tree,geometry,nv,local,receipt)) return
    do p=1,size(tree%nodes)
      if(.not.local(p)) cycle
      if(tree%nodes(p)%level/=clock) cycle
      if(any(receipt%fields(p)%values/=0)) return
    enddo
    call clone_queue(queue,candidate,good);if(.not.good) return
    do p=1,size(tree%nodes)
      if(.not.local(p)) cycle
      b=min(clock,tree%nodes(p)%level)
      call add_packet(candidate(b)%fields(p)%values,receipt%fields(p)%values,good)
      if(.not.good) return
    enddo
    call move_alloc(candidate,result);ok=.true.
  end subroutine queue_eb_held_funding_receipt_3d

  subroutine read_eb_held_funding_stock_3d(tree,geometry,nv,local,clock,queue,stock,ok,open_debt)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv,clock
    logical, intent(in) :: local(:)
    type(amr_tree_state_3d), allocatable, intent(in) :: queue(:)
    type(amr_tree_state_3d), intent(in), optional :: open_debt
    type(amr_tree_state_3d), intent(inout) :: stock
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    integer :: p,b,status,n(3)
    logical :: good
    ok=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    if(clock<0.or.clock>=tree%level_count()) return
    if(.not.valid_queue(tree,geometry,nv,local,queue)) return
    if(present(open_debt)) then
      if(.not.valid_packet(tree,geometry,nv,local,open_debt)) return
    endif
    allocate(candidate%fields(size(tree%nodes)),stat=status)
    if(status/=0) return
    do p=1,size(tree%nodes)
      if(.not.local(p)) cycle
      n=tree%nodes(p)%shape
      allocate(candidate%fields(p)%values(nv+1,n(1),n(2),n(3)),source=0._dp,stat=status)
      if(status/=0) return
      ! The active operator already owns its physical debt. Charge only
      ! inactive leaves, and never spend a virtual covered coarse cell.
      if(tree%nodes(p)%level==clock) cycle
      do b=0,ubound(queue,1)
        call add_packet(candidate%fields(p)%values,queue(b)%fields(p)%values,good)
        if(.not.good) return
      enddo
      if(present(open_debt)) then
        call add_packet(candidate%fields(p)%values,open_debt%fields(p)%values,good)
        if(.not.good) return
      endif
    enddo
    call move_alloc(candidate%fields,stock%fields);ok=.true.
  end subroutine read_eb_held_funding_stock_3d

  subroutine take_eb_held_funding_bucket_3d(tree,geometry,nv,local,clock,queue,receipt,remainder,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv,clock
    logical, intent(in) :: local(:)
    type(amr_tree_state_3d), allocatable, intent(in) :: queue(:)
    type(amr_tree_state_3d), intent(inout) :: receipt
    type(amr_tree_state_3d), allocatable, intent(inout) :: remainder(:)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), allocatable :: candidate(:)
    type(amr_tree_state_3d) :: packet
    integer :: p
    logical :: good
    ok=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    if(clock<0.or.clock>=tree%level_count()) return
    if(.not.valid_queue(tree,geometry,nv,local,queue)) return
    call clone_queue(queue,candidate,good);if(.not.good) return
    ! Move the selected bucket into a private receipt. The caller commits
    ! both outputs only after gas/H application and its own closure accept.
    call move_alloc(candidate(clock)%fields,packet%fields)
    call clone_packet(queue(clock),candidate(clock),good);if(.not.good) return
    do p=1,size(tree%nodes)
      if(local(p)) candidate(clock)%fields(p)%values=0
    enddo
    call move_alloc(packet%fields,receipt%fields)
    call move_alloc(candidate,remainder);ok=.true.
  end subroutine take_eb_held_funding_bucket_3d

  ! A parent may borrow from a physical descendant before that descendant's
  ! native entry. Pay its already queued debt before forming its own old state;
  ! the native inactive-stock protocol must not hide debt on the active level.
  ! Other targets remain in their original min(borrower,target) buckets.
  subroutine take_eb_held_funding_target_3d(tree,geometry,nv,local,clock,queue,receipt,remainder,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv,clock
    logical, intent(in) :: local(:)
    type(amr_tree_state_3d), allocatable, intent(in) :: queue(:)
    type(amr_tree_state_3d), intent(inout) :: receipt
    type(amr_tree_state_3d), allocatable, intent(inout) :: remainder(:)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), allocatable :: candidate(:)
    type(amr_tree_state_3d) :: packet
    integer :: p,b
    logical :: good
    ok=.false.
    if(.not.valid_layout(tree,geometry,nv,local)) return
    if(clock<0.or.clock>=tree%level_count()) return
    if(.not.valid_queue(tree,geometry,nv,local,queue)) return
    call clone_queue(queue,candidate,good);if(.not.good) return
    call clone_packet(queue(0),packet,good);if(.not.good) return
    do p=1,size(tree%nodes)
      if(.not.local(p)) cycle
      packet%fields(p)%values=0
      if(tree%nodes(p)%level/=clock) cycle
      do b=0,clock
        call add_packet(packet%fields(p)%values,queue(b)%fields(p)%values,good)
        if(.not.good) return
        candidate(b)%fields(p)%values=0
      enddo
    enddo
    call move_alloc(packet%fields,receipt%fields)
    call move_alloc(candidate,remainder);ok=.true.
  end subroutine take_eb_held_funding_target_3d

  logical function valid_layout(tree,geometry,nv,local) result(good)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv
    logical, intent(in) :: local(:)
    integer :: p
    good=.false.
    if(.not.tree%is_valid().or.nv<5.or.nv>=huge(1)) return
    if(size(geometry)/=size(tree%nodes).or.size(local)/=size(tree%nodes)) return
    do p=1,size(tree%nodes)
      if(.not.geometry(p)%is_valid()) return
      if(any([geometry(p)%nx,geometry(p)%ny,geometry(p)%nz]/=tree%nodes(p)%shape)) return
    enddo
    good=.true.
  end function valid_layout

  logical function valid_packet(tree,geometry,nv,local,packet,physical_only) result(good)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv
    logical, intent(in) :: local(:)
    type(amr_tree_state_3d), intent(in) :: packet
    logical, optional, intent(in) :: physical_only
    logical :: leaf_only
    integer :: p,n(3),x,y,z
    good=.false.;leaf_only=.true.
    if(present(physical_only)) leaf_only=physical_only
    if(.not.allocated(packet%fields)) return
    if(size(packet%fields)/=size(tree%nodes).or.lbound(packet%fields,1)/=1) return
    do p=1,size(tree%nodes)
      if(allocated(packet%fields(p)%values).neqv.local(p)) return
      if(.not.local(p)) cycle
      n=tree%nodes(p)%shape
      if(any(shape(packet%fields(p)%values)/=[nv+1,n]).or.any(lbound(packet%fields(p)%values)/=1)) return
      if(.not.all(ieee_is_finite(packet%fields(p)%values))) return
      if(.not.leaf_only) cycle
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        if(geometry(p)%volume_fraction(x,y,z)>0.and..not.covered(tree,p,[x,y,z])) cycle
        if(any(packet%fields(p)%values(:,x,y,z)/=0)) return
      enddo
      enddo
      enddo
    enddo
    good=.true.
  end function valid_packet

  logical function valid_queue(tree,geometry,nv,local,queue) result(good)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nv
    logical, intent(in) :: local(:)
    type(amr_tree_state_3d), allocatable, intent(in) :: queue(:)
    integer :: b,p
    good=.false.
    if(.not.allocated(queue)) return
    if(lbound(queue,1)/=0.or.size(queue)/=tree%level_count()) return
    do b=0,ubound(queue,1)
      if(.not.valid_packet(tree,geometry,nv,local,queue(b))) return
      do p=1,size(tree%nodes)
        if(.not.local(p).or.tree%nodes(p)%level>=b) cycle
        if(any(queue(b)%fields(p)%values/=0)) return
      enddo
    enddo
    good=.true.
  end function valid_queue

  subroutine clone_packet(source,destination,ok)
    type(amr_tree_state_3d), intent(in) :: source
    type(amr_tree_state_3d), intent(inout) :: destination
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    integer :: p,status
    ok=.false.
    allocate(candidate%fields(size(source%fields)),stat=status)
    if(status/=0) return
    do p=1,size(source%fields)
      if(.not.allocated(source%fields(p)%values)) cycle
      allocate(candidate%fields(p)%values,source=source%fields(p)%values,stat=status)
      if(status/=0) return
    enddo
    call move_alloc(candidate%fields,destination%fields);ok=.true.
  end subroutine clone_packet

  subroutine clone_queue(source,destination,ok)
    type(amr_tree_state_3d), allocatable, intent(in) :: source(:)
    type(amr_tree_state_3d), allocatable, intent(inout) :: destination(:)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), allocatable :: candidate(:)
    integer :: b,status
    logical :: good
    ok=.false.
    allocate(candidate(0:ubound(source,1)),stat=status)
    if(status/=0) return
    do b=0,ubound(source,1)
      call clone_packet(source(b),candidate(b),good);if(.not.good) return
    enddo
    call move_alloc(candidate,destination);ok=.true.
  end subroutine clone_queue

  subroutine add_packet(value,increment,ok)
    real(dp), intent(inout) :: value(:,:,:,:)
    real(dp), intent(in) :: increment(:,:,:,:)
    logical, intent(out) :: ok
    integer :: c,x,y,z,n(4)
    real(dp) :: a,b
    ok=.false.;n=shape(value)
    ! Check before adding, including builds trapping overflow. Do not create
    ! infinity first and turn a transactional refusal into a floating trap.
    do z=1,n(4)
    do y=1,n(3)
    do x=1,n(2)
    do c=1,n(1)
      a=value(c,x,y,z);b=increment(c,x,y,z)
      if(b>0) then
        if(a>huge(1._dp)-b) return
      else if(b<0) then
        if(a< -huge(1._dp)-b) return
      endif
    enddo
    enddo
    enddo
    enddo
    value=value+increment;ok=.true.
  end subroutine add_packet
end module eb_held_funding_queue_3d_mod
