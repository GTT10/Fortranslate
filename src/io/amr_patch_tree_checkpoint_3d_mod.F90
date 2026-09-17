module amr_patch_tree_checkpoint_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,valid_nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction,valid_elementary_reaction
  use gas_transport_mod, only: gas_transport_species,compatible_transport_database
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use amr_reactive_3d_checkpoint_mod, only: write_selected_checkpoint_context,checkpoint_temperature_is_consistent
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    average_down_amr_patch_tree_3d
  implicit none
  private
  type, public :: amr_tree_checkpoint_3d
    type(amr_patch_tree_3d) :: tree
    type(amr_tree_state_3d) :: state,temperature
    real(dp), allocatable :: initial(:),boundary(:)
    integer(int64), allocatable :: record_offsets(:)
    real(dp) :: time=0,maximum_reflux=0,minimum_theta=1,maximum_remap_error=0
    integer :: steps=0,regrids=0
  end type
  public :: checkpoint_record_offsets_3d
  public :: read_amr_tree_checkpoint_header_3d
  public :: build_amr_tree_checkpoint_context_3d,write_amr_tree_checkpoint_3d,read_amr_tree_checkpoint_3d
contains
  subroutine build_amr_tree_checkpoint_context_3d(config,species,reactions,transport,mole_fractions,identity,integrator, &
      schedule_identity,context,ok)
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mole_fractions(:)
    character(len=*), intent(in) :: identity,integrator,schedule_identity
    character(len=:), allocatable, intent(inout) :: context
    logical, intent(out) :: ok
    type(reactive_3d_config) :: canonical
    character(len=65536) :: line
    character(len=:), allocatable :: buffer
    integer :: unit,status,bytes,cursor,length,i
    ok=.false.
    if (size(species)<1.or.size(mole_fractions)/=size(species)) return
    if (.not.all(ieee_is_finite(mole_fractions))) return
    if (any(mole_fractions<0).or.abs(sum(mole_fractions)-1)>5e-10_dp) return
    do i=1,size(species)
      if (.not.valid_nasa7_species(species(i))) return
    end do
    do i=1,size(reactions)
      if (.not.valid_elementary_reaction(reactions(i),size(species))) return
    end do
    if (.not.compatible_transport_database(species,transport)) return
    if (len_trim(identity)==0.or.(integrator/='explicit'.and.integrator/='implicit')) return
    if (.not.valid_context(trim(schedule_identity)//achar(10))) return
    ! Capture all configuration fields except continuation limits/output path.
    canonical=config;canonical%final_time=0;canonical%maximum_steps=0;canonical%output_file=''
    open(newunit=unit,status='scratch',access='stream',form='formatted',delim='quote',iostat=status)
    if (status/=0) return
    write(unit,*,iostat=status) canonical
    if (status/=0) go to 900
    do i=1,size(species)
      write(unit,*,iostat=status) species(i)
      if (status/=0) go to 900
      write(unit,*,iostat=status) transport(i)
      if (status/=0) go to 900
    end do
    call write_selected_checkpoint_context(unit,reactions,identity,integrator,mole_fractions, &
      config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance,status)
    if (status/=0) go to 900
    write(unit,'(a)',iostat=status) schedule_identity
    if (status/=0) go to 900
    inquire(unit=unit,size=bytes,iostat=status)
    if (status/=0.or.bytes<1) go to 900
    allocate(character(len=bytes)::buffer);cursor=0
    rewind(unit,iostat=status)
    if (status/=0) go to 900
    do
      read(unit,'(a)',iostat=status) line
      if (status<0) exit
      if (status/=0.or.len_trim(line)==len(line)) go to 900
      length=len_trim(line)
      if (length+1>bytes-cursor) go to 900
      buffer(cursor+1:cursor+length+1)=trim(line)//achar(10);cursor=cursor+length+1
    end do
    context=buffer(:cursor);ok=.true.
900 close(unit)
  end subroutine

  subroutine write_amr_tree_checkpoint_3d(path,species,context,snapshot,ok,indexed)
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_checkpoint_3d), intent(in) :: snapshot
    logical, intent(out) :: ok
    logical, optional, intent(in) :: indexed
    logical :: with_index,valid
    integer(int64), allocatable :: offsets(:)
    integer :: unit,status,node,i,j,k,n(3),start,finish
    ok=.false.
    if (.not.valid_snapshot(snapshot,species)) return
    if (.not.valid_context(context)) return
    with_index=.false.
    if (present(indexed)) with_index=indexed
    if (with_index) then
      call checkpoint_record_offsets_3d(snapshot%tree,size(snapshot%initial),offsets,valid)
      if (.not.valid) return
    end if
    open(newunit=unit,file=path,status='replace',action='write',iostat=status)
    if (status/=0) return
    if (with_index) then
      write(unit,'(a)',iostat=status) 'PELEF_AMR_TREE_3D_CHECKPOINT_V2'
    else
      write(unit,'(a)',iostat=status) 'PELEF_AMR_TREE_3D_CHECKPOINT_V1'
    end if
    if (status/=0) go to 900
    write(unit,'(i0)',iostat=status) len(context)
    if (status/=0) go to 900
    start=1
    do while(start<=len(context))
      finish=start+index(context(start:),achar(10))-2
      write(unit,'(a)',iostat=status) context(start:finish)
      if (status/=0) go to 900
      start=finish+2
    end do
    write(unit,'(*(i0,1x))',iostat=status) size(snapshot%tree%nodes),size(snapshot%initial),snapshot%tree%nodes(1)%shape
    if (status/=0) go to 900
    write(unit,'(4(es27.18e3,1x),2(i0,1x))',iostat=status) snapshot%time,snapshot%maximum_reflux, &
      snapshot%minimum_theta,snapshot%maximum_remap_error,snapshot%steps,snapshot%regrids
    if (status/=0) go to 900
    write(unit,'(*(es27.18e3,1x))',iostat=status) snapshot%initial,snapshot%boundary
    if (status/=0) go to 900
    do node=2,size(snapshot%tree%nodes)
      associate(p=>snapshot%tree%nodes(node)%patch)
        write(unit,'(*(i0,1x))',iostat=status) snapshot%tree%nodes(node)%parent,p%coarse_i_lower,p%coarse_i_upper, &
          p%coarse_j_lower,p%coarse_j_upper,p%coarse_k_lower,p%coarse_k_upper,p%refinement_ratio
      end associate
      if (status/=0) go to 900
    end do
    if (with_index) then
      do node=1,size(offsets)
        write(unit,'(i0)',iostat=status) offsets(node)
        if (status/=0) go to 900
      end do
    end if
    do node=1,size(snapshot%tree%nodes)
      n=snapshot%tree%nodes(node)%shape
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            write(unit,'(*(es27.18e3,1x))',iostat=status) snapshot%state%fields(node)%values(:,i,j,k), &
              snapshot%temperature%fields(node)%values(1,i,j,k)
            if (status/=0) go to 900
          end do
        end do
      end do
    end do
    write(unit,'(a)',iostat=status) 'END_AMR_TREE_CHECKPOINT'
900 close(unit,iostat=node)
    ok=status==0.and.node==0
  end subroutine

  subroutine read_amr_tree_checkpoint_3d(path,species,context,root_shape,snapshot,ok)
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: root_shape(3)
    type(amr_tree_checkpoint_3d), intent(inout) :: snapshot
    logical, intent(out) :: ok
    type(amr_tree_checkpoint_3d) :: candidate
    character(len=65536) :: line
    integer :: unit,status,count,nvar,n(3),node,i,j,k
    integer(int64) :: payload_start,position
    logical :: valid
    ok=.false.
    if (.not.valid_context(context)) return
    open(newunit=unit,file=path,status='old',action='read',access='stream',form='formatted',iostat=status)
    if (status/=0) return
    call read_amr_tree_checkpoint_header_3d(unit,species,context,root_shape,candidate,valid)
    if (.not.valid) go to 900
    count=size(candidate%tree%nodes);nvar=size(candidate%initial)
    inquire(unit=unit,pos=payload_start,iostat=status)
    if (status/=0) go to 900
    allocate(candidate%state%fields(count),candidate%temperature%fields(count))
    do node=1,count
      n=candidate%tree%nodes(node)%shape
      allocate(candidate%state%fields(node)%values(nvar,n(1),n(2),n(3)), &
        candidate%temperature%fields(node)%values(1,n(1),n(2),n(3)),stat=status)
      if (status/=0) go to 900
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            read(unit,*,iostat=status) candidate%state%fields(node)%values(:,i,j,k), &
              candidate%temperature%fields(node)%values(1,i,j,k)
            if (status/=0) go to 900
          end do
        end do
      end do
      if (allocated(candidate%record_offsets)) then
        inquire(unit=unit,pos=position,iostat=status)
        if (status/=0) go to 900
        if (position-payload_start/=candidate%record_offsets(node+1)) go to 900
      end if
    end do
    if (allocated(candidate%record_offsets)) then
      inquire(unit=unit,pos=position,iostat=status)
      if (status/=0) go to 900
      if (position-payload_start/=candidate%record_offsets(count+1)) go to 900
    end if
    read(unit,'(a)',iostat=status) line
    if (status/=0) go to 900
    if (trim(line)/='END_AMR_TREE_CHECKPOINT') go to 900
    read(unit,'(a)',iostat=status) line
    if (status>=0) go to 900
    if (.not.valid_snapshot(candidate,species)) go to 900
    snapshot=candidate;ok=.true.
900 close(unit)
  end subroutine

  ! Read only V1/V2 metadata from an already-open formatted unit.
  subroutine read_amr_tree_checkpoint_header_3d(unit,species,context,root_shape,snapshot,ok)
    integer, intent(in) :: unit,root_shape(3)
    type(nasa7_species), intent(in) :: species(:)
    character(len=*), intent(in) :: context
    type(amr_tree_checkpoint_3d), intent(inout) :: snapshot
    logical, intent(out) :: ok
    type(amr_tree_checkpoint_3d) :: candidate
    type(amr_patch_3d), allocatable :: patches(:)
    integer, allocatable :: parents(:),shapes(:,:)
    character(len=65536) :: line
    integer :: status,count,nvar,n(3),node,start,finish,context_length,lo(3),hi(3),ratio,parent
    integer(int64), allocatable :: expected_offsets(:)
    logical :: valid,indexed
    ok=.false.
    if (.not.valid_context(context)) return
    read(unit,'(a)',iostat=status) line
    if (status/=0) return
    indexed=trim(line)=='PELEF_AMR_TREE_3D_CHECKPOINT_V2'
    if (.not.indexed.and.trim(line)/='PELEF_AMR_TREE_3D_CHECKPOINT_V1') return
    read(unit,*,iostat=status) context_length
    if (status/=0) return
    if (context_length/=len(context)) return
    start=1
    do while(start<=len(context))
      finish=start+index(context(start:),achar(10))-2
      read(unit,'(a)',iostat=status) line
      if (status/=0) return
      if (trim(line)/=context(start:finish)) return
      start=finish+2
    end do
    read(unit,*,iostat=status) count,nvar,n
    if (status/=0) return
    if (count<1.or.nvar/=reactive_nvar(size(species)).or.any(n/=root_shape)) return
    allocate(patches(count-1),parents(count-1),shapes(3,count),stat=status)
    if (status/=0) return
    shapes(:,1)=n
    read(unit,*,iostat=status) candidate%time,candidate%maximum_reflux,candidate%minimum_theta, &
      candidate%maximum_remap_error,candidate%steps,candidate%regrids
    if (status/=0) return
    allocate(candidate%initial(nvar),candidate%boundary(nvar))
    read(unit,*,iostat=status) candidate%initial,candidate%boundary
    if (status/=0) return
    do node=2,count
      read(unit,*,iostat=status) parent,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3),ratio
      if (status/=0) return
      if (parent<1.or.parent>=node) return
      n=shapes(:,parent);parents(node-1)=parent
      call initialize_amr_patch_3d(n(1),n(2),n(3),lo(1),hi(1),lo(2),hi(2),lo(3),hi(3),ratio,patches(node-1),valid)
      if (.not.valid) return
      if (any(int(hi-lo+1,int64)*int(ratio,int64)>int(huge(1),int64))) return
      shapes(:,node)=int(int(hi-lo+1,int64)*int(ratio,int64))
    end do
    call initialize_amr_patch_tree_3d(root_shape,parents,patches,candidate%tree,valid)
    if (.not.valid) return
    if (indexed) then
      call checkpoint_record_offsets_3d(candidate%tree,nvar,expected_offsets,valid)
      if (.not.valid) return
      allocate(candidate%record_offsets(size(expected_offsets)),stat=status)
      if (status/=0) return
      do node=1,size(expected_offsets)
        read(unit,*,iostat=status) candidate%record_offsets(node)
        if (status/=0) return
      end do
      if (any(candidate%record_offsets/=expected_offsets)) return
    end if
    snapshot=candidate;ok=.true.
  end subroutine

  subroutine checkpoint_record_offsets_3d(tree,nvar,offsets,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: nvar
    integer(int64), allocatable, intent(out) :: offsets(:)
    logical, intent(out) :: ok
    integer(int64) :: cells,width
    integer :: node,axis,status
    ok=.false.
    if (.not.tree%is_valid().or.nvar<1) return
    width=28_int64*(int(nvar,int64)+1_int64)
    allocate(offsets(size(tree%nodes)+1),stat=status)
    if (status/=0) return
    offsets(1)=0
    do node=1,size(tree%nodes)
      cells=1
      do axis=1,3
        if (cells>huge(cells)/int(tree%nodes(node)%shape(axis),int64)) return
        cells=cells*int(tree%nodes(node)%shape(axis),int64)
      end do
      if (cells>(huge(cells)-offsets(node))/width) return
      offsets(node+1)=offsets(node)+cells*width
    end do
    ok=.true.
  end subroutine

  logical function valid_context(context) result(valid)
    character(len=*), intent(in) :: context
    integer :: start,finish
    valid=.false.
    if (len(context)<1) return
    if (context(len(context):)/=achar(10)) return
    start=1
    do while(start<=len(context))
      finish=start+index(context(start:),achar(10))-2
      if (finish-start+1>=65536) return
      start=finish+2
    end do
    valid=.true.
  end function

  logical function valid_snapshot(s,species) result(valid)
    type(amr_tree_checkpoint_3d), intent(in) :: s
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d) :: synchronized
    integer :: node,nvar
    logical :: consistent
    valid=.false.;nvar=reactive_nvar(size(species))
    if (.not.s%state%is_valid(s%tree).or..not.s%temperature%is_valid(s%tree)) return
    if (size(s%state%fields(1)%values,1)/=nvar.or.size(s%temperature%fields(1)%values,1)/=1) return
    if (.not.allocated(s%initial).or..not.allocated(s%boundary)) return
    if (size(s%initial)/=nvar.or.size(s%boundary)/=nvar) return
    if (.not.all(ieee_is_finite(s%initial)).or..not.all(ieee_is_finite(s%boundary))) return
    if (.not.all(ieee_is_finite([s%time,s%maximum_reflux,s%minimum_theta,s%maximum_remap_error]))) return
    if (min(s%time,s%maximum_reflux,s%minimum_theta,s%maximum_remap_error)<0.or.s%minimum_theta>1) return
    if (min(s%steps,s%regrids)<0.or.s%regrids>s%steps) return
    if ((s%steps==0.and.s%time/=0).or.(s%steps>0.and.s%time<=0)) return
    synchronized=s%state
    call average_down_amr_patch_tree_3d(s%tree,synchronized,consistent)
    if (.not.consistent) return
    do node=1,size(s%tree%nodes)
      if (any(abs(synchronized%fields(node)%values-s%state%fields(node)%values)> &
          64*epsilon(1.0_dp)*max(1.0_dp,abs(s%state%fields(node)%values)))) return
      call checkpoint_temperature_is_consistent(species,s%state%fields(node)%values, &
        s%temperature%fields(node)%values(1,:,:,:),consistent)
      if (.not.consistent) return
    end do
    valid=.true.
  end function
end module amr_patch_tree_checkpoint_3d_mod
