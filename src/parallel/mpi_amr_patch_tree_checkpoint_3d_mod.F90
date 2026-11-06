module mpi_amr_patch_tree_checkpoint_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: iso_c_binding, only: c_char,c_int,c_null_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use reactive_mol_3d_mod, only: convert_regular_energy_3d
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_checkpoint_3d_mod, only: &
    amr_tree_checkpoint_3d,read_amr_tree_checkpoint_header_3d,checkpoint_record_offsets_3d, &
    valid_amr_tree_timestep_history_3d,checkpoint_tree_held_matches
  use amr_reactive_3d_checkpoint_mod, only: checkpoint_temperature_is_consistent
  use mpi_amr_patch_tree_chemistry_3d_mod, only: clone_mpi_amr_tree_inputs_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_distribution_3d_mod, only: partition_amr_tree_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_metadata_3d,validate_mpi_amr_tree_local_3d
  implicit none
  private
  public :: write_mpi_amr_tree_checkpoint_3d,read_mpi_amr_tree_checkpoint_3d
  ! Publish only fully closed checkpoints; retain the previous file on failure.
  interface
    function checkpoint_rename(old,new) bind(C,name='rename') result(status)
      import c_int,c_char
      character(kind=c_char), intent(in) :: old(*),new(*)
      integer(c_int) :: status
    end function
    function checkpoint_unlink(path) bind(C,name='unlink') result(status)
      import c_int,c_char
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int) :: status
    end function
  end interface
contains
  subroutine read_mpi_amr_tree_checkpoint_3d(path,species,context,root_shape,partitioning,exponent, &
      snapshot,owners,comm,ok,allow_internal_energy,frame,domain_lower,root_spacing,mol_energy_reset)
    logical, optional, intent(in) :: mol_energy_reset
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3),root_spacing(3)
    character(len=:), allocatable :: rotation_context
    character(len=*), intent(in) :: path,context,partitioning
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: root_shape(3),exponent
    type(amr_tree_checkpoint_3d), intent(inout) :: snapshot
    integer, allocatable, intent(inout) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    logical, optional, intent(in) :: allow_internal_energy
    type(amr_tree_checkpoint_3d) :: candidate
    integer, allocatable :: assigned(:)
    integer(MPI_OFFSET_KIND), allocatable :: offsets(:),reaction_offsets(:),heat_offsets(:)
    integer(MPI_OFFSET_KIND) :: payload_start,position,footer_position
    real(dp) :: scratch(reactive_nvar(size(species))+2)
    character(len=65536) :: line
    character(len=:), allocatable :: header
    integer :: rank,nranks,status,ios,unit,node,nvar,n(3),i,j,k,low,high,payload_size
    logical :: valid,opened
    ok=.false.;nvar=reactive_nvar(size(species))
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,status)
    if (status/=MPI_SUCCESS) return
    call checkpoint_rotation_context(comm,context,frame,domain_lower,root_spacing,rotation_context,valid)
    if(.not.valid) return
    call same_checkpoint_string(comm,trim(path),valid)
    if (.not.valid) return
    call same_checkpoint_string(comm,context,valid)
    if (.not.valid) return
    call same_checkpoint_string(comm,partitioning,valid)
    if (.not.valid) return
    if (partitioning/='weighted'.and.partitioning/='round_robin') return
    call MPI_Allreduce(exponent,low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(exponent,high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.low/=high.or.low<0.or.low>2) return
    open(newunit=unit,file=path,status='old',action='read',access='stream',form='formatted',iostat=ios)
    opened=ios==0;valid=opened;call agree_checkpoint(comm,valid)
    if (.not.valid) then
      if (opened) close(unit)
      return
    end if
    call read_amr_tree_checkpoint_header_3d(unit,species,rotation_context,root_shape,candidate,valid, &
      allow_sdc=.true.,allow_internal_energy=allow_internal_energy)
    call agree_checkpoint(comm,valid)
    if (.not.valid) go to 900
    payload_size=nvar+1+merge(1,0,allocated(candidate%internal_energy%fields))
    write(line,'(4l1)') allocated(candidate%reaction_source%fields),allocated(candidate%heat_release%fields), &
      allocated(candidate%internal_energy%fields),allocated(candidate%pressure_work_change)
    call same_checkpoint_string(comm,trim(line),valid)
    if(.not.valid) go to 900
    if (partitioning=='weighted') then
      call partition_amr_tree_3d(candidate%tree,nranks,exponent,assigned,valid)
    else
      allocate(assigned(size(candidate%tree%nodes)),stat=ios)
      valid=ios==0
      if (valid) assigned=[(mod(node-1,nranks),node=1,size(assigned))]
    end if
    call agree_checkpoint(comm,valid)
    if (.not.valid) go to 900
    call validate_mpi_amr_tree_metadata_3d(candidate%tree,assigned,nvar,comm,valid)
    if (.not.valid) go to 900
    allocate(offsets(size(assigned)),reaction_offsets(size(assigned)),heat_offsets(size(assigned)),stat=ios)
    valid=ios==0;call agree_checkpoint(comm,valid)
    if (.not.valid) go to 900
    offsets=0
    ! V1 has no on-disk index. Stream through one cell record at a time on rank
    ! zero to index node positions, preserving the serial list-directed syntax.
    ! Only offsets are broadcast; no complete state is retained on rank zero.
    if (rank==0) then
      if (allocated(candidate%record_offsets)) then
        inquire(unit=unit,pos=payload_start,iostat=ios)
        valid=ios==0
        if (valid) valid=candidate%record_offsets(size(assigned)+1)<=huge(payload_start)-payload_start
        if (valid) then
          offsets=payload_start+candidate%record_offsets(:size(assigned))
          footer_position=payload_start+candidate%record_offsets(size(assigned)+1)
        end if
      else
        scan_nodes: do node=1,size(assigned)
          inquire(unit=unit,pos=offsets(node),iostat=ios)
          if (ios/=0) then
            valid=.false.;exit
          end if
          n=candidate%tree%nodes(node)%shape
          do k=1,n(3)
            do j=1,n(2)
              do i=1,n(1)
                read(unit,*,iostat=ios) scratch(:payload_size)
                if (ios/=0) then
                  valid=.false.;exit scan_nodes
                end if
              end do
            end do
          end do
        end do scan_nodes
      end if
      if(valid.and.allocated(candidate%reaction_source%fields)) then
        if(.not.allocated(candidate%record_offsets)) then
          inquire(unit=unit,pos=footer_position,iostat=ios)
          valid=ios==0
        end if
        scan_reactions: do node=1,size(assigned)
          if(.not.valid) exit
          reaction_offsets(node)=footer_position;n=candidate%tree%nodes(node)%shape
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            if(i==1.and.j==1.and.k==1) then
              read(unit,*,pos=reaction_offsets(node),iostat=ios) scratch(:nvar-4)
            else
              read(unit,*,iostat=ios) scratch(:nvar-4)
            end if
            if(ios/=0) then
              valid=.false.;exit scan_reactions
            end if
          end do
          end do
          end do
          inquire(unit=unit,pos=footer_position,iostat=ios)
          valid=ios==0
        end do scan_reactions
      end if
      if(valid.and.allocated(candidate%heat_release%fields)) then
        if(.not.allocated(candidate%record_offsets)) then
          inquire(unit=unit,pos=footer_position,iostat=ios)
          valid=ios==0
        end if
        scan_heat: do node=1,size(assigned)
          if(.not.valid) exit
          heat_offsets(node)=footer_position;n=candidate%tree%nodes(node)%shape
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            if(i==1.and.j==1.and.k==1) then
              read(unit,*,pos=heat_offsets(node),iostat=ios) scratch(:1)
            else
              read(unit,*,iostat=ios) scratch(:1)
            end if
            if(ios/=0) then
              valid=.false.;exit scan_heat
            end if
          end do
          end do
          end do
          inquire(unit=unit,pos=footer_position,iostat=ios)
          valid=ios==0
        end do scan_heat
      end if
      if (valid) then
        if (allocated(candidate%record_offsets)) then
          read(unit,'(a)',pos=footer_position,iostat=ios) line
        else
          read(unit,'(a)',iostat=ios) line
        end if
        valid=ios==0.and.trim(line)=='END_AMR_TREE_CHECKPOINT'
        if (valid) then
          read(unit,'(a)',iostat=ios) line
          valid=ios<0
        end if
      end if
    end if
    call agree_checkpoint(comm,valid)
    if (.not.valid) go to 900
    call MPI_Bcast(offsets,size(offsets),MPI_OFFSET,0,comm,status)
    valid=status==MPI_SUCCESS;call agree_checkpoint(comm,valid)
    if (.not.valid) go to 900
    if(allocated(candidate%reaction_source%fields)) then
      call MPI_Bcast(reaction_offsets,size(reaction_offsets),MPI_OFFSET,0,comm,status)
      valid=status==MPI_SUCCESS;call agree_checkpoint(comm,valid)
      if(.not.valid) go to 900
    end if
    if(allocated(candidate%heat_release%fields)) then
      call MPI_Bcast(heat_offsets,size(heat_offsets),MPI_OFFSET,0,comm,status)
      valid=status==MPI_SUCCESS;call agree_checkpoint(comm,valid)
      if(.not.valid) go to 900
    end if
    allocate(candidate%state%fields(size(assigned)),candidate%temperature%fields(size(assigned)),stat=ios)
    valid=ios==0;call agree_checkpoint(comm,valid)
    if (.not.valid) go to 900
    nodes: do node=1,size(assigned)
      if (assigned(node)/=rank) cycle
      n=candidate%tree%nodes(node)%shape
      allocate(candidate%state%fields(node)%values(nvar,n(1),n(2),n(3)), &
        candidate%temperature%fields(node)%values(1,n(1),n(2),n(3)),stat=ios)
      if (ios/=0) then
        valid=.false.;exit
      end if
      if(allocated(candidate%internal_energy%fields)) then
        allocate(candidate%internal_energy%fields(node)%values(1,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        endif
      endif
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            if (i==1.and.j==1.and.k==1) then
              read(unit,*,pos=offsets(node),iostat=ios) scratch(:payload_size)
            else
              read(unit,*,iostat=ios) scratch(:payload_size)
            end if
            if (ios/=0) then
              valid=.false.;exit nodes
            end if
            candidate%state%fields(node)%values(:,i,j,k)=scratch(:nvar)
            candidate%temperature%fields(node)%values(1,i,j,k)=scratch(nvar+1)
            if(allocated(candidate%internal_energy%fields)) &
              candidate%internal_energy%fields(node)%values(1,i,j,k)=scratch(nvar+2)
          end do
        end do
      end do
      if (allocated(candidate%record_offsets)) then
        inquire(unit=unit,pos=position,iostat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
        if (position-offsets(node)/=candidate%record_offsets(node+1)-candidate%record_offsets(node)) then
          valid=.false.;exit
        end if
      end if
    end do nodes
    call agree_checkpoint(comm,valid)
    if(.not.valid) go to 900
    if(allocated(candidate%reaction_source%fields)) then
      reaction_nodes: do node=1,size(assigned)
        if(assigned(node)/=rank) cycle
        n=candidate%tree%nodes(node)%shape
        allocate(candidate%reaction_source%fields(node)%values(nvar-4,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        end if
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(i==1.and.j==1.and.k==1) then
            read(unit,*,pos=reaction_offsets(node),iostat=ios) candidate%reaction_source%fields(node)%values(:,i,j,k)
          else
            read(unit,*,iostat=ios) candidate%reaction_source%fields(node)%values(:,i,j,k)
          end if
          if(ios/=0) then
            valid=.false.;exit reaction_nodes
          end if
        end do
        end do
        end do
      end do reaction_nodes
      call agree_checkpoint(comm,valid)
    end if
    if(allocated(candidate%heat_release%fields)) then
      heat_nodes: do node=1,size(assigned)
        if(assigned(node)/=rank) cycle
        n=candidate%tree%nodes(node)%shape
        allocate(candidate%heat_release%fields(node)%values(1,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        end if
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(i==1.and.j==1.and.k==1) then
            read(unit,*,pos=heat_offsets(node),iostat=ios) candidate%heat_release%fields(node)%values(:,i,j,k)
          else
            read(unit,*,iostat=ios) candidate%heat_release%fields(node)%values(:,i,j,k)
          end if
          if(ios/=0) then
            valid=.false.;exit heat_nodes
          end if
        end do
        end do
        end do
      end do heat_nodes
      call agree_checkpoint(comm,valid)
    end if
900 close(unit,iostat=ios)
    valid=valid.and.ios==0;call agree_checkpoint(comm,valid)
    if (.not.valid) return
    call prepare_checkpoint(path,species,context,candidate,assigned,candidate%state,candidate%temperature, &
      comm,header,rank,valid,frame,domain_lower,root_spacing,mol_energy_reset)
    if (.not.valid) return
    snapshot=candidate
    call move_alloc(assigned,owners)
    ok=.true.
  end subroutine

  ! Metadata carries tree/ledger/counters; its state/temperature components are unused.
  subroutine write_mpi_amr_tree_checkpoint_3d(path,species,context,metadata,owners,state,temperature,comm,ok,indexed, &
      frame,domain_lower,root_spacing,mol_energy_reset)
    logical, optional, intent(in) :: mol_energy_reset
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3),root_spacing(3)
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_checkpoint_3d), intent(in) :: metadata
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    logical, optional, intent(in) :: indexed
    integer(int64), allocatable :: record_offsets(:)
    character(len=64) :: index_line
    logical :: with_index
    type(MPI_File) :: file
    type(MPI_Status) :: io_status
    integer(MPI_OFFSET_KIND), allocatable :: offsets(:),reaction_offsets(:),heat_offsets(:)
    integer(MPI_OFFSET_KIND) :: total,cells,position
    character(len=:), allocatable :: header,row,buffer,source_row,heat_row,temporary
    integer(c_int) :: publish_status
    character(len=*), parameter :: footer='END_AMR_TREE_CHECKPOINT'//achar(10)
    integer :: rank,status,node,nvar,i,j,k,n(3),axis,ios,used,written
    logical :: valid
    ok=.false.;nvar=reactive_nvar(size(species))
    call prepare_checkpoint(path,species,context,metadata,owners,state,temperature,comm,header,rank,valid, &
      frame,domain_lower,root_spacing,mol_energy_reset)
    if (.not.valid) return
    with_index=.false.
    if (present(indexed)) with_index=indexed
    index_line='V1'
    if (with_index) index_line='V2'
    call same_checkpoint_string(comm,index_line,valid)
    if (.not.valid) return
    if (with_index) then
      call checkpoint_record_offsets_3d(metadata%tree,nvar,record_offsets,valid,allocated(metadata%internal_energy%fields))
      call agree_checkpoint(comm,valid)
      if (.not.valid) return
      if(allocated(metadata%heat_release%fields)) then
        header='PELEF_AMR_TREE_3D_CHECKPOINT_V6'//header(len('PELEF_AMR_TREE_3D_CHECKPOINT_V1')+1:)
      else if(allocated(metadata%source)) then
        header='PELEF_AMR_TREE_3D_CHECKPOINT_V4'//header(len('PELEF_AMR_TREE_3D_CHECKPOINT_V1')+1:)
      else
        header='PELEF_AMR_TREE_3D_CHECKPOINT_V2'//header(len('PELEF_AMR_TREE_3D_CHECKPOINT_V1')+1:)
      end if
      do node=1,size(record_offsets)
        write(index_line,'(i0)') record_offsets(node)
        header=header//trim(index_line)//achar(10)
      end do
    end if
    if(allocated(metadata%previous_timestep)) &
      header='PELEF_AMR_TREE_TIMESTEP_HISTORY_V1'//achar(10)//header
    if(allocated(metadata%internal_energy%fields)) then
      if(allocated(metadata%pressure_work_change)) then
        write(index_line,'(2(es27.18e3,1x))') metadata%energy_correction,metadata%pressure_work_change
        header='PELEF_AMR_TREE_INTERNAL_ENERGY_V2'//achar(10)//trim(index_line)//achar(10)//header
      else
      write(index_line,'(es27.18e3)') metadata%energy_correction
      header='PELEF_AMR_TREE_INTERNAL_ENERGY_V1'//achar(10)//trim(index_line)//achar(10)//header
      endif
    endif
    call same_checkpoint_string(comm,header,valid)
    if (.not.valid) return
    allocate(offsets(size(owners)),reaction_offsets(size(owners)),heat_offsets(size(owners)),stat=ios)
    valid=ios==0;call agree_checkpoint(comm,valid)
    if (.not.valid) return
    allocate(character(len=28*(nvar+1+merge(1,0,allocated(metadata%internal_energy%fields)))-1)::row,stat=ios)
    valid=ios==0;call agree_checkpoint(comm,valid)
    if (.not.valid) return
    allocate(character(len=32*(len(row)+1))::buffer,stat=ios)
    valid=ios==0;call agree_checkpoint(comm,valid)
    if (.not.valid) return
    total=int(len(header),MPI_OFFSET_KIND)
    do node=1,size(owners)
      cells=1
      do axis=1,3
        if (cells>huge(total)/int(metadata%tree%nodes(node)%shape(axis),MPI_OFFSET_KIND)) return
        cells=cells*int(metadata%tree%nodes(node)%shape(axis),MPI_OFFSET_KIND)
      end do
      if (cells>(huge(total)-total)/int(len(row)+1,MPI_OFFSET_KIND)) return
      offsets(node)=total;total=total+cells*int(len(row)+1,MPI_OFFSET_KIND)
    end do
    if(allocated(metadata%reaction_source%fields)) then
      allocate(character(len=28*(nvar-4)-1)::source_row,stat=ios)
      valid=ios==0;call agree_checkpoint(comm,valid)
      if(.not.valid) return
      do node=1,size(owners)
        cells=product(int(metadata%tree%nodes(node)%shape,MPI_OFFSET_KIND))
        if(cells>(huge(total)-total)/int(len(source_row)+1,MPI_OFFSET_KIND)) return
        reaction_offsets(node)=total;total=total+cells*int(len(source_row)+1,MPI_OFFSET_KIND)
      end do
    end if
    if(allocated(metadata%heat_release%fields)) then
      allocate(character(len=27)::heat_row,stat=ios)
      valid=ios==0;call agree_checkpoint(comm,valid)
      if(.not.valid) return
      do node=1,size(owners)
        cells=product(int(metadata%tree%nodes(node)%shape,MPI_OFFSET_KIND))
        if(cells>(huge(total)-total)/int(len(heat_row)+1,MPI_OFFSET_KIND)) return
        heat_offsets(node)=total;total=total+cells*int(len(heat_row)+1,MPI_OFFSET_KIND)
      end do
    end if
    if (total>huge(total)-int(len(footer),MPI_OFFSET_KIND)) return
    temporary=trim(path)//'.tmp'
    call MPI_File_open(comm,temporary,ior(MPI_MODE_WRONLY,ior(MPI_MODE_CREATE,MPI_MODE_EXCL)), &
      MPI_INFO_NULL,file,status)
    if (status/=MPI_SUCCESS) return
    call MPI_File_set_errhandler(file,MPI_ERRORS_RETURN,status)
    valid=status==MPI_SUCCESS;call agree_checkpoint(comm,valid)
    if (valid) then
      call MPI_File_set_size(file,total+int(len(footer),MPI_OFFSET_KIND),status)
      valid=status==MPI_SUCCESS;call agree_checkpoint(comm,valid)
    end if
    if (valid.and.rank==0) then
      call write_text(0_MPI_OFFSET_KIND,header,valid)
      if (valid) call write_text(total,footer,valid)
    end if
    call agree_checkpoint(comm,valid)
    if (valid) then
      nodes: do node=1,size(owners)
        if (owners(node)/=rank) cycle
        n=metadata%tree%nodes(node)%shape;position=offsets(node);used=0
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              if(allocated(metadata%internal_energy%fields)) then
              write(row,'(*(es27.18e3,1x))',iostat=ios) state%fields(node)%values(:,i,j,k), &
                temperature%fields(node)%values(1,i,j,k),metadata%internal_energy%fields(node)%values(1,i,j,k)
              else
              write(row,'(*(es27.18e3,1x))',iostat=ios) state%fields(node)%values(:,i,j,k), &
                temperature%fields(node)%values(1,i,j,k)
              endif
              if (ios/=0) then
                valid=.false.;exit nodes
              end if
              if (used+len(row)+1>len(buffer)) then
                call flush_chunk()
                if (.not.valid) exit nodes
              end if
              buffer(used+1:used+len(row)+1)=row//achar(10);used=used+len(row)+1
            end do
          end do
        end do
        call flush_chunk()
        if (.not.valid) exit
      end do nodes
    end if
    call agree_checkpoint(comm,valid)
    if(valid.and.allocated(metadata%reaction_source%fields)) then
      source_nodes: do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=metadata%tree%nodes(node)%shape;position=reaction_offsets(node);used=0
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          write(source_row,'(*(es27.18e3,1x))',iostat=ios) metadata%reaction_source%fields(node)%values(:,i,j,k)
          if(ios/=0) then
            valid=.false.;exit source_nodes
          end if
          if(used+len(source_row)+1>len(buffer)) then
            call flush_chunk()
            if(.not.valid) exit source_nodes
          end if
          buffer(used+1:used+len(source_row)+1)=source_row//achar(10);used=used+len(source_row)+1
        end do
        end do
        end do
        call flush_chunk()
        if(.not.valid) exit
      end do source_nodes
    end if
    if(valid.and.allocated(metadata%heat_release%fields)) then
      heat_nodes: do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=metadata%tree%nodes(node)%shape;position=heat_offsets(node);used=0
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          write(heat_row,'(es27.18e3)',iostat=ios) metadata%heat_release%fields(node)%values(:,i,j,k)
          if(ios/=0) then
            valid=.false.;exit heat_nodes
          end if
          if(used+len(heat_row)+1>len(buffer)) then
            call flush_chunk()
            if(.not.valid) exit heat_nodes
          end if
          buffer(used+1:used+len(heat_row)+1)=heat_row//achar(10);used=used+len(heat_row)+1
        end do
        end do
        end do
        call flush_chunk()
        if(.not.valid) exit
      end do heat_nodes
    end if
    call agree_checkpoint(comm,valid)
    call MPI_File_close(file,status)
    ok=valid.and.status==MPI_SUCCESS;call agree_checkpoint(comm,ok)
    if(rank==0) then
      if(ok) then
        publish_status=checkpoint_rename(temporary//c_null_char,trim(path)//c_null_char)
        ok=publish_status==0
      endif
      if(.not.ok) publish_status=checkpoint_unlink(temporary//c_null_char)
    endif
    call agree_checkpoint(comm,ok)
  contains
    subroutine write_text(at,text,success)
      integer(MPI_OFFSET_KIND), intent(in) :: at
      character(len=*), intent(in) :: text
      logical, intent(out) :: success
      call MPI_File_write_at(file,at,text,len(text),MPI_CHARACTER,io_status,status)
      success=status==MPI_SUCCESS
      if (.not.success) return
      call MPI_Get_count(io_status,MPI_CHARACTER,written,status)
      success=status==MPI_SUCCESS.and.written==len(text)
    end subroutine
    subroutine flush_chunk()
      if (used==0) return
      call write_text(position,buffer(:used),valid)
      if (.not.valid) return
      position=position+int(used,MPI_OFFSET_KIND);used=0
    end subroutine

  end subroutine
  subroutine prepare_checkpoint(path,species,context,metadata,owners,state,temperature,comm,header,rank,ok, &
      frame,domain_lower,root_spacing,mol_energy_reset)
    logical, optional, intent(in) :: mol_energy_reset
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3),root_spacing(3)
    character(len=:), allocatable :: rotation_context
    real(dp), allocatable :: eos_state(:,:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: h(3),lower(3)
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_checkpoint_3d), intent(in) :: metadata
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(MPI_Comm), intent(in) :: comm
    character(len=:), allocatable, intent(out) :: header
    integer, intent(out) :: rank
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: synchronized,temps
    character(len=65536) :: line
    integer :: nvar,node,start,finish,ios
    logical :: valid,consistent,reset
    ok=.false.;nvar=reactive_nvar(size(species))
    call checkpoint_rotation_context(comm,context,frame,domain_lower,root_spacing,rotation_context,valid)
    if(.not.valid) return
    reset=.false.
    if(present(mol_energy_reset)) reset=mol_energy_reset
    write(line,'(l1)') reset
    call same_checkpoint_string(comm,trim(line),valid)
    if(.not.valid) return
    valid=.not.(reset.and.present(frame))
    call agree_checkpoint(comm,valid)
    if(.not.valid) return
    if(present(frame)) then
      call amr_tree_integer_layout_3d(metadata%tree,origin,domain,ratios,valid)
      call agree_checkpoint(comm,valid);if(.not.valid) return
    endif
    call clone_mpi_amr_tree_inputs_3d(metadata%tree,owners,size(species),state,temperature, &
      synchronized,temps,comm,rank,valid)
    if (.not.valid) return
    valid=allocated(metadata%initial).and.allocated(metadata%boundary).and.len_trim(path)>0
    valid=valid.and.valid_amr_tree_timestep_history_3d(metadata)
    if (valid) valid=size(metadata%initial)==nvar.and.size(metadata%boundary)==nvar
    if (valid) valid=all(ieee_is_finite(metadata%initial)).and.all(ieee_is_finite(metadata%boundary))
    valid=valid.and.all(ieee_is_finite([metadata%time,metadata%maximum_reflux, &
      metadata%minimum_theta,metadata%maximum_remap_error]))
    if (valid) valid=min(metadata%time,metadata%maximum_reflux,metadata%minimum_theta,metadata%maximum_remap_error)>=0 &
      .and.metadata%minimum_theta<=1
    valid=valid.and.min(metadata%steps,metadata%regrids)>=0.and.metadata%regrids<=metadata%steps
    if (valid) valid=.not.((metadata%steps==0.and.metadata%time/=0).or.(metadata%steps>0.and.metadata%time<=0))
    if (len(context)<1) then
      valid=.false.
    else if (context(len(context):)/=achar(10)) then
      valid=.false.
    else
      start=1
      do while(start<=len(context))
        finish=start+index(context(start:),achar(10))-2
        if (finish-start+1>=65536) valid=.false.
        start=finish+2
      end do
    end if
    call agree_checkpoint(comm,valid)
    if (.not.valid) return
    write(line,'(6l1)') allocated(metadata%source),allocated(metadata%reaction_source%fields), &
      allocated(metadata%heat_release%fields),allocated(metadata%internal_energy%fields),allocated(metadata%energy_correction), &
      allocated(metadata%pressure_work_change)
    call same_checkpoint_string(comm,trim(line),valid)
    if(.not.valid) return
    valid=allocated(metadata%internal_energy%fields).eqv.allocated(metadata%energy_correction)
    if(allocated(metadata%energy_correction)) valid=valid.and.ieee_is_finite(metadata%energy_correction)
    if(allocated(metadata%pressure_work_change)) then
      valid=valid.and.allocated(metadata%internal_energy%fields)
      valid=valid.and.ieee_is_finite(metadata%pressure_work_change)
    endif
    call agree_checkpoint(comm,valid)
    if(.not.valid) return
    if(allocated(metadata%pressure_work_change)) then
      write(line,'(es27.18e3)') metadata%pressure_work_change
      call same_checkpoint_string(comm,trim(line),valid)
      if(.not.valid) return
    endif
    if(allocated(metadata%internal_energy%fields)) then
      write(line,'(es27.18e3)') metadata%energy_correction
      call same_checkpoint_string(comm,trim(line),valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(metadata%tree,owners,metadata%internal_energy,1,comm,valid)
      if(.not.valid) return
      temps=metadata%internal_energy
      call average_down_mpi_amr_tree_3d(metadata%tree,owners,temps,1,comm,valid)
      if(.not.valid) return
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        if(.not.all(ieee_is_finite(metadata%internal_energy%fields(node)%values))) valid=.false.
        if(.not.checkpoint_tree_held_matches(state%fields(node)%values, &
          metadata%internal_energy%fields(node)%values(1,:,:,:),temps%fields(node)%values(1,:,:,:),reset)) valid=.false.
      enddo
      call agree_checkpoint(comm,valid)
      if(.not.valid) return
    endif
    if(allocated(metadata%source)) then
      valid=size(metadata%source)==nvar
      if(valid) valid=all(ieee_is_finite(metadata%source))
    else if(allocated(metadata%reaction_source%fields).or.allocated(metadata%heat_release%fields)) then
      valid=.false.
    end if
    call agree_checkpoint(comm,valid)
    if(.not.valid) return
    if(allocated(metadata%reaction_source%fields)) then
      call validate_mpi_amr_tree_local_3d(metadata%tree,owners,metadata%reaction_source,nvar-4,comm,valid)
      if(.not.valid) return
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        if(.not.all(ieee_is_finite(metadata%reaction_source%fields(node)%values))) valid=.false.
      end do
      call agree_checkpoint(comm,valid)
      if(.not.valid) return
    end if
    if(allocated(metadata%heat_release%fields)) then
      call validate_mpi_amr_tree_local_3d(metadata%tree,owners,metadata%heat_release,1,comm,valid)
      if(.not.valid) return
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        if(.not.all(ieee_is_finite(metadata%heat_release%fields(node)%values))) valid=.false.
      end do
      call agree_checkpoint(comm,valid)
      if(.not.valid) return
    end if
    call average_down_mpi_amr_tree_3d(metadata%tree,owners,synchronized,nvar,comm,valid)
    if (.not.valid) return
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      if (any(abs(synchronized%fields(node)%values-state%fields(node)%values)> &
        64*epsilon(1._dp)*max(1._dp,abs(state%fields(node)%values)))) valid=.false.
      eos_state=state%fields(node)%values
      if(present(frame)) then
        h=root_spacing*real(domain(:,1),dp)/real(domain(:,metadata%tree%nodes(node)%level+1),dp)
        lower=domain_lower+real(origin(:,node),dp)*h
        call convert_regular_energy_3d(frame,lower,h,eos_state,0,.false.,consistent)
        if(.not.consistent) then
          valid=.false.;exit
        endif
      endif
      if(allocated(metadata%internal_energy%fields)) then
      call checkpoint_temperature_is_consistent(species,eos_state, &
        temperature%fields(node)%values(1,:,:,:),consistent,metadata%internal_energy%fields(node)%values(1,:,:,:))
      else
      call checkpoint_temperature_is_consistent(species,eos_state, &
        temperature%fields(node)%values(1,:,:,:),consistent)
      endif
      valid=valid.and.consistent
    end do
    call agree_checkpoint(comm,valid)
    if (.not.valid) return
    deallocate(synchronized%fields,temps%fields)
    header='PELEF_AMR_TREE_3D_CHECKPOINT_V1'//achar(10)
    if(allocated(metadata%source)) header='PELEF_AMR_TREE_3D_CHECKPOINT_V3'//achar(10)
    if(allocated(metadata%heat_release%fields)) header='PELEF_AMR_TREE_3D_CHECKPOINT_V5'//achar(10)
    write(line,'(i0)') len(rotation_context)
    header=header//trim(line)//achar(10)//rotation_context
    write(line,'(*(i0,1x))') size(owners),nvar,metadata%tree%nodes(1)%shape
    header=header//trim(line)//achar(10)
    write(line,'(4(es27.18e3,1x),2(i0,1x))') metadata%time,metadata%maximum_reflux, &
      metadata%minimum_theta,metadata%maximum_remap_error,metadata%steps,metadata%regrids
    header=header//trim(line)//achar(10)
    if(allocated(metadata%previous_timestep)) then
      write(line,'(es27.18e3,1x,l1)') metadata%previous_timestep,metadata%post_regrid
      header=header//trim(line)//achar(10)
    endif
    write(line,'(*(es27.18e3,1x))',iostat=ios) metadata%initial,metadata%boundary
    valid=ios==0;call agree_checkpoint(comm,valid)
    if (.not.valid) return
    header=header//trim(line)//achar(10)
    if(allocated(metadata%source)) then
      write(line,'(*(es27.18e3,1x))',iostat=ios) metadata%source
      valid=ios==0;call agree_checkpoint(comm,valid)
      if(.not.valid) return
      header=header//trim(line)//achar(10)
      write(line,'(i0)') merge(1,0,allocated(metadata%reaction_source%fields))
      header=header//trim(line)//achar(10)
    end if
    do node=2,size(owners)
      associate(p=>metadata%tree%nodes(node)%patch)
        write(line,'(*(i0,1x))') metadata%tree%nodes(node)%parent,p%coarse_i_lower,p%coarse_i_upper, &
          p%coarse_j_lower,p%coarse_j_upper,p%coarse_k_lower,p%coarse_k_upper,p%refinement_ratio
      end associate
      header=header//trim(line)//achar(10)
    end do
    call same_checkpoint_string(comm,trim(path),valid)
    if (.not.valid) return
    call same_checkpoint_string(comm,header,valid)
    if (.not.valid) return
    ok=.true.
  end subroutine
    subroutine agree_checkpoint(comm,condition)
      type(MPI_Comm), intent(in) :: comm
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if (ierr==MPI_SUCCESS) condition=result
    end subroutine
    ! Bind storage convention and physical cell coordinates into the archive.
  ! Every caller validates parameters collectively before inspecting payloads.
  subroutine checkpoint_rotation_context(comm,context,frame,domain_lower,root_spacing,bound,ok)
    type(MPI_Comm), intent(in) :: comm
    character(len=*), intent(in) :: context
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3),root_spacing(3)
    character(len=:), allocatable, intent(out) :: bound
    logical, intent(out) :: ok
    integer :: flags(3),low(3),high(3),status
    character(len=512) :: line
    logical :: valid
    ok=.false.
    flags=[merge(1,0,present(frame)),merge(1,0,present(domain_lower)),merge(1,0,present(root_spacing))]
    call MPI_Allreduce(flags,low,3,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,3,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    bound=context
    if(present(frame)) then
      valid=present(domain_lower).and.present(root_spacing).and.frame%is_valid()
      call agree_checkpoint(comm,valid);if(.not.valid) return
      valid=all(ieee_is_finite([domain_lower,root_spacing])).and.all(root_spacing>0)
      call agree_checkpoint(comm,valid);if(.not.valid) return
      write(line,'(i0,10(1x,es27.18e3))') frame%axis,frame%angular_speed,frame%axis_point,domain_lower,root_spacing
      call same_checkpoint_string(comm,trim(line),valid);if(.not.valid) return
      bound=bound//'CHECKPOINT_ROTATION_V1 '//trim(line)//achar(10)
    endif
    call same_checkpoint_string(comm,bound,valid);ok=valid
  end subroutine checkpoint_rotation_context

  subroutine same_checkpoint_string(comm,value,success)
      type(MPI_Comm), intent(in) :: comm
      character(len=*), intent(in) :: value
      logical, intent(out) :: success
      character(len=:), allocatable :: reference
      integer :: minimum,maximum,ierr,allocation_status
      success=.false.
      call MPI_Allreduce(len(value),minimum,1,MPI_INTEGER,MPI_MIN,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Allreduce(len(value),maximum,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if (ierr/=MPI_SUCCESS.or.minimum/=maximum) return
      allocate(character(len=minimum)::reference,stat=allocation_status)
      success=allocation_status==0;call agree_checkpoint(comm,success)
      if (.not.success) return
      reference=value
      call MPI_Bcast(reference,minimum,MPI_CHARACTER,0,comm,ierr)
      success=ierr==MPI_SUCCESS.and.reference==value;call agree_checkpoint(comm,success)
    end subroutine
end module mpi_amr_patch_tree_checkpoint_3d_mod
