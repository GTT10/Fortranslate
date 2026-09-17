module mpi_amr_patch_tree_checkpoint_3d_mod
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_tree_state_3d
  use amr_patch_tree_checkpoint_3d_mod, only: amr_tree_checkpoint_3d
  use amr_reactive_3d_checkpoint_mod, only: checkpoint_temperature_is_consistent
  use mpi_amr_patch_tree_chemistry_3d_mod, only: clone_mpi_amr_tree_inputs_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  implicit none
  private
  public :: write_mpi_amr_tree_checkpoint_3d
contains
  ! Metadata carries tree/ledger/counters; its state/temperature components are unused.
  subroutine write_mpi_amr_tree_checkpoint_3d(path,species,context,metadata,owners,state,temperature,comm,ok)
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_checkpoint_3d), intent(in) :: metadata
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: synchronized,temps
    type(MPI_File) :: file
    type(MPI_Status) :: io_status
    integer(MPI_OFFSET_KIND), allocatable :: offsets(:)
    integer(MPI_OFFSET_KIND) :: total,cells,position
    character(len=:), allocatable :: header,row,buffer
    character(len=65536) :: line
    character(len=*), parameter :: footer='END_AMR_TREE_CHECKPOINT'//achar(10)
    integer :: rank,status,node,nvar,i,j,k,n(3),axis,ios,used,written,start,finish
    logical :: valid,consistent
    ok=.false.;nvar=reactive_nvar(size(species))
    call clone_mpi_amr_tree_inputs_3d(metadata%tree,owners,size(species),state,temperature, &
      synchronized,temps,comm,rank,valid)
    if (.not.valid) return
    valid=allocated(metadata%initial).and.allocated(metadata%boundary).and.len_trim(path)>0
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
    call agree(valid)
    if (.not.valid) return
    call average_down_mpi_amr_tree_3d(metadata%tree,owners,synchronized,nvar,comm,valid)
    if (.not.valid) return
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      if (any(abs(synchronized%fields(node)%values-state%fields(node)%values)> &
        64*epsilon(1._dp)*max(1._dp,abs(state%fields(node)%values)))) valid=.false.
      call checkpoint_temperature_is_consistent(species,state%fields(node)%values, &
        temperature%fields(node)%values(1,:,:,:),consistent)
      valid=valid.and.consistent
    end do
    call agree(valid)
    if (.not.valid) return
    deallocate(synchronized%fields,temps%fields)
    header='PELEF_AMR_TREE_3D_CHECKPOINT_V1'//achar(10)
    write(line,'(i0)') len(context)
    header=header//trim(line)//achar(10)//context
    write(line,'(*(i0,1x))') size(owners),nvar,metadata%tree%nodes(1)%shape
    header=header//trim(line)//achar(10)
    write(line,'(4(es27.18e3,1x),2(i0,1x))') metadata%time,metadata%maximum_reflux, &
      metadata%minimum_theta,metadata%maximum_remap_error,metadata%steps,metadata%regrids
    header=header//trim(line)//achar(10)
    write(line,'(*(es27.18e3,1x))',iostat=ios) metadata%initial,metadata%boundary
    valid=ios==0;call agree(valid)
    if (.not.valid) return
    header=header//trim(line)//achar(10)
    do node=2,size(owners)
      associate(p=>metadata%tree%nodes(node)%patch)
        write(line,'(*(i0,1x))') metadata%tree%nodes(node)%parent,p%coarse_i_lower,p%coarse_i_upper, &
          p%coarse_j_lower,p%coarse_j_upper,p%coarse_k_lower,p%coarse_k_upper,p%refinement_ratio
      end associate
      header=header//trim(line)//achar(10)
    end do
    call same_string(trim(path),valid)
    if (.not.valid) return
    call same_string(header,valid)
    if (.not.valid) return
    allocate(offsets(size(owners)),stat=ios)
    valid=ios==0;call agree(valid)
    if (.not.valid) return
    allocate(character(len=28*(nvar+1)-1)::row,stat=ios)
    valid=ios==0;call agree(valid)
    if (.not.valid) return
    allocate(character(len=32*(len(row)+1))::buffer,stat=ios)
    valid=ios==0;call agree(valid)
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
    if (total>huge(total)-int(len(footer),MPI_OFFSET_KIND)) return
    call MPI_File_open(comm,trim(path),ior(MPI_MODE_WRONLY,MPI_MODE_CREATE),MPI_INFO_NULL,file,status)
    if (status/=MPI_SUCCESS) return
    call MPI_File_set_errhandler(file,MPI_ERRORS_RETURN,status)
    valid=status==MPI_SUCCESS;call agree(valid)
    if (valid) then
      call MPI_File_set_size(file,total+int(len(footer),MPI_OFFSET_KIND),status)
      valid=status==MPI_SUCCESS;call agree(valid)
    end if
    if (valid.and.rank==0) then
      call write_text(0_MPI_OFFSET_KIND,header,valid)
      if (valid) call write_text(total,footer,valid)
    end if
    call agree(valid)
    if (valid) then
      nodes: do node=1,size(owners)
        if (owners(node)/=rank) cycle
        n=metadata%tree%nodes(node)%shape;position=offsets(node);used=0
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              write(row,'(*(es27.18e3,1x))',iostat=ios) state%fields(node)%values(:,i,j,k), &
                temperature%fields(node)%values(1,i,j,k)
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
    call agree(valid)
    call MPI_File_close(file,status)
    ok=valid.and.status==MPI_SUCCESS;call agree(ok)
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
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if (ierr==MPI_SUCCESS) condition=result
    end subroutine
    subroutine same_string(value,success)
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
      success=allocation_status==0;call agree(success)
      if (.not.success) return
      reference=value
      call MPI_Bcast(reference,minimum,MPI_CHARACTER,0,comm,ierr)
      success=ierr==MPI_SUCCESS.and.reference==value;call agree(success)
    end subroutine
  end subroutine
end module mpi_amr_patch_tree_checkpoint_3d_mod
