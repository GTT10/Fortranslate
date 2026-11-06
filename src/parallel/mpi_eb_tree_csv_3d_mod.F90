module mpi_eb_tree_csv_3d_mod
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d
  implicit none
  private
  public :: write_mpi_eb_tree_csv_3d
contains
  subroutine write_mpi_eb_tree_csv_3d(path,tree,owners,species,state,temperature,geometry,comm,ok,heat_release,internal_energy)
    character(len=*), intent(in) :: path
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(amr_tree_state_3d), optional, intent(in) :: heat_release,internal_energy
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(MPI_File) :: file
    type(MPI_Status) :: io_status
    integer(MPI_OFFSET_KIND), allocatable :: local_bytes(:),bytes(:),offset(:)
    integer(MPI_OFFSET_KIND) :: total,position
    character(len=:), allocatable :: header,row,buffer
    real(dp) :: volume
    integer :: rank,status,ios,node,pass,i,j,k,n(3),used,length,written,c,nv,heat_flag,heat_min,heat_max,held_flag,flags
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    held_flag=merge(1,0,present(internal_energy))
    heat_flag=merge(1,0,present(heat_release))
    flags=heat_flag+2*held_flag
    call MPI_Allreduce(flags,heat_min,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,heat_max,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.heat_min/=heat_max) return
    if(present(heat_release)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,heat_release,1,comm,valid)
      if(.not.valid) return
    end if
    if(present(internal_energy)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,valid)
      if(.not.valid) return
    endif
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,valid)
    if (.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,reactive_nvar(size(species)),comm,valid)
    if (.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,valid)
    if (.not.valid) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    valid=len_trim(path)>0
    call agree(valid)
    if (.not.valid) return
    call same_string(trim(path),valid)
    if (.not.valid) return
    header='node,level,i,j,k,fluid_volume,temperature,rho,rhou,rhov,rhow,rhoE'
    do c=1,size(species)
      header=header//',rhoY_'//trim(species(c)%name)
    end do
    if(present(heat_release)) header=header//',heat_release_W_m3'
    if(present(internal_energy)) header=header//',rho_internal_energy_J_m3'
    header=header//new_line('a')
    call same_string(header,valid)
    if (.not.valid) return
    allocate(local_bytes(size(owners)),bytes(size(owners)),offset(size(owners)),stat=ios)
    valid=ios==0;call agree(valid)
    if (.not.valid) return
    allocate(character(len=60+27*(nv+2+heat_flag+held_flag))::row,stat=ios)
    valid=ios==0;call agree(valid)
    if (.not.valid) return
    allocate(character(len=32*(len(row)+1))::buffer,stat=ios)
    valid=ios==0;call agree(valid)
    if (.not.valid) return
    local_bytes=0
    ! Validate/measure all rows before opening or truncating the destination.
    ! The second pass regenerates the same rows into bounded owner-local chunks.
    do pass=1,2
      valid=.true.
      nodes: do node=1,size(owners)
        if (owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        used=0
        if (pass==2) position=offset(node)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              if (covered(tree,node,[i,j,k])) cycle
              if (geometry(node)%volume_fraction(i,j,k)<=0._dp) cycle
              associate(g=>geometry(node))
              volume=g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz
              end associate
              valid=ieee_is_finite(volume).and.ieee_is_finite(temperature%fields(node)%values(1,i,j,k)).and. &
                all(ieee_is_finite(state%fields(node)%values(:,i,j,k)))
              if(present(heat_release)) valid=valid.and.ieee_is_finite(heat_release%fields(node)%values(1,i,j,k))
              if(present(internal_energy)) valid=valid.and.ieee_is_finite(internal_energy%fields(node)%values(1,i,j,k))
              if (.not.valid) exit nodes
              write(row,'(i0,4(",",i0),*(:,",",es26.17e3))',iostat=ios) node,tree%nodes(node)%level,i,j,k, &
                volume,temperature%fields(node)%values(1,i,j,k),state%fields(node)%values(:,i,j,k)
              valid=ios==0
              if (.not.valid) exit nodes
              if(present(heat_release)) then
                length=len_trim(row)
                write(row(length+1:),'(",",es26.17e3)',iostat=ios) heat_release%fields(node)%values(1,i,j,k)
                if(ios/=0) then
                  valid=.false.;exit nodes
                end if
              end if
              if(present(internal_energy)) then
                length=len_trim(row)
                write(row(length+1:),'(",",es26.17e3)',iostat=ios) internal_energy%fields(node)%values(1,i,j,k)
                if(ios/=0) then
                  valid=.false.;exit nodes
                endif
              endif
              length=len_trim(row)+1
              if (pass==1) then
                if (local_bytes(node)>huge(total)-int(length,MPI_OFFSET_KIND)) then
                  valid=.false.;exit nodes
                end if
                local_bytes(node)=local_bytes(node)+int(length,MPI_OFFSET_KIND)
              else
                if (used+length>len(buffer)) then
                  call flush_chunk()
                  if (.not.valid) exit nodes
                end if
                buffer(used+1:used+length)=trim(row)//new_line('a')
                used=used+length
              end if
            end do
          end do
        end do
        if (pass==2) then
          call flush_chunk()
          if (.not.valid) exit
        end if
      end do nodes
      call agree(valid)
      if (pass==2) exit
      if (.not.valid) return
      call MPI_Allreduce(local_bytes,bytes,size(owners),MPI_OFFSET,MPI_SUM,comm,status)
      if (status/=MPI_SUCCESS) return
      total=int(len(header),MPI_OFFSET_KIND)
      do node=1,size(owners)
        if (bytes(node)>huge(total)-total) return
        offset(node)=total;total=total+bytes(node)
      end do
      call MPI_File_open(comm,trim(path),ior(MPI_MODE_WRONLY,MPI_MODE_CREATE),MPI_INFO_NULL,file,status)
      if (status/=MPI_SUCCESS) return
      call MPI_File_set_errhandler(file,MPI_ERRORS_RETURN,status)
      valid=status==MPI_SUCCESS;call agree(valid)
      if (valid) then
        call MPI_File_set_size(file,total,status)
        valid=status==MPI_SUCCESS;call agree(valid)
      end if
      if (valid.and.rank==0) then
        call MPI_File_write_at(file,0_MPI_OFFSET_KIND,header,len(header),MPI_CHARACTER,io_status,status)
        valid=status==MPI_SUCCESS
        if (valid) then
          call MPI_Get_count(io_status,MPI_CHARACTER,written,status)
          valid=status==MPI_SUCCESS.and.written==len(header)
        end if
      end if
      call agree(valid)
      if (.not.valid) exit
    end do
    call MPI_File_close(file,status)
    ok=valid.and.status==MPI_SUCCESS
    call agree(ok)
  contains
    subroutine flush_chunk()
      if (used==0) return
      call MPI_File_write_at(file,position,buffer(:used),used,MPI_CHARACTER,io_status,status)
      valid=status==MPI_SUCCESS
      if (.not.valid) return
      call MPI_Get_count(io_status,MPI_CHARACTER,written,status)
      valid=status==MPI_SUCCESS.and.written==used
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
      success=ierr==MPI_SUCCESS.and.reference==value
      call agree(success)
    end subroutine
  end subroutine
end module mpi_eb_tree_csv_3d_mod
