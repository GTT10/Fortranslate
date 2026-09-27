module mpi_eb_tree_checkpoint_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: iostat_end,iostat_eor
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use eb_reactive_tree_checkpoint_3d_mod, only: eb_tree_checkpoint_3d
  use reactive_eb_3d_checkpoint_mod, only: write_reactive_eb_3d_checkpoint
  use reactive_eb_3d_driver_mod, only: reactive_eb_integrals_3d,reactive_eb_element_integrals_3d, &
    reactive_eb_composition_integrals_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d,integrate_mpi_eb_tree_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  implicit none
  private
  public :: write_mpi_eb_tree_checkpoint_3d
contains
  ! Each owner stages its existing serial patch records in a private scratch file.
  ! Only record sizes are reduced; payloads are written directly in node order.
  subroutine write_mpi_eb_tree_checkpoint_3d(path,prefix_unit,species,reactions,transport,config,tree,owners, &
      geometry,state,temperature,record,comm,ok,bundle_sha256,chemistry_integrator,base_mole_fractions, &
      element_names,element_atoms)
    character(len=*), intent(in) :: path
    integer, intent(in) :: prefix_unit,owners(:)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_eb_3d_config), intent(in) :: config
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(eb_tree_checkpoint_3d), intent(in) :: record
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: bundle_sha256,chemistry_integrator,element_names(:)
    real(dp), optional, intent(in) :: base_mole_fractions(:),element_atoms(:,:)
    type(amr_tree_state_3d) :: averaged,at
    type(reactive_eb_3d_config) :: cfg
    type(MPI_File) :: file
    type(MPI_Status) :: io_status
    integer(MPI_OFFSET_KIND), allocatable :: local_bytes(:),bytes(:),offset(:),local_start(:)
    integer(MPI_OFFSET_KIND) :: cursor,header_bytes,total,position,remaining,source_position
    real(dp), allocatable :: current(:),integral(:),l1(:),elements(:)
    integer :: rank,status,ios,stage,patch,node,nv,ne,descriptor(12),amount,written,close_status,i,j,k,n(3)
    character(len=4096) :: chunk
    character(len=512) :: message
    character(len=:), allocatable :: header,index_text
    character(len=21) :: index_line
    character(len=*), parameter :: ending='END_EB_TREE_CHECKPOINT'
    logical :: valid,opened,sdc
    ok=.false.;opened=.false.;nv=reactive_nvar(size(species));ne=3
    sdc=config%time_integrator=='sdc'
    if (present(element_names)) ne=size(element_names)
    allocate(current(nv),integral(nv),l1(nv),elements(ne))
    averaged=state;at=temperature
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,averaged,at,comm,valid)
    if (.not.valid) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      if (any(abs(averaged%fields(node)%values-state%fields(node)%values)> &
        2.e-12_dp*max(1._dp,abs(averaged%fields(node)%values)))) valid=.false.
    end do
    call agree(valid)
    if (.not.valid) return
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,state,current,comm,valid)
    if (.not.valid) return
    valid=valid_record().and.len_trim(path)>0
    call agree(valid)
    if (.not.valid) return
    call same_string(trim(path),valid)
    if (.not.valid) return
    block
      character(len=128) :: control_text
      write(control_text,*) config%time_integrator,config%chemistry_enabled,config%sdc_iterations,config%sdc_reactor
      call same_string(control_text,valid)
      if(.not.valid) return
    end block
    if(sdc.and.config%chemistry_enabled) then
      call validate_mpi_amr_tree_local_3d(tree,owners,record%reaction_source,nv-4,comm,valid)
      if(.not.valid) return
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        if(.not.all(ieee_is_finite(record%reaction_source%fields(node)%values))) valid=.false.
      end do
      call agree(valid)
      if(.not.valid) return
    end if
    allocate(local_bytes(size(owners)),bytes(size(owners)),offset(size(owners)),local_start(size(owners)))
    local_bytes=0
    open(newunit=stage,status='scratch',access='stream',form='unformatted',iostat=ios)
    valid=ios==0;call agree(valid)
    if (.not.valid) then
      if (ios==0) close(stage)
      return
    end if
    ! Caller supplied only the application policy envelope; append tree metadata.
    if(sdc) then
      write(prefix_unit,'(a)',iostat=ios) 'PELEF_EB_TREE_CHECKPOINT_V4'
    else
      write(prefix_unit,'(a)',iostat=ios) 'PELEF_EB_TREE_CHECKPOINT_V2'
    end if
    if (ios==0) write(prefix_unit,*,iostat=ios) size(owners)
    do node=1,size(owners)
      if (ios/=0) exit
      associate(n=>tree%nodes(node),p=>tree%nodes(node)%patch)
      descriptor=[n%parent,n%level,n%shape,p%refinement_ratio,p%coarse_i_lower,p%coarse_i_upper, &
        p%coarse_j_lower,p%coarse_j_upper,p%coarse_k_lower,p%coarse_k_upper]
      end associate
      write(prefix_unit,*,iostat=ios) descriptor
    end do
    if (ios==0) write(prefix_unit,'(es27.18e3,1x,i0,2(1x,es27.18e3))',iostat=ios) &
      record%time,record%steps,record%minimum_dt,record%minimum_theta
    if (ios==0) write(prefix_unit,'(*(es27.18e3,1x))',iostat=ios) &
      record%initial,record%scale,record%boundary,record%wall,record%chemical
    if(sdc) then
      if(ios==0) write(prefix_unit,*,iostat=ios) config%sdc_iterations, &
          merge(1,0,config%chemistry_enabled)+merge(2,0,config%sdc_reactor=="cvode")
      if(ios==0) write(prefix_unit,'(*(es27.18e3,1x))',iostat=ios) config%external_forcing
      if(ios==0) write(prefix_unit,'(*(es27.18e3,1x))',iostat=ios) record%forcing
    end if
    valid=ios==0;cursor=0
    if (valid) call stage_formatted(prefix_unit,valid)
    call agree(valid)
    if (.not.valid) go to 900
    header_bytes=cursor
    valid=header_bytes<=int(huge(1),MPI_OFFSET_KIND)
    call agree(valid)
    if (.not.valid) go to 900
    allocate(character(len=int(header_bytes))::header,stat=ios)
    valid=ios==0;call agree(valid)
    if (.not.valid) go to 900
    read(stage,pos=1,iostat=ios) header
    valid=ios==0;call agree(valid)
    if (.not.valid) go to 900
    call same_string(header,valid)
    if (.not.valid) go to 900
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      local_start(node)=cursor
      open(newunit=patch,status='scratch',form='formatted',iostat=ios)
      valid=ios==0
      if (.not.valid) exit
      cfg=config
      if(sdc) then
        cfg%time_integrator='split';cfg%external_forcing=0;cfg%sdc_reactor='native'
      end if
      associate(g=>geometry(node),u=>state%fields(node)%values,t=>temperature%fields(node)%values(1,:,:,:))
      cfg%nx=g%nx;cfg%ny=g%ny;cfg%nz=g%nz
      cfg%x_lower=g%x_lower;cfg%x_upper=g%x_upper
      cfg%y_lower=g%y_lower;cfg%y_upper=g%y_upper
      cfg%z_lower=g%z_lower;cfg%z_upper=g%z_upper
      call reactive_eb_integrals_3d(u,g,integral,valid,l1)
      if (valid) then
        if (present(element_atoms)) then
          call reactive_eb_composition_integrals_3d(species,u,g,element_atoms,elements,valid)
        else
          call reactive_eb_element_integrals_3d(species,u,g,elements,valid)
        end if
      end if
      if (valid) call write_reactive_eb_3d_checkpoint('',species,reactions,transport,cfg,g,u,t, &
        record%time,record%steps,integral,l1,elements,record%minimum_dt,0._dp,record%minimum_theta,valid,message, &
        bundle_sha256,chemistry_integrator,base_mole_fractions,patch,element_names,element_atoms,tree_patch=.true.)
      end associate
      if(valid.and.sdc.and.config%chemistry_enabled) then
        write(patch,'(a)',iostat=ios) 'SDC_REACTION_SOURCE'
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(ios==0) write(patch,'(*(es27.18e3,1x))',iostat=ios) record%reaction_source%fields(node)%values(:,i,j,k)
        end do
        end do
        end do
        valid=ios==0
      end if
      if (valid) call stage_formatted(patch,valid)
      close(patch,iostat=ios)
      valid=valid.and.ios==0
      if (.not.valid) exit
      local_bytes(node)=cursor-local_start(node)
    end do
    call agree(valid)
    if (.not.valid) go to 900
    call MPI_Allreduce(local_bytes,bytes,size(owners),MPI_OFFSET,MPI_SUM,comm,status)
    valid=status==MPI_SUCCESS
    if (.not.valid) go to 900
    ! V2/V4 lengths index each gas record plus its optional V4 reaction source.
    valid=all(bytes>0).and.size(owners)<=(huge(1)-len(header))/21
    call agree(valid)
    if (.not.valid) go to 900
    allocate(character(len=21*size(owners))::index_text,stat=ios)
    valid=ios==0;call agree(valid)
    if (.not.valid) go to 900
    do node=1,size(owners)
      write(index_line,'(i20,a)',iostat=ios) bytes(node),new_line('a')
      if (ios/=0) then
        valid=.false.;exit
      end if
      index_text(21*(node-1)+1:21*node)=index_line
    end do
    call agree(valid)
    if (.not.valid) go to 900
    header=header//index_text;header_bytes=len(header)
    total=header_bytes
    do node=1,size(owners)
      if (bytes(node)>huge(total)-total) then
        valid=.false.;exit
      end if
      offset(node)=total;total=total+bytes(node)
    end do
    if (.not.valid) go to 900
    if (total>huge(total)-len(ending)-1) then
      valid=.false.;go to 900
    end if
    call MPI_File_open(comm,trim(path),ior(MPI_MODE_WRONLY,MPI_MODE_CREATE),MPI_INFO_NULL,file,status)
    valid=status==MPI_SUCCESS
    if (.not.valid) go to 900
    opened=.true.
    call MPI_File_set_errhandler(file,MPI_ERRORS_RETURN,status)
    valid=status==MPI_SUCCESS;call agree(valid)
    if (.not.valid) go to 900
    call MPI_File_set_size(file,total+len(ending)+1,status)
    valid=status==MPI_SUCCESS;call agree(valid)
    if (.not.valid) go to 900
    if (rank==0) then
      call write_text(0_MPI_OFFSET_KIND,header,valid)
      if (valid) call write_text(total,ending//new_line('a'),valid)
    end if
    do node=1,size(owners)
      if (.not.valid) exit
      if (owners(node)/=rank) cycle
      remaining=bytes(node);source_position=local_start(node)+1;position=offset(node)
      do while(remaining>0)
        amount=int(min(remaining,int(len(chunk),MPI_OFFSET_KIND)))
        read(stage,pos=source_position,iostat=ios) chunk(:amount)
        valid=ios==0
        if (.not.valid) exit
        call write_text(position,chunk(:amount),valid)
        if (.not.valid) exit
        remaining=remaining-amount;source_position=source_position+amount;position=position+amount
      end do
    end do
    call agree(valid)
900 continue
    if (opened) then
      call MPI_File_close(file,status)
      valid=valid.and.status==MPI_SUCCESS
    end if
    close(stage,iostat=close_status)
    valid=valid.and.close_status==0
    call agree(valid)
    ok=valid
  contains
    subroutine stage_formatted(input,good)
      integer, intent(in) :: input
      logical, intent(out) :: good
      integer :: code,outcode,got
      character(len=4096) :: buffer
      good=.false.
      rewind(input,iostat=code)
      if (code/=0) return
      do
        read(input,'(a)',advance='no',size=got,iostat=code) buffer
        if (code==iostat_end) exit
        if (code/=0.and.code/=iostat_eor) return
        if (cursor>huge(cursor)-got-1) return
        if (got>0) then
          write(stage,pos=cursor+1,iostat=outcode) buffer(:got)
          if (outcode/=0) return
          cursor=cursor+got
        end if
        if (code==iostat_eor) then
          write(stage,pos=cursor+1,iostat=outcode) new_line('a')
          if (outcode/=0) return
          cursor=cursor+1
        end if
      end do
      good=.true.
    end subroutine
    subroutine write_text(at_offset,value,good)
      integer(MPI_OFFSET_KIND), intent(in) :: at_offset
      character(len=*), intent(in) :: value
      logical, intent(out) :: good
      integer :: code
      call MPI_File_write_at(file,at_offset,value,len(value),MPI_CHARACTER,io_status,code)
      good=code==MPI_SUCCESS
      if (.not.good) return
      call MPI_Get_count(io_status,MPI_CHARACTER,written,code)
      good=code==MPI_SUCCESS.and.written==len(value)
    end subroutine
    logical function valid_record() result(good)
      real(dp) :: forcing(nv)
      good=.false.;forcing=0
      if(.not.sdc.and.config%time_integrator/='split') return
      if(.not.all(ieee_is_finite(config%external_forcing))) return
      if(sdc) then
        if(config%sdc_iterations<1.or..not.allocated(record%forcing)) return
        if(size(record%forcing)/=nv) return
        if(.not.all(ieee_is_finite(record%forcing))) return
        forcing=record%forcing
        if(.not.config%chemistry_enabled.and.allocated(record%reaction_source%fields)) return
      else
        if(any(config%external_forcing/=0)) return
        if(allocated(record%forcing).or.allocated(record%reaction_source%fields)) return
      end if
      if (.not.allocated(record%initial).or..not.allocated(record%scale).or..not.allocated(record%boundary).or. &
        .not.allocated(record%wall).or..not.allocated(record%chemical)) return
      if (any([size(record%initial),size(record%scale),size(record%boundary),size(record%wall),size(record%chemical)]/=nv)) return
      if (.not.all(ieee_is_finite([record%time,record%minimum_dt,record%minimum_theta,record%initial,record%scale, &
        record%boundary,record%wall,record%chemical]))) return
      if (record%time<0._dp.or.record%steps<0.or.record%minimum_dt<0._dp.or. &
        record%minimum_theta<0._dp.or.record%minimum_theta>1._dp.or.any(record%scale<0._dp)) return
      if (any(abs(current-record%initial-record%boundary-record%wall-record%chemical-forcing)> &
        5.e-10_dp*max(1._dp,record%scale,abs(current), &
          abs(record%boundary)+abs(record%wall)+abs(record%chemical)+abs(forcing)))) return
      good=.true.
    end function valid_record
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
end module mpi_eb_tree_checkpoint_3d_mod
