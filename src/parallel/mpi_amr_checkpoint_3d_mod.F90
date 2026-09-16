module mpi_amr_checkpoint_3d_mod
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_hierarchy_3d_mod, only: amr_patch_3d
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use simulation_config_amr_reactive_3d_mod, only: amr_reactive_3d_config
  use mpi_amr_sparse_reactive_3d_mod, only: mpi_amr_sparse_distribution_3d,mpi_amr_sparse_hierarchy_3d, &
    initialize_mpi_amr_sparse_distribution_3d
  use amr_reactive_3d_checkpoint_mod, only: amr_checkpoint_envelope,checkpoint_temperature_is_consistent, &
    read_amr_checkpoint_metadata
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use simulation_config_reactive_3d_mod, only: resolved_reactive_boundaries_3d
  implicit none
  private
  public :: write_mpi_amr_checkpoint_3d,read_mpi_amr_checkpoint_3d
contains
  subroutine write_mpi_amr_checkpoint_3d(distribution,hierarchy,path,species,flow_config,amr_config,patch, &
      time,steps,initial_integrals,maximum_reflux,ok,message,bundle_sha256,chemistry_integrator, &
      base_mole_fractions,reactions,transport,maximum_transport_diffusivity,minimum_transport_theta, &
      regrid_count,boundary_change)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(mpi_amr_sparse_hierarchy_3d), intent(in) :: hierarchy
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_3d_config), intent(in) :: flow_config
    type(amr_reactive_3d_config), intent(in) :: amr_config
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: time,initial_integrals(:),maximum_reflux
    integer, intent(in) :: steps
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    character(len=*), optional, intent(in) :: bundle_sha256,chemistry_integrator
    real(dp), optional, intent(in) :: base_mole_fractions(:),maximum_transport_diffusivity,minimum_transport_theta
    real(dp), optional, intent(in) :: boundary_change(:)
    integer, optional, intent(in) :: regrid_count
    type(elementary_reaction), optional, intent(in) :: reactions(:)
    type(gas_transport_species), optional, intent(in) :: transport(:)
    character(len=:), allocatable :: header,middle,footer
    type(MPI_File) :: file
    type(MPI_Status) :: file_status
    integer :: ierr,close_error,nvar,row_bytes
    integer, allocatable :: root_counts(:)
    integer(kind=MPI_OFFSET_KIND) :: coarse_offset,middle_offset,fine_offset,footer_offset,total_bytes,limit
    logical :: valid,all_ok
    ok=.false.;message='Invalid distributed AMR checkpoint'
    ! Each value has width 27; separators plus LF make 28 bytes per value.
    nvar=reactive_nvar(size(species));row_bytes=28*(nvar+1)
    valid=hierarchy%is_valid(distribution,patch).and.hierarchy%nvar==nvar.and.len_trim(path)>0
    call agree(valid,all_ok)
    if (.not.all_ok) return
    call amr_checkpoint_envelope(header,middle,footer,species,flow_config,amr_config,patch, &
      time,steps,initial_integrals,maximum_reflux,valid,message,bundle_sha256,chemistry_integrator, &
      base_mole_fractions,reactions,transport,maximum_transport_diffusivity,minimum_transport_theta,regrid_count,boundary_change)
    call agree(valid,all_ok)
    if (.not.all_ok) return
    message='Checkpoint metadata differs across ranks'
    call same_text(trim(path),all_ok)
    if (.not.all_ok) return
    call same_text(header,all_ok)
    if (.not.all_ok) return
    call same_text(middle,all_ok)
    if (.not.all_ok) return
    call same_text(footer,all_ok)
    if (.not.all_ok) return
    allocate(root_counts(distribution%nranks));root_counts=distribution%coarse_counts
    call MPI_Bcast(root_counts,size(root_counts),MPI_INTEGER,0,distribution%comm,ierr)
    call agree(ierr==MPI_SUCCESS,all_ok)
    if (.not.all_ok) return
    call agree(all(root_counts==distribution%coarse_counts),all_ok)
    if (.not.all_ok) return
    call validate_checkpoint_fields(distribution,hierarchy,species,patch,all_ok,message)
    if (.not.all_ok) return
    ! Guard MPI default-integer row counts and wide file offsets before opening.
    valid=int(max(distribution%coarse_count,distribution%fine_count),int64)*row_bytes<=huge(1)
    limit=huge(total_bytes)-int(len(header),MPI_OFFSET_KIND)-int(len(middle),MPI_OFFSET_KIND)-len(footer)
    coarse_offset=int(len(header),MPI_OFFSET_KIND)
    middle_offset=level_bytes([patch%coarse_nx,patch%coarse_ny,patch%coarse_nz],limit)
    valid=valid.and.middle_offset>=0
    fine_offset=level_bytes([patch%fine_nx(),patch%fine_ny(),patch%fine_nz()],max(0_MPI_OFFSET_KIND,limit-middle_offset))
    valid=valid.and.fine_offset>=0
    call agree(valid,all_ok)
    if (.not.all_ok) return
    footer_offset=coarse_offset+middle_offset+len(middle)+fine_offset
    middle_offset=coarse_offset+middle_offset
    fine_offset=middle_offset+len(middle)
    total_bytes=footer_offset+len(footer)
    message='Failed distributed checkpoint file write'
    call MPI_File_open(distribution%comm,trim(path),ior(MPI_MODE_CREATE,MPI_MODE_WRONLY),MPI_INFO_NULL,file,ierr)
    call agree(ierr==MPI_SUCCESS,all_ok)
    if (.not.all_ok) return
    call MPI_File_set_errhandler(file,MPI_ERRORS_RETURN,ierr)
    call agree(ierr==MPI_SUCCESS,all_ok)
    if (.not.all_ok) goto 900
    call MPI_File_set_size(file,total_bytes,ierr)
    call agree(ierr==MPI_SUCCESS,all_ok)
    if (.not.all_ok) goto 900
    call write_text(header,0_MPI_OFFSET_KIND,all_ok)
    if (.not.all_ok) goto 900
    call write_text(middle,middle_offset,all_ok)
    if (.not.all_ok) goto 900
    call write_text(footer,footer_offset,all_ok)
    if (.not.all_ok) goto 900
    call write_level(hierarchy%coarse_state,hierarchy%coarse_temperature,distribution%coarse_first, &
      patch%coarse_nx,coarse_offset,all_ok)
    if (.not.all_ok) goto 900
    call write_level(hierarchy%fine_state,hierarchy%fine_temperature,distribution%fine_first, &
      patch%fine_nx(),fine_offset,all_ok)
    if (.not.all_ok) goto 900
    ok=.true.
900 call MPI_File_close(file,close_error)
    call agree(close_error==MPI_SUCCESS,all_ok)
    ok=ok.and.all_ok
    if (ok) message=''
  contains
    subroutine agree(value,result)
      logical, intent(in) :: value
      logical, intent(out) :: result
      integer :: status
      call MPI_Allreduce(value,result,1,MPI_LOGICAL,MPI_LAND,distribution%comm,status)
      result=result.and.status==MPI_SUCCESS
    end subroutine
    subroutine same_text(value,result)
      character(len=*), intent(in) :: value
      logical, intent(out) :: result
      character(len=:), allocatable :: reference
      integer :: count,status
      count=len(value)
      call MPI_Bcast(count,1,MPI_INTEGER,0,distribution%comm,status)
      call agree(status==MPI_SUCCESS.and.count==len(value),result)
      if (.not.result) return
      reference=value
      call MPI_Bcast(reference,count,MPI_CHARACTER,0,distribution%comm,status)
      call agree(status==MPI_SUCCESS.and.reference==value,result)
    end subroutine
    function level_bytes(shape,maximum) result(bytes)
      integer, intent(in) :: shape(3)
      integer(kind=MPI_OFFSET_KIND), intent(in) :: maximum
      integer(kind=MPI_OFFSET_KIND) :: bytes
      integer :: axis
      bytes=row_bytes
      do axis=1,3
        if (bytes>maximum/int(shape(axis),MPI_OFFSET_KIND)) then
          bytes=-1;return
        end if
        bytes=bytes*shape(axis)
      end do
    end function
    subroutine write_text(value,offset,result)
      character(len=*), intent(in) :: value
      integer(kind=MPI_OFFSET_KIND), intent(in) :: offset
      logical, intent(out) :: result
      integer :: status
      status=MPI_SUCCESS
      if (distribution%rank==0) call MPI_File_write_at(file,offset,value,len(value),MPI_CHARACTER,file_status,status)
      call agree(status==MPI_SUCCESS,result)
    end subroutine
    subroutine write_level(u,t,first,nx,base,result)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:)
      integer, intent(in) :: first,nx
      integer(kind=MPI_OFFSET_KIND), intent(in) :: base
      logical, intent(out) :: result
      character(len=:), allocatable :: block
      character(len=28*(nvar+1)-1) :: row
      integer :: i,j,k,position,count,status
      integer(kind=MPI_OFFSET_KIND) :: offset
      allocate(character(len=max(1,row_bytes*size(t,1)))::block)
      count=row_bytes*size(t,1);result=.false.
      do k=1,size(t,3)
        do j=1,size(t,2)
          position=0;status=0
          do i=1,size(t,1)
            write(row,'(*(es27.18e3,1x))',iostat=status) u(:,i,j,k),t(i,j,k)
            if (status/=0) exit
            block(position+1:position+row_bytes)=row//new_line('a');position=position+row_bytes
          end do
          call agree(status==0,result)
          if (.not.result) return
          offset=base+((int(k-1,MPI_OFFSET_KIND)*size(t,2)+j-1)*nx+max(1,first)-1)*row_bytes
          call MPI_File_write_at_all(file,offset,block,count,MPI_CHARACTER,file_status,status)
          call agree(status==MPI_SUCCESS,result)
          if (.not.result) return
        end do
      end do
    end subroutine
  end subroutine write_mpi_amr_checkpoint_3d
  subroutine read_mpi_amr_checkpoint_3d(distribution,hierarchy,path,species,flow_config,amr_config,patch, &
      time,steps,initial_integrals,maximum_reflux,ok,message,bundle_sha256,chemistry_integrator, &
      base_mole_fractions,reactions,transport,maximum_transport_diffusivity,minimum_transport_theta, &
      regrid_count,boundary_change)
    type(mpi_amr_sparse_distribution_3d), intent(inout) :: distribution
    type(mpi_amr_sparse_hierarchy_3d), intent(inout) :: hierarchy
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_3d_config), intent(in) :: flow_config
    type(amr_reactive_3d_config), intent(in) :: amr_config
    type(amr_patch_3d), intent(inout) :: patch
    real(dp), intent(inout) :: time,initial_integrals(:),maximum_reflux
    integer, intent(inout) :: steps
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    character(len=*), optional, intent(in) :: bundle_sha256,chemistry_integrator
    real(dp), optional, intent(in) :: base_mole_fractions(:)
    real(dp), optional, intent(inout) :: maximum_transport_diffusivity,minimum_transport_theta,boundary_change(:)
    integer, optional, intent(inout) :: regrid_count
    type(elementary_reaction), optional, intent(in) :: reactions(:)
    type(gas_transport_species), optional, intent(in) :: transport(:)
    type(mpi_amr_sparse_distribution_3d) :: candidate_distribution
    type(mpi_amr_sparse_hierarchy_3d) :: candidate
    type(amr_patch_3d) :: candidate_patch
    real(dp) :: candidate_time,candidate_reflux,diffusivity,theta
    real(dp) :: integrals(reactive_nvar(size(species))),boundary(reactive_nvar(size(species)))
    real(dp), allocatable :: recovered(:,:,:)
    integer :: candidate_steps,count,row_bytes,nvar,ierr,close_error
    integer, allocatable :: root_counts(:)
    integer(int64) :: coarse_offset,fine_offset
    integer(int64) :: root_offsets(2)
    character(len=:), allocatable :: header,middle,footer
    type(MPI_File) :: file
    type(MPI_Status) :: file_status
    logical :: valid,all_ok,physical,scan_mode
    ok=.false.;message='Invalid distributed restart arguments'
    nvar=reactive_nvar(size(species))
    physical=any(resolved_reactive_boundaries_3d(flow_config)/='periodic')
    valid=distribution%is_valid(patch).and.size(initial_integrals)==nvar.and.len_trim(path)>0
    if (physical) then
      valid=valid.and.present(boundary_change)
      if (present(boundary_change)) valid=valid.and.size(boundary_change)==nvar
    end if
    if (flow_config%transport_enabled) then
      valid=valid.and.present(transport).and.present(maximum_transport_diffusivity).and.present(minimum_transport_theta)
    else
      valid=valid.and..not.(present(transport).or.present(maximum_transport_diffusivity).or.present(minimum_transport_theta))
    end if
    call agree(valid,all_ok)
    if (.not.all_ok) return
    call same_text(trim(path),all_ok)
    if (.not.all_ok) return
    scan_mode=.false.
100 candidate=mpi_amr_sparse_hierarchy_3d()
    if (allocated(root_counts)) deallocate(root_counts)
    candidate_patch=patch;candidate_time=0;candidate_reflux=0;candidate_steps=0;count=0
    integrals=0;boundary=0;diffusivity=0;theta=1
    if (flow_config%transport_enabled) then
      call read_amr_checkpoint_metadata(path,species,flow_config,amr_config,candidate_patch, &
        coarse_offset,fine_offset,row_bytes,candidate_time,candidate_steps,integrals,candidate_reflux,valid,message, &
        bundle_sha256,chemistry_integrator,base_mole_fractions,reactions,transport,diffusivity,theta,count,boundary,scan_mode)
    else
      call read_amr_checkpoint_metadata(path,species,flow_config,amr_config,candidate_patch, &
        coarse_offset,fine_offset,row_bytes,candidate_time,candidate_steps,integrals,candidate_reflux,valid,message, &
        bundle_sha256,chemistry_integrator,base_mole_fractions,reactions,regrid_count=count,boundary_change=boundary, &
        scan_rows=scan_mode)
    end if
    call agree(valid,all_ok)
    if (.not.all_ok) then
      if (.not.scan_mode) then
        scan_mode=.true.;goto 100
      end if
      return
    end if
    ! Compare the complete validated context and stored diagnostics across ranks.
    if (flow_config%transport_enabled) then
      call amr_checkpoint_envelope(header,middle,footer,species,flow_config,amr_config,candidate_patch, &
        candidate_time,candidate_steps,integrals,candidate_reflux,valid,message,bundle_sha256,chemistry_integrator, &
        base_mole_fractions,reactions,transport,diffusivity,theta,count,boundary)
    else
      call amr_checkpoint_envelope(header,middle,footer,species,flow_config,amr_config,candidate_patch, &
        candidate_time,candidate_steps,integrals,candidate_reflux,valid,message,bundle_sha256,chemistry_integrator, &
        base_mole_fractions,reactions,regrid_count=count,boundary_change=boundary)
    end if
    call agree(valid,all_ok)
    if (.not.all_ok) return
    message='Distributed restart metadata differs across ranks'
    call same_text(header,all_ok)
    if (.not.all_ok) return
    call same_text(middle,all_ok)
    if (.not.all_ok) return
    call same_text(footer,all_ok)
    if (.not.all_ok) return
    root_offsets=[coarse_offset,fine_offset]
    call MPI_Bcast(root_offsets,2,MPI_INTEGER8,0,distribution%comm,ierr)
    call agree(ierr==MPI_SUCCESS.and.all(root_offsets==[coarse_offset,fine_offset]),all_ok)
    if (.not.all_ok) return
    allocate(root_counts(distribution%nranks));root_counts=distribution%coarse_counts
    call MPI_Bcast(root_counts,size(root_counts),MPI_INTEGER,0,distribution%comm,ierr)
    call agree(ierr==MPI_SUCCESS.and.all(root_counts==distribution%coarse_counts),all_ok)
    if (.not.all_ok) return
    call initialize_mpi_amr_sparse_distribution_3d(candidate_patch,distribution%comm,candidate_distribution,valid,root_counts)
    call agree(valid,all_ok)
    if (.not.all_ok) return
    candidate%nvar=nvar
    allocate(candidate%coarse_state(nvar,candidate_distribution%coarse_count,candidate_patch%coarse_ny, &
      candidate_patch%coarse_nz))
    allocate(candidate%coarse_temperature(candidate_distribution%coarse_count,candidate_patch%coarse_ny, &
      candidate_patch%coarse_nz))
    allocate(candidate%fine_state(nvar,candidate_distribution%fine_count,candidate_patch%fine_ny(),candidate_patch%fine_nz()))
    allocate(candidate%fine_temperature(candidate_distribution%fine_count,candidate_patch%fine_ny(),candidate_patch%fine_nz()))
    valid=int(max(candidate_distribution%coarse_count,candidate_distribution%fine_count),int64)*row_bytes<=huge(1)
    call agree(valid,all_ok)
    if (.not.all_ok) return
    if (scan_mode) then
      call scan_level(candidate%coarse_state,candidate%coarse_temperature,candidate_distribution%coarse_first, &
        candidate_patch%coarse_nx,coarse_offset,valid)
      call agree(valid,all_ok)
      if (.not.all_ok) return
      call scan_level(candidate%fine_state,candidate%fine_temperature,candidate_distribution%fine_first, &
        candidate_patch%fine_nx(),fine_offset,valid)
      call agree(valid,all_ok)
      if (.not.all_ok) return
      goto 910
    end if
    message='Failed distributed checkpoint data read'
    call MPI_File_open(distribution%comm,trim(path),MPI_MODE_RDONLY,MPI_INFO_NULL,file,ierr)
    call agree(ierr==MPI_SUCCESS,all_ok)
    if (.not.all_ok) return
    call MPI_File_set_errhandler(file,MPI_ERRORS_RETURN,ierr)
    call agree(ierr==MPI_SUCCESS,all_ok)
    if (.not.all_ok) goto 900
    call read_level(candidate%coarse_state,candidate%coarse_temperature,candidate_distribution%coarse_first, &
      candidate_patch%coarse_nx,coarse_offset,all_ok)
    if (.not.all_ok) goto 900
    call read_level(candidate%fine_state,candidate%fine_temperature,candidate_distribution%fine_first, &
      candidate_patch%fine_nx(),fine_offset,all_ok)
900 call MPI_File_close(file,close_error)
    call agree(close_error==MPI_SUCCESS.and.all_ok,valid)
    if (.not.valid) then
      scan_mode=.true.;goto 100
    end if
910 call validate_checkpoint_fields(candidate_distribution,candidate,species,candidate_patch,valid,message)
    if (.not.valid) return
    ! Match the serial reader's recovered-temperature commit, not just its tolerance check.
    allocate(recovered,mold=candidate%coarse_temperature)
    call recover_reactive_temperatures_3d(species,candidate%coarse_state,candidate%coarse_temperature, &
      candidate_distribution%coarse_count,candidate_patch%coarse_ny,candidate_patch%coarse_nz,recovered,valid)
    call agree(valid,all_ok)
    if (.not.all_ok) return
    candidate%coarse_temperature=recovered;deallocate(recovered)
    if (candidate_distribution%fine_count>0) then
      allocate(recovered,mold=candidate%fine_temperature)
      call recover_reactive_temperatures_3d(species,candidate%fine_state,candidate%fine_temperature, &
        candidate_distribution%fine_count,candidate_patch%fine_ny(),candidate_patch%fine_nz(),recovered,valid)
      if (valid) candidate%fine_temperature=recovered
    end if
    call agree(valid,all_ok)
    if (.not.all_ok) return
    hierarchy=candidate;distribution=candidate_distribution;patch=candidate_patch
    time=candidate_time;steps=candidate_steps;initial_integrals=integrals;maximum_reflux=candidate_reflux
    if (present(regrid_count)) regrid_count=count
    if (physical) boundary_change=boundary
    if (flow_config%transport_enabled) then
      maximum_transport_diffusivity=diffusivity;minimum_transport_theta=theta
    end if
    ok=.true.;message=''
  contains
    subroutine agree(value,result)
      logical, intent(in) :: value
      logical, intent(out) :: result
      integer :: status
      call MPI_Allreduce(value,result,1,MPI_LOGICAL,MPI_LAND,distribution%comm,status)
      result=result.and.status==MPI_SUCCESS
    end subroutine
    subroutine same_text(value,result)
      character(len=*), intent(in) :: value
      logical, intent(out) :: result
      character(len=:), allocatable :: reference
      integer :: length,status
      length=len(value)
      call MPI_Bcast(length,1,MPI_INTEGER,0,distribution%comm,status)
      call agree(status==MPI_SUCCESS.and.length==len(value),result)
      if (.not.result) return
      reference=value
      call MPI_Bcast(reference,length,MPI_CHARACTER,0,distribution%comm,status)
      call agree(status==MPI_SUCCESS.and.reference==value,result)
    end subroutine
    subroutine scan_level(u,t,first,nx,base,result)
      real(dp), intent(out) :: u(:,:,:,:),t(:,:,:)
      integer, intent(in) :: first,nx
      integer(int64), intent(in) :: base
      logical, intent(out) :: result
      real(dp) :: cell(nvar+1)
      integer :: unit,status,close_status,i,j,k
      result=.false.
      open(newunit=unit,file=path,status='old',action='read',form='formatted',access='stream',iostat=status)
      if (status/=0) return
      do k=1,size(t,3)
        do j=1,size(t,2)
          do i=1,nx
            cell=ieee_value(0.0_dp,ieee_quiet_nan)
            if (i==1.and.j==1.and.k==1) then
              read(unit,*,pos=base+1,iostat=status) cell
            else
              read(unit,*,iostat=status) cell
            end if
            if (status/=0) goto 920
            if (i>=first.and.i<first+size(t,1)) then
              u(:,i-first+1,j,k)=cell(:nvar);t(i-first+1,j,k)=cell(nvar+1)
            end if
          end do
        end do
      end do
      result=.true.
920   close(unit,iostat=close_status)
      result=result.and.close_status==0
    end subroutine scan_level
    subroutine read_level(u,t,first,nx,base,result)
      real(dp), intent(out) :: u(:,:,:,:),t(:,:,:)
      integer, intent(in) :: first,nx
      integer(int64), intent(in) :: base
      logical, intent(out) :: result
      character(len=:), allocatable :: block
      integer :: i,j,k,position,bytes,got,status
      integer(kind=MPI_OFFSET_KIND) :: offset
      logical :: parsed
      allocate(character(len=max(1,row_bytes*size(t,1)))::block)
      bytes=row_bytes*size(t,1);result=.false.
      u=ieee_value(0.0_dp,ieee_quiet_nan);t=ieee_value(0.0_dp,ieee_quiet_nan)
      do k=1,size(t,3)
        do j=1,size(t,2)
          offset=int(base,MPI_OFFSET_KIND)+ &
            ((int(k-1,MPI_OFFSET_KIND)*size(t,2)+j-1)*nx+max(1,first)-1)*row_bytes
          call MPI_File_read_at_all(file,offset,block,bytes,MPI_CHARACTER,file_status,status)
          call agree(status==MPI_SUCCESS,result)
          if (.not.result) return
          call MPI_Get_count(file_status,MPI_CHARACTER,got,status)
          call agree(status==MPI_SUCCESS.and.got==bytes,result)
          if (.not.result) return
          position=0;parsed=.true.
          do i=1,size(t,1)
            if (block(position+row_bytes:position+row_bytes)/=new_line('a')) then
              parsed=.false.;exit
            end if
            read(block(position+1:position+row_bytes-1),*,iostat=status) u(:,i,j,k),t(i,j,k)
            if (status/=0) then
              parsed=.false.;exit
            end if
            position=position+row_bytes
          end do
          call agree(parsed,result)
          if (.not.result) return
        end do
      end do
    end subroutine
  end subroutine read_mpi_amr_checkpoint_3d

  subroutine validate_checkpoint_fields(distribution,hierarchy,species,patch,ok,message)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(mpi_amr_sparse_hierarchy_3d), intent(in) :: hierarchy
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    logical :: valid,all_ok
    integer :: ierr,ci,cj,ck,gi,fi,fj,fk,ii,jj,kk,r
    real(dp) :: local_max(2),global_max(2),weight,restricted(reactive_nvar(size(species)))
    ok=.false.
    message='Invalid distributed checkpoint state or temperature'
    valid=all(ieee_is_finite(hierarchy%coarse_state)).and.all(ieee_is_finite(hierarchy%fine_state)).and. &
      all(ieee_is_finite(hierarchy%coarse_temperature)).and.all(ieee_is_finite(hierarchy%fine_temperature))
    call agree(valid,all_ok)
    if (.not.all_ok) return
    valid=all(hierarchy%coarse_temperature>0).and.all(hierarchy%fine_temperature>0)
    call agree(valid,all_ok)
    if (.not.all_ok) return
    call checkpoint_temperature_is_consistent(species,hierarchy%coarse_state,hierarchy%coarse_temperature,valid)
    call agree(valid,all_ok)
    if (.not.all_ok) return
    if (distribution%fine_count>0) &
      call checkpoint_temperature_is_consistent(species,hierarchy%fine_state,hierarchy%fine_temperature,valid)
    call agree(valid,all_ok)
    if (.not.all_ok) return
    r=patch%refinement_ratio;weight=1.0_dp/real(r**3,dp);local_max=0
    do ck=patch%coarse_k_lower,patch%coarse_k_upper
      fk=(ck-patch%coarse_k_lower)*r+1
      do cj=patch%coarse_j_lower,patch%coarse_j_upper
        fj=(cj-patch%coarse_j_lower)*r+1
        do gi=max(distribution%coarse_first,patch%coarse_i_lower),min(distribution%coarse_last,patch%coarse_i_upper)
          ci=distribution%coarse_local_index(gi);fi=(gi-patch%coarse_i_lower)*r+2-distribution%fine_first
          restricted=0
          do kk=fk,fk+r-1
            do jj=fj,fj+r-1
              do ii=fi,fi+r-1
                restricted=restricted+weight*hierarchy%fine_state(:,ii,jj,kk)
              end do
            end do
          end do
          local_max(1)=max(local_max(1),maxval(abs(restricted-hierarchy%coarse_state(:,ci,cj,ck))))
          local_max(2)=max(local_max(2),maxval(abs(restricted)))
        end do
      end do
    end do
    call agree(all(ieee_is_finite(local_max)),all_ok)
    if (.not.all_ok) return
    call MPI_Allreduce(local_max,global_max,2,MPI_DOUBLE_PRECISION,MPI_MAX,distribution%comm,ierr)
    call agree(ierr==MPI_SUCCESS,all_ok)
    if (.not.all_ok) return
    message='Distributed checkpoint levels are not synchronized'
    if (global_max(1)/max(1.0_dp,global_max(2))>64*epsilon(1.0_dp)) return
    ok=.true.;message=''
  contains
    subroutine agree(value,result)
      logical, intent(in) :: value
      logical, intent(out) :: result
      integer :: status
      call MPI_Allreduce(value,result,1,MPI_LOGICAL,MPI_LAND,distribution%comm,status)
      result=result.and.status==MPI_SUCCESS
    end subroutine
  end subroutine validate_checkpoint_fields
end module mpi_amr_checkpoint_3d_mod
