module spray_mpi_amr_checkpoint_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: iso_c_binding, only: c_int,c_char,c_null_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use spray_parcel_mod, only: spray_parcel,valid_spray_parcel
  use spray_mixture_mod, only: valid_component_inventory
  use spray_mpi_exchange_mod, only: exchange_spray_parcel_bins
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use amr_patch_tree_checkpoint_3d_mod, only: amr_tree_checkpoint_3d
  use mpi_amr_patch_tree_checkpoint_3d_mod, only: read_mpi_amr_tree_checkpoint_3d,write_mpi_amr_tree_checkpoint_3d
  implicit none
  private
  public :: mpi_amr_spray_checkpoint
  interface
    function link_file(old,new) bind(C,name='link') result(status)
      import c_char,c_int
      character(kind=c_char), intent(in) :: old(*),new(*)
      integer(c_int) :: status
    end function
  end interface
contains
  ! Immutable prefix: .pending reserves it, .gas/.particles are payloads, and
  ! the prefix itself is published last with POSIX link. Reads publish atomically.
  ! Shared path/context/model tables are required; at most one million total parcels.
  subroutine mpi_amr_spray_checkpoint(comm,writing,path,context,species,tree,owners,state,temperature, &
      parcels,component_mass,time,steps,event,regridded,last_id,injected,lower,spacing,maximum_parcels,ok,periodic,boundary)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: writing
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(inout) :: tree
    integer, allocatable, intent(inout) :: owners(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(inout) :: time,injected(:)
    integer, intent(inout) :: steps,event
    logical, intent(inout) :: regridded
    integer(int64), intent(inout) :: last_id
    real(dp), intent(in) :: lower(3),spacing(3)
    integer, intent(in) :: maximum_parcels
    logical, intent(out) :: ok
    logical, optional, intent(in) :: periodic(3)
    real(dp), optional, intent(inout) :: boundary(:)
    logical :: wraps(3)
    real(dp), allocatable :: boundary_ledger(:),reference_boundary(:)
    character(len=3) :: boundary_key
    type(amr_tree_checkpoint_3d) :: snapshot
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:),ledger(:)
    integer, allocatable :: assigned(:),counts(:),offsets(:)
    integer :: rank,nranks,ierr,unit,ios,p,nc,np,nvar,info(3),date(8),low(3),high(3),header(3)
    integer(int64) :: mark,tick
    real(dp) :: at_time
    logical :: valid,agreed,exists
    character(len=256) :: token,line
    character(len=:), allocatable :: bound_context
    ok=.false.
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=allocated(owners).and.allocated(parcels).and.allocated(component_mass)
    valid=valid.and.tree%is_valid().and.len_trim(path)>0.and.len(context)>0
    call agree();if (.not.agreed) return
    nc=size(component_mass,1);nvar=5+size(species)
    valid=nc>=1.and.nc<=16.and.all(shape(component_mass)==[nc,size(parcels)]).and.size(injected)==nvar
    if (present(boundary)) valid=valid.and.size(boundary)==nvar
    call agree();if (.not.agreed) return
    wraps=.true.
    if (present(periodic)) wraps=periodic
    write(boundary_key,'(3l1)') wraps
    call same_string(boundary_key,valid);if (.not.valid) return
    header=[nc,nvar,merge(1,0,writing)]
    call MPI_Allreduce(header,low,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(header,high,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.any(low/=high)) return
    ! Verify shared strings before root performs any I/O.
    call same_string(path,valid);if (.not.valid) return
    call same_string(context,valid);if (.not.valid) return
    allocate(ledger(nvar),boundary_ledger(nvar),reference_boundary(nvar));token=''
    boundary_ledger=0
    if (present(boundary)) boundary_ledger=boundary
    if (writing) then
      valid=ieee_is_finite(time).and.all(ieee_is_finite(injected)).and.all(ieee_is_finite(boundary_ledger))
      valid=valid.and.time>=0.and.steps>=0.and.event>=0.and.last_id>=0
      call agree();if (.not.agreed) return
      info=[steps,event,merge(1,0,regridded)];at_time=time;mark=last_id;ledger=injected
      call MPI_Bcast(info,3,MPI_INTEGER,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Bcast(at_time,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Bcast(mark,1,MPI_INTEGER8,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Bcast(ledger,nvar,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      valid=all(info==[steps,event,merge(1,0,regridded)]).and.at_time==time.and.mark==last_id.and.all(ledger==injected)
      call agree();if (.not.agreed) return
      reference_boundary=boundary_ledger
      call MPI_Bcast(reference_boundary,nvar,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      valid=all(reference_boundary==boundary_ledger)
      call agree();if (.not.agreed) return
      drops=parcels;masses=component_mass
      allocate(counts(nranks),offsets(nranks));counts=0;counts(1)=size(drops);offsets=counts(1);offsets(1)=0
      call exchange_spray_parcel_bins(comm,drops,counts,offsets,1000000,valid,masses)
      if (.not.valid) return
      valid=.true.
      if (rank==0) then
        do p=1,size(drops)
          if (drops(p)%id>mark) valid=.false.
        end do
        inquire(file=path,exist=exists);valid=valid.and..not.exists
        inquire(file=path//'.gas',exist=exists);valid=valid.and..not.exists
        inquire(file=path//'.particles',exist=exists);valid=valid.and..not.exists
        if (valid) then
          open(newunit=unit,file=path//'.pending',status='new',action='write',iostat=ios)
          valid=ios==0
          if (valid) then
            call date_and_time(values=date);call system_clock(tick)
            write(token,'(8(i0,":"),i0)') date,tick
            write(unit,'(a)',iostat=ios) 'PELEF_AMR_SPRAY_V1',trim(token),'END_AMR_SPRAY_MANIFEST'
            valid=ios==0
            close(unit,iostat=ios);valid=valid.and.ios==0
          end if
        end if
      end if
      call agree();if (.not.agreed) return
      call MPI_Bcast(token,len(token),MPI_CHARACTER,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      bound_context=context//'SPRAY_GENERATION '//trim(token)//achar(10)
      if (.not.all(wraps)) bound_context=bound_context//'SPRAY_PERIODIC '//boundary_key//achar(10)
      snapshot%tree=tree;snapshot%time=time;snapshot%steps=steps;snapshot%regrids=info(3)
      snapshot%initial=injected;snapshot%boundary=boundary_ledger
      call write_mpi_amr_tree_checkpoint_3d(path//'.gas',species,bound_context,snapshot,owners,state,temperature, &
        comm,valid,indexed=.true.)
      if (.not.valid) return
      if (rank==0) then
        open(newunit=unit,file=path//'.particles',status='new',action='write',iostat=ios)
        valid=ios==0
        if (valid) then
          write(unit,'(a)',iostat=ios) 'PELEF_AMR_PARTICLES_V1',trim(token)
          if (ios==0) write(unit,*,iostat=ios) nc,size(drops),mark,info
          if (ios==0) write(unit,'(*(es27.18e3,1x))',iostat=ios) time,ledger
          do p=1,size(drops)
            if (ios/=0) exit
            write(unit,'(3(i0,1x),*(es27.18e3,1x))',iostat=ios) drops(p)%id,drops(p)%parent_id,drops(p)%material_id, &
              drops(p)%position,drops(p)%velocity,drops(p)%mass,drops(p)%temperature,drops(p)%multiplicity, &
              drops(p)%deformation,drops(p)%deformation_rate,drops(p)%rt_time,drops(p)%khrt_initial_count, &
              drops(p)%khrt_shed_volume,masses(:,p)
          end do
          if (ios==0) write(unit,'(a)',iostat=ios) 'END_AMR_PARTICLES'
          valid=ios==0;close(unit,iostat=ios);valid=valid.and.ios==0
        end if
        if (valid) valid=link_file(path//'.pending'//c_null_char,path//c_null_char)==0
      end if
      call agree();ok=agreed;return
    end if
    valid=.true.
    if (rank==0) then
      open(newunit=unit,file=path,status='old',action='read',iostat=ios)
      valid=ios==0
      if (valid) then
        read(unit,'(a)',iostat=ios) line
        valid=ios==0.and.line=='PELEF_AMR_SPRAY_V1'
        if (valid) read(unit,'(a)',iostat=ios) token
        valid=valid.and.ios==0.and.len_trim(token)>0
        if (valid) read(unit,'(a)',iostat=ios) line
        valid=valid.and.ios==0.and.line=='END_AMR_SPRAY_MANIFEST'
        if (valid) read(unit,'(a)',iostat=ios) line
        valid=valid.and.ios<0
        close(unit)
      end if
    end if
    call agree();if (.not.agreed) return
    call MPI_Bcast(token,len(token),MPI_CHARACTER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    bound_context=context//'SPRAY_GENERATION '//trim(token)//achar(10)
    if (.not.all(wraps)) bound_context=bound_context//'SPRAY_PERIODIC '//boundary_key//achar(10)
    call read_mpi_amr_tree_checkpoint_3d(path//'.gas',species,bound_context,tree%nodes(1)%shape, &
      'round_robin',0,snapshot,assigned,comm,valid)
    if (.not.valid) return
    valid=present(boundary).or.all(snapshot%boundary==0)
    call agree();if (.not.agreed) return
    allocate(drops(0),masses(nc,0));info=0;mark=0
    if (rank==0) then
      open(newunit=unit,file=path//'.particles',status='old',action='read',iostat=ios)
      valid=ios==0
      if (valid) then
        call read_particles()
        close(unit)
      end if
    end if
    call agree();if (.not.agreed) return
    call MPI_Bcast(info,3,MPI_INTEGER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Bcast(mark,1,MPI_INTEGER8,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Bcast(ledger,nvar,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call migrate_amr_tree_spray_parcels(comm,snapshot%tree,assigned,drops,lower,spacing,wraps, &
      maximum_parcels,valid,masses)
    if (.not.valid) return
    tree=snapshot%tree;call move_alloc(assigned,owners)
    call move_alloc(snapshot%state%fields,state%fields);call move_alloc(snapshot%temperature%fields,temperature%fields)
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    time=snapshot%time;steps=info(1);event=info(2);regridded=info(3)==1;last_id=mark;injected=ledger
    if (present(boundary)) boundary=snapshot%boundary
    ok=.true.
  contains
    subroutine read_particles()
      integer :: file_nc
      valid=.false.
      read(unit,'(a)',iostat=ios) line
      if (ios/=0.or.line/='PELEF_AMR_PARTICLES_V1') return
      read(unit,'(a)',iostat=ios) line
      if (ios/=0.or.line/=token) return
      read(unit,*,iostat=ios) file_nc,np,mark,info
      if (ios/=0) return
      if (file_nc/=nc.or.np<0.or.np>1000000.or.mark<0.or.any(info<0).or.info(3)>1) return
      read(unit,*,iostat=ios) at_time,ledger
      if (ios/=0) return
      if (.not.ieee_is_finite(at_time).or..not.all(ieee_is_finite(ledger))) return
      if (at_time/=snapshot%time.or.info(1)/=snapshot%steps.or.info(3)/=snapshot%regrids) return
      if (any(ledger/=snapshot%initial)) return
      deallocate(drops,masses);allocate(drops(np),masses(nc,np),stat=ios)
      if (ios/=0) return
      do p=1,np
        read(unit,*,iostat=ios) drops(p)%id,drops(p)%parent_id,drops(p)%material_id, &
          drops(p)%position,drops(p)%velocity,drops(p)%mass,drops(p)%temperature,drops(p)%multiplicity, &
          drops(p)%deformation,drops(p)%deformation_rate,drops(p)%rt_time,drops(p)%khrt_initial_count, &
          drops(p)%khrt_shed_volume,masses(:,p)
        if (ios/=0) return
        if (.not.valid_spray_parcel(drops(p)).or..not.valid_component_inventory(drops(p),masses(:,p))) return
        if (drops(p)%id>mark) return
        if (p>1) then
          if (drops(p)%id<=drops(p-1)%id) return
        end if
      end do
      read(unit,'(a)',iostat=ios) line
      if (ios/=0.or.line/='END_AMR_PARTICLES') return
      read(unit,'(a)',iostat=ios) line
      valid=ios<0
    end subroutine
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
    subroutine same_string(value,success)
      character(len=*), intent(in) :: value
      logical, intent(out) :: success
      character(len=:), allocatable :: reference
      integer :: length
      length=len(value)
      call MPI_Bcast(length,1,MPI_INTEGER,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) then
        success=.false.;return
      end if
      allocate(character(len=length)::reference);reference(:)=value
      call MPI_Bcast(reference,length,MPI_CHARACTER,0,comm,ierr)
      valid=ierr==MPI_SUCCESS.and.value==reference.and.len(value)==length
      call agree();success=agreed
    end subroutine
  end subroutine
end module
