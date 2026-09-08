module spray_checkpoint_mod
  use, intrinsic :: iso_c_binding, only: c_char,c_int,c_null_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use spray_parcel_mod
  implicit none
  private
  character(len=*), parameter :: magic='PELEF_SPRAY_3D_CHECKPOINT_V1'
  public :: write_spray_checkpoint, read_spray_checkpoint
  ! POSIX link publishes a completely closed temporary file without replacing an
  ! existing destination. This runtime is currently qualified on Linux only.
  interface
    function posix_link(old,new) bind(C,name='link') result(status)
      import c_int,c_char
      character(kind=c_char), intent(in) :: old(*),new(*)
      integer(c_int) :: status
    end function
    function posix_unlink(path) bind(C,name='unlink') result(status)
      import c_int,c_char
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int) :: status
    end function
  end interface
contains
  subroutine write_spray_checkpoint(path,context,species,state,temperature,parcels,time,steps,event, &
      reference_totals,injected_totals,ok,message)
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),time,reference_totals(:),injected_totals(:)
    type(spray_parcel), intent(in) :: parcels(:)
    integer, intent(in) :: steps,event
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), allocatable :: check(:,:,:)
    character(len=:), allocatable :: temporary
    integer :: unit,ios,close_status,i,j,k,p,dims(3)
    integer(c_int) :: status
    logical :: exists,local_ok
    ok=.false.; message='Invalid spray checkpoint state'
    if (len_trim(path)==0 .or. len_trim(context)==0 .or. len_trim(context)>16384) return
    if (size(state,1)/=reactive_nvar(size(species)) .or. size(species)<2) return
    if (size(reference_totals) /= size(state,1) .or. size(injected_totals) /= size(state,1)) return
    if (.not. all(ieee_is_finite([time,reference_totals,injected_totals]))) return
    if (time < 0.0_dp .or. steps < 0 .or. event < 0 .or. size(parcels)>1000000) return
    do p=1,size(parcels)
      if (.not. valid_spray_parcel(parcels(p))) return
      if (p>1) then
        if (parcels(p)%id<=parcels(p-1)%id) return
      end if
    end do
    dims=shape(temperature)
    if (any(dims<2) .or. any(shape(state(1,:,:,:))/=dims)) return
    allocate(check,mold=temperature)
    call recover_reactive_temperatures_3d(species,state,temperature,dims(1),dims(2),dims(3),check,local_ok)
    if (.not. local_ok) return
    if (any(abs(check-temperature)>1.0e-9_dp*max(1.0_dp,abs(temperature)))) return
    inquire(file=trim(path),exist=exists)
    if (exists) then
      message='Checkpoint destination already exists; refusing to overwrite'; return
    end if
    temporary=trim(path)//'.part'
    open(newunit=unit,file=temporary,status='new',action='write',iostat=ios)
    if (ios/=0) then
      message='Cannot create checkpoint temporary file'; return
    end if
    write(unit,'(a)',iostat=ios) magic
    if (ios==0) write(unit,'(a)',iostat=ios) trim(context)
    if (ios==0) write(unit,*,iostat=ios) size(state,1),dims,size(parcels)
    if (ios==0) write(unit,'(es26.17e3,2(1x,i0))',iostat=ios) time,steps,event
    if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) reference_totals
    if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) injected_totals
    do k=1,dims(3)
      do j=1,dims(2)
        do i=1,dims(1)
          if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) temperature(i,j,k),state(:,i,j,k)
        end do
      end do
    end do
    do p=1,size(parcels)
      if (ios/=0) exit
      write(unit,'(i0,1x,*(es26.17e3,1x))',iostat=ios) parcels(p)%id,parcels(p)%position,parcels(p)%velocity, &
        parcels(p)%mass,parcels(p)%temperature,parcels(p)%multiplicity
    end do
    if (ios==0) write(unit,'(a)',iostat=ios) 'END_PELEF_SPRAY_CHECKPOINT'
    close(unit,iostat=close_status)
    if (ios/=0 .or. close_status/=0) then
      status=posix_unlink(temporary//c_null_char)
      message='Checkpoint write failed; destination unchanged'; return
    end if
    status=posix_link(temporary//c_null_char,trim(path)//c_null_char)
    if (status/=0) then
      status=posix_unlink(temporary//c_null_char)
      message='Atomic checkpoint publication failed; destination unchanged'; return
    end if
    status=posix_unlink(temporary//c_null_char)
    ! File visibility is atomic, but durable storage after a power failure is not
    ! promised: this portable Fortran layer does not fsync the directory.
    ok=.true.; message=''
  end subroutine

  subroutine read_spray_checkpoint(path,context,species,state,temperature,parcels,time,steps,event, &
      reference_totals,injected_totals,ok,message)
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),time,reference_totals(:),injected_totals(:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    integer, intent(inout) :: steps,event
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), allocatable :: candidate(:,:,:,:),temp(:,:,:),check(:,:,:)
    type(spray_parcel), allocatable :: drops(:)
    real(dp) :: next_time,ref(size(state,1)),injected(size(state,1)),row(max(9,size(state,1)+1))
    character(len=64) :: tokens(10)
    integer :: unit,ios,i,j,k,p,dims(3),got(5),next_steps,next_event
    character(len=16384) :: line
    logical :: local_ok
    ok=.false.; message='Invalid, truncated or incompatible spray checkpoint'
    if (len_trim(context)==0 .or. len_trim(context)>len(line)) return
    if (size(state,1)/=reactive_nvar(size(species)) .or. size(species)<2) return
    if (size(reference_totals)/=size(state,1) .or. size(injected_totals)/=size(state,1)) return
    dims=shape(temperature)
    if (any(shape(state(1,:,:,:))/=dims)) return
    open(newunit=unit,file=trim(path),status='old',action='read',iostat=ios)
    if (ios/=0) return
    read(unit,'(a)',iostat=ios) line
    if (ios/=0) goto 900
    if (trim(line)/=magic) goto 900
    read(unit,'(a)',iostat=ios) line
    if (ios/=0) goto 900
    if (trim(line)/=trim(context)) goto 900
    call read_integer_row(unit,got,ios)
    if (ios/=0) goto 900
    if (any(got(1:4)/=[size(state,1),dims])) goto 900
    if (got(5)<0 .or. got(5)>1000000) goto 900
    call read_tokens(unit,tokens(1:3),ios)
    if (ios/=0) goto 900
    read(tokens(1),*,iostat=ios) next_time
    if (ios/=0) goto 900
    if (verify(trim(tokens(2)), '0123456789+-')/=0 .or. verify(trim(tokens(3)), '0123456789+-')/=0) goto 900
    read(tokens(2),*,iostat=ios) next_steps
    if (ios/=0) goto 900
    read(tokens(3),*,iostat=ios) next_event
    if (ios/=0) goto 900
    if (.not. ieee_is_finite(next_time)) goto 900
    if (next_time<0.0_dp .or. next_steps<0 .or. next_event<0) goto 900
    call read_real_row(unit,ref,ios)
    if (ios/=0) goto 900
    call read_real_row(unit,injected,ios)
    if (ios/=0) goto 900
    if (.not. all(ieee_is_finite([ref,injected]))) goto 900
    allocate(candidate,mold=state); allocate(temp,mold=temperature); allocate(check,mold=temperature)
    allocate(drops(got(5)))
    do k=1,dims(3)
      do j=1,dims(2)
        do i=1,dims(1)
          call read_real_row(unit,row(1:size(state,1)+1),ios)
          if (ios/=0) goto 900
          temp(i,j,k)=row(1); candidate(:,i,j,k)=row(2:size(state,1)+1)
        end do
      end do
    end do
    do p=1,size(drops)
      call read_tokens(unit,tokens,ios)
      if (ios/=0) goto 900
      if (verify(trim(tokens(1)), '0123456789+-')/=0) goto 900
      read(tokens(1),*,iostat=ios) drops(p)%id
      if (ios/=0) goto 900
      do i=1,9
        read(tokens(i+1),*,iostat=ios) row(i)
        if (ios/=0) goto 900
      end do
      drops(p)%position=row(1:3); drops(p)%velocity=row(4:6)
      drops(p)%mass=row(7); drops(p)%temperature=row(8); drops(p)%multiplicity=row(9)
      if (ios/=0) goto 900
      if (.not. valid_spray_parcel(drops(p))) goto 900
      if (p>1) then
        if (drops(p)%id<=drops(p-1)%id) goto 900
      end if
    end do
    read(unit,'(a)',iostat=ios) line
    if (ios/=0) goto 900
    if (trim(line)/='END_PELEF_SPRAY_CHECKPOINT') goto 900
    read(unit,'(a)',iostat=ios) line
    if (ios>=0) goto 900
    close(unit)
    if (.not. all(ieee_is_finite(temp))) return
    call recover_reactive_temperatures_3d(species,candidate,temp,dims(1),dims(2),dims(3),check,local_ok)
    if (.not. local_ok) return
    if (any(abs(check-temp)>1.0e-9_dp*max(1.0_dp,abs(temp)))) return
    state=candidate; temperature=temp; parcels=drops; time=next_time; steps=next_steps; event=next_event
    reference_totals=ref; injected_totals=injected; ok=.true.; message=''
    return
900 continue
    close(unit)
  end subroutine
  ! Read one complete numeric record, rejecting null values, repeat syntax,
  ! trailing tokens and NaN/Inf. List-directed reads alone can silently accept
  ! "/" or missing values and leave uninitialized storage behind.
  subroutine read_tokens(unit,tokens,status)
    integer, intent(in) :: unit
    character(len=*), intent(out) :: tokens(:)
    integer, intent(out) :: status
    character(len=16384) :: record
    integer :: n,start,finish,length
    tokens=''; status=1
    read(unit,'(a)',iostat=status) record
    if (status/=0) return
    status=1; length=len_trim(record)
    if (length>=len(record)) return
    finish=1
    do n=1,size(tokens)
      do while(finish<=length)
        if (record(finish:finish)/=' ' .and. record(finish:finish)/=achar(9)) exit
        finish=finish+1
      end do
      start=finish
      do while(finish<=length)
        if (record(finish:finish)==' ' .or. record(finish:finish)==achar(9)) exit
        finish=finish+1
      end do
      if (finish==start .or. finish-start>len(tokens)) return
      tokens(n)=record(start:finish-1)
      if (verify(trim(tokens(n)),'0123456789+-.eEdD')/=0) return
    end do
    if (finish<=length) then
      if (verify(record(finish:length),' '//achar(9))/=0) return
    end if
    status=0
  end subroutine
  subroutine read_real_row(unit,values,status)
    integer, intent(in) :: unit
    real(dp), intent(out) :: values(:)
    integer, intent(out) :: status
    character(len=64) :: tokens(size(values))
    integer :: n
    values=ieee_value(0.0_dp,ieee_quiet_nan)
    call read_tokens(unit,tokens,status)
    if (status/=0) return
    do n=1,size(values)
      read(tokens(n),*,iostat=status) values(n)
      if (status/=0) return
    end do
    if (.not. all(ieee_is_finite(values))) status=1
  end subroutine
  subroutine read_integer_row(unit,values,status)
    integer, intent(in) :: unit
    integer, intent(out) :: values(:)
    integer, intent(out) :: status
    character(len=64) :: tokens(size(values))
    integer :: n
    values=-1
    call read_tokens(unit,tokens,status)
    if (status/=0) return
    do n=1,size(values)
      if (verify(trim(tokens(n)),'0123456789+-')/=0) then
        status=1; return
      end if
      read(tokens(n),*,iostat=status) values(n)
      if (status/=0) return
    end do
  end subroutine
end module spray_checkpoint_mod
