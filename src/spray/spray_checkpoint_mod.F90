module spray_checkpoint_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: iso_c_binding, only: c_char,c_int,c_null_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use spray_wall_impact_mod, only: spray_wall_film,valid_wall_film,valid_wall_component_inventory
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory
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
  logical function parent_exists(parcels,index) result(found)
    type(spray_parcel), intent(in) :: parcels(:)
    integer, intent(in) :: index
    integer :: lo,hi,mid
    found=.false.;lo=1;hi=index-1
    do while(lo<=hi)
      mid=lo+(hi-lo)/2
      if (parcels(mid)%id==parcels(index)%parent_id) then
        found=.true.;return
      else if (parcels(mid)%id<parcels(index)%parent_id) then
        lo=mid+1
      else
        hi=mid-1
      end if
    end do
  end function
  function spray_checkpoint_magic(physical,breakup,split,khrt) result(value)
    logical, intent(in) :: physical,breakup,split,khrt
    character(len=40) :: value
    if (khrt) then
      if (split) then
        if (physical) then
          value='PELEF_SPRAY_3D_CHECKPOINT_V10'
        else
          value='PELEF_SPRAY_3D_CHECKPOINT_V9'
        end if
      else
        if (physical) then
          value='PELEF_SPRAY_3D_CHECKPOINT_V8'
        else
          value='PELEF_SPRAY_3D_CHECKPOINT_V7'
        end if
      end if
    else if (split) then
      if (physical) then
        value='PELEF_SPRAY_3D_CHECKPOINT_V6'
      else
        value='PELEF_SPRAY_3D_CHECKPOINT_V5'
      end if
    else if (breakup) then
      if (physical) then
        value='PELEF_SPRAY_3D_CHECKPOINT_V4'
      else
        value='PELEF_SPRAY_3D_CHECKPOINT_V3'
      end if
    else
      if (physical) then
        value='PELEF_SPRAY_3D_CHECKPOINT_V2'
      else
        value=magic
      end if
    end if
  end function
  subroutine write_spray_checkpoint(path,context,species,state,temperature,parcels,time,steps,event, &
      reference_totals,injected_totals,ok,message,boundary_totals,breakup_state,split_state,khrt_state, &
      source_totals,material_state,film_state,wall_films,shedding_state,component_mass,film_component_mass,last_id)
    integer(int64), intent(in), optional :: last_id
    logical, intent(in), optional :: breakup_state,split_state,khrt_state,material_state,film_state,shedding_state
    logical :: has_breakup,has_split,has_khrt,has_materials,has_films,has_shedding
    character(len=80) :: expected_magic
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),time,reference_totals(:),injected_totals(:)
    type(spray_parcel), intent(in) :: parcels(:)
    type(spray_wall_film), intent(in), optional :: wall_films(:)
    integer, intent(in) :: steps,event
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), intent(in), optional :: boundary_totals(:),source_totals(:)
    real(dp), intent(in), optional :: component_mass(:,:),film_component_mass(:,:)
    real(dp), allocatable :: check(:,:,:)
    character(len=:), allocatable :: temporary
    integer :: unit,ios,close_status,i,j,k,p,dims(3)
    integer(c_int) :: status
    logical :: exists,local_ok
    has_films=.false.
    if (present(film_state)) has_films=film_state
    has_materials=.false.
    if (present(material_state)) has_materials=material_state
    has_breakup=.false.
    if (present(breakup_state)) has_breakup=breakup_state
    has_split=.false.
    if (present(split_state)) has_split=split_state
    has_shedding=.false.
    if (present(shedding_state)) has_shedding=shedding_state
    has_khrt=.false.
    if (present(khrt_state)) has_khrt=khrt_state
    expected_magic=spray_checkpoint_magic(present(boundary_totals),has_breakup,has_split,has_khrt)
    ok=.false.; message='Invalid spray checkpoint state'
    if (present(wall_films).and..not.has_films) then
      if (size(wall_films)>0) return
    end if
    if (present(source_totals)) then
      expected_magic=trim(expected_magic)//'_FORCE_V1'
      if (size(source_totals)/=size(state,1)) return
    end if
    if (has_materials) expected_magic=trim(expected_magic)//'_MAT_V1'
    if (has_films) then
      if (.not.present(wall_films)) return
      expected_magic=trim(expected_magic)//'_FILM_V1'
    end if
    if (has_shedding) then
      if (.not.has_khrt.or..not.has_split) return
      expected_magic=trim(expected_magic)//'_SHED_V1'
    end if
    if (present(component_mass)) expected_magic=trim(expected_magic)//'_MIX_V1'
    if ((has_films.and.present(component_mass)).neqv.present(film_component_mass)) return
    if (present(film_component_mass)) expected_magic=trim(expected_magic)//'_FMIX_V1'
    if (present(last_id)) expected_magic=trim(expected_magic)//'_ID_V1'
    if ((has_split.or.has_khrt).and..not.has_breakup) return
    if (present(boundary_totals)) then
      if (size(boundary_totals)/=size(state,1)) return
      if (.not.all(ieee_is_finite(boundary_totals))) return
    end if
    if (present(source_totals)) then
      if (.not.all(ieee_is_finite(source_totals))) return
    end if
    if (len_trim(path)==0 .or. len_trim(context)==0 .or. len_trim(context)>16384) return
    if (size(state,1)/=reactive_nvar(size(species)) .or. size(species)<2) return
    if (size(reference_totals) /= size(state,1) .or. size(injected_totals) /= size(state,1)) return
    if (.not. all(ieee_is_finite([time,reference_totals,injected_totals]))) return
    if (time < 0.0_dp .or. steps < 0 .or. event < 0 .or. size(parcels)>1000000) return
    if (has_films) then
      if (size(wall_films)>1000000) return
      do p=1,size(wall_films)
        if (.not.valid_wall_film(wall_films(p))) return
      end do
    end if
    if (present(component_mass)) then
      if (size(component_mass,1)<1.or.size(component_mass,1)>min(16,size(species))) return
      if (size(component_mass,2)/=size(parcels)) return
      do p=1,size(parcels)
        if (.not.valid_component_inventory(parcels(p),component_mass(:,p))) return
      end do
    end if
    if (present(film_component_mass)) then
      if (any(shape(film_component_mass)/=[size(component_mass,1),size(wall_films)])) return
      do p=1,size(wall_films)
        if (.not.valid_wall_component_inventory(wall_films(p),film_component_mass(:,p))) return
      end do
    end if
    if (present(last_id)) then
      if (last_id<0) return
      if (any(parcels%id>last_id)) return
    end if
    do p=1,size(parcels)
      if (.not. valid_spray_parcel(parcels(p))) return
      if (.not.has_materials.and.parcels(p)%material_id/=1) return
      if (.not.(has_split.or.has_films).and.parcels(p)%parent_id/=0) return
      if (.not.has_shedding.and.(parcels(p)%khrt_initial_count/=0.or.parcels(p)%khrt_shed_volume/=0)) return
      if (.not.has_khrt.and.parcels(p)%rt_time/=-1.0_dp) return
      if (parcels(p)%parent_id>0) then
        if (.not.parent_exists(parcels,p)) return
      end if
      if (.not.has_breakup) then
        if (abs(parcels(p)%deformation)>0.or.abs(parcels(p)%deformation_rate)>0) return
      end if
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
    write(unit,'(a)',iostat=ios) trim(expected_magic)
    if (ios==0) write(unit,'(a)',iostat=ios) trim(context)
    if (ios==0) write(unit,*,iostat=ios) size(state,1),dims,size(parcels)
    if (present(component_mass)) then
      if (ios==0) write(unit,'(i0)',iostat=ios) size(component_mass,1)
    end if
    if (present(last_id)) then
      if (ios==0) write(unit,'(i0)',iostat=ios) last_id
    end if
    if (ios==0) write(unit,'(es26.17e3,2(1x,i0))',iostat=ios) time,steps,event
    if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) reference_totals
    if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) injected_totals
    if (present(boundary_totals)) then
      if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) boundary_totals
    end if
    if (present(source_totals)) then
      if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) source_totals
    end if
    do k=1,dims(3)
      do j=1,dims(2)
        do i=1,dims(1)
          if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) temperature(i,j,k),state(:,i,j,k)
        end do
      end do
    end do
    do p=1,size(parcels)
      if (ios/=0) exit
      if (has_breakup) then
      write(unit,'(i0,1x,*(es26.17e3,1x))',iostat=ios) parcels(p)%id,parcels(p)%position,parcels(p)%velocity, &
        parcels(p)%mass,parcels(p)%temperature,parcels(p)%multiplicity,parcels(p)%deformation,parcels(p)%deformation_rate
      else
      write(unit,'(i0,1x,*(es26.17e3,1x))',iostat=ios) parcels(p)%id,parcels(p)%position,parcels(p)%velocity, &
        parcels(p)%mass,parcels(p)%temperature,parcels(p)%multiplicity
      end if
      if (has_materials.and.ios==0) write(unit,'(i0)',iostat=ios) parcels(p)%material_id
      if ((has_split.or.has_films).and.ios==0) write(unit,'(i0)',iostat=ios) parcels(p)%parent_id
      if (has_khrt.and.ios==0) write(unit,'(es26.17e3)',iostat=ios) parcels(p)%rt_time
      if (has_shedding.and.ios==0) write(unit,'(2es26.17e3)',iostat=ios) &
        parcels(p)%khrt_initial_count,parcels(p)%khrt_shed_volume
      if (present(component_mass)) then
        if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) component_mass(:,p)
      end if
    end do
    if (has_films.and.ios==0) then
      write(unit,'(i0)',iostat=ios) size(wall_films)
      do p=1,size(wall_films)
        if (ios/=0) exit
        write(unit,'(i0)',iostat=ios) wall_films(p)%material_id
        if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) wall_films(p)%position,wall_films(p)%normal, &
          wall_films(p)%mass,wall_films(p)%temperature,wall_films(p)%diameter,wall_films(p)%height
        if (present(film_component_mass)) then
          if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) film_component_mass(:,p)
        end if
      end do
    end if
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
      reference_totals,injected_totals,ok,message,boundary_totals,breakup_state,split_state,khrt_state, &
      source_totals,material_state,film_state,wall_films,shedding_state,component_mass,film_component_mass,last_id)
    integer(int64), intent(inout), optional :: last_id
    integer(int64) :: candidate_last_id
    logical, intent(in), optional :: breakup_state,split_state,khrt_state,material_state,film_state,shedding_state
    logical :: has_breakup,has_split,has_khrt,has_materials,has_films,has_shedding
    character(len=80) :: expected_magic
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),time,reference_totals(:),injected_totals(:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    type(spray_wall_film), allocatable, intent(inout), optional :: wall_films(:)
    type(spray_wall_film), allocatable :: patches(:)
    integer :: film_count
    integer, intent(inout) :: steps,event
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), intent(inout), optional :: boundary_totals(:),source_totals(:)
    real(dp), allocatable, intent(inout), optional :: component_mass(:,:),film_component_mass(:,:)
    real(dp), allocatable :: mass_candidate(:,:),film_mass_candidate(:,:)
    integer :: component_count(1)
    real(dp) :: source(size(state,1))
    real(dp) :: boundary(size(state,1))
    real(dp), allocatable :: candidate(:,:,:,:),temp(:,:,:),check(:,:,:)
    type(spray_parcel), allocatable :: drops(:)
    real(dp) :: next_time,ref(size(state,1)),injected(size(state,1)),row(max(11,size(state,1)+1))
    character(len=64) :: tokens(12)
    integer :: unit,ios,i,j,k,p,dims(3),got(5),next_steps,next_event
    character(len=16384) :: line
    logical :: local_ok
    has_films=.false.
    if (present(film_state)) has_films=film_state
    has_materials=.false.
    if (present(material_state)) has_materials=material_state
    has_breakup=.false.
    if (present(breakup_state)) has_breakup=breakup_state
    has_split=.false.
    if (present(split_state)) has_split=split_state
    has_shedding=.false.
    if (present(shedding_state)) has_shedding=shedding_state
    has_khrt=.false.
    if (present(khrt_state)) has_khrt=khrt_state
    expected_magic=spray_checkpoint_magic(present(boundary_totals),has_breakup,has_split,has_khrt)
    ok=.false.; message='Invalid, truncated or incompatible spray checkpoint'
    if (present(wall_films).and..not.has_films) then
      if (allocated(wall_films)) then
        if (size(wall_films)>0) return
      end if
    end if
    if (present(source_totals)) then
      expected_magic=trim(expected_magic)//'_FORCE_V1'
      if (size(source_totals)/=size(state,1)) return
    end if
    if (has_materials) expected_magic=trim(expected_magic)//'_MAT_V1'
    if (has_films) then
      if (.not.present(wall_films)) return
      expected_magic=trim(expected_magic)//'_FILM_V1'
    end if
    if (has_shedding) then
      if (.not.has_khrt.or..not.has_split) return
      expected_magic=trim(expected_magic)//'_SHED_V1'
    end if
    if (present(component_mass)) expected_magic=trim(expected_magic)//'_MIX_V1'
    if ((has_films.and.present(component_mass)).neqv.present(film_component_mass)) return
    if (present(film_component_mass)) expected_magic=trim(expected_magic)//'_FMIX_V1'
    if (present(last_id)) expected_magic=trim(expected_magic)//'_ID_V1'
    if ((has_split.or.has_khrt).and..not.has_breakup) return
    if (len_trim(context)==0 .or. len_trim(context)>len(line)) return
    if (size(state,1)/=reactive_nvar(size(species)) .or. size(species)<2) return
    if (size(reference_totals)/=size(state,1) .or. size(injected_totals)/=size(state,1)) return
    dims=shape(temperature)
    if (any(shape(state(1,:,:,:))/=dims)) return
    open(newunit=unit,file=trim(path),status='old',action='read',iostat=ios)
    if (ios/=0) return
    read(unit,'(a)',iostat=ios) line
    if (ios/=0) goto 900
    if (trim(line)/=trim(expected_magic)) goto 900
    if (present(boundary_totals)) then
      if (size(boundary_totals)/=size(state,1)) goto 900
    end if
    read(unit,'(a)',iostat=ios) line
    if (ios/=0) goto 900
    if (trim(line)/=trim(context)) goto 900
    call read_integer_row(unit,got,ios)
    if (ios/=0) goto 900
    if (any(got(1:4)/=[size(state,1),dims])) goto 900
    if (got(5)<0 .or. got(5)>1000000) goto 900
    if (present(component_mass)) then
      call read_integer_row(unit,component_count,ios)
      if (ios/=0) goto 900
      if (component_count(1)<1.or.component_count(1)>min(16,size(species))) goto 900
      allocate(mass_candidate(component_count(1),got(5)))
    end if
    if (present(last_id)) then
      call read_tokens(unit,tokens(1:1),ios)
      if (ios/=0) goto 900
      if (verify(trim(tokens(1)),'0123456789+-')/=0) goto 900
      read(tokens(1),*,iostat=ios) candidate_last_id
      if (ios/=0) goto 900
      if (candidate_last_id<0) goto 900
    end if
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
    if (present(boundary_totals)) then
      call read_real_row(unit,boundary,ios)
      if (ios/=0) goto 900
    end if
    if (present(source_totals)) then
      call read_real_row(unit,source,ios)
      if (ios/=0) goto 900
      if (.not.all(ieee_is_finite(source))) goto 900
    end if
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
      call read_tokens(unit,tokens(1:10+merge(2,0,has_breakup)),ios)
      if (ios/=0) goto 900
      if (verify(trim(tokens(1)), '0123456789+-')/=0) goto 900
      read(tokens(1),*,iostat=ios) drops(p)%id
      if (ios/=0) goto 900
      do i=1,9+merge(2,0,has_breakup)
        read(tokens(i+1),*,iostat=ios) row(i)
        if (ios/=0) goto 900
      end do
      drops(p)%position=row(1:3); drops(p)%velocity=row(4:6)
      if (has_breakup) then
        drops(p)%deformation=row(10); drops(p)%deformation_rate=row(11)
      end if
      drops(p)%mass=row(7); drops(p)%temperature=row(8); drops(p)%multiplicity=row(9)
      if (ios/=0) goto 900
      if (has_materials) then
        call read_tokens(unit,tokens(1:1),ios)
        if (ios/=0) goto 900
        if (verify(trim(tokens(1)),'0123456789+-')/=0) goto 900
        read(tokens(1),*,iostat=ios) drops(p)%material_id
        if (ios/=0) goto 900
      end if
      if (has_split.or.has_films) then
        call read_tokens(unit,tokens(1:1),ios)
        if (ios/=0) goto 900
        if (verify(trim(tokens(1)),'0123456789+-')/=0) goto 900
        read(tokens(1),*,iostat=ios) drops(p)%parent_id
        if (ios/=0) goto 900
        if (drops(p)%parent_id>0) then
          if (.not.parent_exists(drops,p)) goto 900
        end if
      end if
      if (has_khrt) then
        call read_real_row(unit,row(1:1),ios)
        if (ios/=0) goto 900
        drops(p)%rt_time=row(1)
      end if
      if (has_shedding) then
        call read_real_row(unit,row(1:2),ios)
        if (ios/=0) goto 900
        drops(p)%khrt_initial_count=row(1);drops(p)%khrt_shed_volume=row(2)
      end if
      if (present(component_mass)) then
        call read_real_row(unit,mass_candidate(:,p),ios)
        if (ios/=0) goto 900
        if (.not.valid_component_inventory(drops(p),mass_candidate(:,p))) goto 900
      end if
      if (.not. valid_spray_parcel(drops(p))) goto 900
      if (p>1) then
        if (drops(p)%id<=drops(p-1)%id) goto 900
      end if
    end do
    if (has_films) then
      call read_tokens(unit,tokens(1:1),ios)
      if (ios/=0) goto 900
      if (verify(trim(tokens(1)),'0123456789+-')/=0) goto 900
      read(tokens(1),*,iostat=ios) film_count
      if (ios/=0) goto 900
      if (film_count<0.or.film_count>1000000) goto 900
      allocate(patches(film_count))
      if (present(film_component_mass)) allocate(film_mass_candidate(component_count(1),film_count))
      do p=1,film_count
        call read_tokens(unit,tokens(1:1),ios)
        if (ios/=0) goto 900
        if (verify(trim(tokens(1)),'0123456789+-')/=0) goto 900
        read(tokens(1),*,iostat=ios) patches(p)%material_id
        if (ios/=0) goto 900
        call read_real_row(unit,row(1:10),ios)
        if (ios/=0) goto 900
        patches(p)%position=row(1:3);patches(p)%normal=row(4:6)
        patches(p)%mass=row(7);patches(p)%temperature=row(8);patches(p)%diameter=row(9);patches(p)%height=row(10)
        if (.not.valid_wall_film(patches(p))) goto 900
        if (present(film_component_mass)) then
          call read_real_row(unit,film_mass_candidate(:,p),ios)
          if (ios/=0) goto 900
          if (.not.valid_wall_component_inventory(patches(p),film_mass_candidate(:,p))) goto 900
        end if
      end do
    end if
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
    if (present(last_id)) then
      if (any(drops%id>candidate_last_id)) return
      last_id=candidate_last_id
    end if
    state=candidate; temperature=temp; parcels=drops; time=next_time; steps=next_steps; event=next_event
    if (present(boundary_totals)) boundary_totals=boundary
    if (present(source_totals)) source_totals=source
    if (has_films) call move_alloc(patches,wall_films)
    if (present(component_mass)) call move_alloc(mass_candidate,component_mass)
    if (present(film_component_mass)) call move_alloc(film_mass_candidate,film_component_mass)
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
