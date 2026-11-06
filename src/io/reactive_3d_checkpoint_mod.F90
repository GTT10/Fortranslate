module reactive_3d_checkpoint_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,nasa7_eos_policy_identity
  use reactive_3d_mod, only: recover_reactive_temperatures_3d,reactive_integrals_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_c_binding, only: c_char,c_int,c_null_char
  implicit none
  private
  public :: write_regular_gas_checkpoint,read_regular_gas_checkpoint
  abstract interface
    subroutine regular_checkpoint_identity(unit,status)
      integer, intent(in) :: unit
      integer, intent(out) :: status
    end subroutine
  end interface
  interface
    function rename_file(old,new) bind(C,name='rename') result(status)
      import c_char,c_int
      character(kind=c_char), intent(in) :: old(*),new(*)
      integer(c_int) :: status
    end function
  end interface
contains
  subroutine write_regular_gas_checkpoint(path,identity,species,state,temperature,time,steps,initial,boundary,source, &
      metrics,ok,message,reaction_source,heat_release,rotating_potential,previous_timestep, &
      internal_energy_density,energy_correction)
    character(len=*), intent(in) :: path
    procedure(regular_checkpoint_identity) :: identity
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),time,initial(:),boundary(:),source(:),metrics(2)
    integer, intent(in) :: steps
    real(dp), optional, intent(in) :: reaction_source(:,:,:,:),heat_release(:,:,:),rotating_potential(:,:,:)
    real(dp), optional, intent(in) :: internal_energy_density(:,:,:),energy_correction
    real(dp), optional, intent(in) :: previous_timestep
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    integer :: unit,ios,closed,i,j,k,n(3),nv
    real(dp), allocatable :: recovered(:,:,:)
    character(len=:), allocatable :: temporary
    ok=.false.;message='Regular gas checkpoint write rejected';n=shape(temperature);nv=size(state,1)
    if(len_trim(path)==0.or.any(n<1).or.size(species)+5/=nv) return
    if(any(shape(state)/=[nv,n]).or.size(initial)/=nv.or.size(boundary)/=nv.or.size(source)/=nv) return
    if(.not.all(ieee_is_finite([time,initial,boundary,source,metrics]))) return
    if(time<0.or.steps<0.or.metrics(1)<0.or.metrics(2)<0.or.metrics(2)>1) return
    if(present(internal_energy_density).neqv.present(energy_correction)) return
    if(present(internal_energy_density)) then
      if(any(shape(internal_energy_density)/=n)) return
      if(.not.all(ieee_is_finite(internal_energy_density)).or..not.ieee_is_finite(energy_correction)) return
      if(abs(energy_correction)>huge(1._dp)/8) return
    endif
    if(present(previous_timestep)) then
      if(.not.ieee_is_finite(previous_timestep)) return
      if(previous_timestep<0.or.previous_timestep>time) return
      if((steps==0).neqv.(previous_timestep==0)) return
    endif
    if(.not.all(ieee_is_finite(state)).or..not.all(ieee_is_finite(temperature))) return
    if(any(state(1,:,:,:)<=0).or.any(state(6:,:,:,:)<0).or.any(temperature<=0)) return
    if(maxval(abs(state))>1e30_dp) return
    if(present(reaction_source)) then
      if(any(shape(reaction_source)/=[nv-4,n])) return
      if(.not.all(ieee_is_finite(reaction_source))) return
    endif
    if(present(heat_release)) then
      if(any(shape(heat_release)/=n)) return
      if(.not.all(ieee_is_finite(heat_release))) return
    endif
    allocate(recovered,mold=temperature)
    call recover_reactive_temperatures_3d(species,state,temperature,n(1),n(2),n(3), &
      recovered,ok,rotating_potential,internal_energy_density)
    if(.not.ok) return
    ok=.false.
    if(any(abs(recovered-temperature)>1e-10_dp+1e-12_dp*abs(temperature))) return
    temporary=trim(path)//'.tmp'
    open(newunit=unit,file=temporary,status='new',action='write',iostat=ios)
    if(ios/=0) then
      message='Cannot create regular checkpoint temporary file';return
    endif
    if(present(internal_energy_density)) then
      write(unit,'(a)',iostat=ios) 'PELEF_REGULAR_GAS_V3'
    else if(present(previous_timestep)) then
      write(unit,'(a)',iostat=ios) 'PELEF_REGULAR_GAS_V2'
    else
      write(unit,'(a)',iostat=ios) 'PELEF_REGULAR_GAS_V1'
    endif
    if(ios==0) call identity(unit,ios)
    if(ios==0.and.all(species%pinned_eos_extrapolation)) &
      write(unit,'(a)',iostat=ios) trim(nasa7_eos_policy_identity(species))
    if(ios==0) write(unit,'(a)',iostat=ios) 'END_REGULAR_GAS_IDENTITY'
    if(ios==0) then
      if(present(internal_energy_density)) then
        write(unit,'(7(i0,1x))',iostat=ios) nv,n,merge(1,0,present(reaction_source)), &
          merge(1,0,present(heat_release)),merge(1,0,present(previous_timestep))
      else
        write(unit,'(6(i0,1x))',iostat=ios) nv,n,merge(1,0,present(reaction_source)),merge(1,0,present(heat_release))
      endif
    endif
    if(ios==0) write(unit,'(es27.18e3,1x,i0)',iostat=ios) time,steps
    if(ios==0) write(unit,'(*(es27.18e3,1x))',iostat=ios) initial
    if(ios==0) write(unit,'(*(es27.18e3,1x))',iostat=ios) boundary
    if(ios==0) write(unit,'(*(es27.18e3,1x))',iostat=ios) source
    if(ios==0) write(unit,'(*(es27.18e3,1x))',iostat=ios) metrics
    if(present(previous_timestep).and.ios==0) then
      write(unit,'(a)',iostat=ios) 'PELEF_REGULAR_GAS_TIMESTEP_V1'
      if(ios==0) write(unit,'(es27.18e3)',iostat=ios) previous_timestep
    endif
    if(present(energy_correction).and.ios==0) then
      write(unit,'(es27.18e3)',iostat=ios) energy_correction
    endif
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(ios/=0) exit
      write(unit,'(*(es27.18e3,1x))',iostat=ios) temperature(i,j,k),state(:,i,j,k)
      if(ios==0.and.present(internal_energy_density)) &
        write(unit,'(es27.18e3)',iostat=ios) internal_energy_density(i,j,k)
      if(ios==0.and.present(reaction_source)) write(unit,'(*(es27.18e3,1x))',iostat=ios) reaction_source(:,i,j,k)
      if(ios==0.and.present(heat_release)) write(unit,'(es27.18e3)',iostat=ios) heat_release(i,j,k)
    enddo
    enddo
    enddo
    if(ios==0) write(unit,'(a)',iostat=ios) 'END_PELEF_REGULAR_GAS'
    if(ios/=0) then
      close(unit,status='delete');return
    endif
    close(unit,iostat=closed)
    if(closed/=0) then
      message='Regular checkpoint close failed';return
    endif
    ios=rename_file(temporary//c_null_char,trim(path)//c_null_char)
    if(ios/=0) then
      message='Regular checkpoint atomic rename failed';return
    endif
    ok=.true.;message=''
  end subroutine

  subroutine read_regular_gas_checkpoint(path,identity,species,state,temperature,time,steps,initial,boundary,source, &
      metrics,ok,message,reaction_source,heat_release,spacing,final_time,rotating_potential,previous_timestep, &
      internal_energy_density,energy_correction)
    character(len=*), intent(in) :: path
    procedure(regular_checkpoint_identity) :: identity
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),time,initial(:),boundary(:),source(:),metrics(2)
    integer, intent(inout) :: steps
    real(dp), optional, intent(inout) :: reaction_source(:,:,:,:),heat_release(:,:,:)
    real(dp), optional, intent(in) :: spacing(3),final_time,rotating_potential(:,:,:)
    real(dp), optional, intent(inout) :: internal_energy_density(:,:,:),energy_correction
    real(dp), optional, intent(inout) :: previous_timestep
    real(dp) :: integrated(size(initial)),correction
    real(dp), allocatable :: auxiliary(:,:,:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), allocatable :: u(:,:,:,:),t(:,:,:),ir(:,:,:,:),heat(:,:,:),recovered(:,:,:),row(:)
    real(dp) :: stamp,totals(size(initial)),wall(size(boundary)),force(size(source)),limits(2),one(1),last_interval
    integer :: unit,expected,ios,code,n(3),nv,dimensions(7),count,i,j,k
    character(len=65536) :: line,reference
    ok=.false.;message='Regular gas checkpoint is invalid, truncated or incompatible';n=shape(temperature);nv=size(state,1)
    if(any(n<1).or.nv/=size(species)+5.or.any(shape(state)/=[nv,n])) return
    if(size(initial)/=nv.or.size(boundary)/=nv.or.size(source)/=nv) return
    if(present(internal_energy_density).neqv.present(energy_correction)) return
    if(present(internal_energy_density)) then
      if(any(shape(internal_energy_density)/=n)) return
      allocate(auxiliary,mold=internal_energy_density)
    endif
    correction=0
    if(present(reaction_source)) then
      if(any(shape(reaction_source)/=[nv-4,n])) return
    endif
    if(present(heat_release)) then
      if(any(shape(heat_release)/=n)) return
    endif
    open(newunit=expected,status='scratch',action='readwrite',iostat=ios)
    if(ios/=0) return
    call identity(expected,ios)
    if(ios==0.and.all(species%pinned_eos_extrapolation)) &
      write(expected,'(a)',iostat=ios) trim(nasa7_eos_policy_identity(species))
    if(ios/=0) then
      close(expected);return
    endif
    rewind(expected)
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    if(ios/=0) then
      close(expected);return
    endif
    read(unit,'(a)',iostat=ios) line
    if(ios/=0) goto 900
    if(present(internal_energy_density)) then
      if(trim(line)/='PELEF_REGULAR_GAS_V3') goto 900
    else if(present(previous_timestep)) then
      if(trim(line)/='PELEF_REGULAR_GAS_V2') goto 900
    else
      if(trim(line)/='PELEF_REGULAR_GAS_V1') goto 900
    endif
    do
      read(expected,'(a)',iostat=code) reference
      if(code<0) exit
      if(code/=0) goto 900
      read(unit,'(a)',iostat=ios) line
      if(ios/=0.or.line/=reference) goto 900
    enddo
    read(unit,'(a)',iostat=ios) line
    if(ios/=0.or.trim(line)/='END_REGULAR_GAS_IDENTITY') goto 900
    read(unit,'(a)',iostat=ios) line
    if(ios/=0) goto 900
    if(present(internal_energy_density)) then
      if(words(line)/=7) goto 900
      read(line,*,iostat=ios) dimensions
      if(ios/=0) goto 900
      if(dimensions(7)/=merge(1,0,present(previous_timestep))) goto 900
    else
      if(words(line)/=6) goto 900
      read(line,*,iostat=ios) dimensions(:6)
    endif
    if(ios/=0) goto 900
    if(any(dimensions(:6)/=[nv,n,merge(1,0,present(reaction_source)),merge(1,0,present(heat_release))])) goto 900
    read(unit,'(a)',iostat=ios) line
    if(ios/=0.or.words(line)/=2) goto 900
    read(line,*,iostat=ios) stamp,count
    if(ios/=0) goto 900
    if(.not.ieee_is_finite(stamp)) goto 900
    if(stamp<0.or.count<0) goto 900
    if(present(final_time)) then
      if(stamp>final_time) goto 900
    endif
    call real_row(unit,totals,ios);if(ios/=0) goto 900
    call real_row(unit,wall,ios);if(ios/=0) goto 900
    call real_row(unit,force,ios);if(ios/=0) goto 900
    call real_row(unit,limits,ios);if(ios/=0) goto 900
    if(present(previous_timestep)) then
      read(unit,'(a)',iostat=ios) line
      if(ios/=0.or.trim(line)/='PELEF_REGULAR_GAS_TIMESTEP_V1') goto 900
      call real_row(unit,one,ios);if(ios/=0) goto 900
      last_interval=one(1)
      if(last_interval<0.or.last_interval>stamp) goto 900
      if((count==0).neqv.(last_interval==0)) goto 900
    endif
    if(limits(1)<0.or.limits(2)<0.or.limits(2)>1) goto 900
    if(maxval(abs(totals))>huge(1._dp)/8.or.maxval(abs(wall))>huge(1._dp)/8.or. &
       maxval(abs(force))>huge(1._dp)/8) goto 900
    if(present(energy_correction)) then
      call real_row(unit,one,ios);if(ios/=0) goto 900
      correction=one(1)
      if(abs(correction)>huge(1._dp)/8) goto 900
    endif
    allocate(u,mold=state);allocate(t,mold=temperature);allocate(recovered,mold=t);allocate(row(nv+1))
    if(present(reaction_source)) allocate(ir,mold=reaction_source)
    if(present(heat_release)) allocate(heat,mold=heat_release)
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      call real_row(unit,row,ios);if(ios/=0) goto 900
      t(i,j,k)=row(1);u(:,i,j,k)=row(2:)
      if(present(internal_energy_density)) then
        call real_row(unit,one,ios);if(ios/=0) goto 900
        auxiliary(i,j,k)=one(1)
      endif
      if(present(reaction_source)) then
        call real_row(unit,ir(:,i,j,k),ios);if(ios/=0) goto 900
      endif
      if(present(heat_release)) then
        call real_row(unit,one,ios);if(ios/=0) goto 900
        heat(i,j,k)=one(1)
      endif
    enddo
    enddo
    enddo
    read(unit,'(a)',iostat=ios) line
    if(ios/=0.or.trim(line)/='END_PELEF_REGULAR_GAS') goto 900
    read(unit,'(a)',iostat=ios) line
    if(ios>=0) goto 900
    if(any(u(1,:,:,:)<=0).or.any(u(6:,:,:,:)<0).or.any(t<=0).or.maxval(abs(u))>1e30_dp) goto 900
    call recover_reactive_temperatures_3d(species,u,t,n(1),n(2),n(3),recovered,ok,rotating_potential,auxiliary)
    if(.not.ok) goto 900
    ok=.false.
    if(any(abs(recovered-t)>1e-10_dp+1e-12_dp*abs(t))) goto 900
    if(present(spacing)) then
      if(.not.all(ieee_is_finite(spacing))) goto 900
      if(any(spacing<=0)) goto 900
      call reactive_integrals_3d(u,n(1),n(2),n(3),spacing(1),spacing(2),spacing(3),integrated,ok)
      if(.not.ok) goto 900
      ok=.false.
      if(totals(1)<=0.or.any(totals(6:)<0)) goto 900
      integrated(5)=integrated(5)-correction
      if(maxval(abs(integrated(1:5)-totals(1:5)-wall(1:5)-force(1:5))/max(1._dp,abs(totals(1:5))))>5e-11_dp) goto 900
    endif
    state=u;temperature=t;time=stamp;steps=count;initial=totals;boundary=wall;source=force;metrics=limits
    if(present(internal_energy_density)) internal_energy_density=auxiliary
    if(present(energy_correction)) energy_correction=correction
    if(present(previous_timestep)) previous_timestep=last_interval
    if(present(reaction_source)) reaction_source=ir
    if(present(heat_release)) heat_release=heat
    close(unit);close(expected);ok=.true.;message='';return
900 continue
    close(unit);close(expected)
  end subroutine

  integer function words(line) result(count)
    character(len=*), intent(in) :: line
    integer :: i
    logical :: inside,space
    count=0;inside=.false.
    do i=1,len_trim(line)
      space=line(i:i)==' '.or.iachar(line(i:i))==9
      if(.not.space.and..not.inside) count=count+1
      inside=.not.space
    enddo
  end function
  subroutine real_row(unit,row,status)
    integer, intent(in) :: unit
    real(dp), intent(out) :: row(:)
    integer, intent(out) :: status
    character(len=65536) :: line
    read(unit,'(a)',iostat=status) line
    if(status/=0) return
    if(words(line)/=size(row)) then
      status=1;return
    endif
    read(line,*,iostat=status) row
    if(status==0) then
      if(.not.all(ieee_is_finite(row))) status=1
    endif
  end subroutine
end module reactive_3d_checkpoint_mod
