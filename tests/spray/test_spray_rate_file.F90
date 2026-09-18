program test_spray_rate_file
  use precision_mod, only: dp
  use spray_injection_profile_mod, only: read_spray_injection_profile
  implicit none
  character(len=1024) :: path
  real(dp) :: times(64),rates(64)
  integer :: unit,mode,i,count
  logical :: ok
  call get_command_argument(1,path)
  if (len_trim(path)==0) error stop 'missing fixture path'
  do mode=1,6
    open(newunit=unit,file=trim(path),status='replace',action='write')
    select case(mode)
    case(1)
      write(unit,'(a)') '# SI columns','','0 0.1','1 0.2 ! endpoint'
    case(2)
      write(unit,'(a)') '0 0.1 extra','1 0.2'
    case(3)
      write(unit,'(a)') '0 0.1','0 0.2'
    case(4)
      write(unit,'(a)') '0 0.1','1 NaN'
    case(5)
      write(unit,'(a)') '0 0.1'//repeat(' ',4100)//'extra','1 0.2'
    case(6)
      do i=0,64
        write(unit,*) i,0.1_dp
      end do
    end select
    close(unit)
    times=42;rates=43;count=44
    call read_spray_injection_profile(trim(path),times,rates,count,ok)
    if (mode==1) then
      if (.not.ok.or.count/=2) error stop 'valid profile rejected'
      if (any(times(:2)/=[0._dp,1._dp]).or.any(rates(:2)/=[.1_dp,.2_dp])) error stop 'wrong values'
      if (any(times(3:)/=0).or.any(rates(3:)/=0)) error stop 'padding'
    else
      if (ok.or.count/=44.or.any(times/=42).or.any(rates/=43)) error stop 'failure changed output'
    end if
  end do
  open(newunit=unit,file=trim(path),status='old')
  close(unit,status='delete')
  call read_spray_injection_profile(trim(path),times,rates,count,ok)
  if (ok.or.count/=44.or.any(times/=42).or.any(rates/=43)) error stop 'missing file changed output'
end program
