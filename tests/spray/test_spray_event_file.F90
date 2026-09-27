program test_spray_event_file
  use precision_mod, only: dp
  use spray_injection_profile_mod, only: spray_event_schedule,read_spray_event_schedule
  implicit none
  type(spray_event_schedule) :: schedule
  character(len=1024) :: path
  integer :: unit,i,mode
  logical :: ok
  call get_command_argument(1,path)
  if(len_trim(path)==0) error stop 'missing path'
  open(newunit=unit,file=trim(path),status='replace')
  do i=1,2048
    write(unit,*) real(i,dp)*1e-6_dp,real(i+1,dp)*1e-6_dp,modulo(i-1,16)+1,1._dp
  end do
  close(unit)
  call read_spray_event_schedule(trim(path),schedule,ok)
  if(.not.ok) error stop 'long schedule rejected'
  if(size(schedule%times)/=2048.or.schedule%ports(2048)/=16) error stop 'schedule truncated'
  do mode=1,5
    open(newunit=unit,file=trim(path),status='replace')
    select case(mode)
    case(1);write(unit,'(a)') '1 0 1 1'
    case(2);write(unit,'(a)') '0 1 17 1'
    case(3);write(unit,'(a)') '0 1 1.5 1'
    case(4);write(unit,'(a)') '0 1 1 1 extra'
    case(5);write(unit,'(a)') '1 2 1 1','0 1 1 1'
    end select
    close(unit)
    call read_spray_event_schedule(trim(path),schedule,ok)
    if(ok) error stop 'invalid event accepted'
    if(size(schedule%times)/=2048.or.schedule%ports(2048)/=16) error stop 'failed read changed schedule'
  end do
  open(newunit=unit,file=trim(path),status='old')
  close(unit,status='delete')
end program
