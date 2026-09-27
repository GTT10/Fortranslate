program test_spray_long_rate_file
  use precision_mod, only: dp
  use spray_injection_profile_mod, only: spray_rate_profile,read_spray_rate_profile,integrate_spray_injection_profile
  implicit none
  type(spray_rate_profile) :: profile
  character(len=1024) :: path
  integer :: unit,i
  real(dp) :: mass
  logical :: ok
  call get_command_argument(1,path)
  if(len_trim(path)==0) error stop 'missing path'
  open(newunit=unit,file=trim(path),status='replace')
  do i=0,10000
    write(unit,'(2es27.18e3)') real(i,dp)/10000,real(i,dp)/20000
  end do
  close(unit)
  call read_spray_rate_profile(trim(path),profile,ok)
  if(.not.ok) error stop 'long profile rejected'
  if(size(profile%times)/=10001) error stop 'long profile truncated'
  call integrate_spray_injection_profile(profile%times,profile%rates,.125_dp,.875_dp,mass,ok)
  if(.not.ok.or.abs(mass-.1875_dp)>1e-13_dp) error stop 'linear rate integral'
  open(newunit=unit,file=trim(path),status='replace')
  write(unit,'(a)') '0 0','1 0.1','0.5 0.2'
  close(unit)
  call read_spray_rate_profile(trim(path),profile,ok)
  if(ok) error stop 'nonmonotone profile accepted'
  if(size(profile%times)/=10001.or.profile%rates(10001)/=.5_dp) error stop 'failure changed profile'
  open(newunit=unit,file=trim(path),status='old')
  close(unit,status='delete')
end program
