program dual_energy_probe
  use precision_mod, only: dp
  use dual_energy_mod, only: dual_energy_policy,reset_dual_energy_cell
  implicit none
  type(dual_energy_policy) :: policy
  real(dp) :: state(6),aux,temp,eta,phi,formation,delta,cv
  integer :: flags(3),ios
  logical :: ok
  do
    read(*,*,iostat=ios) state(:5),aux,temp,flags,eta,phi,formation,cv
    if(ios<0) exit
    if(ios/=0) error stop 'Invalid probe input'
    state(6)=state(1);policy=dual_energy_policy(flags(1)==1,flags(2)==1,flags(3)==1,eta)
    delta=0
    call reset_dual_energy_cell(state,aux,temp,policy,delta,ok,eos,phi)
    if(.not.ok) error stop 'Reset failed'
    write(*,'(3(es25.17e3,1x))') state(5),aux,delta
  enddo
contains
  subroutine eos(t,y,e,success)
    real(dp), intent(in) :: t,y(:)
    real(dp), intent(out) :: e
    logical, intent(out) :: success
    e=cv*t+formation;success=size(y)==1.and.all(y==1)
  end subroutine
end program dual_energy_probe
