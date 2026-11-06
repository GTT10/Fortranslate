program test_constant_particle_transport
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,nasa7_specific_gas_constant
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use transport_database_mod, only: load_h2o2_full_transport
  use gas_transport_mod, only: gas_transport_species,set_constant_transport_coefficients,constant_transport_values
  use mixture_transport_mod, only: mixture_transport_coefficients
  use spray_mixture_mod, only: film_fuel_diffusion_coefficients
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: source(:),runtime(:),saved(:)
  real(dp), allocatable :: y(:),diffusion(:),fuel_diffusion(:)
  real(dp) :: values(4),controls(4),mu,lambda,temperature,pressure,gas_constant,rho,mw
  integer :: i,stage
  logical :: ok,enabled
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_full_transport(source,ok)
  call require(ok,'transport')
  runtime=source;controls=[1.e-5_dp,.02_dp,3.e-5_dp,2.e-5_dp]
  call set_constant_transport_coefficients(runtime,controls,ok)
  call require(ok,'set constant model')
  allocate(y(size(species)),diffusion(size(species)),fuel_diffusion(size(species)))
  y=0
  do i=1,size(species)
    if(species(i)%name=='N2') y(i)=.8_dp
    if(species(i)%name=='H2O') y(i)=.2_dp
    call require(.not.allocated(source(i)%constant_coefficients),'source records preserved')
  enddo
  call require(sum(y)==1,'composition')
  gas_constant=sum(y*[(nasa7_specific_gas_constant(species(i)),i=1,size(species))])
  mw=1/sum(y/species%molecular_weight)
  ! Full constant transport must not depend on molecular collision data.
  runtime%diameter=0
  do stage=1,2
    temperature=300._dp*stage;pressure=101325._dp*stage**2
    rho=pressure/(gas_constant*temperature)
    call mixture_transport_coefficients(species,runtime,y,temperature,pressure,mu,lambda,diffusion,ok)
    call require(ok.and.mu==controls(1).and.lambda==controls(2),'constant skin mu/lambda')
    call require(maxval(abs(diffusion-controls(4)/rho))<1.e-14_dp*controls(4)/rho,'skin rhoD conversion')
    call film_fuel_diffusion_coefficients(species,runtime,y,temperature,pressure,diffusion,fuel_diffusion,ok)
    call require(ok,'fuel diffusion')
    call require(maxval(abs(fuel_diffusion-controls(4)/rho*mw/species%molecular_weight))< &
      1.e-14_dp*maxval(fuel_diffusion),'Drag H molecular-weight conversion')
  enddo
  controls(4)=0
  call set_constant_transport_coefficients(runtime,controls,ok)
  call require(ok,'set zero rhoD')
  call mixture_transport_coefficients(species,runtime,y,600._dp,101325._dp,mu,lambda,diffusion,ok)
  call require(ok.and.all(diffusion==0),'zero rhoD skin evaluation')
  call film_fuel_diffusion_coefficients(species,runtime,y,600._dp,101325._dp,diffusion,fuel_diffusion,ok)
  call require(ok.and.all(fuel_diffusion==0),'zero rhoD fuel conversion')
  controls(4)=2.e-5_dp
  call set_constant_transport_coefficients(runtime,controls,ok)
  call require(ok,'restore positive rhoD')
  saved=runtime
  call set_constant_transport_coefficients(runtime,[1._dp,2._dp,3._dp],ok)
  call require(.not.ok,'invalid control size rejected')
  call constant_transport_values(runtime,values,enabled,ok)
  call require(ok.and.enabled.and.all(values==controls),'invalid setter rollback')
  runtime(2)%constant_coefficients(4)=3.e-5_dp
  call mixture_transport_coefficients(species,runtime,y,600._dp,101325._dp,mu,lambda,diffusion,ok)
  call require(.not.ok.and.mu==0.and.lambda==0.and.all(diffusion==0),'mixed model controls rejected')
  runtime=source
  call set_constant_transport_coefficients(runtime,controls(:2),ok)
  call require(ok,'partial constants')
  call mixture_transport_coefficients(species,source,y,600._dp,101325._dp,mu,lambda,fuel_diffusion,ok)
  call require(ok,'molecular diffusion reference')
  call mixture_transport_coefficients(species,runtime,y,600._dp,101325._dp,mu,lambda,diffusion,ok)
  call require(ok.and.mu==controls(1).and.lambda==controls(2),'partial skin overrides')
  call require(all(diffusion==fuel_diffusion),'partial constants retain molecular diffusion')
  print *, 'PASS: constant skin transport, rhoD, fuel molecular-weight conversion, source preservation and rollback'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    endif
  end subroutine require
end program
