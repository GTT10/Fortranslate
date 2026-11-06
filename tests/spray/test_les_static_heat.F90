! Fixed PeleC LES.H evaluates cp at arithmetic face T,Y for all static models.
! A pure-N2 temperature jump isolates this from PeleF's species-diffusion policy.
program test_les_static_heat
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,universal_gas_constant
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use les_smagorinsky_3d_mod, only: compute_les_ghosted_rhs_3d,smagorinsky_stress,wale_stress,vreman_stress
  implicit none
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: u(15,6,6,6),t(6,6,6),rhs(15,2,2,2),q(15),sound
  real(dp) :: fx(15,0:2,2,2),fy(15,2,0:2,2),fz(15,2,2,0:2),flux(15)
  real(dp) :: gradient(3,3),stress(3,3),mu,diffusion,cp,gas_constant,temperature,heat,expected,velocity(3),error
  real(dp), parameter :: spacing(3)=[.001_dp,.001_dp,.001_dp],cs=.17_dp,prt=.7_dp
  real(dp), parameter :: diagonal(3)=[100._dp,200._dp,300._dp]
  character(len=16), parameter :: models(3)=[character(len=16)::'smagorinsky','wale','vreman']
  real(dp), parameter :: left_t(3)=[300._dp,900._dp,600._dp],right_t(3)=[900._dp,300._dp,600._dp]
  integer :: i,j,k,d,m,b,nitrogen,index(3)
  logical :: ok
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'load pinned thermo')
  nitrogen=0
  do i=1,size(species)
    if(trim(species(i)%name)=='N2') nitrogen=i
  enddo
  call require(nitrogen>0,'N2 species')
  gas_constant=universal_gas_constant/species(nitrogen)%molecular_weight
  ! Independent NASA7 polynomial at the face temperature, rather than an
  ! average of endpoint heat capacities. No rate/transport fits are changed.
  temperature=600
  cp=gas_constant*sum(species(nitrogen)%low_coefficients(1:5)* &
    [1._dp,temperature,temperature**2,temperature**3,temperature**4])
  gradient=0
  do d=1,3
    gradient(d,d)=diagonal(d)
  enddo
  do d=1,3
    do b=1,3
    do k=1,6
      do j=1,6
        do i=1,6
          index=[i,j,k];temperature=left_t(b)
          if(index(d)>3) temperature=right_t(b)
          q=0;q(1)=1;q(2:4)=diagonal*(index-3)*spacing
          q(5+nitrogen)=1
          q(5)=mixture_specific_gas_constant(species,q(6:),ok)*temperature
          call require(ok,'mixture gas constant')
          call reactive_primitive_to_conserved(species,q,u(:,i,j,k),t(i,j,k),sound,ok)
          call require(ok,'admissible supplied halo')
        enddo
      enddo
    enddo
    do m=1,3
      select case(m)
      case(1)
        call smagorinsky_stress(gradient,1._dp,spacing,cs,mu,stress,ok)
      case(2)
        call wale_stress(gradient,1._dp,spacing,cs,mu,stress,ok)
      case(3)
        call vreman_stress(gradient,1._dp,spacing,cs,mu,stress,ok)
      end select
      call require(ok.and.mu>0,'active static closure')
      call compute_les_ghosted_rhs_3d(species,u,t,[2,2,2],spacing,cs,prt,.7_dp,rhs,diffusion,ok, &
        model=trim(models(m)),fx=fx,fy=fy,fz=fz)
      call require(ok.and.diffusion>0,'static face evaluation')
      select case(d)
      case(1)
        flux=fx(:,1,1,1)
      case(2)
        flux=fy(:,1,1,1)
      case(3)
        flux=fz(:,1,1,1)
      end select
      velocity=0;velocity(d)=.5_dp*diagonal(d)*spacing(d)
      heat=flux(5)-dot_product(velocity,flux(2:4))
      expected=-mu*cp/prt*(right_t(b)-left_t(b))/spacing(d)
      error=abs(heat-expected)/max(1._dp,abs(expected))
      write(*,'(a,2(1x,i0),4(1x,es24.16))') trim(models(m)),d,b,cp,heat,expected,error
      call require(error<2e-12_dp,'heat uses cp at face temperature')
      call require(all(flux(6:)==0),'pure composition has no species flux')
      call require(all(abs(flux(2:4)+stress(:,d))<2e-12_dp*max(1._dp,abs(stress(:,d)))), &
        'static stress unchanged')
    enddo
    enddo
  enddo
  call compute_les_ghosted_rhs_3d(species,u,t,[2,2,2],spacing,0._dp,prt,.7_dp,rhs,diffusion,ok, &
    fx=fx,fy=fy,fz=fz)
  call require(ok.and.diffusion==0.and.all(rhs==0),'zero coefficient RHS no-op')
  call require(all(fx==0).and.all(fy==0).and.all(fz==0),'zero coefficient face no-op')
  call check_face_eos_rejection()
  print *, 'PASS: static nonlinear-cp face heat, three closures and three directions'
contains
  subroutine check_face_eos_rejection()
    type(nasa7_species) :: synthetic(size(species))
    ! A test-only caloric curve has valid endpoints but cv<0 at their mean.
    ! This exercises late face-EOS rejection, separately from physical parity.
    synthetic=species
    synthetic(nitrogen)%low_coefficients=[4.5_dp,-1._dp/75,1._dp/90000,0._dp,0._dp,0._dp,0._dp]
    synthetic(nitrogen)%high_coefficients=synthetic(nitrogen)%low_coefficients
    do k=1,6
      do j=1,6
        do i=1,6
          temperature=300;if(i>3) temperature=900
          q=0;q(1)=1;q(2:4)=diagonal*([i,j,k]-3)*spacing;q(5+nitrogen)=1
          q(5)=mixture_specific_gas_constant(synthetic,q(6:),ok)*temperature
          call require(ok,'synthetic gas constant')
          call reactive_primitive_to_conserved(synthetic,q,u(:,i,j,k),t(i,j,k),sound,ok)
          call require(ok,'synthetic endpoint admissible')
        enddo
      enddo
    enddo
    fx=91;fy=92;fz=93;rhs=94;diffusion=95
    call compute_les_ghosted_rhs_3d(synthetic,u,t,[2,2,2],spacing,cs,prt,.7_dp,rhs,diffusion,ok, &
      fx=fx,fy=fy,fz=fz)
    call require(.not.ok.and.all(rhs==0).and.diffusion==0,'invalid face EOS rejects raw rates')
    call require(all(fx==91).and.all(fy==92).and.all(fz==93),'late face EOS rejection preserves public faces')
  end subroutine
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition) then
      print *, 'FAIL: ',message
      error stop 1
    endif
  end subroutine
end program
