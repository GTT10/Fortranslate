module les_dynamic_smagorinsky_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use les_smagorinsky_3d_mod, only: smagorinsky_stress
  implicit none
  private
  public :: dynamic_smagorinsky_quantities,dynamic_smagorinsky_coefficients
  ! Pinned PeleC Constants.H small_num(), including its coefficient floor.
  real(dp), parameter :: regularizer=1e-8_dp
contains
  ! Supplied cell-centered derivatives: gradient(velocity component, derivative).
  ! Filter Kij, RUT, alphaij, alpha and flux_T before coefficient evaluation.
  subroutine dynamic_smagorinsky_quantities(gradient,temperature_gradient,density,velocity,temperature, &
      spacing,filter_ratio,kij,rut,alphaij,alpha,flux_t,ok)
    real(dp), intent(in) :: gradient(3,3),temperature_gradient(3),density,velocity(3),temperature,spacing(3)
    integer, intent(in) :: filter_ratio
    real(dp), intent(out) :: kij(3,3),rut(3),alphaij(3,3),alpha,flux_t(3)
    logical, intent(out) :: ok
    real(dp) :: viscosity,strain(3,3),stress(3,3),magnitude
    integer :: i,j
    kij=0;rut=0;alphaij=0;alpha=0;flux_t=0;ok=.false.
    if (.not.all(ieee_is_finite([temperature_gradient,velocity,temperature]))) return
    if (maxval(abs(temperature_gradient))>1e30_dp.or.maxval(abs(velocity))>1e30_dp) return
    if (temperature<=0.or.temperature>1e8_dp.or.filter_ratio<1.or.filter_ratio>64) return
    call smagorinsky_stress(gradient,density,spacing,1._dp,viscosity,stress,ok)
    if (.not.ok) return
    viscosity=viscosity*real(filter_ratio,dp)**2
    alphaij=stress*real(filter_ratio,dp)**2
    strain=.5_dp*(gradient+transpose(gradient));magnitude=sqrt(2*sum(strain**2))
    alpha=2*viscosity*magnitude;flux_t=viscosity*temperature_gradient
    do j=1,3
      do i=1,3
        kij(i,j)=density*velocity(i)*velocity(j)
      end do
    end do
    rut=density*velocity*temperature
  end subroutine

  ! Inputs alpha* and flux_T are test-filtered grid quantities; beta* and
  ! test_flux_T come from the test-filtered primitive state at the test width.
  ! Output ordering: Cs^2, CI, Cs^2/PrT. This routine performs no spatial filter.
  subroutine dynamic_smagorinsky_coefficients(density,velocity,temperature,kij,rut,alphaij,alpha,flux_t, &
      betaij,beta,test_flux_t,coefficients,ok)
    real(dp), intent(in) :: density,velocity(3),temperature,kij(3,3),rut(3),alphaij(3,3),alpha,flux_t(3)
    real(dp), intent(in) :: betaij(3,3),beta,test_flux_t(3)
    real(dp), intent(out) :: coefficients(3)
    logical, intent(out) :: ok
    real(dp) :: l(3,3),m(3,3),t(3),ke(3),lkk,bma,candidate(3)
    integer :: i,j
    coefficients=0;ok=.false.
    if (.not.all(ieee_is_finite([density,velocity,temperature,rut,alpha,flux_t,beta,test_flux_t]))) return
    if (.not.all(ieee_is_finite(kij)).or..not.all(ieee_is_finite(alphaij))) return
    if (.not.all(ieee_is_finite(betaij))) return
    if (density<=0.or.density>1e8_dp.or.temperature<=0.or.temperature>1e8_dp) return
    if (maxval(abs(velocity))>1e30_dp) return
    ! Bound contractions before evaluating products under floating-point traps.
    if (maxval(abs([rut,alpha,flux_t,beta,test_flux_t]))>1e100_dp) return
    if (maxval(abs(kij))>1e100_dp.or.maxval(abs(alphaij))>1e100_dp.or.maxval(abs(betaij))>1e100_dp) return
    do j=1,3
      do i=1,3
        l(i,j)=kij(i,j)-density*velocity(i)*velocity(j)
      end do
    end do
    m=betaij-alphaij;t=flux_t-test_flux_t;ke=rut-density*velocity*temperature
    lkk=l(1,1)+l(2,2)+l(3,3);bma=beta-alpha+regularizer
    if (bma==0) return
    candidate(1)=max(sum(l*m)/(sum(m*m)+regularizer),regularizer)
    ! Apply the pinned [floor,1] clipping before division to avoid overflow
    ! when the regularized isotropic denominator is nonzero but very small.
    candidate(2)=regularizer
    if ((lkk>0.and.bma>0).or.(lkk<0.and.bma<0)) then
      candidate(2)=1
      if (abs(lkk)<abs(bma)) candidate(2)=max(lkk/bma,regularizer)
    end if
    candidate(3)=max(sum(ke*t)/(sum(t*t)+regularizer),regularizer)
    if (.not.all(ieee_is_finite(candidate))) return
    coefficients=candidate;ok=.true.
  end subroutine
end module
