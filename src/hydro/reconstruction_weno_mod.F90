module reconstruction_weno_mod
  use precision_mod, only: dp
  implicit none
  private

  public :: weno_reconstruct_5js
  public :: weno_reconstruct_5z

contains

  pure subroutine weno_reconstruct_5js(stencil, left_edge, right_edge)
    real(dp), intent(in) :: stencil(5)
    real(dp), intent(out) :: left_edge, right_edge

    real(dp), parameter :: epsilon_weno = 1.0e-6_dp
    real(dp) :: candidate(3), beta(3), alpha(3), normalization

    call weno5_smoothness(stencil, beta)
    beta = 1.0_dp / (epsilon_weno + beta)**2
    alpha = [3.0_dp * beta(1), 6.0_dp * beta(2), beta(3)]
    normalization = 1.0_dp / sum(alpha)
    candidate(1) = &
      (2.0_dp * stencil(3) + 5.0_dp * stencil(4) - stencil(5)) / 6.0_dp
    candidate(2) = &
      (-stencil(2) + 5.0_dp * stencil(3) + 2.0_dp * stencil(4)) / 6.0_dp
    candidate(3) = &
      (2.0_dp * stencil(1) - 7.0_dp * stencil(2) + &
        11.0_dp * stencil(3)) / 6.0_dp
    right_edge = normalization * sum(alpha * candidate)

    call weno5_smoothness(stencil(5:1:-1), beta)
    beta = 1.0_dp / (epsilon_weno + beta)**2
    alpha = [3.0_dp * beta(1), 6.0_dp * beta(2), beta(3)]
    normalization = 1.0_dp / sum(alpha)
    candidate(1) = &
      (2.0_dp * stencil(3) + 5.0_dp * stencil(2) - stencil(1)) / 6.0_dp
    candidate(2) = &
      (-stencil(4) + 5.0_dp * stencil(3) + 2.0_dp * stencil(2)) / 6.0_dp
    candidate(3) = &
      (2.0_dp * stencil(5) - 7.0_dp * stencil(4) + &
        11.0_dp * stencil(3)) / 6.0_dp
    left_edge = normalization * sum(alpha * candidate)
  end subroutine weno_reconstruct_5js

  pure subroutine weno_reconstruct_5z(stencil, left_edge, right_edge)
    real(dp), intent(in) :: stencil(5)
    real(dp), intent(out) :: left_edge, right_edge

    real(dp), parameter :: epsilon_weno = 1.0e-6_dp
    real(dp) :: candidate(3), beta(3), alpha(3), normalization, tau

    call weno5_smoothness(stencil, beta)
    tau = abs(beta(3) - beta(1))
    beta = 1.0_dp + (tau / (epsilon_weno + beta))**2
    alpha = [3.0_dp * beta(1), 6.0_dp * beta(2), beta(3)]
    normalization = 1.0_dp / sum(alpha)
    candidate(1) = &
      (2.0_dp * stencil(3) + 5.0_dp * stencil(4) - stencil(5)) / 6.0_dp
    candidate(2) = &
      (-stencil(2) + 5.0_dp * stencil(3) + 2.0_dp * stencil(4)) / 6.0_dp
    candidate(3) = &
      (2.0_dp * stencil(1) - 7.0_dp * stencil(2) + &
        11.0_dp * stencil(3)) / 6.0_dp
    right_edge = normalization * sum(alpha * candidate)

    call weno5_smoothness(stencil(5:1:-1), beta)
    tau = abs(beta(3) - beta(1))
    beta = 1.0_dp + (tau / (epsilon_weno + beta))**2
    alpha = [3.0_dp * beta(1), 6.0_dp * beta(2), beta(3)]
    normalization = 1.0_dp / sum(alpha)
    candidate(1) = &
      (2.0_dp * stencil(3) + 5.0_dp * stencil(2) - stencil(1)) / 6.0_dp
    candidate(2) = &
      (-stencil(4) + 5.0_dp * stencil(3) + 2.0_dp * stencil(2)) / 6.0_dp
    candidate(3) = &
      (2.0_dp * stencil(5) - 7.0_dp * stencil(4) + &
        11.0_dp * stencil(3)) / 6.0_dp
    left_edge = normalization * sum(alpha * candidate)
  end subroutine weno_reconstruct_5z

  pure subroutine weno5_smoothness(stencil, beta)
    real(dp), intent(in) :: stencil(5)
    real(dp), intent(out) :: beta(3)

    beta(3) = 13.0_dp / 12.0_dp * &
      (stencil(1) - 2.0_dp * stencil(2) + stencil(3))**2 + &
      0.25_dp * &
      (stencil(1) - 4.0_dp * stencil(2) + 3.0_dp * stencil(3))**2
    beta(2) = 13.0_dp / 12.0_dp * &
      (stencil(2) - 2.0_dp * stencil(3) + stencil(4))**2 + &
      0.25_dp * (stencil(2) - stencil(4))**2
    beta(1) = 13.0_dp / 12.0_dp * &
      (stencil(3) - 2.0_dp * stencil(4) + stencil(5))**2 + &
      0.25_dp * &
      (3.0_dp * stencil(3) - 4.0_dp * stencil(4) + stencil(5))**2
  end subroutine weno5_smoothness

end module reconstruction_weno_mod
