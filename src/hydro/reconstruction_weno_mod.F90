module reconstruction_weno_mod
  use precision_mod, only: dp
  implicit none
  private

  public :: weno_reconstruct_5js
  public :: weno_reconstruct_5z
  public :: weno_reconstruct_7z
  public :: weno_reconstruct_3z

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

  pure subroutine weno_reconstruct_7z(stencil, left_edge, right_edge)
    real(dp), intent(in) :: stencil(7)
    real(dp), intent(out) :: left_edge, right_edge

    call weno7z_right_edge(stencil, right_edge)
    call weno7z_right_edge(stencil(7:1:-1), left_edge)
  end subroutine weno_reconstruct_7z

  pure subroutine weno_reconstruct_3z(stencil, left_edge, right_edge)
    real(dp), intent(in) :: stencil(3)
    real(dp), intent(out) :: left_edge, right_edge

    call weno3z_right_edge(stencil, right_edge)
    call weno3z_right_edge(stencil(3:1:-1), left_edge)
  end subroutine weno_reconstruct_3z

  pure subroutine weno7z_right_edge(stencil, edge)
    real(dp), intent(in) :: stencil(7)
    real(dp), intent(out) :: edge

    real(dp), parameter :: epsilon_weno = 1.0e-6_dp
    real(dp), parameter :: linear_weights(4) = &
      [4.0_dp / 35.0_dp, 18.0_dp / 35.0_dp, &
        12.0_dp / 35.0_dp, 1.0_dp / 35.0_dp]
    real(dp) :: candidate(4), beta(4), alpha(4), normalization, tau

    beta(4) = stencil(1) * &
      (547.0_dp * stencil(1) - 3882.0_dp * stencil(2) + &
        4642.0_dp * stencil(3) - 1854.0_dp * stencil(4)) + &
      stencil(2) * &
      (7043.0_dp * stencil(2) - 17246.0_dp * stencil(3) + &
        7042.0_dp * stencil(4)) + &
      stencil(3) * &
      (11003.0_dp * stencil(3) - 9402.0_dp * stencil(4)) + &
      2107.0_dp * stencil(4)**2
    beta(3) = stencil(2) * &
      (267.0_dp * stencil(2) - 1642.0_dp * stencil(3) + &
        1602.0_dp * stencil(4) - 494.0_dp * stencil(5)) + &
      stencil(3) * &
      (2843.0_dp * stencil(3) - 5966.0_dp * stencil(4) + &
        1922.0_dp * stencil(5)) + &
      stencil(4) * &
      (3443.0_dp * stencil(4) - 2522.0_dp * stencil(5)) + &
      547.0_dp * stencil(5)**2
    beta(2) = stencil(3) * &
      (547.0_dp * stencil(3) - 2522.0_dp * stencil(4) + &
        1922.0_dp * stencil(5) - 494.0_dp * stencil(6)) + &
      stencil(4) * &
      (3443.0_dp * stencil(4) - 5966.0_dp * stencil(5) + &
        1602.0_dp * stencil(6)) + &
      stencil(5) * &
      (2843.0_dp * stencil(5) - 1642.0_dp * stencil(6)) + &
      267.0_dp * stencil(6)**2
    beta(1) = stencil(4) * &
      (2107.0_dp * stencil(4) - 9402.0_dp * stencil(5) + &
        7042.0_dp * stencil(6) - 1854.0_dp * stencil(7)) + &
      stencil(5) * &
      (11003.0_dp * stencil(5) - 17246.0_dp * stencil(6) + &
        4642.0_dp * stencil(7)) + &
      stencil(6) * &
      (7043.0_dp * stencil(6) - 3882.0_dp * stencil(7)) + &
      547.0_dp * stencil(7)**2

    tau = abs(beta(4) - beta(1))
    beta = 1.0_dp + (tau / (epsilon_weno + beta))**2
    alpha = linear_weights * beta
    normalization = 1.0_dp / sum(alpha)

    candidate(4) = &
      (-3.0_dp * stencil(1) + 13.0_dp * stencil(2) - &
        23.0_dp * stencil(3) + 25.0_dp * stencil(4)) / 12.0_dp
    candidate(3) = &
      (stencil(2) - 5.0_dp * stencil(3) + &
        13.0_dp * stencil(4) + 3.0_dp * stencil(5)) / 12.0_dp
    candidate(2) = &
      (-stencil(3) + 7.0_dp * stencil(4) + &
        7.0_dp * stencil(5) - stencil(6)) / 12.0_dp
    candidate(1) = &
      (3.0_dp * stencil(4) + 13.0_dp * stencil(5) - &
        5.0_dp * stencil(6) + stencil(7)) / 12.0_dp
    edge = normalization * sum(alpha * candidate)
  end subroutine weno7z_right_edge

  pure subroutine weno3z_right_edge(stencil, edge)
    real(dp), intent(in) :: stencil(3)
    real(dp), intent(out) :: edge

    real(dp), parameter :: epsilon_weno = 1.0e-6_dp
    real(dp) :: candidate(2), beta(2), alpha(2), normalization, tau

    beta(2) = (stencil(1) - stencil(2))**2
    beta(1) = (stencil(2) - stencil(3))**2
    tau = abs(beta(2) - beta(1))
    beta = 1.0_dp + (tau / (epsilon_weno + beta))**2
    alpha = [beta(1), 2.0_dp * beta(2)]
    normalization = 1.0_dp / sum(alpha)
    candidate(2) = (-stencil(1) + 3.0_dp * stencil(2)) / 2.0_dp
    candidate(1) = (stencil(2) + stencil(3)) / 2.0_dp
    edge = normalization * sum(alpha * candidate)
  end subroutine weno3z_right_edge

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
