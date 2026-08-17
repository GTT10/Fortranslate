module time_integrator_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qu, qp
  use state_conversion_mod, only: conserved_to_primitive, state_is_physical
  use eos_ideal_mod, only: ideal_gas_sound_speed
  use boundary_conditions_mod, only: apply_outflow_boundaries
  use finite_volume_mod, only: compute_euler_rhs
  implicit none
  private

  public :: compute_cfl_timestep
  public :: advance_ssprk2
  public :: all_cells_physical

contains

  subroutine compute_cfl_timestep(conserved, nx, dx, gamma, cfl, dt, ok)
    integer, intent(in) :: nx
    real(dp), intent(in) :: conserved(ncons, 0:nx + 1)
    real(dp), intent(in) :: dx, gamma, cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok
    real(dp) :: primitive(nprim), sound_speed, signal_speed, max_signal_speed
    logical :: cell_ok
    integer :: i

    max_signal_speed = 0.0_dp
    ok = cfl > 0.0_dp .and. cfl <= 1.0_dp .and. dx > 0.0_dp
    if (.not. ok) then
      dt = 0.0_dp
      return
    end if

    do i = 1, nx
      call conserved_to_primitive(conserved(:, i), gamma, primitive, cell_ok)
      if (.not. cell_ok) then
        ok = .false.
        dt = 0.0_dp
        return
      end if
      sound_speed = ideal_gas_sound_speed(primitive(qrho), primitive(qp), gamma)
      signal_speed = abs(primitive(qu)) + sound_speed
      max_signal_speed = max(max_signal_speed, signal_speed)
    end do

    if (max_signal_speed <= 0.0_dp) then
      ok = .false.
      dt = 0.0_dp
      return
    end if

    dt = cfl * dx / max_signal_speed
  end subroutine compute_cfl_timestep

  subroutine advance_ssprk2(conserved, nx, dx, dt, gamma, ok)
    integer, intent(in) :: nx
    real(dp), intent(inout) :: conserved(ncons, 0:nx + 1)
    real(dp), intent(in) :: dx, dt, gamma
    logical, intent(out) :: ok
    real(dp), allocatable :: old_state(:, :), stage_state(:, :), rhs(:, :)
    logical :: rhs_ok
    integer :: i

    allocate(old_state(ncons, 0:nx + 1))
    allocate(stage_state(ncons, 0:nx + 1))
    allocate(rhs(ncons, nx))

    call apply_outflow_boundaries(conserved, nx)
    old_state = conserved

    call compute_euler_rhs(old_state, nx, dx, gamma, rhs, rhs_ok)
    if (.not. rhs_ok) then
      ok = .false.
      return
    end if

    stage_state = old_state
    do concurrent (i = 1:nx)
      stage_state(:, i) = old_state(:, i) + dt * rhs(:, i)
    end do
    call apply_outflow_boundaries(stage_state, nx)

    if (.not. all_cells_physical(stage_state, nx, gamma)) then
      ok = .false.
      return
    end if

    call compute_euler_rhs(stage_state, nx, dx, gamma, rhs, rhs_ok)
    if (.not. rhs_ok) then
      ok = .false.
      return
    end if

    do concurrent (i = 1:nx)
      conserved(:, i) = 0.5_dp * old_state(:, i) + &
        0.5_dp * (stage_state(:, i) + dt * rhs(:, i))
    end do
    call apply_outflow_boundaries(conserved, nx)
    ok = all_cells_physical(conserved, nx, gamma)
  end subroutine advance_ssprk2

  pure logical function all_cells_physical(conserved, nx, gamma) result(all_physical)
    integer, intent(in) :: nx
    real(dp), intent(in) :: conserved(ncons, 0:nx + 1)
    real(dp), intent(in) :: gamma
    integer :: i

    all_physical = .true.
    do i = 1, nx
      if (.not. state_is_physical(conserved(:, i), gamma)) then
        all_physical = .false.
        return
      end if
    end do
  end function all_cells_physical

end module time_integrator_mod
