module time_integrator_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qu, qp
  use state_conversion_mod, only: conserved_to_primitive, state_is_physical
  use eos_ideal_mod, only: ideal_gas_sound_speed
  use boundary_conditions_mod, only: apply_boundary_conditions
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

  subroutine advance_ssprk2( &
      conserved, nx, dx, dt, gamma, ok, reconstruction, limiter, boundary_condition)
    integer, intent(in) :: nx
    real(dp), intent(inout) :: conserved(ncons, 0:nx + 1)
    real(dp), intent(in) :: dx, dt, gamma
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: reconstruction, limiter, boundary_condition

    real(dp), allocatable :: old_state(:, :), stage_state(:, :), rhs(:, :)
    character(len=32) :: reconstruction_name, limiter_name, boundary_name
    logical :: rhs_ok, boundary_ok
    integer :: i

    reconstruction_name = "pcm"
    limiter_name = "mc"
    boundary_name = "outflow"
    if (present(reconstruction)) reconstruction_name = trim(reconstruction)
    if (present(limiter)) limiter_name = trim(limiter)
    if (present(boundary_condition)) boundary_name = trim(boundary_condition)

    allocate(old_state(ncons, 0:nx + 1))
    allocate(stage_state(ncons, 0:nx + 1))
    allocate(rhs(ncons, nx))

    call apply_boundary_conditions(conserved, nx, boundary_name, boundary_ok)
    if (.not. boundary_ok) then
      ok = .false.
      return
    end if
    old_state = conserved

    call compute_euler_rhs( &
      old_state, nx, dx, gamma, rhs, rhs_ok, &
      reconstruction_name, limiter_name, boundary_name)
    if (.not. rhs_ok) then
      ok = .false.
      return
    end if

    stage_state = old_state
    do concurrent (i = 1:nx)
      stage_state(:, i) = old_state(:, i) + dt * rhs(:, i)
    end do
    call apply_boundary_conditions(stage_state, nx, boundary_name, boundary_ok)
    if (.not. boundary_ok) then
      ok = .false.
      return
    end if

    if (.not. all_cells_physical(stage_state, nx, gamma)) then
      ok = .false.
      return
    end if

    call compute_euler_rhs( &
      stage_state, nx, dx, gamma, rhs, rhs_ok, &
      reconstruction_name, limiter_name, boundary_name)
    if (.not. rhs_ok) then
      ok = .false.
      return
    end if

    do concurrent (i = 1:nx)
      conserved(:, i) = 0.5_dp * old_state(:, i) + &
        0.5_dp * (stage_state(:, i) + dt * rhs(:, i))
    end do
    call apply_boundary_conditions(conserved, nx, boundary_name, boundary_ok)
    if (.not. boundary_ok) then
      ok = .false.
      return
    end if
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
