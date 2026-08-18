module time_integrator_multispecies_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qu, qp
  use multispecies_state_mod, only: &
    multispecies_nvar, multispecies_state_is_physical, species_closure_error, &
    synchronize_multispecies_thermodynamics
  use state_conversion_mod, only: conserved_to_primitive
  use eos_ideal_mod, only: ideal_gas_sound_speed
  use boundary_conditions_multispecies_mod, only: &
    apply_multispecies_boundary_conditions
  use finite_volume_multispecies_mod, only: compute_multispecies_rhs
  implicit none
  private

  public :: compute_multispecies_cfl_timestep
  public :: advance_multispecies_hydro_step
  public :: all_multispecies_cells_physical
  public :: maximum_species_closure_error

contains

  subroutine compute_multispecies_cfl_timestep( &
      conserved, nx, nspecies, dx, gamma, cfl, dt, ok)
    integer, intent(in) :: nx, nspecies
    real(dp), intent(in) :: conserved(:, 0:), dx, gamma, cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp) :: primitive(nprim), sound_speed, signal_speed, maximum_speed
    logical :: cell_ok
    integer :: i, nvar

    dt = 0.0_dp
    ok = .false.
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0 .or. size(conserved, 1) /= nvar) return
    if (ubound(conserved, 2) < nx + 1 .or. nx < 4 .or. dx <= 0.0_dp) return
    if (cfl <= 0.0_dp .or. cfl > 1.0_dp) return

    maximum_speed = 0.0_dp
    do i = 1, nx
      call conserved_to_primitive( &
        conserved(1:ncons, i), gamma, primitive, cell_ok)
      if (.not. cell_ok) return
      sound_speed = ideal_gas_sound_speed( &
        primitive(qrho), primitive(qp), gamma)
      if (sound_speed <= 0.0_dp) return
      signal_speed = abs(primitive(qu)) + sound_speed
      maximum_speed = max(maximum_speed, signal_speed)
    end do
    if (maximum_speed <= 0.0_dp) return

    dt = cfl * dx / maximum_speed
    ok = .true.
  end subroutine compute_multispecies_cfl_timestep

  subroutine advance_multispecies_hydro_step( &
      conserved, nx, nspecies, dx, dt, gamma, ok, reconstruction, limiter, &
      boundary_condition, riemann_solver, plm_order, use_flattening)
    integer, intent(in) :: nx, nspecies
    real(dp), intent(inout) :: conserved(:, 0:)
    real(dp), intent(in) :: dx, dt, gamma
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: reconstruction, limiter
    character(len=*), intent(in), optional :: boundary_condition
    character(len=*), intent(in), optional :: riemann_solver
    integer, intent(in), optional :: plm_order
    logical, intent(in), optional :: use_flattening

    character(len=32) :: reconstruction_name, limiter_name, boundary_name
    character(len=32) :: riemann_name
    integer :: slope_order
    logical :: flattening_enabled

    reconstruction_name = "pcm"
    limiter_name = "mc"
    boundary_name = "outflow"
    riemann_name = "rusanov"
    slope_order = 2
    flattening_enabled = .false.
    if (present(reconstruction)) reconstruction_name = trim(reconstruction)
    if (present(limiter)) limiter_name = trim(limiter)
    if (present(boundary_condition)) boundary_name = trim(boundary_condition)
    if (present(riemann_solver)) riemann_name = trim(riemann_solver)
    if (present(plm_order)) slope_order = plm_order
    if (present(use_flattening)) flattening_enabled = use_flattening

    if (trim(reconstruction_name) == "pelec_plm") then
      call advance_multispecies_godunov( &
        conserved, nx, nspecies, dx, dt, gamma, ok, limiter_name, &
        boundary_name, riemann_name, slope_order, flattening_enabled)
    else
      call advance_multispecies_ssprk2( &
        conserved, nx, nspecies, dx, dt, gamma, ok, reconstruction_name, &
        limiter_name, boundary_name, riemann_name)
    end if
  end subroutine advance_multispecies_hydro_step

  subroutine advance_multispecies_ssprk2( &
      conserved, nx, nspecies, dx, dt, gamma, ok, reconstruction, limiter, &
      boundary_condition, riemann_solver)
    integer, intent(in) :: nx, nspecies
    real(dp), intent(inout) :: conserved(:, 0:)
    real(dp), intent(in) :: dx, dt, gamma
    logical, intent(out) :: ok
    character(len=*), intent(in) :: reconstruction, limiter
    character(len=*), intent(in) :: boundary_condition, riemann_solver

    real(dp), allocatable :: old_state(:, :), stage_state(:, :), rhs(:, :)
    logical :: rhs_ok, boundary_ok
    integer :: nvar, i

    ok = .false.
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0 .or. size(conserved, 1) /= nvar) return
    if (trim(reconstruction) == "pelec_plm") return

    allocate(old_state(nvar, 0:nx + 1))
    allocate(stage_state(nvar, 0:nx + 1))
    allocate(rhs(nvar, 1:nx))

    call apply_multispecies_boundary_conditions( &
      conserved, nx, boundary_condition, boundary_ok)
    if (.not. boundary_ok) return
    old_state = conserved

    call compute_multispecies_rhs( &
      old_state, nx, nspecies, dx, gamma, rhs, rhs_ok, &
      reconstruction=reconstruction, limiter=limiter, &
      boundary_condition=boundary_condition, riemann_solver=riemann_solver)
    if (.not. rhs_ok) return

    stage_state = old_state
    do i = 1, nx
      stage_state(:, i) = old_state(:, i) + dt * rhs(:, i)
    end do
    call synchronize_multispecies_field(stage_state, nx, gamma, boundary_ok)
    if (.not. boundary_ok) return
    call apply_multispecies_boundary_conditions( &
      stage_state, nx, boundary_condition, boundary_ok)
    if (.not. boundary_ok) return
    if (.not. all_multispecies_cells_physical( &
        stage_state, nx, nspecies, gamma)) return

    call compute_multispecies_rhs( &
      stage_state, nx, nspecies, dx, gamma, rhs, rhs_ok, &
      reconstruction=reconstruction, limiter=limiter, &
      boundary_condition=boundary_condition, riemann_solver=riemann_solver)
    if (.not. rhs_ok) return

    do i = 1, nx
      conserved(:, i) = 0.5_dp * old_state(:, i) + &
        0.5_dp * (stage_state(:, i) + dt * rhs(:, i))
    end do
    call synchronize_multispecies_field(conserved, nx, gamma, boundary_ok)
    if (.not. boundary_ok) return
    call apply_multispecies_boundary_conditions( &
      conserved, nx, boundary_condition, boundary_ok)
    if (.not. boundary_ok) return
    ok = all_multispecies_cells_physical(conserved, nx, nspecies, gamma)
  end subroutine advance_multispecies_ssprk2

  subroutine advance_multispecies_godunov( &
      conserved, nx, nspecies, dx, dt, gamma, ok, limiter, &
      boundary_condition, riemann_solver, plm_order, use_flattening)
    integer, intent(in) :: nx, nspecies, plm_order
    real(dp), intent(inout) :: conserved(:, 0:)
    real(dp), intent(in) :: dx, dt, gamma
    logical, intent(out) :: ok
    character(len=*), intent(in) :: limiter, boundary_condition
    character(len=*), intent(in) :: riemann_solver
    logical, intent(in) :: use_flattening

    real(dp), allocatable :: old_state(:, :), rhs(:, :)
    logical :: rhs_ok, boundary_ok
    integer :: nvar, i

    ok = .false.
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0 .or. size(conserved, 1) /= nvar) return

    allocate(old_state(nvar, 0:nx + 1))
    allocate(rhs(nvar, 1:nx))

    call apply_multispecies_boundary_conditions( &
      conserved, nx, boundary_condition, boundary_ok)
    if (.not. boundary_ok) return
    old_state = conserved

    call compute_multispecies_rhs( &
      old_state, nx, nspecies, dx, gamma, rhs, rhs_ok, &
      reconstruction="pelec_plm", limiter=limiter, &
      boundary_condition=boundary_condition, riemann_solver=riemann_solver, &
      dt=dt, plm_order=plm_order, use_flattening=use_flattening)
    if (.not. rhs_ok) return

    do i = 1, nx
      conserved(:, i) = old_state(:, i) + dt * rhs(:, i)
    end do
    call synchronize_multispecies_field(conserved, nx, gamma, boundary_ok)
    if (.not. boundary_ok) return
    call apply_multispecies_boundary_conditions( &
      conserved, nx, boundary_condition, boundary_ok)
    if (.not. boundary_ok) return
    ok = all_multispecies_cells_physical(conserved, nx, nspecies, gamma)
  end subroutine advance_multispecies_godunov

  pure subroutine synchronize_multispecies_field( &
      conserved, nx, gamma, ok)
    integer, intent(in) :: nx
    real(dp), intent(inout) :: conserved(:, 0:)
    real(dp), intent(in) :: gamma
    logical, intent(out) :: ok
    logical :: cell_ok
    integer :: i

    ok = .false.
    if (ubound(conserved, 2) < nx + 1) return
    do i = 1, nx
      call synchronize_multispecies_thermodynamics( &
        conserved(:, i), gamma, cell_ok)
      if (.not. cell_ok) return
    end do
    ok = .true.
  end subroutine synchronize_multispecies_field

  pure logical function all_multispecies_cells_physical( &
      conserved, nx, nspecies, gamma) result(all_physical)
    integer, intent(in) :: nx, nspecies
    real(dp), intent(in) :: conserved(:, 0:), gamma
    integer :: i

    all_physical = .false.
    if (size(conserved, 1) /= multispecies_nvar(nspecies)) return
    if (ubound(conserved, 2) < nx + 1) return

    all_physical = .true.
    do i = 1, nx
      if (.not. multispecies_state_is_physical( &
          conserved(:, i), gamma, nspecies)) then
        all_physical = .false.
        return
      end if
    end do
  end function all_multispecies_cells_physical

  pure real(dp) function maximum_species_closure_error( &
      conserved, nx, nspecies) result(maximum_error)
    integer, intent(in) :: nx, nspecies
    real(dp), intent(in) :: conserved(:, 0:)
    integer :: i

    maximum_error = huge(1.0_dp)
    if (size(conserved, 1) /= multispecies_nvar(nspecies)) return
    if (ubound(conserved, 2) < nx + 1) return

    maximum_error = 0.0_dp
    do i = 1, nx
      maximum_error = max(maximum_error, &
        species_closure_error(conserved(:, i), nspecies))
    end do
  end function maximum_species_closure_error

end module time_integrator_multispecies_mod
