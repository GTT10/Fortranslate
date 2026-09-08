module finite_volume_3d_mod
  use precision_mod, only: dp
  use state_indices_mod, only: &
    ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: conserved_to_primitive, state_is_physical
  use eos_ideal_mod, only: ideal_gas_sound_speed
  use riemann_flux_mod, only: compute_riemann_flux_x
  use directional_flux_mod, only: &
    compute_riemann_flux_y, compute_riemann_flux_z
  implicit none
  private

  public :: compute_euler_cfl_timestep_3d
  public :: compute_euler_rhs_3d
  public :: advance_euler_ssprk2_3d
  public :: all_cells_physical_3d

contains

  subroutine compute_euler_cfl_timestep_3d( &
      conserved, nx, ny, nz, dx, dy, dz, gamma, cfl, dt, ok)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: conserved(ncons, nx, ny, nz)
    real(dp), intent(in) :: dx, dy, dz, gamma, cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp) :: primitive(nprim), sound_speed, signal_rate, maximum_rate
    logical :: cell_ok
    integer :: i, j, k

    dt = 0.0_dp
    ok = nx >= 2 .and. ny >= 2 .and. nz >= 2 .and. &
      dx > 0.0_dp .and. dy > 0.0_dp .and. dz > 0.0_dp .and. &
      gamma > 1.0_dp .and. cfl > 0.0_dp .and. cfl <= 1.0_dp
    if (.not. ok) return

    maximum_rate = 0.0_dp
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call conserved_to_primitive( &
            conserved(:, i, j, k), gamma, primitive, cell_ok)
          if (.not. cell_ok) then
            ok = .false.
            return
          end if
          sound_speed = ideal_gas_sound_speed( &
            primitive(qrho), primitive(qp), gamma)
          signal_rate = (abs(primitive(qu)) + sound_speed) / dx + &
            (abs(primitive(qv)) + sound_speed) / dy + &
            (abs(primitive(qw)) + sound_speed) / dz
          maximum_rate = max(maximum_rate, signal_rate)
        end do
      end do
    end do

    if (maximum_rate <= 0.0_dp) then
      ok = .false.
      return
    end if
    dt = cfl / maximum_rate
  end subroutine compute_euler_cfl_timestep_3d

  subroutine compute_euler_rhs_3d( &
      conserved, nx, ny, nz, dx, dy, dz, gamma, riemann_solver, rhs, ok)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: conserved(ncons, nx, ny, nz)
    real(dp), intent(in) :: dx, dy, dz, gamma
    character(len=*), intent(in) :: riemann_solver
    real(dp), intent(out) :: rhs(ncons, nx, ny, nz)
    logical, intent(out) :: ok

    real(dp), allocatable :: flux_x(:, :, :, :)
    real(dp), allocatable :: flux_y(:, :, :, :)
    real(dp), allocatable :: flux_z(:, :, :, :)
    logical :: face_ok
    integer :: i, j, k, next_i, next_j, next_k
    integer :: previous_i, previous_j, previous_k

    rhs = 0.0_dp
    ok = nx >= 2 .and. ny >= 2 .and. nz >= 2 .and. &
      dx > 0.0_dp .and. dy > 0.0_dp .and. dz > 0.0_dp .and. &
      gamma > 1.0_dp
    if (.not. ok) return

    allocate(flux_x(ncons, nx, ny, nz))
    allocate(flux_y(ncons, nx, ny, nz))
    allocate(flux_z(ncons, nx, ny, nz))

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          next_i = modulo(i, nx) + 1
          call compute_riemann_flux_x( &
            conserved(:, i, j, k), conserved(:, next_i, j, k), gamma, &
            riemann_solver, flux_x(:, i, j, k), face_ok)
          if (.not. face_ok) then
            ok = .false.
            return
          end if
        end do
      end do
    end do

    do k = 1, nz
      do j = 1, ny
        next_j = modulo(j, ny) + 1
        do i = 1, nx
          call compute_riemann_flux_y( &
            conserved(:, i, j, k), conserved(:, i, next_j, k), gamma, &
            riemann_solver, flux_y(:, i, j, k), face_ok)
          if (.not. face_ok) then
            ok = .false.
            return
          end if
        end do
      end do
    end do

    do k = 1, nz
      next_k = modulo(k, nz) + 1
      do j = 1, ny
        do i = 1, nx
          call compute_riemann_flux_z( &
            conserved(:, i, j, k), conserved(:, i, j, next_k), gamma, &
            riemann_solver, flux_z(:, i, j, k), face_ok)
          if (.not. face_ok) then
            ok = .false.
            return
          end if
        end do
      end do
    end do

    do k = 1, nz
      previous_k = modulo(k - 2, nz) + 1
      do j = 1, ny
        previous_j = modulo(j - 2, ny) + 1
        do i = 1, nx
          previous_i = modulo(i - 2, nx) + 1
          rhs(:, i, j, k) = &
            -(flux_x(:, i, j, k) - &
              flux_x(:, previous_i, j, k)) / dx &
            -(flux_y(:, i, j, k) - &
              flux_y(:, i, previous_j, k)) / dy &
            -(flux_z(:, i, j, k) - &
              flux_z(:, i, j, previous_k)) / dz
        end do
      end do
    end do
  end subroutine compute_euler_rhs_3d

  subroutine advance_euler_ssprk2_3d( &
      conserved, nx, ny, nz, dx, dy, dz, dt, gamma, riemann_solver, ok)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(inout) :: conserved(ncons, nx, ny, nz)
    real(dp), intent(in) :: dx, dy, dz, dt, gamma
    character(len=*), intent(in) :: riemann_solver
    logical, intent(out) :: ok

    real(dp), allocatable :: old_state(:, :, :, :)
    real(dp), allocatable :: stage_state(:, :, :, :)
    real(dp), allocatable :: updated_state(:, :, :, :)
    real(dp), allocatable :: rhs(:, :, :, :)
    logical :: rhs_ok

    ok = dt > 0.0_dp .and. nx >= 2 .and. ny >= 2 .and. nz >= 2
    if (.not. ok) return

    allocate(old_state(ncons, nx, ny, nz))
    allocate(stage_state(ncons, nx, ny, nz))
    allocate(updated_state(ncons, nx, ny, nz))
    allocate(rhs(ncons, nx, ny, nz))
    old_state = conserved

    call compute_euler_rhs_3d( &
      old_state, nx, ny, nz, dx, dy, dz, gamma, riemann_solver, rhs, rhs_ok)
    if (.not. rhs_ok) then
      ok = .false.
      return
    end if
    stage_state = old_state + dt * rhs
    if (.not. all_cells_physical_3d(stage_state, nx, ny, nz, gamma)) then
      ok = .false.
      return
    end if

    call compute_euler_rhs_3d( &
      stage_state, nx, ny, nz, dx, dy, dz, gamma, riemann_solver, rhs, &
      rhs_ok)
    if (.not. rhs_ok) then
      ok = .false.
      return
    end if
    updated_state = 0.5_dp * old_state + 0.5_dp * (stage_state + dt * rhs)
    if (.not. all_cells_physical_3d(updated_state, nx, ny, nz, gamma)) then
      ok = .false.
      return
    end if

    conserved = updated_state
    ok = .true.
  end subroutine advance_euler_ssprk2_3d

  pure logical function all_cells_physical_3d( &
      conserved, nx, ny, nz, gamma) result(all_physical)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: conserved(ncons, nx, ny, nz)
    real(dp), intent(in) :: gamma
    integer :: i, j, k

    all_physical = .true.
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (.not. state_is_physical(conserved(:, i, j, k), gamma)) then
            all_physical = .false.
            return
          end if
        end do
      end do
    end do
  end function all_cells_physical_3d

end module finite_volume_3d_mod
