program pelef_mpi_1d
  use iso_fortran_env, only: real64
  use mpi_f08
  use mpi_domain_1d_mod, only: mpi_domain_1d, initialize_mpi_domain_1d, &
    exchange_periodic_halo_1d, global_minimum_1d, global_sum_1d, gather_state_1d
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: nvar = 5
  real(dp), parameter :: gamma_gas = 1.4_dp
  real(dp), parameter :: pi = acos(-1.0_dp)
  type(mpi_domain_1d) :: domain
  real(dp), allocatable :: state(:, :), next_state(:, :), flux(:, :), global_state(:, :)
  real(dp) :: dx, time, final_time, cfl, local_dt, dt
  real(dp) :: local_integrals(nvar), initial_integrals(nvar), final_integrals(nvar)
  integer :: ierr, steps, i, g, unit, arg_count
  logical :: ok
  character(len=256) :: argument, output_file

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop 'MPI_Init failed'
  call initialize_mpi_domain_1d(domain, 257, MPI_COMM_WORLD, ok)
  if (.not. ok) error stop 'domain initialization failed'

  arg_count = command_argument_count()
  if (arg_count >= 1) then
    call get_command_argument(1, argument)
  else
    argument = 'mpi_entropy.csv'
  end if
  if (trim(argument) == '--halo-test') then
    call run_halo_test(domain)
    call MPI_Finalize(ierr)
    stop
  end if
  output_file = trim(argument)

  dx = 1.0_dp / real(domain%global_cells, dp)
  final_time = 0.20_dp
  cfl = 0.40_dp
  allocate(state(nvar, 0:domain%local_cells + 1))
  allocate(next_state(nvar, 0:domain%local_cells + 1))
  allocate(flux(nvar, 0:domain%local_cells))
  state = 0.0_dp
  do i = 1, domain%local_cells
    g = domain%global_first + i - 1
    call entropy_state((real(g, dp) - 0.5_dp) * dx, state(:, i))
  end do
  call exchange_periodic_halo_1d(domain, state, ok)
  if (.not. ok) error stop 'initial halo exchange failed'
  call compute_integrals(domain, state(:, 1:domain%local_cells), dx, initial_integrals)

  time = 0.0_dp
  steps = 0
  do while (time < final_time .and. steps < 100000)
    call local_stable_timestep(state(:, 1:domain%local_cells), dx, cfl, local_dt)
    call global_minimum_1d(domain, local_dt, dt, ok)
    if (.not. ok .or. dt <= 0.0_dp) error stop 'global timestep failed'
    if (time + dt > final_time) dt = final_time - time
    call exchange_periodic_halo_1d(domain, state, ok)
    if (.not. ok) error stop 'halo exchange failed'
    do i = 0, domain%local_cells
      call rusanov_flux(state(:, i), state(:, i + 1), flux(:, i))
    end do
    next_state = state
    do i = 1, domain%local_cells
      next_state(:, i) = state(:, i) - dt / dx * (flux(:, i) - flux(:, i - 1))
      if (next_state(1, i) <= 0.0_dp) error stop 'negative density'
    end do
    state(:, 1:domain%local_cells) = next_state(:, 1:domain%local_cells)
    time = time + dt
    steps = steps + 1
  end do
  if (steps >= 100000) error stop 'maximum steps reached'

  call compute_integrals(domain, state(:, 1:domain%local_cells), dx, final_integrals)
  if (maxval(abs(final_integrals - initial_integrals)) > 5.0e-11_dp) then
    error stop 'global conservation regression'
  end if
  call gather_state_1d(domain, state(:, 1:domain%local_cells), global_state, 0, ok)
  if (.not. ok) error stop 'state gather failed'
  if (domain%rank == 0) then
    open(newunit=unit, file=trim(output_file), status='replace', action='write')
    write(unit, '(a)') 'x,rho,rhou,rhov,rhow,rhoE'
    do i = 1, domain%global_cells
      write(unit, '(6(es24.16,:,","))') (real(i, dp) - 0.5_dp) * dx, global_state(:, i)
    end do
    close(unit)
    write(*, '(a,i0)') 'MPI ranks: ', domain%nranks
    write(*, '(a,i0)') 'Completed steps: ', steps
    write(*, '(a,es24.16)') 'Final time: ', time
  end if
  call MPI_Finalize(ierr)

contains

  subroutine entropy_state(x, conserved)
    real(dp), intent(in) :: x
    real(dp), intent(out) :: conserved(nvar)
    real(dp) :: rho, velocity, pressure
    rho = 1.0_dp + 0.2_dp * sin(2.0_dp * pi * x)
    velocity = 1.0_dp
    pressure = 1.0_dp
    conserved = 0.0_dp
    conserved(1) = rho
    conserved(2) = rho * velocity
    conserved(5) = pressure / (gamma_gas - 1.0_dp) + 0.5_dp * rho * velocity**2
  end subroutine entropy_state

  subroutine primitive(conserved, rho, velocity, pressure, sound_speed)
    real(dp), intent(in) :: conserved(nvar)
    real(dp), intent(out) :: rho, velocity, pressure, sound_speed
    real(dp) :: kinetic
    rho = conserved(1)
    velocity = conserved(2) / rho
    kinetic = 0.5_dp * (conserved(2)**2 + conserved(3)**2 + conserved(4)**2) / rho
    pressure = (gamma_gas - 1.0_dp) * (conserved(5) - kinetic)
    if (rho <= 0.0_dp .or. pressure <= 0.0_dp) error stop 'invalid Euler state'
    sound_speed = sqrt(gamma_gas * pressure / rho)
  end subroutine primitive

  subroutine physical_flux(conserved, result)
    real(dp), intent(in) :: conserved(nvar)
    real(dp), intent(out) :: result(nvar)
    real(dp) :: rho, u, p, c
    call primitive(conserved, rho, u, p, c)
    result(1) = conserved(2)
    result(2) = conserved(2) * u + p
    result(3) = conserved(3) * u
    result(4) = conserved(4) * u
    result(5) = (conserved(5) + p) * u
  end subroutine physical_flux

  subroutine rusanov_flux(left_state, right_state, result)
    real(dp), intent(in) :: left_state(nvar), right_state(nvar)
    real(dp), intent(out) :: result(nvar)
    real(dp) :: left_flux(nvar), right_flux(nvar)
    real(dp) :: rho_l, u_l, p_l, c_l, rho_r, u_r, p_r, c_r, speed
    call primitive(left_state, rho_l, u_l, p_l, c_l)
    call primitive(right_state, rho_r, u_r, p_r, c_r)
    call physical_flux(left_state, left_flux)
    call physical_flux(right_state, right_flux)
    speed = max(abs(u_l) + c_l, abs(u_r) + c_r)
    result = 0.5_dp * (left_flux + right_flux) - 0.5_dp * speed * (right_state - left_state)
  end subroutine rusanov_flux

  subroutine local_stable_timestep(local_state, spacing, courant, result)
    real(dp), intent(in) :: local_state(:, :), spacing, courant
    real(dp), intent(out) :: result
    real(dp) :: rho, u, p, c, speed
    integer :: cell
    result = huge(1.0_dp)
    do cell = 1, size(local_state, 2)
      call primitive(local_state(:, cell), rho, u, p, c)
      speed = abs(u) + c
      result = min(result, courant * spacing / speed)
    end do
  end subroutine local_stable_timestep

  subroutine compute_integrals(dom, local_state, spacing, global_integrals)
    type(mpi_domain_1d), intent(in) :: dom
    real(dp), intent(in) :: local_state(:, :), spacing
    real(dp), intent(out) :: global_integrals(nvar)
    logical :: reduction_ok
    local_integrals = sum(local_state, dim=2) * spacing
    call global_sum_1d(dom, local_integrals, global_integrals, reduction_ok)
    if (.not. reduction_ok) error stop 'global reduction failed'
  end subroutine compute_integrals

  subroutine run_halo_test(dom)
    type(mpi_domain_1d), intent(in) :: dom
    real(dp), allocatable :: values(:, :)
    integer :: cell, component, global_index, expected_left, expected_right
    logical :: exchange_ok
    allocate(values(15, 0:dom%local_cells + 1))
    values = -1.0_dp
    do cell = 1, dom%local_cells
      global_index = dom%global_first + cell - 1
      do component = 1, 15
        values(component, cell) = 1000.0_dp * real(component, dp) + real(global_index, dp)
      end do
    end do
    call exchange_periodic_halo_1d(dom, values, exchange_ok)
    if (.not. exchange_ok) error stop 'halo test exchange failed'
    expected_left = dom%global_first - 1
    if (expected_left < 1) expected_left = dom%global_cells
    expected_right = dom%global_last + 1
    if (expected_right > dom%global_cells) expected_right = 1
    do component = 1, 15
      if (abs(values(component, 0) - (1000.0_dp * component + expected_left)) > 1.0e-12_dp) then
        error stop 'left halo mismatch'
      end if
      if (abs(values(component, dom%local_cells + 1) - &
          (1000.0_dp * component + expected_right)) > 1.0e-12_dp) then
        error stop 'right halo mismatch'
      end if
    end do
    if (dom%rank == 0) write(*, '(a,i0)') '15-component halo test passed on ranks: ', dom%nranks
  end subroutine run_halo_test

end program pelef_mpi_1d
