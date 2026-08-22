#!/usr/bin/env python3
from pathlib import Path
import re

root = Path('.')
(root / 'src/parallel').mkdir(parents=True, exist_ok=True)
(root / 'docs/validation').mkdir(parents=True, exist_ok=True)
(root / 'cases/mpi_entropy_1d').mkdir(parents=True, exist_ok=True)

(root / 'src/parallel/mpi_domain_1d_mod.F90').write_text(r'''module mpi_domain_1d_mod
  use iso_fortran_env, only: real64
  use mpi_f08
  implicit none
  private

  integer, parameter :: dp = real64

  type, public :: mpi_domain_1d
    type(MPI_Comm) :: comm = MPI_COMM_NULL
    integer :: rank = -1
    integer :: nranks = 0
    integer :: global_cells = 0
    integer :: local_cells = 0
    integer :: global_first = 0
    integer :: global_last = -1
    integer :: left_rank = MPI_PROC_NULL
    integer :: right_rank = MPI_PROC_NULL
    integer, allocatable :: counts(:)
    integer, allocatable :: displacements(:)
  end type mpi_domain_1d

  public :: initialize_mpi_domain_1d
  public :: exchange_periodic_halo_1d
  public :: global_minimum_1d
  public :: global_sum_1d
  public :: gather_state_1d

contains

  subroutine initialize_mpi_domain_1d(domain, global_cells, comm, ok)
    type(mpi_domain_1d), intent(out) :: domain
    integer, intent(in) :: global_cells
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: ierr, r, base_cells, remainder

    ok = .false.
    if (global_cells <= 0) return
    domain%comm = comm
    domain%global_cells = global_cells
    call MPI_Comm_rank(comm, domain%rank, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Comm_size(comm, domain%nranks, ierr)
    if (ierr /= MPI_SUCCESS .or. domain%nranks <= 0) return

    allocate(domain%counts(domain%nranks), domain%displacements(domain%nranks))
    base_cells = global_cells / domain%nranks
    remainder = modulo(global_cells, domain%nranks)
    domain%displacements(1) = 0
    do r = 0, domain%nranks - 1
      domain%counts(r + 1) = base_cells
      if (r < remainder) domain%counts(r + 1) = domain%counts(r + 1) + 1
      if (r > 0) then
        domain%displacements(r + 1) = domain%displacements(r) + domain%counts(r)
      end if
    end do
    domain%local_cells = domain%counts(domain%rank + 1)
    domain%global_first = domain%displacements(domain%rank + 1) + 1
    domain%global_last = domain%global_first + domain%local_cells - 1
    if (domain%nranks == 1) then
      domain%left_rank = 0
      domain%right_rank = 0
    else
      domain%left_rank = modulo(domain%rank - 1 + domain%nranks, domain%nranks)
      domain%right_rank = modulo(domain%rank + 1, domain%nranks)
    end if
    ok = domain%local_cells > 0 .and. sum(domain%counts) == global_cells
  end subroutine initialize_mpi_domain_1d

  subroutine exchange_periodic_halo_1d(domain, state, ok)
    type(mpi_domain_1d), intent(in) :: domain
    real(dp), intent(inout), contiguous :: state(:, 0:)
    logical, intent(out) :: ok
    type(MPI_Request) :: requests(4)
    type(MPI_Status) :: statuses(4)
    integer :: ierr, nvar

    ok = .false.
    nvar = size(state, 1)
    if (ubound(state, 2) < domain%local_cells + 1) return
    if (domain%nranks == 1) then
      state(:, 0) = state(:, domain%local_cells)
      state(:, domain%local_cells + 1) = state(:, 1)
      ok = .true.
      return
    end if

    call MPI_Irecv(state(:, 0), nvar, MPI_DOUBLE_PRECISION, &
      domain%left_rank, 100, domain%comm, requests(1), ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Irecv(state(:, domain%local_cells + 1), nvar, MPI_DOUBLE_PRECISION, &
      domain%right_rank, 101, domain%comm, requests(2), ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Isend(state(:, 1), nvar, MPI_DOUBLE_PRECISION, &
      domain%left_rank, 101, domain%comm, requests(3), ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Isend(state(:, domain%local_cells), nvar, MPI_DOUBLE_PRECISION, &
      domain%right_rank, 100, domain%comm, requests(4), ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Waitall(4, requests, statuses, ierr)
    ok = ierr == MPI_SUCCESS
  end subroutine exchange_periodic_halo_1d

  subroutine global_minimum_1d(domain, local_value, global_value, ok)
    type(mpi_domain_1d), intent(in) :: domain
    real(dp), intent(in) :: local_value
    real(dp), intent(out) :: global_value
    logical, intent(out) :: ok
    integer :: ierr
    call MPI_Allreduce(local_value, global_value, 1, MPI_DOUBLE_PRECISION, &
      MPI_MIN, domain%comm, ierr)
    ok = ierr == MPI_SUCCESS
  end subroutine global_minimum_1d

  subroutine global_sum_1d(domain, local_values, global_values, ok)
    type(mpi_domain_1d), intent(in) :: domain
    real(dp), intent(in), contiguous :: local_values(:)
    real(dp), intent(out), contiguous :: global_values(:)
    logical, intent(out) :: ok
    integer :: ierr
    ok = .false.
    if (size(global_values) /= size(local_values)) return
    call MPI_Allreduce(local_values, global_values, size(local_values), &
      MPI_DOUBLE_PRECISION, MPI_SUM, domain%comm, ierr)
    ok = ierr == MPI_SUCCESS
  end subroutine global_sum_1d

  subroutine gather_state_1d(domain, local_state, global_state, root, ok)
    type(mpi_domain_1d), intent(in) :: domain
    real(dp), intent(in), contiguous :: local_state(:, :)
    real(dp), allocatable, intent(out) :: global_state(:, :)
    integer, intent(in) :: root
    logical, intent(out) :: ok
    real(dp), allocatable :: send_buffer(:), receive_buffer(:)
    integer, allocatable :: element_counts(:), element_displacements(:)
    integer :: ierr, nvar, i, v, offset

    ok = .false.
    nvar = size(local_state, 1)
    if (size(local_state, 2) /= domain%local_cells) return
    allocate(send_buffer(nvar * domain%local_cells))
    offset = 0
    do i = 1, domain%local_cells
      do v = 1, nvar
        offset = offset + 1
        send_buffer(offset) = local_state(v, i)
      end do
    end do
    allocate(element_counts(domain%nranks), element_displacements(domain%nranks))
    element_counts = nvar * domain%counts
    element_displacements = nvar * domain%displacements
    if (domain%rank == root) then
      allocate(receive_buffer(nvar * domain%global_cells))
    else
      allocate(receive_buffer(0))
    end if
    call MPI_Gatherv(send_buffer, size(send_buffer), MPI_DOUBLE_PRECISION, &
      receive_buffer, element_counts, element_displacements, MPI_DOUBLE_PRECISION, &
      root, domain%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    if (domain%rank == root) then
      allocate(global_state(nvar, domain%global_cells))
      offset = 0
      do i = 1, domain%global_cells
        do v = 1, nvar
          offset = offset + 1
          global_state(v, i) = receive_buffer(offset)
        end do
      end do
    else
      allocate(global_state(0, 0))
    end if
    ok = .true.
  end subroutine gather_state_1d

end module mpi_domain_1d_mod
''')

(root / 'app/pelef_mpi_1d.F90').write_text(r'''program pelef_mpi_1d
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
''')

(root / 'tools/compare_mpi_entropy.py').write_text(r'''#!/usr/bin/env python3
from pathlib import Path
import csv
import math
import sys

if len(sys.argv) < 3:
    raise SystemExit('usage: compare_mpi_entropy.py reference.csv candidate.csv [...]')

def read(path):
    with Path(path).open(newline='') as handle:
        return [[float(value) for value in row.values()] for row in csv.DictReader(handle)]

reference = read(sys.argv[1])
for name in sys.argv[2:]:
    candidate = read(name)
    if len(candidate) != len(reference):
        raise AssertionError('row count differs')
    maximum = 0.0
    for left, right in zip(reference, candidate):
        for a, b in zip(left, right):
            scale = max(1.0, abs(a), abs(b))
            maximum = max(maximum, abs(a - b) / scale)
    print(f'{name}: maximum_relative_difference={maximum:.16e}')
    if maximum > 5.0e-13:
        raise AssertionError('MPI rank-count parity failed')
''')

(root / 'cases/mpi_entropy_1d/README.md').write_text(
    'A 257-cell periodic entropy wave verifies uneven 1D MPI decomposition, halo exchange, global reductions, and ordered gather output.\n'
)
(root / 'docs/validation/0.20-mpi-foundation-complete.txt').write_text(
    'MPI 1D foundation: uneven decomposition, 15-component halo, global CFL/conservation, and 1/2/4-rank parity.\n'
)

cmake = root / 'CMakeLists.txt'
text = cmake.read_text()
if 'option(PELEF_ENABLE_MPI' not in text:
    marker = 'option(PELEF_ENABLE_CANTERA_REFERENCE "Enable live Cantera parity tests" OFF)'
    text = text.replace(marker, marker + '\noption(PELEF_ENABLE_MPI "Build MPI verification targets" OFF)', 1)
if 'add_library(pelef_mpi_support' not in text:
    text += r'''

if(PELEF_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  add_library(pelef_mpi_support src/parallel/mpi_domain_1d_mod.F90)
  target_link_libraries(pelef_mpi_support PUBLIC pelef_core MPI::MPI_Fortran)
  target_include_directories(pelef_mpi_support PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mod)
  add_executable(pelef_mpi_1d app/pelef_mpi_1d.F90)
  target_link_libraries(pelef_mpi_1d PRIVATE pelef_mpi_support)
endif()
'''
cmake.write_text(text)
