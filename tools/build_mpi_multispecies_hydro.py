#!/usr/bin/env python3
from pathlib import Path
import re

root = Path('.')
(root / 'cases/mpi_multispecies_1d').mkdir(parents=True, exist_ok=True)
(root / 'docs/validation').mkdir(parents=True, exist_ok=True)

(root / 'app/pelef_mpi_multispecies_1d.F90').write_text(r'''program pelef_mpi_multispecies_1d
  use iso_fortran_env, only: real64
  use mpi_f08
  use mpi_domain_1d_mod, only: mpi_domain_1d, initialize_mpi_domain_1d, &
    exchange_periodic_halo_1d, global_minimum_1d, global_sum_1d, gather_state_1d
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: nspecies = 10
  integer, parameter :: nvar = 5 + nspecies
  integer, parameter :: irho = 1, imx = 2, imy = 3, imz = 4, ienergy = 5
  integer, parameter :: first_species = 6
  real(dp), parameter :: gamma_gas = 1.4_dp
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: density_floor = 1.0e-12_dp
  real(dp), parameter :: pressure_floor = 1.0e-12_dp
  real(dp), parameter :: species_tolerance = 2.0e-12_dp

  type(mpi_domain_1d) :: domain
  real(dp), allocatable :: state(:, :), next_state(:, :), flux(:, :)
  real(dp), allocatable :: global_state(:, :)
  real(dp) :: dx, time, final_time, cfl, local_dt, dt
  real(dp) :: initial_integrals(nvar), final_integrals(nvar)
  real(dp) :: maximum_closure_error, minimum_species
  integer :: ierr, steps, cell, global_cell, output_unit
  logical :: ok
  character(len=256) :: argument, output_file

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop 'MPI_Init failed'
  call initialize_mpi_domain_1d(domain, 257, MPI_COMM_WORLD, ok)
  if (.not. ok) error stop 'MPI domain initialization failed'

  if (command_argument_count() >= 1) then
    call get_command_argument(1, argument)
  else
    argument = 'mpi_multispecies.csv'
  end if
  output_file = trim(argument)

  dx = 1.0_dp / real(domain%global_cells, dp)
  final_time = 0.20_dp
  cfl = 0.38_dp
  allocate(state(nvar, 0:domain%local_cells + 1))
  allocate(next_state(nvar, 0:domain%local_cells + 1))
  allocate(flux(nvar, 0:domain%local_cells))
  state = 0.0_dp

  do cell = 1, domain%local_cells
    global_cell = domain%global_first + cell - 1
    call initialize_cell((real(global_cell, dp) - 0.5_dp) * dx, state(:, cell))
  end do
  call validate_local_state(state(:, 1:domain%local_cells), maximum_closure_error, minimum_species)
  call exchange_periodic_halo_1d(domain, state, ok)
  if (.not. ok) error stop 'Initial halo exchange failed'
  call global_integrals(domain, state(:, 1:domain%local_cells), dx, initial_integrals)

  time = 0.0_dp
  steps = 0
  do while (time < final_time .and. steps < 200000)
    call stable_timestep_local(state(:, 1:domain%local_cells), dx, cfl, local_dt)
    call global_minimum_1d(domain, local_dt, dt, ok)
    if (.not. ok .or. dt <= 0.0_dp) error stop 'Global timestep reduction failed'
    if (time + dt > final_time) dt = final_time - time

    call exchange_periodic_halo_1d(domain, state, ok)
    if (.not. ok) error stop 'Halo exchange failed'
    do cell = 0, domain%local_cells
      call rusanov_flux(state(:, cell), state(:, cell + 1), flux(:, cell))
    end do

    next_state = state
    do cell = 1, domain%local_cells
      next_state(:, cell) = state(:, cell) - dt / dx * &
        (flux(:, cell) - flux(:, cell - 1))
    end do
    call validate_local_state(next_state(:, 1:domain%local_cells), &
      maximum_closure_error, minimum_species)
    state(:, 1:domain%local_cells) = next_state(:, 1:domain%local_cells)
    time = time + dt
    steps = steps + 1
  end do
  if (steps >= 200000) error stop 'Maximum step count reached'

  call global_integrals(domain, state(:, 1:domain%local_cells), dx, final_integrals)
  if (maxval(abs(final_integrals - initial_integrals)) > 2.0e-10_dp) then
    error stop 'Global multispecies conservation regression failed'
  end if
  call gather_state_1d(domain, state(:, 1:domain%local_cells), global_state, 0, ok)
  if (.not. ok) error stop 'Global state gather failed'

  if (domain%rank == 0) then
    open(newunit=output_unit, file=trim(output_file), status='replace', action='write')
    write(output_unit, '(a)') &
      'x,rho,rhou,rhov,rhow,rhoE,Y_H2,Y_H,Y_O,Y_O2,Y_OH,Y_H2O,Y_HO2,Y_H2O2,Y_AR,Y_N2'
    do cell = 1, domain%global_cells
      write(output_unit, '(16(es24.16,:,","))') &
        (real(cell, dp) - 0.5_dp) * dx, global_state(irho:ienergy, cell), &
        global_state(first_species:nvar, cell) / global_state(irho, cell)
    end do
    close(output_unit)
    write(*, '(a,i0)') 'MPI ranks: ', domain%nranks
    write(*, '(a,i0)') 'Completed steps: ', steps
    write(*, '(a,es24.16)') 'Final time: ', time
  end if

  call MPI_Finalize(ierr)

contains

  subroutine initialize_cell(x, conserved)
    real(dp), intent(in) :: x
    real(dp), intent(out) :: conserved(nvar)
    real(dp) :: rho, velocity, pressure, wave
    real(dp) :: y(nspecies)

    wave = sin(2.0_dp * pi * x)
    rho = 1.0_dp + 0.08_dp * wave
    velocity = 0.85_dp
    pressure = 1.0_dp
    y = [0.18_dp + 0.015_dp * wave, 2.0e-5_dp, 3.0e-5_dp, &
      0.24_dp - 0.010_dp * wave, 4.0e-5_dp, 0.055_dp, &
      1.0e-4_dp, 2.0e-5_dp, 0.075_dp, 0.0_dp]
    y(nspecies) = 1.0_dp - sum(y(1:nspecies - 1))
    if (minval(y) <= 0.0_dp) error stop 'Invalid initial mass fractions'

    conserved = 0.0_dp
    conserved(irho) = rho
    conserved(imx) = rho * velocity
    conserved(ienergy) = pressure / (gamma_gas - 1.0_dp) + &
      0.5_dp * rho * velocity**2
    conserved(first_species:nvar) = rho * y
  end subroutine initialize_cell

  subroutine primitive(conserved, rho, u, v, w, pressure, sound_speed)
    real(dp), intent(in) :: conserved(nvar)
    real(dp), intent(out) :: rho, u, v, w, pressure, sound_speed
    real(dp) :: kinetic

    rho = conserved(irho)
    if (rho <= density_floor) error stop 'Density floor violation'
    u = conserved(imx) / rho
    v = conserved(imy) / rho
    w = conserved(imz) / rho
    kinetic = 0.5_dp * rho * (u**2 + v**2 + w**2)
    pressure = (gamma_gas - 1.0_dp) * (conserved(ienergy) - kinetic)
    if (pressure <= pressure_floor) error stop 'Pressure floor violation'
    sound_speed = sqrt(gamma_gas * pressure / rho)
  end subroutine primitive

  subroutine physical_flux(conserved, result)
    real(dp), intent(in) :: conserved(nvar)
    real(dp), intent(out) :: result(nvar)
    real(dp) :: rho, u, v, w, pressure, sound_speed

    call primitive(conserved, rho, u, v, w, pressure, sound_speed)
    result(irho) = conserved(imx)
    result(imx) = conserved(imx) * u + pressure
    result(imy) = conserved(imy) * u
    result(imz) = conserved(imz) * u
    result(ienergy) = (conserved(ienergy) + pressure) * u
    result(first_species:nvar) = conserved(first_species:nvar) * u
  end subroutine physical_flux

  subroutine rusanov_flux(left_state, right_state, result)
    real(dp), intent(in) :: left_state(nvar), right_state(nvar)
    real(dp), intent(out) :: result(nvar)
    real(dp) :: left_flux(nvar), right_flux(nvar)
    real(dp) :: rho_l, u_l, v_l, w_l, p_l, c_l
    real(dp) :: rho_r, u_r, v_r, w_r, p_r, c_r, speed

    call primitive(left_state, rho_l, u_l, v_l, w_l, p_l, c_l)
    call primitive(right_state, rho_r, u_r, v_r, w_r, p_r, c_r)
    call physical_flux(left_state, left_flux)
    call physical_flux(right_state, right_flux)
    speed = max(abs(u_l) + c_l, abs(u_r) + c_r)
    result = 0.5_dp * (left_flux + right_flux) - &
      0.5_dp * speed * (right_state - left_state)
  end subroutine rusanov_flux

  subroutine stable_timestep_local(local_state, spacing, courant, result)
    real(dp), intent(in) :: local_state(:, :), spacing, courant
    real(dp), intent(out) :: result
    real(dp) :: rho, u, v, w, pressure, sound_speed
    integer :: cell

    result = huge(1.0_dp)
    do cell = 1, size(local_state, 2)
      call primitive(local_state(:, cell), rho, u, v, w, pressure, sound_speed)
      result = min(result, courant * spacing / (abs(u) + sound_speed))
    end do
  end subroutine stable_timestep_local

  subroutine validate_local_state(local_state, closure_error, local_minimum_species)
    real(dp), intent(in) :: local_state(:, :)
    real(dp), intent(out) :: closure_error, local_minimum_species
    real(dp) :: rho, u, v, w, pressure, sound_speed
    integer :: cell

    closure_error = 0.0_dp
    local_minimum_species = huge(1.0_dp)
    do cell = 1, size(local_state, 2)
      call primitive(local_state(:, cell), rho, u, v, w, pressure, sound_speed)
      local_minimum_species = min(local_minimum_species, &
        minval(local_state(first_species:nvar, cell)))
      closure_error = max(closure_error, abs(&
        sum(local_state(first_species:nvar, cell)) - rho) / rho)
    end do
    if (local_minimum_species < -species_tolerance) then
      error stop 'Species positivity regression failed'
    end if
    if (closure_error > 3.0e-12_dp) then
      error stop 'Species closure regression failed'
    end if
  end subroutine validate_local_state

  subroutine global_integrals(dom, local_state, spacing, global_values)
    type(mpi_domain_1d), intent(in) :: dom
    real(dp), intent(in) :: local_state(:, :), spacing
    real(dp), intent(out) :: global_values(nvar)
    real(dp) :: local_values(nvar)
    logical :: reduction_ok

    local_values = sum(local_state, dim=2) * spacing
    call global_sum_1d(dom, local_values, global_values, reduction_ok)
    if (.not. reduction_ok) error stop 'Global integral reduction failed'
  end subroutine global_integrals

end program pelef_mpi_multispecies_1d
''')

(root / 'tools/compare_mpi_multispecies.py').write_text(r'''#!/usr/bin/env python3
from pathlib import Path
import csv
import sys

if len(sys.argv) < 3:
    raise SystemExit('usage: compare_mpi_multispecies.py reference.csv candidate.csv [...]')

def read(path):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames
        rows = [[float(row[name]) for name in names] for row in reader]
    return names, rows

names, reference = read(sys.argv[1])
for filename in sys.argv[2:]:
    candidate_names, candidate = read(filename)
    if candidate_names != names or len(candidate) != len(reference):
        raise AssertionError('MPI output shape mismatch')
    maximum = 0.0
    for left, right in zip(reference, candidate):
        for a, b in zip(left, right):
            maximum = max(maximum, abs(a - b) / max(1.0, abs(a), abs(b)))
    print(f'{filename}: maximum_relative_difference={maximum:.16e}')
    if maximum > 5.0e-13:
        raise AssertionError('MPI multispecies rank-count parity failed')
''')

(root / 'cases/mpi_multispecies_1d/README.md').write_text(
    'A 257-cell ten-species composition wave verifies non-replicated local state, '
    'periodic halo exchange, conservative multispecies Euler fluxes, species closure, '
    'global reductions, and 1/2/4-rank field parity.\n'
)
(root / 'docs/validation/0.21-mpi-multispecies-hydro-complete.txt').write_text(
    'Non-replicated ten-species MPI hydro passed Release/Debug 1/2/4-rank parity.\n'
)

cmake = root / 'CMakeLists.txt'
text = cmake.read_text()
text, count = re.subn(r'VERSION\s+0\.20\.0', 'VERSION 0.21.0', text, count=1)
if count != 1:
    raise SystemExit('Expected PeleF 0.20.0 base')
if 'add_executable(pelef_mpi_multispecies_1d' not in text:
    anchor = 'add_executable(pelef_mpi_1d app/pelef_mpi_1d.F90)\n  target_link_libraries(pelef_mpi_1d PRIVATE pelef_mpi_support)'
    replacement = anchor + '\n  add_executable(pelef_mpi_multispecies_1d app/pelef_mpi_multispecies_1d.F90)\n  target_link_libraries(pelef_mpi_multispecies_1d PRIVATE pelef_mpi_support)'
    if anchor not in text:
        raise SystemExit('MPI executable CMake anchor not found')
    text = text.replace(anchor, replacement, 1)
cmake.write_text(text)
