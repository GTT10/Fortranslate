program pelef_transport_probe
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use mixture_transport_mod, only: mixture_transport_coefficients
  implicit none

  integer, parameter :: nstates = 4
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp) :: temperatures(nstates), pressures(nstates)
  real(dp) :: x(7, nstates), y(7), mu, lambda, diffusion(7)
  character(len=24) :: labels(nstates)
  character(len=1024) :: output_path
  logical :: ok
  integer :: unit, status, i, k

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load transport-probe thermodynamics"
  call load_h2o2_elementary_transport(transport, ok)
  if (.not. ok) error stop "Failed to load transport-probe database"

  labels = [character(len=24) :: &
    "stoich_600K", "stoich_1000K", "product_1500K", "lean_1800K_5atm"]
  temperatures = [600.0_dp, 1000.0_dp, 1500.0_dp, 1800.0_dp]
  pressures = [101325.0_dp, 101325.0_dp, 101325.0_dp, 5.0_dp * 101325.0_dp]

  x(:, 1) = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  x(:, 2) = x(:, 1)
  x(:, 3) = [0.04_dp, 1.0e-5_dp, 1.0e-5_dp, 0.06_dp, 1.0e-5_dp, &
    0.30_dp, 0.59997_dp]
  x(:, 4) = [0.10_dp, 1.0e-5_dp, 1.0e-5_dp, 0.05_dp, 1.0e-5_dp, &
    0.02_dp, 0.82997_dp]

  output_path = "transport_probe.csv"
  if (command_argument_count() >= 1) call get_command_argument(1, output_path)
  open(newunit=unit, file=trim(output_path), status="replace", action="write", &
    iostat=status)
  if (status /= 0) error stop "Could not open transport-probe output"
  write(unit, '(a)') &
    "label,temperature,pressure,viscosity,thermal_conductivity," // &
    "D_H2,D_H,D_O,D_O2,D_OH,D_H2O,D_N2," // &
    "Y_H2,Y_H,Y_O,Y_O2,Y_OH,Y_H2O,Y_N2"

  do i = 1, nstates
    x(:, i) = x(:, i) / sum(x(:, i))
    call mass_fractions_from_mole_fractions(species, x(:, i), y, ok)
    if (.not. ok) error stop "Invalid transport-probe composition"
    call mixture_transport_coefficients( &
      species, transport, y, temperatures(i), pressures(i), mu, lambda, &
      diffusion, ok)
    if (.not. ok) error stop "Transport-probe coefficient evaluation failed"
    write(unit, '(a,",",*(es25.16e3,:,","))') trim(labels(i)), &
      temperatures(i), pressures(i), mu, lambda, (diffusion(k), k = 1, 7), &
      (y(k), k = 1, 7)
  end do
  close(unit)
end program pelef_transport_probe
