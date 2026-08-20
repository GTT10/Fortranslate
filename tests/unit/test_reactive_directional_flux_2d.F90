program test_reactive_directional_flux_2d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_species_component, reactive_primitive_to_conserved
  use reactive_2d_mod, only: reactive_riemann_flux_y
  use state_indices_mod, only: irho, imx, imy, imz, iet
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: primitive(:), state(:), flux(:), mass_fractions(:)
  real(dp) :: mole_fractions(7), temperature, sound_speed, mass_flux, tolerance
  logical :: ok
  integer :: k

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "thermodynamic database load failed"
  allocate(primitive(reactive_nprim(size(species))))
  allocate(state(reactive_nvar(size(species))), flux(reactive_nvar(size(species))))
  allocate(mass_fractions(size(species)))
  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  if (.not. ok) error stop "composition conversion failed"

  primitive(1:5) = [0.31_dp, 42.0_dp, -27.0_dp, 6.0_dp, 135000.0_dp]
  do k = 1, size(species)
    primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state, temperature, sound_speed, ok)
  if (.not. ok) error stop "state construction failed"
  call reactive_riemann_flux_y( &
    species, state, state, temperature, temperature, "hllc", flux, ok)
  if (.not. ok) error stop "equal-state y HLLC flux failed"

  mass_flux = state(irho) * primitive(3)
  tolerance = 2.0e-11_dp * max(1.0_dp, maxval(abs(flux)))
  if (abs(flux(irho) - mass_flux) > tolerance) error stop "mass flux mismatch"
  if (abs(flux(imx) - mass_flux * primitive(2)) > tolerance) &
    error stop "x-momentum placement mismatch"
  if (abs(flux(imy) - (mass_flux * primitive(3) + primitive(5))) > tolerance) &
    error stop "normal momentum placement mismatch"
  if (abs(flux(imz) - mass_flux * primitive(4)) > tolerance) &
    error stop "z-momentum placement mismatch"
  if (abs(flux(iet) - (state(iet) + primitive(5)) * primitive(3)) > tolerance) &
    error stop "energy flux mismatch"
  if (abs(sum(flux(reactive_species_component(1): &
      reactive_species_component(size(species)))) - flux(irho)) > tolerance) &
    error stop "species flux closure mismatch"
end program test_reactive_directional_flux_2d
