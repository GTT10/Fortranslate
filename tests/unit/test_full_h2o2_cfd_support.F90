program test_full_h2o2_cfd_support
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: &
    h2o2_full_nspecies, h2o2_full_nreactions, load_h2o2_full_mechanism
  use elementary_kinetics_mod, only: &
    elementary_reaction, reaction_kind_three_body, reaction_kind_falloff
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_full_transport, &
    compatible_transport_database
  use mixture_transport_mod, only: &
    mixture_transport_coefficients, standard_atmosphere
  use simulation_config_reactive_1d_mod, only: &
    reactive_1d_config, reactive_1d_mole_fractions
  use simulation_config_reactive_2d_mod, only: &
    reactive_2d_config, reactive_2d_mole_fractions
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, build_reactive_boundary_set_2d
  use reactive_1d_mod, only: reactive_nvar
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_boundary_set_2d) :: boundaries
  type(reactive_1d_config) :: config_1d
  type(reactive_2d_config) :: config_2d
  real(dp) :: x(h2o2_full_nspecies), y(h2o2_full_nspecies)
  real(dp) :: x_2d(h2o2_full_nspecies)
  real(dp) :: diffusion(h2o2_full_nspecies), viscosity, conductivity
  logical :: ok
  integer :: k, three_body_count, falloff_count

  call load_h2o2_full_thermo(species, ok)
  call require(ok .and. size(species) == h2o2_full_nspecies, &
    "full thermodynamics load")
  call load_h2o2_full_mechanism(reactions, ok)
  call require(ok .and. size(reactions) == h2o2_full_nreactions, &
    "full mechanism load")
  call load_h2o2_full_transport(transport, ok)
  call require(ok, "full transport load")
  call require(compatible_transport_database(species, transport), &
    "full thermo/transport ordering")
  call require(reactive_nvar(size(species)) == 15, &
    "ten-species conserved-state width")

  three_body_count = 0
  falloff_count = 0
  do k = 1, size(reactions)
    if (reactions(k)%kind == reaction_kind_three_body) &
      three_body_count = three_body_count + 1
    if (reactions(k)%kind == reaction_kind_falloff) &
      falloff_count = falloff_count + 1
  end do
  call require(three_body_count > 0, "third-body reactions present")
  call require(falloff_count > 0, "falloff reactions present")

  config_1d%chemistry_model = "full_h2o2"
  call reactive_1d_mole_fractions(config_1d, size(species), x, ok)
  call require(ok, "full 1D composition dispatch")
  call mass_fractions_from_mole_fractions(species, x, y, ok)
  call require(ok, "full mole-to-mass conversion")
  call mixture_transport_coefficients( &
    species, transport, y, 1000.0_dp, standard_atmosphere, viscosity, &
    conductivity, diffusion, ok)
  call require(ok, "full mixture transport evaluation")
  call require(viscosity > 0.0_dp .and. conductivity > 0.0_dp, &
    "full viscosity and conductivity positive")
  call require(minval(diffusion) > 0.0_dp, &
    "full mixture-averaged diffusion positive")

  config_2d%chemistry_model = "full_h2o2"
  call reactive_2d_mole_fractions(config_2d, size(species), x_2d, ok)
  call require(ok, "full 2D composition dispatch")
  call require(maxval(abs(x_2d - x)) < 1.0e-15_dp, &
    "1D/2D full composition agreement")
  call build_reactive_boundary_set_2d(species, config_2d, boundaries, ok)
  call require(ok, "full species physical-boundary construction")

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // trim(message)
      error stop 1
    end if
  end subroutine require

end program test_full_h2o2_cfd_support
