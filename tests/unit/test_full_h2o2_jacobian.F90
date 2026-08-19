program test_full_h2o2_jacobian
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions, mixture_density
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_thermo_mod, only: full_nspecies, load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: &
    load_h2o2_full_mechanism, h2o2_full_mass_fraction_jacobian
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: x(full_nspecies), y(full_nspecies)
  real(dp) :: jacobian(full_nspecies, full_nspecies), density
  logical :: ok

  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop "Failed to load full thermo"
  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load full mechanism"
  x = [2.0_dp, 1.0e-5_dp, 2.0e-6_dp, 1.0_dp, 3.0e-6_dp, &
    1.0e-4_dp, 1.0e-5_dp, 1.0e-6_dp, 0.1_dp, 3.0_dp]
  x = x / sum(x)
  call mass_fractions_from_mole_fractions(species, x, y, ok)
  if (.not. ok) error stop "Composition conversion failed"
  density = mixture_density(species, y, 202650.0_dp, 1100.0_dp, ok)
  if (.not. ok) error stop "Density evaluation failed"
  call h2o2_full_mass_fraction_jacobian( &
    species, reactions, 1100.0_dp, density, y, jacobian, ok)
  if (.not. ok) error stop "Full mass-fraction Jacobian failed"
  if (maxval(abs(jacobian)) <= 1.0_dp) then
    error stop "Full Jacobian did not capture stiff scales"
  end if
  if (any(jacobian /= jacobian)) error stop "Full Jacobian contains NaN"
  write(*, '(a,es24.16)') "test_full_h2o2_jacobian: PASS, max|J|=", &
    maxval(abs(jacobian))
end program test_full_h2o2_jacobian
