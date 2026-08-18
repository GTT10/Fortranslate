! This file is generated from mechanisms/h2o2_full.json.
module h2o2_full_jacobian_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use pressure_dependent_kinetics_mod, only: pressure_dependent_reaction
  use full_reactor_rhs_mod, only: evaluate_full_reactor_rhs
  implicit none
  private

  integer, parameter, public :: h2o2_full_active_species = 9
  public :: evaluate_h2o2_full_jacobian

contains

  subroutine evaluate_h2o2_full_jacobian( &
      species, reactions, density, target_internal_energy, mass_fractions, &
      temperature_guess, derivative, jacobian, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(pressure_dependent_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy
    real(dp), intent(in) :: mass_fractions(:), temperature_guess
    real(dp), intent(out) :: derivative(:)
    real(dp), intent(out) :: jacobian(:, :)
    real(dp), intent(out) :: temperature
    logical, intent(out) :: ok
    real(dp), allocatable :: perturbed(:), perturbed_derivative(:)
    real(dp) :: perturbed_temperature, step
    logical :: local_ok
    integer :: column

    derivative = 0.0_dp
    jacobian = 0.0_dp
    temperature = temperature_guess
    ok = size(species) == size(mass_fractions)
    ok = ok .and. size(species) == h2o2_full_active_species + 1
    ok = ok .and. size(derivative) == size(species)
    ok = ok .and. size(jacobian, 1) == h2o2_full_active_species
    ok = ok .and. size(jacobian, 2) == h2o2_full_active_species
    if (.not. ok) return

    allocate(perturbed(size(species)), perturbed_derivative(size(species)))
    call evaluate_full_reactor_rhs( &
      species, reactions, density, target_internal_energy, mass_fractions, &
      temperature_guess, derivative, temperature, ok=local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if

    do column = 1, h2o2_full_active_species
      step = sqrt(epsilon(1.0_dp)) * &
        max(abs(mass_fractions(column)), 1.0e-7_dp)
      step = max(step, 1.0e-12_dp)
      perturbed = mass_fractions
      if (perturbed(size(species)) > 2.0_dp * step) then
        perturbed(column) = perturbed(column) + step
        perturbed(size(species)) = perturbed(size(species)) - step
        call evaluate_full_reactor_rhs( &
          species, reactions, density, target_internal_energy, perturbed, &
          temperature, perturbed_derivative, perturbed_temperature, ok=local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        jacobian(:, column) = &
          (perturbed_derivative(1:h2o2_full_active_species) - &
           derivative(1:h2o2_full_active_species)) / step
      else if (perturbed(column) > 2.0_dp * step) then
        perturbed(column) = perturbed(column) - step
        perturbed(size(species)) = perturbed(size(species)) + step
        call evaluate_full_reactor_rhs( &
          species, reactions, density, target_internal_energy, perturbed, &
          temperature, perturbed_derivative, perturbed_temperature, ok=local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        jacobian(:, column) = &
          (derivative(1:h2o2_full_active_species) - &
           perturbed_derivative(1:h2o2_full_active_species)) / step
      else
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine evaluate_h2o2_full_jacobian

end module h2o2_full_jacobian_mod
