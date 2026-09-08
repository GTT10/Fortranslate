module selected_mechanism_runtime_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species, valid_nasa7_species
  use elementary_kinetics_mod, only: &
    elementary_reaction, valid_elementary_reaction
  use gas_transport_mod, only: &
    gas_transport_species, load_gas_transport_data, &
    compatible_transport_database
  implicit none
  private

  public :: prepare_selected_mechanism
  public :: validate_selected_temperature_range

contains

  subroutine prepare_selected_mechanism( &
      context, expected_species_count, expected_reaction_count, species, &
      reactions, transport_names, transport_geometries, transport_values, &
      transport, common_temperature_minimum, common_temperature_maximum, &
      ok, message)
    character(len=*), intent(in) :: context
    integer, intent(in) :: expected_species_count, expected_reaction_count
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    character(len=*), intent(in) :: transport_names(:)
    integer, intent(in) :: transport_geometries(:)
    real(dp), intent(in) :: transport_values(:, :)
    type(gas_transport_species), allocatable, intent(out) :: transport(:)
    real(dp), intent(out) :: common_temperature_minimum
    real(dp), intent(out) :: common_temperature_maximum
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp) :: reactant_mass, product_mass, mass_scale
    integer :: first_index, second_index, reaction_index

    ok = .false.
    message = ""
    common_temperature_minimum = 0.0_dp
    common_temperature_maximum = 0.0_dp
    if (expected_species_count < 2 .or. expected_species_count > 32 .or. &
        size(species) /= expected_species_count) then
      message = trim(context) // &
        " mechanism returned an unsupported species count"
      return
    end if
    if (expected_reaction_count < 1 .or. &
        size(reactions) /= expected_reaction_count) then
      message = trim(context) // " mechanism returned the wrong reaction count"
      return
    end if
    do first_index = 1, size(species)
      if (.not. valid_nasa7_species(species(first_index)) .or. &
          len_trim(species(first_index)%name) == 0) then
        message = trim(context) // &
          " mechanism returned invalid thermodynamic data"
        return
      end if
      do second_index = first_index + 1, size(species)
        if (trim(species(first_index)%name) == &
            trim(species(second_index)%name)) then
          message = trim(context) // &
            " mechanism returned duplicate species names"
          return
        end if
      end do
    end do
    do reaction_index = 1, size(reactions)
      if (.not. valid_elementary_reaction( &
          reactions(reaction_index), size(species))) then
        message = trim(context) // " mechanism returned an invalid reaction"
        return
      end if
      reactant_mass = sum( &
        species%molecular_weight * &
          reactions(reaction_index)%reactant_stoich)
      product_mass = sum( &
        species%molecular_weight * &
          reactions(reaction_index)%product_stoich)
      mass_scale = max(1.0_dp, abs(reactant_mass), abs(product_mass))
      if (.not. ieee_is_finite(reactant_mass) .or. &
          .not. ieee_is_finite(product_mass) .or. &
          abs(product_mass - reactant_mass) > 2.0e-12_dp * mass_scale) then
        message = trim(context) // " reaction is not mass balanced"
        return
      end if
    end do

    call load_gas_transport_data( &
      transport_names, transport_geometries, transport_values, transport, ok)
    if (.not. ok .or. .not. allocated(transport)) then
      ok = .false.
      message = trim(context) // " transport database construction failed"
      return
    end if
    if (.not. compatible_transport_database(species, transport)) then
      ok = .false.
      message = trim(context) // &
        " transport order does not match thermodynamics"
      return
    end if

    common_temperature_minimum = maxval(species%temperature_min)
    common_temperature_maximum = minval(species%temperature_max)
    if (.not. ieee_is_finite(common_temperature_minimum) .or. &
        .not. ieee_is_finite(common_temperature_maximum) .or. &
        common_temperature_minimum >= common_temperature_maximum) then
      ok = .false.
      message = trim(context) // &
        " mechanism has no common thermodynamic interval"
      return
    end if
    message = ""
    ok = .true.
  end subroutine prepare_selected_mechanism


  pure subroutine validate_selected_temperature_range( &
      context, requested_temperature_minimum, requested_temperature_maximum, &
      common_temperature_minimum, common_temperature_maximum, ok, message)
    character(len=*), intent(in) :: context
    real(dp), intent(in) :: requested_temperature_minimum
    real(dp), intent(in) :: requested_temperature_maximum
    real(dp), intent(in) :: common_temperature_minimum
    real(dp), intent(in) :: common_temperature_maximum
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    ok = all(ieee_is_finite([ &
      requested_temperature_minimum, requested_temperature_maximum, &
      common_temperature_minimum, common_temperature_maximum])) .and. &
      requested_temperature_minimum <= requested_temperature_maximum .and. &
      common_temperature_minimum < common_temperature_maximum .and. &
      requested_temperature_minimum >= common_temperature_minimum .and. &
      requested_temperature_maximum <= common_temperature_maximum
    if (ok) then
      message = ""
    else
      message = trim(context) // &
        " temperature is outside the mechanism range"
    end if
  end subroutine validate_selected_temperature_range

end module selected_mechanism_runtime_mod
