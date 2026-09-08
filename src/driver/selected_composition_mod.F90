module selected_composition_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  implicit none
  private

  integer, parameter, public :: selected_composition_max_species = 32
  integer, parameter, public :: selected_composition_name_length = 24

  public :: validate_selected_composition_fields
  public :: resolve_selected_composition

contains

  pure subroutine validate_selected_composition_fields( &
      context, composition_count, composition_species, &
      composition_mole_fractions, ok, message)
    character(len=*), intent(in) :: context
    integer, intent(in) :: composition_count
    character(len=*), intent(in) :: composition_species(:)
    real(dp), intent(in) :: composition_mole_fractions(:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: first_index, second_index

    ok = .false.
    message = ""
    if (size(composition_species) /= selected_composition_max_species .or. &
        size(composition_mole_fractions) /= &
          selected_composition_max_species) then
      message = trim(context) // &
        " selected composition storage has the wrong size"
      return
    end if
    if (composition_count < 1 .or. &
        composition_count > selected_composition_max_species) then
      message = trim(context) // " composition_count is out of range"
      return
    end if
    if (any(.not. ieee_is_finite(composition_mole_fractions))) then
      message = trim(context) // &
        " composition contains a nonfinite mole fraction"
      return
    end if
    if (any(len_trim(composition_species(1:composition_count)) == 0)) then
      message = trim(context) // " composition contains a blank species name"
      return
    end if
    if (any(composition_mole_fractions(1:composition_count) < 0.0_dp) .or. &
        maxval(composition_mole_fractions(1:composition_count)) <= 0.0_dp) then
      message = trim(context) // &
        " composition must be nonnegative and nonempty"
      return
    end if
    if (composition_count < selected_composition_max_species) then
      if (any(len_trim(composition_species(composition_count + 1:)) /= 0)) then
        message = trim(context) // &
          " composition has names beyond composition_count"
        return
      end if
      if (any(abs(composition_mole_fractions( &
          composition_count + 1:)) > 0.0_dp)) then
        message = trim(context) // &
          " composition has values beyond composition_count"
        return
      end if
    end if
    do first_index = 1, composition_count - 1
      do second_index = first_index + 1, composition_count
        if (trim(composition_species(first_index)) == &
            trim(composition_species(second_index))) then
          message = trim(context) // &
            " composition contains a duplicate species"
          return
        end if
      end do
    end do
    ok = .true.
  end subroutine validate_selected_composition_fields


  subroutine resolve_selected_composition( &
      context, composition_count, composition_species, &
      composition_mole_fractions, species, mole_fractions, ok, message)
    character(len=*), intent(in) :: context
    integer, intent(in) :: composition_count
    character(len=*), intent(in) :: composition_species(:)
    real(dp), intent(in) :: composition_mole_fractions(:)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(out) :: mole_fractions(:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    logical, allocatable :: assigned(:)
    integer :: input_index, species_index, match_index, match_count
    real(dp) :: scale, total

    mole_fractions = 0.0_dp
    call validate_selected_composition_fields( &
      context, composition_count, composition_species, &
      composition_mole_fractions, ok, message)
    if (.not. ok) return
    ok = .false.
    if (size(species) < 2 .or. &
        size(species) > selected_composition_max_species) then
      message = trim(context) // &
        " selected mechanism species count is unsupported"
      return
    end if
    if (size(mole_fractions) /= size(species)) then
      message = trim(context) // &
        " selected composition output has the wrong size"
      return
    end if
    if (composition_count > size(species)) then
      message = trim(context) // &
        " composition_count exceeds the selected mechanism"
      return
    end if

    allocate(assigned(size(species)))
    assigned = .false.
    do input_index = 1, composition_count
      match_index = 0
      match_count = 0
      do species_index = 1, size(species)
        if (trim(composition_species(input_index)) == &
            trim(species(species_index)%name)) then
          match_index = species_index
          match_count = match_count + 1
        end if
      end do
      if (match_count /= 1) then
        message = trim(context) // &
          " species was not found exactly once in the bundle: " // &
          trim(composition_species(input_index))
        return
      end if
      if (assigned(match_index)) then
        message = trim(context) // &
          " composition contains a duplicate species: " // &
          trim(composition_species(input_index))
        return
      end if
      assigned(match_index) = .true.
      mole_fractions(match_index) = composition_mole_fractions(input_index)
    end do

    scale = maxval(mole_fractions)
    if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) then
      message = trim(context) // &
        " selected composition has no positive total"
      return
    end if
    if (scale <= huge(1.0_dp) / real(size(mole_fractions), dp)) then
      total = sum(mole_fractions)
    else
      mole_fractions = mole_fractions / scale
      total = sum(mole_fractions)
    end if
    if (.not. ieee_is_finite(total) .or. total <= 0.0_dp) then
      message = trim(context) // &
        " selected composition normalization failed"
      return
    end if
    mole_fractions = mole_fractions / total
    message = ""
    ok = .true.
  end subroutine resolve_selected_composition

end module selected_composition_mod
