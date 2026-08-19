module h2o2_full_thermo_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species, valid_nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  implicit none
  private

  integer, parameter, public :: full_h2_index = 1
  integer, parameter, public :: full_h_index = 2
  integer, parameter, public :: full_o_index = 3
  integer, parameter, public :: full_o2_index = 4
  integer, parameter, public :: full_oh_index = 5
  integer, parameter, public :: full_h2o_index = 6
  integer, parameter, public :: full_ho2_index = 7
  integer, parameter, public :: full_h2o2_index = 8
  integer, parameter, public :: full_ar_index = 9
  integer, parameter, public :: full_n2_index = 10
  integer, parameter, public :: full_nspecies = 10

  public :: load_h2o2_full_thermo

contains

  subroutine load_h2o2_full_thermo(species, ok)
    type(nasa7_species), allocatable, intent(out) :: species(:)
    logical, intent(out) :: ok
    type(nasa7_species), allocatable :: elementary_species(:)
    integer :: index

    call load_h2o2_elementary_thermo(elementary_species, ok)
    if (.not. ok .or. size(elementary_species) /= 7) then
      ok = .false.
      return
    end if

    allocate(species(full_nspecies))
    species(full_h2_index) = elementary_species(1)
    species(full_h_index) = elementary_species(2)
    species(full_o_index) = elementary_species(3)
    species(full_o2_index) = elementary_species(4)
    species(full_oh_index) = elementary_species(5)
    species(full_h2o_index) = elementary_species(6)
    call set_ho2(species(full_ho2_index))
    call set_h2o2(species(full_h2o2_index))
    call set_ar(species(full_ar_index))
    species(full_n2_index) = elementary_species(7)

    ok = .true.
    do index = 1, size(species)
      if (.not. valid_nasa7_species(species(index))) then
        ok = .false.
        return
      end if
    end do
  end subroutine load_h2o2_full_thermo

  subroutine set_ho2(species)
    type(nasa7_species), intent(out) :: species

    species%name = "HO2"
    species%molecular_weight = 33.006_dp
    species%temperature_min = 200.0_dp
    species%temperature_mid = 1000.0_dp
    species%temperature_max = 3500.0_dp
    species%low_coefficients = [ &
      4.30179801_dp, -4.74912051e-3_dp, 2.11582891e-5_dp, &
      -2.42763894e-8_dp, 9.29225124e-12_dp, 294.808040_dp, &
      3.71666245_dp ]
    species%high_coefficients = [ &
      4.01721090_dp, 2.23982013e-3_dp, -6.33658150e-7_dp, &
      1.14246370e-10_dp, -1.07908535e-14_dp, 111.856713_dp, &
      3.78510215_dp ]
  end subroutine set_ho2

  subroutine set_h2o2(species)
    type(nasa7_species), intent(out) :: species

    species%name = "H2O2"
    species%molecular_weight = 34.014_dp
    species%temperature_min = 200.0_dp
    species%temperature_mid = 1000.0_dp
    species%temperature_max = 3500.0_dp
    species%low_coefficients = [ &
      4.27611269_dp, -5.42822417e-4_dp, 1.67335701e-5_dp, &
      -2.15770813e-8_dp, 8.62454363e-12_dp, -1.77025821e4_dp, &
      3.43505074_dp ]
    species%high_coefficients = [ &
      4.16500285_dp, 4.90831694e-3_dp, -1.90139225e-6_dp, &
      3.71185986e-10_dp, -2.87908305e-14_dp, -1.78617877e4_dp, &
      2.91615662_dp ]
  end subroutine set_h2o2

  subroutine set_ar(species)
    type(nasa7_species), intent(out) :: species

    species%name = "AR"
    species%molecular_weight = 39.950_dp
    species%temperature_min = 300.0_dp
    species%temperature_mid = 1000.0_dp
    species%temperature_max = 5000.0_dp
    species%low_coefficients = [ &
      2.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -745.375_dp, 4.366_dp ]
    species%high_coefficients = species%low_coefficients
  end subroutine set_ar

end module h2o2_full_thermo_mod
