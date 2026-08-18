module thermo_database_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species, valid_nasa7_species
  implicit none
  private

  integer, parameter, public :: gri_h2_index = 1
  integer, parameter, public :: gri_o2_index = 2
  integer, parameter, public :: gri_h2o_index = 3
  integer, parameter, public :: gri_n2_index = 4

  public :: load_gri30_thermo_subset
  public :: load_h2o2_elementary_thermo
  public :: load_h2o2_full_thermo
  public :: load_toy_isomerization_thermo

contains

  subroutine load_gri30_thermo_subset(species, ok)
    type(nasa7_species), allocatable, intent(out) :: species(:)
    logical, intent(out) :: ok

    allocate(species(4))
    call set_h2(species(gri_h2_index))
    call set_o2(species(gri_o2_index))
    call set_h2o(species(gri_h2o_index))
    call set_n2(species(gri_n2_index))
    ok = all_valid(species)
  end subroutine load_gri30_thermo_subset

  subroutine load_h2o2_elementary_thermo(species, ok)
    type(nasa7_species), allocatable, intent(out) :: species(:)
    logical, intent(out) :: ok

    ! Order matches mechanisms/h2o2_elementary.json.
    allocate(species(7))
    call set_h2(species(1))
    call set_h(species(2))
    call set_o(species(3))
    call set_o2(species(4))
    call set_oh(species(5))
    call set_h2o(species(6))
    call set_n2(species(7))
    ok = all_valid(species)
  end subroutine load_h2o2_elementary_thermo

  subroutine load_h2o2_full_thermo(species, ok)
    type(nasa7_species), allocatable, intent(out) :: species(:)
    logical, intent(out) :: ok

    ! Order matches mechanisms/h2o2_full.json and Cantera h2o2.yaml.
    allocate(species(10))
    call set_h2(species(1))
    call set_h(species(2))
    call set_o(species(3))
    call set_o2(species(4))
    call set_oh(species(5))
    call set_h2o(species(6))
    call set_ho2(species(7))
    call set_h2o2(species(8))
    call set_ar(species(9))
    call set_n2(species(10))
    ok = all_valid(species)
  end subroutine load_h2o2_full_thermo

  subroutine load_toy_isomerization_thermo(species, ok)
    type(nasa7_species), allocatable, intent(out) :: species(:)
    logical, intent(out) :: ok
    integer :: i

    allocate(species(2))
    do i = 1, 2
      species(i)%molecular_weight = 28.0_dp
      species(i)%temperature_min = 200.0_dp
      species(i)%temperature_mid = 1000.0_dp
      species(i)%temperature_max = 5000.0_dp
      species(i)%low_coefficients = [ &
        3.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp ]
      species(i)%high_coefficients = species(i)%low_coefficients
    end do

    species(1)%name = "A"
    species(2)%name = "B"
    species(2)%low_coefficients(6) = -3000.0_dp
    species(2)%high_coefficients(6) = -3000.0_dp
    ok = all_valid(species)
  end subroutine load_toy_isomerization_thermo

  logical function all_valid(species) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    integer :: i

    valid = size(species) > 0
    do i = 1, size(species)
      if (.not. valid_nasa7_species(species(i))) then
        valid = .false.
        return
      end if
    end do
  end function all_valid

  subroutine set_h2(species)
    type(nasa7_species), intent(out) :: species

    species%name = "H2"
    species%molecular_weight = 2.016_dp
    species%temperature_min = 200.0_dp
    species%temperature_mid = 1000.0_dp
    species%temperature_max = 3500.0_dp
    species%low_coefficients = [ &
      2.34433112_dp, 7.98052075e-3_dp, -1.94781510e-5_dp, &
      2.01572094e-8_dp, -7.37611761e-12_dp, -917.935173_dp, &
      0.683010238_dp ]
    species%high_coefficients = [ &
      3.33727920_dp, -4.94024731e-5_dp, 4.99456778e-7_dp, &
      -1.79566394e-10_dp, 2.00255376e-14_dp, -950.158922_dp, &
      -3.20502331_dp ]
  end subroutine set_h2

  subroutine set_h(species)
    type(nasa7_species), intent(out) :: species

    species%name = "H"
    species%molecular_weight = 1.008_dp
    species%temperature_min = 200.0_dp
    species%temperature_mid = 1000.0_dp
    species%temperature_max = 3500.0_dp
    species%low_coefficients = [ &
      2.5_dp, 7.05332819e-13_dp, -1.99591964e-15_dp, &
      2.30081632e-18_dp, -9.27732332e-22_dp, 2.54736599e4_dp, &
      -0.446682853_dp ]
    species%high_coefficients = [ &
      2.50000001_dp, -2.30842973e-11_dp, 1.61561948e-14_dp, &
      -4.73515235e-18_dp, 4.98197357e-22_dp, 2.54736599e4_dp, &
      -0.446682914_dp ]
  end subroutine set_h

  subroutine set_o(species)
    type(nasa7_species), intent(out) :: species

    species%name = "O"
    species%molecular_weight = 15.999_dp
    species%temperature_min = 200.0_dp
    species%temperature_mid = 1000.0_dp
    species%temperature_max = 3500.0_dp
    species%low_coefficients = [ &
      3.16826710_dp, -3.27931884e-3_dp, 6.64306396e-6_dp, &
      -6.12806624e-9_dp, 2.11265971e-12_dp, 2.91222592e4_dp, &
      2.05193346_dp ]
    species%high_coefficients = [ &
      2.56942078_dp, -8.59741137e-5_dp, 4.19484589e-8_dp, &
      -1.00177799e-11_dp, 1.22833691e-15_dp, 2.92175791e4_dp, &
      4.78433864_dp ]
  end subroutine set_o

  subroutine set_o2(species)
    type(nasa7_species), intent(out) :: species

    species%name = "O2"
    species%molecular_weight = 31.998_dp
    species%temperature_min = 200.0_dp
    species%temperature_mid = 1000.0_dp
    species%temperature_max = 3500.0_dp
    species%low_coefficients = [ &
      3.78245636_dp, -2.99673416e-3_dp, 9.84730201e-6_dp, &
      -9.68129509e-9_dp, 3.24372837e-12_dp, -1063.94356_dp, &
      3.65767573_dp ]
    species%high_coefficients = [ &
      3.28253784_dp, 1.48308754e-3_dp, -7.57966669e-7_dp, &
      2.09470555e-10_dp, -2.16717794e-14_dp, -1088.45772_dp, &
      5.45323129_dp ]
  end subroutine set_o2

  subroutine set_oh(species)
    type(nasa7_species), intent(out) :: species

    species%name = "OH"
    species%molecular_weight = 17.007_dp
    species%temperature_min = 200.0_dp
    species%temperature_mid = 1000.0_dp
    species%temperature_max = 3500.0_dp
    species%low_coefficients = [ &
      3.99201543_dp, -2.40131752e-3_dp, 4.61793841e-6_dp, &
      -3.88113333e-9_dp, 1.36411470e-12_dp, 3615.08056_dp, &
      -0.103925458_dp ]
    species%high_coefficients = [ &
      3.09288767_dp, 5.48429716e-4_dp, 1.26505228e-7_dp, &
      -8.79461556e-11_dp, 1.17412376e-14_dp, 3858.65700_dp, &
      4.47669610_dp ]
  end subroutine set_oh

  subroutine set_h2o(species)
    type(nasa7_species), intent(out) :: species

    species%name = "H2O"
    species%molecular_weight = 18.015_dp
    species%temperature_min = 200.0_dp
    species%temperature_mid = 1000.0_dp
    species%temperature_max = 3500.0_dp
    species%low_coefficients = [ &
      4.19864056_dp, -2.03643410e-3_dp, 6.52040211e-6_dp, &
      -5.48797062e-9_dp, 1.77197817e-12_dp, -3.02937267e4_dp, &
      -0.849032208_dp ]
    species%high_coefficients = [ &
      3.03399249_dp, 2.17691804e-3_dp, -1.64072518e-7_dp, &
      -9.70419870e-11_dp, 1.68200992e-14_dp, -3.00042971e4_dp, &
      4.96677010_dp ]
  end subroutine set_h2o

  subroutine set_ho2(species)
    type(nasa7_species), intent(out) :: species

    species%name = "HO2"
    species%molecular_weight = 33.006_dp
    species%temperature_min = 200.0_dp
    species%temperature_mid = 1000.0_dp
    species%temperature_max = 3500.0_dp
    species%low_coefficients = [ &
      4.30179801_dp, -4.74912051e-3_dp, 2.11582891e-5_dp, &
      -2.42763894e-8_dp, 9.29225124e-12_dp, 294.80804_dp, &
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

  subroutine set_n2(species)
    type(nasa7_species), intent(out) :: species

    species%name = "N2"
    species%molecular_weight = 28.014_dp
    species%temperature_min = 300.0_dp
    species%temperature_mid = 1000.0_dp
    species%temperature_max = 5000.0_dp
    species%low_coefficients = [ &
      3.29867700_dp, 1.40824040e-3_dp, -3.96322200e-6_dp, &
      5.64151500e-9_dp, -2.44485400e-12_dp, -1020.89990_dp, &
      3.95037200_dp ]
    species%high_coefficients = [ &
      2.92664000_dp, 1.48797680e-3_dp, -5.68476000e-7_dp, &
      1.00970380e-10_dp, -6.75335100e-15_dp, -922.797700_dp, &
      5.98052800_dp ]
  end subroutine set_n2

end module thermo_database_mod
