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
  public :: load_toy_isomerization_thermo

contains

  subroutine load_gri30_thermo_subset(species, ok)
    type(nasa7_species), allocatable, intent(out) :: species(:)
    logical, intent(out) :: ok
    integer :: i

    allocate(species(4))

    ! NASA7 entries from Cantera gri30.yaml and air.yaml.
    species(gri_h2_index)%name = "H2"
    species(gri_h2_index)%molecular_weight = 2.01588_dp
    species(gri_h2_index)%temperature_min = 200.0_dp
    species(gri_h2_index)%temperature_mid = 1000.0_dp
    species(gri_h2_index)%temperature_max = 3500.0_dp
    species(gri_h2_index)%low_coefficients = [ &
      2.34433112_dp, 7.98052075e-3_dp, -1.94781510e-5_dp, &
      2.01572094e-8_dp, -7.37611761e-12_dp, -917.935173_dp, &
      0.683010238_dp ]
    species(gri_h2_index)%high_coefficients = [ &
      3.33727920_dp, -4.94024731e-5_dp, 4.99456778e-7_dp, &
      -1.79566394e-10_dp, 2.00255376e-14_dp, -950.158922_dp, &
      -3.20502331_dp ]

    species(gri_o2_index)%name = "O2"
    species(gri_o2_index)%molecular_weight = 31.9988_dp
    species(gri_o2_index)%temperature_min = 200.0_dp
    species(gri_o2_index)%temperature_mid = 1000.0_dp
    species(gri_o2_index)%temperature_max = 3500.0_dp
    species(gri_o2_index)%low_coefficients = [ &
      3.78245636_dp, -2.99673416e-3_dp, 9.84730201e-6_dp, &
      -9.68129509e-9_dp, 3.24372837e-12_dp, -1063.94356_dp, &
      3.65767573_dp ]
    species(gri_o2_index)%high_coefficients = [ &
      3.28253784_dp, 1.48308754e-3_dp, -7.57966669e-7_dp, &
      2.09470555e-10_dp, -2.16717794e-14_dp, -1088.45772_dp, &
      5.45323129_dp ]

    species(gri_h2o_index)%name = "H2O"
    species(gri_h2o_index)%molecular_weight = 18.01528_dp
    species(gri_h2o_index)%temperature_min = 200.0_dp
    species(gri_h2o_index)%temperature_mid = 1000.0_dp
    species(gri_h2o_index)%temperature_max = 3500.0_dp
    species(gri_h2o_index)%low_coefficients = [ &
      4.19864056_dp, -2.03643410e-3_dp, 6.52040211e-6_dp, &
      -5.48797062e-9_dp, 1.77197817e-12_dp, -3.02937267e4_dp, &
      -0.849032208_dp ]
    species(gri_h2o_index)%high_coefficients = [ &
      3.03399249_dp, 2.17691804e-3_dp, -1.64072518e-7_dp, &
      -9.70419870e-11_dp, 1.68200992e-14_dp, -3.00042971e4_dp, &
      4.96677010_dp ]

    species(gri_n2_index)%name = "N2"
    species(gri_n2_index)%molecular_weight = 28.0134_dp
    species(gri_n2_index)%temperature_min = 300.0_dp
    species(gri_n2_index)%temperature_mid = 1000.0_dp
    species(gri_n2_index)%temperature_max = 5000.0_dp
    species(gri_n2_index)%low_coefficients = [ &
      3.29867700_dp, 1.40824040e-3_dp, -3.96322200e-6_dp, &
      5.64151500e-9_dp, -2.44485400e-12_dp, -1020.89990_dp, &
      3.95037200_dp ]
    species(gri_n2_index)%high_coefficients = [ &
      2.92664000_dp, 1.48797680e-3_dp, -5.68476000e-7_dp, &
      1.00970380e-10_dp, -6.75335100e-15_dp, -922.797700_dp, &
      5.98052800_dp ]

    ok = .true.
    do i = 1, size(species)
      ok = ok .and. valid_nasa7_species(species(i))
    end do
  end subroutine load_gri30_thermo_subset

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

    ok = .true.
    do i = 1, size(species)
      ok = ok .and. valid_nasa7_species(species(i))
    end do
  end subroutine load_toy_isomerization_thermo

end module thermo_database_mod
