module transport_database_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  implicit none
  private

  integer, parameter, public :: transport_geometry_atom = 0
  integer, parameter, public :: transport_geometry_linear = 1
  integer, parameter, public :: transport_geometry_nonlinear = 2

  type, public :: gas_transport_species
    character(len=24) :: name = ""
    integer :: geometry = transport_geometry_atom
    real(dp) :: well_depth = 0.0_dp ! Lennard-Jones epsilon / k_B [K]
    real(dp) :: diameter = 0.0_dp ! Lennard-Jones sigma [angstrom]
    real(dp) :: dipole = 0.0_dp ! Debye; retained for provenance
    real(dp) :: polarizability = 0.0_dp ! angstrom^3; retained for provenance
    real(dp) :: rotational_relaxation = 0.0_dp
  end type gas_transport_species

  public :: valid_gas_transport_species
  public :: compatible_transport_database
  public :: load_h2o2_elementary_transport
  public :: load_h2o2_full_transport

contains

  logical function valid_gas_transport_species(record) result(valid)
    type(gas_transport_species), intent(in) :: record

    valid = len_trim(record%name) > 0 .and. &
      record%geometry >= transport_geometry_atom .and. &
      record%geometry <= transport_geometry_nonlinear .and. &
      record%well_depth > 0.0_dp .and. record%diameter > 0.0_dp .and. &
      record%dipole >= 0.0_dp .and. record%polarizability >= 0.0_dp .and. &
      record%rotational_relaxation >= 0.0_dp .and. &
      all(ieee_is_finite([record%well_depth, record%diameter, &
        record%dipole, record%polarizability, &
        record%rotational_relaxation]))
  end function valid_gas_transport_species

  logical function compatible_transport_database(species, transport) &
      result(compatible)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    integer :: k

    compatible = size(species) == size(transport) .and. size(species) > 0
    if (.not. compatible) return
    do k = 1, size(species)
      compatible = compatible .and. valid_gas_transport_species(transport(k))
      compatible = compatible .and. &
        trim(species(k)%name) == trim(transport(k)%name)
    end do
  end function compatible_transport_database

  subroutine load_h2o2_elementary_transport(transport, ok)
    type(gas_transport_species), allocatable, intent(out) :: transport(:)
    logical, intent(out) :: ok

    ! Order matches load_h2o2_elementary_thermo and h2o2_elementary.json.
    ! Parameters are pinned to Cantera data/h2o2.yaml at commit
    ! 11a2381011cb6d42e61cc4c195e0f920864bf8d3.
    allocate(transport(7))

    call set_record(transport(1), "H2", transport_geometry_linear, &
      38.0_dp, 2.92_dp, 0.0_dp, 0.79_dp, 280.0_dp)
    call set_record(transport(2), "H", transport_geometry_atom, &
      145.0_dp, 2.05_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    call set_record(transport(3), "O", transport_geometry_atom, &
      80.0_dp, 2.75_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    call set_record(transport(4), "O2", transport_geometry_linear, &
      107.4_dp, 3.458_dp, 0.0_dp, 1.6_dp, 3.8_dp)
    call set_record(transport(5), "OH", transport_geometry_linear, &
      80.0_dp, 2.75_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    call set_record(transport(6), "H2O", transport_geometry_nonlinear, &
      572.4_dp, 2.605_dp, 1.844_dp, 0.0_dp, 4.0_dp)
    call set_record(transport(7), "N2", transport_geometry_linear, &
      97.53_dp, 3.621_dp, 0.0_dp, 1.76_dp, 4.0_dp)

    ok = all_valid(transport)
  end subroutine load_h2o2_elementary_transport


  subroutine load_h2o2_full_transport(transport, ok)
    type(gas_transport_species), allocatable, intent(out) :: transport(:)
    logical, intent(out) :: ok

    ! Order matches load_h2o2_full_thermo and mechanisms/h2o2_full.json.
    allocate(transport(10))
    call set_record(transport(1), "H2", transport_geometry_linear, &
      38.0_dp, 2.92_dp, 0.0_dp, 0.79_dp, 280.0_dp)
    call set_record(transport(2), "H", transport_geometry_atom, &
      145.0_dp, 2.05_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    call set_record(transport(3), "O", transport_geometry_atom, &
      80.0_dp, 2.75_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    call set_record(transport(4), "O2", transport_geometry_linear, &
      107.4_dp, 3.458_dp, 0.0_dp, 1.6_dp, 3.8_dp)
    call set_record(transport(5), "OH", transport_geometry_linear, &
      80.0_dp, 2.75_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    call set_record(transport(6), "H2O", transport_geometry_nonlinear, &
      572.4_dp, 2.605_dp, 1.844_dp, 0.0_dp, 4.0_dp)
    call set_record(transport(7), "HO2", transport_geometry_nonlinear, &
      107.4_dp, 3.458_dp, 0.0_dp, 0.0_dp, 1.0_dp)
    call set_record(transport(8), "H2O2", transport_geometry_nonlinear, &
      107.4_dp, 3.458_dp, 0.0_dp, 0.0_dp, 3.8_dp)
    call set_record(transport(9), "AR", transport_geometry_atom, &
      136.5_dp, 3.33_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    call set_record(transport(10), "N2", transport_geometry_linear, &
      97.53_dp, 3.621_dp, 0.0_dp, 1.76_dp, 4.0_dp)
    ok = all_valid(transport)
  end subroutine load_h2o2_full_transport

  subroutine set_record(record, name, geometry, well_depth, diameter, dipole, &
      polarizability, rotational_relaxation)
    type(gas_transport_species), intent(out) :: record
    character(len=*), intent(in) :: name
    integer, intent(in) :: geometry
    real(dp), intent(in) :: well_depth, diameter, dipole
    real(dp), intent(in) :: polarizability, rotational_relaxation

    record%name = trim(name)
    record%geometry = geometry
    record%well_depth = well_depth
    record%diameter = diameter
    record%dipole = dipole
    record%polarizability = polarizability
    record%rotational_relaxation = rotational_relaxation
  end subroutine set_record

  logical function all_valid(transport) result(valid)
    type(gas_transport_species), intent(in) :: transport(:)
    integer :: k

    valid = size(transport) > 0
    do k = 1, size(transport)
      valid = valid .and. valid_gas_transport_species(transport(k))
    end do
  end function all_valid

end module transport_database_mod
