module transport_database_mod
  use precision_mod, only: dp
  use gas_transport_mod, only: &
    gas_transport_species, transport_geometry_atom, &
    transport_geometry_linear, transport_geometry_nonlinear, &
    valid_gas_transport_species, compatible_transport_database, &
    load_gas_transport_data
  use h2o2_full_mechanism_mod, only: &
    h2o2_full_nspecies, load_h2o2_full_transport_data
  implicit none
  private

  public :: gas_transport_species
  public :: transport_geometry_atom
  public :: transport_geometry_linear
  public :: transport_geometry_nonlinear
  public :: valid_gas_transport_species
  public :: compatible_transport_database
  public :: load_gas_transport_data
  public :: load_h2o2_elementary_transport
  public :: load_h2o2_full_transport

contains

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
    character(len=24) :: names(h2o2_full_nspecies)
    integer :: geometries(h2o2_full_nspecies)
    real(dp) :: values(5, h2o2_full_nspecies)

    call load_h2o2_full_transport_data(names, geometries, values, ok)
    if (.not. ok) return
    call load_gas_transport_data(names, geometries, values, transport, ok)
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
