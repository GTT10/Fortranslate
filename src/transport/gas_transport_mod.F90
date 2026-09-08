module gas_transport_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  implicit none
  private

  integer, parameter :: transport_name_length = 24
  integer, parameter, public :: transport_geometry_atom = 0
  integer, parameter, public :: transport_geometry_linear = 1
  integer, parameter, public :: transport_geometry_nonlinear = 2

  type, public :: gas_transport_species
    character(len=transport_name_length) :: name = ""
    integer :: geometry = transport_geometry_atom
    real(dp) :: well_depth = 0.0_dp ! Lennard-Jones epsilon / k_B [K]
    real(dp) :: diameter = 0.0_dp ! Lennard-Jones sigma [angstrom]
    real(dp) :: dipole = 0.0_dp ! Debye; retained for provenance
    real(dp) :: polarizability = 0.0_dp ! angstrom^3; retained for provenance
    real(dp) :: rotational_relaxation = 0.0_dp
  end type gas_transport_species

  public :: valid_gas_transport_species
  public :: compatible_transport_database
  public :: load_gas_transport_data

contains

  pure logical function valid_gas_transport_species(record) result(valid)
    type(gas_transport_species), intent(in) :: record

    valid = .false.
    if (.not. all(ieee_is_finite([record%well_depth, record%diameter, &
        record%dipole, record%polarizability, &
        record%rotational_relaxation]))) return
    if (len_trim(record%name) == 0) return
    if (record%geometry < transport_geometry_atom) return
    if (record%geometry > transport_geometry_nonlinear) return
    if (record%well_depth <= 0.0_dp) return
    if (record%diameter <= 0.0_dp) return
    if (record%dipole < 0.0_dp) return
    if (record%polarizability < 0.0_dp) return
    if (record%rotational_relaxation < 0.0_dp) return
    valid = .true.
  end function valid_gas_transport_species

  pure logical function compatible_transport_database(species, transport) &
      result(compatible)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    integer :: k

    compatible = .false.
    if (size(species) /= size(transport)) return
    if (size(species) == 0) return
    do k = 1, size(species)
      if (.not. valid_gas_transport_species(transport(k))) return
      if (trim(species(k)%name) /= trim(transport(k)%name)) return
    end do
    compatible = .true.
  end function compatible_transport_database

  subroutine load_gas_transport_data( &
      names, geometries, values, transport, ok)
    character(len=*), intent(in) :: names(:)
    integer, intent(in) :: geometries(:)
    real(dp), intent(in) :: values(:, :)
    type(gas_transport_species), allocatable, intent(out) :: transport(:)
    logical, intent(out) :: ok
    integer :: k

    ok = size(names) > 0 .and. size(geometries) == size(names) .and. &
      size(values, 1) == 5 .and. size(values, 2) == size(names)
    if (.not. ok) return
    do k = 1, size(names)
      if (len_trim(names(k)) == 0 .or. &
          len_trim(names(k)) > transport_name_length) then
        ok = .false.
        return
      end if
    end do

    allocate(transport(size(names)))
    do k = 1, size(names)
      transport(k)%name = trim(names(k))
      transport(k)%geometry = geometries(k)
      transport(k)%well_depth = values(1, k)
      transport(k)%diameter = values(2, k)
      transport(k)%dipole = values(3, k)
      transport(k)%polarizability = values(4, k)
      transport(k)%rotational_relaxation = values(5, k)
      if (.not. valid_gas_transport_species(transport(k))) then
        deallocate(transport)
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine load_gas_transport_data

end module gas_transport_mod
