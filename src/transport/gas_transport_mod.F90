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
    real(dp) :: simple_viscosity(4)=0, simple_conductivity(4)=0
    real(dp), allocatable :: simple_binary(:,:) ! Four log-T coefficients per partner.
  end type gas_transport_species

  public :: valid_gas_transport_species
  public :: compatible_transport_database
  public :: load_gas_transport_data
  public :: set_simple_transport_fits,write_simple_transport_fits,check_simple_transport_fits
  public :: write_transport_context
  public :: transport_numeric_values

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
    if(allocated(record%simple_binary)) then
      if(size(record%simple_binary,1)/=4.or.size(record%simple_binary,2)<2) return
      if(.not.all(ieee_is_finite(record%simple_binary))) return
      if(.not.all(ieee_is_finite([record%simple_viscosity,record%simple_conductivity]))) return
      if(maxval(abs(record%simple_binary))>1e6_dp) return
      if(maxval(abs([record%simple_viscosity,record%simple_conductivity]))>1e6_dp) return
    end if
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
      if(allocated(transport(k)%simple_binary).neqv.allocated(transport(1)%simple_binary)) return
      if(allocated(transport(k)%simple_binary)) then
        if(size(transport(k)%simple_binary,2)/=size(species)) return
      end if
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

  subroutine set_simple_transport_fits(transport,coefficients,ok)
    type(gas_transport_species), intent(inout) :: transport(:)
    real(dp), intent(in) :: coefficients(:,:)
    logical, intent(out) :: ok
    integer :: n,i,j,k
    ok=.false.;n=size(transport)
    if(n<2.or.any(shape(coefficients)/=[4,2*n+n*(n+1)/2])) return
    if(.not.all(ieee_is_finite(coefficients))) return
    if(maxval(abs(coefficients))>1e6_dp) return
    do i=1,n
      if(.not.valid_gas_transport_species(transport(i))) return
    end do
    do i=1,n
      transport(i)%simple_viscosity=coefficients(:,i)
      transport(i)%simple_conductivity=coefficients(:,n+i)
      if(allocated(transport(i)%simple_binary)) deallocate(transport(i)%simple_binary)
      allocate(transport(i)%simple_binary(4,n))
    end do
    do j=1,n
      do i=1,j
        k=2*n+j*(j-1)/2+i
        transport(i)%simple_binary(:,j)=coefficients(:,k)
        transport(j)%simple_binary(:,i)=coefficients(:,k)
      end do
    end do
    ok=.true.
  end subroutine

  pure function transport_numeric_values(record,n) result(values)
    type(gas_transport_species), intent(in) :: record
    integer, intent(in) :: n
    real(dp) :: values(14+4*n)
    values=0
    values(:5)=[record%well_depth,record%diameter,record%dipole,record%polarizability,record%rotational_relaxation]
    if(.not.allocated(record%simple_binary)) return
    values(6)=-1
    if(any(shape(record%simple_binary)/=[4,n])) return
    values(6)=1
    values(7:14)=[record%simple_viscosity,record%simple_conductivity]
    values(15:)=reshape(record%simple_binary,[4*n])
  end function

  subroutine write_transport_context(unit,record,status)
    integer, intent(in) :: unit
    type(gas_transport_species), intent(in) :: record
    integer, intent(out) :: status
    ! Preserve the legacy list-directed representation when fits are absent.
    write(unit,*,iostat=status) record%name,record%geometry,record%well_depth,record%diameter, &
      record%dipole,record%polarizability,record%rotational_relaxation
    if(status/=0) return
    call write_simple_transport_fits(unit,record,status)
  end subroutine

  subroutine write_simple_transport_fits(unit,record,status)
    integer, intent(in) :: unit
    type(gas_transport_species), intent(in) :: record
    integer, intent(out) :: status
    integer :: j
    status=0
    if(.not.allocated(record%simple_binary)) return
    write(unit,'(a)',iostat=status) 'PELEF_SIMPLE_TRANSPORT_V1'
    if(status/=0) return
    write(unit,'(*(es27.18e3,1x))',iostat=status) record%simple_viscosity,record%simple_conductivity
    if(status/=0) return
    do j=1,size(record%simple_binary,2)
      write(unit,'(*(es27.18e3,1x))',iostat=status) record%simple_binary(:,j)
      if(status/=0) return
    end do
  end subroutine

  subroutine check_simple_transport_fits(unit,record,status)
    integer, intent(in) :: unit
    type(gas_transport_species), intent(in) :: record
    integer, intent(out) :: status
    real(dp) :: values(8),binary(4)
    character(len=64) :: marker
    integer :: j
    status=0
    if(.not.allocated(record%simple_binary)) return
    read(unit,'(a)',iostat=status) marker
    if(status/=0) return
    status=1
    if(marker/='PELEF_SIMPLE_TRANSPORT_V1') return
    read(unit,*,iostat=status) values
    if(status/=0) return
    status=1
    if(.not.all(ieee_is_finite(values))) return
    if(any(values/=[record%simple_viscosity,record%simple_conductivity])) return
    do j=1,size(record%simple_binary,2)
      read(unit,*,iostat=status) binary
      if(status/=0) return
      status=1
      if(.not.all(ieee_is_finite(binary))) return
      if(any(binary/=record%simple_binary(:,j))) return
    end do
    status=0
  end subroutine

end module gas_transport_mod
