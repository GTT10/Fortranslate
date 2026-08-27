module reactive_boundary_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use simulation_config_reactive_2d_mod, only: &
    reactive_2d_config, reactive_2d_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nprim, reactive_mass_fraction_component
  implicit none
  private

  integer, parameter, public :: boundary_x_lower = 1
  integer, parameter, public :: boundary_x_upper = 2
  integer, parameter, public :: boundary_y_lower = 3
  integer, parameter, public :: boundary_y_upper = 4

  type, public :: reactive_boundary_face_2d
    character(len=24) :: kind = "periodic"
    character(len=24) :: thermal = "adiabatic"
    real(dp) :: wall_temperature = 300.0_dp
    real(dp) :: wall_velocity(3) = 0.0_dp
    character(len=24) :: wall_species = "impermeable"
    real(dp), allocatable :: prescribed_species_flux(:)
    real(dp) :: inflow_temperature = 300.0_dp
    real(dp), allocatable :: inflow_primitive(:)
  end type reactive_boundary_face_2d

  type, public :: reactive_boundary_set_2d
    type(reactive_boundary_face_2d) :: face(4)
    ! A single-valued embedded boundary is a wall, not a fifth Cartesian
    ! domain face.  Keeping its transport data with the existing boundary
    ! set lets serial, AMR, and MPI paths share the same validated control.
    type(reactive_boundary_face_2d) :: embedded_wall
  end type reactive_boundary_set_2d

  public :: build_reactive_boundary_set_2d
  public :: configure_reactive_embedded_wall_2d
  public :: initialize_periodic_boundary_set_2d
  public :: validate_reactive_boundary_set_2d
  public :: sample_reactive_primitive_2d
  public :: reactive_boundary_is_periodic
  public :: reactive_boundary_is_wall
  public :: reactive_boundary_is_inflow
  public :: reactive_boundary_is_outflow
  public :: reactive_boundary_has_prescribed_species_flux

contains

  subroutine configure_reactive_embedded_wall_2d( &
      boundaries, kind, thermal, wall_temperature, wall_velocity, ok)
    type(reactive_boundary_set_2d), intent(inout) :: boundaries
    character(len=*), intent(in) :: kind, thermal
    real(dp), intent(in) :: wall_temperature, wall_velocity(3)
    logical, intent(out) :: ok

    type(reactive_boundary_set_2d) :: candidate

    candidate = boundaries
    candidate%embedded_wall%kind = trim(kind)
    candidate%embedded_wall%thermal = trim(thermal)
    candidate%embedded_wall%wall_temperature = wall_temperature
    candidate%embedded_wall%wall_velocity = wall_velocity
    call validate_reactive_boundary_set_2d(candidate, ok)
    if (.not. ok) return
    boundaries = candidate
  end subroutine configure_reactive_embedded_wall_2d

  pure logical function reactive_boundary_is_periodic(face) result(is_periodic)
    type(reactive_boundary_face_2d), intent(in) :: face
    is_periodic = trim(face%kind) == "periodic"
  end function reactive_boundary_is_periodic

  pure logical function reactive_boundary_is_wall(face) result(is_wall)
    type(reactive_boundary_face_2d), intent(in) :: face
    is_wall = trim(face%kind) == "slip_wall" .or. &
      trim(face%kind) == "no_slip_wall"
  end function reactive_boundary_is_wall

  pure logical function reactive_boundary_is_inflow(face) result(is_inflow)
    type(reactive_boundary_face_2d), intent(in) :: face
    is_inflow = trim(face%kind) == "inflow"
  end function reactive_boundary_is_inflow

  pure logical function reactive_boundary_is_outflow(face) result(is_outflow)
    type(reactive_boundary_face_2d), intent(in) :: face
    is_outflow = trim(face%kind) == "outflow"
  end function reactive_boundary_is_outflow

  pure logical function reactive_boundary_has_prescribed_species_flux(face) &
      result(has_prescribed_flux)
    type(reactive_boundary_face_2d), intent(in) :: face
    has_prescribed_flux = reactive_boundary_is_wall(face) .and. &
      trim(face%wall_species) == "prescribed"
  end function reactive_boundary_has_prescribed_species_flux

  pure logical function valid_boundary_kind(kind) result(valid)
    character(len=*), intent(in) :: kind
    valid = trim(kind) == "periodic" .or. trim(kind) == "slip_wall" .or. &
      trim(kind) == "no_slip_wall" .or. trim(kind) == "inflow" .or. &
      trim(kind) == "outflow"
  end function valid_boundary_kind

  pure logical function valid_thermal_kind(kind) result(valid)
    character(len=*), intent(in) :: kind
    valid = trim(kind) == "adiabatic" .or. trim(kind) == "isothermal"
  end function valid_thermal_kind

  pure logical function valid_wall_species_kind(kind) result(valid)
    character(len=*), intent(in) :: kind
    valid = trim(kind) == "impermeable" .or. trim(kind) == "prescribed"
  end function valid_wall_species_kind

  subroutine validate_reactive_boundary_set_2d(boundaries, ok)
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    logical, intent(out) :: ok
    real(dp) :: scale, tolerance
    integer :: side

    ok = .true.
    do side = 1, 4
      ok = ok .and. valid_boundary_kind(boundaries%face(side)%kind)
      ok = ok .and. valid_thermal_kind(boundaries%face(side)%thermal)
      ok = ok .and. boundaries%face(side)%wall_temperature > 0.0_dp
      ok = ok .and. allocated(boundaries%face(side)%inflow_primitive)
      ok = ok .and. allocated(boundaries%face(side)%prescribed_species_flux)
      ok = ok .and. boundaries%face(side)%inflow_temperature > 0.0_dp
      ok = ok .and. valid_wall_species_kind( &
        boundaries%face(side)%wall_species)
      if (.not. ok) return
      ok = ok .and. size(boundaries%face(side)%inflow_primitive) >= 6
      ok = ok .and. size(boundaries%face(side)%prescribed_species_flux) == &
        size(boundaries%face(side)%inflow_primitive) - 5
      ok = ok .and. all(ieee_is_finite( &
        boundaries%face(side)%prescribed_species_flux))
      if (.not. ok) return
      scale = max(1.0_dp, maxval(abs( &
        boundaries%face(side)%prescribed_species_flux)))
      tolerance = 2.0e3_dp * epsilon(1.0_dp) * scale
      if (trim(boundaries%face(side)%wall_species) == "impermeable") then
        ok = maxval(abs(boundaries%face(side)%prescribed_species_flux)) <= &
          tolerance
      else
        ok = reactive_boundary_is_wall(boundaries%face(side)) .and. &
          abs(sum(boundaries%face(side)%prescribed_species_flux)) <= tolerance
      end if
      if (.not. ok) return
    end do
    ok = reactive_boundary_is_wall(boundaries%embedded_wall) .and. &
      valid_thermal_kind(boundaries%embedded_wall%thermal) .and. &
      boundaries%embedded_wall%wall_temperature > 0.0_dp .and. &
      all(ieee_is_finite([boundaries%embedded_wall%wall_temperature, &
        boundaries%embedded_wall%wall_velocity])) .and. &
      boundaries%embedded_wall%inflow_temperature > 0.0_dp .and. &
      allocated(boundaries%embedded_wall%inflow_primitive) .and. &
      allocated(boundaries%embedded_wall%prescribed_species_flux) .and. &
      trim(boundaries%embedded_wall%wall_species) == "impermeable"
    if (.not. ok) return
    ok = size(boundaries%embedded_wall%inflow_primitive) >= 6 .and. &
      size(boundaries%embedded_wall%prescribed_species_flux) == &
        size(boundaries%embedded_wall%inflow_primitive) - 5 .and. &
      all(ieee_is_finite( &
        boundaries%embedded_wall%prescribed_species_flux))
    if (.not. ok) return
    scale = max(1.0_dp, maxval(abs( &
      boundaries%embedded_wall%prescribed_species_flux)))
    tolerance = 2.0e3_dp * epsilon(1.0_dp) * scale
    ok = maxval(abs( &
      boundaries%embedded_wall%prescribed_species_flux)) <= tolerance
    if (.not. ok) return
    ok = ok .and. (reactive_boundary_is_periodic(boundaries%face(boundary_x_lower)) .eqv. &
      reactive_boundary_is_periodic(boundaries%face(boundary_x_upper)))
    ok = ok .and. (reactive_boundary_is_periodic(boundaries%face(boundary_y_lower)) .eqv. &
      reactive_boundary_is_periodic(boundaries%face(boundary_y_upper)))
  end subroutine validate_reactive_boundary_set_2d

  subroutine initialize_periodic_boundary_set_2d(nprimitive, boundaries)
    integer, intent(in) :: nprimitive
    type(reactive_boundary_set_2d), intent(out) :: boundaries
    integer :: side, nspecies

    nspecies = max(1, nprimitive - 5)

    do side = 1, 4
      boundaries%face(side)%kind = "periodic"
      boundaries%face(side)%thermal = "adiabatic"
      boundaries%face(side)%wall_temperature = 300.0_dp
      boundaries%face(side)%wall_velocity = 0.0_dp
      boundaries%face(side)%wall_species = "impermeable"
      allocate(boundaries%face(side)%prescribed_species_flux(nspecies))
      boundaries%face(side)%prescribed_species_flux = 0.0_dp
      boundaries%face(side)%inflow_temperature = 300.0_dp
      allocate(boundaries%face(side)%inflow_primitive(nprimitive))
      boundaries%face(side)%inflow_primitive = 0.0_dp
    end do
    boundaries%embedded_wall%kind = "slip_wall"
    boundaries%embedded_wall%thermal = "adiabatic"
    boundaries%embedded_wall%wall_temperature = 300.0_dp
    boundaries%embedded_wall%wall_velocity = 0.0_dp
    boundaries%embedded_wall%wall_species = "impermeable"
    allocate(boundaries%embedded_wall%prescribed_species_flux(nspecies))
    boundaries%embedded_wall%prescribed_species_flux = 0.0_dp
    boundaries%embedded_wall%inflow_temperature = 300.0_dp
    allocate(boundaries%embedded_wall%inflow_primitive(nprimitive))
    boundaries%embedded_wall%inflow_primitive = 0.0_dp
  end subroutine initialize_periodic_boundary_set_2d

  subroutine build_reactive_boundary_set_2d(species, config, boundaries, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_2d_config), intent(in) :: config
    type(reactive_boundary_set_2d), intent(out) :: boundaries
    logical, intent(out) :: ok

    real(dp), allocatable :: mass_fractions(:), primitive(:)
    real(dp), allocatable :: mole_fractions(:)
    real(dp) :: density
    logical :: local_ok
    integer :: side, k

    ok = .false.
    allocate(mass_fractions(size(species)), mole_fractions(size(species)), &
      primitive(reactive_nprim(size(species))))
    call reactive_2d_mole_fractions(config, size(species), mole_fractions, local_ok)
    if (.not. local_ok) return
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, local_ok)
    if (.not. local_ok) return
    density = mixture_density( &
      species, mass_fractions, config%initial_pressure, &
      config%initial_temperature, local_ok)
    if (.not. local_ok) return
    primitive = 0.0_dp
    primitive(1:5) = [density, config%initial_velocity_x, &
      config%initial_velocity_y, 0.0_dp, config%initial_pressure]
    do k = 1, size(species)
      primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
    end do

    boundaries%face(boundary_x_lower)%kind = trim(config%boundary_x_lower)
    boundaries%face(boundary_x_upper)%kind = trim(config%boundary_x_upper)
    boundaries%face(boundary_y_lower)%kind = trim(config%boundary_y_lower)
    boundaries%face(boundary_y_upper)%kind = trim(config%boundary_y_upper)
    boundaries%face(boundary_x_lower)%thermal = trim(config%thermal_x_lower)
    boundaries%face(boundary_x_upper)%thermal = trim(config%thermal_x_upper)
    boundaries%face(boundary_y_lower)%thermal = trim(config%thermal_y_lower)
    boundaries%face(boundary_y_upper)%thermal = trim(config%thermal_y_upper)
    boundaries%face(boundary_x_lower)%wall_temperature = config%wall_temperature_x_lower
    boundaries%face(boundary_x_upper)%wall_temperature = config%wall_temperature_x_upper
    boundaries%face(boundary_y_lower)%wall_temperature = config%wall_temperature_y_lower
    boundaries%face(boundary_y_upper)%wall_temperature = config%wall_temperature_y_upper
    boundaries%face(boundary_x_lower)%wall_velocity = config%wall_velocity_x_lower
    boundaries%face(boundary_x_upper)%wall_velocity = config%wall_velocity_x_upper
    boundaries%face(boundary_y_lower)%wall_velocity = config%wall_velocity_y_lower
    boundaries%face(boundary_y_upper)%wall_velocity = config%wall_velocity_y_upper
    boundaries%face(boundary_x_lower)%wall_species = &
      trim(config%wall_species_x_lower)
    boundaries%face(boundary_x_upper)%wall_species = &
      trim(config%wall_species_x_upper)
    boundaries%face(boundary_y_lower)%wall_species = &
      trim(config%wall_species_y_lower)
    boundaries%face(boundary_y_upper)%wall_species = &
      trim(config%wall_species_y_upper)
    do side = 1, 4
      boundaries%face(side)%inflow_temperature = config%initial_temperature
      allocate(boundaries%face(side)%inflow_primitive(size(primitive)))
      boundaries%face(side)%inflow_primitive = primitive
      allocate(boundaries%face(side)%prescribed_species_flux(size(species)))
      boundaries%face(side)%prescribed_species_flux = 0.0_dp
    end do
    boundaries%embedded_wall%kind = "slip_wall"
    boundaries%embedded_wall%thermal = "adiabatic"
    boundaries%embedded_wall%wall_temperature = 300.0_dp
    boundaries%embedded_wall%wall_velocity = 0.0_dp
    boundaries%embedded_wall%wall_species = "impermeable"
    boundaries%embedded_wall%inflow_temperature = config%initial_temperature
    allocate(boundaries%embedded_wall%inflow_primitive(size(primitive)))
    boundaries%embedded_wall%inflow_primitive = primitive
    allocate(boundaries%embedded_wall%prescribed_species_flux(size(species)))
    boundaries%embedded_wall%prescribed_species_flux = 0.0_dp
    boundaries%face(boundary_x_lower)%prescribed_species_flux = &
      config%prescribed_species_flux_x_lower(1:size(species))
    boundaries%face(boundary_x_upper)%prescribed_species_flux = &
      config%prescribed_species_flux_x_upper(1:size(species))
    boundaries%face(boundary_y_lower)%prescribed_species_flux = &
      config%prescribed_species_flux_y_lower(1:size(species))
    boundaries%face(boundary_y_upper)%prescribed_species_flux = &
      config%prescribed_species_flux_y_upper(1:size(species))
    call validate_reactive_boundary_set_2d(boundaries, ok)
    if (.not. ok) return
    if ((reactive_boundary_has_prescribed_species_flux( &
          boundaries%face(boundary_x_lower)) .or. &
        reactive_boundary_has_prescribed_species_flux( &
          boundaries%face(boundary_x_upper)) .or. &
        reactive_boundary_has_prescribed_species_flux( &
          boundaries%face(boundary_y_lower)) .or. &
        reactive_boundary_has_prescribed_species_flux( &
          boundaries%face(boundary_y_upper))) .and. &
        (.not. config%transport_enabled .or. &
         .not. config%species_diffusion_enabled)) ok = .false.
  end subroutine build_reactive_boundary_set_2d

  recursive subroutine sample_reactive_primitive_2d( &
      primitive, temperature, nx, ny, i, j, boundaries, sampled, &
      sampled_temperature, ok)
    real(dp), intent(in) :: primitive(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny, i, j
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), intent(out) :: sampled(:), sampled_temperature
    logical, intent(out) :: ok

    type(reactive_boundary_face_2d) :: face
    real(dp), allocatable :: base(:)
    real(dp) :: base_temperature
    integer :: mapped_i, mapped_j
    logical :: local_ok

    sampled = 0.0_dp
    sampled_temperature = 0.0_dp
    ok = .false.
    if (size(sampled) /= size(primitive, 1) .or. nx < 1 .or. ny < 1) return
    if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
      sampled = primitive(:, i, j)
      sampled_temperature = temperature(i, j)
      ok = .true.
      return
    end if

    allocate(base(size(sampled)))
    if (i < 1) then
      face = boundaries%face(boundary_x_lower)
      if (reactive_boundary_is_periodic(face)) then
        mapped_i = 1 + modulo(i - 1, nx)
        call sample_reactive_primitive_2d(primitive, temperature, nx, ny, mapped_i, j, boundaries, sampled, sampled_temperature, ok)
        return
      else if (reactive_boundary_is_inflow(face)) then
        sampled = face%inflow_primitive
        sampled_temperature = face%inflow_temperature
        ok = .true.
        return
      else if (reactive_boundary_is_outflow(face)) then
        call sample_reactive_primitive_2d(primitive, temperature, nx, ny, 1, j, boundaries, sampled, sampled_temperature, ok)
        return
      end if
      mapped_i = min(nx, max(1, 1 - i))
      call sample_reactive_primitive_2d(primitive, temperature, nx, ny, mapped_i, j, boundaries, base, base_temperature, local_ok)
      if (.not. local_ok) return
      sampled = base
      sampled(2) = 2.0_dp * face%wall_velocity(1) - base(2)
      if (trim(face%kind) == "no_slip_wall") then
        sampled(3) = 2.0_dp * face%wall_velocity(2) - base(3)
        sampled(4) = 2.0_dp * face%wall_velocity(3) - base(4)
      end if
      if (trim(face%thermal) == "isothermal") then
        sampled_temperature = max(tiny(1.0_dp), 2.0_dp * face%wall_temperature - base_temperature)
      else
        sampled_temperature = base_temperature
      end if
      ok = .true.
      return
    else if (i > nx) then
      face = boundaries%face(boundary_x_upper)
      if (reactive_boundary_is_periodic(face)) then
        mapped_i = 1 + modulo(i - 1, nx)
        call sample_reactive_primitive_2d(primitive, temperature, nx, ny, mapped_i, j, boundaries, sampled, sampled_temperature, ok)
        return
      else if (reactive_boundary_is_inflow(face)) then
        sampled = face%inflow_primitive
        sampled_temperature = face%inflow_temperature
        ok = .true.
        return
      else if (reactive_boundary_is_outflow(face)) then
        call sample_reactive_primitive_2d(primitive, temperature, nx, ny, nx, j, boundaries, sampled, sampled_temperature, ok)
        return
      end if
      mapped_i = min(nx, max(1, 2 * nx + 1 - i))
      call sample_reactive_primitive_2d(primitive, temperature, nx, ny, mapped_i, j, boundaries, base, base_temperature, local_ok)
      if (.not. local_ok) return
      sampled = base
      sampled(2) = 2.0_dp * face%wall_velocity(1) - base(2)
      if (trim(face%kind) == "no_slip_wall") then
        sampled(3) = 2.0_dp * face%wall_velocity(2) - base(3)
        sampled(4) = 2.0_dp * face%wall_velocity(3) - base(4)
      end if
      if (trim(face%thermal) == "isothermal") then
        sampled_temperature = max(tiny(1.0_dp), 2.0_dp * face%wall_temperature - base_temperature)
      else
        sampled_temperature = base_temperature
      end if
      ok = .true.
      return
    end if

    if (j < 1) then
      face = boundaries%face(boundary_y_lower)
      if (reactive_boundary_is_periodic(face)) then
        mapped_j = 1 + modulo(j - 1, ny)
        call sample_reactive_primitive_2d(primitive, temperature, nx, ny, i, mapped_j, boundaries, sampled, sampled_temperature, ok)
        return
      else if (reactive_boundary_is_inflow(face)) then
        sampled = face%inflow_primitive
        sampled_temperature = face%inflow_temperature
        ok = .true.
        return
      else if (reactive_boundary_is_outflow(face)) then
        call sample_reactive_primitive_2d(primitive, temperature, nx, ny, i, 1, boundaries, sampled, sampled_temperature, ok)
        return
      end if
      mapped_j = min(ny, max(1, 1 - j))
      call sample_reactive_primitive_2d(primitive, temperature, nx, ny, i, mapped_j, boundaries, base, base_temperature, local_ok)
      if (.not. local_ok) return
      sampled = base
      sampled(3) = 2.0_dp * face%wall_velocity(2) - base(3)
      if (trim(face%kind) == "no_slip_wall") then
        sampled(2) = 2.0_dp * face%wall_velocity(1) - base(2)
        sampled(4) = 2.0_dp * face%wall_velocity(3) - base(4)
      end if
      if (trim(face%thermal) == "isothermal") then
        sampled_temperature = max(tiny(1.0_dp), 2.0_dp * face%wall_temperature - base_temperature)
      else
        sampled_temperature = base_temperature
      end if
      ok = .true.
      return
    else if (j > ny) then
      face = boundaries%face(boundary_y_upper)
      if (reactive_boundary_is_periodic(face)) then
        mapped_j = 1 + modulo(j - 1, ny)
        call sample_reactive_primitive_2d(primitive, temperature, nx, ny, i, mapped_j, boundaries, sampled, sampled_temperature, ok)
        return
      else if (reactive_boundary_is_inflow(face)) then
        sampled = face%inflow_primitive
        sampled_temperature = face%inflow_temperature
        ok = .true.
        return
      else if (reactive_boundary_is_outflow(face)) then
        call sample_reactive_primitive_2d(primitive, temperature, nx, ny, i, ny, boundaries, sampled, sampled_temperature, ok)
        return
      end if
      mapped_j = min(ny, max(1, 2 * ny + 1 - j))
      call sample_reactive_primitive_2d(primitive, temperature, nx, ny, i, mapped_j, boundaries, base, base_temperature, local_ok)
      if (.not. local_ok) return
      sampled = base
      sampled(3) = 2.0_dp * face%wall_velocity(2) - base(3)
      if (trim(face%kind) == "no_slip_wall") then
        sampled(2) = 2.0_dp * face%wall_velocity(1) - base(2)
        sampled(4) = 2.0_dp * face%wall_velocity(3) - base(4)
      end if
      if (trim(face%thermal) == "isothermal") then
        sampled_temperature = max(tiny(1.0_dp), 2.0_dp * face%wall_temperature - base_temperature)
      else
        sampled_temperature = base_temperature
      end if
      ok = .true.
    end if
  end subroutine sample_reactive_primitive_2d

end module reactive_boundary_2d_mod
