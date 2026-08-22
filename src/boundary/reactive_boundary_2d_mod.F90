module reactive_boundary_2d_mod
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
    real(dp) :: inflow_temperature = 300.0_dp
    real(dp), allocatable :: inflow_primitive(:)
  end type reactive_boundary_face_2d

  type, public :: reactive_boundary_set_2d
    type(reactive_boundary_face_2d) :: face(4)
  end type reactive_boundary_set_2d

  public :: build_reactive_boundary_set_2d
  public :: initialize_periodic_boundary_set_2d
  public :: validate_reactive_boundary_set_2d
  public :: sample_reactive_primitive_2d
  public :: reactive_boundary_is_periodic
  public :: reactive_boundary_is_wall
  public :: reactive_boundary_is_inflow
  public :: reactive_boundary_is_outflow

contains

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

  subroutine validate_reactive_boundary_set_2d(boundaries, ok)
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    logical, intent(out) :: ok
    integer :: side

    ok = .true.
    do side = 1, 4
      ok = ok .and. valid_boundary_kind(boundaries%face(side)%kind)
      ok = ok .and. valid_thermal_kind(boundaries%face(side)%thermal)
      ok = ok .and. boundaries%face(side)%wall_temperature > 0.0_dp
      ok = ok .and. allocated(boundaries%face(side)%inflow_primitive)
      ok = ok .and. boundaries%face(side)%inflow_temperature > 0.0_dp
    end do
    ok = ok .and. (reactive_boundary_is_periodic(boundaries%face(boundary_x_lower)) .eqv. &
      reactive_boundary_is_periodic(boundaries%face(boundary_x_upper)))
    ok = ok .and. (reactive_boundary_is_periodic(boundaries%face(boundary_y_lower)) .eqv. &
      reactive_boundary_is_periodic(boundaries%face(boundary_y_upper)))
  end subroutine validate_reactive_boundary_set_2d

  subroutine initialize_periodic_boundary_set_2d(nprimitive, boundaries)
    integer, intent(in) :: nprimitive
    type(reactive_boundary_set_2d), intent(out) :: boundaries
    integer :: side

    do side = 1, 4
      boundaries%face(side)%kind = "periodic"
      boundaries%face(side)%thermal = "adiabatic"
      boundaries%face(side)%wall_temperature = 300.0_dp
      boundaries%face(side)%wall_velocity = 0.0_dp
      boundaries%face(side)%inflow_temperature = 300.0_dp
      allocate(boundaries%face(side)%inflow_primitive(nprimitive))
      boundaries%face(side)%inflow_primitive = 0.0_dp
    end do
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
    do side = 1, 4
      boundaries%face(side)%inflow_temperature = config%initial_temperature
      allocate(boundaries%face(side)%inflow_primitive(size(primitive)))
      boundaries%face(side)%inflow_primitive = primitive
    end do
    call validate_reactive_boundary_set_2d(boundaries, ok)
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
