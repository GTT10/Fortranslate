module reactive_entropy_wave_3d_problem_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved
  use simulation_config_reactive_3d_mod, only: &
    reactive_3d_config, reactive_3d_mole_fractions
  implicit none
  private

  real(dp), parameter :: pi = acos(-1.0_dp)

  public :: initialize_reactive_entropy_wave_3d
  public :: initialize_reactive_problem_3d
  public :: reactive_entropy_wave_density_3d

contains

  subroutine initialize_reactive_problem_3d( &
      species, config, x, y, z, state, temperature, base_density, &
      mass_fractions, ok, base_mole_fractions)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_3d_config), intent(in) :: config
    real(dp), intent(in) :: x(:), y(:), z(:)
    real(dp), intent(out) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(out) :: base_density, mass_fractions(:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: base_mole_fractions(:)

    real(dp), allocatable :: mole_fractions(:), primitive(:)
    real(dp) :: density, sound_speed, local_temperature
    real(dp) :: delta_x, delta_y, delta_z, radius_squared
    logical :: cell_ok
    integer :: i, j, k, species_index, nspecies

    state = 0.0_dp
    temperature = 0.0_dp
    base_density = 0.0_dp
    mass_fractions = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (nspecies < 1 .or. size(mass_fractions) /= nspecies) return
    if (size(x) /= config%nx .or. size(y) /= config%ny .or. &
        size(z) /= config%nz) return
    if (size(state, 1) /= reactive_nvar(nspecies) .or. &
        size(state, 2) /= config%nx .or. size(state, 3) /= config%ny .or. &
        size(state, 4) /= config%nz) return
    if (size(temperature, 1) /= config%nx .or. &
        size(temperature, 2) /= config%ny .or. &
        size(temperature, 3) /= config%nz) return

    allocate(mole_fractions(nspecies), primitive(reactive_nprim(nspecies)))
    if (present(base_mole_fractions)) then
      if (size(base_mole_fractions) /= nspecies) return
      if (.not. all(ieee_is_finite(base_mole_fractions)) .or. &
          minval(base_mole_fractions) < 0.0_dp .or. &
          abs(sum(base_mole_fractions) - 1.0_dp) > 5.0e-10_dp) return
      mole_fractions = base_mole_fractions
    else
      call reactive_3d_mole_fractions( &
        config, nspecies, mole_fractions, cell_ok)
      if (.not. cell_ok) return
    end if
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, cell_ok)
    if (.not. cell_ok) return
    base_density = mixture_density( &
      species, mass_fractions, config%initial_pressure, &
      config%initial_temperature, cell_ok)
    if (.not. cell_ok) return

    do k = 1, config%nz
      do j = 1, config%ny
        do i = 1, config%nx
          select case (trim(config%problem))
          case ("entropy_wave")
            density = reactive_entropy_wave_density_3d( &
              x(i), y(j), z(k), 0.0_dp, config, base_density)
          case ("uniform_reactor")
            density = base_density
          case ("reactive_hotspot")
            delta_x = periodic_displacement( &
              x(i), config%hotspot_center_x, &
              config%x_upper - config%x_lower)
            delta_y = periodic_displacement( &
              y(j), config%hotspot_center_y, &
              config%y_upper - config%y_lower)
            delta_z = periodic_displacement( &
              z(k), config%hotspot_center_z, &
              config%z_upper - config%z_lower)
            radius_squared = delta_x**2 + delta_y**2 + delta_z**2
            local_temperature = config%initial_temperature + &
              config%hotspot_temperature_rise * exp( &
                -0.5_dp * radius_squared / config%hotspot_width**2)
            density = mixture_density( &
              species, mass_fractions, config%initial_pressure, &
              local_temperature, cell_ok)
            if (.not. cell_ok) return
          case default
            return
          end select
          primitive(1:5) = [ &
            density, config%initial_velocity_x, config%initial_velocity_y, &
            config%initial_velocity_z, config%initial_pressure]
          do species_index = 1, nspecies
            primitive(reactive_mass_fraction_component(species_index)) = &
              mass_fractions(species_index)
          end do
          call reactive_primitive_to_conserved( &
            species, primitive, state(:, i, j, k), temperature(i, j, k), &
            sound_speed, cell_ok)
          if (.not. cell_ok) return
        end do
      end do
    end do
    ok = .true.
  end subroutine initialize_reactive_problem_3d

  subroutine initialize_reactive_entropy_wave_3d( &
      species, config, x, y, z, state, temperature, base_density, &
      mass_fractions, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_3d_config), intent(in) :: config
    real(dp), intent(in) :: x(:), y(:), z(:)
    real(dp), intent(out) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(out) :: base_density, mass_fractions(:)
    logical, intent(out) :: ok

    if (trim(config%problem) /= "entropy_wave") then
      state = 0.0_dp
      temperature = 0.0_dp
      base_density = 0.0_dp
      mass_fractions = 0.0_dp
      ok = .false.
      return
    end if
    call initialize_reactive_problem_3d( &
      species, config, x, y, z, state, temperature, base_density, &
      mass_fractions, ok)
  end subroutine initialize_reactive_entropy_wave_3d

  pure real(dp) function periodic_displacement(value, center, length) &
      result(delta)
    real(dp), intent(in) :: value, center, length

    delta = value - center
    delta = delta - length * real(nint(delta / length), dp)
  end function periodic_displacement

  pure real(dp) function reactive_entropy_wave_density_3d( &
      x, y, z, time, config, base_density) result(density)
    real(dp), intent(in) :: x, y, z, time, base_density
    type(reactive_3d_config), intent(in) :: config

    real(dp) :: length_x, length_y, length_z, phase, phase_speed

    length_x = config%x_upper - config%x_lower
    length_y = config%y_upper - config%y_lower
    length_z = config%z_upper - config%z_lower
    phase_speed = &
      real(config%wave_number_x, dp) * config%initial_velocity_x / length_x + &
      real(config%wave_number_y, dp) * config%initial_velocity_y / length_y + &
      real(config%wave_number_z, dp) * config%initial_velocity_z / length_z
    phase = 2.0_dp * pi * ( &
      real(config%wave_number_x, dp) * &
        (x - config%x_lower) / length_x + &
      real(config%wave_number_y, dp) * &
        (y - config%y_lower) / length_y + &
      real(config%wave_number_z, dp) * &
        (z - config%z_lower) / length_z - time * phase_speed)
    density = base_density * &
      (1.0_dp + config%density_wave_amplitude * sin(phase))
  end function reactive_entropy_wave_density_3d

end module reactive_entropy_wave_3d_problem_mod
