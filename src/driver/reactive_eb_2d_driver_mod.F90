module reactive_eb_2d_driver_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_conserved_to_primitive
  use reactive_2d_mod, only: initialize_reactive_2d
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, build_eb_geometry_2d
  use eb_reactive_hydro_2d_mod, only: advance_reactive_eb_hydro_2d
  use simulation_config_reactive_eb_2d_mod, only: reactive_eb_2d_config
  implicit none
  private

  public :: build_configured_eb_geometry_2d
  public :: compute_reactive_eb_cfl_timestep_2d
  public :: reactive_eb_integrals_2d
  public :: reactive_eb_extrema_2d
  public :: simulate_reactive_eb_2d
  public :: write_reactive_eb_2d_csv

contains

  subroutine build_configured_eb_geometry_2d(config, geometry, ok)
    type(reactive_eb_2d_config), intent(in) :: config
    type(eb_geometry_2d), intent(out) :: geometry
    logical, intent(out) :: ok

    real(dp), allocatable :: level_set(:, :)
    real(dp) :: x, y, distance
    integer :: i, j

    ok = .false.
    allocate(level_set(0:config%flow%nx, 0:config%flow%ny))
    do j = 0, config%flow%ny
      y = config%flow%y_lower + real(j, dp) * &
        (config%flow%y_upper - config%flow%y_lower) / &
        real(config%flow%ny, dp)
      do i = 0, config%flow%nx
        x = config%flow%x_lower + real(i, dp) * &
          (config%flow%x_upper - config%flow%x_lower) / &
          real(config%flow%nx, dp)
        select case (trim(config%geometry))
        case ("plane")
          level_set(i, j) = config%plane_normal_x * x + &
            config%plane_normal_y * y - config%plane_offset
        case ("circle")
          distance = sqrt((x - config%circle_center_x)**2 + &
            (y - config%circle_center_y)**2)
          if (config%circle_fluid_inside) then
            level_set(i, j) = config%circle_radius - distance
          else
            level_set(i, j) = distance - config%circle_radius
          end if
        case default
          return
        end select
      end do
    end do
    call build_eb_geometry_2d( &
      level_set, config%flow%x_lower, config%flow%x_upper, &
      config%flow%y_lower, config%flow%y_upper, geometry, ok)
    if (.not. ok) return
    ok = count(geometry%cell_type /= eb_covered_cell) > 0 .and. &
      count(geometry%cell_type == eb_cut_cell) > 0
  end subroutine build_configured_eb_geometry_2d

  subroutine compute_reactive_eb_cfl_timestep_2d( &
      species, state, temperature, geometry, cfl, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed, rate, maximum_rate
    logical :: local_ok
    integer :: i, j, active_cells

    dt = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid() .or. cfl <= 0.0_dp .or. cfl > 1.0_dp) return
    if (size(state, 1) /= reactive_nvar(size(species)) .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny])) return
    allocate(primitive(reactive_nprim(size(species))))
    maximum_rate = 0.0_dp
    active_cells = 0
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        active_cells = active_cells + 1
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          local_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        rate = (abs(primitive(2)) + sound_speed) / geometry%dx + &
          (abs(primitive(3)) + sound_speed) / geometry%dy
        maximum_rate = max(maximum_rate, rate)
      end do
    end do
    if (active_cells == 0 .or. maximum_rate <= 0.0_dp) return
    dt = cfl / maximum_rate
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp
  end subroutine compute_reactive_eb_cfl_timestep_2d

  subroutine reactive_eb_integrals_2d(state, geometry, integrals, ok)
    real(dp), intent(in) :: state(:, :, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(out) :: integrals(:)
    logical, intent(out) :: ok

    integer :: i, j

    integrals = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid()) return
    if (size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(integrals) /= size(state, 1)) return
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        if (any(.not. ieee_is_finite(state(:, i, j)))) return
        integrals = integrals + &
          geometry%volume_fraction(i, j) * state(:, i, j)
      end do
    end do
    integrals = integrals * geometry%dx * geometry%dy
    ok = all(ieee_is_finite(integrals))
  end subroutine reactive_eb_integrals_2d

  subroutine reactive_eb_extrema_2d( &
      species, state, temperature, geometry, minimum_density, &
      maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, &
      maximum_closure_error, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(out) :: minimum_density, maximum_density
    real(dp), intent(out) :: minimum_pressure, maximum_pressure
    real(dp), intent(out) :: minimum_temperature, maximum_temperature
    real(dp), intent(out) :: maximum_speed, maximum_closure_error
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed, closure, speed
    logical :: local_ok
    integer :: i, j, k, active_cells

    minimum_density = huge(1.0_dp)
    maximum_density = -huge(1.0_dp)
    minimum_pressure = huge(1.0_dp)
    maximum_pressure = -huge(1.0_dp)
    minimum_temperature = huge(1.0_dp)
    maximum_temperature = -huge(1.0_dp)
    maximum_speed = 0.0_dp
    maximum_closure_error = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid()) return
    if (size(state, 1) /= reactive_nvar(size(species)) .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny])) return
    allocate(primitive(reactive_nprim(size(species))))
    active_cells = 0
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        active_cells = active_cells + 1
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          local_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        minimum_density = min(minimum_density, primitive(1))
        maximum_density = max(maximum_density, primitive(1))
        minimum_pressure = min(minimum_pressure, primitive(5))
        maximum_pressure = max(maximum_pressure, primitive(5))
        minimum_temperature = min(minimum_temperature, local_temperature)
        maximum_temperature = max(maximum_temperature, local_temperature)
        speed = sqrt(primitive(2)**2 + primitive(3)**2 + primitive(4)**2)
        maximum_speed = max(maximum_speed, speed)
        closure = 0.0_dp
        do k = 1, size(species)
          closure = closure + &
            primitive(reactive_mass_fraction_component(k))
        end do
        maximum_closure_error = max( &
          maximum_closure_error, abs(closure - 1.0_dp))
      end do
    end do
    ok = active_cells > 0
  end subroutine reactive_eb_extrema_2d

  subroutine simulate_reactive_eb_2d( &
      species, config, state, temperature, geometry, time, steps, &
      initial_integrals, final_integrals, minimum_dt, base_density, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_2d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(out) :: geometry
    real(dp), intent(out) :: time, minimum_dt, base_density
    integer, intent(out) :: steps
    real(dp), allocatable, intent(out) :: initial_integrals(:)
    real(dp), allocatable, intent(out) :: final_integrals(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :)
    real(dp) :: dx, dy, dt, remaining, time_tolerance
    logical :: local_ok
    integer :: nvar

    ok = .false.
    time = 0.0_dp
    steps = 0
    minimum_dt = 0.0_dp
    base_density = 0.0_dp
    call build_configured_eb_geometry_2d(config, geometry, local_ok)
    if (.not. local_ok) return
    call initialize_reactive_2d( &
      species, config%flow, state, temperature, dx, dy, base_density, local_ok)
    if (.not. local_ok) return
    if (abs(dx - geometry%dx) > 8.0_dp * epsilon(1.0_dp) * geometry%dx .or. &
        abs(dy - geometry%dy) > 8.0_dp * epsilon(1.0_dp) * geometry%dy) return

    nvar = size(state, 1)
    allocate(initial_integrals(nvar), final_integrals(nvar))
    call reactive_eb_integrals_2d( &
      state, geometry, initial_integrals, local_ok)
    if (.not. local_ok) return
    allocate(candidate_state, mold=state)
    allocate(candidate_temperature, mold=temperature)
    minimum_dt = huge(1.0_dp)
    time_tolerance = 16.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(config%flow%final_time))

    do
      remaining = config%flow%final_time - time
      if (remaining <= time_tolerance) exit
      if (steps >= config%flow%maximum_steps) return
      call compute_reactive_eb_cfl_timestep_2d( &
        species, state, temperature, geometry, config%flow%cfl, dt, local_ok)
      if (.not. local_ok) return
      dt = min(dt, remaining)
      call advance_reactive_eb_hydro_2d( &
        species, state, temperature, geometry, config%flow%riemann_solver, &
        dt, candidate_state, candidate_temperature, local_ok, &
        config%state_redist_target_volume_fraction)
      if (.not. local_ok) return
      state = candidate_state
      temperature = candidate_temperature
      time = time + dt
      minimum_dt = min(minimum_dt, dt)
      steps = steps + 1
    end do
    time = config%flow%final_time
    call reactive_eb_integrals_2d( &
      state, geometry, final_integrals, local_ok)
    if (.not. local_ok) return
    ok = steps > 0 .and. ieee_is_finite(minimum_dt) .and. minimum_dt > 0.0_dp
  end subroutine simulate_reactive_eb_2d

  subroutine write_reactive_eb_2d_csv( &
      path, species, config, state, temperature, geometry, time, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_2d_config), intent(in) :: config
    real(dp), intent(in) :: state(:, :, :), temperature(:, :), time
    type(eb_geometry_2d), intent(in) :: geometry
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: x, y, local_temperature, sound_speed
    logical :: local_ok
    integer :: unit, status, i, j, k

    ok = .false.
    if (.not. geometry%is_valid()) return
    if (size(state, 1) /= reactive_nvar(size(species)) .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny])) return
    allocate(primitive(reactive_nprim(size(species))))
    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=status)
    if (status /= 0) return
    write(unit, '(a)', advance='no') &
      "time,x,y,volume_fraction,cell_type,boundary_length," // &
      "boundary_normal_x,boundary_normal_y,rho,u,v,w,pressure," // &
      "temperature,rhoE"
    do k = 1, size(species)
      write(unit, '(a)', advance='no') ",Y_" // trim(species(k)%name)
    end do
    write(unit, '(a)') ""
    do j = 1, geometry%ny
      y = config%flow%y_lower + (real(j, dp) - 0.5_dp) * geometry%dy
      do i = 1, geometry%nx
        x = config%flow%x_lower + (real(i, dp) - 0.5_dp) * geometry%dx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          local_temperature, sound_speed, local_ok)
        if (.not. local_ok) then
          close(unit)
          return
        end if
        write(unit, '(*(g0,:,","))') &
          time, x, y, geometry%volume_fraction(i, j), &
          geometry%cell_type(i, j), geometry%boundary_length(i, j), &
          geometry%boundary_normal_x(i, j), &
          geometry%boundary_normal_y(i, j), state(irho, i, j), &
          primitive(2), primitive(3), primitive(4), primitive(5), &
          local_temperature, state(iet, i, j), &
          (primitive(reactive_mass_fraction_component(k)), &
            k = 1, size(species))
      end do
    end do
    close(unit)
    ok = .true.
  end subroutine write_reactive_eb_2d_csv

end module reactive_eb_2d_driver_mod
