module simulation_config_multispecies_mod
  use precision_mod, only: dp
  use constants_mod, only: default_gamma
  use multispecies_state_mod, only: &
    max_supported_species, multispecies_nvar, species_closure_tolerance
  implicit none
  private

  type, public :: multispecies_simulation_config
    integer :: nx = 400
    integer :: max_steps = 100000
    integer :: nspecies = 2
    real(dp) :: x_min = 0.0_dp
    real(dp) :: x_max = 1.0_dp
    real(dp) :: final_time = 0.2_dp
    real(dp) :: cfl = 0.35_dp
    real(dp) :: gamma = default_gamma
    character(len=512) :: output_file = "multispec_sod.csv"
    character(len=32) :: problem = "multispec_sod"
    character(len=32) :: reconstruction = "pelec_plm"
    character(len=32) :: limiter = "mc"
    character(len=32) :: boundary_condition = "outflow"
    character(len=32) :: riemann_solver = "pelec"
    integer :: plm_order = 2
    logical :: use_flattening = .false.
  end type multispecies_simulation_config

  type, public :: multispec_sod_config
    real(dp) :: discontinuity = 0.5_dp
    real(dp) :: rho_left = 1.0_dp
    real(dp) :: velocity_left = 0.0_dp
    real(dp) :: pressure_left = 1.0_dp
    real(dp) :: rho_right = 0.125_dp
    real(dp) :: velocity_right = 0.0_dp
    real(dp) :: pressure_right = 0.1_dp
    real(dp) :: mass_fractions_left(max_supported_species) = 0.0_dp
    real(dp) :: mass_fractions_right(max_supported_species) = 0.0_dp
  end type multispec_sod_config

  public :: read_multispecies_configuration
  public :: validate_multispecies_configuration

contains

  subroutine read_multispecies_configuration(path, config, problem, ok, message)
    character(len=*), intent(in) :: path
    type(multispecies_simulation_config), intent(out) :: config
    type(multispec_sod_config), intent(out) :: problem
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: nx, max_steps, nspecies, plm_order, unit, io_status
    real(dp) :: x_min, x_max, final_time, cfl, gamma
    real(dp) :: discontinuity, rho_left, velocity_left, pressure_left
    real(dp) :: rho_right, velocity_right, pressure_right
    real(dp) :: mass_fractions_left(max_supported_species)
    real(dp) :: mass_fractions_right(max_supported_species)
    logical :: use_flattening
    character(len=512) :: output_file
    character(len=32) :: problem_name, reconstruction, limiter
    character(len=32) :: boundary_condition, riemann_solver

    namelist /simulation/ nx, max_steps, nspecies, x_min, x_max, final_time, &
      cfl, gamma, output_file, problem_name, reconstruction, limiter, &
      boundary_condition, riemann_solver, plm_order, use_flattening
    namelist /multispec_sod_problem/ discontinuity, rho_left, velocity_left, &
      pressure_left, rho_right, velocity_right, pressure_right, &
      mass_fractions_left, mass_fractions_right

    config = multispecies_simulation_config()
    problem = multispec_sod_config()
    problem%mass_fractions_left = 0.0_dp
    problem%mass_fractions_right = 0.0_dp
    problem%mass_fractions_left(1) = 1.0_dp
    problem%mass_fractions_right(2) = 1.0_dp

    nx = config%nx
    max_steps = config%max_steps
    nspecies = config%nspecies
    x_min = config%x_min
    x_max = config%x_max
    final_time = config%final_time
    cfl = config%cfl
    gamma = config%gamma
    output_file = config%output_file
    problem_name = config%problem
    reconstruction = config%reconstruction
    limiter = config%limiter
    boundary_condition = config%boundary_condition
    riemann_solver = config%riemann_solver
    plm_order = config%plm_order
    use_flattening = config%use_flattening

    discontinuity = problem%discontinuity
    rho_left = problem%rho_left
    velocity_left = problem%velocity_left
    pressure_left = problem%pressure_left
    rho_right = problem%rho_right
    velocity_right = problem%velocity_right
    pressure_right = problem%pressure_right
    mass_fractions_left = problem%mass_fractions_left
    mass_fractions_right = problem%mass_fractions_right

    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not open input file:", trim(path)
      return
    end if

    read(unit, nml=simulation, iostat=io_status)
    if (io_status /= 0) then
      close(unit)
      ok = .false.
      write(message, '(a,1x,a)') &
        "Could not read &simulation from:", trim(path)
      return
    end if
    if (trim(problem_name) /= "multispec_sod") then
      close(unit)
      ok = .false.
      write(message, '(a,1x,a)') "Unsupported problem:", trim(problem_name)
      return
    end if

    read(unit, nml=multispec_sod_problem, iostat=io_status)
    close(unit)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') &
        "Could not read &multispec_sod_problem from:", trim(path)
      return
    end if

    config%nx = nx
    config%max_steps = max_steps
    config%nspecies = nspecies
    config%x_min = x_min
    config%x_max = x_max
    config%final_time = final_time
    config%cfl = cfl
    config%gamma = gamma
    config%output_file = trim(output_file)
    config%problem = trim(problem_name)
    config%reconstruction = trim(reconstruction)
    config%limiter = trim(limiter)
    config%boundary_condition = trim(boundary_condition)
    config%riemann_solver = trim(riemann_solver)
    config%plm_order = plm_order
    config%use_flattening = use_flattening

    problem%discontinuity = discontinuity
    problem%rho_left = rho_left
    problem%velocity_left = velocity_left
    problem%pressure_left = pressure_left
    problem%rho_right = rho_right
    problem%velocity_right = velocity_right
    problem%pressure_right = pressure_right
    problem%mass_fractions_left = mass_fractions_left
    problem%mass_fractions_right = mass_fractions_right

    call validate_multispecies_configuration(config, problem, ok, message)
  end subroutine read_multispecies_configuration

  pure subroutine validate_multispecies_configuration( &
      config, problem, ok, message)
    type(multispecies_simulation_config), intent(in) :: config
    type(multispec_sod_config), intent(in) :: problem
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp) :: left_sum, right_sum

    ok = .false.
    message = ""
    if (config%nx < 10) then
      message = "nx must be at least 10"
    else if (multispecies_nvar(config%nspecies) == 0) then
      message = "nspecies is outside the supported range"
    else if (config%x_max <= config%x_min) then
      message = "x_max must be greater than x_min"
    else if (config%final_time <= 0.0_dp) then
      message = "final_time must be positive"
    else if (config%cfl <= 0.0_dp .or. config%cfl > 1.0_dp) then
      message = "cfl must be in (0, 1]"
    else if (config%gamma <= 1.0_dp) then
      message = "gamma must be greater than 1"
    else if (config%max_steps <= 0) then
      message = "max_steps must be positive"
    else if (trim(config%problem) /= "multispec_sod") then
      message = "problem must be multispec_sod"
    else if (.not. valid_reconstruction(config%reconstruction)) then
      message = "reconstruction must be pcm, plm, or pelec_plm"
    else if (.not. valid_limiter(config%limiter)) then
      message = "limiter must be minmod or mc"
    else if (.not. valid_boundary_condition(config%boundary_condition)) then
      message = "boundary_condition must be outflow or periodic"
    else if (.not. valid_riemann_solver(config%riemann_solver)) then
      message = "riemann_solver must be rusanov or pelec"
    else if (config%plm_order /= 2 .and. config%plm_order /= 4) then
      message = "plm_order must be 2 or 4"
    else if (problem%discontinuity <= config%x_min .or. &
             problem%discontinuity >= config%x_max) then
      message = "discontinuity must lie inside the domain"
    else if (problem%rho_left <= 0.0_dp .or. problem%rho_right <= 0.0_dp) then
      message = "left and right densities must be positive"
    else if (problem%pressure_left <= 0.0_dp .or. &
             problem%pressure_right <= 0.0_dp) then
      message = "left and right pressures must be positive"
    else if (any(problem%mass_fractions_left(1:config%nspecies) < 0.0_dp) .or. &
             any(problem%mass_fractions_right(1:config%nspecies) < 0.0_dp)) then
      message = "mass fractions must be non-negative"
    else if (any(problem%mass_fractions_left(1:config%nspecies) > 1.0_dp) .or. &
             any(problem%mass_fractions_right(1:config%nspecies) > 1.0_dp)) then
      message = "mass fractions must not exceed one"
    else
      left_sum = sum(problem%mass_fractions_left(1:config%nspecies))
      right_sum = sum(problem%mass_fractions_right(1:config%nspecies))
      if (abs(left_sum - 1.0_dp) > species_closure_tolerance .or. &
          abs(right_sum - 1.0_dp) > species_closure_tolerance) then
        message = "left and right mass fractions must each sum to one"
      else
        ok = .true.
      end if
    end if
  end subroutine validate_multispecies_configuration

  pure logical function valid_reconstruction(name) result(valid)
    character(len=*), intent(in) :: name
    select case (trim(name))
    case ("pcm", "plm", "pelec_plm")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_reconstruction

  pure logical function valid_limiter(name) result(valid)
    character(len=*), intent(in) :: name
    select case (trim(name))
    case ("minmod", "mc")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_limiter

  pure logical function valid_boundary_condition(name) result(valid)
    character(len=*), intent(in) :: name
    select case (trim(name))
    case ("outflow", "periodic")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_boundary_condition

  pure logical function valid_riemann_solver(name) result(valid)
    character(len=*), intent(in) :: name
    select case (trim(name))
    case ("rusanov", "pelec")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_riemann_solver

end module simulation_config_multispecies_mod
