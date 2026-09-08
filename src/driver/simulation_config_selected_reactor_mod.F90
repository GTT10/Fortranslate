module simulation_config_selected_reactor_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  implicit none
  private

  integer, parameter :: species_name_length = 24
  integer, parameter, public :: selected_reactor_max_species = 32

  type, public :: selected_reactor_config
    character(len=16) :: integrator = "implicit"
    real(dp) :: final_time = 1.0e-6_dp
    real(dp) :: output_interval = 1.0e-7_dp
    real(dp) :: initial_time_step = 1.0e-10_dp
    real(dp) :: minimum_time_step = 1.0e-16_dp
    real(dp) :: maximum_time_step = 1.0e-7_dp
    real(dp) :: relative_tolerance = 1.0e-6_dp
    real(dp) :: absolute_tolerance = 1.0e-12_dp
    integer :: maximum_steps = 100000
    integer :: cvode_max_internal_steps = 100000
    real(dp) :: initial_temperature = 1000.0_dp
    real(dp) :: initial_pressure = 101325.0_dp
    integer :: composition_count = 0
    character(len=species_name_length) :: &
      composition_species(selected_reactor_max_species) = ""
    real(dp) :: &
      composition_mole_fractions(selected_reactor_max_species) = 0.0_dp
    character(len=256) :: output_file = "selected_reactor.csv"
  end type selected_reactor_config

  public :: read_selected_reactor_configuration
  public :: validate_selected_reactor_configuration
  public :: resolve_selected_reactor_composition
  public :: selected_reactor_paths_alias_lexically

contains

  subroutine read_selected_reactor_configuration(path, config, ok, message)
    character(len=*), intent(in) :: path
    type(selected_reactor_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp) :: final_time, output_interval
    real(dp) :: initial_time_step, minimum_time_step, maximum_time_step
    real(dp) :: relative_tolerance, absolute_tolerance
    real(dp) :: initial_temperature, initial_pressure
    real(dp) :: composition_mole_fractions(selected_reactor_max_species)
    integer :: maximum_steps, cvode_max_internal_steps, composition_count
    character(len=species_name_length) :: &
      composition_species(selected_reactor_max_species)
    character(len=16) :: integrator
    character(len=256) :: output_file
    integer :: unit, io_status
    namelist /selected_reactor/ &
      final_time, output_interval, initial_time_step, minimum_time_step, &
      maximum_time_step, relative_tolerance, absolute_tolerance, &
      maximum_steps, cvode_max_internal_steps, &
      initial_temperature, initial_pressure, &
      composition_count, composition_species, composition_mole_fractions, &
      integrator, output_file

    config = selected_reactor_config()
    final_time = config%final_time
    output_interval = config%output_interval
    initial_time_step = config%initial_time_step
    minimum_time_step = config%minimum_time_step
    maximum_time_step = config%maximum_time_step
    relative_tolerance = config%relative_tolerance
    absolute_tolerance = config%absolute_tolerance
    maximum_steps = config%maximum_steps
    cvode_max_internal_steps = config%cvode_max_internal_steps
    initial_temperature = config%initial_temperature
    initial_pressure = config%initial_pressure
    composition_count = config%composition_count
    composition_species = config%composition_species
    composition_mole_fractions = config%composition_mole_fractions
    integrator = config%integrator
    output_file = config%output_file

    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      message = "Could not open selected reactor input"
      return
    end if
    read(unit, nml=selected_reactor, iostat=io_status)
    close(unit)
    if (io_status /= 0) then
      ok = .false.
      message = "Could not parse &selected_reactor namelist"
      return
    end if

    config%final_time = final_time
    config%output_interval = output_interval
    config%initial_time_step = initial_time_step
    config%minimum_time_step = minimum_time_step
    config%maximum_time_step = maximum_time_step
    config%relative_tolerance = relative_tolerance
    config%absolute_tolerance = absolute_tolerance
    config%maximum_steps = maximum_steps
    config%cvode_max_internal_steps = cvode_max_internal_steps
    config%initial_temperature = initial_temperature
    config%initial_pressure = initial_pressure
    config%composition_count = composition_count
    config%composition_species = composition_species
    config%composition_mole_fractions = composition_mole_fractions
    config%integrator = integrator
    config%output_file = output_file
    call validate_selected_reactor_configuration(config, ok, message)
  end subroutine read_selected_reactor_configuration

  subroutine validate_selected_reactor_configuration(config, ok, message)
    type(selected_reactor_config), intent(in) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: first_index, second_index

    ok = .false.
    if (trim(config%integrator) /= "implicit" .and. &
        trim(config%integrator) /= "cvode") then
      message = "Selected reactor integrator must be 'implicit' or 'cvode'"
      return
    end if
    if (any(.not. ieee_is_finite([ &
        config%final_time, config%output_interval, &
        config%initial_time_step, config%minimum_time_step, &
        config%maximum_time_step, config%relative_tolerance, &
        config%absolute_tolerance, config%initial_temperature, &
        config%initial_pressure]))) then
      message = "Selected reactor scalar configuration must be finite"
      return
    end if
    if (config%final_time <= 0.0_dp) then
      message = "Selected reactor final_time must be positive"
      return
    end if
    if (config%output_interval <= 0.0_dp) then
      message = "Selected reactor output_interval must be positive"
      return
    end if
    if (config%initial_time_step <= 0.0_dp .or. &
        config%minimum_time_step <= 0.0_dp) then
      message = "Selected reactor time steps must be positive"
      return
    end if
    if (config%initial_time_step < config%minimum_time_step .or. &
        config%initial_time_step > config%maximum_time_step) then
      message = "Selected reactor initial_time_step must be within the time-step bounds"
      return
    end if
    if (config%maximum_time_step < config%minimum_time_step) then
      message = "Selected reactor maximum_time_step is below the minimum"
      return
    end if
    if (config%output_interval < config%minimum_time_step) then
      message = "Selected reactor output_interval is below the minimum time step"
      return
    end if
    if (config%final_time < config%minimum_time_step) then
      message = "Selected reactor final_time is below the minimum time step"
      return
    end if
    if (config%relative_tolerance <= 0.0_dp .or. &
        config%absolute_tolerance <= 0.0_dp) then
      message = "Selected reactor tolerances must be positive"
      return
    end if
    if (config%maximum_steps <= 0) then
      message = "Selected reactor maximum_steps must be positive"
      return
    end if
    if (config%cvode_max_internal_steps <= 0) then
      message = "Selected reactor cvode_max_internal_steps must be positive"
      return
    end if
    if (config%initial_temperature <= 0.0_dp) then
      message = "Selected reactor initial_temperature must be finite and positive"
      return
    end if
    if (config%initial_pressure <= 0.0_dp) then
      message = "Selected reactor initial_pressure must be finite and positive"
      return
    end if
    if (config%composition_count < 1 .or. &
        config%composition_count > selected_reactor_max_species) then
      message = "Selected reactor composition_count is out of range"
      return
    end if
    if (any(.not. ieee_is_finite(config%composition_mole_fractions))) then
      message = "Selected reactor composition contains a nonfinite mole fraction"
      return
    end if
    if (any(len_trim( &
        config%composition_species(1:config%composition_count)) == 0)) then
      message = "Selected reactor composition contains a blank species name"
      return
    end if
    if (any(config%composition_mole_fractions( &
        1:config%composition_count) < 0.0_dp) .or. &
        maxval(config%composition_mole_fractions( &
          1:config%composition_count)) <= 0.0_dp) then
      message = "Selected reactor composition must be nonnegative and nonempty"
      return
    end if
    if (config%composition_count < selected_reactor_max_species) then
      if (any(len_trim(config%composition_species( &
          config%composition_count + 1:)) /= 0)) then
        message = "Selected reactor composition has names beyond composition_count"
        return
      end if
      if (any(abs(config%composition_mole_fractions( &
          config%composition_count + 1:)) > 0.0_dp)) then
        message = "Selected reactor composition has values beyond composition_count"
        return
      end if
    end if
    do first_index = 1, config%composition_count - 1
      do second_index = first_index + 1, config%composition_count
        if (trim(config%composition_species(first_index)) == &
            trim(config%composition_species(second_index))) then
          message = "Selected reactor composition contains a duplicate species"
          return
        end if
      end do
    end do
    if (len_trim(config%output_file) == 0) then
      message = "Selected reactor output_file must not be blank"
      return
    end if

    message = ""
    ok = .true.
  end subroutine validate_selected_reactor_configuration

  subroutine resolve_selected_reactor_composition( &
      config, species, mole_fractions, ok, message)
    type(selected_reactor_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(out) :: mole_fractions(:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    logical, allocatable :: assigned(:)
    integer :: input_index, species_index, match_index, match_count
    real(dp) :: scale, total

    mole_fractions = 0.0_dp
    call validate_selected_reactor_configuration(config, ok, message)
    if (.not. ok) return
    ok = .false.
    if (size(species) < 2 .or. &
        size(species) > selected_reactor_max_species) then
      message = "Selected reactor mechanism species count is unsupported"
      return
    end if
    if (size(mole_fractions) /= size(species)) then
      message = "Selected reactor composition output has the wrong size"
      return
    end if
    if (config%composition_count < 1 .or. &
        config%composition_count > size(species)) then
      message = "Selected reactor composition_count exceeds the mechanism"
      return
    end if
    allocate(assigned(size(species)))
    assigned = .false.
    do input_index = 1, config%composition_count
      if (len_trim(config%composition_species(input_index)) == 0) then
        message = "Selected reactor composition contains a blank species name"
        return
      end if
      match_index = 0
      match_count = 0
      do species_index = 1, size(species)
        if (trim(config%composition_species(input_index)) == &
            trim(species(species_index)%name)) then
          match_index = species_index
          match_count = match_count + 1
        end if
      end do
      if (match_count /= 1) then
        message = "Selected reactor species was not found exactly once in the bundle: " // &
          trim(config%composition_species(input_index))
        return
      end if
      if (assigned(match_index)) then
        message = "Selected reactor composition contains a duplicate species: " // &
          trim(config%composition_species(input_index))
        return
      end if
      assigned(match_index) = .true.
      mole_fractions(match_index) = &
        config%composition_mole_fractions(input_index)
    end do

    scale = maxval(mole_fractions)
    if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) then
      message = "Selected reactor composition has no positive total"
      return
    end if
    if (scale <= huge(1.0_dp) / real(size(mole_fractions), dp)) then
      total = sum(mole_fractions)
    else
      mole_fractions = mole_fractions / scale
      total = sum(mole_fractions)
    end if
    if (.not. ieee_is_finite(total) .or. total <= 0.0_dp) then
      message = "Selected reactor composition normalization failed"
      return
    end if
    mole_fractions = mole_fractions / total
    message = ""
    ok = .true.
  end subroutine resolve_selected_reactor_composition

  pure logical function selected_reactor_paths_alias_lexically( &
      left, right) result(aliases)
    character(len=*), intent(in) :: left, right

    aliases = normalized_path(left) == normalized_path(right)
  end function selected_reactor_paths_alias_lexically

  pure function normalized_path(path) result(normalized)
    character(len=*), intent(in) :: path
    character(len=1024) :: normalized
    character(len=1024) :: input, component
    integer :: component_starts(512), component_ends(512)
    integer :: input_length, start, finish, component_count, index
    logical :: absolute

    normalized = ""
    input = trim(path)
    input_length = len_trim(input)
    if (input_length == 0) return
    absolute = input(1:1) == "/"
    component_count = 0
    start = 1
    do while (start <= input_length)
      do while (start <= input_length)
        if (input(start:start) /= "/") exit
        start = start + 1
      end do
      if (start > input_length) exit
      finish = start
      do while (finish <= input_length)
        if (input(finish:finish) == "/") exit
        finish = finish + 1
      end do
      component = input(start:finish - 1)
      select case (trim(component))
      case ("", ".")
        continue
      case ("..")
        if (component_count > 0) then
          if (trim(input(component_starts(component_count): &
              component_ends(component_count))) /= "..") then
            component_count = component_count - 1
          else if (.not. absolute) then
            component_count = component_count + 1
            component_starts(component_count) = start
            component_ends(component_count) = finish - 1
          end if
        else if (.not. absolute) then
          component_count = component_count + 1
          component_starts(component_count) = start
          component_ends(component_count) = finish - 1
        end if
      case default
        component_count = component_count + 1
        component_starts(component_count) = start
        component_ends(component_count) = finish - 1
      end select
      start = finish + 1
    end do

    if (absolute) normalized = "/"
    do index = 1, component_count
      if (len_trim(normalized) > 0 .and. &
          trim(normalized) /= "/") normalized = trim(normalized) // "/"
      normalized = trim(normalized) // &
        trim(input(component_starts(index):component_ends(index)))
    end do
    if (len_trim(normalized) == 0) then
      if (absolute) then
        normalized = "/"
      else
        normalized = "."
      end if
    end if
  end function normalized_path

end module simulation_config_selected_reactor_mod
