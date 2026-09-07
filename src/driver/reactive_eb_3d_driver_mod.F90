module reactive_eb_3d_driver_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_cut_cell_3d
  use eb_reactive_hydro_3d_mod, only: &
    advance_reactive_eb_redistributed_euler_3d, &
    advance_reactive_eb_state_redistributed_euler_3d
  use reactive_3d_mod, only: advance_reactive_chemistry_3d
  use eb_reactive_transport_3d_mod, only: advance_reactive_eb_transport_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  implicit none
  private

  public :: initialize_reactive_eb_density_sheet_3d
  public :: advance_reactive_eb_strang_3d
  public :: advance_reactive_eb_full_3d
  public :: reactive_eb_integrals_3d
  public :: reactive_eb_element_integrals_3d
  public :: reactive_eb_element_species_supported_3d
  public :: reactive_eb_extrema_3d
  public :: write_reactive_eb_3d_csv

contains

  subroutine initialize_reactive_eb_density_sheet_3d( &
      species, config, geometry, state, temperature, ok, &
      base_mole_fractions)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_3d_config), intent(in) :: config
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: state(:, :, :, :), temperature(:, :, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: base_mole_fractions(:)

    real(dp), parameter :: default_mole_fractions(7) = [ &
      0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    real(dp), allocatable :: primitive(:), mass_fractions(:)
    real(dp), allocatable :: mole_fractions(:)
    real(dp), allocatable :: regular_state(:), cut_state(:)
    real(dp) :: regular_temperature, cut_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, k, nvar, species_index

    state = 0.0_dp
    temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz])) return

    allocate(primitive(reactive_nprim(size(species))))
    allocate(mass_fractions(size(species)))
    allocate(mole_fractions(size(species)))
    allocate(regular_state(nvar), cut_state(nvar))
    if (present(base_mole_fractions)) then
      if (size(base_mole_fractions) /= size(species)) return
      mole_fractions = base_mole_fractions
    else
      if (size(species) /= size(default_mole_fractions)) return
      mole_fractions = default_mole_fractions
    end if
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, local_ok)
    if (.not. local_ok) return
    primitive(1:5) = [ &
      config%regular_density, config%initial_velocity_x, &
      config%initial_velocity_y, config%initial_velocity_z, &
      config%initial_pressure]
    do species_index = 1, size(species)
      primitive(reactive_mass_fraction_component(species_index)) = &
        mass_fractions(species_index)
    end do
    call reactive_primitive_to_conserved( &
      species, primitive, regular_state, regular_temperature, sound_speed, &
      local_ok)
    if (.not. local_ok) return
    primitive(1) = config%cut_density
    call reactive_primitive_to_conserved( &
      species, primitive, cut_state, cut_temperature, sound_speed, local_ok)
    if (.not. local_ok) return

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_cut_cell_3d) then
            state(:, i, j, k) = cut_state
            temperature(i, j, k) = cut_temperature
          else
            state(:, i, j, k) = regular_state
            temperature(i, j, k) = regular_temperature
          end if
        end do
      end do
    end do
    ok = all(ieee_is_finite(state)) .and. &
      all(ieee_is_finite(temperature))
  end subroutine initialize_reactive_eb_density_sheet_3d

  subroutine advance_reactive_eb_strang_3d( &
      species, reactions, state, temperature, geometry, solver, redistribution, &
      dt, chemistry_enabled, rtol, atol, new_state, new_temperature, ok, &
      target_volume_fraction, chemistry_integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: solver, redistribution
    real(dp), intent(in) :: dt, rtol, atol
    logical, intent(in) :: chemistry_enabled
    real(dp), intent(out) :: new_state(:, :, :, :), new_temperature(:, :, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: hydro_state(:, :, :, :)
    real(dp), allocatable :: hydro_temperature(:, :, :)
    logical, allocatable :: active_mask(:, :, :)
    real(dp) :: target
    logical :: local_ok
    integer :: nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny, geometry%nz]) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. ieee_is_finite(rtol) .or. .not. ieee_is_finite(atol) .or. &
        rtol <= 0.0_dp .or. atol <= 0.0_dp) return
    if (chemistry_enabled .and. size(reactions) < 1) return
    target = 0.5_dp
    if (present(target_volume_fraction)) target = target_volume_fraction
    if (.not. ieee_is_finite(target) .or. target <= 0.0_dp .or. &
        target > 1.0_dp) return
    select case (trim(redistribution))
    case ("flux_redist", "state_redist")
    case default
      return
    end select

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    allocate(hydro_state, source=state)
    allocate(hydro_temperature, source=temperature)
    allocate(active_mask(geometry%nx, geometry%ny, geometry%nz))
    active_mask = geometry%cell_type /= eb_covered_cell_3d

    if (chemistry_enabled) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        geometry%nx, geometry%ny, geometry%nz, 0.5_dp * dt, rtol, atol, &
        local_ok, active_mask, chemistry_integrator)
      if (.not. local_ok) return
      call covered_cells_match_3d( &
        state, temperature, candidate_state, candidate_temperature, geometry, &
        local_ok)
      if (.not. local_ok) return
    end if

    select case (trim(redistribution))
    case ("flux_redist")
      call advance_reactive_eb_redistributed_euler_3d( &
        species, candidate_state, candidate_temperature, geometry, solver, dt, &
        hydro_state, hydro_temperature, local_ok)
    case ("state_redist")
      call advance_reactive_eb_state_redistributed_euler_3d( &
        species, candidate_state, candidate_temperature, geometry, solver, dt, &
        hydro_state, hydro_temperature, local_ok, target)
    end select
    if (.not. local_ok) return
    call covered_cells_match_3d( &
      state, temperature, hydro_state, hydro_temperature, geometry, local_ok)
    if (.not. local_ok) return
    candidate_state = hydro_state
    candidate_temperature = hydro_temperature

    if (chemistry_enabled) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        geometry%nx, geometry%ny, geometry%nz, 0.5_dp * dt, rtol, atol, &
        local_ok, active_mask, chemistry_integrator)
      if (.not. local_ok) return
      call covered_cells_match_3d( &
        state, temperature, candidate_state, candidate_temperature, geometry, &
        local_ok)
      if (.not. local_ok) return
    end if
    if (any(.not. ieee_is_finite(candidate_state)) .or. &
        any(.not. ieee_is_finite(candidate_temperature))) return
    new_state = candidate_state
    new_temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_eb_strang_3d

  subroutine advance_reactive_eb_full_3d( &
      species, reactions, transport, state, temperature, geometry, solver, &
      redistribution, dt, chemistry_enabled, rtol, atol, &
      transport_enabled, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, new_state, &
      new_temperature, minimum_transport_theta, ok, target_volume_fraction, &
      chemistry_integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: solver, redistribution
    real(dp), intent(in) :: dt, rtol, atol
    logical, intent(in) :: chemistry_enabled, transport_enabled
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: new_state(:, :, :, :), new_temperature(:, :, :)
    real(dp), intent(out) :: minimum_transport_theta
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: hydro_state(:, :, :, :)
    real(dp), allocatable :: hydro_temperature(:, :, :)
    logical, allocatable :: active_mask(:, :, :)
    real(dp) :: target, theta_one, theta_two
    logical :: local_ok
    integer :: nvar

    minimum_transport_theta = 1.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz]) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. transport_enabled) then
      call advance_reactive_eb_strang_3d( &
        species, reactions, state, temperature, geometry, solver, &
        redistribution, dt, chemistry_enabled, rtol, atol, new_state, &
        new_temperature, ok, target_volume_fraction, chemistry_integrator)
      return
    end if

    if (size(transport) /= size(species)) return
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. ieee_is_finite(rtol) .or. .not. ieee_is_finite(atol) .or. &
        rtol <= 0.0_dp .or. atol <= 0.0_dp) return
    if (trim(redistribution) /= "state_redist") return
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) return
    if (barodiffusion_enabled .and. .not. species_diffusion_enabled) return
    if (chemistry_enabled .and. size(reactions) < 1) return
    target = 0.5_dp
    if (present(target_volume_fraction)) target = target_volume_fraction
    if (.not. ieee_is_finite(target) .or. target <= 0.0_dp .or. &
        target > 1.0_dp) return

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    allocate(hydro_state, source=state)
    allocate(hydro_temperature, source=temperature)
    allocate(active_mask(geometry%nx, geometry%ny, geometry%nz))
    active_mask = geometry%cell_type /= eb_covered_cell_3d

    if (chemistry_enabled) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        geometry%nx, geometry%ny, geometry%nz, 0.5_dp * dt, rtol, atol, &
        local_ok, active_mask, chemistry_integrator)
      if (.not. local_ok) return
      call covered_cells_match_3d( &
        state, temperature, candidate_state, candidate_temperature, &
        geometry, local_ok)
      if (.not. local_ok) return
    end if

    call advance_reactive_eb_transport_3d( &
      species, transport, candidate_state, candidate_temperature, &
      geometry, 0.5_dp * dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, target, theta_one, local_ok)
    if (.not. local_ok) return
    call covered_cells_match_3d( &
      state, temperature, candidate_state, candidate_temperature, &
      geometry, local_ok)
    if (.not. local_ok) return

    call advance_reactive_eb_state_redistributed_euler_3d( &
      species, candidate_state, candidate_temperature, geometry, solver, &
      dt, hydro_state, hydro_temperature, local_ok, target)
    if (.not. local_ok) return
    call covered_cells_match_3d( &
      state, temperature, hydro_state, hydro_temperature, geometry, &
      local_ok)
    if (.not. local_ok) return
    candidate_state = hydro_state
    candidate_temperature = hydro_temperature

    call advance_reactive_eb_transport_3d( &
      species, transport, candidate_state, candidate_temperature, &
      geometry, 0.5_dp * dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, target, theta_two, local_ok)
    if (.not. local_ok) return
    call covered_cells_match_3d( &
      state, temperature, candidate_state, candidate_temperature, &
      geometry, local_ok)
    if (.not. local_ok) return

    if (chemistry_enabled) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        geometry%nx, geometry%ny, geometry%nz, 0.5_dp * dt, rtol, atol, &
        local_ok, active_mask, chemistry_integrator)
      if (.not. local_ok) return
      call covered_cells_match_3d( &
        state, temperature, candidate_state, candidate_temperature, &
        geometry, local_ok)
      if (.not. local_ok) return
    end if
    if (any(.not. ieee_is_finite(candidate_state)) .or. &
        any(.not. ieee_is_finite(candidate_temperature))) return

    new_state = candidate_state
    new_temperature = candidate_temperature
    minimum_transport_theta = min(theta_one, theta_two)
    ok = .true.
  end subroutine advance_reactive_eb_full_3d

  subroutine covered_cells_match_3d( &
      reference_state, reference_temperature, candidate_state, &
      candidate_temperature, geometry, ok)
    real(dp), intent(in) :: reference_state(:, :, :, :)
    real(dp), intent(in) :: reference_temperature(:, :, :)
    real(dp), intent(in) :: candidate_state(:, :, :, :)
    real(dp), intent(in) :: candidate_temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    logical, intent(out) :: ok

    integer :: i, j, k

    ok = .false.
    if (any(shape(candidate_state) /= shape(reference_state)) .or. &
        any(shape(candidate_temperature) /= shape(reference_temperature))) return
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) /= eb_covered_cell_3d) cycle
          if (.not. all(same_real_bits(candidate_state(:, i, j, k), &
              reference_state(:, i, j, k))) .or. &
              .not. same_real_bits(candidate_temperature(i, j, k), &
                reference_temperature(i, j, k))) return
        end do
      end do
    end do
    ok = .true.
  end subroutine covered_cells_match_3d

  pure elemental logical function same_real_bits(left, right)
    real(dp), intent(in) :: left, right

    same_real_bits = transfer(left, 0_int64) == transfer(right, 0_int64)
  end function same_real_bits

  subroutine reactive_eb_integrals_3d( &
      field, geometry, integrals, ok, l1_integrals)
    real(dp), intent(in) :: field(:, :, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: integrals(:)
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: l1_integrals(:)

    real(dp) :: cell_volume
    integer :: component

    integrals = 0.0_dp
    if (present(l1_integrals)) l1_integrals = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid() .or. &
        size(field, 2) /= geometry%nx .or. &
        size(field, 3) /= geometry%ny .or. &
        size(field, 4) /= geometry%nz .or. &
        size(integrals) /= size(field, 1) .or. &
        any(.not. ieee_is_finite(field))) return
    if (present(l1_integrals)) then
      if (size(l1_integrals) /= size(field, 1)) return
    end if
    cell_volume = geometry%dx * geometry%dy * geometry%dz
    do component = 1, size(field, 1)
      integrals(component) = sum(geometry%volume_fraction * &
        field(component, :, :, :)) * cell_volume
      if (present(l1_integrals)) then
        l1_integrals(component) = sum(geometry%volume_fraction * &
          abs(field(component, :, :, :))) * cell_volume
      end if
    end do
    ok = all(ieee_is_finite(integrals))
    if (present(l1_integrals)) then
      ok = ok .and. all(ieee_is_finite(l1_integrals))
    end if
  end subroutine reactive_eb_integrals_3d

  subroutine reactive_eb_element_integrals_3d( &
      species, state, geometry, element_integrals, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: element_integrals(:)
    logical, intent(out) :: ok

    real(dp) :: atoms(3), cell_volume
    logical :: local_ok
    integer :: i, j, k, species_index, nvar

    element_integrals = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. size(element_integrals) /= 3 .or. &
        .not. geometry%is_valid() .or. size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz) return
    cell_volume = geometry%dx * geometry%dy * geometry%dz
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          if (any(.not. ieee_is_finite(state(:, i, j, k)))) return
          do species_index = 1, size(species)
            call element_atom_counts( &
              species(species_index)%name, atoms, local_ok)
            if (.not. local_ok .or. &
                .not. ieee_is_finite(species(species_index)%molecular_weight) .or. &
                species(species_index)%molecular_weight <= 0.0_dp) return
            element_integrals = element_integrals + &
              geometry%volume_fraction(i, j, k) * &
              state(reactive_species_component(species_index), i, j, k) * &
              atoms / species(species_index)%molecular_weight * cell_volume
          end do
        end do
      end do
    end do
    ok = all(ieee_is_finite(element_integrals))
  end subroutine reactive_eb_element_integrals_3d

  pure logical function reactive_eb_element_species_supported_3d( &
      species) result(supported)
    type(nasa7_species), intent(in) :: species(:)

    real(dp) :: atoms(3)
    logical :: ok
    integer :: species_index

    supported = size(species) >= 1
    do species_index = 1, size(species)
      call element_atom_counts(species(species_index)%name, atoms, ok)
      if (.not. ok) then
        supported = .false.
        return
      end if
    end do
  end function reactive_eb_element_species_supported_3d

  pure subroutine element_atom_counts(name, atoms, ok)
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: atoms(3)
    logical, intent(out) :: ok

    atoms = 0.0_dp
    select case (trim(name))
    case ("H2")
      atoms = [2.0_dp, 0.0_dp, 0.0_dp]
    case ("H")
      atoms = [1.0_dp, 0.0_dp, 0.0_dp]
    case ("O")
      atoms = [0.0_dp, 1.0_dp, 0.0_dp]
    case ("O2")
      atoms = [0.0_dp, 2.0_dp, 0.0_dp]
    case ("OH")
      atoms = [1.0_dp, 1.0_dp, 0.0_dp]
    case ("H2O")
      atoms = [2.0_dp, 1.0_dp, 0.0_dp]
    case ("HO2")
      atoms = [1.0_dp, 2.0_dp, 0.0_dp]
    case ("H2O2")
      atoms = [2.0_dp, 2.0_dp, 0.0_dp]
    case ("N2")
      atoms = [0.0_dp, 0.0_dp, 2.0_dp]
    case ("AR")
      atoms = 0.0_dp
    case default
      ok = .false.
      return
    end select
    ok = .true.
  end subroutine element_atom_counts

  subroutine reactive_eb_extrema_3d( &
      species, state, temperature, geometry, minimum_density, &
      maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, &
      maximum_closure_error, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: minimum_density, maximum_density
    real(dp), intent(out) :: minimum_pressure, maximum_pressure
    real(dp), intent(out) :: minimum_temperature, maximum_temperature
    real(dp), intent(out) :: maximum_speed, maximum_closure_error
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed, speed, closure_error
    logical :: local_ok
    integer :: i, j, k, species_index, active_cells, nvar

    minimum_density = huge(1.0_dp)
    maximum_density = -huge(1.0_dp)
    minimum_pressure = huge(1.0_dp)
    maximum_pressure = -huge(1.0_dp)
    minimum_temperature = huge(1.0_dp)
    maximum_temperature = -huge(1.0_dp)
    maximum_speed = 0.0_dp
    maximum_closure_error = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid() .or. &
        size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz])) return

    allocate(primitive(reactive_nprim(size(species))))
    active_cells = 0
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          active_cells = active_cells + 1
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive, local_temperature, sound_speed, local_ok)
          if (.not. local_ok) return
          minimum_density = min(minimum_density, primitive(1))
          maximum_density = max(maximum_density, primitive(1))
          minimum_pressure = min(minimum_pressure, primitive(5))
          maximum_pressure = max(maximum_pressure, primitive(5))
          minimum_temperature = min(minimum_temperature, local_temperature)
          maximum_temperature = max(maximum_temperature, local_temperature)
          speed = sqrt(sum(primitive(2:4)**2))
          maximum_speed = max(maximum_speed, speed)
          closure_error = 0.0_dp
          do species_index = 1, size(species)
            closure_error = closure_error + &
              primitive(reactive_mass_fraction_component(species_index))
          end do
          maximum_closure_error = max( &
            maximum_closure_error, abs(closure_error - 1.0_dp))
        end do
      end do
    end do
    ok = active_cells > 0 .and. all(ieee_is_finite([ &
      minimum_density, maximum_density, minimum_pressure, &
      maximum_pressure, minimum_temperature, maximum_temperature, &
      maximum_speed, maximum_closure_error]))
  end subroutine reactive_eb_extrema_3d

  subroutine write_reactive_eb_3d_csv( &
      path, species, geometry, state, temperature, time, ok, message)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :), time
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp), allocatable :: primitive(:)
    real(dp) :: x, y, z, centroid_x, centroid_y, centroid_z
    real(dp) :: local_temperature, sound_speed
    logical :: cell_ok
    integer :: unit, io_status, i, j, k, species_index, nvar

    ok = .false.
    message = ""
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid() .or. &
        .not. ieee_is_finite(time) .or. &
        size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz])) then
      message = "Reactive EB 3D CSV input contract is invalid"
      return
    end if
    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive, local_temperature, sound_speed, cell_ok)
          if (.not. cell_ok) then
            write(message, '(a,i0,a,i0,a,i0,a)') &
              "Invalid output state at (", i, ",", j, ",", k, ")"
            return
          end if
        end do
      end do
    end do

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=io_status)
    if (io_status /= 0) then
      write(message, '(a,1x,a)') "Could not create output file:", trim(path)
      return
    end if
    write(unit, '(a)', advance='no') &
      "time,i,j,k,x,y,z,cell_type,volume_fraction,fluid_centroid_x," // &
      "fluid_centroid_y,fluid_centroid_z,rho,u,v,w,pressure," // &
      "temperature,rhoE,rhou,rhov,rhow"
    do species_index = 1, size(species)
      write(unit, '(a)', advance='no') &
        ",Y_" // trim(species(species_index)%name)
    end do
    do species_index = 1, size(species)
      write(unit, '(a)', advance='no') &
        ",rhoY_" // trim(species(species_index)%name)
    end do
    write(unit, '(a)') ""

    do k = 1, geometry%nz
      z = geometry%z_lower + (real(k, dp) - 0.5_dp) * geometry%dz
      do j = 1, geometry%ny
        y = geometry%y_lower + (real(j, dp) - 0.5_dp) * geometry%dy
        do i = 1, geometry%nx
          x = geometry%x_lower + (real(i, dp) - 0.5_dp) * geometry%dx
          centroid_x = x + geometry%cell_centroid_x(i, j, k) * geometry%dx
          centroid_y = y + geometry%cell_centroid_y(i, j, k) * geometry%dy
          centroid_z = z + geometry%cell_centroid_z(i, j, k) * geometry%dz
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive, local_temperature, sound_speed, cell_ok)
          if (.not. cell_ok) then
            close(unit)
            message = "Reactive EB 3D state changed during CSV output"
            return
          end if
          write(unit, '(*(es25.16e3,:,","))') &
            time, real(i, dp), real(j, dp), real(k, dp), x, y, z, &
            real(geometry%cell_type(i, j, k), dp), &
            geometry%volume_fraction(i, j, k), &
            centroid_x, centroid_y, centroid_z, primitive(1:5), &
            local_temperature, state(iet, i, j, k), &
            state(imx, i, j, k), state(imy, i, j, k), &
            state(imz, i, j, k), &
            (primitive(reactive_mass_fraction_component(species_index)), &
              species_index = 1, size(species)), &
            (state(reactive_species_component(species_index), i, j, k), &
              species_index = 1, size(species))
        end do
      end do
    end do
    close(unit)
    ok = .true.
  end subroutine write_reactive_eb_3d_csv

end module reactive_eb_3d_driver_mod
