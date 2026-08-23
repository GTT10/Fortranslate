module reactive_1d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use state_indices_mod, only: irho, imx, imy, imz, iet, ncons
  use nasa7_thermo_mod, only: nasa7_species, nasa7_mass_properties
  use mixture_thermo_mod, only: &
    valid_mixture_composition, mixture_mass_properties, &
    mixture_specific_gas_constant, mixture_pressure, mixture_density, &
    mixture_sound_speed, temperature_from_internal_energy, &
    mass_fractions_from_mole_fractions, mole_fractions_from_mass_fractions
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use mixture_transport_mod, only: mixture_transport_coefficients
  use constant_volume_reactor_mod, only: advance_constant_volume_adaptive, &
    advance_constant_volume_implicit_adaptive
  use slope_limiter_mod, only: limited_slope, minmod3
  use reconstruction_weno_mod, only: &
    weno_reconstruct_5js, weno_reconstruct_5z, &
    weno_reconstruct_7z, weno_reconstruct_3z
  use simulation_config_reactive_1d_mod, only: &
    reactive_1d_config, reactive_1d_mole_fractions
  implicit none
  private

  integer, parameter, public :: reactive_nwaves = 5
  integer, parameter, public :: reactive_wave_minus = 1
  integer, parameter, public :: reactive_wave_contact = 2
  integer, parameter, public :: reactive_wave_shear_y = 3
  integer, parameter, public :: reactive_wave_shear_z = 4
  integer, parameter, public :: reactive_wave_plus = 5
  integer, parameter :: max_chemistry_substeps = 100000
  real(dp), parameter :: species_tolerance = 5.0e-11_dp
  real(dp), parameter :: contact_steepening_cap = 0.5_dp
  real(dp), parameter :: pi = acos(-1.0_dp)

  public :: reactive_nvar, reactive_nprim
  public :: reactive_species_component, reactive_mass_fraction_component
  public :: reactive_primitive_to_conserved
  public :: reactive_conserved_to_primitive
  public :: reactive_rusanov_flux_x
  public :: reactive_hllc_flux_x
  public :: reactive_riemann_flux_x
  public :: reactive_ppm_interface_value
  public :: reactive_ppm_monotone_edges
  public :: reactive_ppm_reconstruct_five
  public :: reactive_ppm_integrate_profile
  public :: reactive_ppm_flattening_coefficient
  public :: reactive_ppm_contact_steepening_factor
  public :: reactive_ppm_apply_contact_steepening
  public :: build_characteristic_ppm_states
  public :: reconstruct_characteristic_ppm_faces
  public :: reconstruct_ppm_faces
  public :: reactive_difference_to_characteristics
  public :: reactive_characteristics_to_difference
  public :: characteristic_limited_slope
  public :: trace_reactive_characteristics
  public :: initialize_reactive_1d
  public :: advance_reactive_hydro
  public :: reactive_cfl_timestep
  public :: reactive_diffusive_flux_x
  public :: reactive_transport_timestep
  public :: advance_reactive_transport
  public :: advance_reactive_chemistry
  public :: advance_reactive_strang
  public :: simulate_reactive_1d
  public :: write_reactive_1d_csv
  public :: reactive_integrals
  public :: reactive_entropy_wave_density
  public :: reactive_composition_wave_exact

contains

  pure integer function reactive_nvar(nspecies) result(nvar)
    integer, intent(in) :: nspecies
    if (nspecies < 1 .or. nspecies > 32) then
      nvar = 0
    else
      nvar = ncons + nspecies
    end if
  end function reactive_nvar

  pure integer function reactive_nprim(nspecies) result(nvar)
    integer, intent(in) :: nspecies
    if (nspecies < 1 .or. nspecies > 32) then
      nvar = 0
    else
      nvar = 5 + nspecies
    end if
  end function reactive_nprim

  pure integer function reactive_species_component(k) result(index)
    integer, intent(in) :: k
    if (k < 1 .or. k > 32) then
      index = 0
    else
      index = ncons + k
    end if
  end function reactive_species_component

  pure integer function reactive_mass_fraction_component(k) result(index)
    integer, intent(in) :: k
    if (k < 1 .or. k > 32) then
      index = 0
    else
      index = 5 + k
    end if
  end function reactive_mass_fraction_component

  subroutine reactive_primitive_to_conserved( &
      species, primitive, conserved, temperature, sound_speed, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: primitive(:)
    real(dp), intent(out) :: conserved(:), temperature, sound_speed
    logical, intent(out) :: ok

    real(dp), allocatable :: y(:)
    real(dp) :: rho, u, v, w, pressure, r_mix
    real(dp) :: molecular_weight, cp, cv, gamma, enthalpy, energy, entropy
    real(dp) :: kinetic
    integer :: nspecies, k

    conserved = 0.0_dp
    temperature = 0.0_dp
    sound_speed = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (size(primitive) /= reactive_nprim(nspecies)) return
    if (size(conserved) /= reactive_nvar(nspecies)) return

    rho = primitive(1)
    u = primitive(2)
    v = primitive(3)
    w = primitive(4)
    pressure = primitive(5)
    if (rho <= density_floor .or. pressure <= pressure_floor) return
    allocate(y(nspecies))
    do k = 1, nspecies
      y(k) = primitive(reactive_mass_fraction_component(k))
    end do
    if (.not. valid_mixture_composition(species, y)) return
    r_mix = mixture_specific_gas_constant(species, y, ok)
    if (.not. ok .or. r_mix <= 0.0_dp) return
    temperature = pressure / (rho * r_mix)
    call mixture_mass_properties( &
      species, y, temperature, molecular_weight, r_mix, cp, cv, gamma, &
      enthalpy, energy, entropy, ok)
    if (.not. ok) return
    sound_speed = sqrt(gamma * pressure / rho)
    kinetic = 0.5_dp * (u * u + v * v + w * w)

    conserved(irho) = rho
    conserved(imx) = rho * u
    conserved(imy) = rho * v
    conserved(imz) = rho * w
    conserved(iet) = rho * (energy + kinetic)
    do k = 1, nspecies
      conserved(reactive_species_component(k)) = rho * y(k)
    end do
    ok = all(ieee_is_finite(conserved)) .and. &
      ieee_is_finite(temperature) .and. sound_speed > 0.0_dp
  end subroutine reactive_primitive_to_conserved

  subroutine reactive_conserved_to_primitive( &
      species, conserved, temperature_guess, primitive, temperature, &
      sound_speed, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: conserved(:), temperature_guess
    real(dp), intent(out) :: primitive(:), temperature, sound_speed
    logical, intent(out) :: ok

    real(dp), allocatable :: y(:)
    real(dp) :: rho, u, v, w, kinetic_density, target_energy, pressure
    integer :: nspecies, k

    primitive = 0.0_dp
    temperature = 0.0_dp
    sound_speed = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (size(conserved) /= reactive_nvar(nspecies)) return
    if (size(primitive) /= reactive_nprim(nspecies)) return
    if (any(.not. ieee_is_finite(conserved))) return
    rho = conserved(irho)
    if (rho <= density_floor) return
    u = conserved(imx) / rho
    v = conserved(imy) / rho
    w = conserved(imz) / rho
    kinetic_density = 0.5_dp * &
      (conserved(imx)**2 + conserved(imy)**2 + conserved(imz)**2) / rho
    target_energy = (conserved(iet) - kinetic_density) / rho

    allocate(y(nspecies))
    call mass_fractions_from_state(conserved, nspecies, y, ok)
    if (.not. ok) return
    call temperature_from_internal_energy( &
      species, y, target_energy, temperature_guess, temperature, ok)
    if (.not. ok) return
    pressure = mixture_pressure(species, y, rho, temperature, ok)
    if (.not. ok .or. pressure <= pressure_floor) return
    sound_speed = mixture_sound_speed(species, y, temperature, ok)
    if (.not. ok) return

    primitive(1:5) = [rho, u, v, w, pressure]
    do k = 1, nspecies
      primitive(reactive_mass_fraction_component(k)) = y(k)
    end do
    ok = all(ieee_is_finite(primitive))
  end subroutine reactive_conserved_to_primitive

  subroutine mass_fractions_from_state(conserved, nspecies, y, ok)
    real(dp), intent(in) :: conserved(:)
    integer, intent(in) :: nspecies
    real(dp), intent(out) :: y(:)
    logical, intent(out) :: ok
    real(dp) :: rho, total, closure
    integer :: k

    y = 0.0_dp
    ok = .false.
    if (size(conserved) /= reactive_nvar(nspecies) .or. &
        size(y) /= nspecies) return
    rho = conserved(irho)
    if (rho <= density_floor) return
    total = 0.0_dp
    do k = 1, nspecies
      if (conserved(reactive_species_component(k)) < &
          -species_tolerance * max(1.0_dp, rho)) return
      y(k) = max(0.0_dp, conserved(reactive_species_component(k))) / rho
      total = total + conserved(reactive_species_component(k))
    end do
    closure = abs(total - rho) / max(rho, density_floor)
    if (closure > species_tolerance .or. sum(y) <= 0.0_dp) return
    y = y / sum(y)
    ok = all(ieee_is_finite(y))
  end subroutine mass_fractions_from_state

  subroutine reactive_physical_flux_x( &
      species, conserved, temperature_guess, flux, temperature, &
      sound_speed, primitive, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: conserved(:), temperature_guess
    real(dp), intent(out) :: flux(:), temperature, sound_speed, primitive(:)
    logical, intent(out) :: ok
    real(dp) :: mass_flux, pressure, u, v, w
    integer :: nspecies, k

    flux = 0.0_dp
    call reactive_conserved_to_primitive( &
      species, conserved, temperature_guess, primitive, temperature, &
      sound_speed, ok)
    if (.not. ok) return
    nspecies = size(species)
    u = primitive(2)
    v = primitive(3)
    w = primitive(4)
    pressure = primitive(5)
    mass_flux = conserved(irho) * u
    flux(irho) = mass_flux
    flux(imx) = mass_flux * u + pressure
    flux(imy) = mass_flux * v
    flux(imz) = mass_flux * w
    flux(iet) = (conserved(iet) + pressure) * u
    do k = 1, nspecies - 1
      flux(reactive_species_component(k)) = mass_flux * &
        primitive(reactive_mass_fraction_component(k))
    end do
    flux(reactive_species_component(nspecies)) = mass_flux - &
      sum(flux(reactive_species_component(1): &
        reactive_species_component(nspecies - 1)))
  end subroutine reactive_physical_flux_x

  subroutine reactive_rusanov_flux_x( &
      species, left_state, right_state, left_temperature_guess, &
      right_temperature_guess, flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: left_state(:), right_state(:)
    real(dp), intent(in) :: left_temperature_guess, right_temperature_guess
    real(dp), intent(out) :: flux(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: fl(:), fr(:), ql(:), qr(:)
    real(dp) :: tl, tr, cl, cr, spectral_radius, mass_flux
    logical :: left_ok, right_ok
    integer :: nspecies, nvar

    nspecies = size(species)
    nvar = reactive_nvar(nspecies)
    flux = 0.0_dp
    ok = .false.
    if (size(left_state) /= nvar .or. size(right_state) /= nvar .or. &
        size(flux) /= nvar) return
    allocate(fl(nvar), fr(nvar))
    allocate(ql(reactive_nprim(nspecies)), qr(reactive_nprim(nspecies)))
    call reactive_physical_flux_x( &
      species, left_state, left_temperature_guess, fl, tl, cl, ql, left_ok)
    call reactive_physical_flux_x( &
      species, right_state, right_temperature_guess, fr, tr, cr, qr, right_ok)
    if (.not. (left_ok .and. right_ok)) return
    spectral_radius = max(abs(ql(2)) + cl, abs(qr(2)) + cr)
    flux = 0.5_dp * (fl + fr) - &
      0.5_dp * spectral_radius * (right_state - left_state)
    mass_flux = flux(irho)
    flux(reactive_species_component(nspecies)) = mass_flux - &
      sum(flux(reactive_species_component(1): &
        reactive_species_component(nspecies - 1)))
    ok = .true.
  end subroutine reactive_rusanov_flux_x

  subroutine reactive_hllc_flux_x( &
      species, left_state, right_state, left_temperature_guess, &
      right_temperature_guess, flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: left_state(:), right_state(:)
    real(dp), intent(in) :: left_temperature_guess, right_temperature_guess
    real(dp), intent(out) :: flux(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: fl(:), fr(:), ql(:), qr(:), star(:)
    real(dp) :: tl, tr, cl, cr, sl, sr, sm, denominator
    real(dp) :: pstar_left, pstar_right, pstar, scale, mass_flux
    logical :: left_ok, right_ok, local_ok
    integer :: nspecies, nvar

    flux = 0.0_dp
    ok = .false.
    nspecies = size(species)
    nvar = reactive_nvar(nspecies)
    if (size(left_state) /= nvar .or. size(right_state) /= nvar .or. &
        size(flux) /= nvar) return

    allocate(fl(nvar), fr(nvar), star(nvar))
    allocate(ql(reactive_nprim(nspecies)), qr(reactive_nprim(nspecies)))
    call reactive_physical_flux_x( &
      species, left_state, left_temperature_guess, fl, tl, cl, ql, left_ok)
    call reactive_physical_flux_x( &
      species, right_state, right_temperature_guess, fr, tr, cr, qr, right_ok)
    if (.not. (left_ok .and. right_ok)) return

    ! Davis bounds are deliberately conservative. They remain valid for the
    ! frozen-composition ideal-gas mixture because each side supplies its own
    ! thermodynamic sound speed.
    sl = min(ql(2) - cl, qr(2) - cr)
    sr = max(ql(2) + cl, qr(2) + cr)
    if (sl >= 0.0_dp) then
      flux = fl
      ok = .true.
      return
    else if (sr <= 0.0_dp) then
      flux = fr
      ok = .true.
      return
    end if

    denominator = ql(1) * (sl - ql(2)) - qr(1) * (sr - qr(2))
    scale = max(1.0_dp, abs(ql(1) * (sl - ql(2))), &
      abs(qr(1) * (sr - qr(2))))
    if (abs(denominator) <= 100.0_dp * epsilon(1.0_dp) * scale) return
    sm = (qr(5) - ql(5) + ql(1) * ql(2) * (sl - ql(2)) - &
      qr(1) * qr(2) * (sr - qr(2))) / denominator
    if (.not. ieee_is_finite(sm)) return
    if (sm < sl .or. sm > sr) return

    pstar_left = ql(5) + ql(1) * (sl - ql(2)) * (sm - ql(2))
    pstar_right = qr(5) + qr(1) * (sr - qr(2)) * (sm - qr(2))
    pstar = 0.5_dp * (pstar_left + pstar_right)
    if (.not. ieee_is_finite(pstar) .or. pstar <= pressure_floor) return

    if (sm >= 0.0_dp) then
      call build_hllc_star_state( &
        species, left_state, ql, left_temperature_guess, sl, sm, pstar, &
        star, local_ok)
      if (.not. local_ok) return
      flux = fl + sl * (star - left_state)
    else
      call build_hllc_star_state( &
        species, right_state, qr, right_temperature_guess, sr, sm, pstar, &
        star, local_ok)
      if (.not. local_ok) return
      flux = fr + sr * (star - right_state)
    end if

    mass_flux = flux(irho)
    flux(reactive_species_component(nspecies)) = mass_flux - &
      sum(flux(reactive_species_component(1): &
        reactive_species_component(nspecies - 1)))
    ok = all(ieee_is_finite(flux))
  end subroutine reactive_hllc_flux_x

  subroutine build_hllc_star_state( &
      species, state, primitive, temperature_guess, wave_speed, &
      contact_speed, star_pressure, star, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:), primitive(:), temperature_guess
    real(dp), intent(in) :: wave_speed, contact_speed, star_pressure
    real(dp), intent(out) :: star(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: checked_primitive(:)
    real(dp) :: denominator, factor, rho_star, checked_temperature
    real(dp) :: checked_sound_speed
    logical :: local_ok
    integer :: nspecies, k

    star = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (size(state) /= reactive_nvar(nspecies) .or. &
        size(star) /= reactive_nvar(nspecies) .or. &
        size(primitive) /= reactive_nprim(nspecies)) return
    denominator = wave_speed - contact_speed
    if (abs(denominator) <= 100.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(wave_speed), abs(contact_speed))) return
    factor = (wave_speed - primitive(2)) / denominator
    rho_star = primitive(1) * factor
    if (.not. ieee_is_finite(rho_star) .or. rho_star <= density_floor) return

    star(irho) = rho_star
    star(imx) = rho_star * contact_speed
    star(imy) = rho_star * primitive(3)
    star(imz) = rho_star * primitive(4)
    star(iet) = ((wave_speed - primitive(2)) * state(iet) - &
      primitive(5) * primitive(2) + star_pressure * contact_speed) / &
      denominator
    do k = 1, nspecies
      star(reactive_species_component(k)) = rho_star * &
        primitive(reactive_mass_fraction_component(k))
    end do
    star(reactive_species_component(nspecies)) = rho_star - &
      sum(star(reactive_species_component(1): &
        reactive_species_component(nspecies - 1)))
    if (any(.not. ieee_is_finite(star))) return

    ! The approximate star state must still admit a positive-temperature EOS
    ! recovery. The recovered pressure is not forced to the HLLC p-star;
    ! this check is strictly a physical-state gate.
    allocate(checked_primitive(reactive_nprim(nspecies)))
    call reactive_conserved_to_primitive( &
      species, star, temperature_guess, checked_primitive, &
      checked_temperature, checked_sound_speed, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine build_hllc_star_state

  subroutine reactive_riemann_flux_x( &
      species, left_state, right_state, left_temperature_guess, &
      right_temperature_guess, solver, flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: left_state(:), right_state(:)
    real(dp), intent(in) :: left_temperature_guess, right_temperature_guess
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: flux(:)
    logical, intent(out) :: ok

    select case (trim(solver))
    case ("rusanov")
      call reactive_rusanov_flux_x( &
        species, left_state, right_state, left_temperature_guess, &
        right_temperature_guess, flux, ok)
    case ("hllc")
      call reactive_hllc_flux_x( &
        species, left_state, right_state, left_temperature_guess, &
        right_temperature_guess, flux, ok)
    case default
      flux = 0.0_dp
      ok = .false.
    end select
  end subroutine reactive_riemann_flux_x

  pure subroutine reactive_difference_to_characteristics( &
      center, difference, sound_speed, characteristic, ok)
    real(dp), intent(in) :: center(:), difference(:), sound_speed
    real(dp), intent(out) :: characteristic(reactive_nwaves)
    logical, intent(out) :: ok
    real(dp) :: rho, c2

    characteristic = 0.0_dp
    ok = .false.
    if (size(center) < 5 .or. size(difference) /= size(center)) return
    rho = center(1)
    if (rho <= density_floor .or. center(5) <= pressure_floor .or. &
        sound_speed <= 0.0_dp) return
    c2 = sound_speed * sound_speed
    characteristic(reactive_wave_minus) = 0.5_dp * &
      (difference(5) / (rho * sound_speed) - difference(2)) * &
      rho / sound_speed
    characteristic(reactive_wave_plus) = 0.5_dp * &
      (difference(5) / (rho * sound_speed) + difference(2)) * &
      rho / sound_speed
    characteristic(reactive_wave_contact) = difference(1) - difference(5) / c2
    characteristic(reactive_wave_shear_y) = difference(3)
    characteristic(reactive_wave_shear_z) = difference(4)
    ok = .true.
  end subroutine reactive_difference_to_characteristics

  pure subroutine reactive_characteristics_to_difference( &
      center, characteristic, sound_speed, difference, ok)
    real(dp), intent(in) :: center(:), characteristic(reactive_nwaves)
    real(dp), intent(in) :: sound_speed
    real(dp), intent(out) :: difference(:)
    logical, intent(out) :: ok
    real(dp) :: rho, c2

    difference = 0.0_dp
    ok = .false.
    if (size(center) < 5 .or. sound_speed <= 0.0_dp) return
    rho = center(1)
    if (rho <= density_floor .or. center(5) <= pressure_floor) return
    c2 = sound_speed * sound_speed
    difference(1) = characteristic(reactive_wave_minus) + &
      characteristic(reactive_wave_contact) + characteristic(reactive_wave_plus)
    difference(2) = (characteristic(reactive_wave_plus) - &
      characteristic(reactive_wave_minus)) * sound_speed / rho
    difference(3) = characteristic(reactive_wave_shear_y)
    difference(4) = characteristic(reactive_wave_shear_z)
    difference(5) = (characteristic(reactive_wave_minus) + &
      characteristic(reactive_wave_plus)) * c2
    ok = .true.
  end subroutine reactive_characteristics_to_difference

  subroutine characteristic_limited_slope( &
      center, dl, dr, sound_speed, limiter, slope, ok)
    real(dp), intent(in) :: center(:), dl(:), dr(:), sound_speed
    character(len=*), intent(in) :: limiter
    real(dp), intent(out) :: slope(:)
    logical, intent(out) :: ok
    real(dp) :: al(reactive_nwaves), ar(reactive_nwaves)
    real(dp) :: a(reactive_nwaves)
    logical :: local_ok
    integer :: nspecies, wave, k, component

    slope = 0.0_dp
    ok = .false.
    nspecies = size(center) - 5
    call reactive_difference_to_characteristics(center, dl, sound_speed, al, local_ok)
    if (.not. local_ok) return
    call reactive_difference_to_characteristics(center, dr, sound_speed, ar, local_ok)
    if (.not. local_ok) return
    do wave = 1, reactive_nwaves
      call limited_slope(al(wave), ar(wave), limiter, a(wave), local_ok)
      if (.not. local_ok) return
    end do
    call reactive_characteristics_to_difference(center, a, sound_speed, slope, local_ok)
    if (.not. local_ok) return
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      call limited_slope(dl(component), dr(component), limiter, &
        slope(component), local_ok)
      if (.not. local_ok) return
    end do
    ok = .true.
  end subroutine characteristic_limited_slope

  pure subroutine trace_reactive_characteristics( &
      center, slope, sound_speed, dtdx, left_state, right_state, ok)
    real(dp), intent(in) :: center(:), slope(:), sound_speed, dtdx
    real(dp), intent(out) :: left_state(:), right_state(:)
    logical, intent(out) :: ok
    real(dp) :: a(reactive_nwaves), weighted(reactive_nwaves)
    real(dp), allocatable :: advective(:)
    real(dp) :: u
    logical :: local_ok
    integer :: nspecies, k, component

    left_state = 0.0_dp
    right_state = 0.0_dp
    ok = .false.
    if (size(slope) /= size(center) .or. dtdx < 0.0_dp) return
    nspecies = size(center) - 5
    allocate(advective(size(center)))
    advective = 0.0_dp
    call reactive_difference_to_characteristics(center, slope, sound_speed, a, local_ok)
    if (.not. local_ok) return
    u = center(2)
    weighted(reactive_wave_minus) = (u - sound_speed) * a(reactive_wave_minus)
    weighted(reactive_wave_contact) = u * a(reactive_wave_contact)
    weighted(reactive_wave_shear_y) = u * a(reactive_wave_shear_y)
    weighted(reactive_wave_shear_z) = u * a(reactive_wave_shear_z)
    weighted(reactive_wave_plus) = (u + sound_speed) * a(reactive_wave_plus)
    call reactive_characteristics_to_difference( &
      center, weighted, sound_speed, advective, local_ok)
    if (.not. local_ok) return
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      advective(component) = u * slope(component)
    end do
    left_state = center - 0.5_dp * slope - 0.5_dp * dtdx * advective
    right_state = center + 0.5_dp * slope - 0.5_dp * dtdx * advective
    ok = .true.
  end subroutine trace_reactive_characteristics

  subroutine fill_ghosts(state, temperature, nx, boundary, ok)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: boundary
    logical, intent(out) :: ok
    ok = .true.
    select case (trim(boundary))
    case ("periodic")
      state(:, 0) = state(:, nx)
      state(:, nx + 1) = state(:, 1)
      temperature(0) = temperature(nx)
      temperature(nx + 1) = temperature(1)
    case ("outflow")
      state(:, 0) = state(:, 1)
      state(:, nx + 1) = state(:, nx)
      temperature(0) = temperature(1)
      temperature(nx + 1) = temperature(nx)
    case default
      ok = .false.
    end select
  end subroutine fill_ghosts

  pure real(dp) function reactive_ppm_interface_value( &
      value_im1, value_i, value_ip1, value_ip2) result(face_value)
    real(dp), intent(in) :: value_im1, value_i, value_ip1, value_ip2
    real(dp) :: lower, upper

    face_value = (7.0_dp * (value_i + value_ip1) - &
      (value_im1 + value_ip2)) / 12.0_dp
    lower = min(value_i, value_ip1)
    upper = max(value_i, value_ip1)
    face_value = max(lower, min(upper, face_value))
  end function reactive_ppm_interface_value

  pure subroutine reactive_ppm_monotone_edges( &
      cell_value, left_edge, right_edge)
    real(dp), intent(in) :: cell_value
    real(dp), intent(inout) :: left_edge, right_edge
    real(dp) :: delta, curvature

    if ((right_edge - cell_value) * (cell_value - left_edge) <= 0.0_dp) then
      left_edge = cell_value
      right_edge = cell_value
      return
    end if

    delta = right_edge - left_edge
    curvature = 6.0_dp * cell_value - 3.0_dp * (left_edge + right_edge)
    if (delta * curvature > delta * delta) then
      left_edge = 3.0_dp * cell_value - 2.0_dp * right_edge
    else if (delta * curvature < -delta * delta) then
      right_edge = 3.0_dp * cell_value - 2.0_dp * left_edge
    end if
  end subroutine reactive_ppm_monotone_edges

  pure subroutine reactive_ppm_reconstruct_five( &
      stencil, flattening, left_edge, right_edge)
    real(dp), intent(in) :: stencil(5), flattening
    real(dp), intent(out) :: left_edge, right_edge
    real(dp) :: dsl, dsr, dsc, slope_left, slope_right, flat

    flat = max(0.0_dp, min(1.0_dp, flattening))

    dsl = 2.0_dp * (stencil(2) - stencil(1))
    dsr = 2.0_dp * (stencil(3) - stencil(2))
    slope_left = 0.0_dp
    if (dsl * dsr > 0.0_dp) then
      dsc = 0.5_dp * (stencil(3) - stencil(1))
      slope_left = sign(min(abs(dsc), min(abs(dsl), abs(dsr))), dsc)
    end if

    dsl = 2.0_dp * (stencil(3) - stencil(2))
    dsr = 2.0_dp * (stencil(4) - stencil(3))
    slope_right = 0.0_dp
    if (dsl * dsr > 0.0_dp) then
      dsc = 0.5_dp * (stencil(4) - stencil(2))
      slope_right = sign(min(abs(dsc), min(abs(dsl), abs(dsr))), dsc)
    end if

    left_edge = 0.5_dp * (stencil(3) + stencil(2)) - &
      (slope_right - slope_left) / 6.0_dp
    left_edge = max(min(stencil(3), stencil(2)), &
      min(max(stencil(3), stencil(2)), left_edge))

    dsl = 2.0_dp * (stencil(3) - stencil(2))
    dsr = 2.0_dp * (stencil(4) - stencil(3))
    slope_left = 0.0_dp
    if (dsl * dsr > 0.0_dp) then
      dsc = 0.5_dp * (stencil(4) - stencil(2))
      slope_left = sign(min(abs(dsc), min(abs(dsl), abs(dsr))), dsc)
    end if

    dsl = 2.0_dp * (stencil(4) - stencil(3))
    dsr = 2.0_dp * (stencil(5) - stencil(4))
    slope_right = 0.0_dp
    if (dsl * dsr > 0.0_dp) then
      dsc = 0.5_dp * (stencil(5) - stencil(3))
      slope_right = sign(min(abs(dsc), min(abs(dsl), abs(dsr))), dsc)
    end if

    right_edge = 0.5_dp * (stencil(4) + stencil(3)) - &
      (slope_right - slope_left) / 6.0_dp
    right_edge = max(min(stencil(4), stencil(3)), &
      min(max(stencil(4), stencil(3)), right_edge))

    left_edge = flat * left_edge + (1.0_dp - flat) * stencil(3)
    right_edge = flat * right_edge + (1.0_dp - flat) * stencil(3)
    call reactive_ppm_monotone_edges(stencil(3), left_edge, right_edge)
  end subroutine reactive_ppm_reconstruct_five

  pure subroutine reactive_ppm_integrate_profile( &
      left_edge, right_edge, center, velocity, sound_speed, dtdx, &
      integral_right, integral_left, ok)
    real(dp), intent(in) :: left_edge, right_edge, center
    real(dp), intent(in) :: velocity, sound_speed, dtdx
    real(dp), intent(out) :: integral_right(3), integral_left(3)
    logical, intent(out) :: ok
    real(dp) :: speeds(3), speed, sigma, s6
    integer :: wave

    integral_right = 0.0_dp
    integral_left = 0.0_dp
    ok = .false.
    if (sound_speed <= 0.0_dp .or. dtdx < 0.0_dp) return
    speeds = [velocity - sound_speed, velocity, velocity + sound_speed]
    if (maxval(abs(speeds)) * dtdx > 1.0_dp + 50.0_dp * epsilon(1.0_dp)) &
      return
    s6 = 6.0_dp * center - 3.0_dp * (left_edge + right_edge)
    do wave = 1, 3
      speed = speeds(wave)
      sigma = abs(speed) * dtdx
      if (speed <= 0.0_dp) then
        integral_right(wave) = right_edge
        integral_left(wave) = left_edge + 0.5_dp * sigma * &
          (right_edge - left_edge + &
            (1.0_dp - 2.0_dp * sigma / 3.0_dp) * s6)
      else
        integral_right(wave) = right_edge - 0.5_dp * sigma * &
          (right_edge - left_edge - &
            (1.0_dp - 2.0_dp * sigma / 3.0_dp) * s6)
        integral_left(wave) = left_edge
      end if
    end do
    ok = all(ieee_is_finite(integral_right)) .and. &
      all(ieee_is_finite(integral_left))
  end subroutine reactive_ppm_integrate_profile

  pure real(dp) function reactive_ppm_flattening_coefficient( &
      pressure, velocity) result(flattening)
    real(dp), intent(in) :: pressure(-3:3), velocity(-3:3)
    real(dp), parameter :: shock_threshold = 0.33_dp
    real(dp), parameter :: zcut1 = 0.75_dp
    real(dp), parameter :: zcut2 = 0.85_dp
    real(dp) :: dp_jump, denominator, zeta, z, z2
    real(dp) :: shifted_plus, shifted_minus, shifted_plus2, shifted_minus2
    real(dp) :: shifted_uplus, shifted_uminus, minimum_pressure
    real(dp) :: chi, chi2, compression
    integer :: shift

    flattening = 1.0_dp
    if (any(pressure <= pressure_floor) .or. &
        any(.not. ieee_is_finite(pressure)) .or. &
        any(.not. ieee_is_finite(velocity))) return

    dp_jump = pressure(1) - pressure(-1)
    shift = merge(1, -1, dp_jump > 0.0_dp)
    denominator = max(tiny(1.0_dp), abs(pressure(2) - pressure(-2)))
    zeta = abs(dp_jump) / denominator
    z = max(0.0_dp, min(1.0_dp, &
      (zeta - zcut1) / (zcut2 - zcut1)))
    compression = merge(1.0_dp, 0.0_dp, velocity(-1) - velocity(1) >= 0.0_dp)
    minimum_pressure = max(pressure_floor, min(pressure(1), pressure(-1)))
    chi = merge(compression, 0.0_dp, &
      abs(dp_jump) / minimum_pressure > shock_threshold)

    shifted_plus = pressure(1 - shift)
    shifted_minus = pressure(-(1 + shift))
    shifted_uplus = velocity(1 - shift)
    shifted_uminus = velocity(-(1 + shift))
    shifted_plus2 = pressure(2 - shift)
    shifted_minus2 = pressure(-(2 + shift))
    dp_jump = shifted_plus - shifted_minus
    denominator = max(tiny(1.0_dp), abs(shifted_plus2 - shifted_minus2))
    zeta = abs(dp_jump) / denominator
    z2 = max(0.0_dp, min(1.0_dp, &
      (zeta - zcut1) / (zcut2 - zcut1)))
    compression = merge(1.0_dp, 0.0_dp, &
      shifted_uminus - shifted_uplus >= 0.0_dp)
    minimum_pressure = max(pressure_floor, min(shifted_plus, shifted_minus))
    chi2 = merge(compression, 0.0_dp, &
      abs(dp_jump) / minimum_pressure > shock_threshold)

    flattening = 1.0_dp - max(chi * z, chi2 * z2)
    flattening = max(0.0_dp, min(1.0_dp, flattening))
  end function reactive_ppm_flattening_coefficient

  pure real(dp) function reactive_ppm_contact_steepening_factor( &
      density, pressure, gamma_effective) result(eta)
    real(dp), intent(in) :: density(-2:2), pressure(-2:2)
    real(dp), intent(in) :: gamma_effective
    real(dp), parameter :: contact_k0 = 0.1_dp
    real(dp), parameter :: eta1 = 20.0_dp
    real(dp), parameter :: eta2 = 0.05_dp
    real(dp), parameter :: density_threshold = 0.01_dp
    real(dp) :: d1rho, d2rho_minus, d2rho_plus
    real(dp) :: minimum_density, minimum_pressure, eta_tilde
    logical :: contact_check, curvature_check, change_check

    eta = 0.0_dp
    if (gamma_effective <= 0.0_dp .or. any(density <= density_floor) .or. &
        any(pressure <= pressure_floor)) return
    d1rho = density(1) - density(-1)
    d2rho_minus = density(0) - 2.0_dp * density(-1) + density(-2)
    d2rho_plus = density(2) - 2.0_dp * density(1) + density(0)
    minimum_density = max(density_floor, min(density(1), density(-1)))
    minimum_pressure = max(pressure_floor, min(pressure(1), pressure(-1)))
    contact_check = gamma_effective * contact_k0 * abs(d1rho) * &
      minimum_pressure >= abs(pressure(1) - pressure(-1)) * minimum_density
    curvature_check = d2rho_plus * d2rho_minus <= 0.0_dp
    change_check = abs(d1rho) >= density_threshold * minimum_density
    if (.not. (contact_check .and. curvature_check .and. change_check)) return
    if (abs(d1rho) <= tiny(1.0_dp)) return
    eta_tilde = -(d2rho_plus - d2rho_minus) / (6.0_dp * d1rho)
    eta = max(0.0_dp, min(1.0_dp, eta1 * (eta_tilde - eta2)))
  end function reactive_ppm_contact_steepening_factor

  pure subroutine reactive_ppm_apply_contact_steepening( &
      stencil, eta, left_edge, right_edge)
    real(dp), intent(in) :: stencil(-2:2), eta
    real(dp), intent(inout) :: left_edge, right_edge
    real(dp) :: slope_minus, slope_plus, left_mc, right_mc, strength

    strength = max(0.0_dp, min(1.0_dp, eta))
    slope_minus = minmod3( &
      0.5_dp * ((stencil(-1) - stencil(-2)) + &
        (stencil(0) - stencil(-1))), &
      2.0_dp * (stencil(-1) - stencil(-2)), &
      2.0_dp * (stencil(0) - stencil(-1)))
    slope_plus = minmod3( &
      0.5_dp * ((stencil(1) - stencil(0)) + &
        (stencil(2) - stencil(1))), &
      2.0_dp * (stencil(1) - stencil(0)), &
      2.0_dp * (stencil(2) - stencil(1)))
    left_mc = stencil(-1) + 0.5_dp * slope_minus
    right_mc = stencil(1) - 0.5_dp * slope_plus
    left_mc = max(min(stencil(-1), stencil(0)), &
      min(max(stencil(-1), stencil(0)), left_mc))
    right_mc = max(min(stencil(0), stencil(1)), &
      min(max(stencil(0), stencil(1)), right_mc))
    left_edge = (1.0_dp - strength) * left_edge + strength * left_mc
    right_edge = (1.0_dp - strength) * right_edge + strength * right_mc
    left_edge = max(min(stencil(-1), stencil(0)), &
      min(max(stencil(-1), stencil(0)), left_edge))
    right_edge = max(min(stencil(0), stencil(1)), &
      min(max(stencil(0), stencil(1)), right_edge))
  end subroutine reactive_ppm_apply_contact_steepening

  pure integer function extended_cell_index(index, nx, boundary) result(source)
    integer, intent(in) :: index, nx
    character(len=*), intent(in) :: boundary

    select case (trim(boundary))
    case ("periodic")
      source = modulo(index - 1, nx) + 1
    case ("outflow")
      source = max(1, min(nx, index))
    case default
      source = 0
    end select
  end function extended_cell_index

  subroutine build_characteristic_ppm_states( &
      species, center, sound_speed, integral_right, integral_left, &
      left_state, right_state, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: center(:), sound_speed
    real(dp), intent(in) :: integral_right(:, :), integral_left(:, :)
    real(dp), intent(out) :: left_state(:), right_state(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: reference(:), dummy_conserved(:)
    real(dp) :: rho_ref, u_ref, p_ref, c_ref, t_ref
    real(dp) :: dum, dpm, drho, dp0, dup, dpp
    real(dp) :: alpham, alpha0, alphap, c2, u
    logical :: local_ok
    integer :: nspecies, nprimitive, k, component

    left_state = 0.0_dp
    right_state = 0.0_dp
    ok = .false.
    nspecies = size(species)
    nprimitive = reactive_nprim(nspecies)
    if (size(center) /= nprimitive .or. size(left_state) /= nprimitive .or. &
        size(right_state) /= nprimitive) return
    if (size(integral_right, 1) /= nprimitive .or. &
        size(integral_left, 1) /= nprimitive .or. &
        size(integral_right, 2) /= 3 .or. &
        size(integral_left, 2) /= 3 .or. sound_speed <= 0.0_dp) return

    allocate(reference(nprimitive), dummy_conserved(reactive_nvar(nspecies)))
    u = center(2)

    left_state = integral_left(:, 2)
    reference = integral_left(:, 1)
    call sanitize_primitive(reference, center, nspecies)
    call reactive_primitive_to_conserved( &
      species, reference, dummy_conserved, t_ref, c_ref, local_ok)
    if (.not. local_ok) return
    rho_ref = reference(1)
    u_ref = reference(2)
    p_ref = reference(5)
    c2 = c_ref * c_ref
    dum = u_ref - integral_left(2, 1)
    dpm = p_ref - integral_left(5, 1)
    drho = rho_ref - integral_left(1, 2)
    dp0 = p_ref - integral_left(5, 2)
    dup = u_ref - integral_left(2, 3)
    dpp = p_ref - integral_left(5, 3)
    alpham = 0.5_dp * (dpm / (rho_ref * c_ref) - dum) * rho_ref / c_ref
    alphap = 0.5_dp * (dpp / (rho_ref * c_ref) + dup) * rho_ref / c_ref
    alpha0 = drho - dp0 / c2
    if (u - sound_speed > 0.0_dp) then
      alpham = 0.0_dp
    else
      alpham = -alpham
    end if
    if (u + sound_speed > 0.0_dp) then
      alphap = 0.0_dp
    else
      alphap = -alphap
    end if
    if (u > 0.0_dp) then
      alpha0 = 0.0_dp
    else
      alpha0 = -alpha0
    end if
    left_state(1) = rho_ref + alpham + alpha0 + alphap
    left_state(2) = u_ref + (alphap - alpham) * c_ref / rho_ref
    left_state(3) = integral_left(3, 2)
    left_state(4) = integral_left(4, 2)
    left_state(5) = p_ref + (alpham + alphap) * c2
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      left_state(component) = integral_left(component, 2)
    end do
    call sanitize_primitive(left_state, center, nspecies)

    right_state = integral_right(:, 2)
    reference = integral_right(:, 3)
    call sanitize_primitive(reference, center, nspecies)
    call reactive_primitive_to_conserved( &
      species, reference, dummy_conserved, t_ref, c_ref, local_ok)
    if (.not. local_ok) return
    rho_ref = reference(1)
    u_ref = reference(2)
    p_ref = reference(5)
    c2 = c_ref * c_ref
    dum = u_ref - integral_right(2, 1)
    dpm = p_ref - integral_right(5, 1)
    drho = rho_ref - integral_right(1, 2)
    dp0 = p_ref - integral_right(5, 2)
    dup = u_ref - integral_right(2, 3)
    dpp = p_ref - integral_right(5, 3)
    alpham = 0.5_dp * (dpm / (rho_ref * c_ref) - dum) * rho_ref / c_ref
    alphap = 0.5_dp * (dpp / (rho_ref * c_ref) + dup) * rho_ref / c_ref
    alpha0 = drho - dp0 / c2
    if (u - sound_speed > 0.0_dp) then
      alpham = -alpham
    else
      alpham = 0.0_dp
    end if
    if (u + sound_speed > 0.0_dp) then
      alphap = -alphap
    else
      alphap = 0.0_dp
    end if
    if (u > 0.0_dp) then
      alpha0 = -alpha0
    else
      alpha0 = 0.0_dp
    end if
    right_state(1) = rho_ref + alpham + alpha0 + alphap
    right_state(2) = u_ref + (alphap - alpham) * c_ref / rho_ref
    right_state(3) = integral_right(3, 2)
    right_state(4) = integral_right(4, 2)
    right_state(5) = p_ref + (alpham + alphap) * c2
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      right_state(component) = integral_right(component, 2)
    end do
    call sanitize_primitive(right_state, center, nspecies)
    ok = all(ieee_is_finite(left_state)) .and. &
      all(ieee_is_finite(right_state))
  end subroutine build_characteristic_ppm_states

  subroutine reconstruct_characteristic_ppm_faces( &
      species, state, temperature, nx, boundary, dtdx, &
      use_contact_steepening, use_shock_flattening, left_faces, right_faces, &
      left_t, right_t, ok, left_ghost_state, right_ghost_state, &
      left_ghost_temperature, right_ghost_temperature, hybrid_weno, &
      weno_scheme)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: boundary
    real(dp), intent(in) :: dtdx
    logical, intent(in) :: use_contact_steepening, use_shock_flattening
    real(dp), intent(out) :: left_faces(:, 0:), right_faces(:, 0:)
    real(dp), intent(out) :: left_t(0:), right_t(0:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: left_ghost_state(:, :)
    real(dp), intent(in), optional :: right_ghost_state(:, :)
    real(dp), intent(in), optional :: left_ghost_temperature(:)
    real(dp), intent(in), optional :: right_ghost_temperature(:)
    logical, intent(in), optional :: hybrid_weno
    integer, intent(in), optional :: weno_scheme

    real(dp), allocatable :: q(:, :), cell_left(:, :), cell_right(:, :)
    real(dp), allocatable :: integral_right(:, :), integral_left(:, :)
    real(dp), allocatable :: edge_left(:), edge_right(:), sound_speed(:)
    real(dp) :: stencil(5), stencil7(7), stencil3(3)
    real(dp) :: contact_stencil(-2:2)
    real(dp) :: density_stencil(-2:2), pressure_stencil(-2:2)
    real(dp) :: pressure_wide(-3:3), velocity_wide(-3:3)
    real(dp) :: flattening, eta, gamma_effective, dummy_c, local_t
    real(dp) :: profile_right(3), profile_left(3)
    logical :: local_ok, wide_ghosts, use_weno
    integer :: nspecies, nprimitive, i, j, component, k, source, lc, rc
    integer :: first_cell, last_cell, layer, selected_weno_scheme

    ok = .false.
    if (dtdx < 0.0_dp) return
    use_weno = .false.
    selected_weno_scheme = 1
    if (present(hybrid_weno)) use_weno = hybrid_weno
    if (present(weno_scheme)) selected_weno_scheme = weno_scheme
    if (use_weno .and. &
        (selected_weno_scheme < 0 .or. selected_weno_scheme > 3)) return
    wide_ghosts = present(left_ghost_state) .and. &
      present(right_ghost_state) .and. &
      present(left_ghost_temperature) .and. &
      present(right_ghost_temperature)
    if (wide_ghosts .neqv. (present(left_ghost_state) .or. &
        present(right_ghost_state) .or. &
        present(left_ghost_temperature) .or. &
        present(right_ghost_temperature))) return
    nspecies = size(species)
    nprimitive = reactive_nprim(nspecies)
    if (wide_ghosts) then
      if (size(left_ghost_state, 1) /= size(state, 1) .or. &
          size(right_ghost_state, 1) /= size(state, 1) .or. &
          size(left_ghost_state, 2) /= 4 .or. &
          size(right_ghost_state, 2) /= 4 .or. &
          size(left_ghost_temperature) /= 4 .or. &
          size(right_ghost_temperature) /= 4) return
      first_cell = 0
      last_cell = nx + 1
      allocate(q(nprimitive, -3:nx + 4), sound_speed(0:nx + 1))
      allocate(cell_left(nprimitive, 0:nx + 1))
      allocate(cell_right(nprimitive, 0:nx + 1))
    else
      first_cell = 1
      last_cell = nx
      allocate(q(nprimitive, -3:nx + 3), sound_speed(1:nx))
      allocate(cell_left(nprimitive, 1:nx), cell_right(nprimitive, 1:nx))
    end if
    allocate(integral_right(nprimitive, 3), integral_left(nprimitive, 3))
    allocate(edge_left(nprimitive), edge_right(nprimitive))

    do i = 1, nx
      call reactive_conserved_to_primitive( &
        species, state(:, i), temperature(i), q(:, i), local_t, &
        sound_speed(i), local_ok)
      if (.not. local_ok) return
    end do
    if (wide_ghosts) then
      do layer = 1, 4
        i = 1 - layer
        call reactive_conserved_to_primitive( &
          species, left_ghost_state(:, layer), &
          left_ghost_temperature(layer), q(:, i), local_t, dummy_c, local_ok)
        if (.not. local_ok) return
        i = nx + layer
        call reactive_conserved_to_primitive( &
          species, right_ghost_state(:, layer), &
          right_ghost_temperature(layer), q(:, i), local_t, dummy_c, local_ok)
        if (.not. local_ok) return
      end do
    else
      do i = -3, nx + 3
        if (i >= 1 .and. i <= nx) cycle
        source = extended_cell_index(i, nx, boundary)
        if (source == 0) return
        q(:, i) = q(:, source)
      end do
    end if

    do i = first_cell, last_cell
      if (wide_ghosts .and. (i == 0 .or. i == nx + 1)) then
        call reactive_conserved_to_primitive( &
          species, merge(left_ghost_state(:, 1), &
            right_ghost_state(:, 1), i == 0), &
          merge(left_ghost_temperature(1), right_ghost_temperature(1), &
            i == 0), q(:, i), local_t, sound_speed(i), local_ok)
        if (.not. local_ok) return
      end if
      flattening = 1.0_dp
      if (use_shock_flattening) then
        do j = -3, 3
          pressure_wide(j) = q(5, i + j)
          velocity_wide(j) = q(2, i + j)
        end do
        flattening = reactive_ppm_flattening_coefficient( &
          pressure_wide, velocity_wide)
      end if
      do component = 1, nprimitive
        stencil = q(component, i - 2:i + 2)
        if (use_weno) then
          select case (selected_weno_scheme)
          case (0)
            call weno_reconstruct_5js( &
              stencil, edge_left(component), edge_right(component))
          case (1)
            call weno_reconstruct_5z( &
              stencil, edge_left(component), edge_right(component))
          case (2)
            stencil7 = q(component, i - 3:i + 3)
            call weno_reconstruct_7z( &
              stencil7, edge_left(component), edge_right(component))
          case (3)
            stencil3 = q(component, i - 1:i + 1)
            call weno_reconstruct_3z( &
              stencil3, edge_left(component), edge_right(component))
          end select
        else
          call reactive_ppm_reconstruct_five( &
            stencil, flattening, edge_left(component), edge_right(component))
        end if
      end do

      if (use_contact_steepening .and. flattening > 0.999_dp) then
        do j = -2, 2
          density_stencil(j) = q(1, i + j)
          pressure_stencil(j) = q(5, i + j)
        end do
        gamma_effective = sound_speed(i)**2 * q(1, i) / q(5, i)
        ! This subset bounds the canonical Colella--Woodward detector at
        ! half strength.  Full-strength density and composition steepening
        ! can over-compress a material interface when coupled to the current
        ! frozen-composition HLLC star-state construction.
        eta = min(contact_steepening_cap, &
          reactive_ppm_contact_steepening_factor( &
          density_stencil, pressure_stencil, gamma_effective))
        if (eta > 0.0_dp) then
          call reactive_ppm_apply_contact_steepening( &
            density_stencil, eta, edge_left(1), edge_right(1))
          do k = 1, nspecies
            component = reactive_mass_fraction_component(k)
            do j = -2, 2
              contact_stencil(j) = q(component, i + j)
            end do
            call reactive_ppm_apply_contact_steepening( &
              contact_stencil, eta, edge_left(component), &
              edge_right(component))
          end do
        end if
      end if

      call sanitize_primitive(edge_left, q(:, i), nspecies)
      call sanitize_primitive(edge_right, q(:, i), nspecies)
      do component = 1, nprimitive
        call reactive_ppm_integrate_profile( &
          edge_left(component), edge_right(component), q(component, i), &
          q(2, i), sound_speed(i), dtdx, profile_right, profile_left, local_ok)
        if (.not. local_ok) return
        integral_right(component, :) = profile_right
        integral_left(component, :) = profile_left
      end do
      call build_characteristic_ppm_states( &
        species, q(:, i), sound_speed(i), integral_right, integral_left, &
        cell_left(:, i), cell_right(:, i), local_ok)
      if (.not. local_ok) return
    end do

    do j = 0, nx
      if (wide_ghosts) then
        lc = j
        rc = j + 1
      else if (j == 0) then
        if (trim(boundary) == "periodic") then
          lc = nx
          rc = 1
        else
          call reactive_primitive_to_conserved( &
            species, q(:, 1), left_faces(:, j), left_t(j), dummy_c, local_ok)
          if (.not. local_ok) return
          right_faces(:, j) = left_faces(:, j)
          right_t(j) = left_t(j)
          cycle
        end if
      else if (j == nx) then
        if (trim(boundary) == "periodic") then
          lc = nx
          rc = 1
        else
          call reactive_primitive_to_conserved( &
            species, q(:, nx), left_faces(:, j), left_t(j), dummy_c, local_ok)
          if (.not. local_ok) return
          right_faces(:, j) = left_faces(:, j)
          right_t(j) = left_t(j)
          cycle
        end if
      else
        lc = j
        rc = j + 1
      end if
      call reactive_primitive_to_conserved( &
        species, cell_right(:, lc), left_faces(:, j), left_t(j), dummy_c, &
        local_ok)
      if (.not. local_ok) return
      call reactive_primitive_to_conserved( &
        species, cell_left(:, rc), right_faces(:, j), right_t(j), dummy_c, &
        local_ok)
      if (.not. local_ok) return
    end do
    ok = .true.
  end subroutine reconstruct_characteristic_ppm_faces

  subroutine reconstruct_ppm_faces( &
      species, state, temperature, nx, boundary, left_faces, right_faces, &
      left_t, right_t, ok, left_ghost_state, right_ghost_state, &
      left_ghost_temperature, right_ghost_temperature)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: boundary
    real(dp), intent(out) :: left_faces(:, 0:), right_faces(:, 0:)
    real(dp), intent(out) :: left_t(0:), right_t(0:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: left_ghost_state(:, :)
    real(dp), intent(in), optional :: right_ghost_state(:, :)
    real(dp), intent(in), optional :: left_ghost_temperature(:)
    real(dp), intent(in), optional :: right_ghost_temperature(:)

    real(dp), allocatable :: q(:, :), qleft(:, :), qright(:, :)
    real(dp), allocatable :: face(:, :)
    real(dp) :: dummy_c, local_temperature
    logical :: local_ok, wide_ghosts
    integer :: nspecies, nprimitive, i, j, component, source, lc, rc
    integer :: first_cell, last_cell, layer

    ok = .false.
    wide_ghosts = present(left_ghost_state) .and. &
      present(right_ghost_state) .and. &
      present(left_ghost_temperature) .and. &
      present(right_ghost_temperature)
    if (wide_ghosts .neqv. (present(left_ghost_state) .or. &
        present(right_ghost_state) .or. &
        present(left_ghost_temperature) .or. &
        present(right_ghost_temperature))) return
    nspecies = size(species)
    nprimitive = reactive_nprim(nspecies)
    if (wide_ghosts) then
      if (size(left_ghost_state, 1) /= size(state, 1) .or. &
          size(right_ghost_state, 1) /= size(state, 1) .or. &
          size(left_ghost_state, 2) /= 4 .or. &
          size(right_ghost_state, 2) /= 4 .or. &
          size(left_ghost_temperature) /= 4 .or. &
          size(right_ghost_temperature) /= 4) return
      first_cell = 0
      last_cell = nx + 1
      allocate(q(nprimitive, -3:nx + 4))
      allocate(qleft(nprimitive, 0:nx + 1))
      allocate(qright(nprimitive, 0:nx + 1))
    else
      first_cell = 1
      last_cell = nx
      allocate(q(nprimitive, -2:nx + 3))
      allocate(qleft(nprimitive, 1:nx), qright(nprimitive, 1:nx))
    end if
    allocate(face(nprimitive, -1:nx + 1))

    do i = 1, nx
      call reactive_conserved_to_primitive( &
        species, state(:, i), temperature(i), q(:, i), local_temperature, &
        dummy_c, local_ok)
      if (.not. local_ok) return
    end do
    if (wide_ghosts) then
      do layer = 1, 4
        i = 1 - layer
        call reactive_conserved_to_primitive( &
          species, left_ghost_state(:, layer), &
          left_ghost_temperature(layer), q(:, i), local_temperature, &
          dummy_c, local_ok)
        if (.not. local_ok) return
        i = nx + layer
        call reactive_conserved_to_primitive( &
          species, right_ghost_state(:, layer), &
          right_ghost_temperature(layer), q(:, i), local_temperature, &
          dummy_c, local_ok)
        if (.not. local_ok) return
      end do
    else
      do i = -2, nx + 3
        if (i >= 1 .and. i <= nx) cycle
        source = extended_cell_index(i, nx, boundary)
        if (source == 0) return
        q(:, i) = q(:, source)
      end do
    end if

    do j = -1, nx + 1
      do component = 1, nprimitive
        face(component, j) = reactive_ppm_interface_value( &
          q(component, j - 1), q(component, j), q(component, j + 1), &
          q(component, j + 2))
      end do
    end do

    do i = first_cell, last_cell
      qleft(:, i) = face(:, i - 1)
      qright(:, i) = face(:, i)
      do component = 1, nprimitive
        call reactive_ppm_monotone_edges( &
          q(component, i), qleft(component, i), qright(component, i))
      end do
      call sanitize_primitive(qleft(:, i), q(:, i), nspecies)
      call sanitize_primitive(qright(:, i), q(:, i), nspecies)
    end do

    do j = 0, nx
      if (wide_ghosts) then
        lc = j
        rc = j + 1
      else if (j == 0) then
        if (trim(boundary) == "periodic") then
          lc = nx
          rc = 1
        else
          call reactive_primitive_to_conserved( &
            species, q(:, 1), left_faces(:, j), left_t(j), dummy_c, local_ok)
          if (.not. local_ok) return
          right_faces(:, j) = left_faces(:, j)
          right_t(j) = left_t(j)
          cycle
        end if
      else if (j == nx) then
        if (trim(boundary) == "periodic") then
          lc = nx
          rc = 1
        else
          call reactive_primitive_to_conserved( &
            species, q(:, nx), left_faces(:, j), left_t(j), dummy_c, local_ok)
          if (.not. local_ok) return
          right_faces(:, j) = left_faces(:, j)
          right_t(j) = left_t(j)
          cycle
        end if
      else
        lc = j
        rc = j + 1
      end if
      call reactive_primitive_to_conserved( &
        species, qright(:, lc), left_faces(:, j), left_t(j), dummy_c, local_ok)
      if (.not. local_ok) return
      call reactive_primitive_to_conserved( &
        species, qleft(:, rc), right_faces(:, j), right_t(j), dummy_c, local_ok)
      if (.not. local_ok) return
    end do
    ok = .true.
  end subroutine reconstruct_ppm_faces

  subroutine reconstruct_faces( &
      species, state, temperature, nx, reconstruction, limiter, boundary, &
      dtdx, use_contact_steepening, use_shock_flattening, left_faces, &
      right_faces, left_t, right_t, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: reconstruction, limiter, boundary
    real(dp), intent(in) :: dtdx
    logical, intent(in) :: use_contact_steepening, use_shock_flattening
    real(dp), intent(out) :: left_faces(:, 0:), right_faces(:, 0:)
    real(dp), intent(out) :: left_t(0:), right_t(0:)
    logical, intent(out) :: ok
    real(dp), allocatable :: q(:, :), c(:), slope(:, :)
    real(dp), allocatable :: qleft(:, :), qright(:, :), temp_cell(:)
    real(dp), allocatable :: dl(:), dr(:)
    real(dp) :: theta, dummy_c
    logical :: local_ok
    integer :: nspecies, nprimitive, i, lc, rc

    ok = .false.
    nspecies = size(species)
    nprimitive = reactive_nprim(nspecies)
    if (trim(reconstruction) == "pcm") then
      do i = 0, nx
        left_faces(:, i) = state(:, i)
        right_faces(:, i) = state(:, i + 1)
        left_t(i) = temperature(i)
        right_t(i) = temperature(i + 1)
      end do
      ok = .true.
      return
    end if
    if (trim(reconstruction) == "ppm") then
      call reconstruct_ppm_faces( &
        species, state, temperature, nx, boundary, left_faces, right_faces, &
        left_t, right_t, ok)
      return
    end if
    if (trim(reconstruction) == "characteristic_ppm") then
      call reconstruct_characteristic_ppm_faces( &
        species, state, temperature, nx, boundary, dtdx, &
        use_contact_steepening, use_shock_flattening, left_faces, &
        right_faces, left_t, right_t, ok)
      return
    end if
    if (trim(reconstruction) /= "characteristic_plm") return

    allocate(q(nprimitive, 0:nx + 1), c(0:nx + 1), temp_cell(0:nx + 1))
    allocate(slope(nprimitive, 1:nx), qleft(nprimitive, 1:nx))
    allocate(qright(nprimitive, 1:nx), dl(nprimitive), dr(nprimitive))
    do i = 0, nx + 1
      call reactive_conserved_to_primitive( &
        species, state(:, i), temperature(i), q(:, i), temp_cell(i), c(i), local_ok)
      if (.not. local_ok) return
    end do
    do i = 1, nx
      dl = q(:, i) - q(:, i - 1)
      dr = q(:, i + 1) - q(:, i)
      call characteristic_limited_slope(q(:, i), dl, dr, c(i), limiter, &
        slope(:, i), local_ok)
      if (.not. local_ok) return
      theta = primitive_slope_scale(q(:, i), slope(:, i), nspecies)
      slope(:, i) = theta * slope(:, i)
      call trace_reactive_characteristics(q(:, i), slope(:, i), c(i), dtdx, &
        qleft(:, i), qright(:, i), local_ok)
      if (.not. local_ok) return
      call sanitize_primitive(qleft(:, i), q(:, i), nspecies)
      call sanitize_primitive(qright(:, i), q(:, i), nspecies)
    end do

    do i = 0, nx
      if (i == 0) then
        if (trim(boundary) == "periodic") then
          lc = nx; rc = 1
        else
          lc = 1; rc = 1
        end if
      else if (i == nx) then
        if (trim(boundary) == "periodic") then
          lc = nx; rc = 1
        else
          lc = nx; rc = nx
        end if
      else
        lc = i; rc = i + 1
      end if
      call reactive_primitive_to_conserved( &
        species, qright(:, lc), left_faces(:, i), left_t(i), dummy_c, local_ok)
      if (.not. local_ok) return
      call reactive_primitive_to_conserved( &
        species, qleft(:, rc), right_faces(:, i), right_t(i), dummy_c, local_ok)
      if (.not. local_ok) return
    end do
    ok = .true.
  end subroutine reconstruct_faces

  pure real(dp) function primitive_slope_scale(center, slope, nspecies) result(theta)
    real(dp), intent(in) :: center(:), slope(:)
    integer, intent(in) :: nspecies
    integer :: k, component
    theta = 1.0_dp
    theta = min(theta, lower_scale(center(1), slope(1), density_floor))
    theta = min(theta, lower_scale(center(5), slope(5), pressure_floor))
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      theta = min(theta, lower_scale(center(component), slope(component), 0.0_dp))
      theta = min(theta, upper_scale(center(component), slope(component), 1.0_dp))
    end do
  end function primitive_slope_scale

  pure real(dp) function lower_scale(center, slope, lower) result(theta)
    real(dp), intent(in) :: center, slope, lower
    real(dp) :: magnitude
    magnitude = abs(slope)
    if (magnitude <= tiny(1.0_dp) .or. center - magnitude > lower) then
      theta = 1.0_dp
    else
      theta = max(0.0_dp, min(1.0_dp, (center - lower) / magnitude))
    end if
  end function lower_scale

  pure real(dp) function upper_scale(center, slope, upper) result(theta)
    real(dp), intent(in) :: center, slope, upper
    real(dp) :: magnitude
    magnitude = abs(slope)
    if (magnitude <= tiny(1.0_dp) .or. center + magnitude < upper) then
      theta = 1.0_dp
    else
      theta = max(0.0_dp, min(1.0_dp, (upper - center) / magnitude))
    end if
  end function upper_scale

  pure subroutine sanitize_primitive(q, fallback, nspecies)
    real(dp), intent(inout) :: q(:)
    real(dp), intent(in) :: fallback(:)
    integer, intent(in) :: nspecies
    real(dp) :: total
    integer :: k, component
    if (q(1) <= density_floor .or. q(5) <= pressure_floor) then
      q = fallback
      return
    end if
    total = 0.0_dp
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      if (q(component) < -1.0e-12_dp) then
        q = fallback
        return
      end if
      q(component) = max(0.0_dp, q(component))
      total = total + q(component)
    end do
    if (total <= tiny(1.0_dp)) then
      q = fallback
      return
    end if
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      q(component) = q(component) / total
    end do
  end subroutine sanitize_primitive

  subroutine reactive_cfl_timestep(species, state, temperature, nx, dx, cfl, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok
    real(dp), allocatable :: q(:)
    real(dp) :: local_t, c, maximum_speed
    integer :: i
    logical :: local_ok
    dt = 0.0_dp
    ok = .false.
    allocate(q(reactive_nprim(size(species))))
    maximum_speed = 0.0_dp
    do i = 1, nx
      call reactive_conserved_to_primitive( &
        species, state(:, i), temperature(i), q, local_t, c, local_ok)
      if (.not. local_ok) return
      maximum_speed = max(maximum_speed, abs(q(2)) + c)
    end do
    if (maximum_speed <= 0.0_dp) return
    dt = cfl * dx / maximum_speed
    ok = dt > 0.0_dp
  end subroutine reactive_cfl_timestep


  subroutine reactive_diffusive_flux_x( &
      species, transport, left_state, right_state, left_temperature_guess, &
      right_temperature_guess, dx, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, flux, ok, face_viscosity, face_conductivity, &
      face_diffusion)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: left_state(:), right_state(:)
    real(dp), intent(in) :: left_temperature_guess, right_temperature_guess
    real(dp), intent(in) :: dx
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux(:)
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: face_viscosity, face_conductivity
    real(dp), intent(out), optional :: face_diffusion(:)

    real(dp), allocatable :: qleft(:), qright(:), yleft(:), yright(:)
    real(dp), allocatable :: yface(:), xleft(:), xright(:), xface(:)
    real(dp), allocatable :: diffusion(:), raw_species_flux(:)
    real(dp), allocatable :: species_flux(:), hleft(:), hright(:)
    real(dp) :: tleft, tright, cleft, cright, tface, pface, rhoface
    real(dp) :: viscosity, conductivity, dudx, dvdx, dwdx, dtdx
    real(dp) :: tau_xx, tau_xy, tau_xz, dlnpdx, correction_flux
    real(dp) :: cp, cv, internal_energy, entropy
    logical :: local_ok
    integer :: nspecies, k

    flux = 0.0_dp
    ok = .false.
    if (present(face_viscosity)) face_viscosity = 0.0_dp
    if (present(face_conductivity)) face_conductivity = 0.0_dp
    if (present(face_diffusion)) face_diffusion = 0.0_dp
    nspecies = size(species)
    if (size(transport) /= nspecies .or. &
        size(left_state) /= reactive_nvar(nspecies) .or. &
        size(right_state) /= reactive_nvar(nspecies) .or. &
        size(flux) /= reactive_nvar(nspecies) .or. dx <= 0.0_dp) return
    if (present(face_diffusion)) then
      if (size(face_diffusion) /= nspecies) return
    end if
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    allocate(qleft(reactive_nprim(nspecies)), qright(reactive_nprim(nspecies)))
    allocate(yleft(nspecies), yright(nspecies), yface(nspecies))
    allocate(xleft(nspecies), xright(nspecies), xface(nspecies))
    allocate(diffusion(nspecies), raw_species_flux(nspecies))
    allocate(species_flux(nspecies), hleft(nspecies), hright(nspecies))

    call reactive_conserved_to_primitive( &
      species, left_state, left_temperature_guess, qleft, tleft, cleft, &
      local_ok)
    if (.not. local_ok) return
    call reactive_conserved_to_primitive( &
      species, right_state, right_temperature_guess, qright, tright, cright, &
      local_ok)
    if (.not. local_ok) return
    do k = 1, nspecies
      yleft(k) = qleft(reactive_mass_fraction_component(k))
      yright(k) = qright(reactive_mass_fraction_component(k))
    end do
    yface = 0.5_dp * (yleft + yright)
    if (sum(yface) <= 0.0_dp) return
    yface = max(0.0_dp, yface)
    yface = yface / sum(yface)
    tface = 0.5_dp * (tleft + tright)
    pface = 0.5_dp * (qleft(5) + qright(5))
    rhoface = 0.5_dp * (qleft(1) + qright(1))
    call mixture_transport_coefficients( &
      species, transport, yface, tface, pface, viscosity, conductivity, &
      diffusion, local_ok)
    if (.not. local_ok) return
    if (present(face_viscosity)) face_viscosity = viscosity
    if (present(face_conductivity)) face_conductivity = conductivity
    if (present(face_diffusion)) face_diffusion = diffusion

    dudx = (qright(2) - qleft(2)) / dx
    dvdx = (qright(3) - qleft(3)) / dx
    dwdx = (qright(4) - qleft(4)) / dx
    dtdx = (tright - tleft) / dx
    tau_xx = 0.0_dp
    tau_xy = 0.0_dp
    tau_xz = 0.0_dp
    if (viscosity_enabled) then
      tau_xx = (4.0_dp / 3.0_dp) * viscosity * dudx
      tau_xy = viscosity * dvdx
      tau_xz = viscosity * dwdx
      flux(imx) = -tau_xx
      flux(imy) = -tau_xy
      flux(imz) = -tau_xz
      flux(iet) = -(tau_xx * 0.5_dp * (qleft(2) + qright(2)) + &
        tau_xy * 0.5_dp * (qleft(3) + qright(3)) + &
        tau_xz * 0.5_dp * (qleft(4) + qright(4)))
    end if
    if (thermal_conduction_enabled) flux(iet) = flux(iet) - conductivity * dtdx

    if (species_diffusion_enabled) then
      call mole_fractions_from_mass_fractions( &
        species, yleft, xleft, local_ok)
      if (.not. local_ok) return
      call mole_fractions_from_mass_fractions( &
        species, yright, xright, local_ok)
      if (.not. local_ok) return
      xface = 0.5_dp * (xleft + xright)
      xface = xface / sum(xface)
      dlnpdx = 0.0_dp
      if (barodiffusion_enabled) then
        dlnpdx = (qright(5) - qleft(5)) / (dx * pface)
      end if
      do k = 1, nspecies
        raw_species_flux(k) = -rhoface * diffusion(k) * &
          ((xright(k) - xleft(k)) / dx + &
            (xface(k) - yface(k)) * dlnpdx)
      end do
      correction_flux = sum(raw_species_flux)
      species_flux = raw_species_flux - yface * correction_flux
      if (nspecies > 1) then
        species_flux(nspecies) = -sum(species_flux(1:nspecies - 1))
      else
        species_flux(1) = 0.0_dp
      end if
      do k = 1, nspecies
        flux(reactive_species_component(k)) = species_flux(k)
        call nasa7_mass_properties( &
          species(k), tleft, cp, cv, hleft(k), internal_energy, entropy, &
          local_ok)
        if (.not. local_ok) return
        call nasa7_mass_properties( &
          species(k), tright, cp, cv, hright(k), internal_energy, entropy, &
          local_ok)
        if (.not. local_ok) return
      end do
      flux(iet) = flux(iet) + sum(0.5_dp * (hleft + hright) * species_flux)
    end if
    flux(irho) = 0.0_dp
    ok = all(ieee_is_finite(flux))
  end subroutine reactive_diffusive_flux_x

  subroutine reactive_transport_timestep( &
      species, transport, state, temperature, nx, dx, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, maximum_diffusivity, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt, maximum_diffusivity
    logical, intent(out) :: ok

    real(dp), allocatable :: q(:), y(:), diffusion(:)
    real(dp) :: local_t, sound_speed, viscosity, conductivity
    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, energy, entropy, rho, candidate
    logical :: local_ok
    integer :: i, k, nspecies

    dt = 0.0_dp
    maximum_diffusivity = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (size(transport) /= nspecies .or. nx < 1 .or. dx <= 0.0_dp .or. &
        transport_cfl <= 0.0_dp .or. transport_cfl > 0.5_dp) return
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) then
      dt = huge(1.0_dp)
      ok = .true.
      return
    end if
    allocate(q(reactive_nprim(nspecies)), y(nspecies), diffusion(nspecies))
    do i = 1, nx
      call reactive_conserved_to_primitive( &
        species, state(:, i), temperature(i), q, local_t, sound_speed, local_ok)
      if (.not. local_ok) return
      rho = q(1)
      do k = 1, nspecies
        y(k) = q(reactive_mass_fraction_component(k))
      end do
      call mixture_transport_coefficients( &
        species, transport, y, local_t, q(5), viscosity, conductivity, &
        diffusion, local_ok)
      if (.not. local_ok) return
      call mixture_mass_properties( &
        species, y, local_t, molecular_weight, gas_constant, cp, cv, gamma, &
        enthalpy, energy, entropy, local_ok)
      if (.not. local_ok .or. rho <= 0.0_dp .or. cv <= 0.0_dp) return
      if (viscosity_enabled) then
        candidate = (4.0_dp / 3.0_dp) * viscosity / rho
        maximum_diffusivity = max(maximum_diffusivity, candidate)
      end if
      if (thermal_conduction_enabled) then
        candidate = conductivity / (rho * cv)
        maximum_diffusivity = max(maximum_diffusivity, candidate)
      end if
      if (species_diffusion_enabled) then
        maximum_diffusivity = max(maximum_diffusivity, maxval(diffusion))
      end if
    end do
    if (maximum_diffusivity <= 0.0_dp) then
      dt = huge(1.0_dp)
    else
      dt = transport_cfl * dx * dx / maximum_diffusivity
    end if
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp .and. &
      ieee_is_finite(maximum_diffusivity)
  end subroutine reactive_transport_timestep

  subroutine reactive_transport_euler_update( &
      species, transport, input_state, input_temperature, nx, dx, dt, &
      boundary, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, output_state, &
      output_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: input_state(:, 0:), input_temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt
    character(len=*), intent(in) :: boundary
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: output_state(:, 0:), output_temperature(0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: work_state(:, :), work_temperature(:)
    real(dp), allocatable :: flux(:, :), q(:)
    real(dp) :: local_t, sound_speed
    logical :: local_ok
    integer :: i, nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    allocate(work_state(nvar, 0:nx + 1), work_temperature(0:nx + 1))
    allocate(flux(nvar, 0:nx), q(reactive_nprim(size(species))))
    work_state = input_state
    work_temperature = input_temperature
    call fill_ghosts(work_state, work_temperature, nx, boundary, local_ok)
    if (.not. local_ok) return
    do i = 0, nx
      call reactive_diffusive_flux_x( &
        species, transport, work_state(:, i), work_state(:, i + 1), &
        work_temperature(i), work_temperature(i + 1), dx, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, flux(:, i), &
        local_ok)
      if (.not. local_ok) return
    end do
    output_state = work_state
    output_temperature = work_temperature
    do i = 1, nx
      output_state(:, i) = work_state(:, i) - dt / dx * &
        (flux(:, i) - flux(:, i - 1))
      call reactive_conserved_to_primitive( &
        species, output_state(:, i), work_temperature(i), q, local_t, &
        sound_speed, local_ok)
      if (.not. local_ok) return
      output_temperature(i) = local_t
    end do
    call fill_ghosts(output_state, output_temperature, nx, boundary, ok)
  end subroutine reactive_transport_euler_update

  subroutine advance_reactive_transport( &
      species, transport, state, temperature, nx, dx, interval, boundary, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, interval
    character(len=*), intent(in) :: boundary
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    logical, intent(out) :: ok

    real(dp), allocatable :: initial_state(:, :), initial_temperature(:)
    real(dp), allocatable :: stage1_state(:, :), stage1_temperature(:)
    real(dp), allocatable :: euler2_state(:, :), euler2_temperature(:)
    real(dp), allocatable :: q(:)
    real(dp) :: local_t, sound_speed, temperature_guess
    logical :: local_ok
    integer :: i, nvar

    ok = .false.
    if (interval < 0.0_dp) return
    if (interval <= 0.0_dp .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      return
    end if
    nvar = reactive_nvar(size(species))
    allocate(initial_state(nvar, 0:nx + 1), initial_temperature(0:nx + 1))
    allocate(stage1_state(nvar, 0:nx + 1), stage1_temperature(0:nx + 1))
    allocate(euler2_state(nvar, 0:nx + 1), euler2_temperature(0:nx + 1))
    allocate(q(reactive_nprim(size(species))))
    initial_state = state
    initial_temperature = temperature

    call reactive_transport_euler_update( &
      species, transport, initial_state, initial_temperature, nx, dx, &
      interval, boundary, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, stage1_state, &
      stage1_temperature, local_ok)
    if (.not. local_ok) return
    call reactive_transport_euler_update( &
      species, transport, stage1_state, stage1_temperature, nx, dx, &
      interval, boundary, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, euler2_state, &
      euler2_temperature, local_ok)
    if (.not. local_ok) return

    state = initial_state
    temperature = initial_temperature
    do i = 1, nx
      state(:, i) = 0.5_dp * (initial_state(:, i) + euler2_state(:, i))
      temperature_guess = 0.5_dp * &
        (initial_temperature(i) + euler2_temperature(i))
      call reactive_conserved_to_primitive( &
        species, state(:, i), temperature_guess, q, local_t, sound_speed, &
        local_ok)
      if (.not. local_ok) return
      temperature(i) = local_t
    end do
    call fill_ghosts(state, temperature, nx, boundary, ok)
  end subroutine advance_reactive_transport

  subroutine reactive_hydro_euler_update( &
      species, input_state, input_temperature, nx, dx, dt, reconstruction, &
      limiter, boundary, riemann_solver, use_contact_steepening, &
      use_shock_flattening, output_state, output_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: input_state(:, 0:), input_temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt
    character(len=*), intent(in) :: reconstruction, limiter, boundary
    character(len=*), intent(in) :: riemann_solver
    logical, intent(in) :: use_contact_steepening, use_shock_flattening
    real(dp), intent(out) :: output_state(:, 0:), output_temperature(0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: work_state(:, :), work_temperature(:)
    real(dp), allocatable :: lf(:, :), rf(:, :), flux(:, :), q(:)
    real(dp), allocatable :: lt(:), rt(:)
    real(dp) :: local_t, c, predictor_dtdx
    logical :: local_ok
    integer :: nvar, i

    ok = .false.
    nvar = reactive_nvar(size(species))
    allocate(work_state(nvar, 0:nx + 1), work_temperature(0:nx + 1))
    work_state = input_state
    work_temperature = input_temperature
    call fill_ghosts(work_state, work_temperature, nx, boundary, local_ok)
    if (.not. local_ok) return

    allocate(lf(nvar, 0:nx), rf(nvar, 0:nx), flux(nvar, 0:nx))
    allocate(lt(0:nx), rt(0:nx), q(reactive_nprim(size(species))))
    predictor_dtdx = 0.0_dp
    if (trim(reconstruction) == "characteristic_plm" .or. &
        trim(reconstruction) == "characteristic_ppm") predictor_dtdx = dt / dx
    call reconstruct_faces( &
      species, work_state, work_temperature, nx, reconstruction, limiter, &
      boundary, predictor_dtdx, use_contact_steepening, &
      use_shock_flattening, lf, rf, lt, rt, local_ok)
    if (.not. local_ok) return
    do i = 0, nx
      call reactive_riemann_flux_x( &
        species, lf(:, i), rf(:, i), lt(i), rt(i), riemann_solver, &
        flux(:, i), local_ok)
      if (.not. local_ok) return
    end do

    output_state = work_state
    output_temperature = work_temperature
    do i = 1, nx
      output_state(:, i) = work_state(:, i) - dt / dx * &
        (flux(:, i) - flux(:, i - 1))
      call reactive_conserved_to_primitive( &
        species, output_state(:, i), work_temperature(i), q, local_t, c, &
        local_ok)
      if (.not. local_ok) return
      output_temperature(i) = local_t
    end do
    call fill_ghosts(output_state, output_temperature, nx, boundary, ok)
  end subroutine reactive_hydro_euler_update

  subroutine advance_reactive_hydro( &
      species, state, temperature, nx, dx, dt, reconstruction, limiter, &
      boundary, ok, riemann_solver, ppm_contact_steepening, &
      ppm_shock_flattening)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt
    character(len=*), intent(in) :: reconstruction, limiter, boundary
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: riemann_solver
    logical, intent(in), optional :: ppm_contact_steepening
    logical, intent(in), optional :: ppm_shock_flattening

    real(dp), allocatable :: initial_state(:, :), initial_temperature(:)
    real(dp), allocatable :: stage1_state(:, :), stage1_temperature(:)
    real(dp), allocatable :: euler2_state(:, :), euler2_temperature(:)
    real(dp), allocatable :: stage2_state(:, :), stage2_temperature(:)
    real(dp), allocatable :: euler3_state(:, :), euler3_temperature(:)
    real(dp), allocatable :: q(:)
    real(dp) :: local_t, c, temperature_guess
    logical :: local_ok, use_contact_steepening, use_shock_flattening
    integer :: nvar, i
    character(len=32) :: selected_solver

    ok = .false.
    selected_solver = "rusanov"
    if (present(riemann_solver)) selected_solver = trim(riemann_solver)
    use_contact_steepening = .false.
    use_shock_flattening = .false.
    if (present(ppm_contact_steepening)) &
      use_contact_steepening = ppm_contact_steepening
    if (present(ppm_shock_flattening)) &
      use_shock_flattening = ppm_shock_flattening
    nvar = reactive_nvar(size(species))
    allocate(initial_state(nvar, 0:nx + 1), initial_temperature(0:nx + 1))
    allocate(stage1_state(nvar, 0:nx + 1), stage1_temperature(0:nx + 1))
    initial_state = state
    initial_temperature = temperature

    if (trim(reconstruction) /= "ppm") then
      call reactive_hydro_euler_update( &
        species, initial_state, initial_temperature, nx, dx, dt, &
        reconstruction, limiter, boundary, selected_solver, &
        use_contact_steepening, use_shock_flattening, stage1_state, &
        stage1_temperature, ok)
      if (ok) then
        state = stage1_state
        temperature = stage1_temperature
      end if
      return
    end if

    allocate(euler2_state(nvar, 0:nx + 1), euler2_temperature(0:nx + 1))
    allocate(stage2_state(nvar, 0:nx + 1), stage2_temperature(0:nx + 1))
    allocate(euler3_state(nvar, 0:nx + 1), euler3_temperature(0:nx + 1))
    allocate(q(reactive_nprim(size(species))))

    ! SSPRK3 stage 1: U1 = U^n + dt L(U^n).
    call reactive_hydro_euler_update( &
      species, initial_state, initial_temperature, nx, dx, dt, "ppm", &
      limiter, boundary, selected_solver, use_contact_steepening, &
      use_shock_flattening, stage1_state, stage1_temperature, local_ok)
    if (.not. local_ok) return

    ! SSPRK3 stage 2: U2 = 3/4 U^n + 1/4 (U1 + dt L(U1)).
    call reactive_hydro_euler_update( &
      species, stage1_state, stage1_temperature, nx, dx, dt, "ppm", limiter, &
      boundary, selected_solver, use_contact_steepening, &
      use_shock_flattening, euler2_state, euler2_temperature, local_ok)
    if (.not. local_ok) return
    stage2_state = initial_state
    stage2_temperature = initial_temperature
    do i = 1, nx
      stage2_state(:, i) = 0.75_dp * initial_state(:, i) + &
        0.25_dp * euler2_state(:, i)
      temperature_guess = 0.75_dp * initial_temperature(i) + &
        0.25_dp * euler2_temperature(i)
      call reactive_conserved_to_primitive( &
        species, stage2_state(:, i), temperature_guess, q, local_t, c, &
        local_ok)
      if (.not. local_ok) return
      stage2_temperature(i) = local_t
    end do
    call fill_ghosts(stage2_state, stage2_temperature, nx, boundary, local_ok)
    if (.not. local_ok) return

    ! SSPRK3 stage 3: U^{n+1} = 1/3 U^n + 2/3 (U2 + dt L(U2)).
    call reactive_hydro_euler_update( &
      species, stage2_state, stage2_temperature, nx, dx, dt, "ppm", limiter, &
      boundary, selected_solver, use_contact_steepening, &
      use_shock_flattening, euler3_state, euler3_temperature, local_ok)
    if (.not. local_ok) return
    state = initial_state
    temperature = initial_temperature
    do i = 1, nx
      state(:, i) = initial_state(:, i) / 3.0_dp + &
        2.0_dp * euler3_state(:, i) / 3.0_dp
      temperature_guess = initial_temperature(i) / 3.0_dp + &
        2.0_dp * euler3_temperature(i) / 3.0_dp
      call reactive_conserved_to_primitive( &
        species, state(:, i), temperature_guess, q, local_t, c, local_ok)
      if (.not. local_ok) return
      temperature(i) = local_t
    end do
    call fill_ghosts(state, temperature, nx, boundary, ok)
  end subroutine advance_reactive_hydro

  subroutine advance_reactive_chemistry( &
      species, reactions, state, temperature, nx, interval, rtol, atol, &
      boundary, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: interval, rtol, atol
    character(len=*), intent(in) :: boundary
    logical, intent(out) :: ok
    real(dp), allocatable :: y(:)
    real(dp) :: rho, kinetic_density, target_energy
    real(dp) :: elapsed, request, accepted, next_step, tolerance
    logical :: local_ok
    integer :: i, k, substeps, newton_iterations, rejected_attempts

    ok = .false.
    if (interval <= 0.0_dp) then
      ok = interval >= 0.0_dp
      return
    end if
    allocate(y(size(species)))
    tolerance = 50.0_dp * epsilon(1.0_dp) * max(1.0_dp, interval)
    do i = 1, nx
      rho = state(irho, i)
      call mass_fractions_from_state(state(:, i), size(species), y, local_ok)
      if (.not. local_ok) return
      kinetic_density = 0.5_dp * &
        (state(imx, i)**2 + state(imy, i)**2 + state(imz, i)**2) / rho
      target_energy = (state(iet, i) - kinetic_density) / rho
      elapsed = 0.0_dp
      request = interval
      substeps = 0
      do while (elapsed < interval - tolerance)
        if (substeps >= max_chemistry_substeps) return
        request = min(request, interval - elapsed)
        if (size(species) == 10) then
          call advance_constant_volume_implicit_adaptive( &
            species, reactions, rho, target_energy, request, rtol, atol, y, &
            temperature(i), accepted, next_step, newton_iterations, &
            rejected_attempts, local_ok)
        else
          call advance_constant_volume_adaptive( &
            species, reactions, rho, target_energy, request, rtol, atol, y, &
            temperature(i), accepted, next_step, local_ok)
        end if
        if (.not. local_ok .or. accepted <= 0.0_dp) return
        elapsed = elapsed + accepted
        request = min(next_step, interval - elapsed)
        substeps = substeps + 1
      end do
      do k = 1, size(species)
        state(reactive_species_component(k), i) = rho * y(k)
      end do
    end do
    call fill_ghosts(state, temperature, nx, boundary, ok)
  end subroutine advance_reactive_chemistry

  subroutine advance_reactive_strang( &
      species, reactions, state, temperature, nx, dx, dt, reconstruction, &
      limiter, boundary, chemistry_enabled, rtol, atol, ok, riemann_solver, &
      ppm_contact_steepening, ppm_shock_flattening, transport, &
      transport_enabled, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt, rtol, atol
    character(len=*), intent(in) :: reconstruction, limiter, boundary
    logical, intent(in) :: chemistry_enabled
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: riemann_solver
    logical, intent(in), optional :: ppm_contact_steepening
    logical, intent(in), optional :: ppm_shock_flattening
    type(gas_transport_species), intent(in), optional :: transport(:)
    logical, intent(in), optional :: transport_enabled, viscosity_enabled
    logical, intent(in), optional :: thermal_conduction_enabled
    logical, intent(in), optional :: species_diffusion_enabled
    logical, intent(in), optional :: barodiffusion_enabled
    logical :: local_ok, use_contact_steepening, use_shock_flattening
    logical :: use_transport, use_viscosity, use_conduction
    logical :: use_species_diffusion, use_barodiffusion
    character(len=32) :: selected_solver

    ok = .false.
    selected_solver = "rusanov"
    if (present(riemann_solver)) selected_solver = trim(riemann_solver)
    use_contact_steepening = .false.
    use_shock_flattening = .false.
    if (present(ppm_contact_steepening)) &
      use_contact_steepening = ppm_contact_steepening
    if (present(ppm_shock_flattening)) &
      use_shock_flattening = ppm_shock_flattening
    use_transport = .false.
    use_viscosity = .true.
    use_conduction = .true.
    use_species_diffusion = .true.
    use_barodiffusion = .true.
    if (present(transport_enabled)) use_transport = transport_enabled
    if (present(viscosity_enabled)) use_viscosity = viscosity_enabled
    if (present(thermal_conduction_enabled)) &
      use_conduction = thermal_conduction_enabled
    if (present(species_diffusion_enabled)) &
      use_species_diffusion = species_diffusion_enabled
    if (present(barodiffusion_enabled)) &
      use_barodiffusion = barodiffusion_enabled
    if (use_transport .and. .not. present(transport)) return
    if (chemistry_enabled) then
      call advance_reactive_chemistry(species, reactions, state, temperature, &
        nx, 0.5_dp * dt, rtol, atol, boundary, local_ok)
      if (.not. local_ok) return
    end if
    if (use_transport) then
      call advance_reactive_transport( &
        species, transport, state, temperature, nx, dx, 0.5_dp * dt, &
        boundary, use_viscosity, use_conduction, use_species_diffusion, &
        use_barodiffusion, local_ok)
      if (.not. local_ok) return
    end if
    call advance_reactive_hydro(species, state, temperature, nx, dx, dt, &
      reconstruction, limiter, boundary, local_ok, selected_solver, &
      use_contact_steepening, use_shock_flattening)
    if (.not. local_ok) return
    if (use_transport) then
      call advance_reactive_transport( &
        species, transport, state, temperature, nx, dx, 0.5_dp * dt, &
        boundary, use_viscosity, use_conduction, use_species_diffusion, &
        use_barodiffusion, local_ok)
      if (.not. local_ok) return
    end if
    if (chemistry_enabled) then
      call advance_reactive_chemistry(species, reactions, state, temperature, &
        nx, 0.5_dp * dt, rtol, atol, boundary, local_ok)
      if (.not. local_ok) return
    end if
    ok = .true.
  end subroutine advance_reactive_strang

  subroutine initialize_reactive_1d(species, config, state, temperature, dx, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: state(:, :), temperature(:)
    real(dp), intent(out) :: dx
    logical, intent(out) :: ok
    real(dp), allocatable :: q(:), base_xmol(:), local_xmol(:)
    real(dp), allocatable :: base_y(:), local_y(:)
    real(dp) :: x, rho, base_rho, local_temperature, c, gaussian
    real(dp) :: phase
    logical :: local_ok
    integer :: i, k, nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    allocate(state(nvar, 0:config%nx + 1), temperature(0:config%nx + 1))
    allocate(q(reactive_nprim(size(species))), base_xmol(size(species)), &
      local_xmol(size(species)))
    allocate(base_y(size(species)), local_y(size(species)))
    dx = (config%x_upper - config%x_lower) / real(config%nx, dp)
    call reactive_1d_mole_fractions(config, size(species), base_xmol, local_ok)
    if (.not. local_ok) return
    call mass_fractions_from_mole_fractions(species, base_xmol, base_y, local_ok)
    if (.not. local_ok) return
    base_rho = mixture_density(species, base_y, config%initial_pressure, &
      config%initial_temperature, local_ok)
    if (.not. local_ok) return

    do i = 1, config%nx
      x = config%x_lower + (real(i, dp) - 0.5_dp) * dx
      local_y = base_y
      select case (trim(config%problem))
      case ("entropy_wave")
        rho = base_rho * (1.0_dp + config%density_wave_amplitude * &
          sin(2.0_dp * pi * (x - config%x_lower) / &
            (config%x_upper - config%x_lower)))
      case ("composition_wave")
        phase = sin(2.0_dp * pi * (x - config%x_lower) / &
          (config%x_upper - config%x_lower))
        local_xmol = base_xmol
        local_xmol(1) = local_xmol(1) + &
          config%composition_wave_amplitude * phase
        local_xmol(size(species)) = local_xmol(size(species)) - &
          config%composition_wave_amplitude * phase
        call mass_fractions_from_mole_fractions( &
          species, local_xmol, local_y, local_ok)
        if (.not. local_ok) return
        rho = mixture_density(species, local_y, config%initial_pressure, &
          config%initial_temperature, local_ok)
        if (.not. local_ok) return
      case ("uniform_reactor")
        rho = base_rho
      case ("reactive_hotspot")
        gaussian = exp(-((x - config%hotspot_center) / config%hotspot_width)**2)
        local_temperature = config%initial_temperature + &
          config%hotspot_temperature_rise * gaussian
        rho = mixture_density(species, local_y, config%initial_pressure, &
          local_temperature, local_ok)
        if (.not. local_ok) return
      case default
        return
      end select
      q(1:5) = [rho, config%initial_velocity, 0.0_dp, 0.0_dp, &
        config%initial_pressure]
      do k = 1, size(species)
        q(reactive_mass_fraction_component(k)) = local_y(k)
      end do
      call reactive_primitive_to_conserved( &
        species, q, state(:, i), temperature(i), c, local_ok)
      if (.not. local_ok) return
    end do
    call fill_ghosts(state, temperature, config%nx, &
      config%boundary_condition, ok)
  end subroutine initialize_reactive_1d

  subroutine reactive_integrals(state, nx, dx, integrals)
    real(dp), intent(in) :: state(:, 0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx
    real(dp), intent(out) :: integrals(5)
    integrals(1) = dx * sum(state(irho, 1:nx))
    integrals(2) = dx * sum(state(imx, 1:nx))
    integrals(3) = dx * sum(state(imy, 1:nx))
    integrals(4) = dx * sum(state(imz, 1:nx))
    integrals(5) = dx * sum(state(iet, 1:nx))
  end subroutine reactive_integrals

  subroutine simulate_reactive_1d( &
      species, reactions, config, state, temperature, dx, time, steps, &
      initial_integrals, final_integrals, ok, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: state(:, :), temperature(:)
    real(dp), intent(out) :: dx, time
    integer, intent(out) :: steps
    real(dp), intent(out) :: initial_integrals(5), final_integrals(5)
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)
    real(dp) :: dt, hydro_dt, transport_dt, maximum_diffusivity, tolerance
    logical :: local_ok

    time = 0.0_dp
    steps = 0
    call initialize_reactive_1d(species, config, state, temperature, dx, local_ok)
    if (.not. local_ok) then
      ok = .false.; return
    end if
    call reactive_integrals(state, config%nx, dx, initial_integrals)
    tolerance = 50.0_dp * epsilon(1.0_dp) * max(1.0_dp, config%final_time)
    do while (time < config%final_time - tolerance)
      if (steps >= config%maximum_steps) then
        ok = .false.; return
      end if
      call reactive_cfl_timestep(species, state, temperature, config%nx, dx, &
        config%cfl, hydro_dt, local_ok)
      if (.not. local_ok) then
        ok = .false.; return
      end if
      dt = hydro_dt
      if (config%transport_enabled) then
        if (.not. present(transport)) then
          ok = .false.; return
        end if
        call reactive_transport_timestep( &
          species, transport, state, temperature, config%nx, dx, &
          config%transport_cfl, config%viscosity_enabled, &
          config%thermal_conduction_enabled, config%species_diffusion_enabled, &
          transport_dt, maximum_diffusivity, local_ok)
        if (.not. local_ok) then
          ok = .false.; return
        end if
        dt = min(dt, transport_dt)
      end if
      dt = min(dt, config%final_time - time)
      if (config%transport_enabled) then
        call advance_reactive_strang( &
          species, reactions, state, temperature, config%nx, dx, dt, &
          config%reconstruction, config%limiter, config%boundary_condition, &
          config%chemistry_enabled, config%chemistry_relative_tolerance, &
          config%chemistry_absolute_tolerance, local_ok, config%riemann_solver, &
          config%ppm_contact_steepening, config%ppm_shock_flattening, &
          transport, .true., config%viscosity_enabled, &
          config%thermal_conduction_enabled, config%species_diffusion_enabled, &
          config%barodiffusion_enabled)
      else
        call advance_reactive_strang( &
          species, reactions, state, temperature, config%nx, dx, dt, &
          config%reconstruction, config%limiter, config%boundary_condition, &
          config%chemistry_enabled, config%chemistry_relative_tolerance, &
          config%chemistry_absolute_tolerance, local_ok, config%riemann_solver, &
          config%ppm_contact_steepening, config%ppm_shock_flattening)
      end if
      if (.not. local_ok) then
        ok = .false.; return
      end if
      time = time + dt
      steps = steps + 1
    end do
    time = config%final_time
    call reactive_integrals(state, config%nx, dx, final_integrals)
    ok = .true.
  end subroutine simulate_reactive_1d

  subroutine write_reactive_1d_csv(path, species, config, state, temperature, &
      dx, time, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: state(:, 0:), temperature(0:), dx, time
    logical, intent(out) :: ok
    real(dp), allocatable :: q(:)
    real(dp) :: local_t, c, x
    logical :: local_ok
    integer :: unit, status, i, k

    ok = .false.
    allocate(q(reactive_nprim(size(species))))
    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=status)
    if (status /= 0) return
    write(unit, '(a)', advance='no') &
      "time,x,rho,u,v,w,pressure,temperature,rhoE"
    do k = 1, size(species)
      write(unit, '(a)', advance='no') ",Y_" // trim(species(k)%name)
    end do
    write(unit, '(a)') ""
    do i = 1, config%nx
      call reactive_conserved_to_primitive(species, state(:, i), temperature(i), &
        q, local_t, c, local_ok)
      if (.not. local_ok) then
        close(unit); return
      end if
      x = config%x_lower + (real(i, dp) - 0.5_dp) * dx
      write(unit, '(*(es25.16e3,:,","))') time, x, state(irho, i), &
        q(2), q(3), q(4), q(5), local_t, state(iet, i), &
        (q(reactive_mass_fraction_component(k)), k = 1, size(species))
    end do
    close(unit)
    ok = .true.
  end subroutine write_reactive_1d_csv

  pure real(dp) function reactive_entropy_wave_density(x, time, config, base_rho) &
      result(rho)
    real(dp), intent(in) :: x, time, base_rho
    type(reactive_1d_config), intent(in) :: config
    real(dp) :: length, shifted
    length = config%x_upper - config%x_lower
    shifted = config%x_lower + modulo(x - config%x_lower - &
      config%initial_velocity * time, length)
    rho = base_rho * (1.0_dp + config%density_wave_amplitude * &
      sin(2.0_dp * pi * (shifted - config%x_lower) / length))
  end function reactive_entropy_wave_density

  subroutine reactive_composition_wave_exact( &
      species, x, time, config, density, mass_fractions, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: x, time
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(out) :: density, mass_fractions(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: mole_fractions(:)
    real(dp) :: length, shifted, phase

    density = 0.0_dp
    mass_fractions = 0.0_dp
    ok = .false.
    if (size(species) < 1 .or. size(mass_fractions) /= size(species)) return
    allocate(mole_fractions(size(species)))
    length = config%x_upper - config%x_lower
    if (length <= 0.0_dp) return
    shifted = config%x_lower + modulo(x - config%x_lower - &
      config%initial_velocity * time, length)
    phase = sin(2.0_dp * pi * (shifted - config%x_lower) / length)
    call reactive_1d_mole_fractions(config, size(species), mole_fractions, ok)
    if (.not. ok) return
    mole_fractions(1) = mole_fractions(1) + &
      config%composition_wave_amplitude * phase
    mole_fractions(size(species)) = mole_fractions(size(species)) - &
      config%composition_wave_amplitude * phase
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, ok)
    if (.not. ok) return
    density = mixture_density(species, mass_fractions, &
      config%initial_pressure, config%initial_temperature, ok)
  end subroutine reactive_composition_wave_exact

end module reactive_1d_mod
