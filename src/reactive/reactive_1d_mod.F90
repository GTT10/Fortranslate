module reactive_1d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use state_indices_mod, only: irho, imx, imy, imz, iet, ncons
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    valid_mixture_composition, mixture_mass_properties, &
    mixture_specific_gas_constant, mixture_pressure, mixture_density, &
    mixture_sound_speed, temperature_from_internal_energy, &
    mass_fractions_from_mole_fractions
  use elementary_kinetics_mod, only: elementary_reaction
  use constant_volume_reactor_mod, only: advance_constant_volume_adaptive
  use slope_limiter_mod, only: limited_slope
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
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
  real(dp), parameter :: pi = acos(-1.0_dp)

  public :: reactive_nvar, reactive_nprim
  public :: reactive_species_component, reactive_mass_fraction_component
  public :: reactive_primitive_to_conserved
  public :: reactive_conserved_to_primitive
  public :: reactive_physical_flux_x
  public :: reactive_rusanov_flux_x
  public :: reactive_hllc_flux_x
  public :: reactive_pelec_flux_x
  public :: compute_reactive_riemann_flux_x
  public :: reactive_difference_to_characteristics
  public :: reactive_characteristics_to_difference
  public :: trace_reactive_characteristics
  public :: initialize_reactive_1d
  public :: reactive_cfl_timestep
  public :: advance_reactive_hydro
  public :: advance_reactive_chemistry
  public :: advance_reactive_strang
  public :: simulate_reactive_1d
  public :: write_reactive_1d_csv
  public :: reactive_integrals
  public :: reactive_entropy_wave_density

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

    real(dp), allocatable :: fl(:), fr(:), ql(:), qr(:), star(:), qstar(:)
    real(dp) :: tl, tr, cl, cr, tstar, cstar
    real(dp) :: rho_l, rho_r, u_l, u_r, p_l, p_r
    real(dp) :: s_l, s_r, s_m, denominator, p_star_l, p_star_r, p_star
    logical :: left_ok, right_ok, star_ok
    integer :: nspecies, nvar

    nspecies = size(species)
    nvar = reactive_nvar(nspecies)
    flux = 0.0_dp
    ok = .false.
    if (size(left_state) /= nvar .or. size(right_state) /= nvar .or. &
        size(flux) /= nvar) return

    allocate(fl(nvar), fr(nvar), star(nvar))
    allocate(ql(reactive_nprim(nspecies)), qr(reactive_nprim(nspecies)))
    allocate(qstar(reactive_nprim(nspecies)))
    call reactive_physical_flux_x( &
      species, left_state, left_temperature_guess, fl, tl, cl, ql, left_ok)
    call reactive_physical_flux_x( &
      species, right_state, right_temperature_guess, fr, tr, cr, qr, right_ok)
    if (.not. (left_ok .and. right_ok)) return

    rho_l = ql(1)
    rho_r = qr(1)
    u_l = ql(2)
    u_r = qr(2)
    p_l = ql(5)
    p_r = qr(5)
    s_l = min(u_l - cl, u_r - cr)
    s_r = max(u_l + cl, u_r + cr)
    if (.not. (ieee_is_finite(s_l) .and. ieee_is_finite(s_r))) return
    if (s_l >= 0.0_dp) then
      flux = fl
      ok = .true.
      return
    else if (s_r <= 0.0_dp) then
      flux = fr
      ok = .true.
      return
    end if

    denominator = rho_l * (s_l - u_l) - rho_r * (s_r - u_r)
    if (abs(denominator) <= 100.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, rho_l * abs(s_l - u_l), rho_r * abs(s_r - u_r))) return
    s_m = (p_r - p_l + rho_l * u_l * (s_l - u_l) - &
      rho_r * u_r * (s_r - u_r)) / denominator
    if (.not. ieee_is_finite(s_m)) return
    if (s_m <= s_l .or. s_m >= s_r) return

    p_star_l = p_l + rho_l * (s_l - u_l) * (s_m - u_l)
    p_star_r = p_r + rho_r * (s_r - u_r) * (s_m - u_r)
    p_star = 0.5_dp * (p_star_l + p_star_r)
    if (.not. ieee_is_finite(p_star) .or. p_star <= pressure_floor) return

    if (s_m >= 0.0_dp) then
      call build_reactive_hllc_star( &
        left_state, ql, s_l, s_m, p_star, nspecies, star, star_ok)
      if (.not. star_ok) return
      call reactive_conserved_to_primitive( &
        species, star, tl, qstar, tstar, cstar, star_ok)
      if (.not. star_ok) return
      flux = fl + s_l * (star - left_state)
    else
      call build_reactive_hllc_star( &
        right_state, qr, s_r, s_m, p_star, nspecies, star, star_ok)
      if (.not. star_ok) return
      call reactive_conserved_to_primitive( &
        species, star, tr, qstar, tstar, cstar, star_ok)
      if (.not. star_ok) return
      flux = fr + s_r * (star - right_state)
    end if
    flux(reactive_species_component(nspecies)) = flux(irho) - &
      sum(flux(reactive_species_component(1): &
        reactive_species_component(nspecies - 1)))
    ok = all(ieee_is_finite(flux))
  end subroutine reactive_hllc_flux_x

  subroutine reactive_pelec_flux_x( &
      species, left_state, right_state, left_temperature_guess, &
      right_temperature_guess, flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: left_state(:), right_state(:)
    real(dp), intent(in) :: left_temperature_guess, right_temperature_guess
    real(dp), intent(out) :: flux(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: ql(:), qr(:), qorigin(:), qstar(:), qinterface(:)
    real(dp), allocatable :: origin_state(:), star_state(:), interface_state(:)
    real(dp), allocatable :: checked_primitive(:)
    real(dp) :: tl, tr, cl, cr, torigin, tstar, tinterface
    real(dp) :: corigin, cstar, cinterface
    real(dp) :: impedance_l, impedance_r, pressure_star, velocity_star
    real(dp) :: density_increment, wave_out, wave_in, shock_speed
    real(dp) :: speed_difference, interpolation_fraction
    real(dp) :: stationary_threshold, regularization, average_sound_speed
    real(dp) :: dummy_temperature, dummy_sound
    logical :: left_ok, right_ok, local_ok, stationary_interface
    integer :: nspecies, nvar, nprimitive, k, component

    nspecies = size(species)
    nvar = reactive_nvar(nspecies)
    nprimitive = reactive_nprim(nspecies)
    flux = 0.0_dp
    ok = .false.
    if (size(left_state) /= nvar .or. size(right_state) /= nvar .or. &
        size(flux) /= nvar) return

    allocate(ql(nprimitive), qr(nprimitive), qorigin(nprimitive))
    allocate(qstar(nprimitive), qinterface(nprimitive), checked_primitive(nprimitive))
    allocate(origin_state(nvar), star_state(nvar), interface_state(nvar))
    call reactive_conserved_to_primitive( &
      species, left_state, left_temperature_guess, ql, tl, cl, left_ok)
    call reactive_conserved_to_primitive( &
      species, right_state, right_temperature_guess, qr, tr, cr, right_ok)
    if (.not. (left_ok .and. right_ok)) return

    impedance_l = max(tiny(1.0_dp), ql(1) * cl)
    impedance_r = max(tiny(1.0_dp), qr(1) * cr)
    pressure_star = max(pressure_floor, &
      ((impedance_r * ql(5) + impedance_l * qr(5)) + &
       impedance_l * impedance_r * (ql(2) - qr(2))) / &
      (impedance_l + impedance_r))
    velocity_star = &
      ((impedance_l * ql(2) + impedance_r * qr(2)) + &
       (ql(5) - qr(5))) / (impedance_l + impedance_r)
    if (.not. (ieee_is_finite(pressure_star) .and. &
        ieee_is_finite(velocity_star))) return

    stationary_threshold = sqrt(epsilon(1.0_dp)) * &
      max(1.0_dp, abs(ql(2)), abs(qr(2)), cl, cr)
    stationary_interface = abs(velocity_star) <= stationary_threshold
    if (stationary_interface) then
      velocity_star = 0.0_dp
      qorigin(1:5) = 0.5_dp * (ql(1:5) + qr(1:5))
      do k = 1, nspecies
        component = reactive_mass_fraction_component(k)
        qorigin(component) = 0.5_dp * (ql(component) + qr(component))
      end do
      call normalize_reactive_mass_fractions(qorigin, nspecies, local_ok)
      if (.not. local_ok) return
    else if (velocity_star > 0.0_dp) then
      qorigin = ql
    else
      qorigin = qr
    end if

    call reactive_primitive_to_conserved( &
      species, qorigin, origin_state, torigin, corigin, local_ok)
    if (.not. local_ok) return
    density_increment = (pressure_star - qorigin(5)) / (corigin * corigin)
    qstar = qorigin
    qstar(1) = qorigin(1) + density_increment
    qstar(2) = velocity_star
    qstar(5) = pressure_star
    if (qstar(1) <= density_floor) return
    call reactive_primitive_to_conserved( &
      species, qstar, star_state, tstar, cstar, local_ok)
    if (.not. local_ok) return

    wave_out = corigin - sign(1.0_dp, velocity_star) * qorigin(2)
    wave_in = cstar - sign(1.0_dp, velocity_star) * velocity_star
    shock_speed = 0.5_dp * (wave_in + wave_out)
    if (pressure_star >= qorigin(5)) then
      wave_out = shock_speed
      wave_in = shock_speed
    end if

    average_sound_speed = 0.5_dp * (cl + cr)
    regularization = sqrt(epsilon(1.0_dp)) * max(1.0_dp, average_sound_speed)
    if (abs(wave_out - wave_in) < regularization) then
      speed_difference = sign(regularization, wave_out - wave_in)
    else
      speed_difference = wave_out - wave_in
    end if
    interpolation_fraction = max(0.0_dp, min(1.0_dp, &
      0.5_dp * (1.0_dp + (wave_out + wave_in) / speed_difference)))
    qinterface = interpolation_fraction * qstar + &
      (1.0_dp - interpolation_fraction) * qorigin
    if (wave_out < 0.0_dp) qinterface = qorigin
    if (wave_in >= 0.0_dp) qinterface = qstar
    call normalize_reactive_mass_fractions(qinterface, nspecies, local_ok)
    if (.not. local_ok) return
    if (qinterface(1) <= density_floor .or. &
        qinterface(5) <= pressure_floor) return

    call reactive_primitive_to_conserved( &
      species, qinterface, interface_state, tinterface, cinterface, local_ok)
    if (.not. local_ok) return
    call reactive_physical_flux_x( &
      species, interface_state, tinterface, flux, dummy_temperature, &
      dummy_sound, checked_primitive, local_ok)
    if (.not. local_ok) return
    ok = all(ieee_is_finite(flux))
  end subroutine reactive_pelec_flux_x

  pure subroutine normalize_reactive_mass_fractions(q, nspecies, ok)
    real(dp), intent(inout) :: q(:)
    integer, intent(in) :: nspecies
    logical, intent(out) :: ok
    real(dp) :: total
    integer :: k, component

    total = 0.0_dp
    ok = .false.
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      if (.not. ieee_is_finite(q(component)) .or. q(component) < -1.0e-12_dp) return
      q(component) = max(0.0_dp, q(component))
      total = total + q(component)
    end do
    if (total <= tiny(1.0_dp)) return
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      q(component) = q(component) / total
    end do
    ok = .true.
  end subroutine normalize_reactive_mass_fractions

  subroutine build_reactive_hllc_star( &
      state, primitive, wave_speed, contact_speed, star_pressure, &
      nspecies, star, ok)
    real(dp), intent(in) :: state(:), primitive(:)
    real(dp), intent(in) :: wave_speed, contact_speed, star_pressure
    integer, intent(in) :: nspecies
    real(dp), intent(out) :: star(:)
    logical, intent(out) :: ok

    real(dp) :: rho, u, v, w, pressure, denominator, factor, rho_star
    integer :: k

    star = 0.0_dp
    ok = .false.
    rho = primitive(1)
    u = primitive(2)
    v = primitive(3)
    w = primitive(4)
    pressure = primitive(5)
    denominator = wave_speed - contact_speed
    if (abs(denominator) <= 100.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(wave_speed), abs(contact_speed))) return
    factor = (wave_speed - u) / denominator
    rho_star = rho * factor
    if (.not. ieee_is_finite(rho_star) .or. rho_star <= density_floor) return

    star(irho) = rho_star
    star(imx) = rho_star * contact_speed
    star(imy) = rho_star * v
    star(imz) = rho_star * w
    star(iet) = ((wave_speed - u) * state(iet) - pressure * u + &
      star_pressure * contact_speed) / denominator
    do k = 1, nspecies - 1
      star(reactive_species_component(k)) = rho_star * &
        primitive(reactive_mass_fraction_component(k))
    end do
    star(reactive_species_component(nspecies)) = rho_star - &
      sum(star(reactive_species_component(1): &
        reactive_species_component(nspecies - 1)))
    ok = all(ieee_is_finite(star)) .and. star(iet) > 0.0_dp
  end subroutine build_reactive_hllc_star

  subroutine compute_reactive_riemann_flux_x( &
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
    case ("pelec")
      call reactive_pelec_flux_x( &
        species, left_state, right_state, left_temperature_guess, &
        right_temperature_guess, flux, ok)
    case default
      flux = 0.0_dp
      ok = .false.
    end select
  end subroutine compute_reactive_riemann_flux_x

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

  subroutine reconstruct_faces( &
      species, state, temperature, nx, reconstruction, limiter, boundary, &
      dtdx, left_faces, right_faces, left_t, right_t, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: reconstruction, limiter, boundary
    real(dp), intent(in) :: dtdx
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

  subroutine advance_reactive_hydro( &
      species, state, temperature, nx, dx, dt, reconstruction, limiter, &
      riemann_solver, boundary, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt
    character(len=*), intent(in) :: reconstruction, limiter, riemann_solver, boundary
    logical, intent(out) :: ok
    real(dp), allocatable :: lf(:, :), rf(:, :), flux(:, :), updated(:, :)
    real(dp), allocatable :: lt(:), rt(:), q(:)
    real(dp) :: local_t, c
    logical :: local_ok
    integer :: nvar, i

    ok = .false.
    nvar = reactive_nvar(size(species))
    call fill_ghosts(state, temperature, nx, boundary, local_ok)
    if (.not. local_ok) return
    allocate(lf(nvar, 0:nx), rf(nvar, 0:nx), flux(nvar, 0:nx))
    allocate(updated(nvar, 1:nx), lt(0:nx), rt(0:nx))
    allocate(q(reactive_nprim(size(species))))
    call reconstruct_faces(species, state, temperature, nx, reconstruction, &
      limiter, boundary, dt / dx, lf, rf, lt, rt, local_ok)
    if (.not. local_ok) return
    do i = 0, nx
      call compute_reactive_riemann_flux_x( &
        species, lf(:, i), rf(:, i), lt(i), rt(i), riemann_solver, &
        flux(:, i), local_ok)
      if (.not. local_ok) return
    end do
    do i = 1, nx
      updated(:, i) = state(:, i) - dt / dx * (flux(:, i) - flux(:, i - 1))
      call reactive_conserved_to_primitive( &
        species, updated(:, i), temperature(i), q, local_t, c, local_ok)
      if (.not. local_ok) return
      temperature(i) = local_t
    end do
    state(:, 1:nx) = updated
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
    integer :: i, k, substeps

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
        call advance_constant_volume_adaptive( &
          species, reactions, rho, target_energy, request, rtol, atol, y, &
          temperature(i), accepted, next_step, local_ok)
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
      limiter, riemann_solver, boundary, chemistry_enabled, rtol, atol, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt, rtol, atol
    character(len=*), intent(in) :: reconstruction, limiter, riemann_solver, boundary
    logical, intent(in) :: chemistry_enabled
    logical, intent(out) :: ok
    logical :: local_ok

    ok = .false.
    if (chemistry_enabled) then
      call advance_reactive_chemistry(species, reactions, state, temperature, &
        nx, 0.5_dp * dt, rtol, atol, boundary, local_ok)
      if (.not. local_ok) return
    end if
    call advance_reactive_hydro(species, state, temperature, nx, dx, dt, &
      reconstruction, limiter, riemann_solver, boundary, local_ok)
    if (.not. local_ok) return
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
    real(dp), allocatable :: q(:), xmol(:), y(:)
    real(dp) :: x, rho, base_rho, local_temperature, c, gaussian
    logical :: local_ok
    integer :: i, k, nvar

    ok = .false.
    if (size(species) /= 7) return
    nvar = reactive_nvar(size(species))
    allocate(state(nvar, 0:config%nx + 1), temperature(0:config%nx + 1))
    allocate(q(reactive_nprim(size(species))), xmol(7), y(7))
    dx = (config%x_upper - config%x_lower) / real(config%nx, dp)
    xmol = [config%x_h2, config%x_h, config%x_o, config%x_o2, config%x_oh, &
      config%x_h2o, config%x_n2]
    call mass_fractions_from_mole_fractions(species, xmol, y, local_ok)
    if (.not. local_ok) return
    base_rho = mixture_density(species, y, config%initial_pressure, &
      config%initial_temperature, local_ok)
    if (.not. local_ok) return

    do i = 1, config%nx
      x = config%x_lower + (real(i, dp) - 0.5_dp) * dx
      select case (trim(config%problem))
      case ("entropy_wave")
        rho = base_rho * (1.0_dp + config%density_wave_amplitude * &
          sin(2.0_dp * pi * (x - config%x_lower) / &
            (config%x_upper - config%x_lower)))
      case ("uniform_reactor")
        rho = base_rho
      case ("reactive_hotspot")
        gaussian = exp(-((x - config%hotspot_center) / config%hotspot_width)**2)
        local_temperature = config%initial_temperature + &
          config%hotspot_temperature_rise * gaussian
        rho = mixture_density(species, y, config%initial_pressure, &
          local_temperature, local_ok)
        if (.not. local_ok) return
      case default
        return
      end select
      q(1:5) = [rho, config%initial_velocity, 0.0_dp, 0.0_dp, &
        config%initial_pressure]
      do k = 1, size(species)
        q(reactive_mass_fraction_component(k)) = y(k)
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
      initial_integrals, final_integrals, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: state(:, :), temperature(:)
    real(dp), intent(out) :: dx, time
    integer, intent(out) :: steps
    real(dp), intent(out) :: initial_integrals(5), final_integrals(5)
    logical, intent(out) :: ok
    real(dp) :: dt, tolerance
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
        config%cfl, dt, local_ok)
      if (.not. local_ok) then
        ok = .false.; return
      end if
      dt = min(dt, config%final_time - time)
      call advance_reactive_strang(species, reactions, state, temperature, &
        config%nx, dx, dt, config%reconstruction, config%limiter, &
        config%riemann_solver, config%boundary_condition, &
        config%chemistry_enabled, &
        config%chemistry_relative_tolerance, config%chemistry_absolute_tolerance, &
        local_ok)
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
    write(unit, '(a)') &
      "time,x,rho,u,v,w,pressure,temperature,rhoE," // &
      "Y_H2,Y_H,Y_O,Y_O2,Y_OH,Y_H2O,Y_N2"
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

end module reactive_1d_mod
