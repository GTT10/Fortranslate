module reactive_transport_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species, nasa7_mass_properties
  use mixture_thermo_mod, only: &
    mole_fractions_from_mass_fractions, mixture_mass_properties
  use transport_database_mod, only: gas_transport_species
  use mixture_transport_mod, only: mixture_transport_coefficients
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_conserved_to_primitive
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, initialize_periodic_boundary_set_2d, &
    sample_reactive_primitive_2d, reactive_boundary_is_periodic, &
    reactive_boundary_is_wall, reactive_boundary_is_inflow
  implicit none
  private

  real(dp), parameter :: species_safety = 0.90_dp

  public :: reactive_transport_timestep_2d
  public :: reactive_transport_fluxes_2d_faces
  public :: reactive_transport_fluxes_2d
  public :: reactive_transport_euler_update_2d
  public :: advance_reactive_transport_2d

contains

  pure integer function periodic_index(index, extent) result(wrapped)
    integer, intent(in) :: index, extent
    wrapped = 1 + modulo(index - 1, extent)
  end function periodic_index

  subroutine recover_primitives_2d( &
      species, state, temperature, nx, ny, primitive, checked_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(out) :: primitive(:, :, :)
    real(dp), intent(out) :: checked_temperature(:, :)
    logical, intent(out) :: ok

    real(dp) :: sound_speed
    logical :: local_ok
    integer :: i, j

    ok = .false.
    if (size(state, 1) /= reactive_nvar(size(species)) .or. &
        size(state, 2) < nx .or. size(state, 3) < ny .or. &
        size(temperature, 1) < nx .or. size(temperature, 2) < ny .or. &
        size(primitive, 1) /= reactive_nprim(size(species)) .or. &
        size(primitive, 2) < nx .or. size(primitive, 3) < ny .or. &
        size(checked_temperature, 1) < nx .or. &
        size(checked_temperature, 2) < ny) return
    do j = 1, ny
      do i = 1, nx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive(:, i, j), &
          checked_temperature(i, j), sound_speed, local_ok)
        if (.not. local_ok) return
      end do
    end do
    ok = .true.
  end subroutine recover_primitives_2d

  subroutine face_transport_data( &
      species, transport, left_primitive, right_primitive, left_temperature, &
      right_temperature, viscosity, conductivity, diffusion, yleft, yright, &
      yface, xleft, xright, hface, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: left_primitive(:), right_primitive(:)
    real(dp), intent(in) :: left_temperature, right_temperature
    real(dp), intent(out) :: viscosity, conductivity, diffusion(:)
    real(dp), intent(out) :: yleft(:), yright(:), yface(:)
    real(dp), intent(out) :: xleft(:), xright(:), hface(:)
    logical, intent(out) :: ok

    real(dp) :: tface, pface, cp, cv, hleft, hright, eint, entropy
    logical :: local_ok
    integer :: k, nspecies

    viscosity = 0.0_dp
    conductivity = 0.0_dp
    diffusion = 0.0_dp
    yleft = 0.0_dp
    yright = 0.0_dp
    yface = 0.0_dp
    xleft = 0.0_dp
    xright = 0.0_dp
    hface = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (size(transport) /= nspecies .or. &
        size(left_primitive) /= reactive_nprim(nspecies) .or. &
        size(right_primitive) /= reactive_nprim(nspecies) .or. &
        size(diffusion) /= nspecies .or. size(yleft) /= nspecies .or. &
        size(yright) /= nspecies .or. size(yface) /= nspecies .or. &
        size(xleft) /= nspecies .or. size(xright) /= nspecies .or. &
        size(hface) /= nspecies) return

    do k = 1, nspecies
      yleft(k) = left_primitive(reactive_mass_fraction_component(k))
      yright(k) = right_primitive(reactive_mass_fraction_component(k))
    end do
    yface = max(0.0_dp, 0.5_dp * (yleft + yright))
    if (sum(yface) <= tiny(1.0_dp)) return
    yface = yface / sum(yface)
    tface = 0.5_dp * (left_temperature + right_temperature)
    pface = 0.5_dp * (left_primitive(5) + right_primitive(5))
    call mixture_transport_coefficients( &
      species, transport, yface, tface, pface, viscosity, conductivity, &
      diffusion, local_ok)
    if (.not. local_ok) return
    call mole_fractions_from_mass_fractions( &
      species, yleft, xleft, local_ok)
    if (.not. local_ok) return
    call mole_fractions_from_mass_fractions( &
      species, yright, xright, local_ok)
    if (.not. local_ok) return
    do k = 1, nspecies
      call nasa7_mass_properties( &
        species(k), left_temperature, cp, cv, hleft, eint, entropy, local_ok)
      if (.not. local_ok) return
      call nasa7_mass_properties( &
        species(k), right_temperature, cp, cv, hright, eint, entropy, local_ok)
      if (.not. local_ok) return
      hface(k) = 0.5_dp * (hleft + hright)
    end do
    ok = ieee_is_finite(viscosity) .and. viscosity > 0.0_dp .and. &
      ieee_is_finite(conductivity) .and. conductivity > 0.0_dp .and. &
      all(ieee_is_finite(diffusion)) .and. all(diffusion >= 0.0_dp) .and. &
      all(ieee_is_finite(hface))
  end subroutine face_transport_data

  subroutine species_face_flux( &
      species, diffusion, yleft, yright, yface, xleft, xright, hface, &
      density_face, pressure_left, pressure_right, pressure_face, spacing, &
      barodiffusion_enabled, species_flux, enthalpy_flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: diffusion(:), yleft(:), yright(:), yface(:)
    real(dp), intent(in) :: xleft(:), xright(:), hface(:)
    real(dp), intent(in) :: density_face, pressure_left, pressure_right
    real(dp), intent(in) :: pressure_face, spacing
    logical, intent(in) :: barodiffusion_enabled
    real(dp), intent(out) :: species_flux(:), enthalpy_flux
    logical, intent(out) :: ok

    real(dp), allocatable :: raw_flux(:)
    real(dp) :: dlnp, correction
    integer :: k, nspecies

    species_flux = 0.0_dp
    enthalpy_flux = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (size(diffusion) /= nspecies .or. size(yleft) /= nspecies .or. &
        size(yright) /= nspecies .or. size(yface) /= nspecies .or. &
        size(xleft) /= nspecies .or. size(xright) /= nspecies .or. &
        size(hface) /= nspecies .or. size(species_flux) /= nspecies .or. &
        density_face <= 0.0_dp .or. pressure_face <= 0.0_dp .or. &
        spacing <= 0.0_dp) return
    allocate(raw_flux(nspecies))
    dlnp = 0.0_dp
    if (barodiffusion_enabled) then
      dlnp = (pressure_right - pressure_left) / (spacing * pressure_face)
    end if
    do k = 1, nspecies
      raw_flux(k) = -density_face * diffusion(k) * &
        ((xright(k) - xleft(k)) / spacing + &
          (0.5_dp * (xleft(k) + xright(k)) - yface(k)) * dlnp)
    end do
    correction = sum(raw_flux)
    species_flux = raw_flux - yface * correction
    if (nspecies > 1) then
      species_flux(nspecies) = -sum(species_flux(1:nspecies - 1))
    else
      species_flux(1) = 0.0_dp
    end if
    enthalpy_flux = sum(hface * species_flux)
    ok = all(ieee_is_finite(species_flux)) .and. &
      ieee_is_finite(enthalpy_flux) .and. &
      abs(sum(species_flux)) <= 2.0e3_dp * epsilon(1.0_dp) * &
        max(1.0_dp, maxval(abs(species_flux)))
  end subroutine species_face_flux

  subroutine transport_face_flux_x( &
      species, transport, primitive, checked_temperature, nx, ny, face_i, j, &
      dx, dy, boundaries, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, flux, &
      species_energy, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: primitive(:, :, :), checked_temperature(:, :)
    integer, intent(in) :: nx, ny, face_i, j
    real(dp), intent(in) :: dx, dy
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux(:), species_energy
    logical, intent(out) :: ok

    real(dp), allocatable :: qleft(:), qright(:), qtmp(:), qtmp2(:)
    real(dp), allocatable :: yleft(:), yright(:), yface(:)
    real(dp), allocatable :: xleft(:), xright(:), diffusion(:), hface(:)
    real(dp), allocatable :: species_flux(:)
    real(dp) :: tleft, tright, ttmp, ttmp2, spacing
    real(dp) :: viscosity, conductivity, density_face, pressure_face
    real(dp) :: dudx, dvdx, dwdx, dudy, dvdy, divu
    real(dp) :: tau_xx, tau_xy, tau_xz, dtdx, uface, vface, wface
    logical :: local_ok, wall_face, slip_face
    integer :: nspecies, nprim, k, side

    flux = 0.0_dp
    species_energy = 0.0_dp
    ok = .false.
    nspecies = size(species)
    nprim = reactive_nprim(nspecies)
    if (face_i < 0 .or. face_i > nx .or. j < 1 .or. j > ny .or. &
        size(flux) /= reactive_nvar(nspecies)) return
    allocate(qleft(nprim), qright(nprim), qtmp(nprim), qtmp2(nprim))
    allocate(yleft(nspecies), yright(nspecies), yface(nspecies))
    allocate(xleft(nspecies), xright(nspecies), diffusion(nspecies))
    allocate(hface(nspecies), species_flux(nspecies))

    call sample_reactive_primitive_2d( &
      primitive, checked_temperature, nx, ny, face_i, j, boundaries, &
      qleft, tleft, local_ok)
    if (.not. local_ok) return
    call sample_reactive_primitive_2d( &
      primitive, checked_temperature, nx, ny, face_i + 1, j, boundaries, &
      qright, tright, local_ok)
    if (.not. local_ok) return

    spacing = dx
    wall_face = .false.
    slip_face = .false.
    side = 0
    if (face_i == 0 .and. .not. &
        reactive_boundary_is_periodic(boundaries%face(1))) side = 1
    if (face_i == nx .and. .not. &
        reactive_boundary_is_periodic(boundaries%face(2))) side = 2
    if (side > 0) then
      wall_face = reactive_boundary_is_wall(boundaries%face(side))
      slip_face = trim(boundaries%face(side)%kind) == 'slip_wall'
      if (reactive_boundary_is_inflow(boundaries%face(side))) spacing = 0.5_dp * dx
    end if

    call face_transport_data( &
      species, transport, qleft, qright, tleft, tright, viscosity, &
      conductivity, diffusion, yleft, yright, yface, xleft, xright, hface, &
      local_ok)
    if (.not. local_ok) return
    density_face = 0.5_dp * (qleft(1) + qright(1))
    pressure_face = 0.5_dp * (qleft(5) + qright(5))
    dudx = (qright(2) - qleft(2)) / spacing
    dvdx = (qright(3) - qleft(3)) / spacing
    dwdx = (qright(4) - qleft(4)) / spacing

    call sample_reactive_primitive_2d( &
      primitive, checked_temperature, nx, ny, face_i, j + 1, boundaries, &
      qtmp, ttmp, local_ok)
    if (.not. local_ok) return
    call sample_reactive_primitive_2d( &
      primitive, checked_temperature, nx, ny, face_i, j - 1, boundaries, &
      qtmp2, ttmp2, local_ok)
    if (.not. local_ok) return
    dudy = 0.5_dp * (qtmp(2) - qtmp2(2)) / (2.0_dp * dy)
    dvdy = 0.5_dp * (qtmp(3) - qtmp2(3)) / (2.0_dp * dy)
    call sample_reactive_primitive_2d( &
      primitive, checked_temperature, nx, ny, face_i + 1, j + 1, &
      boundaries, qtmp, ttmp, local_ok)
    if (.not. local_ok) return
    call sample_reactive_primitive_2d( &
      primitive, checked_temperature, nx, ny, face_i + 1, j - 1, &
      boundaries, qtmp2, ttmp2, local_ok)
    if (.not. local_ok) return
    dudy = dudy + 0.5_dp * (qtmp(2) - qtmp2(2)) / (2.0_dp * dy)
    dvdy = dvdy + 0.5_dp * (qtmp(3) - qtmp2(3)) / (2.0_dp * dy)
    divu = dudx + dvdy

    if (viscosity_enabled) then
      tau_xx = viscosity * (2.0_dp * dudx - (2.0_dp / 3.0_dp) * divu)
      tau_xy = viscosity * (dudy + dvdx)
      tau_xz = viscosity * dwdx
      if (slip_face) then
        tau_xy = 0.0_dp
        tau_xz = 0.0_dp
      end if
      flux(imx) = -tau_xx
      flux(imy) = -tau_xy
      flux(imz) = -tau_xz
      uface = 0.5_dp * (qleft(2) + qright(2))
      vface = 0.5_dp * (qleft(3) + qright(3))
      wface = 0.5_dp * (qleft(4) + qright(4))
      flux(iet) = -(tau_xx * uface + tau_xy * vface + tau_xz * wface)
    end if
    if (thermal_conduction_enabled) then
      dtdx = (tright - tleft) / spacing
      flux(iet) = flux(iet) - conductivity * dtdx
    end if
    if (species_diffusion_enabled .and. .not. wall_face) then
      call species_face_flux( &
        species, diffusion, yleft, yright, yface, xleft, xright, hface, &
        density_face, qleft(5), qright(5), pressure_face, spacing, &
        barodiffusion_enabled, species_flux, species_energy, local_ok)
      if (.not. local_ok) return
      do k = 1, nspecies
        flux(reactive_species_component(k)) = species_flux(k)
      end do
      flux(iet) = flux(iet) + species_energy
    end if
    flux(irho) = 0.0_dp
    ok = all(ieee_is_finite(flux)) .and. ieee_is_finite(species_energy)
  end subroutine transport_face_flux_x

  subroutine transport_face_flux_y( &
      species, transport, primitive, checked_temperature, nx, ny, i, face_j, &
      dx, dy, boundaries, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, flux, &
      species_energy, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: primitive(:, :, :), checked_temperature(:, :)
    integer, intent(in) :: nx, ny, i, face_j
    real(dp), intent(in) :: dx, dy
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux(:), species_energy
    logical, intent(out) :: ok

    real(dp), allocatable :: qlower(:), qupper(:), qtmp(:), qtmp2(:)
    real(dp), allocatable :: ylower(:), yupper(:), yface(:)
    real(dp), allocatable :: xlower(:), xupper(:), diffusion(:), hface(:)
    real(dp), allocatable :: species_flux(:)
    real(dp) :: tlower, tupper, ttmp, ttmp2, spacing
    real(dp) :: viscosity, conductivity, density_face, pressure_face
    real(dp) :: dudy, dvdy, dwdy, dudx, dvdx, divu
    real(dp) :: tau_yx, tau_yy, tau_yz, dtdy, uface, vface, wface
    logical :: local_ok, wall_face, slip_face
    integer :: nspecies, nprim, k, side

    flux = 0.0_dp
    species_energy = 0.0_dp
    ok = .false.
    nspecies = size(species)
    nprim = reactive_nprim(nspecies)
    if (i < 1 .or. i > nx .or. face_j < 0 .or. face_j > ny .or. &
        size(flux) /= reactive_nvar(nspecies)) return
    allocate(qlower(nprim), qupper(nprim), qtmp(nprim), qtmp2(nprim))
    allocate(ylower(nspecies), yupper(nspecies), yface(nspecies))
    allocate(xlower(nspecies), xupper(nspecies), diffusion(nspecies))
    allocate(hface(nspecies), species_flux(nspecies))

    call sample_reactive_primitive_2d( &
      primitive, checked_temperature, nx, ny, i, face_j, boundaries, &
      qlower, tlower, local_ok)
    if (.not. local_ok) return
    call sample_reactive_primitive_2d( &
      primitive, checked_temperature, nx, ny, i, face_j + 1, boundaries, &
      qupper, tupper, local_ok)
    if (.not. local_ok) return

    spacing = dy
    wall_face = .false.
    slip_face = .false.
    side = 0
    if (face_j == 0 .and. .not. &
        reactive_boundary_is_periodic(boundaries%face(3))) side = 3
    if (face_j == ny .and. .not. &
        reactive_boundary_is_periodic(boundaries%face(4))) side = 4
    if (side > 0) then
      wall_face = reactive_boundary_is_wall(boundaries%face(side))
      slip_face = trim(boundaries%face(side)%kind) == 'slip_wall'
      if (reactive_boundary_is_inflow(boundaries%face(side))) spacing = 0.5_dp * dy
    end if

    call face_transport_data( &
      species, transport, qlower, qupper, tlower, tupper, viscosity, &
      conductivity, diffusion, ylower, yupper, yface, xlower, xupper, hface, &
      local_ok)
    if (.not. local_ok) return
    density_face = 0.5_dp * (qlower(1) + qupper(1))
    pressure_face = 0.5_dp * (qlower(5) + qupper(5))
    dudy = (qupper(2) - qlower(2)) / spacing
    dvdy = (qupper(3) - qlower(3)) / spacing
    dwdy = (qupper(4) - qlower(4)) / spacing

    call sample_reactive_primitive_2d( &
      primitive, checked_temperature, nx, ny, i + 1, face_j, boundaries, &
      qtmp, ttmp, local_ok)
    if (.not. local_ok) return
    call sample_reactive_primitive_2d( &
      primitive, checked_temperature, nx, ny, i - 1, face_j, boundaries, &
      qtmp2, ttmp2, local_ok)
    if (.not. local_ok) return
    dudx = 0.5_dp * (qtmp(2) - qtmp2(2)) / (2.0_dp * dx)
    dvdx = 0.5_dp * (qtmp(3) - qtmp2(3)) / (2.0_dp * dx)
    call sample_reactive_primitive_2d( &
      primitive, checked_temperature, nx, ny, i + 1, face_j + 1, &
      boundaries, qtmp, ttmp, local_ok)
    if (.not. local_ok) return
    call sample_reactive_primitive_2d( &
      primitive, checked_temperature, nx, ny, i - 1, face_j + 1, &
      boundaries, qtmp2, ttmp2, local_ok)
    if (.not. local_ok) return
    dudx = dudx + 0.5_dp * (qtmp(2) - qtmp2(2)) / (2.0_dp * dx)
    dvdx = dvdx + 0.5_dp * (qtmp(3) - qtmp2(3)) / (2.0_dp * dx)
    divu = dudx + dvdy

    if (viscosity_enabled) then
      tau_yx = viscosity * (dudy + dvdx)
      tau_yy = viscosity * (2.0_dp * dvdy - (2.0_dp / 3.0_dp) * divu)
      tau_yz = viscosity * dwdy
      if (slip_face) then
        tau_yx = 0.0_dp
        tau_yz = 0.0_dp
      end if
      flux(imx) = -tau_yx
      flux(imy) = -tau_yy
      flux(imz) = -tau_yz
      uface = 0.5_dp * (qlower(2) + qupper(2))
      vface = 0.5_dp * (qlower(3) + qupper(3))
      wface = 0.5_dp * (qlower(4) + qupper(4))
      flux(iet) = -(tau_yx * uface + tau_yy * vface + tau_yz * wface)
    end if
    if (thermal_conduction_enabled) then
      dtdy = (tupper - tlower) / spacing
      flux(iet) = flux(iet) - conductivity * dtdy
    end if
    if (species_diffusion_enabled .and. .not. wall_face) then
      call species_face_flux( &
        species, diffusion, ylower, yupper, yface, xlower, xupper, hface, &
        density_face, qlower(5), qupper(5), pressure_face, spacing, &
        barodiffusion_enabled, species_flux, species_energy, local_ok)
      if (.not. local_ok) return
      do k = 1, nspecies
        flux(reactive_species_component(k)) = species_flux(k)
      end do
      flux(iet) = flux(iet) + species_energy
    end if
    flux(irho) = 0.0_dp
    ok = all(ieee_is_finite(flux)) .and. ieee_is_finite(species_energy)
  end subroutine transport_face_flux_y

  subroutine reactive_transport_fluxes_2d_faces( &
      species, transport, state, temperature, nx, ny, dx, dy, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, flux_x, flux_y, &
      minimum_theta, ok, boundaries)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux_x(:, 0:, :), flux_y(:, :, 0:)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries

    type(reactive_boundary_set_2d) :: active_boundaries
    real(dp), allocatable :: primitive(:, :, :), checked_temperature(:, :)
    real(dp), allocatable :: species_energy_x(:, :), species_energy_y(:, :)
    real(dp), allocatable :: theta_cell(:, :)
    real(dp) :: outgoing, mass, candidate, theta_face
    logical :: local_ok, periodic_x, periodic_y
    integer :: i, j, face_i, face_j, k, component
    integer :: nspecies, nvar, nprim, left_i, right_i, lower_j, upper_j

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    minimum_theta = 1.0_dp
    ok = .false.
    nspecies = size(species)
    nvar = reactive_nvar(nspecies)
    nprim = reactive_nprim(nspecies)
    if (size(transport) /= nspecies .or. nx < 2 .or. ny < 2 .or. &
        dx <= 0.0_dp .or. dy <= 0.0_dp .or. dt < 0.0_dp .or. &
        size(state, 1) /= nvar .or. size(state, 2) < nx .or. &
        size(state, 3) < ny .or. size(temperature, 1) < nx .or. &
        size(temperature, 2) < ny .or. size(flux_x, 1) /= nvar .or. &
        ubound(flux_x, 2) < nx .or. size(flux_x, 3) < ny .or. &
        size(flux_y, 1) /= nvar .or. size(flux_y, 2) < nx .or. &
        ubound(flux_y, 3) < ny) return
    if (present(boundaries)) then
      active_boundaries = boundaries
    else
      call initialize_periodic_boundary_set_2d(nprim, active_boundaries)
    end if
    periodic_x = reactive_boundary_is_periodic(active_boundaries%face(1))
    periodic_y = reactive_boundary_is_periodic(active_boundaries%face(3))
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    allocate(primitive(nprim, nx, ny), checked_temperature(nx, ny))
    allocate(species_energy_x(0:nx, ny), species_energy_y(nx, 0:ny))
    allocate(theta_cell(nx, ny))
    species_energy_x = 0.0_dp
    species_energy_y = 0.0_dp
    call recover_primitives_2d( &
      species, state, temperature, nx, ny, primitive, checked_temperature, &
      local_ok)
    if (.not. local_ok) return

    do j = 1, ny
      do face_i = 0, nx
        call transport_face_flux_x( &
          species, transport, primitive, checked_temperature, nx, ny, face_i, &
          j, dx, dy, active_boundaries, viscosity_enabled, &
          thermal_conduction_enabled, species_diffusion_enabled, &
          barodiffusion_enabled, flux_x(:, face_i, j), &
          species_energy_x(face_i, j), local_ok)
        if (.not. local_ok) return
      end do
    end do
    do face_j = 0, ny
      do i = 1, nx
        call transport_face_flux_y( &
          species, transport, primitive, checked_temperature, nx, ny, i, &
          face_j, dx, dy, active_boundaries, viscosity_enabled, &
          thermal_conduction_enabled, species_diffusion_enabled, &
          barodiffusion_enabled, flux_y(:, i, face_j), &
          species_energy_y(i, face_j), local_ok)
        if (.not. local_ok) return
      end do
    end do

    if (species_diffusion_enabled .and. dt > 0.0_dp) then
      theta_cell = 1.0_dp
      do j = 1, ny
        do i = 1, nx
          do k = 1, nspecies
            component = reactive_species_component(k)
            outgoing = max(flux_x(component, i, j), 0.0_dp) / dx + &
              max(-flux_x(component, i - 1, j), 0.0_dp) / dx + &
              max(flux_y(component, i, j), 0.0_dp) / dy + &
              max(-flux_y(component, i, j - 1), 0.0_dp) / dy
            mass = max(0.0_dp, state(component, i, j))
            if (outgoing > 0.0_dp) then
              candidate = species_safety * mass / (dt * outgoing)
              theta_cell(i, j) = min(theta_cell(i, j), &
                max(0.0_dp, min(1.0_dp, candidate)))
            end if
          end do
        end do
      end do

      do j = 1, ny
        do face_i = 0, nx
          theta_face = 1.0_dp
          left_i = face_i
          right_i = face_i + 1
          if (face_i == 0) then
            if (periodic_x) left_i = nx
          else if (face_i == nx) then
            if (periodic_x) right_i = 1
          end if
          if (left_i >= 1 .and. left_i <= nx) &
            theta_face = min(theta_face, theta_cell(left_i, j))
          if (right_i >= 1 .and. right_i <= nx) &
            theta_face = min(theta_face, theta_cell(right_i, j))
          minimum_theta = min(minimum_theta, theta_face)
          do k = 1, nspecies
            component = reactive_species_component(k)
            flux_x(component, face_i, j) = &
              theta_face * flux_x(component, face_i, j)
          end do
          flux_x(iet, face_i, j) = flux_x(iet, face_i, j) + &
            (theta_face - 1.0_dp) * species_energy_x(face_i, j)
        end do
      end do
      do face_j = 0, ny
        do i = 1, nx
          theta_face = 1.0_dp
          lower_j = face_j
          upper_j = face_j + 1
          if (face_j == 0) then
            if (periodic_y) lower_j = ny
          else if (face_j == ny) then
            if (periodic_y) upper_j = 1
          end if
          if (lower_j >= 1 .and. lower_j <= ny) &
            theta_face = min(theta_face, theta_cell(i, lower_j))
          if (upper_j >= 1 .and. upper_j <= ny) &
            theta_face = min(theta_face, theta_cell(i, upper_j))
          minimum_theta = min(minimum_theta, theta_face)
          do k = 1, nspecies
            component = reactive_species_component(k)
            flux_y(component, i, face_j) = &
              theta_face * flux_y(component, i, face_j)
          end do
          flux_y(iet, i, face_j) = flux_y(iet, i, face_j) + &
            (theta_face - 1.0_dp) * species_energy_y(i, face_j)
        end do
      end do
    end if

    ok = all(ieee_is_finite(flux_x)) .and. all(ieee_is_finite(flux_y)) .and. &
      ieee_is_finite(minimum_theta) .and. minimum_theta >= 0.0_dp .and. &
      minimum_theta <= 1.0_dp
  end subroutine reactive_transport_fluxes_2d_faces

  subroutine reactive_transport_fluxes_2d( &
      species, transport, state, temperature, nx, ny, dx, dy, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, flux_x, flux_y, &
      minimum_theta, ok, boundaries)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux_x(:, :, :), flux_y(:, :, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries

    real(dp), allocatable :: face_x(:, :, :), face_y(:, :, :)
    logical :: local_ok
    integer :: nvar, i, j

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (size(flux_x, 1) /= nvar .or. size(flux_x, 2) < nx .or. &
        size(flux_x, 3) < ny .or. size(flux_y, 1) /= nvar .or. &
        size(flux_y, 2) < nx .or. size(flux_y, 3) < ny) return
    allocate(face_x(nvar, 0:nx, ny), face_y(nvar, nx, 0:ny))
    if (present(boundaries)) then
      call reactive_transport_fluxes_2d_faces( &
        species, transport, state, temperature, nx, ny, dx, dy, dt, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, face_x, face_y, &
        minimum_theta, local_ok, boundaries)
    else
      call reactive_transport_fluxes_2d_faces( &
        species, transport, state, temperature, nx, ny, dx, dy, dt, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, face_x, face_y, &
        minimum_theta, local_ok)
    end if
    if (.not. local_ok) return
    do j = 1, ny
      do i = 1, nx
        flux_x(:, i, j) = face_x(:, i, j)
        flux_y(:, i, j) = face_y(:, i, j)
      end do
    end do
    ok = .true.
  end subroutine reactive_transport_fluxes_2d

  subroutine reactive_transport_timestep_2d( &
      species, transport, state, temperature, nx, ny, dx, dy, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, maximum_diffusivity, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt, maximum_diffusivity
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:), y(:), diffusion(:)
    real(dp) :: checked_temperature, sound_speed, viscosity, conductivity
    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, internal_energy, entropy, candidate, denominator
    logical :: local_ok
    integer :: i, j, k, nspecies

    dt = 0.0_dp
    maximum_diffusivity = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (size(transport) /= nspecies .or. nx < 1 .or. ny < 1 .or. &
        dx <= 0.0_dp .or. dy <= 0.0_dp .or. transport_cfl <= 0.0_dp .or. &
        transport_cfl > 0.5_dp) return
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) then
      dt = huge(1.0_dp)
      ok = .true.
      return
    end if
    allocate(primitive(reactive_nprim(nspecies)), y(nspecies))
    allocate(diffusion(nspecies))
    do j = 1, ny
      do i = 1, nx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          checked_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        do k = 1, nspecies
          y(k) = primitive(reactive_mass_fraction_component(k))
        end do
        call mixture_transport_coefficients( &
          species, transport, y, checked_temperature, primitive(5), viscosity, &
          conductivity, diffusion, local_ok)
        if (.not. local_ok) return
        call mixture_mass_properties( &
          species, y, checked_temperature, molecular_weight, gas_constant, &
          cp, cv, gamma, enthalpy, internal_energy, entropy, local_ok)
        if (.not. local_ok .or. primitive(1) <= 0.0_dp .or. cv <= 0.0_dp) return
        if (viscosity_enabled) then
          candidate = (4.0_dp / 3.0_dp) * viscosity / primitive(1)
          maximum_diffusivity = max(maximum_diffusivity, candidate)
        end if
        if (thermal_conduction_enabled) then
          candidate = conductivity / (primitive(1) * cv)
          maximum_diffusivity = max(maximum_diffusivity, candidate)
        end if
        if (species_diffusion_enabled) then
          maximum_diffusivity = max(maximum_diffusivity, maxval(diffusion))
        end if
      end do
    end do
    if (maximum_diffusivity <= 0.0_dp) then
      dt = huge(1.0_dp)
    else
      denominator = maximum_diffusivity * (1.0_dp / dx**2 + 1.0_dp / dy**2)
      dt = transport_cfl / denominator
    end if
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp .and. &
      ieee_is_finite(maximum_diffusivity)
  end subroutine reactive_transport_timestep_2d

  subroutine reactive_transport_euler_update_2d( &
      species, transport, input_state, input_temperature, nx, ny, dx, dy, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, output_state, &
      output_temperature, minimum_theta, ok, boundaries)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: input_state(:, :, :), input_temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: output_state(:, :, :), output_temperature(:, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries

    real(dp), allocatable :: flux_x(:, :, :), flux_y(:, :, :), primitive(:)
    real(dp) :: checked_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    if (size(output_state, 1) /= nvar .or. size(output_state, 2) < nx .or. &
        size(output_state, 3) < ny .or. size(output_temperature, 1) < nx .or. &
        size(output_temperature, 2) < ny) return
    allocate(flux_x(nvar, 0:nx, ny), flux_y(nvar, nx, 0:ny))
    allocate(primitive(reactive_nprim(size(species))))
    if (present(boundaries)) then
      call reactive_transport_fluxes_2d_faces( &
        species, transport, input_state, input_temperature, nx, ny, dx, dy, &
        dt, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, flux_x, flux_y, &
        minimum_theta, local_ok, boundaries)
    else
      call reactive_transport_fluxes_2d_faces( &
        species, transport, input_state, input_temperature, nx, ny, dx, dy, &
        dt, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, flux_x, flux_y, &
        minimum_theta, local_ok)
    end if
    if (.not. local_ok) return
    output_state = input_state
    output_temperature = input_temperature
    do j = 1, ny
      do i = 1, nx
        output_state(:, i, j) = input_state(:, i, j) - &
          dt / dx * (flux_x(:, i, j) - flux_x(:, i - 1, j)) - &
          dt / dy * (flux_y(:, i, j) - flux_y(:, i, j - 1))
        call reactive_conserved_to_primitive( &
          species, output_state(:, i, j), input_temperature(i, j), primitive, &
          checked_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        output_temperature(i, j) = checked_temperature
      end do
    end do
    ok = .true.
  end subroutine reactive_transport_euler_update_2d

  subroutine advance_reactive_transport_2d( &
      species, transport, state, temperature, nx, ny, dx, dy, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, minimum_theta, ok, &
      boundaries)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, interval
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries

    real(dp), allocatable :: initial_state(:, :, :), initial_temperature(:, :)
    real(dp), allocatable :: stage1_state(:, :, :), stage1_temperature(:, :)
    real(dp), allocatable :: euler2_state(:, :, :), euler2_temperature(:, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: theta1, theta2, checked_temperature, sound_speed
    real(dp) :: temperature_guess
    logical :: local_ok
    integer :: i, j, nvar

    ok = .false.
    minimum_theta = 1.0_dp
    if (interval < 0.0_dp) return
    if (interval <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      return
    end if
    nvar = reactive_nvar(size(species))
    allocate(initial_state(nvar, nx, ny), initial_temperature(nx, ny))
    allocate(stage1_state(nvar, nx, ny), stage1_temperature(nx, ny))
    allocate(euler2_state(nvar, nx, ny), euler2_temperature(nx, ny))
    allocate(primitive(reactive_nprim(size(species))))
    initial_state = state
    initial_temperature = temperature

    if (present(boundaries)) then
      call reactive_transport_euler_update_2d( &
        species, transport, initial_state, initial_temperature, nx, ny, dx, &
        dy, interval, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, stage1_state, &
        stage1_temperature, theta1, local_ok, boundaries)
    else
      call reactive_transport_euler_update_2d( &
        species, transport, initial_state, initial_temperature, nx, ny, dx, &
        dy, interval, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, stage1_state, &
        stage1_temperature, theta1, local_ok)
    end if
    if (.not. local_ok) return
    if (present(boundaries)) then
      call reactive_transport_euler_update_2d( &
        species, transport, stage1_state, stage1_temperature, nx, ny, dx, dy, &
        interval, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, euler2_state, &
        euler2_temperature, theta2, local_ok, boundaries)
    else
      call reactive_transport_euler_update_2d( &
        species, transport, stage1_state, stage1_temperature, nx, ny, dx, dy, &
        interval, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, euler2_state, &
        euler2_temperature, theta2, local_ok)
    end if
    if (.not. local_ok) return

    do j = 1, ny
      do i = 1, nx
        state(:, i, j) = 0.5_dp * &
          (initial_state(:, i, j) + euler2_state(:, i, j))
        temperature_guess = 0.5_dp * &
          (initial_temperature(i, j) + euler2_temperature(i, j))
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature_guess, primitive, &
          checked_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        temperature(i, j) = checked_temperature
      end do
    end do
    minimum_theta = min(theta1, theta2)
    ok = .true.
  end subroutine advance_reactive_transport_2d


end module reactive_transport_2d_mod
