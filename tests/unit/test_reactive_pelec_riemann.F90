program test_reactive_pelec_riemann
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_species_component, reactive_primitive_to_conserved, &
    reactive_pelec_flux_x, reactive_pelec_primitive_flux_x, reactive_riemann_flux_x
  use reactive_2d_mod, only: reactive_riemann_flux_y
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: ql(:), qr(:), ul(:), ur(:), bad_state(:)
  real(dp), allocatable :: flux(:), selected_flux(:), expected(:)
  real(dp) :: xl(7), xr(7), yl(7), yr(7)
  real(dp) :: pressure, rho_l, rho_r, tl, tr, cl, cr
  real(dp) :: interface_density, interface_velocity, interface_pressure
  real(dp) :: velocity, transverse_y, transverse_z
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load PeleC Riemann thermodynamics"
  allocate(ql(reactive_nprim(7)), qr(reactive_nprim(7)))
  allocate(ul(reactive_nvar(7)), ur(reactive_nvar(7)))
  allocate(bad_state(reactive_nvar(7)))
  allocate(flux(reactive_nvar(7)), selected_flux(reactive_nvar(7)))
  allocate(expected(reactive_nvar(7)))

  xl = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  xr = xl
  xr(1) = xr(1) + 0.05_dp
  xr(7) = xr(7) - 0.05_dp
  call mass_fractions_from_mole_fractions(species, xl, yl, ok)
  if (.not. ok) error stop "Failed to construct left PeleC composition"
  call mass_fractions_from_mole_fractions(species, xr, yr, ok)
  if (.not. ok) error stop "Failed to construct right PeleC composition"

  pressure = 1.5e5_dp
  rho_l = mixture_density(species, yl, pressure, 1000.0_dp, ok)
  if (.not. ok) error stop "Failed to construct left PeleC density"
  rho_r = mixture_density(species, yr, pressure, 1000.0_dp, ok)
  if (.not. ok) error stop "Failed to construct right PeleC density"

  call make_state(rho_l, 0.0_dp, 0.0_dp, 0.0_dp, pressure, yl, ql, ul, tl, cl)
  call make_state(rho_r, 0.0_dp, 0.0_dp, 0.0_dp, pressure, yr, qr, ur, tr, cr)
  call reactive_pelec_flux_x( &
    species, ul, ur, tl, tr, flux, ok, interface_density, &
    interface_velocity, interface_pressure)
  if (.not. ok) error stop "Stationary material-contact PeleC flux failed"
  expected = 0.0_dp
  expected(imx) = pressure
  call assert_vector_close(flux, expected, 2.0e-10_dp, &
    "stationary material contact")
  call assert_close(interface_density, 0.5_dp * (rho_l + rho_r), &
    2.0e-12_dp, "stationary interface density")
  call assert_close(interface_velocity, 0.0_dp, 2.0e-12_dp, &
    "stationary interface velocity")
  call assert_close(interface_pressure, pressure, 2.0e-12_dp, &
    "stationary interface pressure")

  velocity = 80.0_dp
  transverse_y = -15.0_dp
  transverse_z = 7.0_dp
  call make_state( &
    rho_l, velocity, transverse_y, transverse_z, pressure, yl, ql, ul, tl, cl)
  call make_state( &
    rho_r, velocity, transverse_y, transverse_z, pressure, yr, qr, ur, tr, cr)
  call reactive_pelec_flux_x(species, ul, ur, tl, tr, flux, ok)
  if (.not. ok) error stop "Moving material-contact PeleC flux failed"
  call physical_flux(ul, pressure, velocity, transverse_y, transverse_z, yl, expected)
  call assert_vector_close(flux, expected, 2.0e-10_dp, &
    "moving material contact")
  call assert_close(sum(flux(6:12)), flux(irho), 2.0e-12_dp, &
    "moving-contact species closure")

  call reactive_pelec_flux_x( &
    species, ul, ul, tl, tl, flux, ok, interface_density, &
    interface_velocity, interface_pressure)
  if (.not. ok) error stop "Equal-state PeleC flux failed"
  call assert_vector_close(flux, expected, 2.0e-10_dp, &
    "equal-state consistency")
  call assert_close(interface_density, rho_l, 2.0e-12_dp, &
    "equal-state interface density")
  call assert_close(interface_velocity, velocity, 2.0e-12_dp, &
    "equal-state interface velocity")
  call assert_close(interface_pressure, pressure, 2.0e-12_dp, &
    "equal-state interface pressure")

  rho_l = mixture_density(species, yl, 5.0e5_dp, 1400.0_dp, ok)
  if (.not. ok) error stop "Failed to construct left shock density"
  rho_r = mixture_density(species, yr, 8.0e4_dp, 700.0_dp, ok)
  if (.not. ok) error stop "Failed to construct right shock density"
  call make_state(rho_l, 150.0_dp, 4.0_dp, -3.0_dp, &
    5.0e5_dp, yl, ql, ul, tl, cl)
  call make_state(rho_r, -20.0_dp, -6.0_dp, 2.0_dp, &
    8.0e4_dp, yr, qr, ur, tr, cr)
  call reactive_pelec_flux_x( &
    species, ul, ur, tl, tr, flux, ok, interface_density, &
    interface_velocity, interface_pressure)
  if (.not. ok) error stop "General-EOS PeleC shock flux failed"
  if (interface_density <= 0.0_dp .or. interface_pressure <= 0.0_dp) &
    error stop "General-EOS PeleC shock interface is non-physical"
  if (.not. all(ieee_is_finite(flux))) &
    error stop "General-EOS PeleC shock flux is non-finite"
  call assert_close(sum(flux(6:12)), flux(irho), 2.0e-12_dp, &
    "shock species closure")
  call reactive_riemann_flux_x( &
    species, ul, ur, tl, tr, "pelec", selected_flux, ok)
  if (.not. ok) error stop "Reactive PeleC dispatcher failed"
  call assert_vector_close(selected_flux, flux, 2.0e-12_dp, &
    "PeleC dispatcher")

  velocity = 12.0_dp
  transverse_y = -25.0_dp
  transverse_z = 3.0_dp
  call make_state( &
    rho_l, velocity, transverse_y, transverse_z, 5.0e5_dp, &
    yl, ql, ul, tl, cl)
  call reactive_riemann_flux_y( &
    species, ul, ul, tl, tl, "pelec", flux, ok)
  if (.not. ok) error stop "Equal-state y-direction PeleC flux failed"
  call physical_flux_y( &
    ul, 5.0e5_dp, velocity, transverse_y, transverse_z, yl, expected)
  call assert_vector_close(flux, expected, 2.0e-10_dp, &
    "y-direction equal-state consistency")

  bad_state = ul
  bad_state(reactive_species_component(7)) = 0.0_dp
  call reactive_pelec_flux_x( &
    species, bad_state, ur, tl, tr, flux, ok, interface_density, &
    interface_velocity, interface_pressure)
  if (ok) error stop "Invalid reactive state was accepted by PeleC solver"
  call assert_close(interface_density, 0.0_dp, 0.0_dp, &
    "failed interface density reset")
  call assert_close(interface_velocity, 0.0_dp, 0.0_dp, &
    "failed interface velocity reset")
  call assert_close(interface_pressure, 0.0_dp, 0.0_dp, &
    "failed interface pressure reset")

  call check_primitive_entry()
  write(*, '(a)') "test_reactive_pelec_riemann: PASS"

contains

  subroutine check_primitive_entry()
    real(dp) :: q(12),u(12),f(12),exact(12),t,c,rf,vf,pf,reference(12),right(12),uright(12),t2,c2
    logical :: good
    ! A tiny moving contact must retain its supplied speed: the conserved API's
    ! EOS-roundtrip stationary band must not affect pressure-traced primitives.
    call make_state(rho_l,1e-8_dp,-15._dp,7._dp,1.5e5_dp,yl,q,u,t,c)
    call reactive_pelec_primitive_flux_x(species,q,q,c,f,good,rf,vf,pf)
    if (.not.good) error stop 'Primitive tiny-velocity Riemann flux rejected'
    call physical_flux(u,q(5),q(2),q(3),q(4),yl,exact)
    call assert_vector_close(f,exact,2e-12_dp,'primitive equal-state flux')
    if (abs(vf-q(2))>1e-14_dp*abs(q(2))) error stop 'Primitive contact velocity clipped'
    if (abs(f(1)-rho_l*q(2))>1e-14_dp*abs(rho_l*q(2))) error stop 'Primitive mass flux clipped'
    if (pf/=q(5)) error stop 'Primitive input pressure changed'
    ! A separately traced pressure changes the momentum flux without an energy
    ! inversion of a stale conserved normal predictor.
    q(2)=0;q(5)=1.6e5_dp
    call reactive_pelec_primitive_flux_x(species,q,q,c,f,good,rf,vf,pf)
    if (.not.good.or.pf/=q(5).or.f(imx)/=q(5)) error stop 'Traced pressure not preserved'
    if (any(f([irho,imy,imz,iet])/=0)) error stop 'Stationary primitive pressure work'
    ! Nondegenerate unequal-state solve retains the existing general-EOS result.
    call make_state(rho_l,80._dp,-15._dp,7._dp,1.5e5_dp,yl,q,u,t,c)
    call make_state(rho_r,-20._dp,12._dp,-5._dp,1.2e5_dp,yr,right,uright,t2,c2)
    call reactive_pelec_primitive_flux_x(species,q,right,.5_dp*(c+c2),f,good)
    if (.not.good) error stop 'Primitive unequal-state flux rejected'
    call reactive_pelec_flux_x(species,u,uright,t,t2,reference,good)
    if (.not.good) error stop 'Conserved unequal-state reference rejected'
    call assert_vector_close(f,reference,2e-10_dp,'primitive general-EOS reference')
    q(5)=-1
    call reactive_pelec_primitive_flux_x(species,q,q,c,f,good,rf,vf,pf)
    if (good.or.any(f/=0).or.rf/=0.or.vf/=0.or.pf/=0) error stop 'Invalid primitive outputs not zero'
    call reactive_pelec_primitive_flux_x(species,ql,qr,0._dp,f,good)
    if (good.or.any(f/=0)) error stop 'Invalid cell sound speed accepted'
  end subroutine check_primitive_entry

  subroutine make_state(rho, u, v, w, p, y, q, conserved, t, c)
    real(dp), intent(in) :: rho, u, v, w, p, y(:)
    real(dp), intent(out) :: q(:), conserved(:), t, c
    logical :: local_ok
    integer :: species_index

    q(1:5) = [rho, u, v, w, p]
    do species_index = 1, size(y)
      q(reactive_mass_fraction_component(species_index)) = y(species_index)
    end do
    call reactive_primitive_to_conserved( &
      species, q, conserved, t, c, local_ok)
    if (.not. local_ok) error stop "Failed to construct PeleC test state"
  end subroutine make_state

  subroutine physical_flux(state, p, u, v, w, y, result_flux)
    real(dp), intent(in) :: state(:), p, u, v, w, y(:)
    real(dp), intent(out) :: result_flux(:)
    real(dp) :: local_mass_flux
    integer :: species_index

    local_mass_flux = state(irho) * u
    result_flux = 0.0_dp
    result_flux(irho) = local_mass_flux
    result_flux(imx) = local_mass_flux * u + p
    result_flux(imy) = local_mass_flux * v
    result_flux(imz) = local_mass_flux * w
    result_flux(iet) = (state(iet) + p) * u
    do species_index = 1, size(y) - 1
      result_flux(reactive_species_component(species_index)) = &
        local_mass_flux * y(species_index)
    end do
    result_flux(reactive_species_component(size(y))) = local_mass_flux - &
      sum(result_flux(reactive_species_component(1): &
        reactive_species_component(size(y) - 1)))
  end subroutine physical_flux

  subroutine physical_flux_y(state, p, u, v, w, y, result_flux)
    real(dp), intent(in) :: state(:), p, u, v, w, y(:)
    real(dp), intent(out) :: result_flux(:)
    real(dp) :: local_mass_flux
    integer :: species_index

    local_mass_flux = state(irho) * v
    result_flux = 0.0_dp
    result_flux(irho) = local_mass_flux
    result_flux(imx) = local_mass_flux * u
    result_flux(imy) = local_mass_flux * v + p
    result_flux(imz) = local_mass_flux * w
    result_flux(iet) = (state(iet) + p) * v
    do species_index = 1, size(y) - 1
      result_flux(reactive_species_component(species_index)) = &
        local_mass_flux * y(species_index)
    end do
    result_flux(reactive_species_component(size(y))) = local_mass_flux - &
      sum(result_flux(reactive_species_component(1): &
        reactive_species_component(size(y) - 1)))
  end subroutine physical_flux_y

  subroutine assert_close(actual, expected_value, relative_tolerance, label)
    real(dp), intent(in) :: actual, expected_value, relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: error

    error = abs(actual - expected_value) / max(1.0_dp, abs(expected_value))
    if (error > relative_tolerance) then
      write(*, '(a,3(1x,es24.16))') &
        trim(label), actual, expected_value, error
      error stop "Reactive PeleC scalar mismatch"
    end if
  end subroutine assert_close

  subroutine assert_vector_close(actual, expected_value, relative_tolerance, label)
    real(dp), intent(in) :: actual(:), expected_value(:), relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: error

    error = maxval(abs(actual - expected_value) / max(1.0_dp, abs(expected_value)))
    if (error > relative_tolerance) then
      write(*, '(a,1x,es24.16)') trim(label), error
      error stop "Reactive PeleC vector mismatch"
    end if
  end subroutine assert_vector_close

end program test_reactive_pelec_riemann
