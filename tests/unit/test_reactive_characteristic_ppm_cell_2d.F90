program test_reactive_characteristic_ppm_cell_2d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nprim, reactive_nvar, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved
  use reactive_2d_mod, only: reconstruct_reactive_characteristic_ppm_cell_2d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: stencil(:, :), minus_plain(:), plus_plain(:)
  real(dp), allocatable :: minus_modified(:), plus_modified(:), conserved(:)
  real(dp) :: base_mole(7), high_mole(7), low_mole(7)
  real(dp) :: base_y(7), high_y(7), low_y(7)
  real(dp) :: temperature, sound_speed, density, difference
  logical :: ok
  integer :: offset, k, component

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "2D PPM-cell thermodynamic database load failed"
  allocate(stencil(reactive_nprim(7), -3:3))
  allocate(minus_plain(reactive_nprim(7)), plus_plain(reactive_nprim(7)))
  allocate(minus_modified(reactive_nprim(7)), plus_modified(reactive_nprim(7)))
  allocate(conserved(reactive_nvar(7)))
  base_mole = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions(species, base_mole, base_y, ok)
  if (.not. ok) error stop "2D PPM-cell base composition failed"

  ! Strong compressive pressure profile: shock flattening must change the
  ! time-traced cell-edge states without violating positivity.
  do offset = -3, 3
    if (offset < 0) then
      stencil(5, offset) = 3.0e5_dp
      stencil(2, offset) = 100.0_dp
    else if (offset > 0) then
      stencil(5, offset) = 1.0e5_dp
      stencil(2, offset) = -100.0_dp
    else
      stencil(5, offset) = 2.0e5_dp
      stencil(2, offset) = 0.0_dp
    end if
    density = mixture_density( &
      species, base_y, stencil(5, offset), 1000.0_dp, ok)
    if (.not. ok) error stop "2D PPM-cell shock density failed"
    stencil(1, offset) = density
    stencil(3, offset) = 0.0_dp
    stencil(4, offset) = 0.0_dp
    do k = 1, size(species)
      stencil(reactive_mass_fraction_component(k), offset) = base_y(k)
    end do
  end do
  call reactive_primitive_to_conserved( &
    species, stencil(:, 0), conserved, temperature, sound_speed, ok)
  if (.not. ok) error stop "2D PPM-cell shock center failed"
  call reconstruct_reactive_characteristic_ppm_cell_2d( &
    species, stencil, sound_speed, 1.0e-7_dp, .false., .false., &
    minus_plain, plus_plain, ok)
  if (.not. ok) error stop "2D unflattened PPM cell failed"
  call reconstruct_reactive_characteristic_ppm_cell_2d( &
    species, stencil, sound_speed, 1.0e-7_dp, .false., .true., &
    minus_modified, plus_modified, ok)
  if (.not. ok) error stop "2D flattened PPM cell failed"
  difference = max(maxval(abs(minus_modified - minus_plain)), &
    maxval(abs(plus_modified - plus_plain)))
  if (difference <= 1.0e-8_dp) &
    error stop "2D shock flattening has no face-state signature"
  if (min(minus_modified(1), plus_modified(1), minus_modified(5), &
      plus_modified(5)) <= 0.0_dp) &
    error stop "2D shock flattening lost positivity"

  ! Constant-pressure material interface: the bounded contact detector must
  ! alter density/species edges while keeping normalized mass fractions.
  high_mole = [0.35570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.49643_dp]
  low_mole = [0.23570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.61643_dp]
  call mass_fractions_from_mole_fractions(species, high_mole, high_y, ok)
  if (.not. ok) error stop "2D high contact composition failed"
  call mass_fractions_from_mole_fractions(species, low_mole, low_y, ok)
  if (.not. ok) error stop "2D low contact composition failed"
  do offset = -3, 3
    stencil(5, offset) = 101325.0_dp
    stencil(2:4, offset) = [150.0_dp, 0.0_dp, 0.0_dp]
    if (offset <= 0) then
      density = mixture_density( &
        species, high_y, 101325.0_dp, 1000.0_dp, ok)
      stencil(1, offset) = density
      do k = 1, size(species)
        stencil(reactive_mass_fraction_component(k), offset) = high_y(k)
      end do
    else
      density = mixture_density( &
        species, low_y, 101325.0_dp, 1000.0_dp, ok)
      stencil(1, offset) = density
      do k = 1, size(species)
        stencil(reactive_mass_fraction_component(k), offset) = low_y(k)
      end do
    end if
    if (.not. ok) error stop "2D contact density failed"
  end do
  call reactive_primitive_to_conserved( &
    species, stencil(:, 0), conserved, temperature, sound_speed, ok)
  if (.not. ok) error stop "2D PPM-cell contact center failed"
  call reconstruct_reactive_characteristic_ppm_cell_2d( &
    species, stencil, sound_speed, 1.0e-7_dp, .false., .false., &
    minus_plain, plus_plain, ok)
  if (.not. ok) error stop "2D plain contact PPM cell failed"
  call reconstruct_reactive_characteristic_ppm_cell_2d( &
    species, stencil, sound_speed, 1.0e-7_dp, .true., .false., &
    minus_modified, plus_modified, ok)
  if (.not. ok) error stop "2D steepened contact PPM cell failed"
  component = reactive_mass_fraction_component(1)
  if (abs(plus_modified(component) - plus_plain(component)) <= 1.0e-9_dp) &
    error stop "2D contact steepening has no H2 face-state signature"
  if (abs(sum(minus_modified(6:)) - 1.0_dp) > 5.0e-13_dp .or. &
      abs(sum(plus_modified(6:)) - 1.0_dp) > 5.0e-13_dp) &
    error stop "2D contact steepening lost composition closure"

  call reconstruct_reactive_characteristic_ppm_cell_2d( &
    species, stencil, sound_speed, 2.0_dp / sound_speed, .false., .false., &
    minus_modified, plus_modified, ok)
  if (ok) error stop "2D characteristic PPM accepted a super-CFL profile"
end program test_reactive_characteristic_ppm_cell_2d
