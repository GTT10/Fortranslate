module reactive_eb_cfl_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  implicit none
  private

  public :: compute_reactive_eb_cfl_timestep_2d

contains

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
    if (.not. geometry%is_valid() .or. .not. ieee_is_finite(cfl) .or. &
        cfl <= 0.0_dp .or. cfl > 1.0_dp) return
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

end module reactive_eb_cfl_2d_mod
