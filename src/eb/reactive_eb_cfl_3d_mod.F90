module reactive_eb_cfl_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use eb_geometry_3d_mod, only: eb_geometry_3d, eb_covered_cell_3d
  implicit none
  private

  public :: compute_reactive_eb_cfl_timestep_3d

contains

  subroutine compute_reactive_eb_cfl_timestep_3d( &
      species, state, temperature, geometry, cfl, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed, rate, maximum_rate
    logical :: local_ok
    integer :: i, j, k, active_cells, nvar

    dt = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid() .or. &
        .not. ieee_is_finite(cfl) .or. cfl <= 0.0_dp .or. &
        cfl > 1.0_dp) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz])) return

    allocate(primitive(reactive_nprim(size(species))))
    maximum_rate = 0.0_dp
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
          rate = (abs(primitive(2)) + sound_speed) / geometry%dx + &
            (abs(primitive(3)) + sound_speed) / geometry%dy + &
            (abs(primitive(4)) + sound_speed) / geometry%dz
          maximum_rate = max(maximum_rate, rate)
        end do
      end do
    end do
    if (active_cells == 0 .or. .not. ieee_is_finite(maximum_rate) .or. &
        maximum_rate <= 0.0_dp) return

    dt = cfl / maximum_rate
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp
    if (.not. ok) dt = 0.0_dp
  end subroutine compute_reactive_eb_cfl_timestep_3d

end module reactive_eb_cfl_3d_mod
