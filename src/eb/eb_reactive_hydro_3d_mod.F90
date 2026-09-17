module eb_reactive_hydro_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive, &
    reactive_riemann_flux_x
  use reactive_2d_mod, only: reactive_riemann_flux_y
  use reactive_directional_flux_3d_mod, only: reactive_riemann_flux_z
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_wall_flux_3d_mod, only: &
    reactive_eb_flux_divergence_3d
  use eb_reactive_redistribution_3d_mod, only: &
    advance_reactive_eb_redistributed_3d, &
    advance_reactive_eb_state_redistributed_3d
  implicit none
  private

  public :: reactive_eb_outflow_riemann_fluxes_3d,reactive_eb_riemann_fluxes_3d
  public :: advance_reactive_eb_euler_3d
  public :: advance_reactive_eb_redistributed_euler_3d
  public :: advance_reactive_eb_state_redistributed_euler_3d

  interface reactive_eb_riemann_fluxes_3d
    module procedure reactive_eb_outflow_riemann_fluxes_3d
  end interface

contains

  subroutine reactive_eb_outflow_riemann_fluxes_3d( &
      species, state, temperature, geometry, solver, x_flux, y_flux, &
      z_flux, ok, ghost_state, ghost_temperature)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: x_flux(:, 0:, :, :)
    real(dp), intent(out) :: y_flux(:, :, 0:, :)
    real(dp), intent(out) :: z_flux(:, :, :, 0:)
    logical, intent(out) :: ok
    ! Optional face halo, indexed 0:n+1 in each spatial direction. Only open
    ! boundary-face donors are read; interior, edge and corner slots are unused.
    ! With no halo, retain the historical extrapolated-outflow boundary states.
    real(dp), optional, intent(in) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    real(dp), allocatable :: left_state(:),right_state(:)
    real(dp) :: left_temperature,right_temperature

    real(dp), allocatable :: candidate_x(:, :, :, :)
    real(dp), allocatable :: candidate_y(:, :, :, :)
    real(dp), allocatable :: candidate_z(:, :, :, :)
    logical :: face_ok
    integer :: i, j, k, lower_i, upper_i, lower_j, upper_j
    integer :: lower_k, upper_k, nvar

    x_flux = 0.0_dp
    y_flux = 0.0_dp
    z_flux = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (.not. valid_riemann_solver(solver)) return
    if (.not. valid_state_shapes(state, temperature, geometry, nvar) .or. &
        size(x_flux, 1) /= nvar .or. &
        size(x_flux, 2) /= geometry%nx + 1 .or. &
        size(x_flux, 3) /= geometry%ny .or. &
        size(x_flux, 4) /= geometry%nz .or. &
        size(y_flux, 1) /= nvar .or. &
        size(y_flux, 2) /= geometry%nx .or. &
        size(y_flux, 3) /= geometry%ny + 1 .or. &
        size(y_flux, 4) /= geometry%nz .or. &
        size(z_flux, 1) /= nvar .or. &
        size(z_flux, 2) /= geometry%nx .or. &
        size(z_flux, 3) /= geometry%ny .or. &
        size(z_flux, 4) /= geometry%nz + 1) return

    if (present(ghost_state).neqv.present(ghost_temperature)) return
    if (present(ghost_state)) then
      if (any(shape(ghost_state)/=[nvar,geometry%nx+2,geometry%ny+2,geometry%nz+2])) return
      if (any(shape(ghost_temperature)/=[geometry%nx+2,geometry%ny+2,geometry%nz+2])) return
    end if
    allocate(left_state(nvar),right_state(nvar))

    allocate(candidate_x(nvar, 0:geometry%nx, geometry%ny, geometry%nz))
    allocate(candidate_y(nvar, geometry%nx, 0:geometry%ny, geometry%nz))
    allocate(candidate_z(nvar, geometry%nx, geometry%ny, 0:geometry%nz))
    candidate_x = 0.0_dp
    candidate_y = 0.0_dp
    candidate_z = 0.0_dp

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 0, geometry%nx
          if (geometry%x_face_fraction(i, j, k) <= 0.0_dp) cycle
          lower_i = max(1, i)
          upper_i = min(geometry%nx, i + 1)
          if (geometry%cell_type(lower_i, j, k) == eb_covered_cell_3d &
              .or. geometry%cell_type(upper_i, j, k) == &
                eb_covered_cell_3d) return
          left_state=state(:,lower_i, j, k);right_state=state(:,upper_i, j, k)
          left_temperature=temperature(lower_i, j, k);right_temperature=temperature(upper_i, j, k)
          if (present(ghost_state)) then
            if (i==0) then
              left_state=ghost_state(:,0,j,k);left_temperature=ghost_temperature(0,j,k)
            end if
            if (i==geometry%nx) then
              right_state=ghost_state(:,geometry%nx+1,j,k);right_temperature=ghost_temperature(geometry%nx+1,j,k)
            end if
          end if
          call reactive_riemann_flux_x( &
            species, left_state, right_state, &
            left_temperature, right_temperature, &
            solver, candidate_x(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do

    do k = 1, geometry%nz
      do j = 0, geometry%ny
        lower_j = max(1, j)
        upper_j = min(geometry%ny, j + 1)
        do i = 1, geometry%nx
          if (geometry%y_face_fraction(i, j, k) <= 0.0_dp) cycle
          if (geometry%cell_type(i, lower_j, k) == eb_covered_cell_3d &
              .or. geometry%cell_type(i, upper_j, k) == &
                eb_covered_cell_3d) return
          left_state=state(:,i, lower_j, k);right_state=state(:,i, upper_j, k)
          left_temperature=temperature(i, lower_j, k);right_temperature=temperature(i, upper_j, k)
          if (present(ghost_state)) then
            if (j==0) then
              left_state=ghost_state(:,i,0,k);left_temperature=ghost_temperature(i,0,k)
            end if
            if (j==geometry%ny) then
              right_state=ghost_state(:,i,geometry%ny+1,k);right_temperature=ghost_temperature(i,geometry%ny+1,k)
            end if
          end if
          call reactive_riemann_flux_y( &
            species, left_state, right_state, &
            left_temperature, right_temperature, &
            solver, candidate_y(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do

    do k = 0, geometry%nz
      lower_k = max(1, k)
      upper_k = min(geometry%nz, k + 1)
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%z_face_fraction(i, j, k) <= 0.0_dp) cycle
          if (geometry%cell_type(i, j, lower_k) == eb_covered_cell_3d &
              .or. geometry%cell_type(i, j, upper_k) == &
                eb_covered_cell_3d) return
          left_state=state(:,i, j, lower_k);right_state=state(:,i, j, upper_k)
          left_temperature=temperature(i, j, lower_k);right_temperature=temperature(i, j, upper_k)
          if (present(ghost_state)) then
            if (k==0) then
              left_state=ghost_state(:,i,j,0);left_temperature=ghost_temperature(i,j,0)
            end if
            if (k==geometry%nz) then
              right_state=ghost_state(:,i,j,geometry%nz+1);right_temperature=ghost_temperature(i,j,geometry%nz+1)
            end if
          end if
          call reactive_riemann_flux_z( &
            species, left_state, right_state, &
            left_temperature, right_temperature, &
            solver, candidate_z(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_x)) .or. &
        any(.not. ieee_is_finite(candidate_y)) .or. &
        any(.not. ieee_is_finite(candidate_z))) return

    x_flux = candidate_x
    y_flux = candidate_y
    z_flux = candidate_z
    ok = .true.
  end subroutine reactive_eb_outflow_riemann_fluxes_3d

  subroutine build_reactive_eb_conservative_rhs_3d( &
      species, state, temperature, geometry, solver, rhs, ok, ghost_state, ghost_temperature, flux)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: rhs(:, :, :, :)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    type(eb_integrated_flux_3d), optional, intent(inout) :: flux

    real(dp), allocatable :: x_flux(:, :, :, :), y_flux(:, :, :, :)
    real(dp), allocatable :: z_flux(:, :, :, :)
    logical :: local_ok
    integer :: nvar

    rhs = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (.not. valid_state_shapes(state, temperature, geometry, nvar) .or. &
        size(rhs, 1) /= nvar .or. &
        size(rhs, 2) /= geometry%nx .or. &
        size(rhs, 3) /= geometry%ny .or. &
        size(rhs, 4) /= geometry%nz) return

    allocate(x_flux(nvar, 0:geometry%nx, geometry%ny, geometry%nz))
    allocate(y_flux(nvar, geometry%nx, 0:geometry%ny, geometry%nz))
    allocate(z_flux(nvar, geometry%nx, geometry%ny, 0:geometry%nz))
    call reactive_eb_outflow_riemann_fluxes_3d( &
      species, state, temperature, geometry, solver, x_flux, y_flux, &
      z_flux, local_ok, ghost_state, ghost_temperature)
    if (.not. local_ok) return
    call reactive_eb_flux_divergence_3d( &
      species, state, temperature, geometry, x_flux, y_flux, z_flux, &
      rhs, local_ok)
    if (.not. local_ok) then
      rhs = 0.0_dp
      return
    end if
    if (present(flux)) then
      flux%x=x_flux;flux%y=y_flux;flux%z=z_flux
    end if
    ok = .true.
  end subroutine build_reactive_eb_conservative_rhs_3d

  subroutine advance_reactive_eb_euler_3d( &
      species, state, temperature, geometry, solver, dt, new_state, &
      new_temperature, ok, ghost_state, ghost_temperature, integrated_flux)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: new_state(:, :, :, :)
    real(dp), intent(out) :: new_temperature(:, :, :)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    type(eb_integrated_flux_3d), optional, intent(inout) :: integrated_flux
    type(eb_integrated_flux_3d) :: flux

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: rhs(:, :, :, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, k, nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (.not. valid_state_shapes(state, temperature, geometry, nvar) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) return

    allocate(rhs(nvar, geometry%nx, geometry%ny, geometry%nz))
    if (present(integrated_flux)) then
      call build_reactive_eb_conservative_rhs_3d( &
        species,state,temperature,geometry,solver,rhs,local_ok,ghost_state,ghost_temperature,flux)
    else
      call build_reactive_eb_conservative_rhs_3d( &
        species,state,temperature,geometry,solver,rhs,local_ok,ghost_state,ghost_temperature)
    end if
    if (.not. local_ok) return

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          candidate_state(:, i, j, k) = &
            state(:, i, j, k) + dt * rhs(:, i, j, k)
          call reactive_conserved_to_primitive( &
            species, candidate_state(:, i, j, k), temperature(i, j, k), &
            primitive, recovered_temperature, sound_speed, local_ok)
          if (.not. local_ok) return
          candidate_temperature(i, j, k) = recovered_temperature
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_state)) .or. &
        any(.not. ieee_is_finite(candidate_temperature))) return

    if (present(integrated_flux)) then
      flux%x=dt*flux%x;flux%y=dt*flux%y;flux%z=dt*flux%z
      if (.not.all(ieee_is_finite(flux%x)).or..not.all(ieee_is_finite(flux%y)).or. &
        .not.all(ieee_is_finite(flux%z))) return
      integrated_flux=flux
    end if
    new_state = candidate_state
    new_temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_eb_euler_3d

  subroutine advance_reactive_eb_redistributed_euler_3d( &
      species, state, temperature, geometry, solver, dt, new_state, &
      new_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: new_state(:, :, :, :)
    real(dp), intent(out) :: new_temperature(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: rhs(:, :, :, :)
    logical :: local_ok
    integer :: nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (.not. valid_state_shapes(state, temperature, geometry, nvar) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) return

    allocate(rhs(nvar, geometry%nx, geometry%ny, geometry%nz))
    call build_reactive_eb_conservative_rhs_3d( &
      species, state, temperature, geometry, solver, rhs, local_ok)
    if (.not. local_ok) return

    call advance_reactive_eb_redistributed_3d( &
      species, state, temperature, geometry, rhs, dt, new_state, &
      new_temperature, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_reactive_eb_redistributed_euler_3d

  subroutine advance_reactive_eb_state_redistributed_euler_3d( &
      species, state, temperature, geometry, solver, dt, new_state, &
      new_temperature, ok, target_volume_fraction)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: new_state(:, :, :, :)
    real(dp), intent(out) :: new_temperature(:, :, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction

    real(dp), allocatable :: rhs(:, :, :, :)
    real(dp) :: target
    logical :: local_ok
    integer :: nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (.not. valid_state_shapes(state, temperature, geometry, nvar) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) return

    allocate(rhs(nvar, geometry%nx, geometry%ny, geometry%nz))
    call build_reactive_eb_conservative_rhs_3d( &
      species, state, temperature, geometry, solver, rhs, local_ok)
    if (.not. local_ok) return

    target = 0.5_dp
    if (present(target_volume_fraction)) target = target_volume_fraction
    call advance_reactive_eb_state_redistributed_3d( &
      species, state, temperature, geometry, rhs, dt, new_state, &
      new_temperature, local_ok, target)
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_reactive_eb_state_redistributed_euler_3d

  pure logical function valid_state_shapes( &
      state, temperature, geometry, nvar) result(valid)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: nvar

    valid = size(state, 1) == nvar .and. &
      size(state, 2) == geometry%nx .and. &
      size(state, 3) == geometry%ny .and. &
      size(state, 4) == geometry%nz .and. &
      all(shape(temperature) == &
        [geometry%nx, geometry%ny, geometry%nz])
  end function valid_state_shapes

  pure logical function valid_riemann_solver(solver) result(valid)
    character(len=*), intent(in) :: solver

    select case (trim(solver))
    case ("rusanov", "hllc", "pelec")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_riemann_solver

end module eb_reactive_hydro_3d_mod
