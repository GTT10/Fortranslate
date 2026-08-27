module amr_eb_reactive_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use slope_limiter_mod, only: limited_slope
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, eb_regular_cell
  use eb_reactive_reconstruction_2d_mod, only: &
    reactive_eb_exterior_state_2d, &
    build_reactive_eb_face_center_fluxes_2d, &
    interpolate_reactive_eb_face_centroid_fluxes_2d
  use eb_reactive_wall_flux_2d_mod, only: reactive_eb_flux_divergence_2d
  use eb_reactive_redistribution_2d_mod, only: &
    advance_reactive_eb_state_redistributed_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, average_down_reactive_eb_state_patch_2d
  use amr_eb_flux_register_2d_mod, only: &
    amr_eb_flux_register_2d, initialize_amr_eb_flux_register_2d, &
    accumulate_coarse_eb_fluxes_2d, accumulate_fine_eb_fluxes_2d, &
    reflux_reactive_eb_state_patch_2d
  implicit none
  private

  type, public :: reactive_eb_patch_exterior_context_2d
    type(reactive_eb_exterior_state_2d) :: start
    type(reactive_eb_exterior_state_2d) :: end
  contains
    procedure :: is_valid => reactive_eb_patch_exterior_context_is_valid
  end type reactive_eb_patch_exterior_context_2d

  public :: prolong_reactive_eb_patch_pcm_2d
  public :: prolong_reactive_eb_patch_linear_2d
  public :: prolong_reactive_eb_patch_2d
  public :: extract_reactive_eb_patch_exterior_context_support_2d
  public :: extract_reactive_eb_patch_exterior_context_2d
  public :: build_reactive_eb_patch_exterior_from_context_2d
  public :: build_reactive_eb_patch_exterior_2d
  public :: advance_reactive_eb_level_2d
  public :: advance_two_level_reactive_eb_hydro_2d

contains

  subroutine prolong_reactive_eb_patch_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      fine_geometry, patch, method, fine_state, fine_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    character(len=*), intent(in) :: method
    real(dp), intent(out) :: fine_state(:, :, :), fine_temperature(:, :)
    logical, intent(out) :: ok

    select case (trim(method))
    case ("pcm")
      call prolong_reactive_eb_patch_pcm_2d( &
        species, coarse_state, coarse_temperature, coarse_geometry, &
        fine_geometry, patch, fine_state, fine_temperature, ok)
    case ("linear")
      call prolong_reactive_eb_patch_linear_2d( &
        species, coarse_state, coarse_temperature, coarse_geometry, &
        fine_geometry, patch, fine_state, fine_temperature, ok)
    case default
      fine_state = 0.0_dp
      fine_temperature = 0.0_dp
      ok = .false.
    end select
  end subroutine prolong_reactive_eb_patch_2d

  pure logical function reactive_eb_patch_exterior_context_is_valid( &
      self, fine_geometry, component_count) result(valid)
    class(reactive_eb_patch_exterior_context_2d), intent(in) :: self
    type(eb_geometry_2d), intent(in) :: fine_geometry
    integer, intent(in) :: component_count

    valid = self%start%is_valid(fine_geometry, component_count) .and. &
      self%end%is_valid(fine_geometry, component_count)
  end function reactive_eb_patch_exterior_context_is_valid

  subroutine prolong_reactive_eb_patch_pcm_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      fine_geometry, patch, fine_state, fine_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    real(dp), intent(out) :: fine_state(:, :, :), fine_temperature(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :), primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: coarse_i, coarse_j, fine_i, fine_j, nvar, ratio

    fine_state = 0.0_dp
    fine_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. size(coarse_state, 1) /= nvar .or. &
        size(coarse_state, 2) /= coarse_geometry%nx .or. &
        size(coarse_state, 3) /= coarse_geometry%ny .or. &
        any(shape(coarse_temperature) /= &
          [coarse_geometry%nx, coarse_geometry%ny]) .or. &
        any(shape(fine_state) /= &
          [nvar, fine_geometry%nx, fine_geometry%ny]) .or. &
        any(shape(fine_temperature) /= &
          [fine_geometry%nx, fine_geometry%ny]) .or. &
        .not. patch%is_valid(coarse_geometry, fine_geometry) .or. &
        any(.not. ieee_is_finite(coarse_state)) .or. &
        any(.not. ieee_is_finite(coarse_temperature))) return

    allocate(candidate_state, mold=fine_state)
    allocate(candidate_temperature, mold=fine_temperature)
    allocate(primitive(reactive_nprim(size(species))))
    candidate_state = 0.0_dp
    candidate_temperature = 0.0_dp
    ratio = patch%refinement_ratio
    do fine_j = 1, fine_geometry%ny
      coarse_j = patch%coarse_j_lower + (fine_j - 1) / ratio
      do fine_i = 1, fine_geometry%nx
        coarse_i = patch%coarse_i_lower + (fine_i - 1) / ratio
        candidate_state(:, fine_i, fine_j) = &
          coarse_state(:, coarse_i, coarse_j)
        candidate_temperature(fine_i, fine_j) = &
          coarse_temperature(coarse_i, coarse_j)
        if (fine_geometry%cell_type(fine_i, fine_j) == eb_covered_cell) cycle
        if (candidate_temperature(fine_i, fine_j) <= 0.0_dp) return
        call reactive_conserved_to_primitive( &
          species, candidate_state(:, fine_i, fine_j), &
          candidate_temperature(fine_i, fine_j), primitive, &
          recovered_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        candidate_temperature(fine_i, fine_j) = recovered_temperature
      end do
    end do
    fine_state = candidate_state
    fine_temperature = candidate_temperature
    ok = .true.
  end subroutine prolong_reactive_eb_patch_pcm_2d

  pure logical function cut_parent_neighbor_connected_2d( &
      geometry, parent_i, parent_j, neighbor_i, neighbor_j) result(connected)
    type(eb_geometry_2d), intent(in) :: geometry
    integer, intent(in) :: parent_i, parent_j, neighbor_i, neighbor_j

    integer :: delta_i, delta_j, x_face_i, y_face_j
    logical :: horizontal_path, vertical_path

    connected = .false.
    if (parent_i < 1 .or. parent_i > geometry%nx .or. &
        parent_j < 1 .or. parent_j > geometry%ny .or. &
        neighbor_i < 1 .or. neighbor_i > geometry%nx .or. &
        neighbor_j < 1 .or. neighbor_j > geometry%ny .or. &
        geometry%cell_type(parent_i, parent_j) == eb_covered_cell .or. &
        geometry%cell_type(neighbor_i, neighbor_j) == eb_covered_cell) return
    delta_i = neighbor_i - parent_i
    delta_j = neighbor_j - parent_j
    if (abs(delta_i) > 1 .or. abs(delta_j) > 1 .or. &
        (delta_i == 0 .and. delta_j == 0)) return
    x_face_i = min(parent_i, neighbor_i)
    y_face_j = min(parent_j, neighbor_j)
    if (delta_j == 0) then
      connected = geometry%x_face_fraction(x_face_i, parent_j) > 0.0_dp
      return
    end if
    if (delta_i == 0) then
      connected = geometry%y_face_fraction(parent_i, y_face_j) > 0.0_dp
      return
    end if

    horizontal_path = &
      geometry%cell_type(neighbor_i, parent_j) /= eb_covered_cell .and. &
      geometry%x_face_fraction(x_face_i, parent_j) > 0.0_dp .and. &
      geometry%y_face_fraction(neighbor_i, y_face_j) > 0.0_dp
    vertical_path = &
      geometry%cell_type(parent_i, neighbor_j) /= eb_covered_cell .and. &
      geometry%y_face_fraction(parent_i, y_face_j) > 0.0_dp .and. &
      geometry%x_face_fraction(x_face_i, neighbor_j) > 0.0_dp
    connected = horizontal_path .or. vertical_path
  end function cut_parent_neighbor_connected_2d

  pure logical function cut_parent_neighbor_connected_within_radius_2d( &
      geometry, parent_i, parent_j, neighbor_i, neighbor_j, radius) &
      result(connected)
    type(eb_geometry_2d), intent(in) :: geometry
    integer, intent(in) :: parent_i, parent_j, neighbor_i, neighbor_j, radius

    integer, parameter :: maximum_radius = 2
    integer, parameter :: maximum_cells = &
      (2 * maximum_radius + 1) * (2 * maximum_radius + 1)
    integer, parameter :: direction_i(4) = [-1, 1, 0, 0]
    integer, parameter :: direction_j(4) = [0, 0, -1, 1]
    integer :: current_i, current_j, direction, head, next_i, next_j, tail
    integer :: offset_i, offset_j, queue_i(maximum_cells)
    integer :: queue_j(maximum_cells)
    logical :: face_open
    logical :: visited(-maximum_radius:maximum_radius, &
      -maximum_radius:maximum_radius)

    connected = .false.
    if (radius == 1) then
      connected = cut_parent_neighbor_connected_2d( &
        geometry, parent_i, parent_j, neighbor_i, neighbor_j)
      return
    end if
    if (radius < 1 .or. radius > maximum_radius) return
    if (parent_i < 1 .or. parent_i > geometry%nx .or. &
        parent_j < 1 .or. parent_j > geometry%ny .or. &
        neighbor_i < 1 .or. neighbor_i > geometry%nx .or. &
        neighbor_j < 1 .or. neighbor_j > geometry%ny) return
    if (geometry%cell_type(parent_i, parent_j) == eb_covered_cell .or. &
        geometry%cell_type(neighbor_i, neighbor_j) == eb_covered_cell .or. &
        (parent_i == neighbor_i .and. parent_j == neighbor_j) .or. &
        abs(neighbor_i - parent_i) > radius .or. &
        abs(neighbor_j - parent_j) > radius) return

    visited = .false.
    visited(0, 0) = .true.
    head = 1
    tail = 1
    queue_i(1) = parent_i
    queue_j(1) = parent_j
    do while (head <= tail)
      current_i = queue_i(head)
      current_j = queue_j(head)
      head = head + 1
      do direction = 1, size(direction_i)
        next_i = current_i + direction_i(direction)
        next_j = current_j + direction_j(direction)
        offset_i = next_i - parent_i
        offset_j = next_j - parent_j
        if (abs(offset_i) > radius .or. abs(offset_j) > radius) cycle
        if (next_i < 1 .or. next_i > geometry%nx .or. &
            next_j < 1 .or. next_j > geometry%ny) cycle
        if (visited(offset_i, offset_j)) cycle
        if (geometry%cell_type(next_i, next_j) == eb_covered_cell) cycle
        if (direction_i(direction) /= 0) then
          face_open = geometry%x_face_fraction( &
            min(current_i, next_i), current_j) > 0.0_dp
        else
          face_open = geometry%y_face_fraction( &
            current_i, min(current_j, next_j)) > 0.0_dp
        end if
        if (.not. face_open) cycle
        if (next_i == neighbor_i .and. next_j == neighbor_j) then
          connected = .true.
          return
        end if
        visited(offset_i, offset_j) = .true.
        tail = tail + 1
        queue_i(tail) = next_i
        queue_j(tail) = next_j
      end do
    end do
  end function cut_parent_neighbor_connected_within_radius_2d

  subroutine build_cut_parent_limited_slopes_2d( &
      geometry, state, parent_i, parent_j, slope_x, slope_y, &
      minimum_state, maximum_state, ok)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: state(:, :, :)
    integer, intent(in) :: parent_i, parent_j
    real(dp), intent(out) :: slope_x(:), slope_y(:)
    real(dp), intent(out) :: minimum_state(:), maximum_state(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: normal_rhs_x(:), normal_rhs_y(:), limiter(:)
    real(dp) :: delta_x, delta_y, determinant, determinant_scale
    real(dp) :: normal_xx, normal_xy, normal_yy, normal_trace
    real(dp) :: predicted_delta, rank_tolerance
    integer :: component, neighbor_i, neighbor_j, ncomp, stencil_radius

    slope_x = 0.0_dp
    slope_y = 0.0_dp
    minimum_state = 0.0_dp
    maximum_state = 0.0_dp
    ok = .false.
    ncomp = size(state, 1)
    if (ncomp < 1 .or. size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(slope_x) /= ncomp .or. size(slope_y) /= ncomp .or. &
        size(minimum_state) /= ncomp .or. &
        size(maximum_state) /= ncomp .or. &
        parent_i < 1 .or. parent_i > geometry%nx .or. &
        parent_j < 1 .or. parent_j > geometry%ny .or. &
        geometry%cell_type(parent_i, parent_j) /= eb_cut_cell) return

    allocate(normal_rhs_x(ncomp), normal_rhs_y(ncomp), limiter(ncomp))
    rank_tolerance = 4096.0_dp * epsilon(1.0_dp)
    stencil_radius = 1
    call build_normal_system(stencil_radius)
    determinant = normal_xx * normal_yy - normal_xy * normal_xy
    determinant_scale = max(1.0_dp, normal_xx * normal_yy)
    if (abs(determinant) <= rank_tolerance * determinant_scale) then
      stencil_radius = 2
      call build_normal_system(stencil_radius)
      determinant = normal_xx * normal_yy - normal_xy * normal_xy
      determinant_scale = max(1.0_dp, normal_xx * normal_yy)
    end if
    normal_trace = normal_xx + normal_yy
    if (abs(determinant) > rank_tolerance * determinant_scale) then
      slope_x = (normal_rhs_x * normal_yy - &
        normal_xy * normal_rhs_y) / determinant
      slope_y = (normal_xx * normal_rhs_y - &
        normal_xy * normal_rhs_x) / determinant
    else if (normal_trace > rank_tolerance) then
      ! The minimum-norm rank-one fit retains variation along the only
      ! connected fluid-centroid direction without inventing a normal slope.
      slope_x = normal_rhs_x / normal_trace
      slope_y = normal_rhs_y / normal_trace
    end if

    limiter = 1.0_dp
    do neighbor_j = max(1, parent_j - stencil_radius), &
        min(geometry%ny, parent_j + stencil_radius)
      do neighbor_i = max(1, parent_i - stencil_radius), &
          min(geometry%nx, parent_i + stencil_radius)
        if (.not. cut_parent_neighbor_connected_within_radius_2d( &
            geometry, parent_i, parent_j, neighbor_i, neighbor_j, &
            stencil_radius)) cycle
        delta_x = real(neighbor_i - parent_i, dp) + &
          geometry%cell_centroid_x(neighbor_i, neighbor_j) - &
          geometry%cell_centroid_x(parent_i, parent_j)
        delta_y = real(neighbor_j - parent_j, dp) + &
          geometry%cell_centroid_y(neighbor_i, neighbor_j) - &
          geometry%cell_centroid_y(parent_i, parent_j)
        do component = 1, ncomp
          predicted_delta = delta_x * slope_x(component) + &
            delta_y * slope_y(component)
          if (predicted_delta > 0.0_dp) then
            limiter(component) = min(limiter(component), &
              (maximum_state(component) - &
               state(component, parent_i, parent_j)) / predicted_delta)
          else if (predicted_delta < 0.0_dp) then
            limiter(component) = min(limiter(component), &
              (minimum_state(component) - &
               state(component, parent_i, parent_j)) / predicted_delta)
          end if
        end do
      end do
    end do
    limiter = max(0.0_dp, min(1.0_dp, limiter))
    slope_x = limiter * slope_x
    slope_y = limiter * slope_y
    ok = all(ieee_is_finite(slope_x)) .and. &
      all(ieee_is_finite(slope_y))

  contains

    subroutine build_normal_system(radius)
      integer, intent(in) :: radius

      normal_rhs_x = 0.0_dp
      normal_rhs_y = 0.0_dp
      normal_xx = 0.0_dp
      normal_xy = 0.0_dp
      normal_yy = 0.0_dp
      minimum_state = state(:, parent_i, parent_j)
      maximum_state = minimum_state
      do neighbor_j = max(1, parent_j - radius), &
          min(geometry%ny, parent_j + radius)
        do neighbor_i = max(1, parent_i - radius), &
            min(geometry%nx, parent_i + radius)
          if (.not. cut_parent_neighbor_connected_within_radius_2d( &
              geometry, parent_i, parent_j, neighbor_i, neighbor_j, &
              radius)) cycle
          minimum_state = min( &
            minimum_state, state(:, neighbor_i, neighbor_j))
          maximum_state = max( &
            maximum_state, state(:, neighbor_i, neighbor_j))
          delta_x = real(neighbor_i - parent_i, dp) + &
            geometry%cell_centroid_x(neighbor_i, neighbor_j) - &
            geometry%cell_centroid_x(parent_i, parent_j)
          delta_y = real(neighbor_j - parent_j, dp) + &
            geometry%cell_centroid_y(neighbor_i, neighbor_j) - &
            geometry%cell_centroid_y(parent_i, parent_j)
          normal_xx = normal_xx + delta_x * delta_x
          normal_xy = normal_xy + delta_x * delta_y
          normal_yy = normal_yy + delta_y * delta_y
          normal_rhs_x = normal_rhs_x + delta_x * &
            (state(:, neighbor_i, neighbor_j) - &
             state(:, parent_i, parent_j))
          normal_rhs_y = normal_rhs_y + delta_y * &
            (state(:, neighbor_i, neighbor_j) - &
             state(:, parent_i, parent_j))
        end do
      end do
    end subroutine build_normal_system

  end subroutine build_cut_parent_limited_slopes_2d

  subroutine prolong_reactive_eb_patch_linear_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      fine_geometry, patch, fine_state, fine_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    real(dp), intent(out) :: fine_state(:, :, :), fine_temperature(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :), primitive(:)
    real(dp), allocatable :: maximum_state(:), minimum_state(:)
    real(dp), allocatable :: slope_x(:), slope_y(:), theta(:)
    real(dp) :: delta, delta_minus, delta_plus
    real(dp) :: mean_offset_x, mean_offset_y, offset_x, offset_y
    real(dp) :: raw_offset_x, raw_offset_y, total_weight, weight
    real(dp) :: recovered_temperature, sound_speed
    logical :: is_cut_parent, local_ok, parent_ok, use_linear
    integer :: child_i, child_j, coarse_i, coarse_j, component
    integer :: fine_i, fine_i_lower, fine_i_upper
    integer :: fine_j, fine_j_lower, fine_j_upper, nvar, ratio

    fine_state = 0.0_dp
    fine_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. size(coarse_state, 1) /= nvar .or. &
        size(coarse_state, 2) /= coarse_geometry%nx .or. &
        size(coarse_state, 3) /= coarse_geometry%ny .or. &
        any(shape(coarse_temperature) /= &
          [coarse_geometry%nx, coarse_geometry%ny]) .or. &
        any(shape(fine_state) /= &
          [nvar, fine_geometry%nx, fine_geometry%ny]) .or. &
        any(shape(fine_temperature) /= &
          [fine_geometry%nx, fine_geometry%ny]) .or. &
        .not. patch%is_valid(coarse_geometry, fine_geometry) .or. &
        any(.not. ieee_is_finite(coarse_state)) .or. &
        any(.not. ieee_is_finite(coarse_temperature))) return

    allocate(candidate_state, mold=fine_state)
    allocate(candidate_temperature, mold=fine_temperature)
    allocate(primitive(reactive_nprim(size(species))))
    allocate(maximum_state(nvar), minimum_state(nvar))
    allocate(slope_x(nvar), slope_y(nvar), theta(nvar))
    candidate_state = 0.0_dp
    candidate_temperature = 0.0_dp
    ratio = patch%refinement_ratio
    do coarse_j = patch%coarse_j_lower, patch%coarse_j_upper
      fine_j_lower = (coarse_j - patch%coarse_j_lower) * ratio + 1
      fine_j_upper = fine_j_lower + ratio - 1
      do coarse_i = patch%coarse_i_lower, patch%coarse_i_upper
        fine_i_lower = (coarse_i - patch%coarse_i_lower) * ratio + 1
        fine_i_upper = fine_i_lower + ratio - 1
        is_cut_parent = &
          coarse_geometry%cell_type(coarse_i, coarse_j) == eb_cut_cell
        use_linear = &
          coarse_geometry%cell_type(coarse_i, coarse_j) == &
            eb_regular_cell .and. &
          all(fine_geometry%cell_type( &
            fine_i_lower:fine_i_upper, fine_j_lower:fine_j_upper) == &
            eb_regular_cell)
        slope_x = 0.0_dp
        slope_y = 0.0_dp
        mean_offset_x = 0.0_dp
        mean_offset_y = 0.0_dp
        if (use_linear) then
          do component = 1, nvar
            delta_minus = 0.0_dp
            delta_plus = 0.0_dp
            if (coarse_i > 1 .and. &
                coarse_geometry%cell_type(coarse_i - 1, coarse_j) == &
                  eb_regular_cell) then
              delta_minus = coarse_state(component, coarse_i, coarse_j) - &
                coarse_state(component, coarse_i - 1, coarse_j)
            end if
            if (coarse_i < coarse_geometry%nx .and. &
                coarse_geometry%cell_type(coarse_i + 1, coarse_j) == &
                  eb_regular_cell) then
              delta_plus = coarse_state(component, coarse_i + 1, coarse_j) - &
                coarse_state(component, coarse_i, coarse_j)
            end if
            call limited_slope( &
              delta_minus, delta_plus, "mc", slope_x(component), local_ok)
            if (.not. local_ok) return

            delta_minus = 0.0_dp
            delta_plus = 0.0_dp
            if (coarse_j > 1 .and. &
                coarse_geometry%cell_type(coarse_i, coarse_j - 1) == &
                  eb_regular_cell) then
              delta_minus = coarse_state(component, coarse_i, coarse_j) - &
                coarse_state(component, coarse_i, coarse_j - 1)
            end if
            if (coarse_j < coarse_geometry%ny .and. &
                coarse_geometry%cell_type(coarse_i, coarse_j + 1) == &
                  eb_regular_cell) then
              delta_plus = coarse_state(component, coarse_i, coarse_j + 1) - &
                coarse_state(component, coarse_i, coarse_j)
            end if
            call limited_slope( &
              delta_minus, delta_plus, "mc", slope_y(component), local_ok)
            if (.not. local_ok) return
          end do
        else if (is_cut_parent) then
          ! Fit one multidimensional gradient to the connected active coarse
          ! fluid centroids. Rank-one support retains its resolved tangent;
          ! the coarse and child envelope limiters prevent new component
          ! bounds without perturbing an exactly affine fit.
          use_linear = .true.
          call build_cut_parent_limited_slopes_2d( &
            coarse_geometry, coarse_state, coarse_i, coarse_j, slope_x, &
            slope_y, minimum_state, maximum_state, local_ok)
          if (.not. local_ok) return

          ! Fine fluid-centroid offsets need not have zero mean in a cut
          ! parent. Remove their volume-weighted mean before reconstruction so
          ! average-down returns the parent state exactly.
          total_weight = 0.0_dp
          do child_j = 1, ratio
            fine_j = fine_j_lower + child_j - 1
            do child_i = 1, ratio
              fine_i = fine_i_lower + child_i - 1
              weight = fine_geometry%volume_fraction(fine_i, fine_j)
              if (weight <= 0.0_dp) cycle
              raw_offset_x = ( &
                real(child_i, dp) - 0.5_dp + &
                fine_geometry%cell_centroid_x(fine_i, fine_j)) / &
                real(ratio, dp) - 0.5_dp
              raw_offset_y = ( &
                real(child_j, dp) - 0.5_dp + &
                fine_geometry%cell_centroid_y(fine_i, fine_j)) / &
                real(ratio, dp) - 0.5_dp
              total_weight = total_weight + weight
              mean_offset_x = mean_offset_x + weight * raw_offset_x
              mean_offset_y = mean_offset_y + weight * raw_offset_y
            end do
          end do
          if (total_weight <= 0.0_dp) then
            use_linear = .false.
          else
            mean_offset_x = mean_offset_x / total_weight
            mean_offset_y = mean_offset_y / total_weight
            theta = 1.0_dp
            do child_j = 1, ratio
              fine_j = fine_j_lower + child_j - 1
              do child_i = 1, ratio
                fine_i = fine_i_lower + child_i - 1
                if (fine_geometry%cell_type(fine_i, fine_j) == &
                    eb_covered_cell) cycle
                offset_x = ( &
                  real(child_i, dp) - 0.5_dp + &
                  fine_geometry%cell_centroid_x(fine_i, fine_j)) / &
                  real(ratio, dp) - 0.5_dp - mean_offset_x
                offset_y = ( &
                  real(child_j, dp) - 0.5_dp + &
                  fine_geometry%cell_centroid_y(fine_i, fine_j)) / &
                  real(ratio, dp) - 0.5_dp - mean_offset_y
                do component = 1, nvar
                  delta = offset_x * slope_x(component) + &
                    offset_y * slope_y(component)
                  if (delta > 0.0_dp) then
                    theta(component) = min(theta(component), &
                      (maximum_state(component) - &
                       coarse_state(component, coarse_i, coarse_j)) / delta)
                  else if (delta < 0.0_dp) then
                    theta(component) = min(theta(component), &
                      (minimum_state(component) - &
                       coarse_state(component, coarse_i, coarse_j)) / delta)
                  end if
                end do
              end do
            end do
            theta = max(0.0_dp, min(1.0_dp, theta))
            slope_x = theta * slope_x
            slope_y = theta * slope_y
          end if
        end if

        parent_ok = .true.
        do child_j = 1, ratio
          fine_j = fine_j_lower + child_j - 1
          do child_i = 1, ratio
            fine_i = fine_i_lower + child_i - 1
            if (is_cut_parent .and. use_linear) then
              offset_x = ( &
                real(child_i, dp) - 0.5_dp + &
                fine_geometry%cell_centroid_x(fine_i, fine_j)) / &
                real(ratio, dp) - 0.5_dp - mean_offset_x
              offset_y = ( &
                real(child_j, dp) - 0.5_dp + &
                fine_geometry%cell_centroid_y(fine_i, fine_j)) / &
                real(ratio, dp) - 0.5_dp - mean_offset_y
            else
              offset_x = (real(child_i, dp) - 0.5_dp) / &
                real(ratio, dp) - 0.5_dp
              offset_y = (real(child_j, dp) - 0.5_dp) / &
                real(ratio, dp) - 0.5_dp
            end if
            candidate_state(:, fine_i, fine_j) = &
              coarse_state(:, coarse_i, coarse_j) + &
              offset_x * slope_x + offset_y * slope_y
            candidate_temperature(fine_i, fine_j) = &
              coarse_temperature(coarse_i, coarse_j)
            if (fine_geometry%cell_type(fine_i, fine_j) == &
                eb_covered_cell) cycle
            if (candidate_temperature(fine_i, fine_j) <= 0.0_dp) then
              parent_ok = .false.
              exit
            end if
            call reactive_conserved_to_primitive( &
              species, candidate_state(:, fine_i, fine_j), &
              candidate_temperature(fine_i, fine_j), primitive, &
              recovered_temperature, sound_speed, local_ok)
            if (.not. local_ok) then
              parent_ok = .false.
              exit
            end if
            candidate_temperature(fine_i, fine_j) = recovered_temperature
          end do
          if (.not. parent_ok) exit
        end do
        if (.not. parent_ok .and. use_linear) then
          ! A component-wise conservative slope may still leave the EOS
          ! admissible set; retry this parent with the qualified PCM state.
          do fine_j = fine_j_lower, fine_j_upper
            do fine_i = fine_i_lower, fine_i_upper
              candidate_state(:, fine_i, fine_j) = &
                coarse_state(:, coarse_i, coarse_j)
              candidate_temperature(fine_i, fine_j) = &
                coarse_temperature(coarse_i, coarse_j)
              call reactive_conserved_to_primitive( &
                species, candidate_state(:, fine_i, fine_j), &
                candidate_temperature(fine_i, fine_j), primitive, &
                recovered_temperature, sound_speed, local_ok)
              if (.not. local_ok) return
              candidate_temperature(fine_i, fine_j) = recovered_temperature
            end do
          end do
          parent_ok = .true.
        end if
        if (.not. parent_ok) return
      end do
    end do
    fine_state = candidate_state
    fine_temperature = candidate_temperature
    ok = .true.
  end subroutine prolong_reactive_eb_patch_linear_2d

  subroutine recover_exterior_cell( &
      species, start_state, end_state, start_temperature, end_temperature, &
      alpha, state, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: start_state(:), end_state(:)
    real(dp), intent(in) :: start_temperature, end_temperature, alpha
    real(dp), intent(out) :: state(:), temperature
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: sound_speed, temperature_guess

    state = (1.0_dp - alpha) * start_state + alpha * end_state
    temperature = 0.0_dp
    ok = .false.
    temperature_guess = (1.0_dp - alpha) * start_temperature + &
      alpha * end_temperature
    if (temperature_guess <= 0.0_dp .or. &
        any(.not. ieee_is_finite(state))) return
    allocate(primitive(reactive_nprim(size(species))))
    call reactive_conserved_to_primitive( &
      species, state, temperature_guess, primitive, temperature, &
      sound_speed, ok)
  end subroutine recover_exterior_cell

  subroutine extract_reactive_eb_patch_exterior_context_support_2d( &
      coarse_i_lower, coarse_j_lower, coarse_start, &
      coarse_start_temperature, coarse_end, coarse_end_temperature, &
      coarse_geometry, fine_geometry, patch, component_count, context, ok)
    integer, intent(in) :: coarse_i_lower, coarse_j_lower
    real(dp), intent(in) :: &
      coarse_start(:, coarse_i_lower:, coarse_j_lower:)
    real(dp), intent(in) :: &
      coarse_start_temperature(coarse_i_lower:, coarse_j_lower:)
    real(dp), intent(in) :: &
      coarse_end(:, coarse_i_lower:, coarse_j_lower:)
    real(dp), intent(in) :: &
      coarse_end_temperature(coarse_i_lower:, coarse_j_lower:)
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    integer, intent(in) :: component_count
    type(reactive_eb_patch_exterior_context_2d), intent(out) :: context
    logical, intent(out) :: ok

    integer :: coarse_i, coarse_i_upper, coarse_j, coarse_j_upper
    integer :: expected_i_lower, expected_i_upper
    integer :: expected_j_lower, expected_j_upper, fine_i, fine_j, ratio

    ok = .false.
    coarse_i_upper = coarse_i_lower + size(coarse_start, 2) - 1
    coarse_j_upper = coarse_j_lower + size(coarse_start, 3) - 1
    expected_i_lower = max(1, patch%coarse_i_lower - 1)
    expected_i_upper = min(coarse_geometry%nx, patch%coarse_i_upper + 1)
    expected_j_lower = max(1, patch%coarse_j_lower - 1)
    expected_j_upper = min(coarse_geometry%ny, patch%coarse_j_upper + 1)
    if (component_count < 1 .or. &
        .not. patch%is_valid(coarse_geometry, fine_geometry) .or. &
        size(coarse_start, 1) /= component_count .or. &
        size(coarse_start, 2) < 1 .or. size(coarse_start, 3) < 1 .or. &
        coarse_i_lower < 1 .or. coarse_i_upper > coarse_geometry%nx .or. &
        coarse_j_lower < 1 .or. coarse_j_upper > coarse_geometry%ny .or. &
        coarse_i_lower > expected_i_lower .or. &
        coarse_i_upper < expected_i_upper .or. &
        coarse_j_lower > expected_j_lower .or. &
        coarse_j_upper < expected_j_upper .or. &
        any(shape(coarse_end) /= shape(coarse_start)) .or. &
        any(shape(coarse_start_temperature) /= &
          [size(coarse_start, 2), size(coarse_start, 3)]) .or. &
        any(shape(coarse_end_temperature) /= &
          shape(coarse_start_temperature)) .or. &
        any(.not. ieee_is_finite(coarse_start)) .or. &
        any(.not. ieee_is_finite(coarse_end)) .or. &
        any(.not. ieee_is_finite(coarse_start_temperature)) .or. &
        any(.not. ieee_is_finite(coarse_end_temperature))) return

    allocate(context%start%x_lower_state(component_count, fine_geometry%ny))
    allocate(context%start%x_upper_state(component_count, fine_geometry%ny))
    allocate(context%start%y_lower_state(component_count, fine_geometry%nx))
    allocate(context%start%y_upper_state(component_count, fine_geometry%nx))
    allocate(context%start%x_lower_temperature(fine_geometry%ny))
    allocate(context%start%x_upper_temperature(fine_geometry%ny))
    allocate(context%start%y_lower_temperature(fine_geometry%nx))
    allocate(context%start%y_upper_temperature(fine_geometry%nx))
    context%start%x_lower_state = 0.0_dp
    context%start%x_upper_state = 0.0_dp
    context%start%y_lower_state = 0.0_dp
    context%start%y_upper_state = 0.0_dp
    context%start%x_lower_temperature = 1.0_dp
    context%start%x_upper_temperature = 1.0_dp
    context%start%y_lower_temperature = 1.0_dp
    context%start%y_upper_temperature = 1.0_dp
    context%end = context%start
    ratio = patch%refinement_ratio

    do fine_j = 1, fine_geometry%ny
      coarse_j = patch%coarse_j_lower + (fine_j - 1) / ratio
      if (patch%coarse_i_lower > 1 .and. &
          fine_geometry%x_face_fraction(0, fine_j) > 0.0_dp) then
        coarse_i = patch%coarse_i_lower - 1
        context%start%x_lower_state(:, fine_j) = &
          coarse_start(:, coarse_i, coarse_j)
        context%end%x_lower_state(:, fine_j) = &
          coarse_end(:, coarse_i, coarse_j)
        context%start%x_lower_temperature(fine_j) = &
          coarse_start_temperature(coarse_i, coarse_j)
        context%end%x_lower_temperature(fine_j) = &
          coarse_end_temperature(coarse_i, coarse_j)
      end if
      if (patch%coarse_i_upper < coarse_geometry%nx .and. &
          fine_geometry%x_face_fraction(fine_geometry%nx, fine_j) > &
            0.0_dp) then
        coarse_i = patch%coarse_i_upper + 1
        context%start%x_upper_state(:, fine_j) = &
          coarse_start(:, coarse_i, coarse_j)
        context%end%x_upper_state(:, fine_j) = &
          coarse_end(:, coarse_i, coarse_j)
        context%start%x_upper_temperature(fine_j) = &
          coarse_start_temperature(coarse_i, coarse_j)
        context%end%x_upper_temperature(fine_j) = &
          coarse_end_temperature(coarse_i, coarse_j)
      end if
    end do
    do fine_i = 1, fine_geometry%nx
      coarse_i = patch%coarse_i_lower + (fine_i - 1) / ratio
      if (patch%coarse_j_lower > 1 .and. &
          fine_geometry%y_face_fraction(fine_i, 0) > 0.0_dp) then
        coarse_j = patch%coarse_j_lower - 1
        context%start%y_lower_state(:, fine_i) = &
          coarse_start(:, coarse_i, coarse_j)
        context%end%y_lower_state(:, fine_i) = &
          coarse_end(:, coarse_i, coarse_j)
        context%start%y_lower_temperature(fine_i) = &
          coarse_start_temperature(coarse_i, coarse_j)
        context%end%y_lower_temperature(fine_i) = &
          coarse_end_temperature(coarse_i, coarse_j)
      end if
      if (patch%coarse_j_upper < coarse_geometry%ny .and. &
          fine_geometry%y_face_fraction(fine_i, fine_geometry%ny) > &
            0.0_dp) then
        coarse_j = patch%coarse_j_upper + 1
        context%start%y_upper_state(:, fine_i) = &
          coarse_start(:, coarse_i, coarse_j)
        context%end%y_upper_state(:, fine_i) = &
          coarse_end(:, coarse_i, coarse_j)
        context%start%y_upper_temperature(fine_i) = &
          coarse_start_temperature(coarse_i, coarse_j)
        context%end%y_upper_temperature(fine_i) = &
          coarse_end_temperature(coarse_i, coarse_j)
      end if
    end do
    ok = context%is_valid(fine_geometry, component_count)
  end subroutine extract_reactive_eb_patch_exterior_context_support_2d

  subroutine extract_reactive_eb_patch_exterior_context_2d( &
      coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, coarse_geometry, fine_geometry, patch, &
      component_count, context, ok)
    real(dp), intent(in) :: coarse_start(:, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :)
    real(dp), intent(in) :: coarse_end(:, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    integer, intent(in) :: component_count
    type(reactive_eb_patch_exterior_context_2d), intent(out) :: context
    logical, intent(out) :: ok

    if (any(shape(coarse_start) /= &
        [component_count, coarse_geometry%nx, coarse_geometry%ny]) .or. &
        any(shape(coarse_start_temperature) /= &
          [coarse_geometry%nx, coarse_geometry%ny]) .or. &
        any(shape(coarse_end) /= shape(coarse_start)) .or. &
        any(shape(coarse_end_temperature) /= &
          shape(coarse_start_temperature))) then
      ok = .false.
      return
    end if
    call extract_reactive_eb_patch_exterior_context_support_2d( &
      1, 1, coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, coarse_geometry, fine_geometry, patch, &
      component_count, context, ok)
  end subroutine extract_reactive_eb_patch_exterior_context_2d

  subroutine build_reactive_eb_patch_exterior_from_context_2d( &
      species, context, coarse_geometry, fine_geometry, patch, alpha, &
      exterior, ok, fine_state, fine_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_patch_exterior_context_2d), intent(in) :: context
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    real(dp), intent(in) :: alpha
    type(reactive_eb_exterior_state_2d), intent(out) :: exterior
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: fine_state(:, :, :)
    real(dp), intent(in), optional :: fine_temperature(:, :)

    type(reactive_eb_exterior_state_2d) :: candidate
    logical :: local_ok, needs_physical_state
    integer :: fine_i, fine_j, nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    needs_physical_state = patch%coarse_i_lower == 1 .or. &
      patch%coarse_i_upper == coarse_geometry%nx .or. &
      patch%coarse_j_lower == 1 .or. &
      patch%coarse_j_upper == coarse_geometry%ny
    if (nvar < 1 .or. .not. ieee_is_finite(alpha) .or. &
        alpha < 0.0_dp .or. alpha > 1.0_dp .or. &
        .not. patch%is_valid(coarse_geometry, fine_geometry) .or. &
        .not. context%is_valid(fine_geometry, nvar)) return
    if (needs_physical_state) then
      if (.not. present(fine_state) .or. &
          .not. present(fine_temperature)) return
      if (any(shape(fine_state) /= &
          [nvar, fine_geometry%nx, fine_geometry%ny]) .or. &
          any(shape(fine_temperature) /= &
            [fine_geometry%nx, fine_geometry%ny]) .or. &
          any(.not. ieee_is_finite(fine_state)) .or. &
          any(.not. ieee_is_finite(fine_temperature)) .or. &
          any(fine_temperature <= 0.0_dp)) return
    end if

    allocate(candidate%x_lower_state(nvar, fine_geometry%ny))
    allocate(candidate%x_upper_state(nvar, fine_geometry%ny))
    allocate(candidate%y_lower_state(nvar, fine_geometry%nx))
    allocate(candidate%y_upper_state(nvar, fine_geometry%nx))
    allocate(candidate%x_lower_temperature(fine_geometry%ny))
    allocate(candidate%x_upper_temperature(fine_geometry%ny))
    allocate(candidate%y_lower_temperature(fine_geometry%nx))
    allocate(candidate%y_upper_temperature(fine_geometry%nx))
    candidate%x_lower_state = 0.0_dp
    candidate%x_upper_state = 0.0_dp
    candidate%y_lower_state = 0.0_dp
    candidate%y_upper_state = 0.0_dp
    candidate%x_lower_temperature = 1.0_dp
    candidate%x_upper_temperature = 1.0_dp
    candidate%y_lower_temperature = 1.0_dp
    candidate%y_upper_temperature = 1.0_dp

    do fine_j = 1, fine_geometry%ny
      if (fine_geometry%x_face_fraction(0, fine_j) > 0.0_dp) then
        if (patch%coarse_i_lower == 1) then
          candidate%x_lower_state(:, fine_j) = fine_state(:, 1, fine_j)
          candidate%x_lower_temperature(fine_j) = fine_temperature(1, fine_j)
        else
          call recover_exterior_cell( &
            species, context%start%x_lower_state(:, fine_j), &
            context%end%x_lower_state(:, fine_j), &
            context%start%x_lower_temperature(fine_j), &
            context%end%x_lower_temperature(fine_j), alpha, &
            candidate%x_lower_state(:, fine_j), &
            candidate%x_lower_temperature(fine_j), local_ok)
          if (.not. local_ok) return
        end if
      end if
      if (fine_geometry%x_face_fraction(fine_geometry%nx, fine_j) > 0.0_dp) then
        if (patch%coarse_i_upper == coarse_geometry%nx) then
          candidate%x_upper_state(:, fine_j) = &
            fine_state(:, fine_geometry%nx, fine_j)
          candidate%x_upper_temperature(fine_j) = &
            fine_temperature(fine_geometry%nx, fine_j)
        else
          call recover_exterior_cell( &
            species, context%start%x_upper_state(:, fine_j), &
            context%end%x_upper_state(:, fine_j), &
            context%start%x_upper_temperature(fine_j), &
            context%end%x_upper_temperature(fine_j), alpha, &
            candidate%x_upper_state(:, fine_j), &
            candidate%x_upper_temperature(fine_j), local_ok)
          if (.not. local_ok) return
        end if
      end if
    end do
    do fine_i = 1, fine_geometry%nx
      if (fine_geometry%y_face_fraction(fine_i, 0) > 0.0_dp) then
        if (patch%coarse_j_lower == 1) then
          candidate%y_lower_state(:, fine_i) = fine_state(:, fine_i, 1)
          candidate%y_lower_temperature(fine_i) = fine_temperature(fine_i, 1)
        else
          call recover_exterior_cell( &
            species, context%start%y_lower_state(:, fine_i), &
            context%end%y_lower_state(:, fine_i), &
            context%start%y_lower_temperature(fine_i), &
            context%end%y_lower_temperature(fine_i), alpha, &
            candidate%y_lower_state(:, fine_i), &
            candidate%y_lower_temperature(fine_i), local_ok)
          if (.not. local_ok) return
        end if
      end if
      if (fine_geometry%y_face_fraction(fine_i, fine_geometry%ny) > 0.0_dp) then
        if (patch%coarse_j_upper == coarse_geometry%ny) then
          candidate%y_upper_state(:, fine_i) = &
            fine_state(:, fine_i, fine_geometry%ny)
          candidate%y_upper_temperature(fine_i) = &
            fine_temperature(fine_i, fine_geometry%ny)
        else
          call recover_exterior_cell( &
            species, context%start%y_upper_state(:, fine_i), &
            context%end%y_upper_state(:, fine_i), &
            context%start%y_upper_temperature(fine_i), &
            context%end%y_upper_temperature(fine_i), alpha, &
            candidate%y_upper_state(:, fine_i), &
            candidate%y_upper_temperature(fine_i), local_ok)
          if (.not. local_ok) return
        end if
      end if
    end do
    if (.not. candidate%is_valid(fine_geometry, nvar)) return
    exterior = candidate
    ok = .true.
  end subroutine build_reactive_eb_patch_exterior_from_context_2d

  subroutine build_reactive_eb_patch_exterior_2d( &
      species, coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, coarse_geometry, fine_geometry, patch, &
      alpha, exterior, ok, fine_state, fine_temperature)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_start(:, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :)
    real(dp), intent(in) :: coarse_end(:, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    real(dp), intent(in) :: alpha
    type(reactive_eb_exterior_state_2d), intent(out) :: exterior
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: fine_state(:, :, :)
    real(dp), intent(in), optional :: fine_temperature(:, :)

    type(reactive_eb_patch_exterior_context_2d) :: context
    logical :: local_ok
    integer :: nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    call extract_reactive_eb_patch_exterior_context_2d( &
      coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, coarse_geometry, fine_geometry, patch, &
      nvar, context, local_ok)
    if (.not. local_ok) return
    call build_reactive_eb_patch_exterior_from_context_2d( &
      species, context, coarse_geometry, fine_geometry, patch, alpha, &
      exterior, ok, fine_state, fine_temperature)
  end subroutine build_reactive_eb_patch_exterior_2d

  subroutine advance_reactive_eb_level_2d( &
      species, state, temperature, geometry, solver, reconstruction, limiter, &
      state_redist_target_volume_fraction, state_redist_max_order, dt, &
      new_state, new_temperature, x_flux, y_flux, ok, exterior)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    character(len=*), intent(in) :: solver, reconstruction, limiter
    real(dp), intent(in) :: state_redist_target_volume_fraction
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: new_state(:, :, :), new_temperature(:, :)
    real(dp), intent(out) :: x_flux(:, 0:, :), y_flux(:, :, 0:)
    logical, intent(out) :: ok
    type(reactive_eb_exterior_state_2d), intent(in), optional :: exterior

    real(dp), allocatable :: center_x(:, :, :), center_y(:, :, :)
    real(dp), allocatable :: conservative_rhs(:, :, :)
    logical :: local_ok
    integer :: nvar

    new_state = state
    new_temperature = temperature
    x_flux = 0.0_dp
    y_flux = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. size(x_flux, 1) /= nvar .or. &
        size(x_flux, 2) /= geometry%nx + 1 .or. &
        size(x_flux, 3) /= geometry%ny .or. &
        size(y_flux, 1) /= nvar .or. size(y_flux, 2) /= geometry%nx .or. &
        size(y_flux, 3) /= geometry%ny + 1) return
    allocate(center_x(nvar, 0:geometry%nx, geometry%ny))
    allocate(center_y(nvar, geometry%nx, 0:geometry%ny))
    if (present(exterior)) then
      call build_reactive_eb_face_center_fluxes_2d( &
        species, state, temperature, geometry, solver, reconstruction, limiter, &
        dt, center_x, center_y, local_ok, exterior)
    else
      call build_reactive_eb_face_center_fluxes_2d( &
        species, state, temperature, geometry, solver, reconstruction, limiter, &
        dt, center_x, center_y, local_ok)
    end if
    if (.not. local_ok) return
    call interpolate_reactive_eb_face_centroid_fluxes_2d( &
      geometry, center_x, center_y, x_flux, y_flux, local_ok)
    if (.not. local_ok) return
    allocate(conservative_rhs(nvar, geometry%nx, geometry%ny))
    call reactive_eb_flux_divergence_2d( &
      species, state, temperature, geometry, x_flux, y_flux, &
      conservative_rhs, local_ok)
    if (.not. local_ok) return
    call advance_reactive_eb_state_redistributed_2d( &
      species, state, temperature, geometry, conservative_rhs, dt, &
      new_state, new_temperature, local_ok, &
      state_redist_target_volume_fraction, state_redist_max_order)
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_reactive_eb_level_2d

  subroutine advance_two_level_reactive_eb_hydro_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      fine_state, fine_temperature, fine_geometry, patch, solver, &
      reconstruction, limiter, state_redist_max_order, dt, new_coarse_state, &
      new_coarse_temperature, new_fine_state, new_fine_temperature, ok, &
      state_redist_target_volume_fraction)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), intent(in) :: fine_state(:, :, :), fine_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: new_coarse_state(:, :, :)
    real(dp), intent(out) :: new_coarse_temperature(:, :)
    real(dp), intent(out) :: new_fine_state(:, :, :)
    real(dp), intent(out) :: new_fine_temperature(:, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: state_redist_target_volume_fraction

    type(amr_eb_flux_register_2d) :: flux_register
    type(reactive_eb_exterior_state_2d) :: exterior
    real(dp), allocatable :: coarse_candidate(:, :, :), coarse_work(:, :, :)
    real(dp), allocatable :: coarse_candidate_temperature(:, :)
    real(dp), allocatable :: coarse_work_temperature(:, :)
    real(dp), allocatable :: fine_candidate(:, :, :), fine_work(:, :, :)
    real(dp), allocatable :: fine_candidate_temperature(:, :)
    real(dp), allocatable :: fine_work_temperature(:, :)
    real(dp), allocatable :: coarse_x_flux(:, :, :), coarse_y_flux(:, :, :)
    real(dp), allocatable :: fine_x_flux(:, :, :), fine_y_flux(:, :, :)
    real(dp) :: alpha, fine_dt, selected_target
    logical :: local_ok
    integer :: nvar, ratio, substep

    new_coarse_state = coarse_state
    new_coarse_temperature = coarse_temperature
    new_fine_state = fine_state
    new_fine_temperature = fine_temperature
    ok = .false.
    nvar = reactive_nvar(size(species))
    selected_target = 0.5_dp
    if (present(state_redist_target_volume_fraction)) &
      selected_target = state_redist_target_volume_fraction
    if (nvar < 1 .or. .not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. ieee_is_finite(selected_target) .or. &
        selected_target <= 0.0_dp .or. selected_target > 1.0_dp .or. &
        .not. patch%is_valid(coarse_geometry, fine_geometry) .or. &
        any(shape(new_coarse_state) /= shape(coarse_state)) .or. &
        any(shape(new_coarse_temperature) /= shape(coarse_temperature)) .or. &
        any(shape(new_fine_state) /= shape(fine_state)) .or. &
        any(shape(new_fine_temperature) /= shape(fine_temperature))) return

    allocate(coarse_candidate, mold=coarse_state)
    allocate(coarse_candidate_temperature, mold=coarse_temperature)
    allocate(coarse_x_flux(nvar, 0:coarse_geometry%nx, coarse_geometry%ny))
    allocate(coarse_y_flux(nvar, coarse_geometry%nx, 0:coarse_geometry%ny))
    call advance_reactive_eb_level_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, solver, &
      reconstruction, limiter, selected_target, state_redist_max_order, dt, &
      coarse_candidate, coarse_candidate_temperature, coarse_x_flux, &
      coarse_y_flux, local_ok)
    if (.not. local_ok) return

    call initialize_amr_eb_flux_register_2d( &
      coarse_geometry, fine_geometry, patch, nvar, flux_register, local_ok)
    if (.not. local_ok) return
    call accumulate_coarse_eb_fluxes_2d( &
      flux_register, coarse_geometry, fine_geometry, patch, &
      coarse_x_flux, coarse_y_flux, dt, local_ok)
    if (.not. local_ok) return

    allocate(fine_candidate, source=fine_state)
    allocate(fine_candidate_temperature, source=fine_temperature)
    allocate(fine_work, mold=fine_state)
    allocate(fine_work_temperature, mold=fine_temperature)
    allocate(fine_x_flux(nvar, 0:fine_geometry%nx, fine_geometry%ny))
    allocate(fine_y_flux(nvar, fine_geometry%nx, 0:fine_geometry%ny))
    ratio = patch%refinement_ratio
    fine_dt = dt / real(ratio, dp)
    do substep = 1, ratio
      if (trim(reconstruction) == "characteristic_plm") then
        alpha = (real(substep, dp) - 0.5_dp) / real(ratio, dp)
      else
        alpha = real(substep - 1, dp) / real(ratio, dp)
      end if
      call build_reactive_eb_patch_exterior_2d( &
        species, coarse_state, coarse_temperature, coarse_candidate, &
        coarse_candidate_temperature, coarse_geometry, fine_geometry, &
        patch, alpha, exterior, local_ok, fine_candidate, &
        fine_candidate_temperature)
      if (.not. local_ok) return
      call advance_reactive_eb_level_2d( &
        species, fine_candidate, fine_candidate_temperature, fine_geometry, &
        solver, reconstruction, limiter, selected_target, &
        state_redist_max_order, fine_dt, fine_work, fine_work_temperature, &
        fine_x_flux, fine_y_flux, local_ok, exterior)
      if (.not. local_ok) return
      fine_candidate = fine_work
      fine_candidate_temperature = fine_work_temperature
      call accumulate_fine_eb_fluxes_2d( &
        flux_register, coarse_geometry, fine_geometry, patch, &
        fine_x_flux, fine_y_flux, fine_dt, local_ok)
      if (.not. local_ok) return
    end do

    allocate(coarse_work, mold=coarse_state)
    allocate(coarse_work_temperature, mold=coarse_temperature)
    call reflux_reactive_eb_state_patch_2d( &
      species, coarse_candidate, coarse_candidate_temperature, &
      coarse_geometry, fine_candidate, fine_candidate_temperature, &
      fine_geometry, patch, flux_register, coarse_work, &
      coarse_work_temperature, fine_work, fine_work_temperature, local_ok)
    if (.not. local_ok) return
    call average_down_reactive_eb_state_patch_2d( &
      species, coarse_work, coarse_work_temperature, coarse_geometry, &
      fine_work, fine_geometry, patch, coarse_candidate, &
      coarse_candidate_temperature, local_ok)
    if (.not. local_ok) return

    new_coarse_state = coarse_candidate
    new_coarse_temperature = coarse_candidate_temperature
    new_fine_state = fine_work
    new_fine_temperature = fine_work_temperature
    ok = .true.
  end subroutine advance_two_level_reactive_eb_hydro_2d

end module amr_eb_reactive_2d_mod
