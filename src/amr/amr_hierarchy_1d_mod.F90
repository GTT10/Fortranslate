module amr_hierarchy_1d_mod
  use precision_mod, only: dp
  use slope_limiter_mod, only: limited_slope
  implicit none
  private

  type, public :: amr_box_1d
    integer :: level = -1
    integer :: lower = 1
    integer :: upper = 0
  contains
    procedure :: cell_count => amr_box_cell_count
    procedure :: is_valid => amr_box_is_valid
  end type amr_box_1d

  type, public :: amr_two_level_hierarchy_1d
    type(amr_box_1d) :: coarse
    type(amr_box_1d) :: fine
    integer :: refinement_ratio = 0
    integer :: fine_coarse_lower = 1
    integer :: fine_coarse_upper = 0
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 0.0_dp
    real(dp) :: coarse_dx = 0.0_dp
    real(dp) :: fine_dx = 0.0_dp
  contains
    procedure :: covered_coarse_cells => amr_covered_coarse_cells
    procedure :: is_valid => amr_hierarchy_is_valid
  end type amr_two_level_hierarchy_1d

  type, public :: amr_flux_register_1d
    real(dp), allocatable :: left(:)
    real(dp), allocatable :: right(:)
  contains
    procedure :: reset => reset_flux_register_1d
  end type amr_flux_register_1d

  public :: initialize_two_level_hierarchy_1d
  public :: prolong_conservative_1d
  public :: restrict_average_1d
  public :: average_down_1d
  public :: level_subcycle_time_steps_1d
  public :: initialize_flux_register_1d
  public :: accumulate_coarse_flux_1d
  public :: accumulate_fine_flux_1d
  public :: reflux_1d
  public :: composite_integral_1d

contains

  pure integer function amr_box_cell_count(self) result(count)
    class(amr_box_1d), intent(in) :: self

    count = max(0, self%upper - self%lower + 1)
  end function amr_box_cell_count

  pure logical function amr_box_is_valid(self) result(valid)
    class(amr_box_1d), intent(in) :: self

    valid = self%level >= 0 .and. self%lower >= 1 .and. &
      self%upper >= self%lower
  end function amr_box_is_valid

  pure integer function amr_covered_coarse_cells(self) result(count)
    class(amr_two_level_hierarchy_1d), intent(in) :: self

    count = max(0, self%fine_coarse_upper - self%fine_coarse_lower + 1)
  end function amr_covered_coarse_cells

  pure logical function amr_hierarchy_is_valid(self) result(valid)
    class(amr_two_level_hierarchy_1d), intent(in) :: self
    real(dp) :: spacing_tolerance

    spacing_tolerance = 32.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(self%coarse_dx), abs(self%fine_dx))
    valid = self%coarse%is_valid() .and. self%fine%is_valid() .and. &
      self%coarse%level == 0 .and. self%fine%level == 1 .and. &
      self%refinement_ratio >= 2 .and. self%x_upper > self%x_lower .and. &
      self%coarse_dx > 0.0_dp .and. self%fine_dx > 0.0_dp .and. &
      self%fine_coarse_lower > self%coarse%lower .and. &
      self%fine_coarse_upper < self%coarse%upper .and. &
      self%fine_coarse_upper >= self%fine_coarse_lower .and. &
      self%fine%lower == &
        (self%fine_coarse_lower - 1) * self%refinement_ratio + 1 .and. &
      self%fine%upper == &
        self%fine_coarse_upper * self%refinement_ratio .and. &
      self%fine%cell_count() == &
        self%covered_coarse_cells() * self%refinement_ratio .and. &
      abs(self%fine_dx * real(self%refinement_ratio, dp) - &
        self%coarse_dx) <= spacing_tolerance
  end function amr_hierarchy_is_valid

  pure subroutine initialize_two_level_hierarchy_1d( &
      coarse_cells, fine_coarse_lower, fine_coarse_upper, &
      refinement_ratio, x_lower, x_upper, hierarchy, ok)
    integer, intent(in) :: coarse_cells
    integer, intent(in) :: fine_coarse_lower, fine_coarse_upper
    integer, intent(in) :: refinement_ratio
    real(dp), intent(in) :: x_lower, x_upper
    type(amr_two_level_hierarchy_1d), intent(out) :: hierarchy
    logical, intent(out) :: ok

    hierarchy%coarse%level = 0
    hierarchy%coarse%lower = 1
    hierarchy%coarse%upper = coarse_cells
    hierarchy%fine%level = 1
    hierarchy%fine%lower = &
      (fine_coarse_lower - 1) * refinement_ratio + 1
    hierarchy%fine%upper = fine_coarse_upper * refinement_ratio
    hierarchy%refinement_ratio = refinement_ratio
    hierarchy%fine_coarse_lower = fine_coarse_lower
    hierarchy%fine_coarse_upper = fine_coarse_upper
    hierarchy%x_lower = x_lower
    hierarchy%x_upper = x_upper
    hierarchy%coarse_dx = 0.0_dp
    hierarchy%fine_dx = 0.0_dp
    if (coarse_cells > 0 .and. refinement_ratio > 0 .and. &
        x_upper > x_lower) then
      hierarchy%coarse_dx = &
        (x_upper - x_lower) / real(coarse_cells, dp)
      hierarchy%fine_dx = &
        hierarchy%coarse_dx / real(refinement_ratio, dp)
    end if
    ok = hierarchy%is_valid()
  end subroutine initialize_two_level_hierarchy_1d

  subroutine prolong_conservative_1d(coarse, hierarchy, fine, ok)
    real(dp), intent(in) :: coarse(:, :)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(out) :: fine(:, :)
    logical, intent(out) :: ok

    real(dp) :: delta_minus, delta_plus, slope, offset
    logical :: slope_ok
    integer :: component, coarse_cell, child, fine_cell

    fine = 0.0_dp
    ok = hierarchy%is_valid() .and. size(coarse, 1) == size(fine, 1) .and. &
      size(coarse, 1) >= 1 .and. &
      size(coarse, 2) == hierarchy%coarse%cell_count() .and. &
      size(fine, 2) == hierarchy%fine%cell_count()
    if (.not. ok) return

    do component = 1, size(coarse, 1)
      do coarse_cell = hierarchy%fine_coarse_lower, &
          hierarchy%fine_coarse_upper
        delta_minus = coarse(component, coarse_cell) - &
          coarse(component, coarse_cell - 1)
        delta_plus = coarse(component, coarse_cell + 1) - &
          coarse(component, coarse_cell)
        call limited_slope( &
          delta_minus, delta_plus, "mc", slope, slope_ok)
        if (.not. slope_ok) then
          ok = .false.
          return
        end if
        do child = 1, hierarchy%refinement_ratio
          fine_cell = &
            (coarse_cell - hierarchy%fine_coarse_lower) * &
            hierarchy%refinement_ratio + child
          offset = (real(child, dp) - 0.5_dp) / &
            real(hierarchy%refinement_ratio, dp) - 0.5_dp
          fine(component, fine_cell) = &
            coarse(component, coarse_cell) + slope * offset
        end do
      end do
    end do
  end subroutine prolong_conservative_1d

  pure subroutine restrict_average_1d(fine, hierarchy, restricted, ok)
    real(dp), intent(in) :: fine(:, :)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(out) :: restricted(:, :)
    logical, intent(out) :: ok

    integer :: component, covered_cell, first_child, last_child

    restricted = 0.0_dp
    ok = hierarchy%is_valid() .and. &
      size(fine, 1) == size(restricted, 1) .and. &
      size(fine, 1) >= 1 .and. &
      size(fine, 2) == hierarchy%fine%cell_count() .and. &
      size(restricted, 2) == hierarchy%covered_coarse_cells()
    if (.not. ok) return

    do covered_cell = 1, hierarchy%covered_coarse_cells()
      first_child = &
        (covered_cell - 1) * hierarchy%refinement_ratio + 1
      last_child = first_child + hierarchy%refinement_ratio - 1
      do component = 1, size(fine, 1)
        restricted(component, covered_cell) = &
          sum(fine(component, first_child:last_child)) / &
          real(hierarchy%refinement_ratio, dp)
      end do
    end do
  end subroutine restrict_average_1d

  subroutine average_down_1d(fine, hierarchy, coarse, ok)
    real(dp), intent(in) :: fine(:, :)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(inout) :: coarse(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: restricted(:, :)

    ok = hierarchy%is_valid() .and. size(fine, 1) == size(coarse, 1) .and. &
      size(fine, 1) >= 1 .and. &
      size(fine, 2) == hierarchy%fine%cell_count() .and. &
      size(coarse, 2) == hierarchy%coarse%cell_count()
    if (.not. ok) return
    allocate(restricted(size(coarse, 1), hierarchy%covered_coarse_cells()))
    call restrict_average_1d(fine, hierarchy, restricted, ok)
    if (.not. ok) return
    coarse(:, hierarchy%fine_coarse_lower:hierarchy%fine_coarse_upper) = &
      restricted
  end subroutine average_down_1d

  pure subroutine level_subcycle_time_steps_1d( &
      hierarchy, coarse_time_step, fine_time_steps, ok)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(in) :: coarse_time_step
    real(dp), intent(out) :: fine_time_steps(:)
    logical, intent(out) :: ok

    fine_time_steps = 0.0_dp
    ok = hierarchy%is_valid() .and. coarse_time_step > 0.0_dp .and. &
      size(fine_time_steps) == hierarchy%refinement_ratio
    if (.not. ok) return
    fine_time_steps = coarse_time_step / real(hierarchy%refinement_ratio, dp)
  end subroutine level_subcycle_time_steps_1d

  subroutine initialize_flux_register_1d(register, variable_count, ok)
    type(amr_flux_register_1d), intent(out) :: register
    integer, intent(in) :: variable_count
    logical, intent(out) :: ok

    ok = variable_count >= 1
    if (.not. ok) return
    allocate(register%left(variable_count), register%right(variable_count))
    call register%reset()
  end subroutine initialize_flux_register_1d

  subroutine reset_flux_register_1d(self)
    class(amr_flux_register_1d), intent(inout) :: self

    if (allocated(self%left)) self%left = 0.0_dp
    if (allocated(self%right)) self%right = 0.0_dp
  end subroutine reset_flux_register_1d

  subroutine accumulate_coarse_flux_1d( &
      register, left_flux, right_flux, time_step, ok)
    type(amr_flux_register_1d), intent(inout) :: register
    real(dp), intent(in) :: left_flux(:), right_flux(:), time_step
    logical, intent(out) :: ok

    ok = valid_flux_register(register, size(left_flux)) .and. &
      size(right_flux) == size(left_flux) .and. time_step >= 0.0_dp
    if (.not. ok) return
    register%left = register%left - time_step * left_flux
    register%right = register%right - time_step * right_flux
  end subroutine accumulate_coarse_flux_1d

  subroutine accumulate_fine_flux_1d( &
      register, left_flux, right_flux, time_step, ok)
    type(amr_flux_register_1d), intent(inout) :: register
    real(dp), intent(in) :: left_flux(:), right_flux(:), time_step
    logical, intent(out) :: ok

    ok = valid_flux_register(register, size(left_flux)) .and. &
      size(right_flux) == size(left_flux) .and. time_step >= 0.0_dp
    if (.not. ok) return
    register%left = register%left + time_step * left_flux
    register%right = register%right + time_step * right_flux
  end subroutine accumulate_fine_flux_1d

  subroutine reflux_1d(coarse, hierarchy, register, ok)
    real(dp), intent(inout) :: coarse(:, :)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    type(amr_flux_register_1d), intent(inout) :: register
    logical, intent(out) :: ok

    ok = hierarchy%is_valid() .and. &
      size(coarse, 2) == hierarchy%coarse%cell_count() .and. &
      valid_flux_register(register, size(coarse, 1))
    if (.not. ok) return
    coarse(:, hierarchy%fine_coarse_lower - 1) = &
      coarse(:, hierarchy%fine_coarse_lower - 1) - &
      register%left / hierarchy%coarse_dx
    coarse(:, hierarchy%fine_coarse_upper + 1) = &
      coarse(:, hierarchy%fine_coarse_upper + 1) + &
      register%right / hierarchy%coarse_dx
    call register%reset()
  end subroutine reflux_1d

  pure subroutine composite_integral_1d( &
      coarse, fine, hierarchy, integral, ok)
    real(dp), intent(in) :: coarse(:, :), fine(:, :)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    integer :: component

    integral = 0.0_dp
    ok = hierarchy%is_valid() .and. size(coarse, 1) == size(fine, 1) .and. &
      size(integral) == size(coarse, 1) .and. &
      size(coarse, 2) == hierarchy%coarse%cell_count() .and. &
      size(fine, 2) == hierarchy%fine%cell_count()
    if (.not. ok) return

    do component = 1, size(integral)
      integral(component) = hierarchy%coarse_dx * ( &
        sum(coarse(component, 1:hierarchy%fine_coarse_lower - 1)) + &
        sum(coarse(component, hierarchy%fine_coarse_upper + 1: &
          hierarchy%coarse%upper))) + &
        hierarchy%fine_dx * sum(fine(component, :))
    end do
  end subroutine composite_integral_1d

  pure logical function valid_flux_register(register, variable_count) &
      result(valid)
    type(amr_flux_register_1d), intent(in) :: register
    integer, intent(in) :: variable_count

    valid = allocated(register%left) .and. allocated(register%right)
    if (.not. valid) return
    valid = variable_count >= 1 .and. &
      size(register%left) == variable_count .and. &
      size(register%right) == variable_count
  end function valid_flux_register

end module amr_hierarchy_1d_mod
