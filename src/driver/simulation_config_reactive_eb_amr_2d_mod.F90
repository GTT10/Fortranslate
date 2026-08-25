module simulation_config_reactive_eb_amr_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use simulation_config_reactive_eb_2d_mod, only: &
    reactive_eb_2d_config, read_reactive_eb_2d_configuration
  implicit none
  private

  type, public :: reactive_eb_amr_2d_config
    type(reactive_eb_2d_config) :: eb
    integer :: coarse_i_lower = 2
    integer :: coarse_i_upper = 3
    integer :: coarse_j_lower = 2
    integer :: coarse_j_upper = 3
    integer :: refinement_ratio = 2
    logical :: dynamic_regridding = .false.
    logical :: regrid_at_initialization = .true.
    logical :: remove_fine_patch_when_untagged = .false.
    integer :: regrid_interval = 1
    real(dp) :: regrid_relative_temperature_gradient = 0.10_dp
    real(dp) :: regrid_absolute_temperature_gradient = 0.0_dp
    real(dp) :: regrid_temperature_scale_floor = 1.0_dp
    integer :: regrid_buffer_cells = 1
    integer :: regrid_minimum_patch_cells_x = 2
    integer :: regrid_minimum_patch_cells_y = 2
    character(len=1024) :: fine_output_file = "reactive_eb_amr_fine_2d.csv"
  end type reactive_eb_amr_2d_config

  public :: read_reactive_eb_amr_2d_configuration

contains

  subroutine read_reactive_eb_amr_2d_configuration( &
      path, config, ok, message)
    character(len=*), intent(in) :: path
    type(reactive_eb_amr_2d_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: coarse_i_lower, coarse_i_upper
    integer :: coarse_j_lower, coarse_j_upper, refinement_ratio
    integer :: regrid_interval, regrid_buffer_cells
    integer :: regrid_minimum_patch_cells_x
    integer :: regrid_minimum_patch_cells_y
    integer :: unit, status
    real(dp) :: regrid_relative_temperature_gradient
    real(dp) :: regrid_absolute_temperature_gradient
    real(dp) :: regrid_temperature_scale_floor
    logical :: dynamic_regridding, regrid_at_initialization
    logical :: remove_fine_patch_when_untagged
    character(len=1024) :: fine_output_file
    namelist /eb_amr/ coarse_i_lower, coarse_i_upper, &
      coarse_j_lower, coarse_j_upper, refinement_ratio, &
      dynamic_regridding, regrid_at_initialization, &
      remove_fine_patch_when_untagged, regrid_interval, &
      regrid_relative_temperature_gradient, &
      regrid_absolute_temperature_gradient, &
      regrid_temperature_scale_floor, regrid_buffer_cells, &
      regrid_minimum_patch_cells_x, regrid_minimum_patch_cells_y, &
      fine_output_file

    config = reactive_eb_amr_2d_config()
    call read_reactive_eb_2d_configuration( &
      path, config%eb, ok, message)
    if (.not. ok) return

    coarse_i_lower = config%coarse_i_lower
    coarse_i_upper = config%coarse_i_upper
    coarse_j_lower = config%coarse_j_lower
    coarse_j_upper = config%coarse_j_upper
    refinement_ratio = config%refinement_ratio
    dynamic_regridding = config%dynamic_regridding
    regrid_at_initialization = config%regrid_at_initialization
    remove_fine_patch_when_untagged = &
      config%remove_fine_patch_when_untagged
    regrid_interval = config%regrid_interval
    regrid_relative_temperature_gradient = &
      config%regrid_relative_temperature_gradient
    regrid_absolute_temperature_gradient = &
      config%regrid_absolute_temperature_gradient
    regrid_temperature_scale_floor = &
      config%regrid_temperature_scale_floor
    regrid_buffer_cells = config%regrid_buffer_cells
    regrid_minimum_patch_cells_x = config%regrid_minimum_patch_cells_x
    regrid_minimum_patch_cells_y = config%regrid_minimum_patch_cells_y
    fine_output_file = config%fine_output_file
    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=status)
    if (status /= 0) then
      ok = .false.
      message = "Could not open reactive EB AMR 2D input"
      return
    end if
    read(unit, nml=eb_amr, iostat=status)
    close(unit)
    if (status /= 0) then
      ok = .false.
      message = "Could not parse eb_amr namelist"
      return
    end if

    if (coarse_i_lower <= 1 .or. &
        coarse_i_upper >= config%eb%flow%nx .or. &
        coarse_j_lower <= 1 .or. &
        coarse_j_upper >= config%eb%flow%ny .or. &
        coarse_i_upper < coarse_i_lower .or. &
        coarse_j_upper < coarse_j_lower) then
      ok = .false.
      message = "EB AMR patch must be a strictly internal coarse rectangle"
      return
    end if
    if (refinement_ratio < 2) then
      ok = .false.
      message = "EB AMR refinement ratio must be at least two"
      return
    end if
    if (regrid_interval < 1 .or. regrid_buffer_cells < 0 .or. &
        regrid_minimum_patch_cells_x < 1 .or. &
        regrid_minimum_patch_cells_x > config%eb%flow%nx - 2 .or. &
        regrid_minimum_patch_cells_y < 1 .or. &
        regrid_minimum_patch_cells_y > config%eb%flow%ny - 2 .or. &
        .not. ieee_is_finite(regrid_relative_temperature_gradient) .or. &
        regrid_relative_temperature_gradient < 0.0_dp .or. &
        .not. ieee_is_finite(regrid_absolute_temperature_gradient) .or. &
        regrid_absolute_temperature_gradient < 0.0_dp .or. &
        .not. ieee_is_finite(regrid_temperature_scale_floor) .or. &
        regrid_temperature_scale_floor <= 0.0_dp) then
      ok = .false.
      message = "Invalid EB AMR dynamic-regridding controls"
      return
    end if
    if (remove_fine_patch_when_untagged .and. &
        .not. dynamic_regridding) then
      ok = .false.
      message = "Fine-patch removal requires dynamic EB AMR regridding"
      return
    end if
    if (len_trim(fine_output_file) == 0 .or. &
        trim(fine_output_file) == trim(config%eb%flow%output_file)) then
      ok = .false.
      message = "EB AMR fine output must be nonempty and distinct"
      return
    end if
    if (config%eb%flow%transport_enabled) then
      ok = .false.
      message = "Reactive EB AMR 2D does not yet support molecular transport"
      return
    end if

    config%coarse_i_lower = coarse_i_lower
    config%coarse_i_upper = coarse_i_upper
    config%coarse_j_lower = coarse_j_lower
    config%coarse_j_upper = coarse_j_upper
    config%refinement_ratio = refinement_ratio
    config%dynamic_regridding = dynamic_regridding
    config%regrid_at_initialization = regrid_at_initialization
    config%remove_fine_patch_when_untagged = &
      remove_fine_patch_when_untagged
    config%regrid_interval = regrid_interval
    config%regrid_relative_temperature_gradient = &
      regrid_relative_temperature_gradient
    config%regrid_absolute_temperature_gradient = &
      regrid_absolute_temperature_gradient
    config%regrid_temperature_scale_floor = regrid_temperature_scale_floor
    config%regrid_buffer_cells = regrid_buffer_cells
    config%regrid_minimum_patch_cells_x = regrid_minimum_patch_cells_x
    config%regrid_minimum_patch_cells_y = regrid_minimum_patch_cells_y
    config%fine_output_file = trim(fine_output_file)
    message = ""
    ok = .true.
  end subroutine read_reactive_eb_amr_2d_configuration

end module simulation_config_reactive_eb_amr_2d_mod
