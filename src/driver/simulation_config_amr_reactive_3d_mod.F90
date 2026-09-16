module simulation_config_amr_reactive_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use simulation_config_reactive_3d_mod, only: reactive_3d_config, resolved_reactive_boundaries_3d
  implicit none
  private

  type, public :: amr_reactive_3d_config
    integer :: regrid_interval_steps=0, regrid_buffer_cells=1
    real(dp) :: regrid_temperature_threshold=1100.0_dp
    integer :: refinement_ratio = 2
    integer :: coarse_i_lower = 3
    integer :: coarse_i_upper = 6
    integer :: coarse_j_lower = 3
    integer :: coarse_j_upper = 6
    integer :: coarse_k_lower = 3
    integer :: coarse_k_upper = 6
    character(len=512) :: coarse_output_file = &
      "amr_reactive_3d_coarse.csv"
    character(len=512) :: fine_output_file = &
      "amr_reactive_3d_fine.csv"
    integer :: checkpoint_interval_steps = 0
    logical :: stop_after_checkpoint = .false.
    character(len=512) :: checkpoint_file = &
      "amr_reactive_3d.chk"
    character(len=512) :: restart_file = ""
  end type amr_reactive_3d_config

  public :: read_amr_reactive_3d_configuration
  public :: validate_amr_reactive_3d_configuration

contains

  subroutine read_amr_reactive_3d_configuration( &
      path, base_config, config, ok, message, allow_selected,allow_physical,allow_boundary_patch)
    character(len=*), intent(in) :: path
    type(reactive_3d_config), intent(in) :: base_config
    type(amr_reactive_3d_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    logical, intent(in), optional :: allow_selected,allow_physical,allow_boundary_patch

    integer :: refinement_ratio,regrid_interval_steps,regrid_buffer_cells
    real(dp) :: regrid_temperature_threshold
    integer :: coarse_i_lower, coarse_i_upper
    integer :: coarse_j_lower, coarse_j_upper
    integer :: coarse_k_lower, coarse_k_upper
    integer :: checkpoint_interval_steps
    integer :: unit, io_status
    character(len=512) :: coarse_output_file, fine_output_file
    character(len=512) :: checkpoint_file, restart_file
    logical :: stop_after_checkpoint, selected_allowed

    namelist /amr_reactive_3d/ &
      refinement_ratio, regrid_interval_steps, regrid_buffer_cells, regrid_temperature_threshold, &
      coarse_i_lower, coarse_i_upper, &
      coarse_j_lower, coarse_j_upper, coarse_k_lower, coarse_k_upper, &
      coarse_output_file, fine_output_file, checkpoint_interval_steps, &
      stop_after_checkpoint, checkpoint_file, restart_file

    config = amr_reactive_3d_config()
    regrid_interval_steps=config%regrid_interval_steps
    regrid_buffer_cells=config%regrid_buffer_cells
    regrid_temperature_threshold=config%regrid_temperature_threshold
    refinement_ratio = config%refinement_ratio
    coarse_i_lower = config%coarse_i_lower
    coarse_i_upper = config%coarse_i_upper
    coarse_j_lower = config%coarse_j_lower
    coarse_j_upper = config%coarse_j_upper
    coarse_k_lower = config%coarse_k_lower
    coarse_k_upper = config%coarse_k_upper
    coarse_output_file = config%coarse_output_file
    fine_output_file = config%fine_output_file
    checkpoint_interval_steps = config%checkpoint_interval_steps
    stop_after_checkpoint = config%stop_after_checkpoint
    checkpoint_file = config%checkpoint_file
    restart_file = config%restart_file
    selected_allowed = .false.
    if (present(allow_selected)) selected_allowed = allow_selected
    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not open input file:", trim(path)
      return
    end if
    read(unit, nml=amr_reactive_3d, iostat=io_status)
    close(unit)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') &
        "Could not read &amr_reactive_3d from:", trim(path)
      return
    end if
    config%regrid_interval_steps=regrid_interval_steps
    config%regrid_buffer_cells=regrid_buffer_cells
    config%regrid_temperature_threshold=regrid_temperature_threshold
    config%refinement_ratio = refinement_ratio
    config%coarse_i_lower = coarse_i_lower
    config%coarse_i_upper = coarse_i_upper
    config%coarse_j_lower = coarse_j_lower
    config%coarse_j_upper = coarse_j_upper
    config%coarse_k_lower = coarse_k_lower
    config%coarse_k_upper = coarse_k_upper
    config%coarse_output_file = trim(coarse_output_file)
    config%fine_output_file = trim(fine_output_file)
    config%checkpoint_interval_steps = checkpoint_interval_steps
    config%stop_after_checkpoint = stop_after_checkpoint
    config%checkpoint_file = trim(checkpoint_file)
    config%restart_file = trim(restart_file)
    call validate_amr_reactive_3d_configuration( &
      base_config, config, ok, message, allow_selected=selected_allowed,allow_physical=allow_physical, &
      allow_boundary_patch=allow_boundary_patch)
  end subroutine read_amr_reactive_3d_configuration

  pure subroutine validate_amr_reactive_3d_configuration( &
      base_config, config, ok, message, allow_selected,allow_physical,allow_boundary_patch)
    type(reactive_3d_config), intent(in) :: base_config
    type(amr_reactive_3d_config), intent(in) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    logical, intent(in), optional :: allow_selected,allow_physical,allow_boundary_patch

    logical :: selected_allowed, selected_model,physical_allowed,touched(6)
    character(len=32) :: faces(6)
    integer :: margin

    ok = .false.
    message = ""
    selected_allowed = .false.
    if (present(allow_selected)) selected_allowed = allow_selected
    physical_allowed=.false.
    if (present(allow_physical)) physical_allowed=allow_physical
    margin=1
    if (present(allow_boundary_patch)) then
      if (allow_boundary_patch) margin=0
    end if
    faces=resolved_reactive_boundaries_3d(base_config)
    touched=[config%coarse_i_lower==1,config%coarse_i_upper==base_config%nx, &
      config%coarse_j_lower==1,config%coarse_j_upper==base_config%ny, &
      config%coarse_k_lower==1,config%coarse_k_upper==base_config%nz]
    selected_model = trim(base_config%thermo_model) == "selected"
    if (.not.physical_allowed.and.any(resolved_reactive_boundaries_3d(base_config) /= "periodic")) then
      message = "Static 3D AMR currently requires periodic boundaries"
    else if (trim(base_config%reconstruction) /= "pcm" .and. &
             trim(base_config%reconstruction) /= "characteristic_plm") then
      message = "Static 3D AMR reconstruction must be pcm or characteristic_plm"
    else if (trim(base_config%limiter) /= "minmod" .and. &
             trim(base_config%limiter) /= "mc") then
      message = "Static 3D AMR limiter must be minmod or mc"
    else if (trim(base_config%thermo_model) /= "elementary" .and. &
             trim(base_config%thermo_model) /= "full_h2o2" .and. &
             .not. (selected_model .and. selected_allowed)) then
      message = "Static 3D AMR requires a fixed H2/O2 thermodynamic model"
    else if (config%regrid_interval_steps<0.or.config%regrid_buffer_cells<0.or. &
        config%regrid_buffer_cells>max(base_config%nx,base_config%ny,base_config%nz)) then
      message='Invalid 3D AMR regrid interval or buffer'
    else if (.not.ieee_is_finite(config%regrid_temperature_threshold)) then
      message='Regrid temperature threshold must be finite'
    else if (config%regrid_temperature_threshold<=0.0_dp) then
      message='Regrid temperature threshold must be positive'
    else if (config%refinement_ratio < 2) then
      message = "3D AMR refinement_ratio must be at least two"
    else if (config%coarse_i_lower <= margin .or. &
             config%coarse_i_upper >= base_config%nx+1-margin .or. &
             config%coarse_i_upper < config%coarse_i_lower) then
      message = "3D AMR x patch lies outside the supported coarse-cell bounds"
    else if (config%coarse_j_lower <= margin .or. &
             config%coarse_j_upper >= base_config%ny+1-margin .or. &
             config%coarse_j_upper < config%coarse_j_lower) then
      message = "3D AMR y patch lies outside the supported coarse-cell bounds"
    else if (config%coarse_k_lower <= margin .or. &
             config%coarse_k_upper >= base_config%nz+1-margin .or. &
             config%coarse_k_upper < config%coarse_k_lower) then
      message = "3D AMR z patch lies outside the supported coarse-cell bounds"
    else if (any(touched.and.faces=='periodic')) then
      message='3D AMR patches touching periodic seams are not supported'
    else if (any(touched).and.config%regrid_interval_steps>0) then
      message='Boundary-touching 3D AMR currently requires a static patch'
    else if (len_trim(config%coarse_output_file) == 0 .or. &
             len_trim(config%fine_output_file) == 0) then
      message = "3D AMR output paths must be nonempty"
    else if (trim(config%coarse_output_file) == &
             trim(config%fine_output_file)) then
      message = "3D AMR coarse and fine output paths must differ"
    else if (config%checkpoint_interval_steps < 0) then
      message = "3D AMR checkpoint_interval_steps must be nonnegative"
    else if (config%stop_after_checkpoint .and. &
             config%checkpoint_interval_steps == 0) then
      message = "3D AMR checkpoint stop requires a positive interval"
    else if (config%checkpoint_interval_steps > 0 .and. &
             len_trim(config%checkpoint_file) == 0) then
      message = "3D AMR checkpoint path must be nonempty"
    else if (config%checkpoint_interval_steps > 0 .and. &
             (trim(config%checkpoint_file) == &
                trim(config%coarse_output_file) .or. &
              trim(config%checkpoint_file) == &
                trim(config%fine_output_file))) then
      message = "3D AMR checkpoint and CSV paths must differ"
    else if (len_trim(config%restart_file) > 0 .and. &
             (trim(config%restart_file) == &
                trim(config%coarse_output_file) .or. &
              trim(config%restart_file) == &
                trim(config%fine_output_file))) then
      message = "3D AMR restart input and CSV paths must differ"
    else
      ok = .true.
    end if
  end subroutine validate_amr_reactive_3d_configuration

end module simulation_config_amr_reactive_3d_mod
