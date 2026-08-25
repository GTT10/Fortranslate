module simulation_config_reactive_eb_amr_2d_mod
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
    integer :: unit, status
    character(len=1024) :: fine_output_file
    namelist /eb_amr/ coarse_i_lower, coarse_i_upper, &
      coarse_j_lower, coarse_j_upper, refinement_ratio, fine_output_file

    config = reactive_eb_amr_2d_config()
    call read_reactive_eb_2d_configuration( &
      path, config%eb, ok, message)
    if (.not. ok) return

    coarse_i_lower = config%coarse_i_lower
    coarse_i_upper = config%coarse_i_upper
    coarse_j_lower = config%coarse_j_lower
    coarse_j_upper = config%coarse_j_upper
    refinement_ratio = config%refinement_ratio
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
    if (len_trim(fine_output_file) == 0 .or. &
        trim(fine_output_file) == trim(config%eb%flow%output_file)) then
      ok = .false.
      message = "EB AMR fine output must be nonempty and distinct"
      return
    end if
    if (config%eb%flow%chemistry_enabled .or. &
        config%eb%flow%transport_enabled) then
      ok = .false.
      message = "Reactive EB AMR 2D currently supports hydrodynamics only"
      return
    end if

    config%coarse_i_lower = coarse_i_lower
    config%coarse_i_upper = coarse_i_upper
    config%coarse_j_lower = coarse_j_lower
    config%coarse_j_upper = coarse_j_upper
    config%refinement_ratio = refinement_ratio
    config%fine_output_file = trim(fine_output_file)
    message = ""
    ok = .true.
  end subroutine read_reactive_eb_amr_2d_configuration

end module simulation_config_reactive_eb_amr_2d_mod
