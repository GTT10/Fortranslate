program pelef_reactive_eb_3d
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_elementary_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use simulation_config_reactive_eb_3d_mod, only: &
    reactive_eb_3d_config, read_reactive_eb_3d_configuration
  use reactive_eb_3d_application_mod, only: &
    run_reactive_eb_3d_application
  use constant_volume_reactor_mod, only: forced_reactor_integrator
#ifdef PELEF_SDC_CVODE
  use sundials_forced_reactor_mod, only: advance_forced_reactor_cvode
#endif
  implicit none
#ifdef PELEF_SDC_CVODE
  logical, parameter :: cvode_available=.true.
#else
  logical, parameter :: cvode_available=.false.
#endif

  type(reactive_eb_3d_config) :: config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  character(len=1024) :: input_path, message
  logical :: ok

  if (command_argument_count() /= 1) then
    write(*, '(a)') "Usage: pelef_reactive_eb_3d <input.nml>"
    error stop 2
  end if
  call get_command_argument(1, input_path)
  call read_reactive_eb_3d_configuration( &
    trim(input_path), config, ok, message,allow_difmag=.true.,allow_ctu=.true., &
    allow_boundaries=.true.,allow_inlet_composition=.true.,allow_sdc=.true.,allow_sdc_cvode=cvode_available, &
    allow_pinned_eos=.true.)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if
  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load elementary thermodynamics"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load elementary mechanism"
  call load_h2o2_elementary_transport(transport, ok)
  if (.not. ok) error stop "Failed to load elementary transport"
#ifdef PELEF_SDC_CVODE
  if(config%sdc_reactor=='cvode') then
    call run(advance_forced_reactor_cvode)
  else
    call run()
  endif
#else
  call run()
#endif
contains
  subroutine run(evaluate_reaction)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
  call run_reactive_eb_3d_application( &
    trim(input_path), "reactive planar EB 3D", "", config, species, &
    reactions, transport,evaluate_reaction=evaluate_reaction)
  end subroutine run
end program pelef_reactive_eb_3d
