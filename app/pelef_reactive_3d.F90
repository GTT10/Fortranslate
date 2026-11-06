program pelef_reactive_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: &
    load_h2o2_elementary_transport, load_h2o2_full_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo,load_gamma_law_thermo
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use simulation_config_reactive_3d_mod, only: &
    reactive_3d_config, read_reactive_3d_configuration, &
    reactive_3d_mole_fractions
#ifdef PELEF_REGULAR_LES
  use les_regular_application_mod, only: run_reactive_3d_application=>run_regular_les_application_3d
#else
  use reactive_3d_application_mod, only: run_reactive_3d_application
#endif
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

  type(reactive_3d_config) :: config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp), allocatable :: mole_fractions(:)
  character(len=1024) :: input_path, message
  logical :: ok
  real(dp) :: gamma
  integer :: unit,ios
  namelist /gamma_law/ gamma

  if (command_argument_count() /= 1) then
    write(*, '(a)') "Usage: pelef_reactive_3d <input.nml>"
    error stop 2
  end if
  call get_command_argument(1, input_path)
  call read_reactive_3d_configuration(trim(input_path), config, ok, message, &
    allow_pressure_outlet=.true., allow_mach_inlet=.true., allow_characteristic_ppm=.true., &
    allow_difmag=.true.,allow_external_forcing=.true.,allow_sdc=.true.,allow_mol=.true.,allow_dual_energy=.true., &
    allow_sdc_cvode=cvode_available,allow_gamma_law=.true., &
    allow_pinned_eos=.true.)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if

  select case (trim(config%thermo_model))
  case('gamma_law')
    gamma=1.4_dp
    open(newunit=unit,file=trim(input_path),status='old',action='read',iostat=ios)
    if(ios/=0) error stop 'Cannot open gamma-law input'
    read(unit,nml=gamma_law,iostat=ios)
    close(unit)
    if(ios>0) error stop 'Invalid gamma_law namelist'
    call load_gamma_law_thermo(gamma,species,ok)
    if(.not.ok) error stop 'Invalid gamma-law heat-capacity ratio'
    allocate(reactions(0),transport(1))
    ! No molecular database: only prescribed, nondiffusing transport is allowed.
    transport(1)%name='AIR'
    write(*,'(a,es24.16)') 'Gamma-law ratio: ',gamma
  case ("elementary")
    call load_h2o2_elementary_thermo(species, ok)
    if (.not. ok) error stop "Failed to load elementary thermodynamics"
    call load_h2o2_elementary_mechanism(reactions, ok)
    if (.not. ok) error stop "Failed to load elementary mechanism"
    call load_h2o2_elementary_transport(transport, ok)
  case ("full_h2o2")
    call load_h2o2_full_thermo(species, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 thermodynamics"
    call load_h2o2_full_mechanism(reactions, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 mechanism"
    call load_h2o2_full_transport(transport, ok)
  case default
    ok = .false.
  end select
  if (.not. ok .or. .not. allocated(species) .or. &
      .not. allocated(reactions) .or. .not. allocated(transport)) then
    error stop "Failed to load reactive 3D model"
  end if

  allocate(mole_fractions(size(species)))
  call reactive_3d_mole_fractions( &
    config, size(species), mole_fractions, ok)
  if (.not. ok) error stop "Failed to resolve reactive 3D composition"
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
  call run_reactive_3d_application( &
    trim(input_path), "3D reactive flow", "", config, species, reactions, &
    transport, mole_fractions,evaluate_reaction=evaluate_reaction)
  end subroutine run
end program pelef_reactive_3d
