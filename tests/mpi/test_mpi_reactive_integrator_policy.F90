program test_mpi_reactive_integrator_policy
  use mpi_f08
  use precision_mod, only: dp
  use mpi_domain_1d_mod, only: mpi_domain_1d, initialize_mpi_domain_1d
  use mpi_reactive_1d_mod, only: advance_mpi_reactive_strang_adaptive
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_full_transport
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved
  implicit none

  type(mpi_domain_1d) :: domain
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp), allocatable :: state(:, :), saved_state(:, :)
  real(dp), allocatable :: temperature(:), saved_temperature(:)
  real(dp), allocatable :: primitive(:), mole_fractions(:), mass_fractions(:)
  real(dp), allocatable :: conserved(:)
  real(dp) :: density, cell_temperature, sound_speed, accepted_interval
  integer :: ierr, nvar, i, k, rejected_trials
  logical :: ok, local_pass, global_pass
  character(len=16) :: divergent_policy

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call initialize_mpi_domain_1d(domain, 19, MPI_COMM_WORLD, ok)
  if (.not. ok) error stop "MPI domain initialization failed"
  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop "Thermodynamics load failed"
  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop "Mechanism load failed"
  call load_h2o2_full_transport(transport, ok)
  if (.not. ok) error stop "Transport load failed"

  nvar = reactive_nvar(size(species))
  allocate(state(nvar, 0:domain%local_cells + 1))
  allocate(temperature(0:domain%local_cells + 1))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(mole_fractions(size(species)), mass_fractions(size(species)))
  allocate(conserved(nvar))
  mole_fractions = 0.0_dp
  mole_fractions(1) = 2.0_dp
  mole_fractions(4) = 1.0_dp
  mole_fractions(10) = 3.0_dp
  mole_fractions = mole_fractions / sum(mole_fractions)
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  if (.not. ok) error stop "Composition conversion failed"
  cell_temperature = 1000.0_dp
  density = mixture_density( &
    species, mass_fractions, 101325.0_dp, cell_temperature, ok)
  if (.not. ok) error stop "Mixture density failed"
  primitive = 0.0_dp
  primitive(1:5) = [density, 0.0_dp, 0.0_dp, 0.0_dp, 101325.0_dp]
  do k = 1, size(species)
    primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, conserved, cell_temperature, sound_speed, ok)
  if (.not. ok) error stop "Reactive state initialization failed"
  do i = 0, domain%local_cells + 1
    state(:, i) = conserved
    temperature(i) = cell_temperature
  end do
  saved_state = state
  saved_temperature = temperature

  call advance_mpi_reactive_strang_adaptive( &
    domain, species, reactions, transport, state, temperature, 1.0e-4_dp, &
    1.0e-8_dp, 1.0e-12_dp, .true., 1.0e-7_dp, 1.0e-13_dp, &
    .false., .false., .false., .false., .false., "rusanov", &
    accepted_interval, rejected_trials, ok, "implicit")
  local_pass = ok .and. accepted_interval > 0.0_dp .and. &
    rejected_trials == 0 .and. &
    any(state(:, 1:domain%local_cells) /= &
      saved_state(:, 1:domain%local_cells))
  call MPI_Allreduce( &
    local_pass, global_pass, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS .or. .not. global_pass) then
    error stop "MPI chemistry-integrator control advance failed"
  end if

  state = saved_state
  temperature = saved_temperature
  if (domain%nranks == 1) then
    divergent_policy = "invalid"
  else if (domain%rank == 0) then
    divergent_policy = "explicit"
  else
    divergent_policy = "implicit"
  end if
  call advance_mpi_reactive_strang_adaptive( &
    domain, species, reactions, transport, state, temperature, 1.0e-4_dp, &
    1.0e-8_dp, 1.0e-12_dp, .true., 1.0e-7_dp, 1.0e-13_dp, &
    .false., .false., .false., .false., .false., "rusanov", &
    accepted_interval, rejected_trials, ok, divergent_policy)
  local_pass = .not. ok .and. accepted_interval == 0.0_dp .and. &
    rejected_trials > 0 .and. all(state == saved_state) .and. &
    all(temperature == saved_temperature)
  call MPI_Allreduce( &
    local_pass, global_pass, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS .or. .not. global_pass) then
    error stop "MPI chemistry-integrator consensus rollback failed"
  end if
  if (domain%rank == 0) then
    write(*, '(a)') "test_mpi_reactive_integrator_policy: PASS"
  end if

  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"
end program test_mpi_reactive_integrator_policy
