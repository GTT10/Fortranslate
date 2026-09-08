program test_mpi_amr_selected_context_mismatch
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use simulation_config_reactive_1d_mod, only: &
    reactive_1d_config, reactive_1d_mole_fractions
  use mpi_amr_reactive_1d_application_mod, only: &
    run_mpi_amr_reactive_1d_application
  implicit none

  character(len=*), parameter :: bundle_sha256 = &
    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_1d_config) :: config
  real(dp), allocatable :: composition(:)
  logical :: ok
  integer :: ierr, rank

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"
  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Context mismatch thermodynamics load failed"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "Context mismatch mechanism load failed"
  call load_h2o2_elementary_transport(transport, ok)
  if (.not. ok) error stop "Context mismatch transport load failed"
  call configure_case(config)
  allocate(composition(size(species)))
  call reactive_1d_mole_fractions(config, size(species), composition, ok)
  if (.not. ok) error stop "Context mismatch composition failed"
  if (rank == 1) then
    composition(1) = composition(1) + 1.0e-6_dp
    composition(2) = composition(2) - 1.0e-6_dp
  end if

  call run_mpi_amr_reactive_1d_application( &
    MPI_COMM_WORLD, "selected context mismatch test", config, species, &
    reactions, transport, "mpi_amr_context_mismatch.csv", bundle_sha256, &
    composition, "explicit")
  error stop "Rank-dependent selected context was unexpectedly accepted"

contains

  subroutine configure_case(local_config)
    type(reactive_1d_config), intent(out) :: local_config

    local_config = reactive_1d_config()
    local_config%nx = 8
    local_config%x_lower = 0.0_dp
    local_config%x_upper = 0.008_dp
    local_config%final_time = 1.0e-9_dp
    local_config%problem = "uniform_reactor"
    local_config%boundary_condition = "periodic"
    local_config%chemistry_model = "elementary"
    local_config%chemistry_enabled = .true.
    local_config%transport_enabled = .false.
    local_config%initial_temperature = 1200.0_dp
    local_config%initial_pressure = 101325.0_dp
    local_config%x_h2 = 0.29570_dp
    local_config%x_h = 1.0e-5_dp
    local_config%x_o = 1.0e-5_dp
    local_config%x_o2 = 0.14784_dp
    local_config%x_oh = 1.0e-5_dp
    local_config%x_h2o = 0.0_dp
    local_config%x_n2 = 0.55643_dp
    local_config%amr_enabled = .true.
    local_config%amr_max_levels = 2
    local_config%amr_minimum_patch_cells = 2
  end subroutine configure_case

end program test_mpi_amr_selected_context_mismatch
