program test_mpi_amr_transport_context_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: load_h2o2_elementary_transport
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use mpi_amr_reactive_3d_application_mod, only: &
    validate_mpi_amr_transport_context_3d
  implicit none

  type(nasa7_species), allocatable :: baseline_species(:), species(:)
  type(gas_transport_species), allocatable :: baseline_transport(:), transport(:)
  type(reactive_3d_config) :: baseline_config, config
  logical :: ok, local_pass, global_pass
  integer :: ierr, rank, nranks

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_size failed"
  if (nranks /= 2) then
    if (rank == 0) write(*, '(a,i0)') &
      "test_mpi_amr_transport_context_3d requires 2 ranks, got ", nranks
    call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
    error stop "Invalid MPI rank count"
  end if

  call load_h2o2_elementary_thermo(baseline_species, ok)
  if (.not. ok) then
    call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
    error stop "Failed to load elementary thermodynamics"
  end if
  call load_h2o2_elementary_transport(baseline_transport, ok)
  if (.not. ok) then
    call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
    error stop "Failed to load elementary transport"
  end if

  baseline_config = reactive_3d_config()
  baseline_config%transport_enabled = .true.
  baseline_config%viscosity_enabled = .true.
  baseline_config%thermal_conduction_enabled = .true.
  baseline_config%species_diffusion_enabled = .true.
  baseline_config%barodiffusion_enabled = .false.
  baseline_config%transport_cfl = 0.35_dp
  baseline_config%reconstruction = "characteristic_plm"
  baseline_config%limiter = "mc"

  local_pass = .true.

  call reset_context()
  call expect_context(.true., "matching context")

  call reset_context()
  if (rank == 1) transport(1)%well_depth = &
    nearest(transport(1)%well_depth, 1.0_dp)
  call expect_context(.false., "nonroot transport value mismatch")

  call reset_context()
  if (rank == 1) transport(1)%name(1:1) = "X"
  call expect_context(.false., "nonroot transport name mismatch")

  call reset_context()
  if (rank == 1) transport(1)%geometry = mod(transport(1)%geometry + 1, 3)
  call expect_context(.false., "nonroot transport geometry mismatch")

  call reset_context()
  if (rank == 1) config%transport_enabled = .false.
  call expect_context(.false., "transport-enabled policy mismatch")

  call reset_context()
  if (rank == 0) config%transport_enabled = .false.
  call expect_context(.false., "root transport-enabled policy mismatch")

  call reset_context()
  if (rank == 1) config%viscosity_enabled = .false.
  call expect_context(.false., "viscosity policy mismatch")

  call reset_context()
  if (rank == 1) config%thermal_conduction_enabled = .false.
  call expect_context(.false., "thermal-conduction policy mismatch")

  call reset_context()
  if (rank == 1) config%species_diffusion_enabled = .false.
  call expect_context(.false., "species-diffusion policy mismatch")

  call reset_context()
  if (rank == 1) config%barodiffusion_enabled = .true.
  call expect_context(.false., "barodiffusion policy mismatch")

  call reset_context()
  if (rank == 1) config%transport_cfl = &
    nearest(config%transport_cfl, 1.0_dp)
  call expect_context(.false., "transport CFL mismatch")

  call reset_context()
  config%transport_cfl = ieee_value(0.0_dp, ieee_quiet_nan)
  call expect_context(.false., "non-finite transport CFL")

  call reset_context()
  config%transport_cfl = 0.0_dp
  call expect_context(.false., "non-positive transport CFL")

  call reset_context()
  if (rank == 1) config%reconstruction = "pcm"
  call expect_context(.false., "reconstruction mismatch")

  call reset_context()
  if (rank == 1) config%limiter = "minmod"
  call expect_context(.false., "limiter mismatch")

  call MPI_Allreduce( &
    local_pass, global_pass, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS) then
    if (rank == 0) write(*, '(a)') &
      "test_mpi_amr_transport_context_3d: FAIL"
    call MPI_Finalize(ierr)
    error stop "MPI transport context validation reduction failed"
  end if
  if (.not. global_pass) then
    if (rank == 0) write(*, '(a)') &
      "test_mpi_amr_transport_context_3d: FAIL"
    call MPI_Finalize(ierr)
    error stop "MPI transport context validation test failed"
  end if
  if (rank == 0) write(*, '(a)') &
    "test_mpi_amr_transport_context_3d: PASS"

  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

contains

  subroutine reset_context()
    config = baseline_config
    species = baseline_species
    transport = baseline_transport
  end subroutine reset_context

  subroutine expect_context(expected, label)
    logical, intent(in) :: expected
    character(len=*), intent(in) :: label

    logical :: context_ok

    call validate_mpi_amr_transport_context_3d( &
      MPI_COMM_WORLD, config, species, transport, context_ok)
    if (context_ok .neqv. expected) then
      write(*, '(a,i0,2a,l2)') &
        "rank ", rank, " unexpected transport context result for ", &
        trim(label), context_ok
      local_pass = .false.
    end if
  end subroutine expect_context

end program test_mpi_amr_transport_context_3d
