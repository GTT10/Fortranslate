program test_mpi_amr_selected_context_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: &
    elementary_reaction, reaction_kind_elementary
  use gas_transport_mod, only: gas_transport_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use transport_database_mod, only: load_h2o2_full_transport
  use mpi_amr_reactive_3d_application_mod, only: &
    validate_mpi_amr_selected_context_3d
  implicit none

  type(nasa7_species), allocatable :: species(:), baseline_species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(elementary_reaction), allocatable :: baseline_reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(gas_transport_species), allocatable :: baseline_transport(:)
  real(dp), allocatable :: composition(:), baseline_composition(:)
  real(dp) :: relative_tolerance, absolute_tolerance
  character(len=64) :: bundle_sha256
  character(len=32) :: chemistry_integrator
  logical :: ok, local_pass, global_pass
  integer :: ierr, rank, reaction_index

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"

  call load_h2o2_full_thermo(baseline_species, ok)
  if (ok) call load_h2o2_full_mechanism(baseline_reactions, ok)
  if (ok) call load_h2o2_full_transport(baseline_transport, ok)
  if (.not. ok) then
    call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
    error stop "Failed to load selected-context fixture"
  end if
  allocate(baseline_composition(size(baseline_species)))
  baseline_composition = 1.0_dp / real(size(baseline_species), dp)
  local_pass = .true.

  call reset_context()
  call expect_context(.true., "baseline")

  call reset_context()
  if (rank == 1) bundle_sha256(64:64) = "b"
  call expect_context(.false., "bundle SHA")

  call reset_context()
  if (rank == 1) chemistry_integrator = "implicit"
  call expect_context(.false., "integrator")

  call reset_context()
  if (rank == 1) relative_tolerance = 2.0_dp * relative_tolerance
  call expect_context(.false., "chemistry tolerance")

  call reset_context()
  if (rank == 1) then
    composition(1) = composition(1) + 1.0e-12_dp
    composition(2) = composition(2) - 1.0e-12_dp
  end if
  call expect_context(.false., "composition")

  call reset_context()
  if (rank == 1) then
    species(1)%low_coefficients(1) = &
      nearest(species(1)%low_coefficients(1), 1.0_dp)
  end if
  call expect_context(.false., "NASA7")

  call reset_context()
  if (rank == 1) reactions(1)%equation(128:128) = "x"
  call expect_context(.false., "reaction equation")

  call reset_context()
  if (rank == 1) then
    reactions(1)%forward_rate%pre_exponential = nearest( &
      reactions(1)%forward_rate%pre_exponential, 1.0_dp)
  end if
  call expect_context(.false., "reaction rate")

  call reset_context()
  if (rank == 1) then
    transport(1)%diameter = nearest(transport(1)%diameter, 1.0_dp)
  end if
  call expect_context(.false., "transport")

  call reset_context()
  if (rank == 1) then
    do reaction_index = 1, size(reactions)
      if (reactions(reaction_index)%kind == reaction_kind_elementary) then
        allocate(reactions(reaction_index)%third_body_efficiencies(1))
        reactions(reaction_index)%third_body_efficiencies = 1.0_dp
        exit
      end if
    end do
  end if
  call expect_context(.false., "malformed third-body shape")

  call MPI_Allreduce( &
    local_pass, global_pass, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS .or. .not. global_pass) then
    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
    error stop "Selected MPI 3D AMR context validation failed"
  end if
  if (rank == 0) write(*, '(a)') &
    "test_mpi_amr_selected_context_3d: PASS"

  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

contains

  subroutine reset_context()
    bundle_sha256 = repeat("a", len(bundle_sha256))
    chemistry_integrator = "explicit"
    relative_tolerance = 2.0e-7_dp
    absolute_tolerance = 1.0e-12_dp
    composition = baseline_composition
    species = baseline_species
    reactions = baseline_reactions
    transport = baseline_transport
  end subroutine reset_context

  subroutine expect_context(expected, label)
    logical, intent(in) :: expected
    character(len=*), intent(in) :: label

    logical :: context_ok

    call validate_mpi_amr_selected_context_3d( &
      MPI_COMM_WORLD, bundle_sha256, chemistry_integrator, &
      relative_tolerance, absolute_tolerance, composition, species, &
      reactions, transport, context_ok)
    if (context_ok .neqv. expected) then
      write(*, '(a,i0,2a,l2)') &
        "rank ", rank, " unexpected selected context result for ", &
        trim(label), context_ok
      local_pass = .false.
    end if
  end subroutine expect_context

end program test_mpi_amr_selected_context_3d
