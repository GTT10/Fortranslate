module mpi_reactive_eos_policy_mod
  use mpi_f08
  use nasa7_thermo_mod, only: nasa7_species, nasa7_eos_policy_identity
  implicit none
  private
  public :: agree_mpi_reactive_eos_1d

contains

  subroutine agree_mpi_reactive_eos_1d(species, comm, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok

    character(len=128) :: identity, reference
    logical :: local_ok
    integer :: ierr, rank

    identity = nasa7_eos_policy_identity(species)
    local_ok = trim(identity) /= "PELEF_INVALID_EOS_POLICY"
    ok = .false.
    call MPI_Comm_rank(comm, rank, ierr)
    if (ierr /= MPI_SUCCESS) return
    reference = ""
    if (rank == 0) reference = identity
    call MPI_Bcast(reference, len(reference), MPI_CHARACTER, 0, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    local_ok = local_ok .and. identity == reference
    call MPI_Allreduce(local_ok, ok, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    ok = ierr == MPI_SUCCESS .and. ok
  end subroutine agree_mpi_reactive_eos_1d
end module mpi_reactive_eos_policy_mod
