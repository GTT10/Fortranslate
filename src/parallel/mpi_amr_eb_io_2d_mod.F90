module mpi_amr_eb_io_2d_mod
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use eb_geometry_2d_mod, only: eb_geometry_2d
  use amr_eb_regrid_2d_mod, only: reactive_eb_patch_set_2d
  use simulation_config_reactive_eb_2d_mod, only: reactive_eb_2d_config
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config
  use reactive_eb_2d_driver_mod, only: write_reactive_eb_2d_csv
  use reactive_eb_amr_2d_driver_mod, only: &
    write_reactive_eb_amr_patch_set_2d_checkpoint
  use mpi_amr_eb_patch_2d_mod, only: &
    mpi_amr_eb_patch_distribution_2d, &
    mpi_amr_eb_sparse_patch_set_2d, &
    gather_sparse_owned_reactive_eb_patch_set_to_root_2d
  implicit none
  private

  public :: write_sparse_owned_reactive_eb_patch_set_2d_checkpoint
  public :: write_sparse_owned_reactive_eb_patch_set_2d_csv

contains

  subroutine write_sparse_owned_reactive_eb_patch_set_2d_checkpoint( &
      path, species, config, distribution, sparse_patch_set, &
      coarse_geometry, patch_set_template, root, time, steps, regrids, &
      minimum_dt, base_density, ok, local_transfers)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    integer, intent(in) :: root
    real(dp), intent(in) :: time, minimum_dt, base_density
    integer, intent(in) :: steps, regrids
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_transfers

    type(reactive_eb_patch_set_2d) :: materialized_patch_set
    real(dp), allocatable :: coarse_state(:, :, :)
    real(dp), allocatable :: coarse_temperature(:, :)
    logical :: materialized, write_ok
    integer :: ierr, transfers

    ok = .false.
    transfers = 0
    if (present(local_transfers)) local_transfers = 0
    call gather_sparse_owned_reactive_eb_patch_set_to_root_2d( &
      distribution, sparse_patch_set, coarse_geometry, patch_set_template, &
      root, coarse_state, coarse_temperature, materialized_patch_set, &
      materialized, transfers)
    if (.not. materialized) return

    write_ok = .true.
    if (distribution%rank == root) then
      call write_reactive_eb_amr_patch_set_2d_checkpoint( &
        path, species, config, coarse_state, coarse_temperature, &
        coarse_geometry, materialized_patch_set, time, steps, regrids, &
        minimum_dt, base_density, write_ok)
    end if
    call MPI_Bcast( &
      write_ok, 1, MPI_LOGICAL, root, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. write_ok) return

    ok = .true.
    if (present(local_transfers)) local_transfers = transfers
  end subroutine write_sparse_owned_reactive_eb_patch_set_2d_checkpoint

  subroutine write_sparse_owned_reactive_eb_patch_set_2d_csv( &
      species, config, distribution, sparse_patch_set, coarse_geometry, &
      patch_set_template, root, time, ok, local_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    integer, intent(in) :: root
    real(dp), intent(in) :: time
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_transfers

    type(reactive_eb_2d_config) :: child_config
    type(reactive_eb_patch_set_2d) :: materialized_patch_set
    real(dp), allocatable :: coarse_state(:, :, :)
    real(dp), allocatable :: coarse_temperature(:, :)
    character(len=1024) :: child_path
    logical :: materialized, path_ok, write_ok
    integer :: child, ierr, transfers

    ok = .false.
    transfers = 0
    if (present(local_transfers)) local_transfers = 0
    call gather_sparse_owned_reactive_eb_patch_set_to_root_2d( &
      distribution, sparse_patch_set, coarse_geometry, patch_set_template, &
      root, coarse_state, coarse_temperature, materialized_patch_set, &
      materialized, transfers)
    if (.not. materialized) return

    write_ok = .true.
    if (distribution%rank == root) then
      write_ok = config%multipatch_enabled .and. &
        .not. config%three_level_enabled .and. &
        len_trim(config%eb%flow%output_file) > 0 .and. &
        (materialized_patch_set%patch_count() == 0 .or. &
         (len_trim(config%fine_output_file) > 0 .and. &
          trim(config%fine_output_file) /= &
            trim(config%eb%flow%output_file)))
      if (write_ok) call write_reactive_eb_2d_csv( &
        config%eb%flow%output_file, species, config%eb, coarse_state, &
        coarse_temperature, coarse_geometry, time, write_ok)
      do child = 1, materialized_patch_set%patch_count()
        if (.not. write_ok) exit
        call make_patch_output_path_2d( &
          config%fine_output_file, child, child_path, path_ok)
        if (.not. path_ok) then
          write_ok = .false.
          exit
        end if
        child_config = config%eb
        child_config%flow%nx = &
          materialized_patch_set%children(child)%geometry%nx
        child_config%flow%ny = &
          materialized_patch_set%children(child)%geometry%ny
        child_config%flow%x_lower = &
          materialized_patch_set%children(child)%geometry%x_lower
        child_config%flow%x_upper = &
          materialized_patch_set%children(child)%geometry%x_upper
        child_config%flow%y_lower = &
          materialized_patch_set%children(child)%geometry%y_lower
        child_config%flow%y_upper = &
          materialized_patch_set%children(child)%geometry%y_upper
        child_config%flow%output_file = trim(child_path)
        call write_reactive_eb_2d_csv( &
          child_path, species, child_config, &
          materialized_patch_set%children(child)%state, &
          materialized_patch_set%children(child)%temperature, &
          materialized_patch_set%children(child)%geometry, time, write_ok)
      end do
    end if
    call MPI_Bcast( &
      write_ok, 1, MPI_LOGICAL, root, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. write_ok) return

    ok = .true.
    if (present(local_transfers)) local_transfers = transfers
  end subroutine write_sparse_owned_reactive_eb_patch_set_2d_csv

  subroutine make_patch_output_path_2d( &
      base_path, patch_index, output_path, ok)
    character(len=*), intent(in) :: base_path
    integer, intent(in) :: patch_index
    character(len=*), intent(out) :: output_path
    logical, intent(out) :: ok

    character(len=:), allocatable :: candidate
    character(len=16) :: index_text
    integer :: dot

    output_path = ""
    ok = .false.
    if (len_trim(base_path) == 0 .or. patch_index < 1) return
    write(index_text, '(i4.4)') patch_index
    dot = scan(trim(base_path), ".", back=.true.)
    if (dot > 1) then
      candidate = trim(base_path(:dot - 1)) // "_patch" // &
        trim(index_text) // trim(base_path(dot:))
    else
      candidate = trim(base_path) // "_patch" // trim(index_text) // ".csv"
    end if
    if (len(candidate) > len(output_path)) return
    output_path = candidate
    ok = .true.
  end subroutine make_patch_output_path_2d

end module mpi_amr_eb_io_2d_mod
