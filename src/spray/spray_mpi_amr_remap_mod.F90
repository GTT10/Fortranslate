module spray_mpi_amr_remap_mod
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use spray_parcel_mod, only: spray_parcel
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use mpi_amr_patch_tree_remap_3d_mod, only: remap_mpi_amr_tree_reactive_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels
  implicit none
  private
  public :: remap_mpi_amr_spray
contains
  ! Same physical root domain; refinement/coarsening and owner changes are atomic
  ! across gas, particles and the published topology. Common species are borrowed.
  subroutine remap_mpi_amr_spray(comm,species,tree,owners,state,temperature,parcels,component_mass, &
      new_tree,new_owners,lower,spacing,periodic,maximum_parcels,ok)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(inout) :: tree
    integer, allocatable, intent(inout) :: owners(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    type(amr_patch_tree_3d), intent(in) :: new_tree
    integer, intent(in) :: new_owners(:),maximum_parcels
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: gas,temp,checked
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    logical :: valid,agreed
    integer :: ierr
    ok=.false.
    valid=allocated(owners).and.allocated(parcels).and.allocated(component_mass)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    checked=temperature
    call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,state,checked,comm,valid)
    if (.not.valid) return
    drops=parcels;masses=component_mass
    call migrate_amr_tree_spray_parcels(comm,new_tree,new_owners,drops,lower,spacing,periodic, &
      maximum_parcels,valid,masses)
    if (.not.valid) return
    call remap_mpi_amr_tree_reactive_3d(tree,owners,state,temperature,new_tree,new_owners, &
      species,gas,temp,comm,valid)
    if (.not.valid) return
    tree=new_tree;owners=new_owners
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    ok=.true.
  end subroutine
end module
