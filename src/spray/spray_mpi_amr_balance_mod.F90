module spray_mpi_amr_balance_mod
  use mpi_f08
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use mpi_amr_patch_tree_distribution_3d_mod, only: rebalance_mpi_amr_tree_3d
  use spray_parcel_mod, only: spray_parcel
  use spray_mixture_mod, only: valid_component_inventory
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels
  implicit none
  private
  public :: rebalance_mpi_amr_spray
contains
  subroutine rebalance_mpi_amr_spray(comm,tree,owners,nvar,state,temperature,parcels,component_mass, &
      lower,spacing,periodic,local_cost,minimum_gain,maximum_parcels,changed,ok)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(inout) :: owners(:)
    integer, intent(in) :: nvar,maximum_parcels
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),local_cost(:),minimum_gain
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: changed,ok
    type(amr_tree_state_3d) :: gas,temp
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    integer, allocatable :: assigned(:)
    logical :: valid,agreed,moved
    integer :: ierr,p
    ok=.false.;changed=.false.
    valid=allocated(parcels).and.allocated(component_mass).and.maximum_parcels>=0.and.maximum_parcels<=1000000
    call agree();if (.not.agreed) return
    valid=size(component_mass,2)==size(parcels).and.size(component_mass,1)>=1.and.size(component_mass,1)<=16
    call agree();if (.not.agreed) return
    do p=1,size(parcels)
      valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p))
    end do
    call agree();if (.not.agreed) return
    gas=state;temp=temperature;assigned=owners
    call rebalance_mpi_amr_tree_3d(tree,assigned,nvar,gas,temp,local_cost,minimum_gain,comm,moved,valid)
    if (.not.valid) return
    if (.not.moved) then
      ok=.true.;return
    end if
    drops=parcels;masses=component_mass
    call migrate_amr_tree_spray_parcels(comm,tree,assigned,drops,lower,spacing,periodic,maximum_parcels,valid,masses)
    if (.not.valid) return
    owners=assigned
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    changed=.true.;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
end module
