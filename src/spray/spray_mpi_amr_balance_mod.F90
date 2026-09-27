module spray_mpi_amr_balance_mod
  use mpi_f08
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use mpi_amr_patch_tree_distribution_3d_mod, only: rebalance_mpi_amr_tree_3d, &
    rebalance_mpi_amr_tree_with_field_3d
  use spray_parcel_mod, only: spray_parcel
  use spray_wall_impact_mod, only: spray_wall_film,valid_wall_component_inventory
  use spray_mpi_amr_film_mod, only: migrate_mpi_amr_wall_films
  use spray_mixture_mod, only: valid_component_inventory
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels,migrate_eb_tree_spray_parcels
  implicit none
  private
  public :: rebalance_mpi_amr_spray
contains
  subroutine rebalance_mpi_amr_spray(comm,tree,owners,nvar,state,temperature,parcels,component_mass, &
      lower,spacing,periodic,local_cost,minimum_gain,maximum_parcels,changed,ok,wall_films,film_mass, &
      maximum_films,nodal_phi,reaction_source)
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
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
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    integer, optional, intent(in) :: maximum_films
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    integer :: film_flag,low_flag,high_flag
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    type(amr_tree_state_3d) :: gas,temp,history
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    integer, allocatable :: assigned(:)
    logical :: valid,agreed,moved
    integer :: ierr,p
    ok=.false.;changed=.false.
    valid=allocated(parcels).and.allocated(component_mass).and.maximum_parcels>=0.and.maximum_parcels<=1000000
    valid=valid.and.(present(wall_films).eqv.present(film_mass))
    valid=valid.and.(present(wall_films).eqv.present(maximum_films))
    if (present(wall_films)) valid=valid.and.allocated(wall_films)
    if (present(film_mass)) valid=valid.and.allocated(film_mass)
    call agree();if (.not.agreed) return
    film_flag=merge(1,0,present(wall_films))+2*merge(1,0,present(nodal_phi))+4*merge(1,0,present(reaction_source))
    call MPI_Allreduce(film_flag,low_flag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,high_flag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_flag/=high_flag) return
    if (present(wall_films)) then
      valid=size(film_mass,2)==size(wall_films).and.size(film_mass,1)>=1.and.size(film_mass,1)<=16
      valid=valid.and.maximum_films>=0.and.maximum_films<=1000000
      call agree();if (.not.agreed) return
      do p=1,size(wall_films)
        valid=valid.and.valid_wall_component_inventory(wall_films(p),film_mass(:,p))
      end do
      call agree();if (.not.agreed) return
    end if
    valid=size(component_mass,2)==size(parcels).and.size(component_mass,1)>=1.and.size(component_mass,1)<=16
    call agree();if (.not.agreed) return
    do p=1,size(parcels)
      valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p))
    end do
    call agree();if (.not.agreed) return
    gas=state;temp=temperature;assigned=owners
    if(present(reaction_source)) then
      history=reaction_source
      call rebalance_mpi_amr_tree_with_field_3d(tree,assigned,nvar,gas,temp,nvar-4,history, &
        local_cost,minimum_gain,comm,moved,valid)
    else
    call rebalance_mpi_amr_tree_3d(tree,assigned,nvar,gas,temp,local_cost,minimum_gain,comm,moved,valid)
    end if
    if (.not.valid) return
    if (.not.moved.and..not.present(nodal_phi)) then
      ok=.true.;return
    end if
    drops=parcels;masses=component_mass
    if(present(nodal_phi)) then
      call migrate_eb_tree_spray_parcels(comm,tree,assigned,nodal_phi,drops,lower,spacing, &
        maximum_parcels,valid,masses,periodic=periodic)
    else
      call migrate_amr_tree_spray_parcels(comm,tree,assigned,drops,lower,spacing,periodic,maximum_parcels,valid,masses)
    end if
    if (.not.valid) return
    if (present(wall_films)) then
      patches=wall_films;patch_mass=film_mass
      call migrate_mpi_amr_wall_films(comm,tree,assigned,patches,patch_mass,lower,spacing,maximum_films,valid,nodal_phi, &
        periodic=periodic)
      if (.not.valid) return
    end if
    if(.not.moved) then
      ok=.true.;return
    end if
    if(present(reaction_source)) call move_alloc(history%fields,reaction_source%fields)
    owners=assigned
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    if (present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    changed=.true.;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
end module
