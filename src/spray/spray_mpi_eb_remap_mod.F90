module spray_mpi_eb_remap_mod
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_parcel_mod, only: spray_parcel
  use spray_wall_impact_mod, only: spray_wall_film
  use spray_mpi_amr_film_mod, only: migrate_mpi_amr_wall_films
  use mpi_eb_tree_remap_3d_mod, only: remap_mpi_eb_tree_3d
  use spray_mpi_amr_tree_mod, only: migrate_eb_tree_spray_parcels
  implicit none
  private
  public :: remap_mpi_eb_spray
contains
  ! Immutable physical geometry, common root/refinement lattice. All state commits together.
  subroutine remap_mpi_eb_spray(comm,species,tree,owners,geometry,state,temperature,parcels,components, &
      new_tree,new_owners,finest,phi,lower,spacing,capacity,ok,wall_films,film_mass,maximum_films,reaction_source)
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    type(amr_tree_state_3d), allocatable :: ir
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(inout) :: tree
    integer, allocatable, intent(inout) :: owners(:)
    type(eb_geometry_3d), allocatable, intent(inout) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: components(:,:)
    type(amr_patch_tree_3d), intent(in) :: new_tree
    integer, intent(in) :: new_owners(:),capacity
    type(eb_geometry_3d), intent(in) :: finest
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),spacing(3)
    logical, intent(out) :: ok
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    integer, optional, intent(in) :: maximum_films
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    type(amr_tree_state_3d) :: gas,temp
    type(eb_geometry_3d), allocatable :: candidate_geometry(:)
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    logical :: valid,agreed
    integer :: ierr,flags,lowflags,highflags
    ok=.false.
    flags=merge(1,0,present(wall_films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(maximum_films))
    call MPI_Allreduce(flags,lowflags,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,highflags,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.lowflags/=highflags) return
    if(flags/=0.and.flags/=7) return
    valid=allocated(owners).and.allocated(geometry).and.allocated(parcels).and.allocated(components)
    if(present(wall_films)) valid=valid.and.allocated(wall_films).and.allocated(film_mass)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    gas=state;temp=temperature;candidate_geometry=geometry;drops=parcels;masses=components
    if(present(reaction_source)) then
      allocate(ir);ir=reaction_source
    end if
    call remap_mpi_eb_tree_3d(species,tree,new_tree,owners,new_owners,finest,candidate_geometry,gas,temp,comm,valid,ir)
    if(.not.valid) return
    call migrate_eb_tree_spray_parcels(comm,new_tree,new_owners,phi,drops,lower,spacing,capacity,valid,masses)
    if(.not.valid) return
    if(present(wall_films)) then
      patches=wall_films;patch_mass=film_mass
      call migrate_mpi_amr_wall_films(comm,new_tree,new_owners,patches,patch_mass, &
        lower,spacing,maximum_films,valid,phi)
      if(.not.valid) return
    end if
    if(present(reaction_source)) reaction_source=ir
    tree=new_tree;owners=new_owners
    call move_alloc(candidate_geometry,geometry)
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(drops,parcels);call move_alloc(masses,components)
    if(present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    ok=.true.
  end subroutine
end module
