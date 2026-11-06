module spray_mpi_amr_remap_mod
  use mpi_f08
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use spray_mpi_amr_step_mod, only: convert_mpi_amr_tree_rotating_energy
  use nasa7_thermo_mod, only: nasa7_species
  use spray_parcel_mod, only: spray_parcel
  use spray_wall_impact_mod, only: spray_wall_film
  use spray_mpi_amr_film_mod, only: migrate_mpi_amr_wall_films
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use mpi_amr_patch_tree_remap_3d_mod, only: remap_mpi_amr_tree_reactive_3d,remap_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels
  implicit none
  private
  public :: remap_mpi_amr_spray
contains
  ! Same physical root domain; refinement/coarsening and owner changes are atomic
  ! across gas, particles and the published topology. Common species are borrowed.
  subroutine remap_mpi_amr_spray(comm,species,tree,owners,state,temperature,parcels,component_mass, &
      new_tree,new_owners,lower,spacing,periodic,maximum_parcels,ok,wall_films,film_mass,maximum_films, &
      reaction_source,heat_release,frame,internal_energy)
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    type(rotating_frame_3d), optional, intent(in) :: frame
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
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    integer, optional, intent(in) :: maximum_films
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    integer :: film_flag,low_flag,high_flag
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source,heat_release
    type(amr_tree_state_3d), allocatable :: old_rates,old_heat,rates,heat,old_held,held
    type(amr_tree_state_3d) :: gas,temp,checked,ordinary
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    logical :: valid,agreed
    integer :: ierr
    ok=.false.
    valid=allocated(owners).and.allocated(parcels).and.allocated(component_mass)
    valid=valid.and.(present(wall_films).eqv.present(film_mass))
    valid=valid.and.(present(wall_films).eqv.present(maximum_films))
    if (present(wall_films)) valid=valid.and.allocated(wall_films)
    if (present(film_mass)) valid=valid.and.allocated(film_mass)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    film_flag=merge(1,0,present(wall_films))+2*merge(1,0,present(reaction_source))+ &
      4*merge(1,0,present(heat_release))+8*merge(1,0,present(frame))+ &
      16*merge(1,0,present(internal_energy))
    call MPI_Allreduce(film_flag,low_flag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,high_flag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_flag/=high_flag) return
    if(present(internal_energy)) then
      allocate(old_held,held);old_held=internal_energy
    endif
    checked=temperature
    if(present(frame)) then
      ordinary=state
      call convert_mpi_amr_tree_rotating_energy(comm,tree,owners,ordinary,lower,spacing,frame,.false.,valid)
      if(.not.valid) return
      call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,ordinary,checked,comm,valid,internal_energy=old_held)
    else
      call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,state,checked,comm,valid,internal_energy=old_held)
    endif
    if (.not.valid) return
    drops=parcels;masses=component_mass
    call migrate_amr_tree_spray_parcels(comm,new_tree,new_owners,drops,lower,spacing,periodic, &
      maximum_parcels,valid,masses)
    if (.not.valid) return
    if(present(reaction_source)) then
      allocate(old_rates,rates);old_rates=reaction_source
    end if
    if(present(heat_release)) then
      allocate(old_heat,heat);old_heat=heat_release
    end if
    if(present(frame)) then
      ! Preserve rotating conserved energy through refinement/coarsening.
      ! Converting before remap would conserve a different energy variable.
      call remap_mpi_amr_tree_3d(tree,owners,state,new_tree,new_owners,size(species)+5,gas,comm,valid)
      if(.not.valid) return
      call remap_mpi_amr_tree_3d(tree,owners,checked,new_tree,new_owners,1,temp,comm,valid)
      if(.not.valid) return
      if(present(reaction_source)) then
        call remap_mpi_amr_tree_3d(tree,owners,old_rates,new_tree,new_owners,size(species)+1,rates,comm,valid)
        if(.not.valid) return
      endif
      if(present(heat_release)) then
        call remap_mpi_amr_tree_3d(tree,owners,old_heat,new_tree,new_owners,1,heat,comm,valid)
        if(.not.valid) return
      endif
      if(present(internal_energy)) then
        call remap_mpi_amr_tree_3d(tree,owners,old_held,new_tree,new_owners,1,held,comm,valid)
        if(.not.valid) return
      endif
      ordinary=gas
      call convert_mpi_amr_tree_rotating_energy(comm,new_tree,new_owners,ordinary,lower,spacing,frame,.false.,valid)
      if(.not.valid) return
      call recover_mpi_amr_tree_temperatures_3d(new_tree,new_owners,species,ordinary,temp,comm,valid,internal_energy=held)
    else
      call remap_mpi_amr_tree_reactive_3d(tree,owners,state,temperature,new_tree,new_owners, &
        species,gas,temp,comm,valid,old_rates,rates,old_heat,heat,old_held,held)
    endif
    if (.not.valid) return
    if (present(wall_films)) then
      patches=wall_films;patch_mass=film_mass
      call migrate_mpi_amr_wall_films(comm,new_tree,new_owners,patches,patch_mass,lower,spacing,maximum_films,valid)
      if (.not.valid) return
    end if
    if(present(reaction_source)) reaction_source=rates
    if(present(heat_release)) heat_release=heat
    if(present(internal_energy)) call move_alloc(held%fields,internal_energy%fields)
    tree=new_tree;owners=new_owners
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    if (present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    ok=.true.
  end subroutine
end module
