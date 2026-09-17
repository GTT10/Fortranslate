module spray_mpi_amr_reactive_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options
  use spray_mpi_amr_film_step_mod, only: advance_mpi_amr_film_coupling
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use spray_mpi_amr_step_mod, only: advance_mpi_amr_mixture_coupling
  use mpi_amr_patch_tree_split_3d_mod, only: advance_mpi_amr_tree_split_3d
  implicit none
  private
  public :: advance_mpi_amr_spray_reactive
contains
  ! S/2 - (R/2 T/2 H T/2 R/2) - S/2 with atomic publication.
  ! Common tables/options and globally unique parcel IDs are caller contracts.
  subroutine advance_mpi_amr_spray_reactive(comm,tree,owners,species,reactions,transport,liquids,options, &
      state,temperature,parcels,component_mass,lower,spacing,dt,reconstruction,limiter,solver, &
      chemistry,viscosity,conduction,diffusion,barodiffusion,rtol,atol,particle_steps, &
      maximum_correction,minimum_theta,ok,last_id,chemistry_integrator,subcycling, &
      boundary_faces,boundary_change,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
      wall_films,film_mass,film_settings,film_steps)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),dt,rtol,atol
    character(len=*), intent(in) :: reconstruction,limiter,solver
    logical, intent(in) :: chemistry,viscosity,conduction,diffusion,barodiffusion
    integer, intent(out) :: particle_steps
    real(dp), intent(inout) :: maximum_correction,minimum_theta
    logical, intent(out) :: ok
    integer(int64), optional, intent(inout) :: last_id
    character(len=*), optional, intent(in) :: chemistry_integrator
    logical, optional, intent(in) :: subcycling
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(inout) :: boundary_change(:)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    integer, optional, intent(out) :: film_steps
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    real(dp) :: film_boundary(reactive_nvar(size(species))),film_ledger(reactive_nvar(size(species)))
    integer :: film_first,film_second,film_flag,low_flag,high_flag
    type(amr_tree_state_3d) :: gas,temp
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    real(dp) :: boundary(reactive_nvar(size(species))),correction,theta
    real(dp) :: first_boundary(reactive_nvar(size(species))),second_boundary(reactive_nvar(size(species)))
    integer(int64) :: candidate_id
    integer :: first_steps,second_steps,ierr
    logical :: valid,agreed
    character(len=32) :: faces(6)
    ok=.false.;particle_steps=0;film_first=0;film_second=0;film_ledger=0
    if (present(film_steps)) film_steps=0
    valid=allocated(parcels).and.allocated(component_mass).and.ieee_is_finite(dt)
    if (valid) valid=dt>0.and.dt<=1e6_dp
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    if (present(boundary_change)) valid=valid.and.size(boundary_change)==size(boundary)
    valid=valid.and.(present(wall_films).eqv.present(film_mass))
    valid=valid.and.(present(wall_films).eqv.present(film_settings))
    if (present(wall_films)) valid=valid.and.allocated(wall_films)
    if (present(film_mass)) valid=valid.and.allocated(film_mass)
    faces='periodic'
    if (present(boundary_faces)) then
      valid=valid.and.all(len_trim(boundary_faces)<=len(faces))
      faces=boundary_faces
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    film_flag=merge(1,0,present(wall_films))
    call MPI_Allreduce(film_flag,low_flag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,high_flag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_flag/=high_flag) return
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    if (present(wall_films)) then
      patches=wall_films;patch_mass=film_mass
      call advance_mpi_amr_film_coupling(comm,species,transport,liquids,options,film_settings,gas,temp, &
        tree,owners,patches,patch_mass,lower,spacing,.5_dp*dt,faces,film_first,film_boundary,valid)
      if (.not.valid) return
      film_ledger=film_boundary
    end if
    candidate_id=0
    if (present(last_id)) candidate_id=last_id
    call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
      drops,masses,lower,spacing,.5_dp*dt,first_steps,valid,candidate_id,faces,first_boundary)
    if (.not.valid) return
    boundary=0;correction=0;theta=1
    call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,gas,temp,spacing,dt,faces, &
      reconstruction,limiter,solver,chemistry,viscosity,conduction,diffusion,barodiffusion,rtol,atol, &
      boundary,correction,theta,comm,valid,chemistry_integrator=chemistry_integrator,subcycling=subcycling, &
      inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_temperature=wall_temperature, &
      wall_velocity=wall_velocity)
    if (.not.valid) return
    call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
      drops,masses,lower,spacing,.5_dp*dt,second_steps,valid,candidate_id,faces,second_boundary)
    if (.not.valid) return
    if (present(wall_films)) then
      call advance_mpi_amr_film_coupling(comm,species,transport,liquids,options,film_settings,gas,temp, &
        tree,owners,patches,patch_mass,lower,spacing,.5_dp*dt,faces,film_second,film_boundary,valid)
      if (.not.valid) return
      film_ledger=film_ledger+film_boundary
    end if
    boundary=boundary+first_boundary+second_boundary+film_ledger
    valid=all(ieee_is_finite(boundary))
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    if (present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    if (present(film_steps)) film_steps=film_first+film_second
    if (present(last_id)) last_id=candidate_id
    if (present(boundary_change)) boundary_change=boundary
    maximum_correction=correction;minimum_theta=theta
    particle_steps=first_steps+second_steps;ok=.true.
  end subroutine
end module
