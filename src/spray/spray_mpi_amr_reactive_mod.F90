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
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use spray_mpi_amr_step_mod, only: advance_mpi_amr_mixture_coupling
  use mpi_amr_patch_tree_split_3d_mod, only: advance_mpi_amr_tree_split_3d
  implicit none
  private
  public :: advance_mpi_amr_spray_reactive
contains
  ! Periodic S/2 - (R/2 T/2 H T/2 R/2) - S/2 with atomic publication.
  ! Common tables/options and globally unique parcel IDs are caller contracts.
  subroutine advance_mpi_amr_spray_reactive(comm,tree,owners,species,reactions,transport,liquids,options, &
      state,temperature,parcels,component_mass,lower,spacing,dt,reconstruction,limiter,solver, &
      chemistry,viscosity,conduction,diffusion,barodiffusion,rtol,atol,particle_steps, &
      maximum_correction,minimum_theta,ok,last_id,chemistry_integrator,subcycling)
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
    type(amr_tree_state_3d) :: gas,temp
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    real(dp) :: boundary(reactive_nvar(size(species))),correction,theta
    integer(int64) :: candidate_id
    integer :: first_steps,second_steps,ierr
    logical :: valid,agreed
    character(len=8), parameter :: faces(6)='periodic'
    ok=.false.;particle_steps=0
    valid=allocated(parcels).and.allocated(component_mass).and.ieee_is_finite(dt)
    if (valid) valid=dt>0.and.dt<=1e6_dp
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    candidate_id=0
    if (present(last_id)) candidate_id=last_id
    call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
      drops,masses,lower,spacing,.5_dp*dt,first_steps,valid,candidate_id)
    if (.not.valid) return
    boundary=0;correction=0;theta=1
    call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,gas,temp,spacing,dt,faces, &
      reconstruction,limiter,solver,chemistry,viscosity,conduction,diffusion,barodiffusion,rtol,atol, &
      boundary,correction,theta,comm,valid,chemistry_integrator=chemistry_integrator,subcycling=subcycling)
    if (.not.valid) return
    call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
      drops,masses,lower,spacing,.5_dp*dt,second_steps,valid,candidate_id)
    if (.not.valid) return
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    if (present(last_id)) last_id=candidate_id
    maximum_correction=correction;minimum_theta=theta
    particle_steps=first_steps+second_steps;ok=.true.
  end subroutine
end module
