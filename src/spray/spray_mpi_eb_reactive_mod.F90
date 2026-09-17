module spray_mpi_eb_reactive_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_mpi_amr_step_mod, only: advance_mpi_amr_mixture_coupling
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_full_3d
  implicit none
  private
  public :: advance_mpi_eb_spray_reactive
contains
  ! S/2 - existing EB reaction/transport/hydro - S/2. All publication is atomic.
  ! Replicated tables/options and a consistent geometry/nodal field are caller contracts.
  ! Ledgers use change-in-system signs, including the particle wall impulse.
  subroutine advance_mpi_eb_spray_reactive(comm,tree,owners,geometry,phi,species,reactions,transport, &
      liquids,options,state,temperature,parcels,component_mass,lower,spacing,dt,solver, &
      chemistry,transport_enabled,viscosity,conduction,diffusion,barodiffusion,rtol,atol, &
      particle_steps,minimum_theta,boundary_change,wall_change,chemical_change,ok,integrator,subcycling)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),spacing(3),dt,rtol,atol
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    character(len=*), intent(in) :: solver
    logical, intent(in) :: chemistry,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    integer, intent(out) :: particle_steps
    real(dp), intent(inout) :: minimum_theta,boundary_change(:),wall_change(:),chemical_change(:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    logical, optional, intent(in) :: subcycling
    type(amr_tree_state_3d) :: gas,temp
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    real(dp) :: chemical(reactive_nvar(size(species))),first(reactive_nvar(size(species)))
    real(dp) :: second(reactive_nvar(size(species))),theta
    integer :: first_steps,second_steps,ierr,nv
    logical :: valid,agreed
    ok=.false.;particle_steps=0;nv=size(boundary)
    valid=allocated(parcels).and.allocated(component_mass).and.ieee_is_finite(dt)
    if(valid) valid=dt>0.and.dt<=1e6_dp
    valid=valid.and.size(boundary_change)==nv.and.size(wall_change)==nv.and.size(chemical_change)==nv
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    boundary=0;wall=0;chemical=0;first=0;second=0;theta=1
    call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
      drops,masses,lower,spacing,.5_dp*dt,first_steps,valid,boundary_change=first, &
      eb_geometry=geometry,nodal_phi=phi)
    if(.not.valid) return
    call advance_mpi_eb_tree_full_3d(species,reactions,transport,tree,owners,geometry,solver,dt, &
      chemistry,rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion, &
      gas,temp,boundary,wall,chemical,theta,comm,valid,integrator,subcycling)
    if(.not.valid) return
    call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
      drops,masses,lower,spacing,.5_dp*dt,second_steps,valid,boundary_change=second, &
      eb_geometry=geometry,nodal_phi=phi)
    if(.not.valid) return
    state=gas;temperature=temp
    call move_alloc(drops,parcels)
    call move_alloc(masses,component_mass)
    boundary_change=boundary;wall_change=wall+first+second;chemical_change=chemical
    minimum_theta=theta;particle_steps=first_steps+second_steps;ok=.true.
  end subroutine
end module
