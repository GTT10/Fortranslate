module spray_mpi_eb_reactive_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options,valid_spray_options,valid_spray_parcel
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_mpi_amr_step_mod, only: advance_mpi_amr_mixture_coupling
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_full_3d
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options
  use spray_mpi_amr_film_step_mod, only: advance_mpi_amr_film_coupling
  use spray_mpi_mod, only: split_mpi_spray_parcels
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  implicit none
  private
  public :: advance_mpi_eb_spray_reactive
contains
  ! S/2 - F/2 - EB reaction/transport/hydro - F/2 - S/2; atomic publication.
  ! Supplying last_id with film state enables EB impacts in both particle halves.
  ! Replicated tables/options and a consistent geometry/nodal field are caller contracts.
  ! Ledgers use change-in-system signs, including the particle wall impulse.
  subroutine advance_mpi_eb_spray_reactive(comm,tree,owners,geometry,phi,species,reactions,transport, &
      liquids,options,state,temperature,parcels,component_mass,lower,spacing,dt,solver, &
      chemistry,transport_enabled,viscosity,conduction,diffusion,barodiffusion,rtol,atol, &
      particle_steps,minimum_theta,boundary_change,wall_change,chemical_change,ok,integrator,subcycling, &
      films,film_mass,film_settings,film_steps,last_id,model_rejected,split_created,particle_outflow, &
      acceleration,source_change,inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry, &
      boundary_faces,inflow_state)
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
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    real(dp), optional, intent(in) :: difmag
    real(dp) :: coefficient,low_coefficient,high_coefficient
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    logical, optional, intent(in) :: subcycling
    type(spray_wall_film), allocatable, optional, intent(inout) :: films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    integer, optional, intent(out) :: film_steps
    integer(int64), optional, intent(inout) :: last_id
    logical, optional, intent(out) :: model_rejected
    integer(int64), optional, intent(out) :: split_created
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(inout) :: source_change(:)
    real(dp) :: accel(3),low_accel(3),high_accel(3),force_ledger(reactive_nvar(size(species)))
    logical, optional, intent(in) :: particle_outflow
    character(len=32), allocatable :: particle_faces(:)
    logical :: outflow
    real(dp), allocatable :: half_exterior(:)
    real(dp) :: exterior(reactive_nvar(size(species)))
    integer(int64) :: before_split
    real(dp) :: low_multiplicity,high_multiplicity
    integer(int64) :: watermark
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    real(dp) :: first_film(reactive_nvar(size(species))),second_film(reactive_nvar(size(species)))
    integer :: flags,lowflags,highflags,first_film_steps,second_film_steps
    character(len=32) :: film_faces(6)
    type(amr_tree_state_3d) :: gas,temp
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    real(dp) :: chemical(reactive_nvar(size(species))),first(reactive_nvar(size(species)))
    real(dp) :: second(reactive_nvar(size(species))),theta
    integer :: first_steps,second_steps,ierr,nv
    logical :: valid,agreed
    ok=.false.;particle_steps=0;nv=size(boundary)
    if(present(film_steps)) film_steps=0
    if(present(split_created)) split_created=0
    if(present(model_rejected)) model_rejected=.false.
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    valid=ieee_is_finite(coefficient).and.coefficient>=0
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(coefficient,low_coefficient,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(coefficient,high_coefficient,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.low_coefficient/=high_coefficient) return
    accel=0;force_ledger=0
    if(present(acceleration)) accel=acceleration
    valid=all(ieee_is_finite(accel))
    valid=valid.and.all(abs(accel)<=1e6_dp)
    if(any(accel/=0)) valid=valid.and.present(source_change)
    if(present(source_change)) valid=valid.and.size(source_change)==nv
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(accel,low_accel,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(accel,high_accel,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(low_accel/=high_accel)) return
    outflow=.false.;exterior=0
    if(present(particle_outflow)) outflow=particle_outflow
    valid=.true.
    if(present(boundary_faces)) then
      valid=all(boundary_faces=='inflow'.or.boundary_faces=='outflow'.or.boundary_faces=='periodic')
      if(any(boundary_faces=='periodic')) then
        particle_faces=boundary_faces
        where(particle_faces=='inflow') particle_faces='outflow'
      end if
      valid=valid.and.outflow
    else
      valid=.not.present(inflow_state)
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    watermark=0
    if(present(last_id)) watermark=last_id
    flags=merge(1,0,present(films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))+ &
      8*merge(1,0,present(last_id))+16*merge(1,0,outflow)
    call MPI_Allreduce(flags,lowflags,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,highflags,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.lowflags/=highflags) return
    if(outflow) then
      if(flags/=16.and.flags/=24.and.flags/=31) return
      allocate(half_exterior(nv));half_exterior=0
      flags=flags-16
    end if
    if(flags/=0.and.flags/=7.and.flags/=8.and.flags/=15) return
    valid=valid_spray_options(options)
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    valid=valid.and.(options%maximum_multiplicity==0.or.present(last_id))
    if(present(films)) valid=valid.and.allocated(films).and.allocated(film_mass)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(options%maximum_multiplicity,low_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(options%maximum_multiplicity,high_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.low_multiplicity/=high_multiplicity) return
    if(present(films)) then
      patches=films;patch_mass=film_mass
    end if
    first_film=0;second_film=0;first_film_steps=0;second_film_steps=0;film_faces='slip_wall'
    valid=allocated(parcels).and.allocated(component_mass).and.ieee_is_finite(dt)
    if(valid) valid=dt>0.and.dt<=1e6_dp
    valid=valid.and.size(boundary_change)==nv.and.size(wall_change)==nv.and.size(chemical_change)==nv
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    boundary=0;wall=0;chemical=0;first=0;second=0;theta=1
    if(any(accel/=0)) then
      call force_kick()
      if(.not.valid) return
    end if
    call particle_half(first_steps,first)
    if(.not.valid) return
    if(present(films)) then
      call advance_mpi_amr_film_coupling(comm,species,transport,liquids,options,film_settings,gas,temp, &
        tree,owners,patches,patch_mass,lower,spacing,.5_dp*dt,film_faces,first_film_steps,first_film,valid,geometry,phi)
      if(.not.valid) return
    end if
    call advance_mpi_eb_tree_full_3d(species,reactions,transport,tree,owners,geometry,solver,dt, &
      chemistry,rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion, &
      gas,temp,boundary,wall,chemical,theta,comm,valid,integrator,subcycling, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      difmag=coefficient,extended_geometry=extended_geometry,boundary_faces=boundary_faces,inflow_state=inflow_state)
    if(.not.valid) return
    if(present(films)) then
      call advance_mpi_amr_film_coupling(comm,species,transport,liquids,options,film_settings,gas,temp, &
        tree,owners,patches,patch_mass,lower,spacing,.5_dp*dt,film_faces,second_film_steps,second_film,valid,geometry,phi)
      if(.not.valid) return
    end if
    call particle_half(second_steps,second)
    if(.not.valid) return
    before_split=watermark
    if(options%maximum_multiplicity>0) then
      call split_mpi_spray_parcels(comm,drops,options%maximum_multiplicity,options%maximum_parcels,watermark,valid,masses)
      if(.not.valid) return
    end if
    if(any(accel/=0)) then
      call force_kick()
      if(.not.valid) return
    end if
    state=gas;temperature=temp
    call move_alloc(drops,parcels)
    call move_alloc(masses,component_mass)
    if(present(films)) then
      call move_alloc(patches,films);call move_alloc(patch_mass,film_mass)
    end if
    if(present(last_id)) last_id=watermark
    if(present(split_created)) split_created=watermark-before_split
    if(present(film_steps)) film_steps=first_film_steps+second_film_steps
    boundary_change=boundary+exterior;wall_change=wall+first+second+first_film+second_film;chemical_change=chemical
    if(present(source_change)) source_change=force_ledger
    minimum_theta=theta;particle_steps=first_steps+second_steps;ok=.true.
  contains
    subroutine force_kick()
      real(dp) :: local(nv),global(nv)
      call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,gas,nv,comm,valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temp,1,comm,valid)
      if(.not.valid) return
      call local_force_kick(local)
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      valid=ierr==MPI_SUCCESS.and.agreed
      if(.not.valid) return
      call MPI_Allreduce(local,global,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      valid=ierr==MPI_SUCCESS
      if(valid) force_ledger=force_ledger+global
    end subroutine
    subroutine local_force_kick(local)
      real(dp), intent(out) :: local(nv)
      real(dp) :: dv(3),mass,momentum(3),energy,volume,velocity(3),prim(nv),checked_t,sound
      integer :: rank,status,node,i,j,k,p,n(3)
      valid=.false.;local=0;dv=.5_dp*dt*accel
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          mass=gas%fields(node)%values(1,i,j,k)
          if(.not.ieee_is_finite(mass).or.mass<=0) return
          velocity=gas%fields(node)%values(2:4,i,j,k)/mass
          if(any(abs(velocity+dv)>1e7_dp)) return
          momentum=mass*dv
          energy=dot_product(gas%fields(node)%values(2:4,i,j,k)+.5_dp*momentum,dv)
          gas%fields(node)%values(2:4,i,j,k)=gas%fields(node)%values(2:4,i,j,k)+momentum
          gas%fields(node)%values(5,i,j,k)=gas%fields(node)%values(5,i,j,k)+energy
          call reactive_conserved_to_primitive(species,gas%fields(node)%values(:,i,j,k), &
            temp%fields(node)%values(1,i,j,k),prim,checked_t,sound,valid)
          if(.not.valid) return
          valid=.false.
          if(abs(checked_t-temp%fields(node)%values(1,i,j,k))> &
            1e-9_dp*max(1._dp,abs(temp%fields(node)%values(1,i,j,k)))) return
          if(covered(tree,node,[i,j,k])) cycle
          volume=geometry(node)%volume_fraction(i,j,k)*geometry(node)%dx*geometry(node)%dy*geometry(node)%dz
          local(2:4)=local(2:4)+volume*momentum;local(5)=local(5)+volume*energy
        end do
        end do
        end do
      end do
      do p=1,size(drops)
        if(drops(p)%mass<=0) cycle
        mass=drops(p)%mass*drops(p)%multiplicity;velocity=drops(p)%velocity
        drops(p)%velocity=velocity+dv
        if(.not.valid_spray_parcel(drops(p))) return
        local(2:4)=local(2:4)+mass*dv
        local(5)=local(5)+mass*dot_product(velocity+.5_dp*dv,dv)
      end do
      valid=all(ieee_is_finite(local))
    end subroutine
    subroutine particle_half(steps,exchange)
      integer, intent(out) :: steps
      real(dp), intent(inout) :: exchange(:)
      logical :: rejected
      if(present(last_id).and.present(films)) then
        call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
          drops,masses,lower,spacing,.5_dp*dt,steps,valid,last_id=watermark,boundary_change=exchange, &
          wall_films=patches,film_mass=patch_mass,film_settings=film_settings, &
          eb_geometry=geometry,nodal_phi=phi,model_rejected=rejected,eb_outflow_change=half_exterior,boundary_faces=particle_faces)
        if(present(model_rejected)) model_rejected=rejected
      else if(present(last_id)) then
        call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
          drops,masses,lower,spacing,.5_dp*dt,steps,valid,last_id=watermark,boundary_change=exchange, &
          eb_geometry=geometry,nodal_phi=phi,eb_outflow_change=half_exterior,boundary_faces=particle_faces)
      else
        call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
          drops,masses,lower,spacing,.5_dp*dt,steps,valid,boundary_change=exchange, &
          eb_geometry=geometry,nodal_phi=phi,eb_outflow_change=half_exterior,boundary_faces=particle_faces)
      end if
      if(valid.and.outflow) exterior=exterior+half_exterior
    end subroutine
  end subroutine
end module
