module spray_mpi_amr_reactive_mod
  use les_filter_3d_mod, only: les_filter_weights
  use les_mpi_amr_tree_mod, only: advance_mpi_amr_tree_les_3d
  use spray_mpi_mod, only: split_mpi_spray_parcels
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options,valid_spray_parcel,valid_spray_options
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options
  use spray_mpi_amr_film_step_mod, only: advance_mpi_amr_film_coupling
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use spray_mpi_amr_step_mod, only: advance_mpi_amr_mixture_coupling
  use mpi_amr_patch_tree_split_3d_mod, only: advance_mpi_amr_tree_split_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  implicit none
  private
  public :: advance_mpi_amr_spray_reactive
contains
  ! S/2 - F/2 - L/2 - (R/2 T/2 H T/2 R/2) - L/2 - F/2 - S/2 with atomic publication.
  ! Common tables/options and globally unique parcel IDs are caller contracts.
  subroutine advance_mpi_amr_spray_reactive(comm,tree,owners,species,reactions,transport,liquids,options, &
      state,temperature,parcels,component_mass,lower,spacing,dt,reconstruction,limiter,solver, &
      chemistry,viscosity,conduction,diffusion,barodiffusion,rtol,atol,particle_steps, &
      maximum_correction,minimum_theta,ok,last_id,chemistry_integrator,subcycling, &
      boundary_faces,boundary_change,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
      wall_films,film_mass,film_settings,film_steps,les_settings,les_model,outlet_pressure, &
      acceleration,source_change,split_created, &
      les_test_model,les_test_ratio,les_grid_ratio,les_harmonic_mean)
    integer(int64), optional, intent(out) :: split_created
    integer(int64) :: before_split
    real(dp) :: low_multiplicity,high_multiplicity
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(inout) :: source_change(:)
    real(dp) :: accel(3),low_accel(3),high_accel(3)
    real(dp), optional, intent(in) :: outlet_pressure(6)
    real(dp) :: pressure(6),pressure_min(6),pressure_max(6)
    ! Static LES settings are [eddy coefficient, PrT, ScT, CI]; dynamic is enabled by model.
    real(dp), optional, intent(in) :: les_settings(4)
    character(len=*), optional, intent(in) :: les_model,les_test_model
    integer, optional, intent(in) :: les_test_ratio,les_grid_ratio
    logical, optional, intent(in) :: les_harmonic_mean
    character(len=64) :: filter,reference_filter
    integer :: filter_controls(3),filter_low(3),filter_high(3)
    real(dp), allocatable :: weights(:)
    logical :: filter_ok
    real(dp) :: sgs(4),sgs_min(4),sgs_max(4),les_correction,les_max
    character(len=16) :: model,reference_model
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(nasa7_species), intent(in) :: species(:)
    real(dp) :: les_boundary(reactive_nvar(size(species))),les_ledger(reactive_nvar(size(species)))
    real(dp) :: force_ledger(reactive_nvar(size(species)))
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
    if (present(split_created)) split_created=0
    if (present(film_steps)) film_steps=0
    sgs=[0._dp,.7_dp,.7_dp,0._dp];model='smagorinsky';les_ledger=0;les_max=0
    if (present(les_settings)) sgs=les_settings
    if (present(les_model)) model=les_model
    valid=allocated(parcels).and.allocated(component_mass).and.ieee_is_finite(dt)
    if (valid) valid=dt>0.and.dt<=1e6_dp
    valid=valid.and.all(ieee_is_finite(sgs))
    valid=valid.and.sgs(1)>=0.and.sgs(1)<=2.and.sgs(4)>=0.and.sgs(4)<=1
    valid=valid.and.minval(sgs(2:3))>=1e-6_dp.and.maxval(sgs(2:3))<=1e6_dp
    valid=valid.and.(model=='smagorinsky'.or.model=='wale'.or.model=='vreman'.or.model=='dynamic')
    valid=valid.and.valid_spray_options(options)
    valid=valid.and.(options%maximum_multiplicity==0.or.present(last_id))
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    if (present(boundary_change)) valid=valid.and.size(boundary_change)==size(boundary)
    valid=valid.and.(present(wall_films).eqv.present(film_mass))
    valid=valid.and.(present(wall_films).eqv.present(film_settings))
    if (present(wall_films)) valid=valid.and.allocated(wall_films).and.present(last_id)
    if (present(film_mass)) valid=valid.and.allocated(film_mass)
    faces='periodic'
    if (present(boundary_faces)) then
      valid=valid.and.all(len_trim(boundary_faces)<=len(faces))
      faces=boundary_faces
    end if
    accel=0;force_ledger=0
    if (present(acceleration)) accel=acceleration
    valid=valid.and.all(ieee_is_finite(accel)).and.all(abs(accel)<=1e6_dp)
    if (any(accel/=0)) valid=valid.and.present(source_change)
    if (present(source_change)) valid=valid.and.size(source_change)==size(boundary)
    pressure=0
    if (present(outlet_pressure)) pressure=outlet_pressure
    valid=valid.and.all(ieee_is_finite(pressure))
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    valid=valid.and..not.any(pressure>0.and.faces/='outflow')
    filter='box';filter_controls=[2,1,0]
    if(present(les_test_model)) filter=les_test_model
    if(present(les_test_ratio)) filter_controls(1)=les_test_ratio
    if(present(les_grid_ratio)) filter_controls(2)=les_grid_ratio
    if(present(les_harmonic_mean)) filter_controls(3)=merge(1,0,les_harmonic_mean)
    if(model=='dynamic') then
      call les_filter_weights(filter,filter_controls(1),weights,filter_ok)
      valid=valid.and.filter_ok.and.filter_controls(2)>=1.and.filter_controls(2)<=64
      valid=valid.and.all(faces=='periodic'.or.faces=='inflow'.or.faces=='outflow'.or. &
        faces=='slip_wall'.or.faces=='no_slip_wall'.or.faces=='isothermal_wall')
    else
      if(present(les_test_model).or.present(les_test_ratio).or. &
        present(les_grid_ratio).or.present(les_harmonic_mean)) valid=.false.
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(options%maximum_multiplicity,low_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(options%maximum_multiplicity,high_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_multiplicity/=high_multiplicity) return
    call MPI_Allreduce(accel,low_accel,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(accel,high_accel,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.any(low_accel/=high_accel)) return
    call MPI_Allreduce(pressure,pressure_min,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(pressure,pressure_max,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(pressure_min/=pressure_max)) return
    call MPI_Allreduce(sgs,sgs_min,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(sgs,sgs_max,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(sgs_min/=sgs_max)) return
    call MPI_Allreduce(filter_controls,filter_low,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(filter_controls,filter_high,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(filter_low/=filter_high)) return
    reference_filter=filter
    call MPI_Bcast(reference_filter,len(reference_filter),MPI_CHARACTER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=filter==reference_filter
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    reference_model=model
    call MPI_Bcast(reference_model,len(reference_model),MPI_CHARACTER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=model==reference_model
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
    if (present(wall_films)) then
      call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
        drops,masses,lower,spacing,.5_dp*dt,first_steps,valid,candidate_id,faces,first_boundary, &
        patches,patch_mass,film_settings)
    else
      call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
        drops,masses,lower,spacing,.5_dp*dt,first_steps,valid,candidate_id,faces,first_boundary)
    end if
    if (.not.valid) return
    call force_kick()
    if (.not.valid) return
    call les_half_step()
    if (.not.valid) return
    boundary=0;correction=0;theta=1
    call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,gas,temp,spacing,dt,faces, &
      reconstruction,limiter,solver,chemistry,viscosity,conduction,diffusion,barodiffusion,rtol,atol, &
      boundary,correction,theta,comm,valid,chemistry_integrator=chemistry_integrator,subcycling=subcycling, &
      inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_temperature=wall_temperature, &
      wall_velocity=wall_velocity,outlet_pressure=pressure)
    if (.not.valid) return
    call les_half_step()
    if (.not.valid) return
    call force_kick()
    if (.not.valid) return
    if (present(wall_films)) then
      call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
        drops,masses,lower,spacing,.5_dp*dt,second_steps,valid,candidate_id,faces,second_boundary, &
        patches,patch_mass,film_settings)
    else
      call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
        drops,masses,lower,spacing,.5_dp*dt,second_steps,valid,candidate_id,faces,second_boundary)
    end if
    if (.not.valid) return
    if (present(wall_films)) then
      call advance_mpi_amr_film_coupling(comm,species,transport,liquids,options,film_settings,gas,temp, &
        tree,owners,patches,patch_mass,lower,spacing,.5_dp*dt,faces,film_second,film_boundary,valid)
      if (.not.valid) return
      film_ledger=film_ledger+film_boundary
    end if
    before_split=candidate_id
    if (options%maximum_multiplicity>0) then
      call split_mpi_spray_parcels(comm,drops,options%maximum_multiplicity,options%maximum_parcels, &
        candidate_id,valid,masses)
      if (.not.valid) return
    end if
    boundary=boundary+first_boundary+second_boundary+film_ledger+les_ledger
    valid=all(ieee_is_finite(boundary)).and.all(ieee_is_finite(force_ledger))
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    if (present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    if (present(film_steps)) film_steps=film_first+film_second
    if (present(split_created)) split_created=candidate_id-before_split
    if (present(last_id)) last_id=candidate_id
    if (present(boundary_change)) boundary_change=boundary
    if (present(source_change)) source_change=force_ledger
    maximum_correction=max(correction,les_max);minimum_theta=theta
    particle_steps=first_steps+second_steps;ok=.true.
  contains
    subroutine force_kick()
      real(dp) :: local(size(boundary)),global(size(boundary))
      valid=.true.
      if (all(accel==0)) return
      call validate_mpi_amr_tree_local_3d(tree,owners,gas,size(boundary),comm,valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temp,1,comm,valid)
      if(.not.valid) return
      call local_force_kick(local)
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      valid=ierr==MPI_SUCCESS.and.agreed
      if(.not.valid) return
      call MPI_Allreduce(local,global,size(boundary),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      valid=ierr==MPI_SUCCESS
      if(valid) force_ledger=force_ledger+global
    end subroutine
    subroutine local_force_kick(local)
      real(dp), intent(out) :: local(size(boundary))
      real(dp) :: dv(3),mass,momentum(3),energy,volume,velocity(3),prim(size(boundary)),checked_t,sound
      real(dp) :: h(3,size(tree%nodes))
      integer :: rank,status,node,i,j,k,p,n(3)
      valid=.false.;local=0;dv=.5_dp*dt*accel
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
      h(:,1)=spacing
      do node=2,size(tree%nodes)
        h(:,node)=h(:,tree%nodes(node)%parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      end do
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
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
          volume=product(h(:,node))
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
    subroutine les_half_step()
      valid=.true.
      if (sgs(1)==0.and.sgs(4)==0.and.model/='dynamic') return
      call advance_mpi_amr_tree_les_3d(tree,owners,species,gas,temp,spacing,.5_dp*dt,faces, &
        sgs(1),sgs(2),sgs(3),les_boundary,les_correction,comm,valid, &
        inflow_state=inflow_state,inflow_temperature=inflow_temperature, &
        wall_temperature=wall_temperature,wall_velocity=wall_velocity,model=model,ci=sgs(4),outlet_pressure=pressure, &
        test_model=les_test_model,test_ratio=les_test_ratio,grid_ratio=les_grid_ratio,harmonic_mean=les_harmonic_mean)
      if (.not.valid) return
      les_ledger=les_ledger+les_boundary;les_max=max(les_max,les_correction)
    end subroutine
  end subroutine
end module
