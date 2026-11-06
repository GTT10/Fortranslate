module spray_reactive_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use spray_wall_impact_mod, only: spray_wall_options,spray_wall_film
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use reactive_3d_mod, only: advance_reactive_chemistry_3d, advance_reactive_full_3d, recover_reactive_temperatures_3d
  use spray_parcel_mod
  use spray_coupling_3d_mod, only: advance_spray_coupling_3d, advance_mixture_coupling_3d, valid_spray_boundary_faces
  use les_smagorinsky_3d_mod, only: advance_les_ssprk2_3d
  use les_dynamic_smagorinsky_mod, only: advance_dynamic_les_ssprk2_3d
  use les_filter_3d_mod, only: les_filter_weights
  use les_sdc_3d_mod, only: advance_les_sdc_3d
  use constant_volume_reactor_mod, only: forced_reactor_integrator,advance_forced_reactor_interval
  use spray_sdc_3d_mod, only: spray_sdc_stage_3d,advance_spray_sdc_stage_3d,kick_serial_sdc_films
#ifdef PELEF_SPRAY_CVODE
  use sundials_forced_reactor_mod, only: advance_forced_reactor_cvode
  use sundials_constant_volume_reactor_mod, only: cvode_reactor_context, initialize_constant_volume_cvode, &
    advance_constant_volume_cvode, finalize_constant_volume_cvode
#endif
  use reactive_transport_3d_mod, only: valid_constant_transport_coefficients
  implicit none
  private
  public :: advance_spray_chemistry_3d, advance_spray_reactive_3d, spray_cvode_available
  public :: advance_spray_reactive_sdc_3d
contains
  subroutine advance_spray_reactive_sdc_3d(species,reactions,transport,components,options,state,temperature, &
      parcels,component_mass,lower,spacing,dt,time,iterations,chemistry,transport_flags,force,cs,prt,sct,rtol,atol,faces, &
      boundary_change,source_change,exchange_change,particle_wall_change,particle_boundary_change,particle_force_change,ok, &
      acceleration,evaluate_reaction,last_id,chemistry_method,model,ci,test_filter,test_ratio,grid_ratio,harmonic_mean, &
      inflow_state,inflow_temperature,inflow_mach,inflow_pressure,outlet_pressure,wall_temperature,wall_velocity, &
      difmag,constant_coefficients,breakup_properties,wall_options,wall_films,film_component_mass, &
      reaction_source,heat_release,mol_mode,mol_lax_friedrichs)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_options), optional, intent(in) :: breakup_properties(:)
    type(spray_wall_options), optional, intent(in) :: wall_options
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_component_mass(:,:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), optional, intent(inout) :: reaction_source(:,:,:,:),heat_release(:,:,:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),dt,time,force(3),cs,prt,sct,rtol,atol
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry,transport_flags(4)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(inout) :: boundary_change(:),source_change(:),exchange_change(:),particle_wall_change(5)
    real(dp), intent(inout) :: particle_boundary_change(:),particle_force_change(5)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: acceleration(3)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(forced_reactor_integrator), pointer :: reactor
    character(len=*), optional, intent(in) :: chemistry_method,model,test_filter
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    real(dp), optional, intent(in) :: ci,difmag,constant_coefficients(:)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: outlet_pressure(6),wall_temperature(:),wall_velocity(:,:)

    integer(int64), optional, intent(inout) :: last_id
    integer(int64), allocatable :: watermark
    type(spray_sdc_stage_3d) :: work
    type(spray_parcel), allocatable :: drops(:)
    type(spray_wall_film), allocatable :: films(:)
    real(dp), allocatable :: film_masses(:,:),film_rate(:,:,:,:),old_film_rate(:,:,:,:)
    real(dp) :: film_wall(5),old_film_wall(5)
    real(dp), allocatable :: masses(:,:),gas(:,:,:,:),temp(:,:,:),ir(:,:,:,:),heat(:,:,:)
    real(dp) :: boundary(size(state,1)),source(size(state,1)),exchange(size(state,1)),outside(size(state,1))
    real(dp) :: wall(5),forced(5),diffusion
    integer :: nv
    logical, optional, intent(in) :: mol_mode,mol_lax_friedrichs
    logical :: valid,mol
    integer :: particle_iterations,c
    real(dp), allocatable :: last_particle_rate(:,:,:,:),recovered(:,:,:)
    real(dp) :: new_exchange(size(state,1))
    ok=.false.;nv=size(state,1)
    mol=.false.;if(present(mol_mode)) mol=mol_mode
    particle_iterations=iterations;if(mol) particle_iterations=1
    allocate(last_particle_rate,mold=state);last_particle_rate=0
    if(.not.allocated(parcels).or..not.allocated(component_mass)) return
    if(present(wall_options).neqv.present(wall_films)) return
    if(present(wall_films).neqv.present(film_component_mass)) return
    if(present(wall_films)) then
      if(.not.allocated(wall_films).or..not.allocated(film_component_mass)) return
      films=wall_films;film_masses=film_component_mass
      allocate(film_rate,mold=state);allocate(old_film_rate,mold=state)
      film_rate=0;old_film_rate=0;film_wall=0;old_film_wall=0
    endif
    if(present(last_id)) then
      allocate(watermark);watermark=last_id
    endif
    if(size(boundary_change)/=nv.or.size(source_change)/=nv.or.size(exchange_change)/=nv.or. &
      size(particle_boundary_change)/=nv) return
    reactor=>implicit_reactor
    if(present(chemistry_method)) then
      select case(trim(chemistry_method))
      case('implicit')
      case('cvode')
#ifdef PELEF_SPRAY_CVODE
        reactor=>advance_forced_reactor_cvode
#else
        return
#endif
      case default
        return
      end select
    endif
    if(present(evaluate_reaction)) reactor=>evaluate_reaction
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    if(present(reaction_source)) then
      if(.not.chemistry.and..not.mol) return
      allocate(ir,source=reaction_source)
    endif
    if(present(heat_release)) allocate(heat,mold=heat_release)
    wall=0;forced=0;outside=0;exchange=0
    if(chemistry) then
      call advance_les_sdc_3d(species,gas,temp,spacing,dt,time,iterations,cs,prt,sct,boundary,source,diffusion,valid, &
        external_force=force,transport=transport,transport_flags=transport_flags,reactions=reactions, &
        chemistry_relative_tolerance=rtol,chemistry_absolute_tolerance=atol,boundary_faces=faces, &
        evaluate_reaction=reactor,auxiliary_source=particle_source,auxiliary_change=exchange, &
        model=model,ci=ci,test_filter=test_filter,test_ratio=test_ratio,grid_ratio=grid_ratio,harmonic_mean=harmonic_mean, &
        inflow_state=inflow_state,inflow_temperature=inflow_temperature,inflow_mach=inflow_mach, &
        inflow_pressure=inflow_pressure,outlet_pressure=outlet_pressure,wall_temperature=wall_temperature, &
        wall_velocity=wall_velocity,difmag=difmag,constant_coefficients=constant_coefficients, &
        reaction_source=ir,heat_release=heat,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs)
    else
      call advance_les_sdc_3d(species,gas,temp,spacing,dt,time,iterations,cs,prt,sct,boundary,source,diffusion,valid, &
        external_force=force,transport=transport,transport_flags=transport_flags,boundary_faces=faces, &
        auxiliary_source=particle_source,auxiliary_change=exchange, &
        model=model,ci=ci,test_filter=test_filter,test_ratio=test_ratio,grid_ratio=grid_ratio,harmonic_mean=harmonic_mean, &
        inflow_state=inflow_state,inflow_temperature=inflow_temperature,inflow_mach=inflow_mach, &
        inflow_pressure=inflow_pressure,outlet_pressure=outlet_pressure,wall_temperature=wall_temperature, &
        wall_velocity=wall_velocity,difmag=difmag,constant_coefficients=constant_coefficients, &
        reaction_source=ir,heat_release=heat,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs)
    endif
    if(.not.valid) return
    if(mol.and.chemistry.and.iterations>1) then
      ! Later gas passes omit the provisional new auxiliary family. Physical
      ! particle/film MK exchange has completed and must be applied once.
      gas=gas+.5_dp*dt*last_particle_rate
      if(.not.all(ieee_is_finite(gas)).or.any(gas(1,:,:,:)<=0).or.any(gas(6:,:,:,:)<0)) return
      allocate(recovered,mold=temp)
      call recover_reactive_temperatures_3d(species,gas,temp,size(gas,2),size(gas,3),size(gas,4),recovered,valid)
      if(.not.valid) return
      temp=recovered
      do c=1,nv
        new_exchange(c)=sum(last_particle_rate(c,:,:,:))*product(spacing)
      enddo
      exchange=exchange+.5_dp*dt*new_exchange
      if(.not.all(ieee_is_finite(exchange))) return
    endif
    if(work%next_stage/=0.or.any(gas(6:,:,:,:)<0)) return
    state=gas;temperature=temp;call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    if(present(reaction_source)) reaction_source=ir
    if(present(heat_release)) heat_release=heat
    boundary_change=boundary;source_change=source;exchange_change=exchange
    particle_wall_change=wall;particle_boundary_change=outside;particle_force_change=forced;ok=.true.
    if(present(last_id)) last_id=watermark
    if(present(wall_films)) then
      call move_alloc(films,wall_films);call move_alloc(film_masses,film_component_mass)
    endif
  contains
    subroutine implicit_reactor(species,reactions,initial_energy_density,external_species, &
        external_energy,duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps)
      type(nasa7_species), intent(in) :: species(:)
      type(elementary_reaction), intent(in) :: reactions(:)
      real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
      real(dp), intent(inout) :: species_density(:),temperature
      integer, intent(out) :: completed_steps
      logical, intent(out) :: ok
      integer, optional, intent(in) :: maximum_steps
      call advance_forced_reactor_interval(species,reactions,initial_energy_density,external_species, &
        external_energy,duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps=maximum_steps)
    end subroutine implicit_reactor
    subroutine particle_source(u,t,stage_time,interval,stage,rate,good)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),stage_time,interval
      integer, intent(in) :: stage
      real(dp), intent(out) :: rate(:,:,:,:)
      logical, intent(out) :: good
      good=.false.
      if(.not.ieee_is_finite(stage_time)) return
      if(present(wall_films)) then
        ! Existing films take the old kick before flight creates new films.
        if(stage==0) then
          call kick_serial_sdc_films(species,transport,components,options,wall_options,u,t,lower,spacing, &
            .5_dp*interval,faces,films,film_masses,film_rate,film_wall,good)
          if(.not.good) return
          old_film_rate=film_rate;old_film_wall=film_wall
        endif
        call advance_spray_sdc_stage_3d(species,transport,components,options,u,t,lower,spacing,time,interval,stage, &
          particle_iterations,faces,drops,masses,work,rate,wall,outside,forced,good,acceleration,watermark,breakup_properties, &
          wall_options,films,film_masses)
      else
        call advance_spray_sdc_stage_3d(species,transport,components,options,u,t,lower,spacing,time,interval,stage, &
          particle_iterations,faces,drops,masses,work,rate,wall,outside,forced,good,acceleration,watermark,breakup_properties)
      endif
      if(.not.good) return
      if(.not.present(wall_films)) then
        if(stage==1.and.mol) last_particle_rate=rate
        return
      endif
      if(stage==particle_iterations) then
        call kick_serial_sdc_films(species,transport,components,options,wall_options,u,t,lower,spacing, &
          .5_dp*interval,faces,films,film_masses,film_rate,film_wall,good)
        if(.not.good) return
        film_wall=film_wall+old_film_wall
      else
        film_rate=old_film_rate;film_wall=old_film_wall
      endif
      rate=rate+film_rate;wall=wall+film_wall
      if(stage==1.and.mol) last_particle_rate=rate
    end subroutine
  end subroutine advance_spray_reactive_sdc_3d

  logical function spray_cvode_available() result(available)
#ifdef PELEF_SPRAY_CVODE
    available=.true.
#else
    available=.false.
#endif
  end function

  subroutine advance_spray_chemistry_3d(species,reactions,state,temperature,interval,rtol,atol,method,ok,message)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: interval,rtol,atol
    character(len=*), intent(in) :: method
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    integer :: dims(3)
#ifdef PELEF_SPRAY_CVODE
    type(cvode_reactor_context) :: context
    real(dp), allocatable :: candidate(:,:,:,:),temp(:,:,:),check_t(:,:,:)
    real(dp) :: primitive(reactive_nprim(size(species))),y(size(species)),rho,e,t,sound,initial_step
    integer :: i,j,k
    logical :: local_ok,finalize_ok
    character(len=256) :: finalize_message
#endif
    ok=.false.; message='Invalid spray chemistry request'
    dims=shape(temperature)
    if (size(state,1) /= reactive_nvar(size(species))) return
    if (any(shape(state(1,:,:,:)) /= dims) .or. any(dims < 1)) return
    if (.not. all(ieee_is_finite([interval,rtol,atol]))) return
    if (interval < 0.0_dp .or. interval > 1.0e6_dp .or. min(rtol,atol) <= 0.0_dp) return
    select case(trim(method))
    case('implicit','explicit')
      call advance_reactive_chemistry_3d(species,reactions,state,temperature,dims(1),dims(2),dims(3), &
        interval,rtol,atol,ok,chemistry_integrator=method)
      if (ok) message=''
    case('cvode')
#ifndef PELEF_SPRAY_CVODE
      message='CVODE requested but this spray runtime was built without SUNDIALS'
#else
      allocate(candidate,source=state); allocate(temp,source=temperature); allocate(check_t,mold=temperature)
      call recover_reactive_temperatures_3d(species,state,temperature,dims(1),dims(2),dims(3),check_t,local_ok)
      if (.not. local_ok) return
      if (interval == 0.0_dp) then
        ok=.true.; message=''; return
      end if
      initial_step=min(interval,1.0e-12_dp)
      do k=1,dims(3)
        do j=1,dims(2)
          do i=1,dims(1)
            call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k),primitive,t,sound,local_ok)
            if (.not. local_ok) return
            rho=primitive(1); y=primitive(6:)
            e=state(5,i,j,k)/rho-0.5_dp*sum(primitive(2:4)**2)
            ! A new context after transport/spray uses the current rho and energy.
            ! Finalize every cell, including failed advances: no stale state or slot leak.
            call initialize_constant_volume_cvode(context,species,reactions,rho,e,0.0_dp,initial_step, &
              min(initial_step,interval*1.0e-12_dp),interval,rtol,atol,100000,y,t,local_ok,message)
            if (.not. local_ok) return
            call advance_constant_volume_cvode(context,interval,y,t,local_ok,message)
            call finalize_constant_volume_cvode(context,finalize_ok,finalize_message)
            if (.not. finalize_ok) then
              message=finalize_message; return
            end if
            if (.not. local_ok) return
            candidate(6:,i,j,k)=rho*y; temp(i,j,k)=t
          end do
        end do
      end do
      call recover_reactive_temperatures_3d(species,candidate,temp,dims(1),dims(2),dims(3),check_t,local_ok)
      if (.not. local_ok) return
      ! Chemistry never overwrites density, momentum or total energy.
      state=candidate; temperature=temp; ok=.true.; message=''
#endif
    case default
      message='Unknown chemistry backend (choose implicit, explicit, or cvode)'
    end select
  end subroutine

  subroutine advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,parcels, &
      lower,spacing,dt,chemistry_enabled,transport_enabled,chemistry_method,rtol,atol,cs,prt,sct,ok,message, &
      reconstruction,riemann_solver,boundary_faces,boundary_change,inflow_state,inflow_temperature, &
      wall_temperature,wall_velocity,acceleration,source_change,materials,breakup_options,wall_options,wall_films, &
      mixture_components,component_mass,film_component_mass,les_model,les_test_filter,les_test_ratio,les_grid_ratio, &
      les_harmonic_mean,les_ci,outlet_pressure,inflow_mach,inflow_pressure,difmag, &
      constant_coefficients)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquid
    type(spray_liquid), intent(in), optional :: materials(:)
    type(spray_liquid), intent(in), optional :: mixture_components(:)
    real(dp), allocatable, intent(inout), optional :: component_mass(:,:),film_component_mass(:,:)
    real(dp), allocatable :: mixture_mass(:,:),mixture_film_mass(:,:)
    type(spray_options), intent(in), optional :: breakup_options(:)
    type(spray_wall_options), intent(in), optional :: wall_options
    type(spray_wall_film), allocatable, intent(inout), optional :: wall_films(:)
    type(spray_wall_film), allocatable :: patches(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(in) :: lower(3),spacing(3),dt,rtol,atol,cs,prt,sct
    logical, intent(in) :: chemistry_enabled,transport_enabled
    character(len=*), optional, intent(in) :: les_model,les_test_filter
    integer, optional, intent(in) :: les_test_ratio,les_grid_ratio
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: les_ci,difmag
    real(dp) :: isotropic
    logical, optional, intent(in) :: les_harmonic_mean
    logical :: harmonic
    character(len=32) :: selected_les,selected_filter
    integer :: test_ratio,grid_ratio
    real(dp), allocatable :: filter_weights(:)
    character(len=*), intent(in) :: chemistry_method
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), allocatable :: candidate(:,:,:,:),temp(:,:,:)
    type(spray_parcel), allocatable :: drops(:)
    character(len=*), intent(in), optional :: reconstruction,riemann_solver
    character(len=32) :: selected_reconstruction,selected_riemann
    real(dp), intent(in), optional :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    character(len=*), intent(in), optional :: boundary_faces(6)
    real(dp), intent(out), optional :: boundary_change(:)
    real(dp), intent(in), optional :: acceleration(3)
    real(dp), intent(out), optional :: source_change(:)
    real(dp) :: force_ledger(size(state,1)),accel(3),coupled_accel(3),particle_force(5)
    character(len=32) :: faces(6)
    real(dp) :: ledger(size(state,1)),exchange(size(state,1))
    real(dp) :: theta
    integer :: dims(3),substeps,half
    logical :: local_ok
    ok=.false.; message='Invalid coupled spray step'
    if(present(constant_coefficients)) then
      if(.not.valid_constant_transport_coefficients(constant_coefficients)) return
    endif
    faces='periodic'; ledger=0;force_ledger=0;accel=0
    if (present(source_change)) then
      source_change=0
      if (size(source_change)/=size(state,1)) return
    end if
    if (present(acceleration)) accel=acceleration
    if (.not.all(ieee_is_finite(accel))) return
    if (maxval(abs(accel))>1e6_dp) return
    if ((any(accel/=0).or.options%fixed_parts).and..not.present(source_change)) then
      message='External or constraint source requires a ledger'; return
    end if
    if (present(boundary_faces)) faces=boundary_faces
    if(present(inflow_pressure).and..not.present(inflow_mach)) return
    if(present(inflow_mach)) then
      if(.not.all(ieee_is_finite(inflow_mach))) return
      if(any(inflow_mach>100).or.any(inflow_mach>=0.and.faces/='inflow')) return
      if(present(inflow_pressure)) then
        if(.not.all(ieee_is_finite(inflow_pressure))) return
        if(any(inflow_pressure<0).or.any(inflow_pressure>1e15_dp)) return
        if(any(inflow_pressure>0.and.inflow_mach<0)) return
      end if
    end if
    if (present(outlet_pressure)) then
      if (.not.all(ieee_is_finite(outlet_pressure))) return
      if (any(outlet_pressure<0).or.any(outlet_pressure>1e15_dp)) return
      if (any(outlet_pressure>0.and.faces/='outflow')) return
    end if
    if (present(boundary_change)) boundary_change=0
    selected_les='smagorinsky';selected_filter='box';test_ratio=2;grid_ratio=1;harmonic=.false.
    isotropic=0
    if (present(les_ci)) isotropic=les_ci
    if (.not.ieee_is_finite(isotropic)) return
    if (isotropic<0.or.isotropic>1) return
    if (present(les_harmonic_mean)) harmonic=les_harmonic_mean
    if (present(les_model)) then
      select case(les_model)
      case('smagorinsky','wale','vreman','dynamic')
        selected_les=les_model
      case default
        message='Invalid LES model'; return
      end select
    end if
    if (harmonic.and.selected_les/='dynamic') then
      message='Harmonic LES interpolation requires dynamic model';return
    end if
    if (selected_les=='dynamic'.and.isotropic/=0) then
      message='Fixed LES CI requires a static model';return
    end if
    if (selected_les=='dynamic') then
      if (present(les_test_filter)) selected_filter=les_test_filter
      if (present(les_test_ratio)) test_ratio=les_test_ratio
      if (present(les_grid_ratio)) grid_ratio=les_grid_ratio
      call les_filter_weights(selected_filter,test_ratio,filter_weights,local_ok)
      if (.not.local_ok.or.grid_ratio<1.or.grid_ratio>64) then
        message='Invalid dynamic LES filter';return
      end if
    end if
    if (.not.allocated(parcels)) return
    if (present(mixture_components).neqv.present(component_mass)) then
      message='Mixture components and mass storage must be supplied together';return
    end if
    if (present(component_mass)) then
      if (.not.allocated(component_mass)) return
      if (present(materials).or.present(breakup_options)) then
        message='Mixture and per-port pure-liquid property tables cannot be combined';return
      end if
      mixture_mass=component_mass
    end if
    if (.not.valid_spray_boundary_faces(faces)) return
    if (present(boundary_change)) then
      if (size(boundary_change)/=size(state,1)) return
    end if
    if (any(faces/='periodic')) then
      if (.not.present(boundary_change)) then
        message='Physical spray boundaries require a surface ledger'; return
      end if
    end if
    selected_reconstruction='characteristic_plm'; selected_riemann='hllc'
    if (present(reconstruction)) selected_reconstruction=reconstruction
    if (present(difmag)) then
      if (.not.ieee_is_finite(difmag)) return
      if (difmag<0) return
    end if
    if (present(riemann_solver)) selected_riemann=riemann_solver
    if (.not. all(ieee_is_finite([dt,rtol,atol,cs,prt,sct]))) return
    if (cs<0.or.cs>2.or.min(prt,sct)<1e-6_dp.or.max(prt,sct)>1e6_dp) return
    if (dt <= 0.0_dp .or. dt > 1.0e6_dp) return
    dims=shape(temperature)
    if (size(state,1)/=reactive_nvar(size(species)).or.any(shape(state(1,:,:,:))/=dims)) return
    if (any(spacing<=0).or..not.all(ieee_is_finite(spacing))) return
    allocate(patches(0))
    if (present(wall_options)) then
      if (any(wall_options%enabled).and..not.present(wall_films)) return
      if (any(wall_options%enabled).and.present(wall_velocity)) then
        if (any(wall_velocity/=0)) return
      end if
    end if
    if (present(wall_films)) then
      if (.not.allocated(wall_films)) return
      patches=wall_films
    end if
    if ((present(component_mass).and.present(wall_films)).neqv.present(film_component_mass)) then
      message='Mixed wall films require film component storage';return
    end if
    if (present(component_mass)) then
      allocate(mixture_film_mass(size(mixture_components),0))
      if (present(film_component_mass)) then
        if (.not.allocated(film_component_mass)) return
        mixture_film_mass=film_component_mass
      end if
    end if
    coupled_accel=0
    if(present(component_mass).and.options%mixture_time_scheme=='euler') coupled_accel=accel
    allocate(candidate,source=state); allocate(temp,source=temperature); allocate(drops,source=parcels)
    if (any(accel/=0)) then
      call force_kick(candidate,drops,0.5_dp*dt,local_ok)
      if (.not.local_ok) return
    end if
    ! Symmetric operator ordering; the nonlinear frozen-film parcel integrator
    ! remains first order in its substep. No second-order spray claim is implied.
    do half=1,2
      if (half == 1) then
        if (present(component_mass)) then
          call advance_mixture_coupling_3d(species,transport,mixture_components,options,candidate,temp,lower,spacing, &
            drops,mixture_mass,0.5_dp*dt,substeps,local_ok,message,faces,exchange, &
            wall_options,patches,mixture_film_mass,acceleration=coupled_accel,force_change=particle_force)
        else
        call advance_spray_coupling_3d(species,transport,liquid,options,candidate,temp,lower,spacing, &
          drops,0.5_dp*dt,substeps,local_ok,message,faces,exchange,materials,breakup_options,wall_options,patches, &
          force_change=particle_force)
        end if
        if (.not. local_ok) return
        force_ledger(:5)=force_ledger(:5)+particle_force
        ledger=ledger+exchange
        call les_half_step()
        if (.not. local_ok) then
          message='LES step rejected'; return
        end if
        ledger=ledger+exchange
      end if
      if (chemistry_enabled) then
        call advance_spray_chemistry_3d(species,reactions,candidate,temp,0.5_dp*dt,rtol,atol,chemistry_method,local_ok,message)
        if (.not. local_ok) return
      end if
      if (half == 1) then
        call advance_reactive_full_3d(species,reactions,transport,candidate,temp,dims(1),dims(2),dims(3), &
          spacing(1),spacing(2),spacing(3),dt,trim(selected_riemann),.false.,rtol,atol,transport_enabled, &
          .true.,.true.,.true.,.false.,theta,local_ok,trim(selected_reconstruction),'mc',chemistry_integrator='implicit', &
          boundary_faces=faces,boundary_change=exchange,inflow_state=inflow_state,inflow_temperature=inflow_temperature, &
          wall_temperature=wall_temperature,wall_velocity=wall_velocity,outlet_pressure=outlet_pressure, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure,difmag=difmag, &
          constant_coefficients=constant_coefficients)
        if (.not. local_ok) then
          message='Hydro/molecular transport step rejected'; return
        end if
        ledger=ledger+exchange
      else
        call les_half_step()
        if (.not. local_ok) then
          message='LES step rejected'; return
        end if
        ledger=ledger+exchange
        if (present(component_mass)) then
          call advance_mixture_coupling_3d(species,transport,mixture_components,options,candidate,temp,lower,spacing, &
            drops,mixture_mass,0.5_dp*dt,substeps,local_ok,message,faces,exchange, &
            wall_options,patches,mixture_film_mass,acceleration=coupled_accel,force_change=particle_force)
        else
        call advance_spray_coupling_3d(species,transport,liquid,options,candidate,temp,lower,spacing, &
          drops,0.5_dp*dt,substeps,local_ok,message,faces,exchange,materials,breakup_options,wall_options,patches, &
          force_change=particle_force)
        end if
        if (.not. local_ok) return
        force_ledger(:5)=force_ledger(:5)+particle_force
        ledger=ledger+exchange
      end if
    end do
    if (any(accel/=0)) then
      call force_kick(candidate,drops,0.5_dp*dt,local_ok)
      if (.not.local_ok) return
    end if
    if (present(source_change)) source_change=force_ledger
    if (present(boundary_change)) boundary_change=ledger
    state=candidate; temperature=temp; parcels=drops
    if (present(component_mass)) call move_alloc(mixture_mass,component_mass)
    if (present(film_component_mass)) call move_alloc(mixture_film_mass,film_component_mass)
    if (present(wall_films)) call move_alloc(patches,wall_films); ok=.true.; message=''
  contains
    subroutine les_half_step()
      if (selected_les=='dynamic') then
        exchange=0;local_ok=.true.
        ! cs=0 retains the existing disabled convention. A positive cs enables
        ! dynamic coefficients; its magnitude and fixed prt/sct do not set them.
        if (cs>0) call advance_dynamic_les_ssprk2_3d(species,candidate,temp,spacing,selected_filter, &
          test_ratio,.5_dp*dt,local_ok,grid_ratio,harmonic,boundary_faces=faces,boundary_change=exchange, &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_temperature=wall_temperature, &
          wall_velocity=wall_velocity,outlet_pressure=outlet_pressure, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      else
        call advance_les_ssprk2_3d(species,candidate,temp,spacing,cs,prt,sct,.5_dp*dt,local_ok,faces,exchange, &
          inflow_state,inflow_temperature,wall_temperature,wall_velocity,model=selected_les,ci=isotropic, &
          outlet_pressure=outlet_pressure, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      end if
    end subroutine

    subroutine force_kick(gas,particles,interval,accepted)
      real(dp), intent(inout) :: gas(:,:,:,:)
      type(spray_parcel), intent(inout) :: particles(:)
      real(dp), intent(in) :: interval
      logical, intent(out) :: accepted
      real(dp) :: dv(3),momentum(3),energy,mass,velocity(3),volume
      real(dp) :: checked_t(size(gas,2),size(gas,3),size(gas,4))
      integer :: i,j,k,p
      accepted=.false.;dv=interval*accel;volume=product(spacing)
      do k=1,size(gas,4)
        do j=1,size(gas,3)
          do i=1,size(gas,2)
            mass=gas(1,i,j,k)
            if (mass<=0.or..not.ieee_is_finite(mass)) return
            if (maxval(abs(gas(2:4,i,j,k)/mass+dv))>1e7_dp) return
            momentum=mass*dv
            energy=dot_product(gas(2:4,i,j,k)+0.5_dp*momentum,dv)
            gas(2:4,i,j,k)=gas(2:4,i,j,k)+momentum
            gas(5,i,j,k)=gas(5,i,j,k)+energy
            force_ledger(2:4)=force_ledger(2:4)+volume*momentum
            force_ledger(5)=force_ledger(5)+volume*energy
          end do
        end do
      end do
      do p=1,size(particles)
        if(options%fixed_parts) cycle
        if(present(component_mass).and.options%mixture_time_scheme=='euler') cycle
        if (particles(p)%mass<=0) cycle
        mass=particles(p)%mass*particles(p)%multiplicity
        velocity=particles(p)%velocity
        particles(p)%velocity=velocity+dv
        force_ledger(2:4)=force_ledger(2:4)+mass*dv
        force_ledger(5)=force_ledger(5)+mass*dot_product(velocity+0.5_dp*dv,dv)
        if (.not.valid_spray_parcel(particles(p))) return
      end do
      if (.not.all(ieee_is_finite(gas)).or..not.all(ieee_is_finite(force_ledger))) return
      call recover_reactive_temperatures_3d(species,gas,temp,size(gas,2),size(gas,3),size(gas,4),checked_t,accepted)
      if (.not.accepted) return
      accepted=all(abs(checked_t-temp)<=1e-9_dp*max(1.0_dp,abs(temp)))
    end subroutine
  end subroutine
end module spray_reactive_3d_mod
