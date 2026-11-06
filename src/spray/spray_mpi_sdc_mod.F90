module spray_mpi_sdc_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use reactive_mol_3d_mod, only: advance_reactive_mol_3d,convert_regular_energy_3d
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use reactive_3d_mod, only: advance_reactive_sdc_3d,reactive_sdc_source_evaluator,compute_reactive_external_source_3d
  use reactive_1d_mod, only: reactive_conserved_to_primitive
  use mixture_thermo_mod, only: mixture_pressure
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options,spray_film,valid_spray_options,valid_spray_liquid
  use spray_mixture_mod, only: valid_component_inventory,advance_mixture_frozen_gas,mixture_trace_fraction, &
    liquid_mixture_properties,build_mixture_film,apply_mixture_breakup_viscosity
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,valid_wall_options, &
    valid_wall_component_inventory,wall_liquid_conductivity_at_temperature,mixture_wall_film_rates, &
    advance_mixture_wall_film_exchange
  use spray_coupling_3d_mod, only: apply_parcel_boundaries,valid_spray_boundary_faces
  use spray_breakup_mod, only: advance_spray_breakup,advance_tab_no_birth_khrt
  use spray_mpi_mod, only: sample_mpi_spray_field,deposit_mpi_spray_field,advance_mpi_khrt_shedding, &
    split_mpi_reference_etab_parcels,migrate_spray_parcels,impact_mpi_planar_walls,migrate_mpi_wall_films, &
    split_mpi_spray_parcels
  use les_mpi_boundary_mod, only: fill_mpi_les_boundary_ghosts_3d,fill_mpi_sdc_source_ghosts_3d
  use les_smagorinsky_3d_mod, only: compute_les_ghosted_rhs_3d
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_rhs_ghosted_3d
  use les_filter_3d_mod, only: les_filter_weights,les_filter_code
  use spray_partition_mod, only: spray_cartesian_block
  implicit none
  private
  type, public :: mpi_sdc_film_stage
    integer :: next_stage=0,iterations=0,cells(3)=0,parts(3)=0
    real(dp) :: time=0,interval=0,lower(3)=0,spacing(3)=0,old_wall(5)=0
    character(len=32) :: faces(6)=''
    logical :: births_attached=.false.
    type(spray_wall_film), allocatable :: midpoint(:)
    real(dp), allocatable :: midpoint_mass(:,:),old_source(:,:,:,:),material_values(:)
    character(len=:), allocatable :: material_names
  end type mpi_sdc_film_stage
  type, public :: mpi_sdc_particle_stage
    integer :: next_stage=0,iterations=0,cells(3)=0,parts(3)=0
    integer(int64) :: initial_watermark=0,midpoint_watermark=0
    real(dp) :: time=0,interval=0,lower(3)=0,spacing(3)=0,acceleration(3)=0,flow_cfl=-1
    real(dp) :: old_wall(5)=0,old_force(5)=0
    character(len=32) :: faces(6)=''
    type(spray_parcel), allocatable :: midpoint(:)
    type(spray_wall_film), allocatable :: midpoint_films(:)
    real(dp), allocatable :: midpoint_mass(:,:),midpoint_film_mass(:,:),old_source(:,:,:,:),old_boundary(:)
    real(dp), allocatable :: material_values(:)
    character(len=:), allocatable :: material_names
  end type mpi_sdc_particle_stage
  public :: advance_mpi_cartesian_gas_sdc,kick_mpi_sdc_films,agree_mpi_sdc_materials
  public :: advance_mpi_sdc_film_stage,append_mpi_sdc_film_births
  public :: sample_mpi_sdc_primitive,deposit_mpi_sdc_cell_source,kick_mpi_sdc_transfer,kick_mpi_sdc_breakup
  public :: drift_mpi_sdc_particles
  public :: advance_mpi_sdc_particle_stage
  public :: advance_mpi_cartesian_spray_sdc
contains
  ! Explicit values avoid comparing derived-type padding. The mechanism and
  ! transport identity is owned by the gas/app contract; particle material
  ! tables and optional wall controls are agreed here before stage branching.
  subroutine agree_mpi_sdc_materials(comm,liquids,options,ok,settings,values_snapshot,names_snapshot,breakup_properties)
    type(MPI_Comm), intent(in) :: comm
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    logical, intent(out) :: ok
    type(spray_wall_options), optional, intent(in) :: settings
    type(spray_options), optional, intent(in) :: breakup_properties(:)
    real(dp), allocatable, optional, intent(inout) :: values_snapshot(:)
    character(len=:), allocatable, optional, intent(inout) :: names_snapshot
    real(dp), allocatable :: values(:),reference(:)
    character(len=:), allocatable :: names,reference_names
    integer :: config(3),minimum(3),maximum(3),p,status
    logical :: valid,agreed
    ok=.false.;config=[size(liquids),merge(1,0,present(settings)),-1]
    if(present(breakup_properties)) config(3)=size(breakup_properties)
    call MPI_Allreduce(config,minimum,3,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(config,maximum,3,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(minimum/=maximum)) return
    valid=size(liquids)>=1.and.size(liquids)<=16.and.valid_spray_options(options)
    do p=1,size(liquids)
      valid=valid.and.valid_spray_liquid(liquids(p))
    enddo
    if(present(breakup_properties)) then
      valid=valid.and.size(breakup_properties)>=1.and.size(breakup_properties)<=16
      do p=1,size(breakup_properties)
        valid=valid.and.valid_spray_options(breakup_properties(p))
      enddo
    endif
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    values=[merge(1._dp,0._dp,options%fixed_parts),merge(1._dp,0._dp,options%khrt_children), &
      options%khrt_reference_diameter_cubed,options%khrt_parcel_factor,options%etab_parcel_factor, &
      options%khrt_b0,options%khrt_b1,options%khrt_c3,options%maximum_multiplicity,real(options%maximum_parcels,dp), &
      options%surface_tension,options%liquid_viscosity,merge(1._dp,0._dp,options%drag), &
      merge(1._dp,0._dp,options%heat),merge(1._dp,0._dp,options%evaporation), &
      merge(1._dp,0._dp,options%abramzon_sirignano),options%max_fractional_change,options%particle_cfl, &
      options%minimum_diameter,real(options%maximum_substeps,dp),options%liquid_viscosity_coefficients, &
      options%liquid_viscosity_temperature_range]
    names=options%breakup_model//options%drag_model//options%mixture_time_scheme//options%liquid_viscosity_model
    do p=1,size(liquids)
      associate(l=>liquids(p))
        values=[values,real(l%vapor_index,dp),l%density,l%heat_capacity,l%reference_temperature,l%reference_pressure, &
          l%latent_heat,l%reference_energy,l%antoine_coefficients,l%antoine_temperature_range, &
          l%density_coefficients,l%density_temperature_range,l%critical_temperature, &
          l%heat_capacity_coefficients,l%heat_capacity_temperature_range,l%normal_boiling_temperature, &
          l%molecular_weight,l%viscosity,l%viscosity_coefficients,l%viscosity_temperature_range, &
          l%conductivity,l%conductivity_coefficients,l%conductivity_temperature_range]
        names=names//l%vapor_name//l%vapor_pressure_model//l%density_model//l%boiling_model//l%heat_capacity_model// &
          l%viscosity_model//l%conductivity_model
      end associate
    enddo
    if(present(settings)) then
      values=[values,merge(1._dp,0._dp,settings%enabled),settings%temperature,settings%contact_angle, &
        real(settings%maximum_films,dp),merge(1._dp,0._dp,settings%evolve),settings%liquid_conductivity, &
        settings%minimum_height,settings%conductivity_coefficients,settings%conductivity_temperature_range, &
        settings%port_conductivities,reshape(settings%port_conductivity_coefficients,[64]), &
        reshape(settings%port_conductivity_temperature_ranges,[32])]
      names=names//settings%conductivity_model
      do p=1,16
        names=names//settings%port_conductivity_models(p)
      enddo
    endif
    if(present(breakup_properties)) then
      do p=1,size(breakup_properties)
        associate(b=>breakup_properties(p))
          values=[values,b%surface_tension,b%liquid_viscosity,b%liquid_viscosity_coefficients, &
            b%liquid_viscosity_temperature_range]
          names=names//b%liquid_viscosity_model
        end associate
      enddo
    endif
    reference=values;reference_names=names
    call MPI_Bcast(reference,size(reference),MPI_DOUBLE_PRECISION,0,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Bcast(reference_names,len(reference_names),MPI_CHARACTER,0,comm,status)
    if(status/=MPI_SUCCESS) return
    valid=all(ieee_is_finite(values)).and.all(values==reference).and.names==reference_names
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    if(present(values_snapshot)) call move_alloc(values,values_snapshot)
    if(present(names_snapshot)) call move_alloc(names,names_snapshot)
    ok=.true.
  end subroutine agree_mpi_sdc_materials

  ! Predictor films stay private. Intermediate gas corrections reuse the old
  ! film rate; only the final correction advances midpoint films and commits.
  subroutine advance_mpi_sdc_film_stage(comm,species,transport,liquids,options,settings,state,temperature, &
      lower,spacing,cells,parts,time,interval,stage,iterations,faces,films,film_mass,work,source,wall_change,ok)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(spray_wall_options), intent(in) :: settings
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),lower(3),spacing(3),time,interval
    integer, intent(in) :: cells(3),parts(3),stage,iterations
    character(len=*), intent(in) :: faces(6)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    real(dp), allocatable, intent(inout) :: film_mass(:,:)
    type(mpi_sdc_film_stage), intent(inout) :: work
    real(dp), intent(inout) :: source(:,:,:,:),wall_change(5)
    logical, intent(out) :: ok
    type(mpi_sdc_film_stage) :: candidate
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: masses(:,:),rhs(:,:,:,:),values(:)
    character(len=:), allocatable :: names
    character(len=32) :: reference_faces(6)
    real(dp) :: controls(8),minimum(8),maximum(8),wall(5),half
    integer :: config(8),low_config(8),high_config(8),status
    logical :: valid,agreed
    ok=.false.;controls=[time,interval,lower,spacing];config=[cells,parts,stage,iterations]
    valid=all(ieee_is_finite(controls)).and.time>=0.and.interval>0.and.interval<=1e6_dp
    valid=valid.and.all(spacing>=1e-12_dp).and.all(spacing<=1e9_dp).and.maxval(abs(lower))<=1e9_dp
    valid=valid.and.stage>=0.and.iterations>=1.and.stage<=iterations
    valid=valid.and.allocated(films).and.allocated(film_mass).and.valid_wall_options(settings,faces)
    valid=valid.and.all(shape(source)==shape(state))
    call agree();if(.not.agreed) return
    call MPI_Allreduce(config,low_config,8,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(config,high_config,8,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low_config/=high_config)) return
    call MPI_Allreduce(controls,minimum,8,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,8,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(minimum/=maximum)) return
    reference_faces=faces
    call MPI_Bcast(reference_faces,6*len(reference_faces),MPI_CHARACTER,0,comm,status)
    if(status/=MPI_SUCCESS) return
    valid=all(faces==reference_faces);call agree();if(.not.agreed) return
    call agree_mpi_sdc_materials(comm,liquids,options,valid,settings,values,names)
    if(.not.valid) return
    if(stage==0) then
      valid=work%next_stage==0
    else
      valid=work%next_stage==stage.and.work%iterations==iterations.and.work%time==time.and.work%interval==interval
      valid=valid.and.all(work%cells==cells).and.all(work%parts==parts).and.all(work%lower==lower)
      valid=valid.and.all(work%spacing==spacing).and.all(work%faces==faces)
      valid=valid.and.allocated(work%material_values).and.allocated(work%material_names)
      valid=valid.and.allocated(work%midpoint).and.allocated(work%midpoint_mass).and.allocated(work%old_source)
      if(valid) valid=size(work%material_values)==size(values).and.all(shape(work%old_source)==shape(source))
      if(valid) valid=all(work%material_values==values).and.work%material_names==names
    endif
    call agree();if(.not.agreed) return
    half=.5_dp*interval
    valid=half>0;call agree();if(.not.agreed) return
    candidate=work
    if(stage>0.and.stage<iterations) then
      source=work%old_source;wall_change=work%old_wall;work%next_stage=stage+1;ok=.true.;return
    endif
    if(stage==0) then
      patches=films;masses=film_mass
    else
      patches=work%midpoint;masses=work%midpoint_mass
    endif
    rhs=source;wall=0
    call kick_mpi_sdc_films(comm,species,transport,liquids,options,settings,state,temperature,lower,spacing, &
      cells,parts,half,faces,patches,masses,rhs,wall,valid)
    if(.not.valid) return
    if(stage==0) then
      candidate%midpoint=patches;candidate%midpoint_mass=masses;candidate%old_source=rhs;candidate%old_wall=wall
      candidate%time=time;candidate%interval=interval;candidate%iterations=iterations
      candidate%cells=cells;candidate%parts=parts;candidate%lower=lower;candidate%spacing=spacing;candidate%faces=faces
      candidate%material_values=values;candidate%material_names=names;candidate%births_attached=.false.
    else
      valid=.not.any(wall>huge(1._dp)-max(work%old_wall,0._dp)).and. &
        .not.any(wall< -huge(1._dp)-min(work%old_wall,0._dp))
      call agree();if(.not.agreed) return
      wall=wall+work%old_wall
    endif
    candidate%next_stage=stage+1
    if(stage==iterations) then
      call move_alloc(patches,films);call move_alloc(masses,film_mass);candidate%next_stage=0
    endif
    work=candidate;source=rhs;wall_change=wall;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine advance_mpi_sdc_film_stage

  ! Wall deposits born during predictor drift enter only the final film kick.
  subroutine append_mpi_sdc_film_births(comm,films,mass,maximum_films,work,ok)
    type(MPI_Comm), intent(in) :: comm
    type(spray_wall_film), intent(in) :: films(:)
    real(dp), intent(in) :: mass(:,:)
    integer, intent(in) :: maximum_films
    type(mpi_sdc_film_stage), intent(inout) :: work
    logical, intent(out) :: ok
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: masses(:,:)
    integer :: np,nold,nc,p,status,low_limit,high_limit
    logical :: valid,agreed
    ok=.false.;np=size(films);nc=size(mass,1)
    valid=work%next_stage==1.and.work%iterations>=1.and..not.work%births_attached
    valid=valid.and.allocated(work%midpoint).and.allocated(work%midpoint_mass)
    valid=valid.and.nc>=1.and.nc<=16.and.size(mass,2)==np.and.maximum_films>=1.and.maximum_films<=1000000
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(maximum_films,low_limit,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(maximum_films,high_limit,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low_limit/=high_limit) return
    nold=size(work%midpoint)
    valid=np<=maximum_films.and.nold<=maximum_films-np.and.all(shape(work%midpoint_mass)==[nc,nold])
    do p=1,np
      valid=valid.and.valid_wall_component_inventory(films(p),mass(:,p))
    enddo
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    allocate(patches(nold+np),masses(nc,nold+np))
    patches(:nold)=work%midpoint;patches(nold+1:)=films
    masses(:,:nold)=work%midpoint_mass;masses(:,nold+1:)=mass
    call move_alloc(patches,work%midpoint);call move_alloc(masses,work%midpoint_mass)
    work%births_attached=.true.;ok=.true.
  end subroutine append_mpi_sdc_film_births

  ! Each owned film samples one immutable primitive gas stage. Local adaptive
  ! substeps contain no MPI; acceptance precedes source routing and publication.
  subroutine kick_mpi_sdc_films(comm,species,transport,liquids,options,settings,state,temperature, &
      lower,spacing,cells,parts,duration,faces,films,film_mass,source,wall_change,ok)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(spray_wall_options), intent(in) :: settings
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),lower(3),spacing(3),duration
    integer, intent(in) :: cells(3),parts(3)
    character(len=*), intent(in) :: faces(6)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    real(dp), allocatable, intent(inout) :: film_mass(:,:)
    real(dp), intent(inout) :: source(:,:,:,:),wall_change(5)
    logical, intent(out) :: ok
    type(spray_wall_film), allocatable :: patches(:),kept_films(:)
    type(spray_wall_film) :: trial
    real(dp), allocatable :: masses(:,:),kept_mass(:,:),rhs(:,:,:,:),positions(:,:),primitive(:,:),temps(:),rates(:,:)
    real(dp) :: evaporation(size(liquids)),hv(size(liquids)),trial_mass(size(liquids)),qg,qw,tdot
    real(dp) :: conductivity,elapsed,dt,change(size(state,1)),wall_energy,local_wall(5),total_wall(5)
    real(dp) :: low_duration,high_duration
    integer :: p,c,np,axis,face,total_steps,attempt,status
    logical :: valid,agreed,accepted
    ok=.false.
    call agree_mpi_sdc_materials(comm,liquids,options,valid,settings)
    if(.not.valid) return
    valid=allocated(films).and.allocated(film_mass).and.valid_wall_options(settings,faces)
    valid=valid.and.ieee_is_finite(duration).and.duration>0.and.duration<=1e6_dp
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(duration,low_duration,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(duration,high_duration,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low_duration/=high_duration) return
    np=size(films)
    valid=all(shape(film_mass)==[size(liquids),np]).and.np<=settings%maximum_films
    valid=valid.and.all(shape(source)==shape(state)).and.size(transport)==size(species)
    valid=valid.and.all(liquids%vapor_index<=size(species))
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    do p=1,np
      valid=valid.and.valid_wall_component_inventory(films(p),film_mass(:,p))
    enddo
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    patches=films;masses=film_mass;allocate(positions(3,np),rates(size(state,1),np));rates=0
    local_wall=0;total_steps=0
    do p=1,np
      positions(:,p)=patches(p)%position
    enddo
    if(settings%evolve) then
      call sample_mpi_sdc_primitive(comm,species,state,temperature,positions,lower,spacing,cells,parts,faces, &
        primitive,temps,valid)
      if(.not.valid) return
      film_loop: do p=1,np
        axis=maxloc(abs(patches(p)%normal),dim=1)
        face=2*axis-merge(1,0,patches(p)%normal(axis)>0)
        if(.not.settings%enabled(face)) then
          valid=.false.;exit
        endif
        elapsed=0
        do while(elapsed<duration.and.patches(p)%mass>0)
          if(total_steps>=options%maximum_substeps) then
            valid=.false.;exit film_loop
          endif
          conductivity=wall_liquid_conductivity_at_temperature(settings,patches(p)%material_id,patches(p)%temperature)
          call mixture_wall_film_rates(patches(p),liquids,masses(:,p),species,transport,primitive(:,p),temps(p),options, &
            conductivity,settings%temperature(face),evaporation,hv,qg,qw,tdot,valid)
          if(.not.valid) exit film_loop
          dt=duration-elapsed
          do c=1,size(liquids)
            if(evaporation(c)>0.and.masses(c,p)>mixture_trace_fraction*patches(p)%mass) &
              dt=min(dt,options%max_fractional_change*masses(c,p)/evaporation(c))
          enddo
          if(abs(tdot)>0) dt=min(dt,min(2._dp,options%max_fractional_change*patches(p)%temperature)/abs(tdot))
          accepted=.false.
          do attempt=1,30
            if(dt<=0.or.elapsed+dt<=elapsed) exit
            trial=patches(p);trial_mass=masses(:,p)
            call advance_mixture_wall_film_exchange(trial,liquids,trial_mass,evaporation*dt,hv,qg*dt,qw*dt, &
              settings%minimum_height,change,wall_energy,valid,deplete_traces=.true.,pressure=primitive(5,p), &
              gas_molecular_weights=species(liquids%vapor_index)%molecular_weight, &
              euler_temperature=options%mixture_time_scheme=='euler')
            if(valid.and.trial%mass>0) valid= &
              wall_liquid_conductivity_at_temperature(settings,trial%material_id,trial%temperature)>0
            if(valid) then
              accepted=.true.;exit
            endif
            dt=.5_dp*dt
          enddo
          if(.not.accepted) then
            valid=.false.;exit film_loop
          endif
          patches(p)=trial;masses(:,p)=trial_mass;local_wall(5)=local_wall(5)-wall_energy
          if(duration<1) then
            if(any(abs(change)>huge(1._dp)*duration)) then
              valid=.false.;exit film_loop
            endif
          endif
          rates(:,p)=rates(:,p)+change/duration
          elapsed=min(duration,elapsed+dt);total_steps=total_steps+1
        enddo
      enddo film_loop
    endif
    valid=valid.and.all(ieee_is_finite(rates)).and.all(ieee_is_finite(local_wall))
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    rhs=source;rhs=0
    call deposit_mpi_sdc_cell_source(comm,rhs,positions,rates,lower,spacing,cells,parts,faces,valid)
    if(.not.valid) return
    call MPI_Allreduce(local_wall,total_wall,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if(status/=MPI_SUCCESS) return
    valid=all(ieee_is_finite(rhs)).and.all(ieee_is_finite(total_wall))
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    allocate(kept_mass(size(liquids),count(patches%mass>0)))
    do c=1,size(liquids)
      kept_mass(c,:)=pack(masses(c,:),patches%mass>0)
    enddo
    kept_films=pack(patches,patches%mass>0)
    call move_alloc(kept_films,films);call move_alloc(kept_mass,film_mass)
    source=rhs;wall_change=total_wall;ok=.true.
  end subroutine kick_mpi_sdc_films

  subroutine sample_mpi_sdc_primitive(comm,species,state,temperature,positions,lower,spacing,cells,parts,faces, &
      primitive,gas_temperature,ok)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),positions(:,:),lower(3),spacing(3)
    integer, intent(in) :: cells(3),parts(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), allocatable, intent(inout) :: primitive(:,:),gas_temperature(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: field(:,:,:,:),samples(:,:),candidate(:,:),temps(:)
    real(dp) :: recovered,sound
    integer :: n(3),nv,i,j,k,p,status
    logical :: valid,agreed
    ok=.false.;n=shape(temperature);nv=5+size(species)
    valid=all(shape(state)==[nv,n]).and.size(positions,1)==3
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    allocate(field(nv+1,n(1),n(2),n(3)));valid=.true.
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k),field(:nv,i,j,k), &
        recovered,sound,valid)
      if(.not.valid) exit
      field(nv+1,i,j,k)=recovered
    enddo
    if(.not.valid) exit
    enddo
    if(.not.valid) exit
    enddo
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call sample_mpi_spray_field(comm,field,positions,lower,spacing,cells,parts,faces,samples,valid)
    if(.not.valid) return
    candidate=samples(:nv,:);temps=samples(nv+1,:)
    do p=1,size(temps)
      candidate(5,p)=mixture_pressure(species,candidate(6:,p),candidate(1,p),temps(p),valid)
      if(.not.valid) exit
    enddo
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call move_alloc(candidate,primitive);call move_alloc(temps,gas_temperature);ok=.true.
  end subroutine

  ! Source feedback belongs to the containing cell, unlike interpolation.
  subroutine deposit_mpi_sdc_cell_source(comm,source,positions,changes,lower,spacing,cells,parts,faces,ok)
    type(MPI_Comm), intent(in) :: comm
    real(dp), intent(inout) :: source(:,:,:,:)
    real(dp), intent(in) :: positions(:,:),changes(:,:),lower(3),spacing(3)
    integer, intent(in) :: cells(3),parts(3)
    character(len=*), intent(in) :: faces(6)
    logical, intent(out) :: ok
    call deposit_mpi_spray_field(comm,source,positions,changes,lower,spacing,cells,parts,faces,ok, &
      containing_cell=.true.)
  end subroutine

  ! Transfer-only half kick. The SDC stage owner performs breakup, drift,
  ! wall impact and population migration.
  subroutine kick_mpi_sdc_transfer(comm,species,transport,liquids,options,state,temperature,parcels,masses, &
      lower,spacing,cells,parts,faces,interval,acceleration,source,force_change,ok, &
      sampled_primitive,sampled_temperature,reynolds_number,breakup_properties)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),lower(3),spacing(3),interval,acceleration(3)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: masses(:,:)
    integer, intent(in) :: cells(3),parts(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(inout) :: source(:,:,:,:),force_change(5)
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: sampled_primitive(:,:),sampled_temperature(:),reynolds_number(:)
    type(spray_options), optional, intent(in) :: breakup_properties(:)
    type(spray_options) :: properties
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: component_mass(:,:),positions(:,:),primitive(:,:),temps(:),changes(:,:),rhs(:,:,:,:),reynolds(:)
    real(dp) :: controls(4),low(4),high(4),local_force(5),global_force(5),parcel_force(5),half
    integer :: p,np,status
    logical :: valid,agreed
    ok=.false.;controls=[interval,acceleration]
    call agree_mpi_sdc_materials(comm,liquids,options,valid,breakup_properties=breakup_properties)
    if(.not.valid) return
    valid=allocated(parcels).and.allocated(masses).and.valid_spray_options(options)
    valid=valid.and.all(ieee_is_finite(controls)).and.interval>0.and.all(abs(acceleration)<=1e6_dp)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(controls,low,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    np=size(parcels);half=.5_dp*interval
    valid=half>0.and.all(shape(masses)==[size(liquids),np]).and.all(shape(source)==shape(state))
    valid=valid.and.all(liquids%vapor_index<=size(species)).and.size(transport)==size(species)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    drops=parcels;component_mass=masses
    allocate(positions(3,np),changes(size(state,1),np),reynolds(np));local_force=0;reynolds=0
    do p=1,np
      positions(:,p)=drops(p)%position
      valid=valid.and.valid_component_inventory(drops(p),component_mass(:,p))
      if(present(breakup_properties)) valid=valid.and.drops(p)%material_id>=1.and. &
        drops(p)%material_id<=size(breakup_properties)
    enddo
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call sample_mpi_sdc_primitive(comm,species,state,temperature,positions,lower,spacing,cells,parts,faces, &
      primitive,temps,valid)
    if(.not.valid) return
    changes=0
    do p=1,np
      if(drops(p)%mass<=0) cycle
      properties=options
      if(present(breakup_properties)) call select_mpi_sdc_properties(breakup_properties(drops(p)%material_id),properties)
      if(options%mixture_time_scheme=='euler') then
        call advance_mixture_frozen_gas(species,transport,primitive(:,p),temps(p),drops(p),liquids,component_mass(:,p), &
          properties,interval,changes(:,p),valid,time_modifier=.5_dp,acceleration=acceleration, &
          force_change=parcel_force,reynolds_number=reynolds(p))
      else
        call advance_mixture_frozen_gas(species,transport,primitive(:,p),temps(p),drops(p),liquids,component_mass(:,p), &
          properties,half,changes(:,p),valid,acceleration=acceleration,force_change=parcel_force,reynolds_number=reynolds(p))
      endif
      if(.not.valid) exit
      drops(p)%position=positions(:,p)
      local_force=local_force+parcel_force
    enddo
    valid=valid.and.all(ieee_is_finite(local_force)).and.all(ieee_is_finite(changes))
    if(half<1) valid=valid.and.all(abs(changes)<=huge(1._dp)*half)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    changes=changes/half
    rhs=source;rhs=0
    call deposit_mpi_sdc_cell_source(comm,rhs,positions,changes,lower,spacing,cells,parts,faces,valid)
    if(.not.valid) return
    call MPI_Allreduce(local_force,global_force,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if(status/=MPI_SUCCESS) return
    valid=all(ieee_is_finite(rhs)).and.all(ieee_is_finite(global_force))
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call move_alloc(drops,parcels);call move_alloc(component_mass,masses)
    if(present(sampled_primitive)) call move_alloc(primitive,sampled_primitive)
    if(present(sampled_temperature)) call move_alloc(temps,sampled_temperature)
    if(present(reynolds_number)) call move_alloc(reynolds,reynolds_number)
    source=rhs;force_change=global_force;ok=.true.
  end subroutine kick_mpi_sdc_transfer

  subroutine select_mpi_sdc_properties(material,properties)
    type(spray_options), intent(in) :: material
    type(spray_options), intent(inout) :: properties
    properties%surface_tension=material%surface_tension
    properties%liquid_viscosity=material%liquid_viscosity
    properties%liquid_viscosity_model=material%liquid_viscosity_model
    properties%liquid_viscosity_coefficients=material%liquid_viscosity_coefficients
    properties%liquid_viscosity_temperature_range=material%liquid_viscosity_temperature_range
  end subroutine select_mpi_sdc_properties

  ! Consume the pre-transfer gas sample and transfer Reynolds number, retaining
  ! MKD/MK clocks. Only the old-time MKD caller requests replacement births.
  ! This routine owns breakup and its gas compensation; drift is separate.
  subroutine kick_mpi_sdc_breakup(comm,species,transport,liquids,options,parcels,masses,primitive,gas_temperature, &
      reynolds,lower,spacing,cells,parts,faces,flow_duration,khrt_duration,source_half,make_new_drops, &
      do_khrt,last_id,source,ok,breakup_properties,tangent_accumulator)
    use, intrinsic :: iso_fortran_env, only: int64
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: masses(:,:)
    real(dp), intent(in) :: primitive(:,:),gas_temperature(:),reynolds(:),lower(3),spacing(3)
    integer, intent(in) :: cells(3),parts(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: flow_duration,khrt_duration,source_half
    logical, intent(in) :: make_new_drops,do_khrt
    integer(int64), intent(inout) :: last_id
    real(dp), intent(inout) :: source(:,:,:,:)
    logical, intent(out) :: ok
    type(spray_options), optional, intent(in) :: breakup_properties(:)
    ! Indexed by the pre-birth local parent population. The stage owner moves
    ! this payload with parents between flow subcycles and resets it after birth.
    real(dp), allocatable, optional, intent(inout) :: tangent_accumulator(:)
    type(spray_parcel), allocatable :: drops(:)
    type(spray_options) :: properties
    type(spray_liquid) :: effective
    type(spray_film) :: film
    real(dp), allocatable :: mass(:,:),positions(:,:),tangent(:),changes(:,:),rhs(:,:,:,:),feedback(:,:)
    real(dp), allocatable :: birth_positions(:,:),energy(:),momentum(:,:)
    real(dp) :: fraction(size(liquids)),transfer(size(liquids)),enthalpy(size(liquids))
    real(dp) :: density,cp,e,old_mass,clock,tangent_delta,controls(6),low(6),high(6)
    integer(int64) :: watermark,low_watermark,high_watermark
    integer :: p,np,nv,status
    logical :: valid,agreed
    ok=.false.
    call agree_mpi_sdc_materials(comm,liquids,options,valid,breakup_properties=breakup_properties)
    if(.not.valid) return
    controls=[flow_duration,khrt_duration,source_half,merge(1._dp,0._dp,make_new_drops),merge(1._dp,0._dp,do_khrt), &
      merge(1._dp,0._dp,present(tangent_accumulator))]
    valid=allocated(parcels).and.allocated(masses).and.all(ieee_is_finite(controls)).and.last_id>=0
    valid=valid.and.flow_duration>0.and.khrt_duration>0.and.source_half>0
    valid=valid.and.maxval(controls(:3))<=1e6_dp.and.all(liquids%vapor_index<=size(species))
    call agree();if(.not.agreed) return
    call MPI_Allreduce(controls,low,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    call MPI_Allreduce(last_id,low_watermark,1,MPI_INTEGER8,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(last_id,high_watermark,1,MPI_INTEGER8,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low_watermark/=high_watermark) return
    np=size(parcels);nv=5+size(species)
    valid=all(shape(masses)==[size(liquids),np]).and.all(shape(primitive)==[nv,np])
    valid=valid.and.size(gas_temperature)==np.and.size(reynolds)==np.and.size(transport)==size(species)
    valid=valid.and.size(source,1)==nv.and.np<=options%maximum_parcels
    if(present(tangent_accumulator)) then
      valid=valid.and.allocated(tangent_accumulator)
      if(valid) valid=size(tangent_accumulator)==np
      if(valid) valid=all(ieee_is_finite(tangent_accumulator)).and.all(tangent_accumulator>=0)
    endif
    call agree();if(.not.agreed) return
    valid=all(ieee_is_finite(primitive)).and.all(ieee_is_finite(gas_temperature)).and.all(ieee_is_finite(reynolds))
    valid=valid.and.all(reynolds>=0).and.all(ieee_is_finite(source))
    do p=1,np
      valid=valid.and.valid_component_inventory(parcels(p),masses(:,p)).and.parcels(p)%id<=last_id
      if(p>1) valid=valid.and.parcels(p)%id>parcels(p-1)%id
      if(present(breakup_properties)) valid=valid.and.parcels(p)%material_id>=1.and. &
        parcels(p)%material_id<=size(breakup_properties)
    enddo
    call agree();if(.not.agreed) return
    drops=parcels;mass=masses;watermark=last_id;rhs=source
    allocate(positions(3,np),tangent(np),changes(nv,np));tangent=0;changes=0
    if(present(tangent_accumulator)) tangent=tangent_accumulator
    do p=1,np
      positions(:,p)=drops(p)%position
      if(drops(p)%mass<=0.or.options%breakup_model=='none') cycle
      if(options%breakup_model=='khrt'.and.(.not.do_khrt.or.options%khrt_children)) cycle
      properties=options
      if(present(breakup_properties)) call select_mpi_sdc_properties(breakup_properties(drops(p)%material_id),properties)
      old_mass=drops(p)%mass;fraction=mass(:,p)/old_mass
      call liquid_mixture_properties(liquids,fraction,drops(p)%temperature,density,cp,e,valid, &
        primitive(5,p),species(liquids%vapor_index)%molecular_weight)
      if(.not.valid) exit
      effective=liquids(1);effective%density=density;effective%density_model='constant'
      effective%density_coefficients=0;effective%density_temperature_range=0
      effective%heat_capacity=cp;effective%heat_capacity_model='constant'
      effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
      effective%reference_temperature=drops(p)%temperature;effective%reference_energy=e
      effective%boiling_model='reject';effective%critical_temperature=0;effective%normal_boiling_temperature=0
      call build_mixture_film(species,transport,primitive(:,p),gas_temperature(p),drops(p),liquids,fraction,properties, &
        film,transfer,enthalpy,valid)
      if(.not.valid) exit
      call apply_mixture_breakup_viscosity(liquids,fraction,drops(p)%temperature,properties,valid, &
        primitive(5,p),species(liquids%vapor_index)%molecular_weight)
      if(.not.valid) exit
      clock=khrt_duration
      if(options%breakup_model=='etab') clock=flow_duration
      if(options%etab_parcel_factor>=0) then
        call advance_spray_breakup(drops(p),effective,properties,primitive(1,p),primitive(2:4,p), &
          film%density,film%viscosity,clock,valid,tangent_delta,reynolds_number=reynolds(p))
        if(valid) then
          valid=ieee_is_finite(tangent_delta).and.tangent_delta>=0
          if(valid) valid=tangent(p)<=huge(1._dp)-tangent_delta
          if(valid) tangent(p)=tangent(p)+tangent_delta
        endif
        if(valid.and..not.make_new_drops) call advance_tab_no_birth_khrt(drops(p),effective,properties, &
          primitive(1,p),primitive(2:4,p),reynolds(p),khrt_duration,valid)
      else
        call advance_spray_breakup(drops(p),effective,properties,primitive(1,p),primitive(2:4,p), &
          film%density,film%viscosity,clock,valid,reynolds_number=reynolds(p))
      endif
      if(.not.valid) exit
      mass(:,p)=mass(:,p)*(drops(p)%mass/old_mass)
    enddo
    call agree();if(.not.agreed) return
    if(options%khrt_children.and.do_khrt) then
      call advance_mpi_khrt_shedding(comm,liquids,options,drops,mass,primitive(1,:),primitive(2:4,:),reynolds, &
        khrt_duration,watermark,energy,valid,gas_momentum=momentum,pressure=primitive(5,:), &
        gas_molecular_weights=species(liquids%vapor_index)%molecular_weight,birth_positions=positions, &
        child_drift_dt=0._dp,make_new_drops=make_new_drops,breakup_properties=breakup_properties)
      if(.not.valid) return
      changes(2:4,:)=momentum;changes(5,:)=energy
    endif
    if(options%etab_parcel_factor>=0.and.make_new_drops.and.do_khrt) then
      call split_mpi_reference_etab_parcels(comm,drops,mass,options%maximum_multiplicity,options%maximum_parcels, &
        options%etab_parcel_factor,tangent,0._dp,watermark,birth_positions,feedback,valid)
      if(.not.valid) return
      positions=birth_positions;changes(2:5,:)=feedback
    endif
    valid=all(ieee_is_finite(changes))
    if(source_half<1) valid=valid.and.all(abs(changes)<=huge(1._dp)*source_half)
    call agree();if(.not.agreed) return
    changes=changes/source_half
    call deposit_mpi_sdc_cell_source(comm,rhs,positions,changes,lower,spacing,cells,parts,faces,valid)
    if(.not.valid) return
    if(present(tangent_accumulator)) call move_alloc(tangent,tangent_accumulator)
    call move_alloc(drops,parcels);call move_alloc(mass,masses);last_id=watermark;source=rhs;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine kick_mpi_sdc_breakup

  ! The outer stage owns this private flight result. Particle, material mass,
  ! tangent history, films, IDs and boundary ledgers publish only together.
  ! MKD predicts a private migrated midpoint. MK reuses the old rate until
  ! the final stage, which commits the second kick and population split.
  ! Film thermal stages remain a separate source family; these films hold
  ! wall deposits and migrate with the particle midpoint.
  subroutine advance_mpi_sdc_particle_stage(comm,species,transport,liquids,options,state,temperature, &
      lower,spacing,cells,parts,time,interval,stage,iterations,faces,parcels,masses,work, &
      source,wall_change,boundary_change,force_change,last_id,ok,acceleration,breakup_properties, &
      wall_options,wall_films,film_component_mass,native_flow_cfl)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),lower(3),spacing(3),time,interval
    integer, intent(in) :: cells(3),parts(3),stage,iterations
    character(len=*), intent(in) :: faces(6)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: masses(:,:)
    type(mpi_sdc_particle_stage), intent(inout) :: work
    real(dp), intent(inout) :: source(:,:,:,:),wall_change(5),boundary_change(:),force_change(5)
    integer(int64), intent(inout) :: last_id
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: acceleration(3),native_flow_cfl
    type(spray_options), optional, intent(in) :: breakup_properties(:)
    type(spray_wall_options), optional, intent(in) :: wall_options
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_component_mass(:,:)
    type(mpi_sdc_particle_stage) :: candidate
    type(spray_parcel), allocatable :: drops(:)
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: mass(:,:),patch_mass(:,:),rhs(:,:,:,:),cycle_rhs(:,:,:,:),tangent(:)
    real(dp), allocatable :: primitive(:,:),temps(:),reynolds(:),values(:),kept_mass(:,:)
    character(len=:), allocatable :: names
    character(len=32) :: reference_faces(6)
    logical, allocatable :: alive(:)
    real(dp) :: accel(3),coupled(3),cfl,controls(13),low(13),high(13),half,sub_interval
    real(dp) :: wall(5),outside(size(boundary_change)),forced(5),kick_force(5),drift_wall(5),drift_out(size(boundary_change))
    integer :: config(12),minimum(12),maximum(12),cycles,cycle,c,p,status
    integer(int64) :: watermark,low_id,high_id
    logical :: valid,agreed
    ok=.false.;accel=0;cfl=-1
    if(present(acceleration)) accel=acceleration
    if(present(native_flow_cfl)) cfl=native_flow_cfl
    config=[cells,parts,stage,iterations,size(species),merge(1,0,present(wall_options)), &
      merge(1,0,present(wall_films)),merge(1,0,present(film_component_mass))]
    controls=[lower,spacing,time,interval,accel,cfl,real(work%next_stage,dp)]
    valid=allocated(parcels).and.allocated(masses).and.iterations>=1.and.stage>=0.and.stage<=iterations
    valid=valid.and.all(ieee_is_finite(controls)).and.time>=0.and.interval>0.and.interval<=1e6_dp.and.cfl>=-1
    valid=valid.and.all(abs(accel)<=1e6_dp).and.last_id>=0.and.valid_spray_boundary_faces(faces)
    valid=valid.and.(present(wall_options).eqv.present(wall_films))
    valid=valid.and.(present(wall_films).eqv.present(film_component_mass))
    call agree();if(.not.agreed) return
    call MPI_Allreduce(config,minimum,12,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(config,maximum,12,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(minimum/=maximum)) return
    call MPI_Allreduce(controls,low,13,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,13,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    call MPI_Allreduce(last_id,low_id,1,MPI_INTEGER8,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(last_id,high_id,1,MPI_INTEGER8,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low_id/=high_id) return
    reference_faces=faces
    call MPI_Bcast(reference_faces,6*len(reference_faces),MPI_CHARACTER,0,comm,status)
    if(status/=MPI_SUCCESS) return
    valid=all(reference_faces==faces);call agree();if(.not.agreed) return
    call agree_mpi_sdc_materials(comm,liquids,options,valid,wall_options,values,names,breakup_properties)
    if(.not.valid) return
    valid=all(shape(masses)==[size(liquids),size(parcels)]).and.all(shape(source)==shape(state))
    valid=valid.and.size(boundary_change)==5+size(species).and.size(transport)==size(species)
    valid=valid.and.all(liquids%vapor_index<=size(species)).and.size(parcels)<=options%maximum_parcels
    if(present(wall_films)) then
      valid=valid.and.allocated(wall_films).and.allocated(film_component_mass)
      if(valid) valid=all(shape(film_component_mass)==[size(liquids),size(wall_films)])
      if(valid) valid=valid_wall_options(wall_options,faces).and.size(wall_films)<=wall_options%maximum_films
    endif
    if(stage==0) then
      valid=valid.and.work%next_stage==0
    else
      valid=valid.and.work%next_stage==stage.and.work%iterations==iterations
      valid=valid.and.work%time==time.and.work%interval==interval.and.work%initial_watermark==last_id
      valid=valid.and.all(work%lower==lower).and.all(work%spacing==spacing).and.all(work%cells==cells)
      valid=valid.and.all(work%parts==parts).and.all(work%faces==faces)
      valid=valid.and.all(work%acceleration==accel).and.work%flow_cfl==cfl
      valid=valid.and.allocated(work%midpoint).and.allocated(work%midpoint_mass).and.allocated(work%old_source)
      valid=valid.and.allocated(work%old_boundary)
      valid=valid.and.allocated(work%material_values).and.allocated(work%material_names)
      valid=valid.and.(allocated(work%midpoint_films).eqv.present(wall_films))
      if(valid) valid=size(work%material_values)==size(values).and.all(shape(work%old_source)==shape(source))
      if(valid) valid=all(shape(work%midpoint_mass)==[size(liquids),size(work%midpoint)])
      if(valid) valid=size(work%old_boundary)==size(boundary_change)
      if(valid.and.present(wall_films)) then
        valid=allocated(work%midpoint_film_mass)
        if(valid) valid=all(shape(work%midpoint_film_mass)==[size(liquids),size(work%midpoint_films)])
      endif
      if(valid) valid=all(work%material_values==values).and.work%material_names==names
    endif
    call agree();if(.not.agreed) return
    if(stage>0.and.stage<iterations) then
      source=work%old_source;wall_change=work%old_wall;boundary_change=work%old_boundary
      force_change=work%old_force;work%next_stage=stage+1;ok=.true.;return
    endif
    half=.5_dp*interval;cycles=1
    if(stage==0.and.cfl>.5_dp) then
      valid=cfl<=.5_dp*real(min(options%maximum_substeps,huge(cycles)-1),dp)
      call agree();if(.not.agreed) return
      cycles=ceiling(cfl/.5_dp)
    endif
    sub_interval=interval/real(cycles,dp)
    candidate=work;watermark=last_id
    if(stage==0) then
      drops=parcels;mass=masses
      if(present(wall_films)) then
        patches=wall_films;patch_mass=film_component_mass
      endif
    else
      drops=work%midpoint;mass=work%midpoint_mass;watermark=work%midpoint_watermark
      if(present(wall_films)) then
        patches=work%midpoint_films;patch_mass=work%midpoint_film_mass
      endif
    endif
    rhs=source;rhs=0;cycle_rhs=rhs;wall=0;outside=0;forced=0;coupled=0
    tangent=spread(0._dp,1,size(drops))
    if(options%mixture_time_scheme=='euler'.or.options%fixed_parts) coupled=accel
    if(stage==0.and.options%mixture_time_scheme/='euler'.and..not.options%fixed_parts) then
      call force_half();if(.not.valid) return
    endif
    do cycle=1,cycles
      call kick_mpi_sdc_transfer(comm,species,transport,liquids,options,state,temperature,drops,mass, &
        lower,spacing,cells,parts,faces,sub_interval,coupled,cycle_rhs,kick_force,valid, &
        primitive,temps,reynolds,breakup_properties)
      if(.not.valid) return
      rhs=rhs+cycle_rhs/real(cycles,dp);forced=forced+kick_force
      call kick_mpi_sdc_breakup(comm,species,transport,liquids,options,drops,mass,primitive,temps,reynolds, &
        lower,spacing,cells,parts,faces,sub_interval,half,half,stage==0.and.cycle==cycles,cycle==cycles, &
        watermark,rhs,valid,breakup_properties,tangent)
      if(.not.valid) return
      if(stage==0) then
        if(cycle==cycles) tangent=spread(0._dp,1,size(drops))
        if(present(wall_films)) then
          call drift_mpi_sdc_particles(comm,species,liquids,options,state,temperature,drops,mass, &
            lower,spacing,cells,parts,faces,sub_interval,watermark,drift_wall,drift_out,valid, &
            tangent,wall_options,patches,patch_mass,breakup_properties)
        else
          call drift_mpi_sdc_particles(comm,species,liquids,options,state,temperature,drops,mass, &
            lower,spacing,cells,parts,faces,sub_interval,watermark,drift_wall,drift_out,valid,tangent)
        endif
        if(.not.valid) return
        wall=wall+drift_wall;outside=outside+drift_out
      endif
    enddo
    if(stage==iterations.and.options%mixture_time_scheme/='euler'.and..not.options%fixed_parts) then
      call force_half();if(.not.valid) return
    endif
    if(stage==0) then
      candidate%midpoint=drops;candidate%midpoint_mass=mass;candidate%midpoint_watermark=watermark
      if(present(wall_films)) then
        candidate%midpoint_films=patches;candidate%midpoint_film_mass=patch_mass
      else
        if(allocated(candidate%midpoint_films)) deallocate(candidate%midpoint_films)
        if(allocated(candidate%midpoint_film_mass)) deallocate(candidate%midpoint_film_mass)
      endif
      candidate%old_source=rhs;candidate%old_wall=wall;candidate%old_boundary=outside;candidate%old_force=forced
      candidate%time=time;candidate%interval=interval;candidate%iterations=iterations
      candidate%initial_watermark=last_id;candidate%lower=lower;candidate%spacing=spacing
      candidate%cells=cells;candidate%parts=parts;candidate%faces=faces
      candidate%acceleration=accel;candidate%flow_cfl=cfl
      candidate%material_values=values;candidate%material_names=names
    else
      wall=wall+work%old_wall;outside=outside+work%old_boundary;forced=forced+work%old_force
      alive=drops%mass>0;allocate(kept_mass(size(liquids),count(alive)))
      do c=1,size(liquids)
        kept_mass(c,:)=pack(mass(c,:),alive)
      enddo
      drops=pack(drops,alive);call move_alloc(kept_mass,mass)
      if(options%maximum_multiplicity>0.and.options%etab_parcel_factor<0) then
        call split_mpi_spray_parcels(comm,drops,options%maximum_multiplicity,options%maximum_parcels, &
          watermark,valid,mass)
        if(.not.valid) return
      endif
    endif
    valid=all(ieee_is_finite(rhs)).and.all(ieee_is_finite([wall,outside,forced]))
    call agree();if(.not.agreed) return
    candidate%next_stage=stage+1
    if(stage==iterations) then
      call move_alloc(drops,parcels);call move_alloc(mass,masses);last_id=watermark;candidate%next_stage=0
      if(present(wall_films)) then
        call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_component_mass)
      endif
    endif
    work=candidate;source=rhs;wall_change=wall;boundary_change=outside;force_change=forced;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine
    subroutine force_half()
      real(dp) :: old(5),new(5),local_force(5),global_force(5),rho,cp,energy
      valid=.true.;local_force=0
      do p=1,size(drops)
        if(drops(p)%mass<=0) cycle
        valid=valid_component_inventory(drops(p),mass(:,p));if(.not.valid) exit
        call liquid_mixture_properties(liquids,mass(:,p)/drops(p)%mass,drops(p)%temperature,rho,cp,energy,valid)
        if(.not.valid) exit
        old=drops(p)%mass*drops(p)%multiplicity*[1._dp,drops(p)%velocity,energy+.5_dp*sum(drops(p)%velocity**2)]
        drops(p)%velocity=drops(p)%velocity+half*accel
        valid=all(ieee_is_finite(drops(p)%velocity)).and.all(abs(drops(p)%velocity)<=1e7_dp)
        if(.not.valid) exit
        new=drops(p)%mass*drops(p)%multiplicity*[1._dp,drops(p)%velocity,energy+.5_dp*sum(drops(p)%velocity**2)]
        local_force=local_force+new-old
      enddo
      valid=valid.and.all(ieee_is_finite(local_force));call agree();if(.not.agreed) return
      call MPI_Allreduce(local_force,global_force,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
      valid=status==MPI_SUCCESS
      if(valid) valid=all(ieee_is_finite(global_force))
      call agree();if(.not.agreed) return
      forced=forced+global_force
    end subroutine
  end subroutine advance_mpi_sdc_particle_stage

  subroutine drift_mpi_sdc_particles(comm,species,liquids,options,state,temperature,parcels,masses, &
      lower,spacing,cells,parts,faces,duration,last_id,wall_change,boundary_change,ok, &
      tab_tangent_history,wall_options,wall_films,film_component_mass,breakup_properties)
    use, intrinsic :: iso_fortran_env, only: int64
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),lower(3),spacing(3),duration
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: masses(:,:)
    integer, intent(in) :: cells(3),parts(3)
    character(len=*), intent(in) :: faces(6)
    integer(int64), intent(inout) :: last_id
    real(dp), intent(inout) :: wall_change(5),boundary_change(:)
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: tab_tangent_history(:)
    type(spray_wall_options), optional, intent(in) :: wall_options
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_component_mass(:,:)
    type(spray_options), optional, intent(in) :: breakup_properties(:)
    type(spray_parcel), allocatable :: drops(:)
    type(spray_wall_film), allocatable :: patches(:)
    type(spray_liquid) :: effective
    real(dp), allocatable :: mass(:,:),history(:),patch_mass(:,:),positions(:,:),primitive(:,:),temps(:),before(:,:)
    real(dp), allocatable :: kept_mass(:,:)
    logical, allocatable :: hit(:),alive(:)
    character(len=32) :: reference_faces(6)
    real(dp) :: controls(7),low(7),high(7),upper(3),local_wall(5),global_wall(5),impact_wall(5),total(5),delta(5)
    real(dp) :: exterior(size(boundary_change)),global_exterior(size(boundary_change)),rho,cp,e
    integer :: config(12),minimum(12),maximum(12),p,c,np,status
    integer(int64) :: watermark,low_watermark,high_watermark
    logical :: valid,agreed
    ok=.false.;controls=[lower,spacing,duration]
    config=[cells,parts,size(species),size(liquids),merge(1,0,present(tab_tangent_history)), &
      merge(1,0,present(wall_options)),merge(1,0,present(wall_films)),merge(1,0,present(film_component_mass))]
    valid=allocated(parcels).and.allocated(masses).and.all(ieee_is_finite(controls)).and.duration>0.and.duration<=1e6_dp
    valid=valid.and.size(boundary_change)==5+size(species).and.last_id>=0.and.valid_spray_boundary_faces(faces)
    valid=valid.and.all(spacing>=1e-12_dp).and.all(spacing<=1e9_dp).and.maxval(abs(lower))<=1e9_dp
    valid=valid.and.all(cells>=2).and.all(cells<=100000000).and.all(liquids%vapor_index<=size(species))
    valid=valid.and.(present(wall_options).eqv.present(wall_films))
    valid=valid.and.(present(wall_films).eqv.present(film_component_mass))
    call agree();if(.not.agreed) return
    call MPI_Allreduce(config,minimum,12,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(config,maximum,12,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(minimum/=maximum)) return
    call MPI_Allreduce(controls,low,7,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,7,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    call MPI_Allreduce(last_id,low_watermark,1,MPI_INTEGER8,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(last_id,high_watermark,1,MPI_INTEGER8,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low_watermark/=high_watermark) return
    reference_faces=faces
    call MPI_Bcast(reference_faces,6*len(reference_faces),MPI_CHARACTER,0,comm,status)
    if(status/=MPI_SUCCESS) return
    valid=all(faces==reference_faces);call agree();if(.not.agreed) return
    call agree_mpi_sdc_materials(comm,liquids,options,valid,wall_options,breakup_properties=breakup_properties)
    if(.not.valid) return
    np=size(parcels);upper=lower+spacing*real(cells,dp)
    valid=all(shape(masses)==[size(liquids),np]).and.np<=options%maximum_parcels
    valid=valid.and.all(ieee_is_finite(upper)).and.all(spacing>0).and.all(upper>lower)
    if(present(tab_tangent_history)) then
      valid=valid.and.allocated(tab_tangent_history)
      if(valid) valid=size(tab_tangent_history)==np
      if(valid) valid=all(ieee_is_finite(tab_tangent_history)).and.all(tab_tangent_history>=0)
    endif
    if(present(wall_options)) then
      valid=valid.and.valid_wall_options(wall_options,faces)
      valid=valid.and.allocated(wall_films).and.allocated(film_component_mass)
      if(valid) valid=all(shape(film_component_mass)==[size(liquids),size(wall_films)])
    endif
    call agree();if(.not.agreed) return
    drops=parcels;mass=masses;watermark=last_id
    if(present(tab_tangent_history)) history=tab_tangent_history
    if(present(wall_films)) then
      patches=wall_films;patch_mass=film_component_mass
    endif
    allocate(positions(3,np),before(5,np),hit(np));hit=.false.;local_wall=0;impact_wall=0;exterior=0
    do p=1,np
      valid=valid_component_inventory(drops(p),mass(:,p)).and.drops(p)%id<=watermark
      if(.not.valid) exit
      do c=1,3
        if(faces(2*c)=='periodic') &
          drops(p)%position(c)=lower(c)+modulo(drops(p)%position(c)-lower(c),upper(c)-lower(c))
      enddo
      positions(:,p)=drops(p)%position;before(:,p)=0
      if(drops(p)%mass<=0) cycle
      call inventory(drops(p),mass(:,p),before(:,p),valid)
      if(.not.valid) exit
    enddo
    call agree();if(.not.agreed) return
    if(.not.options%fixed_parts) then
      if(present(wall_options)) then
        if(any(wall_options%enabled)) then
          call sample_mpi_sdc_primitive(comm,species,state,temperature,positions,lower,spacing,cells,parts,faces, &
            primitive,temps,valid)
          if(.not.valid) return
        endif
      endif
      do p=1,np
        if(drops(p)%mass>0) drops(p)%position=drops(p)%position+duration*drops(p)%velocity
      enddo
      valid=.true.
      do p=1,np
        valid=valid.and.all(ieee_is_finite(drops(p)%position)).and.all(abs(drops(p)%position)<=1e9_dp)
      enddo
      call agree();if(.not.agreed) return
      if(present(wall_options)) then
        if(any(wall_options%enabled)) then
          call impact_mpi_planar_walls(comm,liquids,options,wall_options,faces, &
            species(liquids%vapor_index)%molecular_weight,primitive(5,:),positions,lower,upper,drops,mass, &
            patches,patch_mass,watermark,hit,impact_wall,valid,breakup_properties)
          if(.not.valid) return
          if(present(tab_tangent_history)) history=[history,spread(0._dp,1,size(drops)-np)]
        endif
      endif
      do p=1,np
        if(hit(p).or.drops(p)%mass<=0) cycle
        call liquid_mixture_properties(liquids,mass(:,p)/drops(p)%mass,drops(p)%temperature,rho,cp,e,valid)
        if(.not.valid) exit
        effective=liquids(1);effective%reference_temperature=drops(p)%temperature;effective%reference_energy=e
        call apply_parcel_boundaries(drops(p),effective,lower,upper,faces,delta,valid)
        if(.not.valid) exit
        if(drops(p)%mass==0) then
          exterior(:5)=exterior(:5)-before(:,p)
          do c=1,size(liquids)
            exterior(5+liquids(c)%vapor_index)=exterior(5+liquids(c)%vapor_index)-mass(c,p)*drops(p)%multiplicity
          enddo
          mass(:,p)=0
        else
          call inventory(drops(p),mass(:,p),total,valid)
          if(.not.valid) exit
          local_wall=local_wall+total-before(:,p)
        endif
      enddo
      call agree();if(.not.agreed) return
    endif
    alive=drops%mass>0;allocate(kept_mass(size(liquids),count(alive)))
    do c=1,size(liquids)
      kept_mass(c,:)=pack(mass(c,:),alive)
    enddo
    drops=pack(drops,alive);call move_alloc(kept_mass,mass)
    if(present(tab_tangent_history)) then
      history=pack(history,alive)
      call migrate_spray_parcels(comm,drops,lower,spacing,cells,parts, &
        [faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic'],options%maximum_parcels,valid, &
        component_mass=mass,tab_tangent_history=history)
    else
      call migrate_spray_parcels(comm,drops,lower,spacing,cells,parts, &
        [faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic'],options%maximum_parcels,valid,component_mass=mass)
    endif
    if(.not.valid) return
    if(present(wall_films)) then
      call migrate_mpi_wall_films(comm,patches,patch_mass,lower,spacing,cells,parts,wall_options%maximum_films,valid)
      if(.not.valid) return
    endif
    call MPI_Allreduce(local_wall,global_wall,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if(status/=MPI_SUCCESS) return
    global_wall=global_wall-impact_wall
    call MPI_Allreduce(exterior,global_exterior,size(exterior),MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if(status/=MPI_SUCCESS) return
    valid=all(ieee_is_finite(global_wall)).and.all(ieee_is_finite(global_exterior))
    call agree();if(.not.agreed) return
    call move_alloc(drops,parcels);call move_alloc(mass,masses)
    if(present(tab_tangent_history)) call move_alloc(history,tab_tangent_history)
    if(present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_component_mass)
    endif
    last_id=watermark;wall_change=global_wall;boundary_change=global_exterior;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine
    subroutine inventory(drop,components,total,accepted)
      type(spray_parcel), intent(in) :: drop
      real(dp), intent(in) :: components(:)
      real(dp), intent(out) :: total(5)
      logical, intent(out) :: accepted
      real(dp) :: density,cp,energy
      call liquid_mixture_properties(liquids,components/drop%mass,drop%temperature,density,cp,energy,accepted)
      if(.not.accepted) return
      total=drop%mass*drop%multiplicity*[1._dp,drop%velocity,energy+.5_dp*sum(drop%velocity**2)]
      accepted=all(ieee_is_finite(total))
    end subroutine
  end subroutine drift_mpi_sdc_particles

  ! The gas integrator evaluates one old and iterations new source stages.
  ! Particle/film candidates remain local to this interval until every gas,
  ! chemistry, particle and wall stage and the final ledgers are accepted.
  subroutine advance_mpi_cartesian_spray_sdc(comm,species,reactions,transport,liquids,options, &
      state,temperature,parcels,masses,lower,spacing,cells,dt,time,iterations,force,flags,rtol,atol,difmag, &
      faces,reaction_source,heat_release,last_id,boundary_change,source_change,exchange_change, &
      particle_wall_change,particle_boundary_change,particle_force_change,ok,acceleration,evaluate_reaction, &
      breakup_properties,wall_options,wall_films,film_component_mass,native_flow_cfl, &
      inflow_state,inflow_temperature,inflow_mach,inflow_pressure,outlet_pressure,wall_temperature, &
      wall_velocity,constant_coefficients,les_parameters,les_model,les_ci,les_test_filter,les_test_ratio, &
      les_grid_ratio,les_harmonic_mean,mol_mode,mol_lax_friedrichs,frame)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),reaction_source(:,:,:,:),heat_release(:,:,:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: masses(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),dt,time,force(3),rtol,atol,difmag
    integer, intent(in) :: cells(3),iterations
    logical, intent(in) :: flags(4)
    character(len=*), intent(in) :: faces(6)
    integer(int64), intent(inout) :: last_id
    real(dp), intent(inout) :: boundary_change(:),source_change(:),exchange_change(:),particle_wall_change(5)
    real(dp), intent(inout) :: particle_boundary_change(:),particle_force_change(5)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: acceleration(3),native_flow_cfl
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(spray_options), optional, intent(in) :: breakup_properties(:)
    type(spray_wall_options), optional, intent(in) :: wall_options
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_component_mass(:,:)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: outlet_pressure(6),wall_temperature(:),wall_velocity(:,:),constant_coefficients(:)
    real(dp), optional, intent(in) :: les_parameters(3),les_ci
    character(len=*), optional, intent(in) :: les_model,les_test_filter
    integer, optional, intent(in) :: les_test_ratio,les_grid_ratio
    logical, optional, intent(in) :: les_harmonic_mean,mol_mode,mol_lax_friedrichs
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(mpi_sdc_particle_stage) :: particle_work
    type(spray_parcel), allocatable :: drops(:)
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: mass(:,:),patch_mass(:,:),gas(:,:,:,:),temp(:,:,:),ir(:,:,:,:),heat(:,:,:)
    real(dp), allocatable :: film_rate(:,:,:,:),old_film_rate(:,:,:,:),last_particle_rate(:,:,:,:),recovered(:,:,:)
    real(dp) :: wall(5),outside(size(state,1)),forced(5),boundary(size(state,1)),source(size(state,1))
    real(dp) :: film_wall(5),old_film_wall(5),old_exchange(size(state,1)),new_exchange(size(state,1))
    real(dp) :: local_exchange(size(state,1)),exchange(size(state,1)),new_rotating_exchange(size(state,1))
    real(dp) :: local_lower(3),position(3)
    real(dp), allocatable :: potential(:,:,:)
    integer :: rank,nranks,status,nv,stage,c,particle_iterations,config(4),minimum(4),maximum(4),first(3),extent(3),x,y,z
    integer(int64) :: watermark
    logical :: valid,agreed,gas_accepted,mol
    ok=.false.;nv=size(state,1);stage=0
    mol=.false.;if(present(mol_mode)) mol=mol_mode
    particle_iterations=iterations;if(mol) particle_iterations=1
    config=[merge(1,0,present(wall_options)),merge(1,0,present(wall_films)), &
      merge(1,0,present(film_component_mass)),merge(1,0,present(frame))]
    valid=allocated(parcels).and.allocated(masses).and.(present(wall_options).eqv.present(wall_films))
    valid=valid.and.(present(wall_films).eqv.present(film_component_mass))
    valid=valid.and.size(exchange_change)==nv.and.size(particle_boundary_change)==nv
    valid=valid.and.size(boundary_change)==nv.and.size(source_change)==nv
    call agree(valid,agreed);if(.not.agreed) return
    call MPI_Allreduce(config,minimum,4,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(config,maximum,4,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(minimum/=maximum)) return
    if(present(wall_films)) then
      valid=allocated(wall_films).and.allocated(film_component_mass)
      call agree(valid,agreed);if(.not.agreed) return
      patches=wall_films;patch_mass=film_component_mass
    endif
    call MPI_Comm_size(comm,nranks,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    call agree(valid,agreed);if(.not.agreed) return
    local_lower=lower+real(first,dp)*spacing
    if(present(frame)) then
      valid=frame%is_valid().and.mol
    else
      valid=.true.
    endif
    call agree(valid,agreed);if(.not.agreed) return
    if(present(frame)) then
      allocate(potential(size(temperature,1),size(temperature,2),size(temperature,3)))
      do z=1,size(potential,3)
      do y=1,size(potential,2)
      do x=1,size(potential,1)
        position=local_lower+(real([x,y,z],dp)-.5_dp)*spacing
        call rotating_potential_3d(frame,position,potential(x,y,z),valid)
        if(.not.valid) exit
      enddo
        if(.not.valid) exit
      enddo
        if(.not.valid) exit
      enddo
      call agree(valid,agreed);if(.not.agreed) return
    endif
    new_rotating_exchange=0
    gas=state;temp=temperature;ir=reaction_source;heat=heat_release;drops=parcels;mass=masses;watermark=last_id
    allocate(film_rate,mold=state);allocate(old_film_rate,mold=state)
    allocate(last_particle_rate,mold=state);last_particle_rate=0
    film_rate=0;old_film_rate=0;film_wall=0;old_film_wall=0
    old_exchange=0;new_exchange=0;wall=0;outside=0;forced=0
    call advance_mpi_cartesian_gas_sdc(comm,species,reactions,transport,gas,temp,cells,spacing,dt,time,iterations, &
      force,flags,rtol,atol,difmag,faces,ir,heat,boundary,source,gas_accepted,auxiliary_source=particle_source, &
      evaluate_reaction=evaluate_reaction,inflow_state=inflow_state,inflow_temperature=inflow_temperature, &
      inflow_mach=inflow_mach,inflow_pressure=inflow_pressure,outlet_pressure=outlet_pressure, &
      wall_temperature=wall_temperature,wall_velocity=wall_velocity,constant_coefficients=constant_coefficients, &
      les_parameters=les_parameters,les_model=les_model,les_ci=les_ci,les_test_filter=les_test_filter, &
      les_test_ratio=les_test_ratio,les_grid_ratio=les_grid_ratio,les_harmonic_mean=les_harmonic_mean, &
      mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,domain_lower=lower,frame=frame)
    if(.not.gas_accepted) return
    ! MOL completes physical MK at the first new gas stage. Later reacting
    ! passes omit its source; restore that completed exchange exactly once.
    if(mol.and.size(reactions)>0.and.iterations>1) then
      gas=gas+.5_dp*dt*last_particle_rate
      valid=all(ieee_is_finite(gas)).and.all(gas(1,:,:,:)>0).and.all(gas(6:,:,:,:)>=0)
      call agree(valid,agreed);if(.not.agreed) return
      allocate(recovered,mold=temp)
      call recover_reactive_temperatures_3d(species,gas,temp,size(gas,2),size(gas,3),size(gas,4),recovered,valid,potential)
      if(valid) temp=recovered
      call agree(valid,agreed);if(.not.agreed) return
      call MPI_Allreduce(new_rotating_exchange,exchange,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
      if(status/=MPI_SUCCESS) return
      source=source+.5_dp*dt*exchange
    endif
    valid=stage==particle_iterations+1.and.particle_work%next_stage==0
    local_exchange=.5_dp*dt*(old_exchange+new_exchange)
    valid=valid.and.all(ieee_is_finite(local_exchange)).and.all(ieee_is_finite([wall,outside,forced]))
    call agree(valid,agreed);if(.not.agreed) return
    call MPI_Allreduce(local_exchange,exchange,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if(status/=MPI_SUCCESS) return
    source=source-exchange
    valid=all(ieee_is_finite(exchange)).and.all(ieee_is_finite(source))
    call agree(valid,agreed);if(.not.agreed) return
    state=gas;temperature=temp;reaction_source=ir;heat_release=heat
    call move_alloc(drops,parcels);call move_alloc(mass,masses);last_id=watermark
    if(present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_component_mass)
    endif
    boundary_change=boundary;source_change=source;exchange_change=exchange
    particle_wall_change=wall;particle_boundary_change=outside;particle_force_change=forced;ok=.true.
  contains
    subroutine agree(local_ok,all_ok)
      logical, intent(in) :: local_ok
      logical, intent(out) :: all_ok
      integer :: code
      call MPI_Allreduce(local_ok,all_ok,1,MPI_LOGICAL,MPI_LAND,comm,code)
      if(code/=MPI_SUCCESS) all_ok=.false.
    end subroutine
    subroutine particle_source(u,t,stamp,rates,good)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),stamp
      real(dp), intent(out) :: rates(:,:,:,:,:)
      logical, intent(out) :: good
      logical :: accepted
      real(dp) :: inventory(nv)
      real(dp), allocatable :: ordinary_gas(:,:,:,:)
      rates=0
      valid=stage<=particle_iterations.and.size(rates,5)==1.and.ieee_is_finite(stamp)
      call agree(valid,good);if(.not.good) return
      ordinary_gas=u
      if(present(frame)) then
        call convert_regular_energy_3d(frame,local_lower,spacing,ordinary_gas,0,.false.,accepted)
        call agree(accepted,good);if(.not.good) return
      endif
      if(present(wall_films)) then
        if(stage==0) then
          call kick_mpi_sdc_films(comm,species,transport,liquids,options,wall_options,ordinary_gas,t,lower,spacing, &
            cells,[nranks,1,1],.5_dp*dt,faces,patches,patch_mass,film_rate,film_wall,accepted)
          if(.not.accepted) then
            good=.false.;return
          endif
          old_film_rate=film_rate;old_film_wall=film_wall
        endif
        call advance_mpi_sdc_particle_stage(comm,species,transport,liquids,options,ordinary_gas,t,lower,spacing,cells, &
          [nranks,1,1],time,dt,stage,particle_iterations,faces,drops,mass,particle_work,rates(:,:,:,:,1),wall,outside, &
          forced,watermark,accepted,acceleration,breakup_properties,wall_options,patches,patch_mass,native_flow_cfl)
      else
        call advance_mpi_sdc_particle_stage(comm,species,transport,liquids,options,ordinary_gas,t,lower,spacing,cells, &
          [nranks,1,1],time,dt,stage,particle_iterations,faces,drops,mass,particle_work,rates(:,:,:,:,1),wall,outside, &
          forced,watermark,accepted,acceleration,breakup_properties,native_flow_cfl=native_flow_cfl)
      endif
      good=accepted;if(.not.good) return
      if(present(wall_films)) then
        if(stage==particle_iterations) then
          call kick_mpi_sdc_films(comm,species,transport,liquids,options,wall_options,ordinary_gas,t,lower,spacing, &
            cells,[nranks,1,1],.5_dp*dt,faces,patches,patch_mass,film_rate,film_wall,accepted)
          if(.not.accepted) then
            good=.false.;return
          endif
          film_wall=film_wall+old_film_wall
        else
          film_rate=old_film_rate;film_wall=old_film_wall
        endif
        rates(:,:,:,:,1)=rates(:,:,:,:,1)+film_rate;wall=wall+film_wall
      endif
      do c=1,nv
        inventory(c)=sum(rates(c,:,:,:,1))*product(spacing)
      enddo
      valid=all(ieee_is_finite(inventory)).and.all(ieee_is_finite(rates)).and.all(ieee_is_finite(wall))
      call agree(valid,good);if(.not.good) return
      if(stage==0) then
        old_exchange=inventory
      else
        new_exchange=inventory
      endif
      ! Physical particle exchange uses ordinary relative energy. The published
      ! gas source subtracts Phi*drho. Leave that storage-convention change in
      ! source_change after removing physical exchange below; it closes the
      ! mixed gas/particle inventory without pretending it is particle heating.
      if(present(frame)) rates(5,:,:,:,1)=rates(5,:,:,:,1)-potential*rates(1,:,:,:,1)
      valid=all(ieee_is_finite(rates))
      call agree(valid,good);if(.not.good) return
      if(stage>0) then
        last_particle_rate=rates(:,:,:,:,1)
        do c=1,nv
          new_rotating_exchange(c)=sum(rates(c,:,:,:,1))*product(spacing)
        enddo
      endif
      stage=stage+1
    end subroutine
  end subroutine advance_mpi_cartesian_spray_sdc

  ! Owned x slabs share CTU/source halos. All state and ledger publication is
  ! collective; particle stage owners can provide the auxiliary source family.
  subroutine advance_mpi_cartesian_gas_sdc(comm,species,reactions,transport,state,temperature,cells,spacing, &
      dt,time,iterations,force,flags,rtol,atol,difmag,faces,reaction_source,heat_release,boundary_change, &
      source_change,ok,auxiliary_source,evaluate_reaction,inflow_state,inflow_temperature,inflow_mach, &
      inflow_pressure,outlet_pressure,wall_temperature,wall_velocity,constant_coefficients,les_parameters, &
      les_model,les_ci,les_test_filter,les_test_ratio,les_grid_ratio,les_harmonic_mean,mol_mode,mol_lax_friedrichs, &
      domain_lower,frame)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),reaction_source(:,:,:,:),heat_release(:,:,:)
    real(dp), intent(inout) :: boundary_change(:),source_change(:)
    integer, intent(in) :: cells(3),iterations
    real(dp), intent(in) :: spacing(3),dt,time,force(3),rtol,atol,difmag
    logical, intent(in) :: flags(4)
    character(len=*), intent(in) :: faces(6)
    logical, intent(out) :: ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    logical, optional, intent(in) :: mol_mode,mol_lax_friedrichs
    procedure(reactive_sdc_source_evaluator), optional :: auxiliary_source
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: outlet_pressure(6),wall_temperature(:),wall_velocity(:,:),constant_coefficients(:)
    real(dp), optional, intent(in) :: les_parameters(3),les_ci
    character(len=*), optional, intent(in) :: les_model,les_test_filter
    integer, optional, intent(in) :: les_test_ratio,les_grid_ratio
    logical, optional, intent(in) :: les_harmonic_mean
    real(dp) :: parameters(3),isotropic,old_les_rate(size(state,1)),new_les_rate(size(state,1))
    real(dp) :: les_inventory(size(state,1)),global_les_inventory(size(state,1))
    integer :: closure_code,filter_code,test_ratio,grid_ratio,source_calls
    character(len=24) :: closure,filter
    logical :: harmonic,les_enabled,mol,use_lax
    type(elementary_reaction), allocatable :: active_reactions(:)
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),ir(:,:,:,:),heat(:,:,:),reference(:),controls(:)
    real(dp) :: boundary(size(state,1)),source(size(state,1)),total_boundary(size(state,1)),total_source(size(state,1))
    real(dp), allocatable :: work(:,:,:,:),work_t(:,:,:),ordinary(:,:,:,:)
    real(dp) :: global_lower(3),local_lower(3),frame_values(4)
    integer :: frame_axis
    integer :: rank,nranks,status,nv,first(3),extent(3),config(33),low(33),high(33),nc,nfamilies
    character(len=32) :: reference_faces(6)
    logical :: valid,agreed,accepted,physical(6)
    ok=.false.;nv=size(species)+5;nc=0
    mol=.false.;use_lax=.false.
    global_lower=0;frame_values=0;frame_axis=0
    if(present(domain_lower)) global_lower=domain_lower
    if(present(frame)) then
      frame_axis=frame%axis;frame_values=[frame%angular_speed,frame%axis_point]
    endif
    if(present(mol_mode)) mol=mol_mode
    if(present(mol_lax_friedrichs)) use_lax=mol_lax_friedrichs
    if(present(constant_coefficients)) nc=size(constant_coefficients)
    parameters=[0._dp,.7_dp,.7_dp];isotropic=0;closure='smagorinsky';filter='box'
    test_ratio=2;grid_ratio=1;harmonic=.false.;les_enabled=present(les_parameters)
    if(present(les_parameters)) parameters=les_parameters
    if(present(les_ci)) isotropic=les_ci
    if(present(les_model)) closure=les_model
    if(present(les_test_filter)) filter=les_test_filter
    if(present(les_test_ratio)) test_ratio=les_test_ratio
    if(present(les_grid_ratio)) grid_ratio=les_grid_ratio
    if(present(les_harmonic_mean)) harmonic=les_harmonic_mean
    closure_code=0
    select case(closure)
    case('smagorinsky');closure_code=1
    case('wale');closure_code=2
    case('vreman');closure_code=3
    case('dynamic');closure_code=4
    end select
    filter_code=les_filter_code(filter)
    config=[cells,iterations,size(species),size(reactions),size(transport),merge(1,0,flags), &
      merge(1,0,present(auxiliary_source)),merge(1,0,present(evaluate_reaction)), &
      merge(1,0,present(inflow_state)),merge(1,0,present(inflow_temperature)), &
      merge(1,0,present(inflow_mach)),merge(1,0,present(inflow_pressure)),merge(1,0,present(outlet_pressure)), &
      merge(1,0,present(wall_temperature)),merge(1,0,present(wall_velocity)), &
      merge(1,0,present(constant_coefficients)),nc,merge(1,0,les_enabled),closure_code,filter_code, &
      test_ratio,grid_ratio,merge(1,0,harmonic),merge(1,0,mol),merge(1,0,use_lax),merge(1,0,present(frame)),frame_axis, &
      merge(1,0,present(domain_lower))]
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(config,low,size(config),MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(config,high,size(config),MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    call agree(all(low==high),agreed);if(.not.agreed) return
    reference_faces=faces
    call MPI_Bcast(reference_faces,6*len(reference_faces),MPI_CHARACTER,0,comm,status)
    if(status/=MPI_SUCCESS) return
    call agree(all(faces==reference_faces),agreed);if(.not.agreed) return
    allocate(controls(22+nc),reference(22+nc))
    controls(:11)=[spacing,dt,time,force,rtol,atol,difmag]
    controls(12:15)=[parameters,isotropic]
    controls(16:22)=[global_lower,frame_values]
    if(present(constant_coefficients)) controls(23:)=constant_coefficients
    reference=controls
    call MPI_Bcast(reference,size(reference),MPI_DOUBLE_PRECISION,0,comm,status)
    if(status/=MPI_SUCCESS) return
    call agree(all(ieee_is_finite(controls)).and.all(controls==reference),agreed)
    if(.not.agreed) return
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    valid=valid.and.all(cells>=2).and.all(extent>=1).and.iterations>=1.and.size(transport)==size(species)
    local_lower=global_lower+real(first,dp)*spacing
    valid=valid.and.(mol.or..not.use_lax)
    if(present(frame)) valid=valid.and.frame%is_valid().and.mol
    valid=valid.and.all(shape(state)==[nv,extent]).and.all(shape(temperature)==extent)
    valid=valid.and.all(shape(reaction_source)==[nv-4,extent]).and.all(shape(heat_release)==extent)
    valid=valid.and.size(boundary_change)==nv.and.size(source_change)==nv
    valid=valid.and.all(ieee_is_finite(state)).and.all(ieee_is_finite(temperature))
    valid=valid.and.all(ieee_is_finite(reaction_source)).and.all(ieee_is_finite(heat_release))
    valid=valid.and.(size(reactions)>0.or..not.present(evaluate_reaction))
    call agree(valid,agreed);if(.not.agreed) return
    valid=parameters(1)>=0.and.parameters(1)<=2.and.minval(parameters(2:))>=1e-6_dp
    valid=valid.and.maxval(parameters(2:))<=1e6_dp.and.isotropic>=0.and.isotropic<=1
    valid=valid.and.closure_code>0.and.filter_code>0.and.test_ratio>=1.and.grid_ratio>=1
    valid=valid.and.(closure=='dynamic'.or..not.harmonic).and.(closure/='dynamic'.or.isotropic==0)
    valid=valid.and.(les_enabled.or.(parameters(1)==0.and.isotropic==0))
    call agree(valid,agreed);if(.not.agreed) return
    ! Validate physical controls, halos and rank-consistent inlet/wall settings
    ! before entering the shared core's matched optional branches.
    ordinary=state
    if(present(frame)) then
      call convert_regular_energy_3d(frame,local_lower,spacing,ordinary,0,.false.,valid)
      call agree(valid,agreed);if(.not.agreed) return
    endif
    call fill_mpi_les_boundary_ghosts_3d(comm,species,ordinary,temperature,cells,faces,4,work,work_t,valid, &
      inflow_state,inflow_temperature,wall_velocity,outlet_pressure,inflow_mach,inflow_pressure,viscous_wall=.false.)
    if(.not.valid) return
    if(present(wall_temperature)) then
      valid=size(wall_temperature)==6
      call agree(valid,agreed);if(.not.agreed) return
      block
        real(dp) :: wall_reference(6)
        wall_reference=wall_temperature
        call MPI_Bcast(wall_reference,6,MPI_DOUBLE_PRECISION,0,comm,status)
        if(status/=MPI_SUCCESS) return
        call agree(all(ieee_is_finite(wall_temperature)).and.all(wall_temperature==wall_reference),agreed)
        if(.not.agreed) return
      end block
    endif
    gas=state;temp=temperature;heat=heat_release
    if(size(reactions)>0) then
      active_reactions=reactions;ir=reaction_source
    endif
    physical=.true.;physical(1)=rank==0;physical(2)=rank==nranks-1
    old_les_rate=0;new_les_rate=0;source_calls=0
    nfamilies=2+merge(1,0,present(auxiliary_source))+merge(1,0,les_enabled)
    if(mol) then
      call advance_reactive_mol_3d(species,gas,temp,spacing,dt,time,iterations,size(reactions)>0,rtol,atol,accepted, &
        reactions=active_reactions,reaction_source=ir,heat_release=heat,transport=transport,transport_flags=flags, &
        source_count=nfamilies-1,evaluate_sources=evaluate_sources,source_change=source,boundary_change=boundary, &
        boundary_faces=faces,inflow_state=inflow_state,inflow_temperature=inflow_temperature,inflow_mach=inflow_mach, &
        inflow_pressure=inflow_pressure,outlet_pressure=outlet_pressure,wall_temperature=wall_temperature, &
        wall_velocity=wall_velocity,constant_coefficients=constant_coefficients,evaluate_reaction=evaluate_reaction, &
        fill_halos=state_halo,agree_stage=agree,physical_faces=physical,lax_friedrichs=use_lax,domain_lower=local_lower, &
        frame=frame,agree_flux_blend=agree_window_blend)
    else
    call advance_reactive_sdc_3d(species,gas,temp,spacing,dt,time,iterations,nfamilies,evaluate_sources,accepted, &
      difmag=difmag,transport=transport,transport_flags=flags,source_change=source,boundary_change=boundary, &
      reactions=active_reactions,reaction_source=ir,heat_release=heat, &
      chemistry_relative_tolerance=rtol,chemistry_absolute_tolerance=atol,boundary_faces=faces, &
      inflow_state=inflow_state,inflow_temperature=inflow_temperature,inflow_mach=inflow_mach, &
      inflow_pressure=inflow_pressure,outlet_pressure=outlet_pressure,wall_temperature=wall_temperature, &
      wall_velocity=wall_velocity,constant_coefficients=constant_coefficients,evaluate_reaction=evaluate_reaction, &
      fill_halos=state_halo,fill_source_halos=source_halo,agree_stage=agree,physical_faces=physical)
    endif
    if(.not.accepted) return
    ! Accepted owned cells must be writable as a conservative restart state.
    valid=all(ieee_is_finite(gas)).and.all(ieee_is_finite(temp)).and.all(ieee_is_finite(heat))
    valid=valid.and.all(gas(1,:,:,:)>0).and.all(gas(6:,:,:,:)>=0).and.all(temp>0)
    if(allocated(ir)) valid=valid.and.all(ieee_is_finite(ir))
    call agree(valid,agreed);if(.not.agreed) return
    call MPI_Allreduce(boundary,total_boundary,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(source,total_source,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if(status/=MPI_SUCCESS) return
    call agree(all(ieee_is_finite(total_boundary)).and.all(ieee_is_finite(total_source)),agreed)
    if(.not.agreed) return
    if(les_enabled) then
      if(mol.and.iterations>1.and.size(reactions)>0) new_les_rate=0
      les_inventory=.5_dp*dt*(old_les_rate+new_les_rate)
      call agree(all(ieee_is_finite(les_inventory)),agreed);if(.not.agreed) return
      call MPI_Allreduce(les_inventory,global_les_inventory,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
      if(status/=MPI_SUCCESS) return
      total_boundary=total_boundary+global_les_inventory;total_source=total_source-global_les_inventory
      call agree(all(ieee_is_finite(total_boundary)).and.all(ieee_is_finite(total_source)),agreed)
      if(.not.agreed) return
    endif
    state=gas;temperature=temp;heat_release=heat
    reaction_source=0
    if(allocated(ir)) reaction_source=ir
    boundary_change=total_boundary;source_change=total_source;ok=.true.
  contains
    subroutine agree(local_ok,all_ok)
      logical, intent(in) :: local_ok
      logical, intent(out) :: all_ok
      integer :: code
      call MPI_Allreduce(local_ok,all_ok,1,MPI_LOGICAL,MPI_LAND,comm,code)
      if(code/=MPI_SUCCESS) all_ok=.false.
    end subroutine
    subroutine agree_window_blend(local_blend,global_blend,good)
      real(dp), intent(in) :: local_blend
      real(dp), intent(out) :: global_blend
      logical, intent(out) :: good
      integer :: code
      call MPI_Allreduce(local_blend,global_blend,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,code)
      good=code==MPI_SUCCESS
    end subroutine agree_window_blend
    subroutine evaluate_sources(u,t,stamp,rates,good)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),stamp
      real(dp), intent(out) :: rates(:,:,:,:,:)
      logical, intent(out) :: good
      logical :: local_ok
      rates=0
      call compute_reactive_external_source_3d(u,force,rates(:,:,:,:,1),local_ok)
      call agree(local_ok,good);if(.not.good) return
      if(les_enabled) then
        call evaluate_les(u,t,rates(:,:,:,:,2),good)
        if(.not.good) return
      endif
      if(present(auxiliary_source)) then
        if(les_enabled) then
          call auxiliary_source(u,t,stamp,rates(:,:,:,:,3:3),good)
        else
          call auxiliary_source(u,t,stamp,rates(:,:,:,:,2:2),good)
        endif
      endif
      source_calls=source_calls+1
    end subroutine
    subroutine evaluate_les(u,t,rate,good)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:)
      real(dp), intent(out) :: rate(:,:,:,:)
      logical, intent(out) :: good
      real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),weights(:),lx(:,:,:,:),ly(:,:,:,:),lz(:,:,:,:),ordinary_u(:,:,:,:)
      real(dp) :: diffusion,inventory(nv)
      character(len=32) :: local_faces(6)
      integer :: width,c
      logical :: local_ok,physical_boundary(6)
      local_faces=faces
      if(rank/=0) local_faces(1)='outflow'
      if(rank/=nranks-1) local_faces(2)='outflow'
      physical_boundary=physical.and.faces/='periodic'
      width=2
      if(closure=='dynamic') then
        call les_filter_weights(filter,test_ratio,weights,local_ok)
        call agree(local_ok,good);if(.not.good) return
        width=(size(weights)-1)/2+5
      endif
      call agree(width<=128,good);if(.not.good) return
      ordinary_u=u
      if(present(frame)) then
        call convert_regular_energy_3d(frame,local_lower,spacing,ordinary_u,0,.false.,local_ok)
        call agree(local_ok,good);if(.not.good) return
      endif
      call fill_mpi_les_boundary_ghosts_3d(comm,species,ordinary_u,t,cells,faces,width,ghost,gt,good, &
        inflow_state,inflow_temperature,wall_velocity,outlet_pressure,inflow_mach,inflow_pressure,viscous_wall=.true.)
      if(.not.good) return
      allocate(lx(nv,extent(1)+1,extent(2),extent(3)),ly(nv,extent(1),extent(2)+1,extent(3)), &
        lz(nv,extent(1),extent(2),extent(3)+1))
      inventory=0
      if(closure=='dynamic') then
        call compute_dynamic_les_rhs_ghosted_3d(species,ghost,gt,spacing,filter,test_ratio,width, &
          rate,lx,ly,lz,local_ok,grid_ratio,harmonic,diffusion,local_faces,wall_temperature)
        call agree(local_ok,good);if(.not.good) return
        do c=1,nv
          if(physical_boundary(1)) inventory(c)=inventory(c)+sum(lx(c,1,:,:))*spacing(2)*spacing(3)
          if(physical_boundary(2)) inventory(c)=inventory(c)-sum(lx(c,extent(1)+1,:,:))*spacing(2)*spacing(3)
          if(physical_boundary(3)) inventory(c)=inventory(c)+sum(ly(c,:,1,:))*spacing(1)*spacing(3)
          if(physical_boundary(4)) inventory(c)=inventory(c)-sum(ly(c,:,extent(2)+1,:))*spacing(1)*spacing(3)
          if(physical_boundary(5)) inventory(c)=inventory(c)+sum(lz(c,:,:,1))*spacing(1)*spacing(2)
          if(physical_boundary(6)) inventory(c)=inventory(c)-sum(lz(c,:,:,extent(3)+1))*spacing(1)*spacing(2)
        enddo
      else
        call compute_les_ghosted_rhs_3d(species,ghost,gt,extent,spacing,parameters(1),parameters(2),parameters(3), &
          rate,diffusion,local_ok,boundary_faces=local_faces,boundary_rate=inventory,wall_temperature=wall_temperature, &
          boundary_face_mask=physical_boundary,model=closure,ci=isotropic,fx=lx,fy=ly,fz=lz)
        call agree(local_ok,good);if(.not.good) return
      endif
      call agree(all(ieee_is_finite(inventory)),good);if(.not.good) return
      new_les_rate=inventory
      if(source_calls==0) old_les_rate=inventory
    end subroutine evaluate_les
    subroutine state_halo(u,t,stamp,width,viscous,ghost,gt,good, &
      internal_energy_density,internal_energy_ghost)
    real(dp), optional, intent(in) :: internal_energy_density(:,:,:)
    real(dp), optional, intent(out) :: internal_energy_ghost(:,:,:)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),stamp
      integer, intent(in) :: width
      logical, intent(in) :: viscous
      real(dp), intent(out) :: ghost(:,:,:,:),gt(:,:,:)
      logical, intent(out) :: good
      real(dp), allocatable :: work_h(:,:,:)
      logical :: local_ok
      ordinary=u
      if(present(frame)) then
        call convert_regular_energy_3d(frame,local_lower,spacing,ordinary,0,.false.,good)
        call agree(good,local_ok);good=local_ok
        if(.not.good) return
      endif
      if(present(internal_energy_density)) then
        if(.not.present(internal_energy_ghost)) then
          good=.false.;return
        endif
      call fill_mpi_les_boundary_ghosts_3d(comm,species,ordinary,t,cells,faces,width,work,work_t,good, &
        inflow_state,inflow_temperature,wall_velocity,outlet_pressure,inflow_mach,inflow_pressure,viscous, &
        internal_energy_density=internal_energy_density,internal_energy_ghost=work_h)
      else
      call fill_mpi_les_boundary_ghosts_3d(comm,species,ordinary,t,cells,faces,width,work,work_t,good, &
        inflow_state,inflow_temperature,wall_velocity,outlet_pressure,inflow_mach,inflow_pressure,viscous)
      endif
      if(good) then
        if(present(frame)) then
          call convert_regular_energy_3d(frame,local_lower,spacing,work,width,.true.,good)
          call agree(good,local_ok);good=local_ok
          if(.not.good) return
        endif
        if(present(internal_energy_ghost)) internal_energy_ghost=work_h
        ghost=work;gt=work_t
      endif
    end subroutine
    subroutine source_halo(rate,stamp,ghost,good)
      real(dp), intent(in) :: rate(:,:,:,:),stamp
      real(dp), intent(out) :: ghost(:,:,:,:)
      logical, intent(out) :: good
      call fill_mpi_sdc_source_ghosts_3d(comm,rate,cells,faces,4,work,good)
      if(good) ghost=work
    end subroutine
  end subroutine advance_mpi_cartesian_gas_sdc
end module spray_mpi_sdc_mod
