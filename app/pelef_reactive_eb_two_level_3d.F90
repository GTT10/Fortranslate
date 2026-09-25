program pelef_reactive_eb_two_level_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  use, intrinsic :: iso_fortran_env, only: iostat_end, iostat_eor
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
#ifndef PELEF_SELECTED_EB_TWO_LEVEL
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use transport_database_mod, only: load_h2o2_elementary_transport
#else
  use @PELEF_BUNDLE_MODULE_NAME@, only: &
    load_selected_thermo => @PELEF_BUNDLE_THERMO_LOADER_NAME@, &
    load_selected_reactions => @PELEF_BUNDLE_LOADER_NAME@, &
    load_selected_transport => @PELEF_BUNDLE_TRANSPORT_LOADER_NAME@, &
    selected_nspecies => @PELEF_BUNDLE_NSPECIES_SYMBOL@, &
    selected_nreactions => @PELEF_BUNDLE_NREACTIONS_SYMBOL@, &
    selected_integrator => @PELEF_BUNDLE_CHEMISTRY_INTEGRATOR_SYMBOL@, &
    element_count => @PELEF_BUNDLE_NELEMENTS_SYMBOL@, &
    selected_element_names => @PELEF_BUNDLE_ELEMENT_NAMES_SYMBOL@, &
    selected_element_atoms => @PELEF_BUNDLE_ELEMENT_ATOMS_SYMBOL@
  use selected_mechanism_runtime_mod, only: prepare_selected_mechanism,validate_selected_temperature_range
  use simulation_config_reactive_eb_3d_mod, only: resolve_reactive_eb_3d_selected_composition
  use reactive_eb_3d_driver_mod, only: reactive_eb_composition_integrals_3d
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions,mixture_specific_gas_constant
#endif
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use reactive_eb_3d_driver_mod, only: reactive_eb_integrals_3d,reactive_eb_element_integrals_3d
  use reactive_eb_3d_checkpoint_mod, only: write_reactive_eb_3d_checkpoint,read_reactive_eb_3d_checkpoint
  use eb_reactive_transfer_3d_mod, only: restrict_reactive_eb_state_3d
  use reactive_eb_cfl_3d_mod, only: reactive_eb_signal_speeds_3d,compute_reactive_eb_euler_cfl_timestep_3d
  use eb_reactive_transport_3d_mod, only: reactive_eb_transport_patch_timestep_3d
  use eb_reactive_two_level_3d_mod, only: advance_reactive_eb_two_level_full_3d, &
    eb_reactive_patch_3d,advance_reactive_eb_multipatch_full_3d,remap_reactive_eb_multipatch_3d, &
    tag_reactive_eb_temperature_3d
  use selected_composition_mod, only: selected_composition_max_species,selected_composition_name_length
  implicit none
#ifndef PELEF_SELECTED_EB_TWO_LEVEL
  integer, parameter :: element_count=3
#endif
  character(len=:), allocatable :: element_names(:)
  real(dp), allocatable :: element_atoms(:,:)
  real(dp), allocatable :: base_mole_fractions(:)
  character(len=:), allocatable :: bundle_identity,checkpoint_integrator
  integer :: composition_count=0
  character(len=selected_composition_name_length) :: composition_species(selected_composition_max_species)=''
  real(dp) :: composition_mole_fractions(selected_composition_max_species)=0._dp
#ifdef PELEF_SELECTED_EB_TWO_LEVEL
  character(len=24) :: transport_names(selected_nspecies)
  integer :: transport_geometries(selected_nspecies)
  real(dp) :: transport_values(5,selected_nspecies),mass_fractions(selected_nspecies)
  real(dp) :: common_tmin,common_tmax,gas_constant,regular_t,cut_t
#endif
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_eb_3d_config) :: config
  type(eb_geometry_3d) :: full,coarse
  type(eb_geometry_3d), pointer :: fine
  type(eb_reactive_patch_3d), allocatable, target :: patches(:)
  real(dp), pointer :: uf(:,:,:,:),tf(:,:,:)
  logical :: temperature_adaptive=.false.
  real(dp) :: tag_temperature=1000._dp
  integer :: tag_tile_shape(3)=[2,2,2],tag_buffer=1,tag_maximum_patches=16
  integer :: regrid_count=0,regrid_stage=0,regrid_patch_count(16)=0
  integer :: regrid_lower(3,16,16)=0,regrid_upper(3,16,16)=0
  integer :: original_count,original_lower(3,16),original_upper(3,16)
  real(dp) :: regrid_times(16)=0._dp,next_time
  integer :: patch_count=1,patch_lower(3,16)=0,patch_upper(3,16)=0
  integer :: nx=8,ny=8,nz=8,refinement_ratio=2,fine_lower(3)=[2,2,2],fine_upper(3)=[3,3,3],maximum_steps=100
  character(len=16) :: geometry_type='plane',plane_axis='z',extrusion_axis='z',riemann_solver='rusanov'
#ifdef PELEF_SELECTED_EB_TWO_LEVEL
  character(len=16) :: chemistry_integrator=selected_integrator
#else
  character(len=16) :: chemistry_integrator='explicit'
#endif
  character(len=512) :: level_set_file='',section_file='',output_file='eb_two_level.csv'
  real(dp) :: domain_lower(3)=0._dp,domain_upper(3)=1._dp,plane_position=.31_dp,time_step=1.e-8_dp,final_time=3.e-8_dp
  character(len=512) :: checkpoint_file='',restart_file=''
  integer :: checkpoint_interval_steps=0,maximum_step_retries=12
  logical :: stop_after_checkpoint=.false.
  logical :: automatic_timestep=.false.,subcycling=.false.
  real(dp) :: hydro_cfl=.5_dp,transport_cfl=.35_dp,difmag=.1_dp
  real(dp) :: rtol=2.e-7_dp,atol=1.e-12_dp
  logical :: chemistry_enabled=.true.,transport_enabled=.true.,viscosity_enabled=.true.
  logical :: thermal_conduction_enabled=.true.,species_diffusion_enabled=.true.,barodiffusion_enabled=.true.
  namelist /eb_two_level/ nx,ny,nz,refinement_ratio,fine_lower,fine_upper,maximum_steps,geometry_type,plane_axis, &
    extrusion_axis,level_set_file,section_file,output_file,domain_lower,domain_upper,plane_position,time_step,final_time, &
    rtol,atol,chemistry_enabled,transport_enabled,viscosity_enabled,thermal_conduction_enabled, &
    species_diffusion_enabled,barodiffusion_enabled,chemistry_integrator,riemann_solver,difmag, &
    checkpoint_file,restart_file,checkpoint_interval_steps,stop_after_checkpoint,maximum_step_retries, &
    automatic_timestep,hydro_cfl,transport_cfl,subcycling,patch_count,patch_lower,patch_upper, &
    composition_count,composition_species,composition_mole_fractions, &
    regrid_count,regrid_times,regrid_patch_count,regrid_lower,regrid_upper, &
    temperature_adaptive,tag_temperature,tag_tile_shape,tag_buffer,tag_maximum_patches
  real(dp), allocatable :: uc(:,:,:,:),tc(:,:,:),initial(:),current(:),scale(:), &
    boundary(:),wall(:),chemical(:),total_boundary(:),total_wall(:),total_chemical(:)
  logical, allocatable :: refined(:,:,:)
  real(dp) :: time,dt,theta,error,minimum_theta,minimum_dt
  integer :: unit,status,close_status,nvar,step,dims(3),nf(3),c,p,retry,rejected_steps
  character(len=1024) :: input_path,message,budget_path
  logical :: ok
  if (command_argument_count()/=1) error stop 'Usage: pelef_reactive_eb_two_level_3d <input.nml>'
  call get_command_argument(1,input_path)
  open(newunit=unit,file=trim(input_path),status='old',action='read',iostat=status)
  if (status/=0) error stop 'Cannot open two-level EB input'
  message=""
  read(unit,nml=eb_two_level,iostat=status,iomsg=message)
  close(unit,iostat=close_status)
  if (status/=0.or.close_status/=0) then
    write(*,'(a)') trim(message)
    error stop 'Invalid two-level EB input'
  end if
  if (min(nx,ny,nz)<2.or.refinement_ratio<2.or.maximum_steps<1) error stop 'Invalid grid or step count'
  if (any(mod([nx,ny,nz],refinement_ratio)/=0)) error stop 'Finest dimensions must divide by refinement ratio'
  dims=[nx,ny,nz]/refinement_ratio
  if (patch_count<0.or.patch_count>16) error stop 'Invalid patch count (0:16)'
  if (patch_count==0) then
    fine_lower=0;fine_upper=0;patch_lower=0;patch_upper=0
  else if (patch_count==1) then
    patch_lower(:,1)=fine_lower;patch_upper(:,1)=fine_upper
  else
    fine_lower=patch_lower(:,1);fine_upper=patch_upper(:,1)
  end if
  do p=1,patch_count
    if (any(patch_lower(:,p)<1).or.any(patch_upper(:,p)<patch_lower(:,p)).or.any(patch_upper(:,p)>dims)) &
      error stop 'Invalid patch box'
  end do
  original_count=patch_count;original_lower=patch_lower;original_upper=patch_upper
  if (temperature_adaptive) then
    if (any(regrid_patch_count/=0).or.any(regrid_lower/=0).or.any(regrid_upper/=0)) &
      error stop "Temperature adaptation cannot use prescribed boxes"
  end if
  call validate_schedule()
  allocate(patches(patch_count))
  nullify(fine,uf,tf)
  if (patch_count>0) fine=>patches(1)%geometry
  if (.not.all(ieee_is_finite([domain_lower,domain_upper,time_step,final_time,rtol,atol]))) &
    error stop 'Nonfinite two-level controls'
  if (any(domain_upper<=domain_lower).or.min(time_step,final_time,rtol,atol)<=0._dp) &
    error stop 'Invalid two-level physical controls'
  if (barodiffusion_enabled.and..not.species_diffusion_enabled.and.transport_enabled) &
    error stop 'Barodiffusion requires species diffusion'
  if(.not.ieee_is_finite(difmag).or.difmag<0) error stop 'Invalid artificial-viscosity coefficient'
  select case(trim(riemann_solver))
  case('rusanov','hllc','pelec')
  case default;error stop 'Invalid Riemann solver'
  end select
  if (chemistry_integrator/='explicit'.and.chemistry_integrator/='implicit') error stop 'Invalid chemistry integrator'
  if (len_trim(output_file)==0) error stop 'Empty output path'
  budget_path=trim(output_file)//'.budget.csv'
  if (output_file==input_path.or.budget_path==input_path.or.output_file==level_set_file.or. &
    budget_path==level_set_file.or.output_file==section_file.or.budget_path==section_file) &
    error stop 'Output aliases an input file'
  if (.not.all(ieee_is_finite([hydro_cfl,transport_cfl]))) error stop 'Nonfinite CFL controls'
  if (hydro_cfl<=0._dp.or.hydro_cfl>1._dp.or.transport_cfl<=0._dp.or.transport_cfl>.5_dp) &
    error stop 'Invalid CFL controls'
  if (maximum_step_retries<0.or.maximum_step_retries>60) error stop 'Invalid step retry limit (0:60)'
  if (checkpoint_interval_steps<0) error stop 'Invalid checkpoint interval'
  if ((checkpoint_interval_steps>0.and.len_trim(checkpoint_file)==0).or. &
    (stop_after_checkpoint.and.checkpoint_interval_steps==0)) error stop 'Invalid checkpoint policy'
  if (len_trim(restart_file)>0) then
    if (restart_file==output_file.or.restart_file==budget_path) error stop 'Output aliases restart'
  end if
  if (len_trim(checkpoint_file)>0) then
    if (checkpoint_file==input_path.or.checkpoint_file==output_file.or.checkpoint_file==budget_path.or. &
      checkpoint_file==level_set_file.or.checkpoint_file==section_file.or.checkpoint_file==restart_file) &
      error stop 'Checkpoint aliases another file'
  end if
  config%cfl=hydro_cfl;config%transport_cfl=transport_cfl
  config%redistribution='none';config%difmag=difmag
  config%maximum_steps=maximum_steps;config%final_time=final_time;config%riemann_solver=riemann_solver
  config%chemistry_enabled=chemistry_enabled;config%transport_enabled=transport_enabled
  config%chemistry_relative_tolerance=rtol;config%chemistry_absolute_tolerance=atol
  config%viscosity_enabled=viscosity_enabled;config%thermal_conduction_enabled=thermal_conduction_enabled
  config%species_diffusion_enabled=species_diffusion_enabled;config%barodiffusion_enabled=barodiffusion_enabled
  config%nx=nx;config%ny=ny;config%nz=nz
  config%x_lower=domain_lower(1);config%y_lower=domain_lower(2);config%z_lower=domain_lower(3)
  config%x_upper=domain_upper(1);config%y_upper=domain_upper(2);config%z_upper=domain_upper(3)
  config%geometry_type=geometry_type;config%plane_axis=plane_axis;config%plane_position=plane_position
  config%extrusion_axis=extrusion_axis;config%section_file=section_file;config%level_set_file=level_set_file
  call build_configured_eb_geometry_3d(config,full,ok)
  if (.not.ok) error stop 'Cannot construct finest EB geometry'
  call coarsen_eb_geometry_3d(full,refinement_ratio,coarse,ok)
  if (.not.ok) error stop 'Cannot represent coarse EB geometry'
  do p=1,patch_count
    call extract_eb_geometry_3d(full,(patch_lower(:,p)-1)*refinement_ratio+1, &
      patch_upper(:,p)*refinement_ratio,patches(p)%geometry,ok)
    if (.not.ok) error stop 'Cannot extract fine EB patch'
  end do
#ifdef PELEF_SELECTED_EB_TWO_LEVEL
  call load_selected_thermo(species,ok)
  if (.not.ok) error stop 'Selected thermodynamic load failed'
  call load_selected_reactions(reactions,ok)
  if (.not.ok) error stop 'Selected reaction load failed'
  call load_selected_transport(transport_names,transport_geometries,transport_values,ok)
  if (.not.ok) error stop 'Selected transport load failed'
  call prepare_selected_mechanism('Two-level EB',selected_nspecies,selected_nreactions,species,reactions, &
    transport_names,transport_geometries,transport_values,transport,common_tmin,common_tmax,ok,message)
  if (.not.ok) then
    write(*,'(a)') trim(message)
    error stop 'Invalid selected mechanism'
  end if
  element_names=selected_element_names; element_atoms=selected_element_atoms
  config%thermo_model='selected'
  config%composition_count=composition_count;config%composition_species=composition_species
  config%composition_mole_fractions=composition_mole_fractions
  allocate(base_mole_fractions(size(species)))
  call resolve_reactive_eb_3d_selected_composition(config,species,base_mole_fractions,ok,message)
  if (.not.ok) then
    write(*,'(a)') trim(message)
    error stop 'Invalid selected composition'
  end if
  call mass_fractions_from_mole_fractions(species,base_mole_fractions,mass_fractions,ok)
  if (.not.ok) error stop 'Invalid selected mass fractions'
  gas_constant=mixture_specific_gas_constant(species,mass_fractions,ok)
  if (.not.ok) error stop 'Invalid selected gas constant'
  regular_t=config%initial_pressure/(config%regular_density*gas_constant)
  cut_t=config%initial_pressure/(config%cut_density*gas_constant)
  call validate_selected_temperature_range('Two-level EB',min(regular_t,cut_t),max(regular_t,cut_t), &
    common_tmin,common_tmax,ok,message)
  if (.not.ok) then
    write(*,'(a)') trim(message)
    error stop 'Selected initial temperature outside mechanism range'
  end if
  bundle_identity='@PELEF_BUNDLE_BUNDLE_SHA256@'
  checkpoint_integrator=chemistry_integrator
#else
  if (composition_count/=0.or.any(composition_species/='').or.any(composition_mole_fractions/=0._dp)) &
    error stop 'Composition inputs require selected executable'
  call load_h2o2_elementary_thermo(species,ok)
  if (.not.ok) error stop 'Thermodynamic load failed'
  call load_h2o2_elementary_mechanism(reactions,ok)
  if (.not.ok) error stop 'Mechanism load failed'
  call load_h2o2_elementary_transport(transport,ok)
  if (.not.ok) error stop 'Transport load failed'
#endif
  nvar=reactive_nvar(size(species))
  allocate(uc(nvar,dims(1),dims(2),dims(3)),tc(dims(1),dims(2),dims(3)), &
    refined(dims(1),dims(2),dims(3)))
  config%nx=dims(1);config%ny=dims(2);config%nz=dims(3)
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,uc,tc,ok,base_mole_fractions)
  if (.not.ok) error stop 'Coarse initialization failed'
  refined=.false.
  do p=1,patch_count
    nf=[patches(p)%geometry%nx,patches(p)%geometry%ny,patches(p)%geometry%nz]
    allocate(patches(p)%state(nvar,nf(1),nf(2),nf(3)),patches(p)%temperature(nf(1),nf(2),nf(3)))
    config%nx=nf(1);config%ny=nf(2);config%nz=nf(3)
    call initialize_reactive_eb_density_sheet_3d(species,config,patches(p)%geometry, &
      patches(p)%state,patches(p)%temperature,ok,base_mole_fractions)
    if (.not.ok) error stop 'Fine initialization failed'
    if (any(refined(patch_lower(1,p):patch_upper(1,p),patch_lower(2,p):patch_upper(2,p), &
      patch_lower(3,p):patch_upper(3,p)))) error stop 'Overlapping patch boxes'
    refined(patch_lower(1,p):patch_upper(1,p),patch_lower(2,p):patch_upper(2,p), &
      patch_lower(3,p):patch_upper(3,p))=.true.
  end do
  call bind_first_patch()
  allocate(initial(nvar),current(nvar),scale(nvar),boundary(nvar),wall(nvar),chemical(nvar), &
    total_boundary(nvar),total_wall(nvar),total_chemical(nvar))
  call inventory(initial)
  do c=1,nvar
    scale(c)=coarse%dx*coarse%dy*coarse%dz*sum(abs(uc(c,:,:,:))*coarse%volume_fraction,mask=.not.refined)
    do p=1,patch_count
      scale(c)=scale(c)+patches(p)%geometry%dx*patches(p)%geometry%dy*patches(p)%geometry%dz* &
        sum(abs(patches(p)%state(c,:,:,:))*patches(p)%geometry%volume_fraction)
    end do
    scale(c)=max(1.e-30_dp,scale(c))
  end do
  total_boundary=0._dp;total_wall=0._dp;total_chemical=0._dp
  time=0._dp;step=0;minimum_theta=1._dp;minimum_dt=0._dp;rejected_steps=0
  if (len_trim(restart_file)>0) call restore_checkpoint()
  do while(time<final_time)
    if (step>=maximum_steps) error stop 'Maximum steps reached before final time'
    dt=min(time_step,final_time-time)
    if (regrid_stage<regrid_count) dt=min(dt,regrid_times(regrid_stage+1)-time)
    if (automatic_timestep) call limit_timestep(dt)
    do retry=0,maximum_step_retries
      if (time+dt<=time) error stop 'Time step cannot advance time'
      ! The full update is transactional: a rejected trial leaves both grids,
      ! per-step ledgers and the limiter untouched. Count only accepted steps.
      if (patch_count==1) then
      call advance_reactive_eb_two_level_full_3d(species,reactions,transport,coarse,fine,fine_lower,refinement_ratio, &
        riemann_solver,dt,chemistry_enabled,rtol,atol,transport_enabled,viscosity_enabled,thermal_conduction_enabled, &
        species_diffusion_enabled,barodiffusion_enabled,uc,tc,uf,tf,boundary,wall,chemical,theta,ok, &
        chemistry_integrator,subcycling,difmag=difmag)
      else
        call advance_reactive_eb_multipatch_full_3d(species,reactions,transport,coarse,patch_lower(:,1:patch_count), &
          refinement_ratio,riemann_solver,dt,chemistry_enabled,rtol,atol,transport_enabled,viscosity_enabled, &
          thermal_conduction_enabled,species_diffusion_enabled,barodiffusion_enabled,uc,tc,patches, &
          boundary,wall,chemical,theta,ok,chemistry_integrator,subcycling,difmag=difmag)
        call bind_first_patch()
      end if
      if (ok) exit
      rejected_steps=rejected_steps+1
      if (retry==maximum_step_retries) error stop 'Two-level EB step retry limit exhausted'
      dt=0.5_dp*dt
    end do
    total_boundary=total_boundary+boundary;total_wall=total_wall+wall;total_chemical=total_chemical+chemical
    next_time=time+dt
    if (regrid_stage<regrid_count) then
      if (dt==regrid_times(regrid_stage+1)-time) next_time=regrid_times(regrid_stage+1)
    end if
    time=next_time;step=step+1;minimum_theta=min(minimum_theta,theta)
    if (step==1) minimum_dt=dt
    minimum_dt=min(minimum_dt,dt)
    if (regrid_stage<regrid_count) then
      if (time>=regrid_times(regrid_stage+1)) then
        if (temperature_adaptive) call select_temperature_layout(regrid_stage+1)
        call activate_layout(regrid_stage+1,.true.)
        regrid_stage=regrid_stage+1
        call check_balance()
      end if
    end if
    if (checkpoint_interval_steps>0) then
      if (mod(step,checkpoint_interval_steps)==0) then
        call save_checkpoint()
        if (stop_after_checkpoint) exit
      end if
    end if
  end do
  call check_balance()
  call write_output()
  write(*,'(a,i0)') 'Rejected trials this invocation: ',rejected_steps
  write(*,'(a,i0)') 'Completed steps: ',step
  if (regrid_count>0) write(*,'(a,i0)') 'Completed regrids: ',regrid_stage
  write(*,'(a,es24.16)') 'Final time: ',time
  write(*,'(a,es24.16)') 'Composite inventory error: ',error
  write(*,'(a,es24.16)') 'Minimum transport theta: ',minimum_theta
  write(*,'(a,es24.16)') 'Minimum accepted dt: ',minimum_dt
contains
  subroutine bind_first_patch()
    nullify(fine,uf,tf)
    if (patch_count>0) then
      fine=>patches(1)%geometry;uf=>patches(1)%state;tf=>patches(1)%temperature
    end if
  end subroutine bind_first_patch

  subroutine validate_schedule()
    integer :: event,a,b,n
    if (regrid_count<0.or.regrid_count>16) error stop 'Invalid regrid_count (0:16)'
    if (temperature_adaptive) then
      if (regrid_count==0) error stop 'Temperature adaptation needs event times'
      if (.not.ieee_is_finite(tag_temperature)) error stop 'Nonfinite tag temperature'
      if (tag_temperature<=0._dp.or.any(tag_tile_shape<1).or.tag_buffer<0.or. &
        tag_maximum_patches<1.or.tag_maximum_patches>16) error stop 'Invalid temperature tagging controls'
    end if
    do event=1,regrid_count
      if (.not.ieee_is_finite(regrid_times(event))) error stop 'Nonfinite regrid time'
      if (regrid_times(event)<=0._dp) error stop 'Nonpositive regrid time'
      if (event>1) then
        if (regrid_times(event)<=regrid_times(event-1)) error stop 'Regrid times must increase'
      end if
      n=regrid_patch_count(event)
      if (n<0.or.n>16) error stop 'Invalid scheduled patch count (0:16)'
      do a=1,n
        if (any(regrid_lower(:,a,event)<1).or.any(regrid_upper(:,a,event)>dims).or. &
          any(regrid_upper(:,a,event)<regrid_lower(:,a,event))) error stop 'Invalid scheduled patch box'
        do b=1,a-1
          if (all(regrid_lower(:,a,event)<=regrid_upper(:,b,event)).and. &
            all(regrid_lower(:,b,event)<=regrid_upper(:,a,event))) error stop 'Overlapping scheduled patches'
        end do
      end do
    end do
  end subroutine validate_schedule

  subroutine select_temperature_layout(event)
    integer, intent(in) :: event
    integer, allocatable :: lower(:,:),upper(:,:)
    integer :: n
    logical :: valid
    call tag_reactive_eb_temperature_3d(coarse,tc,patch_lower(:,1:patch_count),refinement_ratio,patches, &
      tag_temperature,tag_tile_shape,tag_buffer,tag_maximum_patches,lower,upper,valid)
    if (.not.valid) error stop 'Temperature tagging failed or exceeded patch capacity'
    n=size(lower,2)
    regrid_patch_count(event)=n;regrid_lower(:,:,event)=0;regrid_upper(:,:,event)=0
    regrid_lower(:,1:n,event)=lower;regrid_upper(:,1:n,event)=upper
    write(*,'(a,i0,a,i0)') 'Temperature regrid event ',event,': patches=',n
  end subroutine select_temperature_layout

  subroutine validate_tagged_layout(lower,upper)
    integer, intent(in) :: lower(:,:),upper(:,:)
    integer :: a,axis,tile(3)
    tile=min(tag_tile_shape,dims)
    do a=1,size(lower,2)
      if (any(lower(:,a)<1).or.any(lower(:,a)>dims)) error stop 'Tagged checkpoint box outside domain'
      if (any(mod(lower(:,a)-1,tile)/=0)) error stop 'Tagged checkpoint box is not tile aligned'
      if (any(upper(:,a)/=lower(:,a)+min(tile-1,dims-lower(:,a)))) error stop 'Tagged checkpoint extent mismatch'
      if (a>1) then
        do axis=3,1,-1
          if (lower(axis,a)==lower(axis,a-1)) cycle
          if (lower(axis,a)<lower(axis,a-1)) error stop 'Tagged checkpoint order mismatch'
          exit
        end do
        if (all(lower(:,a)==lower(:,a-1))) error stop 'Duplicate tagged checkpoint tile'
      end if
    end do
  end subroutine validate_tagged_layout

  subroutine activate_layout(event,remapping)
    integer, intent(in) :: event
    logical, intent(in) :: remapping
    type(eb_reactive_patch_3d), allocatable, target :: candidates(:)
    integer :: a,n,d(3)
    logical :: valid
    n=regrid_patch_count(event)
    allocate(candidates(n))
    do a=1,n
      call extract_eb_geometry_3d(full,(regrid_lower(:,a,event)-1)*refinement_ratio+1, &
        regrid_upper(:,a,event)*refinement_ratio,candidates(a)%geometry,valid)
      if (.not.valid) error stop 'Cannot extract scheduled patch'
      d=[candidates(a)%geometry%nx,candidates(a)%geometry%ny,candidates(a)%geometry%nz]
      allocate(candidates(a)%state(nvar,d(1),d(2),d(3)),candidates(a)%temperature(d(1),d(2),d(3)))
      candidates(a)%state=0._dp;candidates(a)%temperature=0._dp
    end do
    if (remapping) then
      call remap_reactive_eb_multipatch_3d(species,coarse,patch_lower(:,1:patch_count),refinement_ratio, &
        patches,regrid_lower(:,1:n,event),refinement_ratio,candidates,uc,tc,valid)
      if (.not.valid) error stop 'Scheduled EB remap failed'
    end if
    call move_alloc(candidates,patches)
    patch_count=n;patch_lower=0;patch_upper=0
    patch_lower(:,1:n)=regrid_lower(:,1:n,event);patch_upper(:,1:n)=regrid_upper(:,1:n,event)
    fine_lower=patch_lower(:,1);fine_upper=patch_upper(:,1)
    call bind_first_patch()
    refined=.false.
    do a=1,n
      refined(patch_lower(1,a):patch_upper(1,a),patch_lower(2,a):patch_upper(2,a), &
        patch_lower(3,a):patch_upper(3,a))=.true.
    end do
  end subroutine activate_layout

  subroutine schedule_record(stream,reading)
    integer, intent(in) :: stream
    logical, intent(in) :: reading
    integer :: ios,n,stage,a,count_saved,lower_saved(3,16),upper_saved(3,16)
    real(dp) :: times(16),stored_threshold
    integer :: stored_tile(3),stored_buffer,stored_capacity
    character(len=64) :: marker
    if (.not.reading) then
      write(stream,'(a)',iostat=ios) 'REGRID_SCHEDULE'
      if (ios/=0) error stop 'Schedule marker write failed'
      write(stream,*,iostat=ios) regrid_count,regrid_stage,original_count
      if (ios/=0) error stop 'Schedule header write failed'
      if (temperature_adaptive) then
        write(stream,'(es26.17e3,5(1x,i0))',iostat=ios) tag_temperature,tag_tile_shape,tag_buffer,tag_maximum_patches
        if (ios/=0) error stop 'Tag policy write failed'
      end if
      write(stream,*,iostat=ios) original_lower(:,1:original_count),original_upper(:,1:original_count)
      if (ios/=0) error stop 'Initial topology write failed'
      write(stream,'(*(es26.17e3,1x))',iostat=ios) regrid_times(1:regrid_count)
      if (ios/=0) error stop 'Schedule times write failed'
      do a=1,regrid_count
        n=regrid_patch_count(a)
        write(stream,*,iostat=ios) n,regrid_lower(:,1:n,a),regrid_upper(:,1:n,a)
        if (ios/=0) error stop 'Schedule topology write failed'
      end do
      return
    end if
    read(stream,'(a)',iostat=ios) marker
    if (ios/=0.or.marker/='REGRID_SCHEDULE') error stop 'Invalid schedule marker'
    read(stream,*,iostat=ios) n,stage,count_saved
    if (ios/=0) error stop 'Missing schedule header'
    if (n/=regrid_count.or.stage<0.or.stage>n.or.count_saved/=original_count) error stop 'Schedule header mismatch'
    if (temperature_adaptive) then
      read(stream,*,iostat=ios) stored_threshold,stored_tile,stored_buffer,stored_capacity
      if (ios/=0) error stop 'Missing tag policy'
      if (.not.ieee_is_finite(stored_threshold)) error stop 'Nonfinite stored tag threshold'
      if (stored_threshold/=tag_temperature.or.any(stored_tile/=tag_tile_shape).or. &
        stored_buffer/=tag_buffer.or.stored_capacity/=tag_maximum_patches) error stop 'Temperature tag policy mismatch'
    end if
    read(stream,*,iostat=ios) lower_saved(:,1:count_saved),upper_saved(:,1:count_saved)
    if (ios/=0) error stop 'Missing initial topology'
    if (any(lower_saved(:,1:count_saved)/=original_lower(:,1:count_saved)).or. &
      any(upper_saved(:,1:count_saved)/=original_upper(:,1:count_saved))) error stop 'Initial topology mismatch'
    read(stream,*,iostat=ios) times(1:n)
    if (ios/=0) error stop 'Missing schedule times'
    if (any(times(1:n)/=regrid_times(1:n))) error stop 'Schedule times mismatch'
    do a=1,regrid_count
      ! Read the count first without letting an untrusted value set bounds.
      read(stream,*,iostat=ios) count_saved
      if (ios/=0) error stop 'Missing scheduled topology'
      if (temperature_adaptive) then
        if (count_saved<0.or.count_saved>tag_maximum_patches) error stop 'Tagged checkpoint capacity mismatch'
        if (a>stage.and.count_saved/=0) error stop 'Future tagged topology must be empty'
      else
        if (count_saved/=regrid_patch_count(a)) error stop 'Scheduled patch count mismatch'
      end if
      backspace(stream,iostat=ios)
      if (ios/=0) error stop 'Cannot reread scheduled topology'
      read(stream,*,iostat=ios) count_saved,lower_saved(:,1:count_saved),upper_saved(:,1:count_saved)
      if (ios/=0) error stop 'Invalid scheduled topology'
      if (temperature_adaptive) then
        call validate_tagged_layout(lower_saved(:,1:count_saved),upper_saved(:,1:count_saved))
        regrid_patch_count(a)=count_saved
        regrid_lower(:,1:count_saved,a)=lower_saved(:,1:count_saved)
        regrid_upper(:,1:count_saved,a)=upper_saved(:,1:count_saved)
      else
        if (any(lower_saved(:,1:count_saved)/=regrid_lower(:,1:count_saved,a)).or. &
          any(upper_saved(:,1:count_saved)/=regrid_upper(:,1:count_saved,a))) error stop 'Scheduled topology mismatch'
      end if
    end do
    call validate_schedule()
    regrid_stage=stage
    if (stage>0) call activate_layout(stage,.false.)
  end subroutine schedule_record

  function checkpoint_magic() result(magic)
    character(len=64) :: magic
    ! V1/V2 retain common-step fixed/automatic compatibility. V3/V4 select
    ! subcycling with the same layout; patch configs bind the CFL controls.
    magic='PELEF_EB_TWO_LEVEL_CHECKPOINT_V1'
    if (automatic_timestep) magic='PELEF_EB_TWO_LEVEL_CHECKPOINT_V2'
    if (subcycling) then
      magic='PELEF_EB_TWO_LEVEL_CHECKPOINT_V3'
      if (automatic_timestep) magic='PELEF_EB_TWO_LEVEL_CHECKPOINT_V4'
    end if
    if (patch_count>1) then
      magic='PELEF_EB_MULTIPATCH_CHECKPOINT_V1'
      if (automatic_timestep) magic='PELEF_EB_MULTIPATCH_CHECKPOINT_V2'
      if (subcycling) then
        magic='PELEF_EB_MULTIPATCH_CHECKPOINT_V3'
        if (automatic_timestep) magic='PELEF_EB_MULTIPATCH_CHECKPOINT_V4'
      end if
    end if
    if (patch_count==0) then
      magic='PELEF_EB_COARSE_CHECKPOINT_V1'
      if (automatic_timestep) magic='PELEF_EB_COARSE_CHECKPOINT_V2'
      if (subcycling) then
        magic='PELEF_EB_COARSE_CHECKPOINT_V3'
        if (automatic_timestep) magic='PELEF_EB_COARSE_CHECKPOINT_V4'
      end if
    end if
    if (regrid_count>0) then
      magic='PELEF_EB_SCHEDULED_CHECKPOINT_V1'
      if (automatic_timestep) magic='PELEF_EB_SCHEDULED_CHECKPOINT_V2'
      if (subcycling) then
        magic='PELEF_EB_SCHEDULED_CHECKPOINT_V3'
        if (automatic_timestep) magic='PELEF_EB_SCHEDULED_CHECKPOINT_V4'
      end if
    end if
    if (temperature_adaptive) then
      magic='PELEF_EB_TEMPERATURE_CHECKPOINT_V1'
      if (automatic_timestep) magic='PELEF_EB_TEMPERATURE_CHECKPOINT_V2'
      if (subcycling) then
        magic='PELEF_EB_TEMPERATURE_CHECKPOINT_V3'
        if (automatic_timestep) magic='PELEF_EB_TEMPERATURE_CHECKPOINT_V4'
      end if
    end if
  end function checkpoint_magic

  subroutine limit_timestep(interval)
    real(dp), intent(inout) :: interval
    real(dp), allocatable :: parent(:,:,:,:),parent_t(:,:,:)
    real(dp) :: coarse_speeds(4),fine_speeds(4),donors(4),coarse_dt,fine_dt,coarse_d,fine_d,donor_d,unused_d,factor
    logical :: valid,fine_active
    integer :: patch_index
    real(dp) :: extra_speeds(4),extra_dt,extra_d
    ! Match the parent averaging used at the start of the coupled update;
    ! estimation must not publish state changes of its own.
    parent=uc;parent_t=tc
    if (patch_count>0) then
    call restrict_reactive_eb_state_3d(species,fine,refinement_ratio,uf,tf, &
      parent(:,fine_lower(1):fine_upper(1),fine_lower(2):fine_upper(2),fine_lower(3):fine_upper(3)), &
      parent_t(fine_lower(1):fine_upper(1),fine_lower(2):fine_upper(2),fine_lower(3):fine_upper(3)),valid)
    if (.not.valid) error stop 'Cannot average parents for timestep selection'
    end if
    do patch_index=2,patch_count
      call restrict_reactive_eb_state_3d(species,patches(patch_index)%geometry,refinement_ratio, &
        patches(patch_index)%state,patches(patch_index)%temperature, &
        parent(:,patch_lower(1,patch_index):patch_upper(1,patch_index), &
          patch_lower(2,patch_index):patch_upper(2,patch_index),patch_lower(3,patch_index):patch_upper(3,patch_index)), &
        parent_t(patch_lower(1,patch_index):patch_upper(1,patch_index), &
          patch_lower(2,patch_index):patch_upper(2,patch_index),patch_lower(3,patch_index):patch_upper(3,patch_index)),valid)
      if (.not.valid) error stop 'Cannot average additional timestep parent'
    end do
    call reactive_eb_signal_speeds_3d(species,parent,parent_t,coarse,coarse_speeds,valid)
    if (.not.valid) error stop 'Invalid coarse timestep state'
    fine_active=.false.
    if (patch_count>0) fine_active=any(fine%volume_fraction>0._dp)
    fine_speeds=0._dp
    if (fine_active) then
      call reactive_eb_signal_speeds_3d(species,uf,tf,fine,fine_speeds,valid)
      if (.not.valid) error stop 'Invalid fine timestep state'
    end if
    donors=max(coarse_speeds,fine_speeds)
    do patch_index=2,patch_count
      if (.not.any(patches(patch_index)%geometry%volume_fraction>0._dp)) cycle
      call reactive_eb_signal_speeds_3d(species,patches(patch_index)%state,patches(patch_index)%temperature, &
        patches(patch_index)%geometry,extra_speeds,valid)
      if (.not.valid) error stop 'Invalid additional timestep donor'
      donors=max(donors,extra_speeds)
    end do
    call compute_reactive_eb_euler_cfl_timestep_3d(species,parent,parent_t,coarse,hydro_cfl,coarse_dt,valid,donors)
    if (.not.valid) error stop 'Cannot select coarse hydro interval'
    fine_dt=huge(1._dp)
    if (fine_active) then
      call compute_reactive_eb_euler_cfl_timestep_3d(species,uf,tf,fine,hydro_cfl,fine_dt,valid,donors)
      if (.not.valid) error stop 'Cannot select fine hydro interval'
    end if
    do patch_index=2,patch_count
      if (.not.any(patches(patch_index)%geometry%volume_fraction>0._dp)) cycle
      call compute_reactive_eb_euler_cfl_timestep_3d(species,patches(patch_index)%state, &
        patches(patch_index)%temperature,patches(patch_index)%geometry,hydro_cfl,extra_dt,valid,donors)
      if (.not.valid) error stop 'Cannot select additional hydro interval'
      fine_dt=min(fine_dt,extra_dt)
    end do
    interval=min(interval,coarse_dt)
    factor=1._dp
    if (subcycling) factor=real(refinement_ratio,dp)
    if (fine_dt<interval/factor) interval=fine_dt*factor
    if (.not.transport_enabled) return
    call reactive_eb_transport_patch_timestep_3d(species,transport,parent,parent_t,coarse,transport_cfl, &
      viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,coarse_dt,coarse_d,valid)
    if (.not.valid) error stop 'Cannot evaluate coarse diffusivity'
    fine_d=0._dp
    if (fine_active) then
      call reactive_eb_transport_patch_timestep_3d(species,transport,uf,tf,fine,transport_cfl, &
        viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,fine_dt,fine_d,valid)
      if (.not.valid) error stop 'Cannot evaluate fine diffusivity'
    end if
    donor_d=max(coarse_d,fine_d)
    do patch_index=2,patch_count
      if (.not.any(patches(patch_index)%geometry%volume_fraction>0._dp)) cycle
      call reactive_eb_transport_patch_timestep_3d(species,transport,patches(patch_index)%state, &
        patches(patch_index)%temperature,patches(patch_index)%geometry,transport_cfl, &
        viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,extra_dt,extra_d,valid)
      if (.not.valid) error stop 'Cannot evaluate additional diffusivity'
      donor_d=max(donor_d,extra_d)
    end do
    call reactive_eb_transport_patch_timestep_3d(species,transport,parent,parent_t,coarse,transport_cfl, &
      viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,coarse_dt,unused_d,valid,donor_d)
    if (.not.valid) error stop 'Cannot select coarse diffusion interval'
    fine_dt=huge(1._dp)
    if (fine_active) then
      call reactive_eb_transport_patch_timestep_3d(species,transport,uf,tf,fine,transport_cfl, &
        viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,fine_dt,unused_d,valid,donor_d)
      if (.not.valid) error stop 'Cannot select fine diffusion interval'
    end if
    do patch_index=2,patch_count
      if (.not.any(patches(patch_index)%geometry%volume_fraction>0._dp)) cycle
      call reactive_eb_transport_patch_timestep_3d(species,transport,patches(patch_index)%state, &
        patches(patch_index)%temperature,patches(patch_index)%geometry,transport_cfl, &
        viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,extra_dt,extra_d,valid,donor_d)
      if (.not.valid) error stop 'Cannot select additional diffusion interval'
      fine_dt=min(fine_dt,extra_dt)
    end do
    ! Bound the full interval, conservatively also bounding each transport half.
    interval=min(interval,coarse_dt)
    factor=1._dp
    if (subcycling) factor=real(refinement_ratio,dp)**2
    if (fine_dt<interval/factor) interval=fine_dt*factor
  end subroutine limit_timestep

#ifdef PELEF_SELECTED_EB_TWO_LEVEL
  logical function chemical_elements_valid(initial_values,changes,accepted_steps) result(valid)
    real(dp), intent(in) :: initial_values(:),changes(:)
    integer, intent(in) :: accepted_steps
    real(dp) :: element_scale(element_count),element_change(element_count),roundoff(element_count)
    integer :: e
    element_scale=matmul(element_atoms,abs(initial_values(6:))/species%molecular_weight)
    element_change=matmul(element_atoms,changes(6:)/species%molecular_weight)
    ! A zero initial element still sees floating-point species mass closure.
    ! Bound that absolute error by total species mass and atoms per unit mass;
    ! retain the relative tolerance for each populated elemental inventory.
    do e=1,element_count
      roundoff(e)=64._dp*epsilon(1._dp)*real(max(1,accepted_steps),dp)*sum(abs(initial_values(6:)))* &
        maxval(element_atoms(e,:)/species%molecular_weight)
    end do
    valid=all(ieee_is_finite(element_change)).and.all(ieee_is_finite(element_scale)).and.all(ieee_is_finite(roundoff))
    if (valid) valid=all(abs(element_change)<=2.e-10_dp*element_scale+roundoff)
    if (.not.valid) then
      do e=1,element_count
        write(*,'(a,1x,2(es24.16,1x))') trim(element_names(e)),element_change(e),element_scale(e)
      end do
    end if
  end function chemical_elements_valid
#endif

  subroutine check_balance()
#ifdef PELEF_SELECTED_EB_TWO_LEVEL
    if (.not.chemical_elements_valid(initial,total_chemical,step)) error stop 'Chemical element conservation failed'
#endif
    call inventory(current)
    ! Use initial L1, final and exchange scales with an absolute floor of one.
    error=maxval(abs(current-initial-total_boundary-total_wall-total_chemical)/ &
      max(1._dp,scale,abs(current),abs(total_boundary)+abs(total_wall)+abs(total_chemical)))
    if (.not.ieee_is_finite(error)) error stop 'Nonfinite composite inventory residual'
    if (error>5.e-10_dp) error stop 'Composite inventory balance failed'
  end subroutine check_balance

  function patch_config(geometry) result(result_config)
    type(eb_geometry_3d), intent(in) :: geometry
    type(reactive_eb_3d_config) :: result_config
    result_config=config
    result_config%nx=geometry%nx;result_config%ny=geometry%ny;result_config%nz=geometry%nz
    result_config%x_lower=geometry%x_lower;result_config%x_upper=geometry%x_upper
    result_config%y_lower=geometry%y_lower;result_config%y_upper=geometry%y_upper
    result_config%z_lower=geometry%z_lower;result_config%z_upper=geometry%z_upper
  end function patch_config

  subroutine patch_record(stream,reading,geometry,state,temperature,record_time,record_step,record_dt,record_theta)
    integer, intent(in) :: stream
    logical, intent(in) :: reading
    type(eb_geometry_3d), intent(inout) :: geometry
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: record_time,record_dt,record_theta
    integer, intent(in) :: record_step
    type(reactive_eb_3d_config) :: cfg
    real(dp) :: integrals(nvar),l1(nvar),elements(element_count),got_time,got_dt,got_theta,diffusivity
    real(dp) :: expected(nvar),expected_l1(nvar),expected_elements(element_count)
    integer :: got_step
    logical :: valid
    character(len=512) :: reason
    cfg=patch_config(geometry)
    if (reading) then
      call read_reactive_eb_3d_checkpoint('',species,reactions,transport,cfg,geometry,state,temperature, &
        got_time,got_step,integrals,l1,elements,got_dt,diffusivity,got_theta,valid,reason,bundle_sha256=bundle_identity, &
        chemistry_integrator=checkpoint_integrator,base_mole_fractions=base_mole_fractions,checkpoint_unit=stream, &
        element_names=element_names,element_atoms=element_atoms)
      if (.not.valid) then
        write(*,'(a)') trim(reason)
        error stop 'Invalid two-level patch record'
      end if
      if (got_time/=record_time.or.got_step/=record_step.or.got_dt/=record_dt.or.got_theta/=record_theta.or. &
        diffusivity/=0._dp) error stop 'Inconsistent patch metadata'
    end if
    ! Embedded diagnostic inventories describe this snapshot, not the composite
    ! initial condition. The outer record owns initial and cumulative budgets.
    call reactive_eb_integrals_3d(state,geometry,expected,valid,expected_l1)
    if (.not.valid) error stop 'Cannot integrate checkpoint patch'
#ifdef PELEF_SELECTED_EB_TWO_LEVEL
    call reactive_eb_composition_integrals_3d(species,state,geometry,element_atoms,expected_elements,valid)
#else
    call reactive_eb_element_integrals_3d(species,state,geometry,expected_elements,valid)
#endif
    if (.not.valid) error stop 'Cannot integrate checkpoint elements'
    if (reading) then
      if (any(abs(integrals-expected)>2.e-12_dp*max(1._dp,abs(expected))).or. &
        any(abs(l1-expected_l1)>2.e-12_dp*max(1._dp,expected_l1)).or. &
        any(abs(elements-expected_elements)>2.e-12_dp*max(1._dp,expected_elements))) &
        error stop 'Inconsistent patch inventories'
    else
      call write_reactive_eb_3d_checkpoint('',species,reactions,transport,cfg,geometry,state,temperature, &
        record_time,record_step,expected,expected_l1,expected_elements,record_dt,0._dp,record_theta, &
        valid,reason,bundle_sha256=bundle_identity, &
        chemistry_integrator=checkpoint_integrator,base_mole_fractions=base_mole_fractions,checkpoint_unit=stream, &
        element_names=element_names,element_atoms=element_atoms)
      if (.not.valid) then
        write(*,'(a)') trim(reason)
        error stop 'Cannot serialize two-level patch'
      end if
    end if
  end subroutine patch_record

  subroutine save_checkpoint()
    integer :: scratch,destination,ios,amount
    character(len=4096) :: chunk
    call check_balance()
    open(newunit=scratch,status='scratch',form='formatted',action='readwrite',iostat=ios)
    if (ios/=0) error stop 'Cannot stage two-level checkpoint'
    write(scratch,'(a)',iostat=ios) trim(checkpoint_magic())
    if (ios/=0) error stop 'Checkpoint magic write failed'
    if (regrid_count>0) call schedule_record(scratch,.false.)
    write(scratch,*,iostat=ios) refinement_ratio,fine_lower,fine_upper
    if (ios/=0) error stop 'Checkpoint topology write failed'
    if (patch_count>1) then
      write(scratch,*,iostat=ios) patch_count,patch_lower(:,2:patch_count),patch_upper(:,2:patch_count)
      if (ios/=0) error stop 'Additional checkpoint topology write failed'
    end if
    write(scratch,'(a)',iostat=ios) trim(chemistry_integrator)
    if (ios/=0) error stop 'Checkpoint integrator write failed'
    write(scratch,'(4(es26.17e3,1x),i0)',iostat=ios) time_step,time,minimum_dt,minimum_theta,step
    if (ios/=0) error stop 'Checkpoint clock write failed'
    write(scratch,'(*(es26.17e3,1x))',iostat=ios) initial,scale,total_boundary,total_wall,total_chemical
    if (ios/=0) error stop 'Checkpoint inventory write failed'
    call patch_record(scratch,.false.,coarse,uc,tc,time,step,minimum_dt,minimum_theta)
    if (patch_count>0) call patch_record(scratch,.false.,fine,uf,tf,time,step,minimum_dt,minimum_theta)
    do p=2,patch_count
      call patch_record(scratch,.false.,patches(p)%geometry,patches(p)%state,patches(p)%temperature, &
        time,step,minimum_dt,minimum_theta)
    end do
    write(scratch,'(a)',iostat=ios) 'END_TWO_LEVEL_CHECKPOINT'
    if (ios/=0) error stop 'Checkpoint terminator write failed'
    rewind(scratch,iostat=ios)
    if (ios/=0) error stop 'Checkpoint rewind failed'
    ! Only publish after both records have passed their full preflight. Like
    ! the standalone writer this is not an atomic replacement on I/O failure.
    open(newunit=destination,file=trim(checkpoint_file),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open checkpoint destination'
    do
      read(scratch,'(a)',advance='no',size=amount,iostat=ios) chunk
      if (ios==iostat_end) exit
      if (ios/=0.and.ios/=iostat_eor) error stop 'Checkpoint staged read failed'
      write(destination,'(a)',advance='no',iostat=status) chunk(:amount)
      if (status/=0) error stop 'Checkpoint publication failed'
      if (ios==iostat_eor) then
        write(destination,'(a)',iostat=status) ''
        if (status/=0) error stop 'Checkpoint newline failed'
      end if
    end do
    close(destination,iostat=ios)
    if (ios/=0) error stop 'Checkpoint destination close failed'
    close(scratch,iostat=ios)
    if (ios/=0) error stop 'Checkpoint staging close failed'
  end subroutine save_checkpoint

  subroutine restore_checkpoint()
    integer :: stream,ios,ratio,lower(3),upper(3),saved_step,component,patch_index,saved_count
    integer, allocatable :: stored_lower(:,:),stored_upper(:,:)
    type(eb_reactive_patch_3d), allocatable :: extras(:)
    real(dp) :: saved_time,saved_dt,saved_theta,saved_time_step
    real(dp) :: saved_initial(nvar),saved_scale(nvar),saved_boundary(nvar),saved_wall(nvar),saved_chemical(nvar)
    real(dp) :: values(nvar),residual
    real(dp), allocatable :: candidate_c(:,:,:,:),candidate_f(:,:,:,:),candidate_tc(:,:,:),candidate_tf(:,:,:)
    real(dp), allocatable :: averaged(:,:,:,:),averaged_t(:,:,:)
    logical :: valid
    type(eb_geometry_3d) :: geometry_c,geometry_f
    character(len=128) :: marker,integrator
    open(newunit=stream,file=trim(restart_file),status='old',action='read',iostat=ios)
    if (ios/=0) error stop 'Cannot open two-level restart'
    read(stream,'(a)',iostat=ios) marker
    if (ios/=0) error stop 'Missing checkpoint magic'
    if (marker/=checkpoint_magic()) error stop 'Invalid checkpoint magic or timestep mode'
    if (regrid_count>0) call schedule_record(stream,.true.)
    allocate(stored_lower(3,max(0,patch_count-1)),stored_upper(3,max(0,patch_count-1)))
    ratio=-1;lower=-1;upper=-1
    read(stream,*,iostat=ios) ratio,lower,upper
    if (ios/=0) error stop 'Missing checkpoint topology'
    if (ratio/=refinement_ratio.or.any(lower/=fine_lower).or.any(upper/=fine_upper)) &
      error stop 'Checkpoint topology mismatch'
    if (patch_count>1) then
      saved_count=-1;stored_lower=-1;stored_upper=-1
      read(stream,*,iostat=ios) saved_count,stored_lower,stored_upper
      if (ios/=0) error stop 'Missing additional checkpoint topology'
      if (saved_count/=patch_count.or.any(stored_lower/=patch_lower(:,2:patch_count)).or. &
        any(stored_upper/=patch_upper(:,2:patch_count))) error stop 'Additional checkpoint topology mismatch'
    end if
    read(stream,'(a)',iostat=ios) integrator
    if (ios/=0) error stop 'Missing checkpoint integrator'
    if (integrator/=chemistry_integrator) error stop 'Checkpoint integrator mismatch'
    saved_time_step=ieee_value(0._dp,ieee_quiet_nan)
    saved_time=saved_time_step;saved_dt=saved_time_step;saved_theta=saved_time_step;saved_step=-1
    saved_initial=saved_time_step;saved_scale=saved_time_step
    saved_boundary=saved_time_step;saved_wall=saved_time_step;saved_chemical=saved_time_step
    read(stream,*,iostat=ios) saved_time_step,saved_time,saved_dt,saved_theta,saved_step
    if (ios/=0) error stop 'Missing checkpoint clock'
    if (.not.all(ieee_is_finite([saved_time_step,saved_time,saved_dt,saved_theta]))) &
      error stop 'Nonfinite checkpoint clock'
    if (saved_step<1) error stop 'Invalid checkpoint step'
    if (regrid_count>0) then
      if (count(regrid_times(1:regrid_count)<=saved_time)/=regrid_stage) &
        error stop 'Checkpoint regrid stage disagrees with time'
    end if
    if (saved_time_step/=time_step) error stop 'Checkpoint time_step mismatch'
    read(stream,*,iostat=ios) saved_initial,saved_scale,saved_boundary,saved_wall,saved_chemical
    if (ios/=0) error stop 'Missing checkpoint inventories'
    if (.not.all(ieee_is_finite([saved_initial,saved_scale,saved_boundary,saved_wall,saved_chemical]))) &
      error stop 'Nonfinite checkpoint inventories'
    if (any(saved_scale<0._dp).or.any(abs(saved_initial-initial)>2.e-12_dp*max(1._dp,scale)).or. &
      any(abs(saved_scale-scale)>2.e-12_dp*max(1._dp,scale))) error stop 'Checkpoint initial inventory mismatch'
    candidate_c=uc;candidate_tc=tc
    geometry_c=coarse
    call patch_record(stream,.true.,geometry_c,candidate_c,candidate_tc,saved_time,saved_step,saved_dt,saved_theta)
    if (patch_count>0) then
      candidate_f=uf;candidate_tf=tf;geometry_f=fine
      call patch_record(stream,.true.,geometry_f,candidate_f,candidate_tf,saved_time,saved_step,saved_dt,saved_theta)
    end if
    extras=patches
    do patch_index=2,patch_count
      call patch_record(stream,.true.,extras(patch_index)%geometry,extras(patch_index)%state, &
        extras(patch_index)%temperature,saved_time,saved_step,saved_dt,saved_theta)
    end do
    read(stream,'(a)',iostat=ios) marker
    if (ios/=0) error stop 'Missing composite checkpoint terminator'
    if (marker/='END_TWO_LEVEL_CHECKPOINT') error stop 'Invalid composite checkpoint terminator'
    read(stream,'(a)',iostat=ios) marker
    if (ios/=iostat_end) error stop 'Trailing composite checkpoint data'
    close(stream,iostat=ios)
    if (ios/=0) error stop 'Checkpoint read close failed'
    do component=1,nvar
      values(component)=coarse%dx*coarse%dy*coarse%dz* &
        sum(candidate_c(component,:,:,:)*coarse%volume_fraction,mask=.not.refined)
      if (patch_count>0) values(component)=values(component)+ &
        fine%dx*fine%dy*fine%dz*sum(candidate_f(component,:,:,:)*fine%volume_fraction)
      do patch_index=2,patch_count
        values(component)=values(component)+patches(patch_index)%geometry%dx*patches(patch_index)%geometry%dy* &
          patches(patch_index)%geometry%dz*sum(extras(patch_index)%state(component,:,:,:)* &
            patches(patch_index)%geometry%volume_fraction)
      end do
    end do
    residual=maxval(abs(values-saved_initial-saved_boundary-saved_wall-saved_chemical)/ &
      max(1._dp,saved_scale,abs(values),abs(saved_boundary)+abs(saved_wall)+abs(saved_chemical)))
    if (.not.ieee_is_finite(residual)) error stop 'Nonfinite checkpoint balance'
    if (residual>5.e-10_dp) error stop 'Checkpoint composite balance mismatch'
#ifdef PELEF_SELECTED_EB_TWO_LEVEL
    if (.not.chemical_elements_valid(saved_initial,saved_chemical,saved_step)) error stop 'Checkpoint chemical element imbalance'
#endif
    if (patch_count>0) then
    averaged=candidate_c(:,fine_lower(1):fine_upper(1),fine_lower(2):fine_upper(2),fine_lower(3):fine_upper(3))
    averaged_t=candidate_tc(fine_lower(1):fine_upper(1),fine_lower(2):fine_upper(2),fine_lower(3):fine_upper(3))
    call restrict_reactive_eb_state_3d(species,fine,refinement_ratio,candidate_f,candidate_tf,averaged,averaged_t,valid)
    if (.not.valid) error stop 'Checkpoint restriction failed'
    if (any(abs(averaged-candidate_c(:,fine_lower(1):fine_upper(1),fine_lower(2):fine_upper(2), &
      fine_lower(3):fine_upper(3)))>2.e-12_dp*max(1._dp,abs(averaged)))) &
      error stop 'Checkpoint parent averages mismatch'
    end if
    do patch_index=2,patch_count
      lower=patch_lower(:,patch_index);upper=patch_upper(:,patch_index)
      averaged=candidate_c(:,lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))
      averaged_t=candidate_tc(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))
      call restrict_reactive_eb_state_3d(species,patches(patch_index)%geometry,refinement_ratio, &
        extras(patch_index)%state,extras(patch_index)%temperature,averaged,averaged_t,valid)
      if (.not.valid) error stop 'Additional checkpoint restriction failed'
      if (any(abs(averaged-candidate_c(:,lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)))> &
        2.e-12_dp*max(1._dp,abs(averaged)))) error stop 'Additional checkpoint parent mismatch'
    end do
    ! Publish both patches and the outer diagnostics only after the complete
    ! file, physics, geometry, EOS and budget have passed validation.
    uc=candidate_c;tc=candidate_tc
    if (patch_count>0) then
      uf=candidate_f;tf=candidate_tf
    end if
    do patch_index=2,patch_count
      patches(patch_index)%state=extras(patch_index)%state
      patches(patch_index)%temperature=extras(patch_index)%temperature
    end do
    initial=saved_initial;scale=saved_scale
    total_boundary=saved_boundary;total_wall=saved_wall;total_chemical=saved_chemical
    time=saved_time;step=saved_step;minimum_dt=saved_dt;minimum_theta=saved_theta
  end subroutine restore_checkpoint

  subroutine inventory(values)
    real(dp), intent(out) :: values(:)
    integer :: component,patch_index
    do component=1,size(values)
      values(component)=coarse%dx*coarse%dy*coarse%dz*sum(uc(component,:,:,:)*coarse%volume_fraction,mask=.not.refined)
      do patch_index=1,patch_count
        values(component)=values(component)+patches(patch_index)%geometry%dx*patches(patch_index)%geometry%dy* &
          patches(patch_index)%geometry%dz*sum(patches(patch_index)%state(component,:,:,:)* &
          patches(patch_index)%geometry%volume_fraction)
      end do
    end do
  end subroutine inventory

  subroutine write_output()
    integer :: stream,ios,component,level,i,j,k,patch_index
    open(newunit=stream,file=trim(output_file),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open composite CSV'
    if (patch_count>1) then
      write(stream,'(a)',advance='no',iostat=ios) 'patch,'
      if (ios/=0) error stop 'CSV patch header failure'
    end if
    write(stream,'(a)',advance='no',iostat=ios) 'level,i,j,k,x,y,z,kappa,temperature,rho,rhou,rhov,rhow,rhoE'
    if (ios/=0) error stop 'CSV header failure'
    do component=1,size(species)
      write(stream,'(a)',advance='no',iostat=ios) ',rhoY_'//trim(species(component)%name)
      if (ios/=0) error stop 'CSV species header failure'
    end do
    write(stream,*,iostat=ios)
    if (ios/=0) error stop 'CSV header terminator failure'
    do level=0,1
      if (level==0) then
        do k=1,coarse%nz
          do j=1,coarse%ny
            do i=1,coarse%nx
              if (refined(i,j,k).or.coarse%volume_fraction(i,j,k)<=0._dp) cycle
              call row(stream,level,i,j,k,coarse,uc(:,i,j,k),tc(i,j,k),0)
            end do
          end do
        end do
      else
        do patch_index=1,patch_count
          do k=1,patches(patch_index)%geometry%nz
            do j=1,patches(patch_index)%geometry%ny
              do i=1,patches(patch_index)%geometry%nx
                if (patches(patch_index)%geometry%volume_fraction(i,j,k)<=0._dp) cycle
                call row(stream,level,i,j,k,patches(patch_index)%geometry,patches(patch_index)%state(:,i,j,k), &
                  patches(patch_index)%temperature(i,j,k),patch_index)
              end do
            end do
          end do
        end do
      end if
    end do
    close(stream,iostat=ios)
    if (ios/=0) error stop 'CSV close failure'
    open(newunit=stream,file=trim(budget_path),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open budget CSV'
    write(stream,'(a)',iostat=ios) 'component,initial,final,boundary,wall,chemistry,residual'
    if (ios/=0) error stop 'Budget header failure'
    do component=1,nvar
      write(stream,'(i0,6(",",es26.17e3))',iostat=ios) component,initial(component),current(component), &
        total_boundary(component),total_wall(component),total_chemical(component), &
        current(component)-initial(component)-total_boundary(component)-total_wall(component)-total_chemical(component)
      if (ios/=0) error stop 'Budget row failure'
    end do
    close(stream,iostat=ios)
    if (ios/=0) error stop 'Budget close failure'
  end subroutine write_output

    subroutine row(stream,level,i,j,k,geometry,state,temperature,patch_index)
      integer, intent(in) :: stream,level,i,j,k,patch_index
      integer :: ios
      type(eb_geometry_3d), intent(in) :: geometry
      real(dp), intent(in) :: state(:),temperature
      if (patch_count>1) then
        write(stream,'(i0,a)',advance='no',iostat=ios) patch_index,','
        if (ios/=0) error stop 'CSV patch index failure'
      end if
      write(stream,'(i0,3(",",i0),*(:,",",es26.17e3))',iostat=ios) level,i,j,k, &
        geometry%x_lower+(real(i,dp)-.5_dp+geometry%cell_centroid_x(i,j,k))*geometry%dx, &
        geometry%y_lower+(real(j,dp)-.5_dp+geometry%cell_centroid_y(i,j,k))*geometry%dy, &
        geometry%z_lower+(real(k,dp)-.5_dp+geometry%cell_centroid_z(i,j,k))*geometry%dz, &
        geometry%volume_fraction(i,j,k),temperature,state
      if (ios/=0) error stop 'Composite CSV row failure'
    end subroutine row
end program pelef_reactive_eb_two_level_3d
