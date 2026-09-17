program pelef_reactive_eb_two_level_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  use, intrinsic :: iso_fortran_env, only: iostat_end, iostat_eor
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use reactive_eb_3d_driver_mod, only: reactive_eb_integrals_3d,reactive_eb_element_integrals_3d
  use reactive_eb_3d_checkpoint_mod, only: write_reactive_eb_3d_checkpoint,read_reactive_eb_3d_checkpoint
  use eb_reactive_transfer_3d_mod, only: restrict_reactive_eb_state_3d
  use reactive_eb_cfl_3d_mod, only: reactive_eb_signal_speeds_3d,compute_reactive_eb_euler_cfl_timestep_3d
  use eb_reactive_transport_3d_mod, only: reactive_eb_transport_patch_timestep_3d
  use eb_reactive_two_level_3d_mod, only: advance_reactive_eb_two_level_full_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_eb_3d_config) :: config
  type(eb_geometry_3d) :: full,coarse,fine
  integer :: nx=8,ny=8,nz=8,refinement_ratio=2,fine_lower(3)=[2,2,2],fine_upper(3)=[3,3,3],maximum_steps=100
  character(len=16) :: geometry_type='plane',plane_axis='z',extrusion_axis='z',riemann_solver='rusanov'
  character(len=16) :: chemistry_integrator='explicit'
  character(len=512) :: level_set_file='',section_file='',output_file='eb_two_level.csv'
  real(dp) :: domain_lower(3)=0._dp,domain_upper(3)=1._dp,plane_position=.31_dp,time_step=1.e-8_dp,final_time=3.e-8_dp
  character(len=512) :: checkpoint_file='',restart_file=''
  integer :: checkpoint_interval_steps=0,maximum_step_retries=12
  logical :: stop_after_checkpoint=.false.
  logical :: automatic_timestep=.false.,subcycling=.false.
  real(dp) :: hydro_cfl=.5_dp,transport_cfl=.35_dp
  real(dp) :: rtol=2.e-7_dp,atol=1.e-12_dp
  logical :: chemistry_enabled=.true.,transport_enabled=.true.,viscosity_enabled=.true.
  logical :: thermal_conduction_enabled=.true.,species_diffusion_enabled=.true.,barodiffusion_enabled=.true.
  namelist /eb_two_level/ nx,ny,nz,refinement_ratio,fine_lower,fine_upper,maximum_steps,geometry_type,plane_axis, &
    extrusion_axis,level_set_file,section_file,output_file,domain_lower,domain_upper,plane_position,time_step,final_time, &
    rtol,atol,chemistry_enabled,transport_enabled,viscosity_enabled,thermal_conduction_enabled, &
    species_diffusion_enabled,barodiffusion_enabled,chemistry_integrator,riemann_solver, &
    checkpoint_file,restart_file,checkpoint_interval_steps,stop_after_checkpoint,maximum_step_retries, &
    automatic_timestep,hydro_cfl,transport_cfl,subcycling
  real(dp), allocatable :: uc(:,:,:,:),tc(:,:,:),uf(:,:,:,:),tf(:,:,:),initial(:),current(:),scale(:), &
    boundary(:),wall(:),chemical(:),total_boundary(:),total_wall(:),total_chemical(:)
  logical, allocatable :: refined(:,:,:)
  real(dp) :: time,dt,theta,error,minimum_theta,minimum_dt
  integer :: unit,status,close_status,nvar,step,dims(3),nf(3),c,retry,rejected_steps
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
  if (any(fine_lower<1).or.any(fine_upper<fine_lower).or.any(fine_upper>dims)) error stop 'Invalid fine patch box'
  if (.not.all(ieee_is_finite([domain_lower,domain_upper,time_step,final_time,rtol,atol]))) &
    error stop 'Nonfinite two-level controls'
  if (any(domain_upper<=domain_lower).or.min(time_step,final_time,rtol,atol)<=0._dp) &
    error stop 'Invalid two-level physical controls'
  if (barodiffusion_enabled.and..not.species_diffusion_enabled.and.transport_enabled) &
    error stop 'Barodiffusion requires species diffusion'
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
  config%redistribution='none'
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
  call extract_eb_geometry_3d(full,(fine_lower-1)*refinement_ratio+1,fine_upper*refinement_ratio,fine,ok)
  if (.not.ok) error stop 'Cannot extract fine EB patch'
  call load_h2o2_elementary_thermo(species,ok)
  if (.not.ok) error stop 'Thermodynamic load failed'
  call load_h2o2_elementary_mechanism(reactions,ok)
  if (.not.ok) error stop 'Mechanism load failed'
  call load_h2o2_elementary_transport(transport,ok)
  if (.not.ok) error stop 'Transport load failed'
  nvar=reactive_nvar(size(species));nf=[fine%nx,fine%ny,fine%nz]
  allocate(uc(nvar,dims(1),dims(2),dims(3)),tc(dims(1),dims(2),dims(3)), &
    uf(nvar,nf(1),nf(2),nf(3)),tf(nf(1),nf(2),nf(3)),refined(dims(1),dims(2),dims(3)))
  config%nx=dims(1);config%ny=dims(2);config%nz=dims(3)
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,uc,tc,ok)
  if (.not.ok) error stop 'Coarse initialization failed'
  config%nx=nf(1);config%ny=nf(2);config%nz=nf(3)
  call initialize_reactive_eb_density_sheet_3d(species,config,fine,uf,tf,ok)
  if (.not.ok) error stop 'Fine initialization failed'
  refined=.false.;refined(fine_lower(1):fine_upper(1),fine_lower(2):fine_upper(2),fine_lower(3):fine_upper(3))=.true.
  allocate(initial(nvar),current(nvar),scale(nvar),boundary(nvar),wall(nvar),chemical(nvar), &
    total_boundary(nvar),total_wall(nvar),total_chemical(nvar))
  call inventory(initial)
  do c=1,nvar
    scale(c)=max(1.e-30_dp,coarse%dx*coarse%dy*coarse%dz*sum(abs(uc(c,:,:,:))*coarse%volume_fraction, &
      mask=.not.refined)+fine%dx*fine%dy*fine%dz*sum(abs(uf(c,:,:,:))*fine%volume_fraction))
  end do
  total_boundary=0._dp;total_wall=0._dp;total_chemical=0._dp
  time=0._dp;step=0;minimum_theta=1._dp;minimum_dt=0._dp;rejected_steps=0
  if (len_trim(restart_file)>0) call restore_checkpoint()
  do while(time<final_time)
    if (step>=maximum_steps) error stop 'Maximum steps reached before final time'
    dt=min(time_step,final_time-time)
    if (automatic_timestep) call limit_timestep(dt)
    do retry=0,maximum_step_retries
      if (time+dt<=time) error stop 'Time step cannot advance time'
      ! The full update is transactional: a rejected trial leaves both grids,
      ! per-step ledgers and the limiter untouched. Count only accepted steps.
      call advance_reactive_eb_two_level_full_3d(species,reactions,transport,coarse,fine,fine_lower,refinement_ratio, &
        riemann_solver,dt,chemistry_enabled,rtol,atol,transport_enabled,viscosity_enabled,thermal_conduction_enabled, &
        species_diffusion_enabled,barodiffusion_enabled,uc,tc,uf,tf,boundary,wall,chemical,theta,ok,chemistry_integrator,subcycling)
      if (ok) exit
      rejected_steps=rejected_steps+1
      if (retry==maximum_step_retries) error stop 'Two-level EB step retry limit exhausted'
      dt=0.5_dp*dt
    end do
    total_boundary=total_boundary+boundary;total_wall=total_wall+wall;total_chemical=total_chemical+chemical
    time=time+dt;step=step+1;minimum_theta=min(minimum_theta,theta)
    if (step==1) minimum_dt=dt
    minimum_dt=min(minimum_dt,dt)
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
  write(*,'(a,es24.16)') 'Final time: ',time
  write(*,'(a,es24.16)') 'Composite inventory error: ',error
  write(*,'(a,es24.16)') 'Minimum transport theta: ',minimum_theta
  write(*,'(a,es24.16)') 'Minimum accepted dt: ',minimum_dt
contains
  function checkpoint_magic() result(magic)
    character(len=32) :: magic
    ! V1/V2 retain common-step fixed/automatic compatibility. V3/V4 select
    ! subcycling with the same layout; patch configs bind the CFL controls.
    magic='PELEF_EB_TWO_LEVEL_CHECKPOINT_V1'
    if (automatic_timestep) magic='PELEF_EB_TWO_LEVEL_CHECKPOINT_V2'
    if (subcycling) then
      magic='PELEF_EB_TWO_LEVEL_CHECKPOINT_V3'
      if (automatic_timestep) magic='PELEF_EB_TWO_LEVEL_CHECKPOINT_V4'
    end if
  end function checkpoint_magic

  subroutine limit_timestep(interval)
    real(dp), intent(inout) :: interval
    real(dp), allocatable :: parent(:,:,:,:),parent_t(:,:,:)
    real(dp) :: coarse_speeds(4),fine_speeds(4),donors(4),coarse_dt,fine_dt,coarse_d,fine_d,donor_d,unused_d,factor
    logical :: valid,fine_active
    ! Match the parent averaging used at the start of the coupled update;
    ! estimation must not publish state changes of its own.
    parent=uc;parent_t=tc
    call restrict_reactive_eb_state_3d(species,fine,refinement_ratio,uf,tf, &
      parent(:,fine_lower(1):fine_upper(1),fine_lower(2):fine_upper(2),fine_lower(3):fine_upper(3)), &
      parent_t(fine_lower(1):fine_upper(1),fine_lower(2):fine_upper(2),fine_lower(3):fine_upper(3)),valid)
    if (.not.valid) error stop 'Cannot average parents for timestep selection'
    call reactive_eb_signal_speeds_3d(species,parent,parent_t,coarse,coarse_speeds,valid)
    if (.not.valid) error stop 'Invalid coarse timestep state'
    fine_active=any(fine%volume_fraction>0._dp)
    fine_speeds=0._dp
    if (fine_active) then
      call reactive_eb_signal_speeds_3d(species,uf,tf,fine,fine_speeds,valid)
      if (.not.valid) error stop 'Invalid fine timestep state'
    end if
    donors=max(coarse_speeds,fine_speeds)
    call compute_reactive_eb_euler_cfl_timestep_3d(species,parent,parent_t,coarse,hydro_cfl,coarse_dt,valid,donors)
    if (.not.valid) error stop 'Cannot select coarse hydro interval'
    fine_dt=huge(1._dp)
    if (fine_active) then
      call compute_reactive_eb_euler_cfl_timestep_3d(species,uf,tf,fine,hydro_cfl,fine_dt,valid,donors)
      if (.not.valid) error stop 'Cannot select fine hydro interval'
    end if
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
    call reactive_eb_transport_patch_timestep_3d(species,transport,parent,parent_t,coarse,transport_cfl, &
      viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,coarse_dt,unused_d,valid,donor_d)
    if (.not.valid) error stop 'Cannot select coarse diffusion interval'
    fine_dt=huge(1._dp)
    if (fine_active) then
      call reactive_eb_transport_patch_timestep_3d(species,transport,uf,tf,fine,transport_cfl, &
        viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,fine_dt,unused_d,valid,donor_d)
      if (.not.valid) error stop 'Cannot select fine diffusion interval'
    end if
    ! Bound the full interval, conservatively also bounding each transport half.
    interval=min(interval,coarse_dt)
    factor=1._dp
    if (subcycling) factor=real(refinement_ratio,dp)**2
    if (fine_dt<interval/factor) interval=fine_dt*factor
  end subroutine limit_timestep

  subroutine check_balance()
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
    real(dp) :: integrals(nvar),l1(nvar),elements(3),got_time,got_dt,got_theta,diffusivity
    real(dp) :: expected(nvar),expected_l1(nvar),expected_elements(3)
    integer :: got_step
    logical :: valid
    character(len=512) :: reason
    cfg=patch_config(geometry)
    if (reading) then
      call read_reactive_eb_3d_checkpoint('',species,reactions,transport,cfg,geometry,state,temperature, &
        got_time,got_step,integrals,l1,elements,got_dt,diffusivity,got_theta,valid,reason,checkpoint_unit=stream)
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
    call reactive_eb_element_integrals_3d(species,state,geometry,expected_elements,valid)
    if (.not.valid) error stop 'Cannot integrate checkpoint elements'
    if (reading) then
      if (any(abs(integrals-expected)>2.e-12_dp*max(1._dp,abs(expected))).or. &
        any(abs(l1-expected_l1)>2.e-12_dp*max(1._dp,expected_l1)).or. &
        any(abs(elements-expected_elements)>2.e-12_dp*max(1._dp,expected_elements))) &
        error stop 'Inconsistent patch inventories'
    else
      call write_reactive_eb_3d_checkpoint('',species,reactions,transport,cfg,geometry,state,temperature, &
        record_time,record_step,expected,expected_l1,expected_elements,record_dt,0._dp,record_theta, &
        valid,reason,checkpoint_unit=stream)
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
    write(scratch,'(a)',iostat=ios) checkpoint_magic()
    if (ios/=0) error stop 'Checkpoint magic write failed'
    write(scratch,*,iostat=ios) refinement_ratio,fine_lower,fine_upper
    if (ios/=0) error stop 'Checkpoint topology write failed'
    write(scratch,'(a)',iostat=ios) trim(chemistry_integrator)
    if (ios/=0) error stop 'Checkpoint integrator write failed'
    write(scratch,'(4(es26.17e3,1x),i0)',iostat=ios) time_step,time,minimum_dt,minimum_theta,step
    if (ios/=0) error stop 'Checkpoint clock write failed'
    write(scratch,'(*(es26.17e3,1x))',iostat=ios) initial,scale,total_boundary,total_wall,total_chemical
    if (ios/=0) error stop 'Checkpoint inventory write failed'
    call patch_record(scratch,.false.,coarse,uc,tc,time,step,minimum_dt,minimum_theta)
    call patch_record(scratch,.false.,fine,uf,tf,time,step,minimum_dt,minimum_theta)
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
    integer :: stream,ios,ratio,lower(3),upper(3),saved_step,component
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
    ratio=-1;lower=-1;upper=-1
    read(stream,*,iostat=ios) ratio,lower,upper
    if (ios/=0) error stop 'Missing checkpoint topology'
    if (ratio/=refinement_ratio.or.any(lower/=fine_lower).or.any(upper/=fine_upper)) &
      error stop 'Checkpoint topology mismatch'
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
    if (saved_time_step/=time_step) error stop 'Checkpoint time_step mismatch'
    read(stream,*,iostat=ios) saved_initial,saved_scale,saved_boundary,saved_wall,saved_chemical
    if (ios/=0) error stop 'Missing checkpoint inventories'
    if (.not.all(ieee_is_finite([saved_initial,saved_scale,saved_boundary,saved_wall,saved_chemical]))) &
      error stop 'Nonfinite checkpoint inventories'
    if (any(saved_scale<0._dp).or.any(abs(saved_initial-initial)>2.e-12_dp*max(1._dp,scale)).or. &
      any(abs(saved_scale-scale)>2.e-12_dp*max(1._dp,scale))) error stop 'Checkpoint initial inventory mismatch'
    candidate_c=uc;candidate_tc=tc;candidate_f=uf;candidate_tf=tf
    geometry_c=coarse;geometry_f=fine
    call patch_record(stream,.true.,geometry_c,candidate_c,candidate_tc,saved_time,saved_step,saved_dt,saved_theta)
    call patch_record(stream,.true.,geometry_f,candidate_f,candidate_tf,saved_time,saved_step,saved_dt,saved_theta)
    read(stream,'(a)',iostat=ios) marker
    if (ios/=0) error stop 'Missing composite checkpoint terminator'
    if (marker/='END_TWO_LEVEL_CHECKPOINT') error stop 'Invalid composite checkpoint terminator'
    read(stream,'(a)',iostat=ios) marker
    if (ios/=iostat_end) error stop 'Trailing composite checkpoint data'
    close(stream,iostat=ios)
    if (ios/=0) error stop 'Checkpoint read close failed'
    do component=1,nvar
      values(component)=coarse%dx*coarse%dy*coarse%dz* &
        sum(candidate_c(component,:,:,:)*coarse%volume_fraction,mask=.not.refined) &
        +fine%dx*fine%dy*fine%dz*sum(candidate_f(component,:,:,:)*fine%volume_fraction)
    end do
    residual=maxval(abs(values-saved_initial-saved_boundary-saved_wall-saved_chemical)/ &
      max(1._dp,saved_scale,abs(values),abs(saved_boundary)+abs(saved_wall)+abs(saved_chemical)))
    if (.not.ieee_is_finite(residual)) error stop 'Nonfinite checkpoint balance'
    if (residual>5.e-10_dp) error stop 'Checkpoint composite balance mismatch'
    averaged=candidate_c(:,fine_lower(1):fine_upper(1),fine_lower(2):fine_upper(2),fine_lower(3):fine_upper(3))
    averaged_t=candidate_tc(fine_lower(1):fine_upper(1),fine_lower(2):fine_upper(2),fine_lower(3):fine_upper(3))
    call restrict_reactive_eb_state_3d(species,fine,refinement_ratio,candidate_f,candidate_tf,averaged,averaged_t,valid)
    if (.not.valid) error stop 'Checkpoint restriction failed'
    if (any(abs(averaged-candidate_c(:,fine_lower(1):fine_upper(1),fine_lower(2):fine_upper(2), &
      fine_lower(3):fine_upper(3)))>2.e-12_dp*max(1._dp,abs(averaged)))) &
      error stop 'Checkpoint parent averages mismatch'
    ! Publish both patches and the outer diagnostics only after the complete
    ! file, physics, geometry, EOS and budget have passed validation.
    uc=candidate_c;tc=candidate_tc;uf=candidate_f;tf=candidate_tf
    initial=saved_initial;scale=saved_scale
    total_boundary=saved_boundary;total_wall=saved_wall;total_chemical=saved_chemical
    time=saved_time;step=saved_step;minimum_dt=saved_dt;minimum_theta=saved_theta
  end subroutine restore_checkpoint

  subroutine inventory(values)
    real(dp), intent(out) :: values(:)
    integer :: component
    do component=1,size(values)
      values(component)=coarse%dx*coarse%dy*coarse%dz*sum(uc(component,:,:,:)*coarse%volume_fraction,mask=.not.refined) &
        +fine%dx*fine%dy*fine%dz*sum(uf(component,:,:,:)*fine%volume_fraction)
    end do
  end subroutine inventory

  subroutine write_output()
    integer :: stream,ios,component,level,i,j,k
    open(newunit=stream,file=trim(output_file),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open composite CSV'
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
              call row(stream,level,i,j,k,coarse,uc(:,i,j,k),tc(i,j,k))
            end do
          end do
        end do
      else
        do k=1,fine%nz
          do j=1,fine%ny
            do i=1,fine%nx
              if (fine%volume_fraction(i,j,k)<=0._dp) cycle
              call row(stream,level,i,j,k,fine,uf(:,i,j,k),tf(i,j,k))
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

    subroutine row(stream,level,i,j,k,geometry,state,temperature)
      integer, intent(in) :: stream,level,i,j,k
      integer :: ios
      type(eb_geometry_3d), intent(in) :: geometry
      real(dp), intent(in) :: state(:),temperature
      write(stream,'(i0,3(",",i0),*(:,",",es26.17e3))',iostat=ios) level,i,j,k, &
        geometry%x_lower+(real(i,dp)-.5_dp+geometry%cell_centroid_x(i,j,k))*geometry%dx, &
        geometry%y_lower+(real(j,dp)-.5_dp+geometry%cell_centroid_y(i,j,k))*geometry%dy, &
        geometry%z_lower+(real(k,dp)-.5_dp+geometry%cell_centroid_z(i,j,k))*geometry%dz, &
        geometry%volume_fraction(i,j,k),temperature,state
      if (ios/=0) error stop 'Composite CSV row failure'
    end subroutine row
end program pelef_reactive_eb_two_level_3d
