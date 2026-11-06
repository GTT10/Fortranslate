! Matched-state probe for the two-level planar EB reacting TG control.
! This exercises the production tree integrator with imported SI cell states.
! LiDryer CVODE chemistry and pinned Simple molecular transport; planar EB.
program eb_tree_reacting_probe
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,configure_nasa7_eos_policy,universal_gas_constant
  use lidryer_mechanism_mod, only: load_lidryer_thermo,load_lidryer_mechanism,lidryer_nspecies,load_lidryer_transport
  use sundials_forced_reactor_mod, only: advance_forced_reactor_cvode
  use dual_energy_mod, only: dual_energy_policy,reset_nasa7_dual_energy_cell
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species,load_gas_transport_data,set_simple_transport_fits
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    build_eb_tree_extended_geometry_3d,initialize_eb_tree_reaction_source_3d, &
    advance_reactive_eb_tree_sdc_3d,advance_reactive_eb_tree_sdc_subcycle_3d
  use eb_reactive_transfer_3d_mod, only: plan_eb_linear_halo_3d,interpolate_eb_linear_halo_3d, &
    reconstruct_eb_pressure_state_3d
  use mixture_thermo_mod, only: mixture_pressure
  use reactive_1d_mod, only: reactive_nprim,reactive_conserved_to_primitive,reactive_primitive_to_conserved
  use eb_mol_hydro_3d_mod, only: compute_mol_eb_hydro_3d,eb_mol_diagnostics_3d
  use eb_reactive_patch_tree_3d_mod, only: synchronize_reactive_eb_tree_flux_3d, &
    synchronize_reactive_eb_tree_scalar_flux_3d,integrate_reactive_eb_tree_3d,average_down_reactive_eb_tree_3d
  use eb_reactive_patch_tree_3d_mod, only: compute_mol_eb_tree_hydro_3d
  implicit none
  integer, parameter :: nv=5+lidryer_nspecies
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: patch(1)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,temperature,seen,held,history
  type(dual_energy_policy) :: policy
  type(eb_reactive_tree_3d) :: eb
  type(eb_geometry_3d) :: finest,padded
  type(eb_geometry_3d) :: physical(2)
  type(eb_geometry_3d), allocatable :: extended(:)
  type(eb_integrated_flux_3d), allocatable :: held_flux(:)
  real(dp) :: plane,wall(nv),chemical(nv),transport_values(5,9),fits(4,63)
  integer :: geometries(9),transport_on
  character(len=24) :: transport_names(9)
  real(dp) :: spacing(3),dt,budget(nv),source(nv),correction,theta,values(2*nv-2),force(3),energy_correction
  real(dp) :: covered_energy_gap,covered_energy_scale
  real(dp) :: reaction_rtol=1e-10_dp,reaction_atol=1e-14_dp
  real(dp) :: cvode_atol_min=huge(1._dp),cvode_atol_max=0
  integer, parameter :: initial_cache_capacity=128
  integer(int64) :: initial_keys(2*lidryer_nspecies+6,initial_cache_capacity)
  real(dp) :: initial_results(lidryer_nspecies+1,initial_cache_capacity)
  integer :: initial_limits(initial_cache_capacity),initial_steps(initial_cache_capacity)
  integer :: initial_entries=0,initial_hits=0,initial_solves=0
  integer :: unit,ios,steps,iterations,subcycle,node,i,j,k,records,row,n(3),step,c,root_size,half_size
  type(amr_tree_state_3d) :: evolution_rate,evolution_held_rate
  type(eb_integrated_flux_3d) :: evolution_flux(2),evolution_held_flux(2)
  logical :: audit_evolving=.false.,audit_reset=.false.,smooth_contact=.false.,consistent_contact=.false.,smooth_cut=.false.
  integer :: stage_unit,stage_counter=0
  integer :: sync_trace_unit,sync_trace_event=0
  integer :: mol_trace_unit,mol_trace_event=0,mol_trace_step=0
  logical :: audit_stage_flux=.false.
  procedure(trace_mol_stage), pointer :: mol_observer=>null()
  procedure(trace_mol_sync), pointer :: sync_observer=>null()
  real(dp) :: reset_energy_signed=0,reset_energy_absolute=0
  logical :: paired_reacting=.false.,velocity_reacting=.false.,momentum_reacting=.false.,reacting_option=.false.
  logical :: audit_shear=.false.,audit_shear_initial=.false.,audit_consistent=.false.,audit_connected=.false.
  logical :: audit_native=.false.,audit_initial=.false.,audit_separate=.false.,audit_pair=.false.,audit_mixing=.false.
  character(len=32) :: option
  logical :: ok
  logical :: density_absolute_tolerance=.false.
  integer :: eos_convergence=0
  character(len=1024) :: input,prefix,path
  character(len=32) :: magic,faces(6)
  if(command_argument_count()<2.or.command_argument_count()>7) &
    error stop 'Usage: probe input.dat output-prefix [--audit-initial [--separate-momentum]]'
  if(command_argument_count()>=3) then
    call get_command_argument(3,option)
    paired_reacting=trim(option)=='--paired-reacting'
    velocity_reacting=trim(option)=='--velocity-reacting'.or.paired_reacting
    momentum_reacting=trim(option)=='--momentum-reacting'.or.velocity_reacting
    reacting_option=momentum_reacting.or.trim(option)=='--default-reacting-trace'
    if(reacting_option.and.command_argument_count()/=3) error stop 'Reacting policy accepts no audit options'
    if(.not.reacting_option.and.trim(option)/='--audit-initial'.and.trim(option)/='--audit-mixing'.and. &
      trim(option)/='--audit-native'.and.trim(option)/='--audit-shear'.and. &
      trim(option)/='--audit-shear-stage-flux'.and.trim(option)/='--audit-shear-connected'.and. &
      trim(option)/='--audit-shear-initial') error stop 'Unknown probe option'
    audit_shear_initial=trim(option)=='--audit-shear-initial'
    audit_stage_flux=trim(option)=='--audit-shear-stage-flux'.or.trim(option)=='--default-reacting-trace'
    if(audit_stage_flux) then
      mol_observer=>trace_mol_stage;sync_observer=>trace_mol_sync
    endif
    audit_connected=trim(option)=='--audit-shear-connected'.or.audit_stage_flux
    audit_shear=trim(option)=='--audit-shear'.or.audit_shear_initial.or.audit_connected
    if(audit_shear.and.command_argument_count()/=3) error stop 'Shear audit accepts no other options'
    audit_native=trim(option)=='--audit-native'.or.audit_shear
    audit_initial=.not.reacting_option;audit_mixing=trim(option)=='--audit-mixing'
    if(audit_mixing.and.command_argument_count()>4) error stop 'Mixing audit accepts one profile option'
  endif
  if(command_argument_count()>=4) then
    call get_command_argument(4,option)
    if(audit_mixing) then
      if(trim(option)/='--smooth-cut-contact') error stop 'Unknown mixing profile'
      smooth_contact=.true.;smooth_cut=.true.
    else
      if(trim(option)/='--separate-momentum') error stop 'Unknown audit policy'
      audit_separate=.true.
    endif
  endif
  if(command_argument_count()>=5) then
    call get_command_argument(5,option)
    if(trim(option)/='--pressure-pair'.and.trim(option)/='--pressure-evolve'.and. &
      trim(option)/='--pressure-evolve-reset') error stop 'Unknown pressure pair option'
    audit_reset=trim(option)=='--pressure-evolve-reset'
    audit_evolving=trim(option)=='--pressure-evolve'.or.audit_reset
    audit_pair=.true.
  endif
  if(command_argument_count()>=6) then
    call get_command_argument(6,option)
    if(.not.audit_evolving) error stop 'Contact control requires evolution'
    select case(trim(option))
    case('--smooth-contact');smooth_contact=.true.
    case('--smooth-cut-contact');smooth_contact=.true.;smooth_cut=.true.
    case('--imported-state')
    case default;error stop 'Invalid contact initial state option'
    end select
  endif
  if(command_argument_count()==7) then
    call get_command_argument(7,option)
    if(trim(option)/='--consistent-velocity') error stop 'Unknown contact policy'
    consistent_contact=.true.
  endif
  call get_command_argument(1,input);call get_command_argument(2,prefix)
  open(newunit=unit,file=trim(input),status='old',action='read',iostat=ios)
  if(ios/=0) error stop 'Cannot open probe input'
  read(unit,*,iostat=ios) magic
  if(ios/=0) error stop 'Invalid probe input magic'
  root_size=8
  select case(magic)
  case('PELEF_REACTING_EB_REFERENCE_V1')
    ! Preserve the frozen 8-cubed input layout and numerical path.
  case('PELEF_REACTING_EB_REFERENCE_V2','PELEF_REACTING_EB_REFERENCE_V3','PELEF_REACTING_EB_REFERENCE_V4')
    read(unit,*,iostat=ios) root_size
    if(ios/=0) error stop 'Invalid root size row'
    if(.not.any(root_size==[8,16,32])) error stop 'Root size must be 8, 16 or 32'
    if(magic=='PELEF_REACTING_EB_REFERENCE_V3'.or.magic=='PELEF_REACTING_EB_REFERENCE_V4') then
      density_absolute_tolerance=.true.
      ! Explicit density tolerances in SI; V1/V2 retain their historical values.
      read(unit,*,iostat=ios) reaction_rtol,reaction_atol
      if(ios/=0) error stop 'Invalid reaction tolerances'
      if(.not.all(ieee_is_finite([reaction_rtol,reaction_atol]))) error stop 'Nonfinite reaction tolerances'
      if(reaction_rtol<=0.or.reaction_rtol>1e-10_dp.or.reaction_atol<=0.or.reaction_atol>1e-15_dp) &
        error stop 'Reaction tolerances require rtol <= 1e-10 and density atol <= 1e-15'
    endif
    if(magic=='PELEF_REACTING_EB_REFERENCE_V4') then
      read(unit,*,iostat=ios) eos_convergence
      if(ios/=0) error stop 'Invalid EOS convergence control'
      if(eos_convergence/=0.and.eos_convergence/=1) error stop 'EOS convergence control must be 0 or 1'
    endif
  case default
    error stop 'Invalid probe input magic'
  end select
  half_size=root_size/2
  read(unit,*,iostat=ios) steps,iterations,subcycle,dt
  if(ios/=0) error stop 'Invalid stepping row'
  if(steps<0.or.steps>100.or.iterations<1.or.(subcycle/=0.and.subcycle/=1).or.dt<=0) error stop 'Invalid stepping controls'
  read(unit,*,iostat=ios) spacing
  if(ios/=0.or.any(spacing<=0)) error stop 'Invalid probe spacing'
  read(unit,*,iostat=ios) plane,transport_on
  if(ios/=0.or.(transport_on/=0.and.transport_on/=1)) error stop 'Invalid EB/transport controls'
  do c=1,63
    read(unit,*,iostat=ios) fits(:,c)
    if(ios/=0) error stop 'Missing pinned Simple fit'
  enddo
  read(unit,*,iostat=ios) records
  if(ios/=0.or.records/=2*root_size**3) error stop 'Probe requires complete coarse and fine cell records'
  call load_lidryer_thermo(species,ok)
  if(.not.ok) error stop 'Cannot load LiDryer thermo'
  call load_lidryer_mechanism(reactions,ok)
  if(.not.ok) error stop "Cannot load LiDryer reactions"
  call configure_nasa7_eos_policy(species,.true.,ok,converged_inverse=eos_convergence==1)
  if(.not.ok) error stop "Invalid EOS policy"
  policy%allow_negative_energy=.false.;policy%eta2=2
  energy_correction=0
  call load_lidryer_transport(transport_names,geometries,transport_values,ok)
  if(.not.ok) error stop 'Transport source load failed'
  call load_gas_transport_data(transport_names,geometries,transport_values,transport,ok)
  if(.not.ok) error stop 'Transport database failed'
  call set_simple_transport_fits(transport,fits,ok)
  if(.not.ok) error stop 'Pinned Simple fits failed'
  call initialize_amr_patch_3d(root_size,root_size,root_size,1,half_size,1,half_size,1,half_size,2,patch(1),ok)
  if(.not.ok) error stop 'Invalid reference patch'
  call initialize_amr_patch_tree_3d([root_size,root_size,root_size],[1],patch,tree,ok)
  if(.not.ok) error stop 'Invalid reference tree'
  allocate(state%fields(2),temperature%fields(2),seen%fields(2),held%fields(2),history%fields(2))
  do node=1,2
    n=tree%nodes(node)%shape
    allocate(state%fields(node)%values(nv,n(1),n(2),n(3)),temperature%fields(node)%values(1,n(1),n(2),n(3)), &
      seen%fields(node)%values(1,n(1),n(2),n(3)))
    allocate(held%fields(node)%values(1,n(1),n(2),n(3)),history%fields(node)%values(nv-4,n(1),n(2),n(3)))
    history%fields(node)%values=0;held%fields(node)%values=0
    state%fields(node)%values=0;temperature%fields(node)%values=0;seen%fields(node)%values=0
  enddo
  do row=1,records
    read(unit,*,iostat=ios) node,i,j,k,values
    if(ios/=0.or.node<1.or.node>2) error stop 'Invalid imported cell'
    n=tree%nodes(node)%shape
    if(any([i,j,k]<1).or.any([i,j,k]>n)) error stop 'Imported cell index outside patch'
    if(seen%fields(node)%values(1,i,j,k)/=0) error stop 'Duplicate imported cell'
    seen%fields(node)%values(1,i,j,k)=1
    state%fields(node)%values(:,i,j,k)=values(:nv)
    temperature%fields(node)%values(1,i,j,k)=values(nv+2)
    held%fields(node)%values(1,i,j,k)=values(nv+1)
    ! Retain native time-zero estimates for diagnostics; the production
    ! initializer below replaces them before the first predictor.
    history%fields(node)%values(:,i,j,k)=values(nv+3:)
  enddo
  read(unit,*,iostat=ios) magic
  if(ios>=0) error stop 'Unexpected records after imported cells'
  close(unit)
  do node=1,2
    if(any(seen%fields(node)%values/=1)) error stop 'Missing imported cell'
  enddo
  call build_axis_plane_eb_geometry_3d(2*root_size,2*root_size,2*root_size, &
    0._dp,root_size*spacing(1),0._dp,root_size*spacing(2),0._dp,root_size*spacing(3),'z',plane,finest,ok)
  if(.not.ok) error stop 'Invalid EB geometry'
  call initialize_reactive_eb_tree_3d(species,tree,finest,state%fields(1)%values, &
    temperature%fields(1)%values(1,:,:,:),eb,ok)
  if(.not.ok) error stop 'EB tree initialization failed'
  do node=1,2
    eb%nodes(node)%state=state%fields(node)%values
    eb%nodes(node)%temperature=temperature%fields(node)%values(1,:,:,:)
    physical(node)=eb%nodes(node)%geometry
  enddo
  call build_axis_plane_eb_geometry_3d(2*(root_size+12),2*(root_size+12),2*(root_size+12), &
    -6*spacing(1),(root_size+6)*spacing(1),-6*spacing(2),(root_size+6)*spacing(2), &
    -6*spacing(3),(root_size+6)*spacing(3),'z',plane,padded,ok)
  if(.not.ok) error stop 'Invalid padded EB geometry'
  call build_eb_tree_extended_geometry_3d(tree,physical,padded,extended,ok,halo_width=6)
  if(.not.ok) error stop 'Invalid EB halo geometry'
  faces='periodic';faces(5:6)='slip_wall';force=0
  if(audit_initial) then
    if(eos_convergence/=1) error stop 'Initial audit requires converged V4 EOS'
    if(audit_shear_initial) then
      call initialize_analytic_shear()
      call audit_shear_spatial()
    else if(audit_native) then
      if(smooth_contact) call initialize_smooth_contact()
      if(audit_shear) call initialize_analytic_shear()
      call audit_native_pressure()
    else if(audit_mixing) then
      if(smooth_contact) call initialize_smooth_contact()
      call audit_conservative_mixing()
    else if(audit_evolving) then
      if(smooth_contact) call initialize_smooth_contact()
      call audit_pressure_evolution()
    else
      call audit_initial_operator()
    endif
    stop
  endif
  call initialize_eb_tree_reaction_source_3d(species,reactions,eb,dt,subcycle==1, &
    reaction_rtol,reaction_atol,history,ok,evaluate_reaction=cached_initial_reaction,internal_energy=held)
  if(.not.ok) error stop 'Initial EB reaction estimate refused'
  write(*,'(a,3(i0,1x))') 'Initial reaction cache: entries/hits/solves ',initial_entries,initial_hits,initial_solves
  if(audit_stage_flux) then
        mol_trace_event=0;sync_trace_event=0
        open(newunit=sync_trace_unit,file=trim(prefix)//'_default'//'_sync_flux.csv',status='replace')
        write(sync_trace_unit,'(a)') 'step,event,phase,depth,time_start_s,interval_s,node,probe,i,j,k,axis,side,'// &
          'vfrac,aperture,pressure_Pa,held_energy,state_mx,state_my,state_mz,flux_mx,flux_my,flux_mz,flux_held,'// &
          'raw_delta_mx,raw_delta_my,raw_delta_mz,applied_delta_mx,applied_delta_my,applied_delta_mz'
        open(newunit=mol_trace_unit,file=trim(prefix)//'_default'//'_stage_flux.csv',status='replace')
        write(mol_trace_unit,'(a)') 'step,event,stage,depth,time_start_s,interval_s,node,probe,i,j,k,axis,side,'// &
          'vfrac,aperture,pressure_cell_Pa,pressure_left_Pa,pressure_right_Pa,'// &
          'raw_mx,raw_my,raw_mz,centroid_mx,centroid_my,centroid_mz,tree_mx,tree_my,tree_mz,'// &
          'rate_mx,rate_my,rate_mz,wall_mx,wall_my,wall_mz,state_mx,state_my,state_mz'
  endif
  do step=0,steps
    mol_trace_step=step+1
    write(path,'(a,"_",i3.3,".csv")') trim(prefix),step
    open(newunit=unit,file=trim(path),status='replace',action='write',iostat=ios)
    if(ios/=0) error stop 'Cannot open probe output'
    write(unit,'(a)',advance='no') 'node,level,i,j,k,time_s,density,xmom,ymom,zmom,rho_E'
    do c=1,size(species)
      write(unit,'(a)',advance='no') ',rho_'//trim(species(c)%name)
    enddo
    write(unit,'(a)',advance='no') ',rho_e,Temp'
    do c=1,size(species)
      write(unit,'(a)',advance='no') ',rho_omega_'//trim(species(c)%name)
    enddo
    write(unit,'(a)') ',rhoe_dot,vfrac'
    do node=1,2
      n=tree%nodes(node)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        write(unit,'(5(i0,","),*(es25.16e3,:,","))') node,tree%nodes(node)%level,i,j,k,step*dt, &
          eb%nodes(node)%state(:,i,j,k),held%fields(node)%values(1,i,j,k),eb%nodes(node)%temperature(i,j,k), &
          history%fields(node)%values(:,i,j,k),eb%nodes(node)%geometry%volume_fraction(i,j,k)
      enddo
      enddo
      enddo
    enddo
    close(unit)
    write(path,'(a,"_",i3.3,".budget.csv")') trim(prefix),step
    open(newunit=unit,file=trim(path),status="replace",action="write")
    write(unit,'(a)') "step,time_s,numerical_energy_J"
    write(unit,'(i0,2(",",es25.16e3))') step,step*dt,energy_correction
    close(unit)
    if(step==steps) exit
    budget=0;wall=0;chemical=0;source=0;theta=1
    if(subcycle==1) then
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,iterations,.true., &
        [transport_on==1,transport_on==1,transport_on==1,transport_on==1],force,reaction_rtol,reaction_atol, &
        eb,extended,faces,budget,wall,chemical,source,theta,ok,history,redistribution='none', &
        evaluate_reaction=controlled_reaction,mol_mode=.true.,internal_energy=held, &
        energy_correction=energy_correction,energy_policy=policy,internal_energy_flux=held_flux, &
        separate_momentum=momentum_reacting,consistent_velocity=velocity_reacting,paired_pressure=paired_reacting, &
        observe_mol=mol_observer,observe_mol_sync=sync_observer)
    else
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,iterations,.true., &
        [transport_on==1,transport_on==1,transport_on==1,transport_on==1],force,reaction_rtol,reaction_atol, &
        eb,extended,faces,budget,wall,chemical,source,theta,ok,history,redistribution='none', &
        evaluate_reaction=controlled_reaction,mol_mode=.true.,internal_energy=held, &
        energy_correction=energy_correction,energy_policy=policy,internal_energy_flux=held_flux, &
        separate_momentum=momentum_reacting,consistent_velocity=velocity_reacting,paired_pressure=paired_reacting, &
        observe_mol=mol_observer)
    endif
    if(.not.ok) error stop 'Reference-state MOL step rejected'
    ! This case closes root level zero with eta2=2. Covered coarse cells
    ! remain interpolation donors: their averaged held energy must be made
    ! consistent with total energy, without resetting the fine-level leaves.
    do k=1,half_size
    do j=1,half_size
    do i=1,half_size
      if(eb%nodes(1)%geometry%volume_fraction(i,j,k)<=0) cycle
      covered_energy_gap=eb%nodes(1)%state(5,i,j,k)-held%fields(1)%values(1,i,j,k)- &
        .5_dp*sum(eb%nodes(1)%state(2:4,i,j,k)**2)/eb%nodes(1)%state(1,i,j,k)
      covered_energy_scale=max(1._dp,abs(eb%nodes(1)%state(5,i,j,k)),abs(held%fields(1)%values(1,i,j,k)))
      if(abs(covered_energy_gap)>64*epsilon(1._dp)*covered_energy_scale) &
        error stop 'Closing covered coarse energy inconsistent'
    enddo
    enddo
    enddo
  enddo
  write(*,'(a,2(es25.16e3,1x))') 'CVODE absolute density tolerance range: ',cvode_atol_min,cvode_atol_max
contains
  ! Same imported control, ordinary held-MOL integration, no chemistry or
  ! molecular transport. Each row independently closes the composite gas
  ! inventory against returned boundary/wall/source and numerical-energy terms.
  subroutine trace_mol_sync(phase,depth,time,interval,snapshot,stored,flux,scalar_flux, &
      raw_increment,applied_increment)
    integer, intent(in) :: phase,depth
    real(dp), intent(in) :: time,interval
    type(eb_reactive_tree_3d), intent(in) :: snapshot
    type(amr_tree_state_3d), intent(in) :: stored
    type(amr_tree_state_3d), optional, intent(in) :: raw_increment,applied_increment
    type(eb_integrated_flux_3d), intent(in) :: flux(:),scalar_flux(:)
    integer :: p,cell(3),f(3),a,side,probe_id,zcut
    real(dp) :: aperture,face(3),held_face,q(nv),t,cs,raw_delta(3),applied_delta(3)
    logical :: good
    sync_trace_event=sync_trace_event+1
    do probe_id=1,13
      p=2
      if(probe_id>=3.and.probe_id<=5) p=1
      zcut=max(1,int(plane/physical(p)%dz)+1)
      select case(probe_id)
      case(1);cell=[1,root_size/2,zcut]
      case(2);cell=[1,root_size/2+1,root_size]
      case(3);cell=[root_size,root_size/4,zcut]
      case(4);cell=[root_size,root_size/4+1,root_size/2]
      case(5);cell=[1,root_size/4+1,root_size/2+1]
      case(6);cell=[1,1,zcut]
      case(7);cell=[1,1,root_size]
      case(8);cell=[1,1,zcut+1]
      case(9);cell=[1,1,root_size-1]
      ! Saved 32-cell original-reference extrema; retain historical probe IDs.
      case(10:13)
        if(root_size/=32) cycle
        select case(probe_id)
        case(10);cell=[25,32,26]
        case(11);cell=[25,32,25]
        case(12);cell=[26,32,30]
        case(13);cell=[15,23,32]
        end select
      end select
      if(any(cell>root_size)) cycle
      if(physical(p)%volume_fraction(cell(1),cell(2),cell(3))<=0) cycle
      call reactive_conserved_to_primitive(species,snapshot%nodes(p)%state(:,cell(1),cell(2),cell(3)), &
        snapshot%nodes(p)%temperature(cell(1),cell(2),cell(3)),q,t,cs,good, &
        internal_energy_density=stored%fields(p)%values(1,cell(1),cell(2),cell(3)))
      if(.not.good) error stop 'Sync diagnostic EOS refused'
      raw_delta=0;applied_delta=0
      if(phase==1) then
        if(.not.present(raw_increment).or..not.present(applied_increment)) error stop 'Missing closure increments'
        raw_delta=raw_increment%fields(p)%values(2:4,cell(1),cell(2),cell(3))
        applied_delta=applied_increment%fields(p)%values(2:4,cell(1),cell(2),cell(3))
      endif
      do a=1,3
      do side=0,1
        f=cell;f(a)=f(a)-1+side
        select case(a)
        case(1)
          aperture=physical(p)%x_face_fraction(f(1),f(2),f(3))
          face=flux(p)%x(2:4,f(1),f(2),f(3));held_face=scalar_flux(p)%x(1,f(1),f(2),f(3))
        case(2)
          aperture=physical(p)%y_face_fraction(f(1),f(2),f(3))
          face=flux(p)%y(2:4,f(1),f(2),f(3));held_face=scalar_flux(p)%y(1,f(1),f(2),f(3))
        case(3)
          aperture=physical(p)%z_face_fraction(f(1),f(2),f(3))
          face=flux(p)%z(2:4,f(1),f(2),f(3));held_face=scalar_flux(p)%z(1,f(1),f(2),f(3))
        end select
        write(sync_trace_unit,'(4(i0,","),2(es25.16e3,","),7(i0,","),*(es25.16e3,:,","))') &
          mol_trace_step,sync_trace_event,phase,depth,(mol_trace_step-1)*dt+time,interval, &
          p,probe_id,cell,a,side,physical(p)%volume_fraction(cell(1),cell(2),cell(3)),aperture,q(5), &
          stored%fields(p)%values(1,cell(1),cell(2),cell(3)), &
          snapshot%nodes(p)%state(2:4,cell(1),cell(2),cell(3)),face,held_face,raw_delta,applied_delta
      enddo
      enddo
    enddo
    flush(sync_trace_unit)
  end subroutine trace_mol_sync

  subroutine trace_mol_stage(stage,depth,time,interval,snapshot,diagnostics,flux)
    integer, intent(in) :: stage,depth
    real(dp), intent(in) :: time,interval
    type(eb_reactive_tree_3d), intent(in) :: snapshot
    type(eb_mol_diagnostics_3d), intent(in) :: diagnostics(:)
    type(eb_integrated_flux_3d), intent(in) :: flux(:)
    integer :: p,h,g,cell(3),f(3),r(3),v(3),a,side,probe_id
    real(dp) :: aperture,pair(2),raw(3),centroid(3),tree_face(3)
    mol_trace_event=mol_trace_event+1
    p=2
    if(.not.allocated(diagnostics(p)%rate)) return
    h=(size(diagnostics(p)%cell_pressure,1)-root_size)/2
    g=h-3
    ! IDs 1:5 retain their historical meaning in the spatial/sync tables.
    ! IDs 6:9 cover the candidate extrema and neighboring cut/full cells.
    do probe_id=1,13
      select case(probe_id)
      case(1);cell=[1,root_size/2,max(1,int(plane/physical(p)%dz)+1)]
      case(2);cell=[1,root_size/2+1,root_size]
      case(3:5);cycle
      case(6);cell=[1,1,max(1,int(plane/physical(p)%dz)+1)]
      case(7);cell=[1,1,root_size]
      case(8);cell=[1,1,max(1,int(plane/physical(p)%dz)+1)+1]
      case(9);cell=[1,1,root_size-1]
      ! Saved 32-cell original-reference extrema; retain historical probe IDs.
      case(10:13)
        if(root_size/=32) cycle
        select case(probe_id)
        case(10);cell=[25,32,26]
        case(11);cell=[25,32,25]
        case(12);cell=[26,32,30]
        case(13);cell=[15,23,32]
        end select
      end select
      if(any(cell>root_size)) cycle
      if(physical(p)%volume_fraction(cell(1),cell(2),cell(3))<=0) cycle
      do a=1,3
      do side=0,1
        v=cell;v(a)=v(a)-1+side
        f=v+g;r=f+1
        select case(a)
        case(1)
          aperture=physical(p)%x_face_fraction(v(1),v(2),v(3))
          pair=diagnostics(p)%pressure%x(:,r(1),r(2),r(3))
          raw=diagnostics(p)%raw_flux%x(2:4,r(1),r(2),r(3))
          centroid=diagnostics(p)%centroid_flux%x(2:4,f(1),f(2),f(3))
          tree_face=flux(p)%x(2:4,v(1),v(2),v(3))/interval
        case(2)
          aperture=physical(p)%y_face_fraction(v(1),v(2),v(3))
          pair=diagnostics(p)%pressure%y(:,r(1),r(2),r(3))
          raw=diagnostics(p)%raw_flux%y(2:4,r(1),r(2),r(3))
          centroid=diagnostics(p)%centroid_flux%y(2:4,f(1),f(2),f(3))
          tree_face=flux(p)%y(2:4,v(1),v(2),v(3))/interval
        case(3)
          aperture=physical(p)%z_face_fraction(v(1),v(2),v(3))
          pair=diagnostics(p)%pressure%z(:,r(1),r(2),r(3))
          raw=diagnostics(p)%raw_flux%z(2:4,r(1),r(2),r(3))
          centroid=diagnostics(p)%centroid_flux%z(2:4,f(1),f(2),f(3))
          tree_face=flux(p)%z(2:4,v(1),v(2),v(3))/interval
        end select
        write(mol_trace_unit,'(4(i0,","),2(es25.16e3,","),7(i0,","),*(es25.16e3,:,","))') &
          mol_trace_step,mol_trace_event,stage,depth,(mol_trace_step-1)*dt+time,interval, &
          p,probe_id,cell,a,side,physical(p)%volume_fraction(cell(1),cell(2),cell(3)),aperture, &
          diagnostics(p)%cell_pressure(cell(1)+h,cell(2)+h,cell(3)+h),pair,raw,centroid,tree_face, &
          diagnostics(p)%rate(2:4,cell(1)+g,cell(2)+g,cell(3)+g), &
          diagnostics(p)%wall_rate(2:4,cell(1)+g,cell(2)+g,cell(3)+g), &
          snapshot%nodes(p)%state(2:4,cell(1),cell(2),cell(3))
      enddo
      enddo
    enddo
    flush(mol_trace_unit)
  end subroutine trace_mol_stage

  subroutine audit_native_pressure()
    type(eb_reactive_tree_3d) :: initial_tree
    type(amr_tree_state_3d) :: initial_held
    real(dp) :: initial(nv),total(nv),ledger(nv),error,maxp,q(nv),t,cs,maxgap,gap
    real(dp) :: native_correction,pressure_work
    integer :: mode,it,p,x,y,z,u,v,c,w,modes
    logical :: good,paired,separate,consistent
    character(len=16) :: label
    initial_tree=eb;initial_held=held
    call integrate_reactive_eb_tree_3d(eb,initial,good)
    if(.not.good) error stop 'Native initial inventory refused'
    modes=2
    if(audit_shear) modes=4
    do mode=1,modes
      separate=.true.;consistent=.true.;paired=mode==2;label='velocity'
      if(paired) label='paired'
      if(audit_shear) then
        separate=mode>=2;consistent=mode>=3;paired=mode==4
        select case(mode)
        case(1);label='default'
        case(2);label='momentum'
        case(3);label='velocity'
        case(4);label='paired'
        end select
      endif
      eb=initial_tree;held=initial_held;ledger=0;native_correction=0
      if(audit_stage_flux) then
        mol_trace_event=0;sync_trace_event=0
        open(newunit=sync_trace_unit,file=trim(prefix)//'_'//trim(label)//'_sync_flux.csv',status='replace')
        write(sync_trace_unit,'(a)') 'step,event,phase,depth,time_start_s,interval_s,node,probe,i,j,k,axis,side,'// &
          'vfrac,aperture,pressure_Pa,held_energy,state_mx,state_my,state_mz,flux_mx,flux_my,flux_mz,flux_held,'// &
          'raw_delta_mx,raw_delta_my,raw_delta_mz,applied_delta_mx,applied_delta_my,applied_delta_mz'
        open(newunit=mol_trace_unit,file=trim(prefix)//'_'//trim(label)//'_stage_flux.csv',status='replace')
        write(mol_trace_unit,'(a)') 'step,event,stage,depth,time_start_s,interval_s,node,probe,i,j,k,axis,side,'// &
          'vfrac,aperture,pressure_cell_Pa,pressure_left_Pa,pressure_right_Pa,'// &
          'raw_mx,raw_my,raw_mz,centroid_mx,centroid_my,centroid_mz,tree_mx,tree_my,tree_mz,'// &
          'rate_mx,rate_my,rate_mz,wall_mx,wall_my,wall_mz,state_mx,state_my,state_mz'
      endif
      open(newunit=u,file=trim(prefix)//'_'//trim(label)//'_native.csv',status='replace')
      open(newunit=v,file=trim(prefix)//'_'//trim(label)//'_native_inventory.csv',status='replace')
      write(u,'(a)') 'step,time_s,max_pressure_offset_Pa,max_budget_residual,numerical_energy_J,max_E_H_KE_J_m3'
      write(v,'(a)') 'step,component,initial,final,physical_change,numerical_change,residual'
      if(audit_shear) then
        open(newunit=w,file=trim(prefix)//'_'//trim(label)//'_analytic.csv',status='replace')
        write(w,'(a)') 'step,time_s,component,L1,L2,Linf'
      endif
      do it=0,steps
        call integrate_reactive_eb_tree_3d(eb,total,good)
        if(.not.good) error stop 'Native inventory refused'
        error=0
        do c=1,nv
          gap=total(c)-initial(c)-ledger(c)
          if(c==5) gap=gap-native_correction
          if(c>=2.and.c<=4) then
            error=max(error,abs(gap)/max(abs(initial(c)),abs(total(c)),abs(initial(1))))
          else
            error=max(error,abs(gap)/max(abs(initial(c)),abs(total(c)),abs(initial(1))*1.e-14_dp))
          endif
          write(v,'(i0,",",i0,5(",",es25.16e3))') it,c,initial(c),total(c),ledger(c), &
            merge(native_correction,0._dp,c==5),gap
        enddo
        ! Transverse momentum may start at zero: use the gas mass inventory
        ! as its scale in the actual gate, retaining raw residuals in the CSV.
        do c=2,4
          gap=total(c)-initial(c)-ledger(c)
          if(abs(gap)>1.e-12_dp*max(abs(initial(c)),abs(initial(1)))) error stop 'Native momentum budget failed'
        enddo
        do c=1,nv
          if(c>=2.and.c<=4) cycle
          gap=total(c)-initial(c)-ledger(c)
          if(c==5) gap=gap-native_correction
          if(abs(gap)>1.e-12_dp*max(abs(initial(c)),abs(total(c)),abs(initial(1))*1.e-14_dp)) &
            error stop 'Native material/energy budget failed'
        enddo
        maxp=0;maxgap=0
        do p=1,2
        do z=1,root_size
        do y=1,root_size
        do x=1,root_size
          if(p==1.and.all([x,y,z]<=half_size)) cycle
          if(physical(p)%volume_fraction(x,y,z)<=0) cycle
          call reactive_conserved_to_primitive(species,eb%nodes(p)%state(:,x,y,z), &
            eb%nodes(p)%temperature(x,y,z),q,t,cs,good,internal_energy_density=held%fields(p)%values(1,x,y,z))
          if(.not.good) error stop 'Native pressure EOS refused'
          maxp=max(maxp,abs(q(5)-101325._dp))
          gap=eb%nodes(p)%state(5,x,y,z)-held%fields(p)%values(1,x,y,z)- &
            .5_dp*sum(eb%nodes(p)%state(2:4,x,y,z)**2)/eb%nodes(p)%state(1,x,y,z)
          maxgap=max(maxgap,abs(gap))
        enddo
        enddo
        enddo
        enddo
        write(u,'(i0,5(",",es25.16e3))') it,it*dt,maxp,error,native_correction,maxgap
        flush(u);flush(v)
        if(audit_shear) call write_analytic_norms(initial_tree,initial_held,it,w,trim(label))
        if(it==steps) exit
        mol_trace_step=it+1
        budget=0;wall=0;chemical=0;source=0;theta=1;pressure_work=0
        if(subcycle==1) then
          call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,iterations,.false., &
            [.false.,.false.,.false.,.false.],force,reaction_rtol,reaction_atol,eb,extended,faces, &
            budget,wall,chemical,source,theta,good,mol_mode=.true.,internal_energy=held, &
            energy_correction=native_correction,energy_policy=policy,internal_energy_flux=held_flux, &
            internal_energy_pressure_work_change=pressure_work,separate_momentum=separate, &
            consistent_velocity=consistent,paired_pressure=paired,connected_directions=audit_connected.and.separate, &
            observe_mol=mol_observer,observe_mol_sync=sync_observer)
        else
          call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,iterations,.false., &
            [.false.,.false.,.false.,.false.],force,reaction_rtol,reaction_atol,eb,extended,faces, &
            budget,wall,chemical,source,theta,good,mol_mode=.true.,internal_energy=held, &
            energy_correction=native_correction,energy_policy=policy,internal_energy_flux=held_flux, &
            internal_energy_pressure_work_change=pressure_work,separate_momentum=separate, &
            consistent_velocity=consistent,paired_pressure=paired,connected_directions=audit_connected.and.separate, &
            observe_mol=mol_observer)
        endif
        if(.not.good) error stop 'Native pressure evolution refused'
        ledger=ledger+budget+wall+chemical+source
      enddo
      if(audit_stage_flux) then
        close(mol_trace_unit);close(sync_trace_unit)
      endif
      close(u);close(v)
      if(audit_shear) close(w)
    enddo
    eb=initial_tree;held=initial_held
  end subroutine audit_native_pressure

  ! Observe the production spatial operator before any reaction callback.
  ! Face values are positive-coordinate, intensive SI fluxes (divide out dt).
  ! Output every face so downstream analysis can area-average fine faces.
  ! Initial spatial operator only: no chemistry, advance, reset or reflux.
  subroutine audit_shear_spatial()
    character(len=512) :: saved_prefix
    character(len=8), parameter :: labels(4)=[character(len=8) :: 'default','momentum','velocity','paired']
    integer :: mode
    saved_prefix=prefix
    do mode=1,4
      audit_separate=mode>=2;audit_consistent=mode>=3;audit_pair=mode==4
      prefix=trim(saved_prefix)//'_'//trim(labels(mode))
      call audit_initial_operator()
    enddo
    prefix=saved_prefix
  end subroutine audit_shear_spatial

  subroutine audit_initial_operator()
    type(amr_tree_state_3d) :: rates,old_halos,aux_rates,aux_halos
    type(eb_integrated_flux_3d), allocatable :: fluxes(:),aux_fluxes(:)
    real(dp) :: q(reactive_nprim(size(species))),t,cs,f(nv),aperture
    integer :: u,p,a,x,y,z,c,idx(3),lo(3),hi(3),ncell(3)
    logical :: good
    call compute_mol_eb_tree_hydro_3d(species,eb,extended,dt,2,rates,fluxes,good, &
      boundary_faces=faces,old_state_halos=old_halos,internal_energy_rate=aux_rates, &
      internal_energy_flux=aux_fluxes,internal_energy_state=held, &
      internal_energy_old_state_halos=aux_halos,separate_momentum=audit_separate, &
      consistent_velocity=audit_consistent,paired_pressure=audit_pair.and.audit_shear_initial)
    if(.not.good) error stop 'Initial spatial audit refused'
    open(newunit=u,file=trim(prefix)//'_initial_pressure.csv',status='replace')
    if(audit_shear_initial) then
      write(u,'(a)') 'node,i,j,k,vfrac,pressure_Pa,temperature_K,rho_e_J_m3,u_m_s,v_m_s,w_m_s'
    else
      write(u,'(a)') 'node,i,j,k,vfrac,pressure_Pa,temperature_K,rho_e_J_m3'
    endif
    do p=1,size(eb%nodes)
      ncell=tree%nodes(p)%shape
      do z=1,ncell(3)+12
      do y=1,ncell(2)+12
      do x=1,ncell(1)+12
        if(extended(p)%volume_fraction(x,y,z)<=0) cycle
        call reactive_conserved_to_primitive(species,old_halos%fields(p)%values(:,x,y,z), &
          1000._dp,q,t,cs,good,internal_energy_density=aux_halos%fields(p)%values(1,x,y,z))
        if(.not.good) error stop 'Audit EOS refused'
        if(audit_shear_initial) then
        write(u,'(4(i0,","),*(es25.16e3,:,","))') p,x-6,y-6,z-6, &
          extended(p)%volume_fraction(x,y,z),q(5),t,aux_halos%fields(p)%values(1,x,y,z),q(2:4)
        else
        write(u,'(4(i0,","),*(es25.16e3,:,","))') p,x-6,y-6,z-6, &
          extended(p)%volume_fraction(x,y,z),q(5),t,aux_halos%fields(p)%values(1,x,y,z)
        endif
      enddo
      enddo
      enddo
    enddo
    close(u)
    open(newunit=u,file=trim(prefix)//'_initial_flux.csv',status='replace')
    if(audit_shear_initial) then
    write(u,'(a)') 'node,axis,i,j,k,aperture,mass_kg_m2_s,normal_momentum_Pa,energy_W_m2'// &
        ',x_momentum_Pa,y_momentum_Pa,z_momentum_Pa'
    else
    write(u,'(a)') 'node,axis,i,j,k,aperture,mass_kg_m2_s,normal_momentum_Pa,energy_W_m2'
    endif
    do p=1,size(eb%nodes)
      ncell=tree%nodes(p)%shape
      do a=1,3
        lo=1;lo(a)=0;hi=ncell
        do z=lo(3),hi(3)
        do y=lo(2),hi(2)
        do x=lo(1),hi(1)
          idx=[x,y,z]
          select case(a)
          case(1)
            f=fluxes(p)%x(:,x,y,z)/dt
            aperture=eb%nodes(p)%geometry%x_face_fraction(x,y,z)
          case(2)
            f=fluxes(p)%y(:,x,y,z)/dt
            aperture=eb%nodes(p)%geometry%y_face_fraction(x,y,z)
          case(3)
            f=fluxes(p)%z(:,x,y,z)/dt
            aperture=eb%nodes(p)%geometry%z_face_fraction(x,y,z)
          end select
          if(audit_shear_initial) then
          write(u,'(5(i0,","),*(es25.16e3,:,","))') p,a,idx,aperture,f(1),f(a+1),f(5),f(2:4)
          else
          write(u,'(5(i0,","),*(es25.16e3,:,","))') p,a,idx,aperture,f(1),f(a+1),f(5)
          endif
        enddo
        enddo
        enddo
      enddo
    enddo
    close(u)
    if(audit_shear_initial) then
      open(newunit=u,file=trim(prefix)//'_initial_rates.csv',status='replace')
      write(u,'(a)') 'node,i,j,k,component,rate_SI_per_s'
      do p=1,size(eb%nodes)
      do z=1,root_size
      do y=1,root_size
      do x=1,root_size
        if(physical(p)%volume_fraction(x,y,z)<=0) cycle
        if(p==1.and.all([x,y,z]<=half_size)) cycle
        do c=1,nv
          write(u,'(5(i0,","),es25.16e3)') p,x,y,z,c,rates%fields(p)%values(c,x,y,z)
        enddo
      enddo
      enddo
      enddo
      enddo
      close(u)
      return
    endif
    if(audit_pair) call audit_pressure_step(rates,aux_rates,fluxes,aux_fluxes,'_baseline')
    call audit_initial_stencils(old_halos,aux_halos)
    call audit_pressure_ghosts(old_halos,aux_halos,fluxes,.false.)
    if(audit_pair) call audit_pressure_ghosts(old_halos,aux_halos,fluxes,.true.)
    write(*,'(a,l1)') 'Initial spatial audit complete; separate_momentum=',audit_separate
  end subroutine audit_initial_operator
  ! Experimental ghost-only primitive pressure reconstruction. Conserved
  ! valid cells/restriction are untouched; ghost energy is NOT a conservative
  ! prolongation. Reflux/trajectory qualification remains separate.
  subroutine audit_pressure_ghosts(halos,aux,baseline_flux,paired)
    type(amr_tree_state_3d), intent(in) :: halos,aux
    type(eb_integrated_flux_3d), intent(in) :: baseline_flux(:)
    logical, intent(in) :: paired
    character(len=40) :: suffix
    type(eb_integrated_flux_3d) :: coarse_flux,step_flux(2),step_held_flux(2)
    type(amr_tree_state_3d) :: step_rate,step_held_rate
    type(eb_integrated_flux_3d) :: flux,held_face
    real(dp), allocatable :: gas(:,:,:,:),rhoe(:,:,:),temp(:,:,:),rate(:,:,:,:),erate(:,:,:),rate_halo(:,:,:,:)
    real(dp) :: stencil_values(nv+2,6)
    real(dp) :: value(nv+2),neighbor(nv+2),q(nv),t,cs,ke,delta,target(1),pn(1,6),offset(3)
    real(dp) :: boundary(nv),wall_flux(nv),f(nv),aperture,max_delta,max_pressure_error
    integer(int64) :: origins(3,2),domains(3,2),global(3)
    integer :: reflected(3)
    integer :: shapes(3,2),requests(4,6),signs(3,6),idx(3),src(3),n,a,x,y,z,j,u,p,lo(3),hi(3),changes
    logical :: good,linear
    origins=0;domains(:,1)=root_size;domains(:,2)=2*root_size;shapes=root_size;n=root_size
    suffix='_pressure_candidate'
    if(paired) suffix='_pressure_pair_candidate'
    gas=halos%fields(2)%values;rhoe=aux%fields(2)%values(1,:,:,:)
    allocate(temp(n+12,n+12,n+12),source=1000._dp)
    changes=0;max_delta=0;max_pressure_error=0
    do z=1,n+12
    do y=1,n+12
    do x=1,n+12
      if(extended(2)%volume_fraction(x,y,z)<=0) cycle
      idx=[x,y,z]-6;reflected=1
      if(idx(3)<1.or.idx(3)>2*n) reflected(3)=-1
      global=int(idx-1,int64);global(1:2)=modulo(global(1:2),domains(1:2,2))
      if(global(3)<0) global(3)=-global(3)-1
      if(global(3)>=2*n) global(3)=4*n-global(3)-1
      if(all(global>=0).and.all(global<n)) cycle ! same-level/reflected fine donor
      src=int(global/2)+1
      call pack(1,src,value);target(1)=donor_pressure(1,src,paired)
      call plan_eb_linear_halo_3d([0,1],shapes,origins,domains,1,src,2,global, &
        extended(1),6,faces,requests,signs,offset,linear,good)
      if(.not.good) error stop 'Pressure ghost planner refused'
      if(linear) then
        do j=1,6
          call pack(requests(1,j),requests(2:4,j),neighbor)
          stencil_values(:,j)=neighbor
          stencil_values(2:4,j)=stencil_values(2:4,j)*real(signs(:,j),dp)
          pn(1,j)=donor_pressure(requests(1,j),requests(2:4,j),paired)
        enddo
        call interpolate_eb_linear_halo_3d(target,pn,offset,good)
        if(.not.good) error stop 'Pressure ghost interpolation refused'
      endif
      if(consistent_contact.and.paired) then
        if(linear) then
          call interpolate_eb_linear_halo_3d(value,stencil_values,offset,good, &
            separate_momentum=.true.,consistent_velocity=.true.)
          if(.not.good) error stop 'Consistent velocity interpolation refused'
        endif
        gas(2:4,x,y,z)=value(2:4)*real(reflected,dp)
        ke=.5_dp*sum(gas(2:4,x,y,z)**2)/gas(1,x,y,z)
        gas(5,x,y,z)=rhoe(x,y,z)+ke
      endif
      temp(x,y,z)=1000._dp
      call reconstruct_eb_pressure_state_3d(species,gas(:,x,y,z),rhoe(x,y,z),temp(x,y,z),target(1),delta,good)
      if(.not.good) error stop 'Pressure ghost reconstruction refused'
      changes=changes+1;max_delta=max(max_delta,abs(delta))
      call reactive_conserved_to_primitive(species,gas(:,x,y,z),temp(x,y,z),q,t,cs,good, &
        internal_energy_density=rhoe(x,y,z))
      if(.not.good) error stop 'Pressure ghost readback refused'
      max_pressure_error=max(max_pressure_error,abs(q(5)-target(1)))
    enddo
    enddo
    enddo
    if(any(gas(:,7:n+6,7:n+6,7:n+6)/=halos%fields(2)%values(:,7:n+6,7:n+6,7:n+6))) &
      error stop 'Pressure ghost modified valid gas'
    if(any(rhoe(7:n+6,7:n+6,7:n+6)/=aux%fields(2)%values(1,7:n+6,7:n+6,7:n+6))) &
      error stop 'Pressure ghost modified valid held energy'
    if(any(gas(1,:,:,:)/=halos%fields(2)%values(1,:,:,:)).or. &
      any(gas(6:,:,:,:)/=halos%fields(2)%values(6:,:,:,:))) error stop 'Pressure ghost modified material'
    if(.not.consistent_contact) then
    if(any(gas(:4,:,:,:)/=halos%fields(2)%values(:4,:,:,:)).or. &
      any(gas(6:,:,:,:)/=halos%fields(2)%values(6:,:,:,:))) error stop 'Pressure ghost modified material or momentum'
    endif
    allocate(rate(nv,n,n,n),erate(n,n,n));rate=0;erate=0
    call compute_mol_eb_hydro_3d(species,extended(2),gas,temp,dt,2,rate,flux,boundary,wall_flux,good, &
      spatial_rate_halo=rate_halo,internal_energy_rate=erate,internal_energy_integrated_flux=held_face, &
      internal_energy_density=rhoe)
    if(.not.good) error stop 'Pressure ghost spatial operator refused'
    if(paired) then
      allocate(step_rate%fields(2),step_held_rate%fields(2))
      step_rate%fields(2)%values=rate
      allocate(step_held_rate%fields(2)%values(1,n,n,n))
      step_held_rate%fields(2)%values(1,:,:,:)=erate
      step_flux(2)=flux;step_held_flux(2)=held_face
    endif
    if(.not.audit_evolving) then
    open(newunit=u,file=trim(prefix)//trim(suffix)//'_initial_pressure.csv',status='replace')
    write(u,'(a)') 'node,i,j,k,vfrac,pressure_Pa,temperature_K,rho_e_J_m3'
    do p=1,2
    do z=1,n+12
    do y=1,n+12
    do x=1,n+12
      if(extended(p)%volume_fraction(x,y,z)<=0) cycle
      if(p==2) then
        value(:nv)=gas(:,x,y,z);value(nv+2)=rhoe(x,y,z)
      else
        value(:nv)=halos%fields(p)%values(:,x,y,z);value(nv+2)=aux%fields(p)%values(1,x,y,z)
        if(paired) then
          src=[x,y,z]-6;src(1:2)=modulo(src(1:2)-1,n)+1
          if(src(3)<1) src(3)=1-src(3)
          if(src(3)>n) src(3)=2*n+1-src(3)
          if(all(src<=half_size)) then
            t=1000._dp
            call reconstruct_eb_pressure_state_3d(species,value(:nv),value(nv+2),t, &
              donor_pressure(1,src,.true.),delta,good)
            if(.not.good) error stop 'Paired pressure output reconstruction refused'
          endif
        endif
      endif
      call reactive_conserved_to_primitive(species,value(:nv),1000._dp,q,t,cs,good, &
        internal_energy_density=value(nv+2))
      if(.not.good) error stop 'Pressure candidate output EOS refused'
      write(u,'(4(i0,","),*(es25.16e3,:,","))') p,x-6,y-6,z-6,extended(p)%volume_fraction(x,y,z),q(5),t,value(nv+2)
    enddo
    enddo
    enddo
    enddo
    close(u)
    endif
    coarse_flux=baseline_flux(1)
    if(paired) then
      gas=halos%fields(1)%values;rhoe=aux%fields(1)%values(1,:,:,:);temp=1000
      do z=1,n+12
      do y=1,n+12
      do x=1,n+12
        if(extended(1)%volume_fraction(x,y,z)<=0) cycle
        idx=[x,y,z]-6;src=idx
        src(1:2)=modulo(src(1:2)-1,n)+1
        if(src(3)<1) src(3)=1-src(3)
        if(src(3)>n) src(3)=2*n+1-src(3)
        if(any(src>half_size)) cycle
        target(1)=donor_pressure(1,src,.true.)
        temp(x,y,z)=1000._dp
        call reconstruct_eb_pressure_state_3d(species,gas(:,x,y,z),rhoe(x,y,z),temp(x,y,z),target(1),delta,good)
        if(.not.good) error stop 'Paired coarse reconstruction refused'
      enddo
      enddo
      enddo
      do z=1,n
      do y=1,n
      do x=1,n
        if(all([x,y,z]<=half_size)) cycle
        if(any(gas(:,x+6,y+6,z+6)/=halos%fields(1)%values(:,x+6,y+6,z+6))) &
          error stop 'Paired pressure modified composite coarse cell'
      enddo
      enddo
      enddo
      call compute_mol_eb_hydro_3d(species,extended(1),gas,temp,dt,2,rate,coarse_flux,boundary,wall_flux,good, &
        spatial_rate_halo=rate_halo,internal_energy_rate=erate,internal_energy_integrated_flux=held_face, &
        internal_energy_density=rhoe)
      if(.not.good) error stop 'Paired coarse spatial operator refused'
      step_rate%fields(1)%values=rate
      allocate(step_held_rate%fields(1)%values(1,n,n,n))
      step_held_rate%fields(1)%values(1,:,:,:)=erate
      step_flux(1)=coarse_flux;step_held_flux(1)=held_face
    endif
    if(audit_evolving) then
      evolution_rate=step_rate;evolution_held_rate=step_held_rate
      evolution_flux=step_flux;evolution_held_flux=step_held_flux
      return
    endif
    open(newunit=u,file=trim(prefix)//trim(suffix)//'_initial_flux.csv',status='replace')
    write(u,'(a)') 'node,axis,i,j,k,aperture,mass_kg_m2_s,normal_momentum_Pa,energy_W_m2'
    do p=1,2
      if(p==1) held_face=coarse_flux
      if(p==2) held_face=flux
      do a=1,3
        lo=1;lo(a)=0;hi=n
        do z=lo(3),hi(3)
        do y=lo(2),hi(2)
        do x=lo(1),hi(1)
          select case(a)
          case(1);f=held_face%x(:,x,y,z)/dt;aperture=physical(p)%x_face_fraction(x,y,z)
          case(2);f=held_face%y(:,x,y,z)/dt;aperture=physical(p)%y_face_fraction(x,y,z)
          case(3);f=held_face%z(:,x,y,z)/dt;aperture=physical(p)%z_face_fraction(x,y,z)
          end select
          write(u,'(5(i0,","),*(es25.16e3,:,","))') p,a,x,y,z,aperture,f(1),f(a+1),f(5)
        enddo
        enddo
        enddo
      enddo
    enddo
    close(u)
    if(paired) call audit_pressure_step(step_rate,step_held_rate,step_flux,step_held_flux,'_paired')
    write(*,'(a,i0,2(1x,es25.16e3))') 'Pressure ghost count/max energy correction/EOS residual: ', &
      changes,max_delta,max_pressure_error
  end subroutine audit_pressure_ghosts
  ! SSPRK2 diagnostic: rebuild donor pressures/halos at both stages, and
  ! synchronize each Euler map before its convex combination. Hydro only.
  ! Exact cell averages of a sinusoidal mole-fraction contact at fixed T,p,u.
  ! At fixed temperature all conserved densities are affine in mole fractions.
  ! Exact finite-volume averages of steady inviscid shear:
  ! rho,T,p,Y constant; u=10*sin(2*pi*y/L), v=w=0. The z wall is tangent.
  ! Total energy includes averaged u^2, not the square of averaged u.
  subroutine initialize_analytic_shear()
    type(eb_reactive_tree_3d) :: restricted
    type(amr_tree_state_3d) :: restricted_held
    real(dp) :: q(nv),gas(nv),t,cs,width,angle,mean_u,mean_u2,h,err
    integer :: p,x,y,z,c,n2
    logical :: good
    n2=0
    do c=1,size(species)
      if(trim(species(c)%name)=='N2') n2=c
    enddo
    if(n2==0) error stop 'Shear N2 unavailable'
    q=0;q(1)=101325._dp*species(n2)%molecular_weight/(universal_gas_constant*900._dp)
    q(5)=101325;q(5+n2)=1
    call reactive_primitive_to_conserved(species,q,gas,t,cs,good)
    if(.not.good) error stop 'Shear EOS refused'
    h=gas(5)
    do p=1,2
      width=2*acos(-1._dp)/(root_size*real(2**(p-1),dp))
      do y=1,root_size
        angle=(y-.5_dp)*width
        mean_u=10*sin(angle)*sin(width/2)/(width/2)
        mean_u2=50*(1-cos(2*angle)*sin(width)/width)
        gas(2)=gas(1)*mean_u;gas(5)=h+.5_dp*gas(1)*mean_u2
        do z=1,root_size
        do x=1,root_size
          eb%nodes(p)%state(:,x,y,z)=gas
          eb%nodes(p)%temperature(x,y,z)=t
          held%fields(p)%values(1,x,y,z)=h
          if(physical(p)%volume_fraction(x,y,z)<=0) then
            eb%nodes(p)%state(:,x,y,z)=0;held%fields(p)%values(1,x,y,z)=0
          endif
        enddo
        enddo
      enddo
    enddo
    restricted=eb;restricted_held=held
    call average_down_reactive_eb_tree_3d(species,restricted,good,internal_energy=restricted_held)
    if(.not.good) error stop 'Shear restriction refused'
    err=maxval(abs(restricted%nodes(1)%state-eb%nodes(1)%state)/max(1._dp,abs(eb%nodes(1)%state)))
    write(*,'(a,es24.16)') 'Analytic shear restriction residual: ',err
    if(err>1e-12_dp) error stop 'Shear analytic averages inconsistent'
  end subroutine initialize_analytic_shear

  subroutine write_analytic_norms(reference,reference_held,it,unit,label)
    type(eb_reactive_tree_3d), intent(in) :: reference
    type(amr_tree_state_3d), intent(in) :: reference_held
    integer, intent(in) :: it,unit
    character(len=*), intent(in) :: label
    real(dp) :: l1(nv+3),l2(nv+3),linf(nv+3),difference(nv+3),vol,total,q(nv),t,cs
    integer :: p,x,y,z,c,locations(4,nv+3),v,a,side,idx(3),cf(nv+3),cut(nv+3),wall(nv+3)
    real(dp) :: fractions(nv+3),signed_error(nv+3)
    logical :: near_cf
    logical :: good
    l1=0;l2=0;linf=0;total=0
    locations=0;fractions=0;signed_error=0;cf=0;cut=0;wall=0
    do p=1,2
    do z=1,root_size
    do y=1,root_size
    do x=1,root_size
      if(p==1.and.all([x,y,z]<=half_size)) cycle
      if(physical(p)%volume_fraction(x,y,z)<=0) cycle
      vol=physical(p)%volume_fraction(x,y,z)*product(spacing/real(2**(p-1),dp))
      call reactive_conserved_to_primitive(species,eb%nodes(p)%state(:,x,y,z), &
        eb%nodes(p)%temperature(x,y,z),q,t,cs,good,internal_energy_density=held%fields(p)%values(1,x,y,z))
      if(.not.good) error stop 'Shear norm EOS refused'
      difference(:nv)=eb%nodes(p)%state(:,x,y,z)-reference%nodes(p)%state(:,x,y,z)
      difference(nv+1)=held%fields(p)%values(1,x,y,z)-reference_held%fields(p)%values(1,x,y,z)
      difference(nv+2)=t-reference%nodes(p)%temperature(x,y,z)
      difference(nv+3)=q(5)-101325._dp
      near_cf=.false.
      if(p==2) then
        near_cf=x==1.or.x==root_size.or.y==1.or.y==root_size.or.z==root_size
      else
        do a=1,3
        do side=-1,1,2
          idx=[x,y,z];idx(a)=idx(a)+side
          if(a<=2) idx(a)=modulo(idx(a)-1,root_size)+1
          if(all(idx>=1).and.all(idx<=half_size)) near_cf=.true.
        enddo
        enddo
      endif
      do c=1,nv+3
        if(abs(difference(c))<=linf(c)) cycle
        locations(:,c)=[p,x,y,z];fractions(c)=physical(p)%volume_fraction(x,y,z)
        signed_error(c)=difference(c);cf(c)=merge(1,0,near_cf)
        cut(c)=merge(1,0,fractions(c)<1._dp)
        wall(c)=merge(1,0,(p==1.and.z==root_size).or.(plane<=0.and.z==1))
      enddo
      l1=l1+vol*abs(difference);l2=l2+vol*difference**2;linf=max(linf,abs(difference));total=total+vol
    enddo
    enddo
    enddo
    enddo
    do c=1,nv+3
      write(unit,'(i0,",",es25.16e3,",",i0,3(",",es25.16e3))') it,it*dt,c,l1(c)/total,sqrt(l2(c)/total),linf(c)
    enddo
    flush(unit)
    if(it==0) then
      open(newunit=v,file=trim(prefix)//'_'//label//'_locations.csv',status='replace')
      write(v,'(a)') 'step,time_s,component,node,i,j,k,vfrac,coarse_fine_adjacent,cut_cell,'// &
        'domain_z_wall_adjacent,signed_error'
    else
      open(newunit=v,file=trim(prefix)//'_'//label//'_locations.csv',status='old',position='append')
    endif
    do c=1,nv+3
      write(v,'(i0,",",es25.16e3,5(",",i0),",",es25.16e3,3(",",i0),",",es25.16e3)') &
        it,it*dt,c,locations(:,c),fractions(c),cf(c),cut(c),wall(c),signed_error(c)
    enddo
    close(v)
  end subroutine write_analytic_norms

  subroutine initialize_smooth_contact()
    type(eb_reactive_tree_3d) :: restricted
    type(amr_tree_state_3d) :: restricted_held
    real(dp) :: pi,width,angle,mole(lidryer_nspecies),mw,q(nv),gas(nv),t,cs,err
    integer :: p,x,y,z,h2,n2,c
    logical :: good
    if(.not.smooth_cut.and.(any(physical(1)%volume_fraction/=1).or.any(physical(2)%volume_fraction/=1))) &
      error stop 'Smooth contact requires regular geometry'
    h2=0;n2=0
    do c=1,size(species)
      if(trim(species(c)%name)=='H2') h2=c
      if(trim(species(c)%name)=='N2') n2=c
    enddo
    if(h2==0.or.n2==0) error stop 'Smooth contact species unavailable'
    pi=acos(-1._dp)
    do p=1,2
      width=2*pi/(root_size*real(2**(p-1),dp))
      do x=1,root_size
        angle=(x-.5_dp)*width
        mole=0;mole(h2)=.2_dp+.05_dp*cos(angle)*sin(width/2)/(width/2)
        mole(n2)=1-mole(h2);mw=sum(mole*species%molecular_weight)
        q=0;q(1)=101325._dp*mw/(universal_gas_constant*900._dp)
        q(2)=10;q(5)=101325;q(6:nv)=mole*species%molecular_weight/mw
        call reactive_primitive_to_conserved(species,q,gas,t,cs,good)
        if(.not.good) error stop 'Smooth contact EOS refused'
        do z=1,root_size
        do y=1,root_size
          if(physical(p)%volume_fraction(x,y,z)<=0) then
            eb%nodes(p)%state(:,x,y,z)=0;held%fields(p)%values(1,x,y,z)=0
            eb%nodes(p)%temperature(x,y,z)=t
            cycle
          endif
          eb%nodes(p)%state(:,x,y,z)=gas;eb%nodes(p)%temperature(x,y,z)=t
          held%fields(p)%values(1,x,y,z)=gas(5)-.5_dp*gas(1)*100._dp
        enddo
        enddo
      enddo
    enddo
    restricted=eb;restricted_held=held
    call average_down_reactive_eb_tree_3d(species,restricted,good,internal_energy=restricted_held)
    if(.not.good) error stop 'Smooth contact restriction refused'
    err=maxval(abs(restricted%nodes(1)%state-eb%nodes(1)%state)/ &
      max(abs(eb%nodes(1)%state),1e-14_dp))
    write(*,'(a,es24.16)') 'Smooth contact restriction residual: ',err
    write(*,*) 'Smooth contact temperature range: ',minval(eb%nodes(1)%temperature),maxval(eb%nodes(1)%temperature)
    if(err>1e-12_dp) error stop 'Smooth contact averages inconsistent'
  end subroutine initialize_smooth_contact

  ! Local convex mixing of two unrefined fluid cells, without an AMR operator.
  ! Report the pressure-preserving energy gap but never apply that correction.
  subroutine audit_conservative_mixing()
    real(dp) :: left(nv+2),right(nv+2),mixed(nv+2),q(nv),rebuilt(nv),t,cs
    real(dp) :: alpha,p_left,p_right,p_mix,gap,ke,target,base_alpha,tl,tr,qr(nv),original_left(nv+2)
    integer :: x,y,z,u,level,variant
    logical :: good
    open(newunit=u,file=trim(prefix)//'_mixing.csv',status='replace')
    write(u,'(a)') 'variant,i,j,k,fraction,p_left_Pa,p_right_Pa,p_mixed_Pa,linear_pressure_Pa,'// &
      'temperature_mixed_K,pressure_preserving_energy_gap_J_m3,temperature_left_K,temperature_right_K'
    do z=1,root_size
    do y=1,root_size
    do x=2,root_size
      if(all([x-1,y,z]<=half_size)) cycle
      if(physical(1)%volume_fraction(x,y,z)<=0.or.physical(1)%volume_fraction(x-1,y,z)<=0) cycle
      call pack(1,[x-1,y,z],left);call pack(1,[x,y,z],right)
      base_alpha=abs(right(2)/right(1))*dt/physical(1)%dx
      if(base_alpha<=0.or.base_alpha>=1) error stop 'Mixing audit requires positive subcell displacement'
      p_left=pressure_of(left);p_right=pressure_of(right)
      original_left=left
      call reactive_conserved_to_primitive(species,right(:nv),right(nv+1),qr,tr,cs,good, &
        internal_energy_density=right(nv+2))
      if(.not.good) error stop 'Mixing right donor EOS refused'
      do variant=0,2
        left=original_left
        call reactive_conserved_to_primitive(species,left(:nv),left(nv+1),q,tl,cs,good, &
          internal_energy_density=left(nv+2))
        if(.not.good) error stop 'Mixing left donor EOS refused'
        ! 0: imported donors; 1: equal composition; 2: equal temperature.
        ! Controls reconstruct only the left donor at its original pressure.
        if(variant>0) then
          if(variant==1) q(6:nv)=qr(6:nv)
          if(variant==2) tl=tr
          q(1)=p_left/(universal_gas_constant*sum(q(6:nv)/species%molecular_weight)*tl)
          call reactive_primitive_to_conserved(species,q,rebuilt,t,cs,good)
          if(.not.good) error stop 'Mixing control donor EOS refused'
          left=[rebuilt,t,rebuilt(5)-.5_dp*sum(rebuilt(2:4)**2)/rebuilt(1)]
        endif
      do level=0,2
        alpha=base_alpha/(2._dp**level)
        mixed=(1-alpha)*right+alpha*left
        call reactive_conserved_to_primitive(species,mixed(:nv),mixed(nv+1),q,t,cs,good, &
          internal_energy_density=mixed(nv+2))
        if(.not.good) error stop 'Mixed EOS refused'
        p_mix=q(5);target=(1-alpha)*p_right+alpha*p_left;q(5)=target
        call reactive_primitive_to_conserved(species,q,rebuilt,t,cs,good)
        if(.not.good) error stop 'Mixing target EOS refused'
        ke=.5_dp*sum(mixed(2:4)**2)/mixed(1)
        gap=rebuilt(5)-ke-mixed(nv+2)
        ! t here belongs to the pressure-preserving reconstruction; recover
        ! the actual mixed temperature for the observation column.
        call reactive_conserved_to_primitive(species,mixed(:nv),mixed(nv+1),q,t,cs,good, &
          internal_energy_density=mixed(nv+2))
        if(.not.good) error stop 'Mixed temperature refused'
        write(u,'(4(i0,","),*(es25.16e3,:,","))') variant,x,y,z,alpha,p_left,p_right,p_mix,target,t,gap,tl,tr
      enddo
      enddo
    enddo
    enddo
    enddo
    close(u)
  end subroutine audit_conservative_mixing

  subroutine audit_pressure_evolution()
    type(eb_reactive_tree_3d) :: initial_tree,old
    type(amr_tree_state_3d) :: initial_held,old_held
    real(dp) :: inventory(nv),h0
    integer :: mode,it,p,u,v,w
    logical :: good,paired
    character(len=16) :: label
    initial_tree=eb;initial_held=held
    call integrate_reactive_eb_tree_3d(eb,inventory,good)
    if(.not.good) error stop 'Evolution initial inventory refused'
    h0=held_inventory()
    do mode=1,2
      paired=mode==2;label='baseline'
      if(paired) label='paired'
      eb=initial_tree;held=initial_held;reset_energy_signed=0;reset_energy_absolute=0
      open(newunit=u,file=trim(prefix)//'_'//trim(label)//'_evolution.csv',status='replace')
      open(newunit=v,file=trim(prefix)//'_'//trim(label)//'_evolution_inventory.csv',status='replace')
      write(u,'(a)') 'step,time_s,max_conserved_relative_error,held_relative_error,'// &
        'max_pressure_offset_Pa,max_E_minus_H_minus_KE_J_m3,kinetic_energy_J,signed_E_H_KE_gap_J,L1_E_H_KE_gap_J,'// &
        'numerical_energy_correction_J,absolute_energy_correction_J,raw_total_energy_relative_change'
      write(v,'(a)') 'step,component,value'
      open(newunit=stage_unit,file=trim(prefix)//'_'//trim(label)//'_stages.csv',status='replace')
      write(stage_unit,'(a)') 'stage,node,i,j,k,pre_sync_pressure_offset_Pa,post_sync_pressure_offset_Pa'
      stage_counter=0
      open(newunit=w,file=trim(prefix)//'_'//trim(label)//'_evolution_cells.csv',status='replace')
      write(w,'(a)') 'step,node,i,j,k,vfrac,volume_m3,pressure_offset_Pa,E_minus_H_minus_KE_J_m3'
      call evolution_record(0,u,v,w,inventory,h0)
      do it=1,steps
        old=eb;old_held=held
        call evolution_euler(paired)
        call evolution_euler(paired)
        do p=1,2
          eb%nodes(p)%state=.5_dp*(old%nodes(p)%state+eb%nodes(p)%state)
          held%fields(p)%values=.5_dp*(old_held%fields(p)%values+held%fields(p)%values)
        enddo
        call refresh_evolution_temperature()
        if(audit_reset) call apply_evolution_reset()
        call evolution_record(it,u,v,w,inventory,h0)
      enddo
      close(u);close(v);close(w);close(stage_unit)
    enddo
    eb=initial_tree;held=initial_held
  end subroutine audit_pressure_evolution

  subroutine evolution_euler(paired)
    logical, intent(in) :: paired
    type(amr_tree_state_3d) :: rates,aux_rates,halos,aux_halos
    type(eb_integrated_flux_3d), allocatable :: flux(:),aux_flux(:)
    integer :: p,x,y,z
    real(dp) :: before(root_size,root_size,root_size,2),after(root_size,root_size,root_size,2)
    logical :: good
    call compute_mol_eb_tree_hydro_3d(species,eb,extended,dt,2,rates,flux,good, &
      boundary_faces=faces,old_state_halos=halos,internal_energy_rate=aux_rates, &
      internal_energy_flux=aux_flux,internal_energy_state=held, &
      internal_energy_old_state_halos=aux_halos,separate_momentum=audit_separate)
    if(.not.good) error stop 'Evolution spatial stage refused'
    if(smooth_contact) call audit_contact_halos(halos,aux_halos)
    if(paired) then
      call audit_pressure_ghosts(halos,aux_halos,flux,.true.)
      rates=evolution_rate;aux_rates=evolution_held_rate
      flux=evolution_flux;aux_flux=evolution_held_flux
    endif
    do p=1,2
      eb%nodes(p)%state=eb%nodes(p)%state+dt*rates%fields(p)%values
      held%fields(p)%values=held%fields(p)%values+dt*aux_rates%fields(p)%values
    enddo
    call stage_pressure(before)
    call synchronize_reactive_eb_tree_flux_3d(species,eb,flux,good,boundary_faces=faces)
    if(.not.good) error stop 'Evolution gas reflux refused'
    call synchronize_reactive_eb_tree_scalar_flux_3d(species,eb,held,aux_flux,good,boundary_faces=faces)
    if(.not.good) error stop 'Evolution held reflux refused'
    call stage_pressure(after)
    stage_counter=stage_counter+1
    do p=1,2
    do z=1,root_size
    do y=1,root_size
    do x=1,root_size
      if(p==1.and.all([x,y,z]<=half_size)) cycle
      if(physical(p)%volume_fraction(x,y,z)<=0) cycle
      write(stage_unit,'(5(i0,","),*(es25.16e3,:,","))') stage_counter,p,x,y,z, &
        before(x,y,z,p),after(x,y,z,p)
    enddo
    enddo
    enddo
    enddo
    call refresh_evolution_temperature()
  end subroutine evolution_euler

  subroutine audit_contact_halos(gas,aux)
    type(amr_tree_state_3d), intent(in) :: gas,aux
    real(dp) :: q(nv),t,cs,errors(3),peak(3)
    integer :: p,x,y,z,loc(3,3),c
    logical :: good
    do p=1,2
      peak=0;loc=0
      do z=1,root_size+12
      do y=1,root_size+12
      do x=1,root_size+12
        if(extended(p)%volume_fraction(x,y,z)<=0) cycle
        call reactive_conserved_to_primitive(species,gas%fields(p)%values(:,x,y,z),900._dp, &
          q,t,cs,good,internal_energy_density=aux%fields(p)%values(1,x,y,z))
        if(.not.good) error stop 'Contact halo EOS refused'
        errors=[abs(q(5)-101325._dp),abs(t-900._dp),maxval(abs(q(2:4)-[10._dp,0._dp,0._dp]))]
        do c=1,3
          if(errors(c)>peak(c)) then
            peak(c)=errors(c);loc(:,c)=[x,y,z]-6
          endif
        enddo
      enddo
      enddo
      enddo
      write(*,*) 'CONTACT_HALO stage,node,p,T,u,locations ',stage_counter+1,p,peak,loc
    enddo
  end subroutine audit_contact_halos

  subroutine stage_pressure(offset)
    real(dp), intent(out) :: offset(root_size,root_size,root_size,2)
    real(dp) :: q(nv),t,cs
    integer :: p,x,y,z
    logical :: good
    offset=0
    do p=1,2
    do z=1,root_size
    do y=1,root_size
    do x=1,root_size
      if(p==1.and.all([x,y,z]<=half_size)) cycle
      if(physical(p)%volume_fraction(x,y,z)<=0) cycle
      call reactive_conserved_to_primitive(species,eb%nodes(p)%state(:,x,y,z), &
        eb%nodes(p)%temperature(x,y,z),q,t,cs,good,internal_energy_density=held%fields(p)%values(1,x,y,z))
      if(.not.good) error stop 'Stage pressure EOS refused'
      offset(x,y,z,p)=q(5)-101325._dp
    enddo
    enddo
    enddo
    enddo
  end subroutine stage_pressure

  ! Apply the existing cell policy transactionally on accepted-step copies.
  ! This is not the full prepared held-MOL/re-reflux closure.
  subroutine apply_evolution_reset()
    type(eb_reactive_tree_3d) :: candidate
    type(amr_tree_state_3d) :: stored
    real(dp) :: delta,volume,change,absolute_change,previous_held
    integer :: p,x,y,z
    logical :: good
    candidate=eb;stored=held;change=0;absolute_change=0
    do p=1,2
    do z=1,root_size
    do y=1,root_size
    do x=1,root_size
      if(p==1.and.all([x,y,z]<=half_size)) cycle
      if(physical(p)%volume_fraction(x,y,z)<=0) cycle
      previous_held=stored%fields(p)%values(1,x,y,z);delta=0
      call reset_nasa7_dual_energy_cell(species,candidate%nodes(p)%state(:,x,y,z), &
        stored%fields(p)%values(1,x,y,z),candidate%nodes(p)%temperature(x,y,z),policy,delta,good)
      if(.not.good) error stop 'Evolution cell reset refused'
      if(stored%fields(p)%values(1,x,y,z)/=previous_held) error stop 'Unexpected diagnostic held reset'
      volume=physical(p)%dx*physical(p)%dy*physical(p)%dz*physical(p)%volume_fraction(x,y,z)
      change=change+volume*delta;absolute_change=absolute_change+volume*abs(delta)
    enddo
    enddo
    enddo
    enddo
    call average_down_reactive_eb_tree_3d(species,candidate,good,internal_energy=stored)
    if(.not.good) error stop 'Evolution reset restriction refused'
    eb=candidate;held=stored
    reset_energy_signed=reset_energy_signed+change;reset_energy_absolute=reset_energy_absolute+absolute_change
    call refresh_evolution_temperature()
  end subroutine apply_evolution_reset

  subroutine refresh_evolution_temperature()
    integer :: p,x,y,z
    real(dp) :: q(nv),t,cs
    logical :: good
    do p=1,2
    do z=1,root_size
    do y=1,root_size
    do x=1,root_size
      if(physical(p)%volume_fraction(x,y,z)<=0) cycle
      call reactive_conserved_to_primitive(species,eb%nodes(p)%state(:,x,y,z), &
        eb%nodes(p)%temperature(x,y,z),q,t,cs,good,internal_energy_density=held%fields(p)%values(1,x,y,z))
      if(.not.good) error stop 'Evolution EOS refused'
      eb%nodes(p)%temperature(x,y,z)=t
    enddo
    enddo
    enddo
    enddo
  end subroutine refresh_evolution_temperature

  real(dp) function held_inventory() result(total)
    integer :: p,x,y,z
    total=0
    do p=1,2
    do z=1,root_size
    do y=1,root_size
    do x=1,root_size
      if(p==1.and.all([x,y,z]<=half_size)) cycle
      total=total+physical(p)%dx*physical(p)%dy*physical(p)%dz*physical(p)%volume_fraction(x,y,z)* &
        held%fields(p)%values(1,x,y,z)
    enddo
    enddo
    enddo
    enddo
  end function held_inventory

  subroutine evolution_record(it,u,v,w,initial,h0)
    integer, intent(in) :: it,u,v,w
    real(dp), intent(in) :: initial(nv),h0
    real(dp) :: total(nv),h,error,max_p,max_gap,q(nv),t,cs,ke,volume,kinetic,gap,signed_gap,l1_gap,residual
    integer :: p,x,y,z,c
    logical :: good
    call integrate_reactive_eb_tree_3d(eb,total,good)
    if(.not.good) error stop 'Evolution inventory refused'
    h=held_inventory();error=0
    do c=1,nv
      if(c>=2.and.c<=4) cycle
      residual=total(c)-initial(c)
      if(c==5) residual=residual-reset_energy_signed
      error=max(error,abs(residual)/max(abs(initial(c)),abs(initial(1))*1e-14_dp))
    enddo
    if(error>1e-12_dp.or.abs(h-h0)/abs(h0)>1e-12_dp) error stop 'Evolution budget failed'
    max_p=0;max_gap=0;kinetic=0;signed_gap=0;l1_gap=0
    do p=1,2
    do z=1,root_size
    do y=1,root_size
    do x=1,root_size
      if(p==1.and.all([x,y,z]<=half_size)) cycle
      if(physical(p)%volume_fraction(x,y,z)<=0) cycle
      call reactive_conserved_to_primitive(species,eb%nodes(p)%state(:,x,y,z),eb%nodes(p)%temperature(x,y,z), &
        q,t,cs,good,internal_energy_density=held%fields(p)%values(1,x,y,z))
      if(.not.good) error stop 'Evolution record EOS refused'
      max_p=max(max_p,abs(q(5)-101325._dp))
      ke=.5_dp*sum(eb%nodes(p)%state(2:4,x,y,z)**2)/eb%nodes(p)%state(1,x,y,z)
      gap=eb%nodes(p)%state(5,x,y,z)-ke-held%fields(p)%values(1,x,y,z)
      max_gap=max(max_gap,abs(gap))
      volume=physical(p)%dx*physical(p)%dy*physical(p)%dz*physical(p)%volume_fraction(x,y,z)
      write(w,'(5(i0,","),*(es25.16e3,:,","))') it,p,x,y,z,physical(p)%volume_fraction(x,y,z), &
        volume,q(5)-101325._dp,gap
      kinetic=kinetic+volume*ke;signed_gap=signed_gap+volume*gap;l1_gap=l1_gap+volume*abs(gap)
    enddo
    enddo
    enddo
    enddo
    write(u,'(i0,11(",",es25.16e3))') it,it*dt,error,abs(h-h0)/abs(h0),max_p,max_gap,kinetic,signed_gap,l1_gap, &
      reset_energy_signed,reset_energy_absolute,(total(5)-initial(5))/abs(initial(5))
    do c=1,nv
      write(v,'(i0,",",i0,",",es25.16e3)') it,c,total(c)
    enddo
    write(v,'(i0,",held,",es25.16e3)') it,h
    flush(u);flush(v)
  end subroutine evolution_record

  ! Forward-Euler hydro-only increment with production gas/scalar reflux.
  ! This diagnostic does not claim SDC/chemistry or accepted-app equivalence.
  subroutine audit_pressure_step(rates,held_rates,gas_faces,held_faces,tag)
    type(amr_tree_state_3d), intent(in) :: rates,held_rates
    character(len=*), intent(in) :: tag
    type(eb_integrated_flux_3d), intent(in) :: gas_faces(2),held_faces(2)
    type(eb_reactive_tree_3d) :: updated
    type(amr_tree_state_3d) :: updated_held
    type(eb_integrated_flux_3d) :: flux(2),aux_flux(2)
    real(dp) :: initial(nv),final(nv),h0,h1,volume,ke,q(nv),t,cs,pgas,qgas(nv)
    real(dp) :: max_pressure,max_gap,max_update,max_eos_gap,component_error,held_error
    integer :: p,x,y,z,u,c
    logical :: good
    updated=eb;updated_held=held;flux=gas_faces;aux_flux=held_faces
    call integrate_reactive_eb_tree_3d(eb,initial,good)
    if(.not.good) error stop 'Step initial inventory refused'
    do p=1,2
      updated%nodes(p)%state=eb%nodes(p)%state+dt*rates%fields(p)%values
      updated_held%fields(p)%values=held%fields(p)%values+dt*held_rates%fields(p)%values
    enddo
    call synchronize_reactive_eb_tree_flux_3d(species,updated,flux,good,boundary_faces=faces)
    if(.not.good) error stop 'Step gas reflux refused'
    call synchronize_reactive_eb_tree_scalar_flux_3d(species,updated,updated_held,aux_flux,good,boundary_faces=faces)
    if(.not.good) error stop 'Step held reflux refused'
    call integrate_reactive_eb_tree_3d(updated,final,good)
    if(.not.good) error stop 'Step final inventory refused'
    h0=0;h1=0;max_pressure=0;max_gap=0;max_update=0;max_eos_gap=0
    do p=1,2
    do z=1,root_size
    do y=1,root_size
    do x=1,root_size
      if(p==1.and.all([x,y,z]<=half_size)) cycle
      if(physical(p)%volume_fraction(x,y,z)<=0) cycle
      volume=physical(p)%volume_fraction(x,y,z)*physical(p)%dx*physical(p)%dy*physical(p)%dz
      h0=h0+volume*held%fields(p)%values(1,x,y,z)
      h1=h1+volume*updated_held%fields(p)%values(1,x,y,z)
      call reactive_conserved_to_primitive(species,updated%nodes(p)%state(:,x,y,z), &
        eb%nodes(p)%temperature(x,y,z),q,t,cs,good,internal_energy_density=updated_held%fields(p)%values(1,x,y,z))
      if(.not.good) error stop 'Step held EOS refused'
      max_pressure=max(max_pressure,abs(q(5)-101325._dp))
      call reactive_conserved_to_primitive(species,updated%nodes(p)%state(:,x,y,z), &
        eb%nodes(p)%temperature(x,y,z),qgas,t,cs,good)
      if(.not.good) error stop 'Step total-energy EOS refused'
      pgas=qgas(5);max_eos_gap=max(max_eos_gap,abs(pgas-q(5)))
      ke=.5_dp*sum(updated%nodes(p)%state(2:4,x,y,z)**2)/updated%nodes(p)%state(1,x,y,z)
      max_gap=max(max_gap,abs(updated%nodes(p)%state(5,x,y,z)-ke-updated_held%fields(p)%values(1,x,y,z)))
      max_update=max(max_update,maxval(abs(updated%nodes(p)%state(:,x,y,z)-eb%nodes(p)%state(:,x,y,z))))
    enddo
    enddo
    enddo
    enddo
    ! Closed periodic/slip-wall hydro conserves mass, energy and each species.
    ! Wall-normal momentum has its own wall-force budget and is excluded.
    component_error=max(abs(final(1)-initial(1))/max(abs(initial(1)),tiny(1._dp)), &
      abs(final(5)-initial(5))/max(abs(initial(5)),tiny(1._dp)))
    do c=6,nv
      component_error=max(component_error,abs(final(c)-initial(c))/max(abs(initial(c)),abs(initial(1))*1e-14_dp))
    enddo
    held_error=abs(h1-h0)/max(abs(h0),tiny(1._dp))
    if(component_error>1e-12_dp.or.held_error>1e-12_dp) error stop 'Step composite budget failed'
    open(newunit=u,file=trim(prefix)//tag//'_pressure_step.csv',status='replace')
    write(u,'(a)') 'dt_s,max_conserved_relative_error,held_relative_error,max_pressure_offset_Pa,'// &
      'max_E_minus_H_minus_KE_J_m3,max_state_change,max_total_held_pressure_gap_Pa'
    write(u,'(*(es25.16e3,:,","))') dt,component_error,held_error,max_pressure,max_gap,max_update,max_eos_gap
    close(u)
    open(newunit=u,file=trim(prefix)//tag//'_pressure_step_inventory.csv',status='replace')
    write(u,'(a)') 'component,initial,final'
    do c=1,nv
      write(u,'(i0,2(",",es25.16e3))') c,initial(c),final(c)
    enddo
    write(u,'(a,2(",",es25.16e3))') 'held',h0,h1
    close(u)
  end subroutine audit_pressure_step

  real(dp) function donor_pressure(p,idx,paired) result(pressure)
    integer, intent(in) :: p,idx(3)
    logical, intent(in) :: paired
    real(dp) :: value(nv+2),weight,total
    integer :: x,y,z
    call pack(p,idx,value);pressure=pressure_of(value)
    if(.not.paired.or.p/=1.or.any(idx>half_size)) return
    pressure=0;total=0
    do z=2*idx(3)-1,2*idx(3)
    do y=2*idx(2)-1,2*idx(2)
    do x=2*idx(1)-1,2*idx(1)
      weight=physical(2)%volume_fraction(x,y,z)
      if(weight<=0) cycle
      call pack(2,[x,y,z],value)
      pressure=pressure+weight*pressure_of(value);total=total+weight
    enddo
    enddo
    enddo
    if(total<=0) error stop 'Paired pressure empty covered donor'
    pressure=pressure/total
  end function donor_pressure
  ! Independent reconstruction using the same public planner/kernel. No state
  ! mutation: pressure envelope and child-average checks observe conserved data.
  subroutine audit_initial_stencils(halos,aux)
    type(amr_tree_state_3d), intent(in) :: halos,aux
    integer(int64) :: origins(3,2),domains(3,2),global(3)
    integer :: shapes(3,2),requests(4,6),signs(3,6),idx(3),src(3)
    integer :: u,v,x,y,z,a,b,d,q,ii,jj,kk,c
    real(dp) :: center(nv+2),neighbors(nv+2,6),candidate(nv+2),mean(nv+2),child(nv+2)
    real(dp) :: default_value(nv+2),separate_value(nv+2)
    real(dp) :: offset(3),child_offset(3),pc,pn(6),ph,ps,err,conservation,w,total
    real(dp) :: restricted(nv+2),volume,pmin,pmax,eos_center(nv+2),eos_neighbors(nv+2,6),eos_p
    real(dp) :: envelope_theta,envelope_p
    logical :: good,linear
    origins=0;domains(:,1)=root_size;domains(:,2)=2*root_size;shapes=root_size
    open(newunit=v,file=trim(prefix)//'_initial_limiter.csv',status='replace')
    write(v,'(a)') 'i,j,k,component,center,xminus,xplus,yminus,yplus,zminus,zplus,'// &
      'offset_x,offset_y,offset_z,default_value,separate_value'
    open(newunit=u,file=trim(prefix)//'_initial_stencils.csv',status='replace')
    write(u,'(a)') 'i,j,k,donor_i,donor_j,donor_k,linear,donor_p,neighbor_p_min,neighbor_p_max,'// &
      'ghost_p,stored_T_pressure,halo_scaled_residual,child_mean_scaled_residual,pressure_envelope_excess_Pa,'// &
      'eos_T_ghost_p,envelope_theta,envelope_ghost_p'
    do z=0,root_size+1
    do y=0,root_size+1
    do x=0,root_size+1
      idx=[x,y,z]
      if(count(idx<1.or.idx>root_size)/=1.or.z==0) cycle
      if(extended(2)%volume_fraction(x+6,y+6,z+6)<=0) cycle
      global=int(idx-1,int64)
      global(1:2)=modulo(global(1:2),domains(1:2,2))
      src=int(global/2)+1
      call pack(1,src,center)
      call plan_eb_linear_halo_3d([0,1],shapes,origins,domains,1,src,2,global, &
        extended(1),6,faces,requests,signs,offset,linear,good)
      if(.not.good) error stop 'Audit stencil planner refused'
      candidate=center;pc=pressure_of(center);pn=pc;mean=center
      if(linear) then
        do q=1,6
          call pack(requests(1,q),requests(2:4,q),neighbors(:,q))
          neighbors(2:4,q)=neighbors(2:4,q)*real(signs(:,q),dp)
          pn(q)=pressure_of(neighbors(:,q))
        enddo
        default_value=center;separate_value=center
        call interpolate_eb_linear_halo_3d(default_value,neighbors,offset,good)
        if(.not.good) error stop 'Audit default limiter refused'
        call interpolate_eb_linear_halo_3d(separate_value,neighbors,offset,good,separate_momentum=.true.)
        if(.not.good) error stop 'Audit separate limiter refused'
        do c=1,nv+2
          write(v,'(4(i0,","),*(es25.16e3,:,","))') idx,c,center(c),neighbors(c,:),offset, &
            default_value(c),separate_value(c)
        enddo
        call interpolate_eb_linear_halo_3d(candidate,neighbors,offset,good,audit_separate)
        if(.not.good) error stop 'Audit stencil interpolation refused'
        mean=0
        do d=0,1
        do b=0,1
        do a=0,1
          child=center;child_offset=.5_dp*real([a,b,d],dp)-.25_dp
          call interpolate_eb_linear_halo_3d(child,neighbors,child_offset,good,audit_separate)
          if(.not.good) error stop 'Audit child refused'
          mean=mean+child/8
        enddo
        enddo
        enddo
      endif
      eos_center=center
      call eos_temperature(eos_center)
      if(linear) then
        eos_neighbors=neighbors
        do q=1,6
          call eos_temperature(eos_neighbors(:,q))
        enddo
        call interpolate_eb_linear_halo_3d(eos_center,eos_neighbors,offset,good,audit_separate)
        if(.not.good) error stop 'Audit EOS-temperature stencil refused'
      endif
      eos_p=pressure_of(eos_center)
      ph=pressure_of(candidate)
      ps=mixture_pressure(species,candidate(6:nv)/candidate(1),candidate(1),candidate(nv+1),good)
      if(.not.good) error stop 'Audit stored-temperature pressure refused'
      err=maxval(abs(candidate(:nv)-halos%fields(2)%values(:,x+6,y+6,z+6))/max(1._dp,abs(candidate(:nv))))
      err=max(err,abs(candidate(nv+2)-aux%fields(2)%values(1,x+6,y+6,z+6))/max(1._dp,abs(candidate(nv+2))))
      conservation=maxval(abs(mean-center)/max(1._dp,abs(center)))
      pmin=min(pc,minval(pn));pmax=max(pc,maxval(pn))
      envelope_theta=1
      if(linear) call envelope_trial(center,neighbors,pmin,pmax,envelope_theta)
      envelope_p=pressure_of(center+envelope_theta*(candidate-center))
      write(u,'(7(i0,","),*(es25.16e3,:,","))') idx,src,merge(1,0,linear),pc,minval(pn),maxval(pn), &
        ph,ps,err,conservation,max(0._dp,ph-pmax,pmin-ph),eos_p,envelope_theta,envelope_p
    enddo
    enddo
    enddo
    close(u)
    close(v)
    open(newunit=u,file=trim(prefix)//'_initial_restriction.csv',status='replace')
    write(u,'(a)') 'i,j,k,coarse_p,restricted_p,fine_p_min,fine_p_max,gas_scaled_residual,held_scaled_residual'
    do z=1,half_size
    do y=1,half_size
    do x=1,half_size
      if(physical(1)%volume_fraction(x,y,z)<=0) cycle
      restricted=0;total=0;pmin=huge(1._dp);pmax=-huge(1._dp)
      do kk=2*z-1,2*z
      do jj=2*y-1,2*y
      do ii=2*x-1,2*x
        volume=physical(2)%volume_fraction(ii,jj,kk)
        if(volume<=0) cycle
        call pack(2,[ii,jj,kk],child)
        restricted=restricted+volume*child;total=total+volume
        ph=pressure_of(child);pmin=min(pmin,ph);pmax=max(pmax,ph)
      enddo
      enddo
      enddo
      if(total<=0) error stop 'Audit empty restriction'
      restricted=restricted/total
      call pack(1,[x,y,z],center)
      err=maxval(abs(restricted(:nv)-center(:nv))/max(1._dp,abs(center(:nv))))
      w=abs(restricted(nv+2)-center(nv+2))/max(1._dp,abs(center(nv+2)))
      write(u,'(3(i0,","),*(es25.16e3,:,","))') x,y,z,pressure_of(center),pressure_of(restricted),pmin,pmax,err,w
    enddo
    enddo
    enddo
    close(u)
  end subroutine audit_initial_stencils

  subroutine pack(p,idx,value)
    integer, intent(in) :: p,idx(3)
    real(dp), intent(out) :: value(nv+2)
    value=[eb%nodes(p)%state(:,idx(1),idx(2),idx(3)),eb%nodes(p)%temperature(idx(1),idx(2),idx(3)), &
      held%fields(p)%values(1,idx(1),idx(2),idx(3))]
  end subroutine pack

  ! Diagnostic only: one factor for all eight children's conserved slopes.
  ! This preserves their parent mean, but may collapse accuracy to PCM.
  ! Record that cost before considering any production numerical policy.
  subroutine envelope_trial(center,neighbors,pmin,pmax,theta)
    real(dp), intent(in) :: center(nv+2),neighbors(nv+2,6),pmin,pmax
    real(dp), intent(out) :: theta
    real(dp) :: delta(nv+2,8),child(nv+2),shift(3),lo,hi,trial,p,tolerance
    integer :: a,b,c,q,it,j
    logical :: good,admitted
    q=0
    do c=0,1
    do b=0,1
    do a=0,1
      q=q+1;child=center;shift=.5_dp*real([a,b,c],dp)-.25_dp
      call interpolate_eb_linear_halo_3d(child,neighbors,shift,good,audit_separate)
      if(.not.good) error stop 'Envelope trial interpolation refused'
      delta(:,q)=child-center
    enddo
    enddo
    enddo
    tolerance=128*epsilon(1._dp)*max(1._dp,abs(pmin),abs(pmax))
    lo=0;hi=1;trial=1
    do it=1,42
      admitted=.true.
      do j=1,8
        p=pressure_of(center+trial*delta(:,j))
        if(p<pmin-tolerance.or.p>pmax+tolerance) admitted=.false.
      enddo
      if(admitted) then
        lo=trial
        if(it==1) exit
      else
        hi=trial
      endif
      trial=.5_dp*(lo+hi)
    enddo
    theta=lo
  end subroutine envelope_trial

  subroutine eos_temperature(value)
    real(dp), intent(inout) :: value(nv+2)
    real(dp) :: q(reactive_nprim(size(species))),t,cs
    logical :: good
    call reactive_conserved_to_primitive(species,value(:nv),value(nv+1),q,t,cs,good, &
      internal_energy_density=value(nv+2))
    if(.not.good) error stop 'Audit EOS-temperature recovery refused'
    value(nv+1)=t
  end subroutine eos_temperature

  real(dp) function pressure_of(value) result(p)
    real(dp), intent(in) :: value(nv+2)
    real(dp) :: q(reactive_nprim(size(species))),t,cs
    logical :: good
    call reactive_conserved_to_primitive(species,value(:nv),value(nv+1),q,t,cs,good, &
      internal_energy_density=value(nv+2))
    if(.not.good) error stop 'Audit stencil EOS refused'
    p=q(5)
  end function pressure_of

  ! This probe loads one immutable mechanism. Cache only its initialization
  ! callbacks, using every floating-point input bit and the step-limit option.
  ! No quantization, approximate matching, failed solves or evolved states.
  subroutine cached_initial_reaction(records,rxns,energy,forcing,energy_forcing,duration,rtol,atol, &
      rho_y,temp,completed_steps,accepted,maximum_steps)
    type(nasa7_species), intent(in) :: records(:)
    type(elementary_reaction), intent(in) :: rxns(:)
    real(dp), intent(in) :: energy,forcing(:),energy_forcing,duration,rtol,atol
    real(dp), intent(inout) :: rho_y(:),temp
    integer, intent(out) :: completed_steps
    logical, intent(out) :: accepted
    integer, optional, intent(in) :: maximum_steps
    integer(int64) :: key(2*lidryer_nspecies+6)
    integer :: slot,limit
    limit=-1
    if(present(maximum_steps)) limit=maximum_steps
    key=transfer([energy,energy_forcing,duration,rtol,merge(reaction_atol,atol,density_absolute_tolerance), &
      temp,rho_y,forcing],key)
    do slot=1,initial_entries
      if(initial_limits(slot)/=limit) cycle
      if(.not.all(initial_keys(:,slot)==key)) cycle
      rho_y=initial_results(:lidryer_nspecies,slot);temp=initial_results(lidryer_nspecies+1,slot)
      completed_steps=initial_steps(slot);accepted=.true.;initial_hits=initial_hits+1
      return
    enddo
    initial_solves=initial_solves+1
    call controlled_reaction(records,rxns,energy,forcing,energy_forcing,duration,rtol,atol, &
      rho_y,temp,completed_steps,accepted,maximum_steps)
    if(.not.accepted.or.initial_entries==initial_cache_capacity) return
    initial_entries=initial_entries+1;slot=initial_entries
    initial_keys(:,slot)=key;initial_limits(slot)=limit;initial_steps(slot)=completed_steps
    initial_results(:,slot)=[rho_y,temp]
  end subroutine cached_initial_reaction
  subroutine controlled_reaction(records,rxns,energy,forcing,energy_forcing,duration,rtol,atol, &
      rho_y,temp,completed_steps,accepted,maximum_steps)
    type(nasa7_species), intent(in) :: records(:)
    type(elementary_reaction), intent(in) :: rxns(:)
    real(dp), intent(in) :: energy,forcing(:),energy_forcing,duration,rtol,atol
    real(dp), intent(inout) :: rho_y(:),temp
    integer, intent(out) :: completed_steps
    logical, intent(out) :: accepted
    integer, optional, intent(in) :: maximum_steps
    real(dp) :: effective_atol
    ! react_sdc_state_3d multiplies its mass-fraction tolerance by density.
    ! V3 explicitly requests a density tolerance, so bypass that conversion
    ! here. V1/V2 pass through the historical, already density-scaled value.
    effective_atol=merge(reaction_atol,atol,density_absolute_tolerance)
    cvode_atol_min=min(cvode_atol_min,effective_atol);cvode_atol_max=max(cvode_atol_max,effective_atol)
    call advance_forced_reactor_cvode(records,rxns,energy,forcing,energy_forcing,duration,rtol, &
      effective_atol,rho_y,temp,completed_steps,accepted,maximum_steps)
  end subroutine controlled_reaction
end program eb_tree_reacting_probe
