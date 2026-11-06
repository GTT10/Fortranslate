program test_mpi_eb_tree_full_3d
  use iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_finite
  use mpi_f08
  use precision_mod, only: dp
  use eb_held_funding_queue_3d_mod, only: eb_held_funding_context_3d
  use dual_energy_mod, only: dual_energy_policy,reset_nasa7_dual_energy_cell
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval,forced_reactor_integrator
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: temperature_from_internal_energy
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport,load_h2o2_full_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,average_down_reactive_eb_tree_field_3d, &
    integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_full_3d,react_sdc_eb_tree_3d
  use eb_reactive_patch_tree_3d_mod, only: advance_reactive_eb_tree_sdc_3d,build_eb_tree_extended_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: compute_sdc_eb_tree_transport_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_full_3d,react_sdc_mpi_eb_tree_3d,advance_mpi_eb_tree_sdc_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_subcycle_3d,mpi_eb_sdc_recursive_source, &
    mpi_eb_mol_interval_3d,reconstruct_mpi_eb_mol_interval_3d,reconstruct_mpi_eb_mol_state_redist_interval_3d
  use eb_reactive_patch_tree_3d_mod, only: advance_reactive_eb_tree_sdc_subcycle_3d,eb_sdc_recursive_source, &
    eb_mol_held_interval_3d,replay_reactive_eb_mol_held_group_3d,reconstruct_reactive_eb_mol_held_registration_3d
  use mpi_eb_tree_hydro_3d_mod, only: compute_mol_mpi_eb_tree_hydro_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_sdc_mpi_eb_tree_source_halo_3d
  use eb_reactive_transport_3d_mod, only: compute_sdc_eb_transport_halo_3d
  use mpi_eb_tree_transport_3d_mod, only: compute_sdc_mpi_eb_tree_transport_3d, &
    compute_sdc_mpi_eb_tree_transport_halos_3d
  use eb_reactive_redistribution_3d_mod, only: state_redist_register_3d,post_eb_spatial_density_rate_3d, &
    accumulate_state_redist_mol_registration_stage_3d, &
    initialize_state_redist_tree_register_3d,accumulate_state_redist_register_3d, &
    accumulate_state_redist_coarse_register_3d
  use mpi_eb_tree_flux_3d_mod, only: synchronize_mpi_eb_tree_register_flux_3d,assemble_mpi_eb_tree_reflux_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d
  use mpi_eb_tree_redistribution_3d_mod, only: compute_state_redist_mpi_level_stage_3d, &
    assemble_state_redist_mpi_interface_rereflux_3d
  implicit none
  logical :: scalar_replay_budget_test=.false.,native_funding_replay_test=.false.
  logical :: native_material_replay_test=.false.
  logical :: material_refusal_probe=.false.
  real(dp), allocatable :: material_packet(:)
  real(dp), parameter :: material_rate=1.8e8_dp
  real(dp), parameter :: material_refusal_rate=2.2e8_dp
  integer, parameter :: material_target(3)=[3,3,2]
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),base_mole_fractions(:)
  real(dp), allocatable :: before(:),after(:),boundary(:),primitive(:),wall(:),chemical(:)
  real(dp) :: residual,theta,speed
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k,mode,integrator_index,rank,nranks,status
  integer :: serial_source_calls=0,mpi_source_calls=0,source_failure=0
  integer :: prepared_parent_refusal_call=0
  integer, allocatable :: diffusion_owners(:)
  type(eb_geometry_3d), allocatable :: diffusion_geometry(:),spatial_geometry(:)
  type(amr_patch_tree_3d) :: spatial_topology
  integer :: spatial_failure=0,spatial_calls=0
  logical :: diffusion_failure=.false.
  integer :: history_post_calls=0,history_recursive_calls=0
  integer :: history_calls=0
  integer :: history_failure_depth=0,history_failure_iteration=0
  logical :: history_failure=.false.
  integer :: reactor_calls=0,reactor_failure_rank=0
  logical :: reactor_failure=.false.
  logical :: reactor_replay_failure=.false.
  logical :: replay_refusal_only=.false.
  logical :: reactor_replay_refused=.false.
  real(dp) :: reactor_replay_duration=0
  integer :: reactor_replay_calls=0,reactor_replay_refusal_call=0,reactor_replay_refusal_actual=0
  integer :: diffusion_calls=0
  logical :: held_diffusion_failure=.false.
  integer :: auxiliary_family_calls=0,auxiliary_family_failure=0
  logical :: auxiliary_family_water=.false.,auxiliary_family_reacting=.false.
  integer, allocatable :: auxiliary_family_owners(:)
  type(eb_geometry_3d), allocatable :: auxiliary_family_geometry(:)
  integer :: predictor_calls=0
  integer, allocatable :: predictor_owners(:)
  type(eb_geometry_3d), allocatable :: predictor_geometry(:)
  type(amr_tree_state_3d) :: predictor_metrics
  character(len=32) :: test_mode,order_argument
  integer :: requested_order=0,order_status
  logical :: ok,chemistry_on,transport_on
  character(len=8), parameter :: integrators(2)=[character(len=8)::'explicit','implicit']
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  call get_command_argument(1,test_mode)
  native_material_replay_test=test_mode=='held_native_material'
  native_funding_replay_test=test_mode=='held_native_funding'.or.native_material_replay_test
  if(native_funding_replay_test) test_mode='held_natural_replay'
  scalar_replay_budget_test=test_mode=='held_scalar_replay'.or.test_mode=='held_prepared_scalar_replay'
  if(test_mode=='held_scalar_replay') test_mode='held_natural_replay'
  if(test_mode=='held_prepared_scalar_replay') test_mode='held_prepared_natural'
  if(test_mode=='flux_refusal'.or.test_mode=='routed_refusal'.or.test_mode=='stored_reactor'.or. &
    test_mode=='held_policy'.or.test_mode=='held_chem_policy'.or.test_mode=='held_recursive_policy'.or. &
    test_mode=='held_prepared'.or.test_mode=='held_prepared_recursive'.or. &
    test_mode=='held_recursive_chemistry'.or.test_mode=='held_prepared_serial'.or. &
    test_mode=='held_natural_replay'.or.test_mode=='held_prepared_parent'.or.test_mode=='held_prepared_natural') &
    call require(nranks>1,'collective refusal fixture requires multiple MPI ranks; check launcher ABI')
  if(test_mode=='flux_refusal') then
    test_mode='flux_replay';replay_refusal_only=.true.
  else if(test_mode=='routed_refusal') then
    test_mode='routed_replay';replay_refusal_only=.true.
  endif
  call get_command_argument(2,order_argument)
  material_refusal_probe=native_material_replay_test.and.order_argument=='refusal-probe'
  if(len_trim(order_argument)>0.and..not.material_refusal_probe) then
    read(order_argument,*,iostat=order_status) requested_order
    call require(order_status==0,'reconstruction order argument parses')
    call require(requested_order==0.or.requested_order==2.or.requested_order==4,'supported reconstruction order argument')
    call require(test_mode=='state_redist'.or.test_mode=='state_recursive','high-order fixture mode')
  endif
  if(test_mode=='routed_replay'.or.test_mode=='routed_predictor'.or. &
     test_mode=='flux_replay'.or.test_mode=='state_replay'.or.test_mode=='state_budget'.or.test_mode=='held_replay'.or. &
     test_mode=='held_natural_replay'.or.test_mode=='held_prepared_natural') then
    ! This failure involves peroxide chemistry absent from the legacy seven-
    ! species fixture. Load the pinned full mechanism for this mode alone.
    call load_h2o2_full_thermo(species,ok)
    call require(ok,'routed replay full thermo')
    call load_h2o2_full_transport(transport,ok)
    call require(ok,'routed replay full transport')
    call load_h2o2_full_mechanism(reactions,ok)
    call require(ok,'routed replay full reactions')
    allocate(base_mole_fractions(size(species)),source=0._dp)
    do i=1,size(species)
      if(trim(species(i)%name)=='N2') base_mole_fractions(i)=1
      if(test_mode=='held_replay') then
        select case(trim(species(i)%name))
        case('H2O');base_mole_fractions(i)=.01_dp
        case('N2');base_mole_fractions(i)=.99_dp
        end select
      endif
    enddo
  endif
  if(test_mode=='seed_aux') species%pinned_eos_extrapolation=.true.
  nv=reactive_nvar(size(species))
  if(test_mode=='held_auxiliary') then
    call require(nranks>=3,'paired source fixture requires owner and idle ranks')
    call check_held_auxiliary_family()
    if(rank==0) print *, 'PASS: actual held auxiliary input/source/ledger family'
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='recursive_srd'.or.test_mode=='common_srd'.or.test_mode=='held_common'.or. &
    test_mode=='nonmol_srd'.or.test_mode=='nonmol_rec'.or. &
    test_mode=='spatial_diff'.or.test_mode=='recursive_diff'.or. &
    test_mode=='nonmol_diff'.or.test_mode=='nonmol_rd') then
    call check_recursive_registers()
    if(rank==0) print *, 'PASS: MPI recursive MOL interface registers'
    call MPI_Finalize(status)
    stop
  endif
  allocate(before(nv),after(nv),boundary(nv),primitive(reactive_nprim(size(species))))
  allocate(wall(nv),chemical(nv))
  do axis=3,3
  do count=5,5
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,2,3,2,boxes(1),ok)
  call require(ok,'child left')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,2,3,2,boxes(2),ok)
  call require(ok,'child right')
  call initialize_amr_patch_3d(4,4,4,3,4,3,4,2,3,2,boxes(3),ok)
  call require(ok,'grandchild left')
  call initialize_amr_patch_3d(4,4,4,1,2,3,4,2,3,2,boxes(4),ok)
  call require(ok,'grandchild right')
  if (count==5) then
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3],boxes,topology,ok)
  else
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],boxes(:3),topology,ok)
  end if
  call require(ok,'topology')
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,axes(axis),.31_dp,finest,ok)
  call require(ok,'finest geometry')
  call coarsen_eb_geometry_3d(finest,2,middle,ok)
  call require(ok,'middle geometry')
  call coarsen_eb_geometry_3d(middle,2,coarse,ok)
  call require(ok,'coarse geometry')
  allocate(u(nv,4,4,4),t(4,4,4))
  config%nx=4;config%ny=4;config%nz=4
  if(test_mode=='routed_replay'.or.test_mode=='routed_predictor'.or. &
     test_mode=='flux_replay'.or.test_mode=='state_replay'.or.test_mode=='state_budget'.or.test_mode=='held_replay'.or. &
     test_mode=='held_natural_replay'.or.test_mode=='held_prepared_natural') then
    call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok,base_mole_fractions)
  else
    call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
  endif
  call require(ok,'root state')
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
  call require(ok,'tree initialization')
  do node=2,count
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
    tree%nodes(node)%state(5,:,:,:)=1.1_dp*tree%nodes(node)%state(5,:,:,:)
    do k=1,4
    do j=1,4
    do i=1,4
      if (tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0._dp) cycle
      call reactive_conserved_to_primitive(species,tree%nodes(node)%state(:,i,j,k),t(i,j,k),primitive, &
        tree%nodes(node)%temperature(i,j,k),speed,ok)
      call require(ok,'temperature perturbation EOS')
    end do
    end do
    end do
  end do
  call average_down_reactive_eb_tree_3d(species,tree,ok)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree
  if(test_mode=='rotating_reaction') then
    call check_rotating_sdc_reaction()
    deallocate(u,t)
    cycle
  end if
  if(test_mode=='stored_reactor') then
    call check_sdc_reaction()
    deallocate(u,t)
    cycle
  endif
  call check_sdc_loop()
  if(test_mode=="held_diffusion".or.test_mode=="mol".or.test_mode=="stored_mol".or.test_mode=="held_policy".or. &
     test_mode=="held_chem_policy".or.test_mode=="held_recursive_policy".or. &
     test_mode=="held_replay".or.test_mode=="held_natural_replay".or. &
     test_mode=="held_prepared".or.test_mode=="held_prepared_recursive".or. &
     test_mode=="held_recursive_chemistry".or.test_mode=="held_prepared_serial".or. &
     test_mode=="held_prepared_parent".or.test_mode=="held_prepared_natural".or. &
     test_mode=="rotating_mol".or.test_mode=="subcycle".or. &
     test_mode=="redist_level".or.test_mode=="redist_recursive".or. &
     test_mode=="state_redist".or.test_mode=="state_recursive".or.test_mode=="aux_source".or.test_mode=="seed_aux".or. &
     test_mode=="aux_recursive".or.test_mode=="aux_route".or.test_mode=="diffusion".or.test_mode=="reactor".or. &
     test_mode=="interval_history".or.test_mode=="recursive_history".or. &
     test_mode=="redist_history".or.test_mode=="route_history".or.test_mode=="route_first".or. &
     test_mode=="routed_replay".or.test_mode=="routed_predictor".or. &
     test_mode=="flux_replay".or.test_mode=="state_replay".or.test_mode=="state_budget") then
    deallocate(u,t)
    cycle
  end if
  call check_sdc_reaction()
  do integrator_index=1,2
  do mode=0,3
    if (integrator_index==2.and.mode/=3) cycle
    chemistry_on=btest(mode,0);transport_on=btest(mode,1)
    tree=saved
    call advance_reactive_eb_tree_full_3d(species,reactions,transport,'hllc',1.e-8_dp,chemistry_on, &
      2.e-7_dp,1.e-12_dp,transport_on,.true.,.true.,.true.,.true.,tree,boundary,wall,chemical,theta,ok, &
      integrators(integrator_index))
    call require(ok,'combined step')
    call integrate_reactive_eb_tree_3d(tree,after,ok)
    call require(ok,'final inventory')
    residual=maxval(abs(after-before-boundary-wall-chemical)/max(1._dp,abs(before)))
    print *, 'integrator, mode, residual: ',integrators(integrator_index),mode,residual
    call require(residual<3.e-14_dp,'combined conservation')
    call require(theta>=0._dp.and.theta<=1._dp,'limiter range')
    if (.not.transport_on) call require(theta==1._dp,'disabled transport limiter')
    if (.not.chemistry_on) call require(all(chemical==0._dp),'disabled chemistry ledger')
    if (chemistry_on) then
      call require(any(chemical(6:)/=0._dp),'nonzero chemical species exchange')
      call require(all(abs(chemical(:5))<=3.e-14_dp*max(1._dp,abs(before(:5)))),'chemical invariants')
    end if
    call check_mpi_full()
  end do
  end do
  deallocate(u,t)
  end do
  end do
  if (rank==0) print *, 'test_mpi_eb_tree_full_3d: PASS'
  call MPI_Finalize(status)
contains
  subroutine check_recursive_registers()
    type(amr_patch_tree_3d) :: nested
    type(amr_patch_3d) :: patches(2)
    type(eb_geometry_3d) :: fg,mg,cg,padded,physical(3)
    type(eb_geometry_3d), allocatable :: wide(:)
    type(eb_reactive_tree_3d) :: candidate,seed,held_baseline_gas
    type(eb_reactive_tree_3d) :: nonmol_reference(2)
    type(amr_tree_state_3d) :: local,lt,kept,kt,held_seed,held_serial,held_local,held_baseline
    type(eb_integrated_flux_3d), allocatable :: held_serial_flux(:),held_local_flux(:),held_baseline_flux(:)
    type(elementary_reaction) :: reactions(0)
    type(state_redist_register_3d), allocatable :: serial_regs(:),mpi_regs(:),kept_regs(:)
    real(dp) :: root_state(nv,4,4,4),root_temp(4,4,4),q(reactive_nprim(size(species))),pos(3),tt,cs,kinetic, &
      b(nv),w(nv),ch(nv),sc(nv),aux(nv),lim,initial(nv),final(nv),error, &
      mb(nv),mw(nv),mch(nv),msc(nv),mlim,mpi_final(nv),bb(nv),bw(nv),bch(nv),bsc(nv),blim
    integer(int64) :: counts(3),mpi_counts(3)
    real(dp) :: reference_b(nv,2),reference_w(nv,2),reference_ch(nv,2),reference_sc(nv,2),reference_lim(2)
    integer(int64) :: reference_counts(3,2)
    integer :: p,x,y,z,c,step,n(3),layout,owners(3),callback_failure
    character(len=32) :: faces(6)
    logical :: good,recursive_mol,spatial_nonmol
    integer :: recursive_iterations
    spatial_nonmol=test_mode=='nonmol_diff'.or.test_mode=='nonmol_rd'
    recursive_mol=test_mode/='nonmol_rec'.and..not.spatial_nonmol
    recursive_iterations=merge(1,2,recursive_mol)
    call initialize_amr_patch_3d(4,4,4,2,3,2,3,2,3,2,patches(1),good)
    call require(good,'recursive ledger middle patch')
    call initialize_amr_patch_3d(4,4,4,2,2,2,2,2,2,2,patches(2),good)
    call require(good,'recursive ledger properly nested leaf patch')
    call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,nested,good)
    call require(good,'recursive ledger three-level topology')
    call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.43_dp,fg,good)
    call require(good,'recursive ledger fine geometry')
    call coarsen_eb_geometry_3d(fg,2,mg,good)
    call require(good,'recursive ledger middle geometry')
    call coarsen_eb_geometry_3d(mg,2,cg,good)
    call require(good,'recursive ledger root geometry')
    config%nx=4;config%ny=4;config%nz=4
    call initialize_reactive_eb_density_sheet_3d(species,config,cg,root_state,root_temp,good)
    call require(good,'recursive ledger root gas')
    call initialize_reactive_eb_tree_3d(species,nested,fg,root_state,root_temp,candidate,good)
    call require(good,'recursive ledger initialized tree')
    do p=1,3
      physical(p)=candidate%nodes(p)%geometry;n=nested%nodes(p)%shape
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        if(physical(p)%volume_fraction(x,y,z)<=0) cycle
        pos=[physical(p)%x_lower+(x-.5_dp)*physical(p)%dx,physical(p)%y_lower+(y-.5_dp)*physical(p)%dy, &
          physical(p)%z_lower+(z-.5_dp)*physical(p)%dz]
        q=0;q(1)=1+.01_dp*pos(1);q(2:4)=.1_dp*pos;q(5)=250000+1000*pos(2)
        do c=1,size(species)
          if(species(c)%name=='N2') q(c+5)=.95_dp
          if(species(c)%name=='H2') q(c+5)=.05_dp
        enddo
        call reactive_primitive_to_conserved(species,q,candidate%nodes(p)%state(:,x,y,z),tt,cs,good)
        call require(good,'recursive ledger smooth EOS')
        candidate%nodes(p)%temperature(x,y,z)=tt
      enddo
      enddo
      enddo
    enddo
    call average_down_reactive_eb_tree_3d(species,candidate,good)
    call require(good,'recursive ledger initial restriction')
    if(test_mode=='nonmol_srd'.or.test_mode=='nonmol_rec'.or.spatial_nonmol) then
      call build_axis_plane_eb_geometry_3d(88,88,88,-2.25_dp,3.25_dp,-2.25_dp,3.25_dp,-2.25_dp,3.25_dp, &
        'z',.43_dp,padded,good)
    else
      call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
        'z',.43_dp,padded,good)
    endif
    call require(good,'recursive ledger padded geometry')
    call build_eb_tree_extended_geometry_3d(nested,physical,padded,wide,good, &
      halo_width=merge(9,6,test_mode=='nonmol_srd'.or.test_mode=='nonmol_rec'.or.spatial_nonmol))
    call require(good,'recursive ledger extended metrics')
    faces='outflow';seed=candidate
    if(test_mode=='held_common') then
      allocate(held_seed%fields(3))
      do p=1,3
        n=nested%nodes(p)%shape
        allocate(held_seed%fields(p)%values(1,n(1),n(2),n(3)),source=0._dp)
        do z=1,n(3)
        do y=1,n(2)
        do x=1,n(1)
          if(physical(p)%volume_fraction(x,y,z)<=0) cycle
          kinetic=.5_dp*sum(seed%nodes(p)%state(2:4,x,y,z)**2)/seed%nodes(p)%state(1,x,y,z)
          held_seed%fields(p)%values(1,x,y,z)=(seed%nodes(p)%state(5,x,y,z)-kinetic)* &
            (.9_dp+.01_dp*real(x-1,dp)+.01_dp*real(z-1,dp))
        enddo
        enddo
        enddo
      enddo
      call average_down_reactive_eb_tree_field_3d(seed,held_seed,good)
      call require(good,'properly nested held-state initial restriction')
    endif
    spatial_topology=nested;diffusion_geometry=physical;spatial_geometry=wide
    do layout=1,3
      owners=[0,1,mod(2,nranks)]
      if(layout==2) owners=nranks-1-owners
      if(layout==3) owners=0
      candidate=seed;diffusion_owners=owners;spatial_failure=0;spatial_calls=0
      if(test_mode=='held_common') then
        held_serial=held_seed
        allocate(held_local%fields(3))
        if(allocated(held_serial_flux)) deallocate(held_serial_flux)
        if(allocated(held_local_flux)) deallocate(held_local_flux)
      endif
      allocate(local%fields(3),lt%fields(3))
      do p=1,3
        if(owners(p)/=rank) cycle
        local%fields(p)%values=seed%nodes(p)%state;n=nested%nodes(p)%shape
        allocate(lt%fields(p)%values(1,n(1),n(2),n(3)))
        lt%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
        if(test_mode=='held_common') held_local%fields(p)=held_seed%fields(p)
      enddo
    do step=1,2
      call integrate_reactive_eb_tree_3d(candidate,initial,good)
      call require(good,'recursive ledger initial inventory')
      if(test_mode=='recursive_diff'.or.test_mode=='nonmol_rd') then
        call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,recursive_iterations,.false., &
          [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
          b,w,ch,sc,lim,good,level_steps=counts,redistribution='state_redist',mol_mode=recursive_mol,state_redist_max_order=2)
        call require(good,'recursive callback serial built-in transport control')
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,nested,owners,physical,local,lt, &
          1.e-10_dp,recursive_iterations,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
          1.e-6_dp,1.e-12_dp,wide,faces, &
          mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,level_steps=mpi_counts,redistribution='state_redist',mol_mode=recursive_mol, &
          state_redist_max_order=2,evaluate_spatial_diffusion=spatial_diffusion_source)
        call require(good,'recursive callback MPI transport accepted')
        call require(all(mpi_counts==counts),'recursive callback counts agree')
      elseif(test_mode=='spatial_diff'.or.test_mode=='nonmol_diff') then
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,recursive_iterations,.false., &
          [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
          b,w,ch,sc,lim,good,redistribution='state_redist',mol_mode=recursive_mol,state_redist_max_order=2)
        call require(good,'spatial callback serial built-in transport control')
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,nested,owners,physical,local,lt,wide, &
          1.e-10_dp,recursive_iterations,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,faces, &
          mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,redistribution='state_redist',mol_mode=recursive_mol, &
          state_redist_max_order=2,state_redist_registers=mpi_regs,evaluate_spatial_diffusion=spatial_diffusion_source)
        call require(good,'MPI spatial callback reproduces built-in transport')
        do p=1,3
          if(owners(p)/=rank) cycle
          if(allocated(mpi_regs(p)%coarse_delta)) call require(all(mpi_regs(p)%coarse_delta==0),'callback coarse consumed')
          if(allocated(mpi_regs(p)%fine_delta)) call require(all(mpi_regs(p)%fine_delta==0),'callback fine consumed')
        enddo
      elseif(test_mode=='nonmol_srd') then
        if(step==1) then
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.false., &
          [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
          b,w,ch,sc,lim,good,redistribution='state_redist',mol_mode=.false.,state_redist_max_order=2, &
          embedded_wall_temperature=700._dp,embedded_noslip=.true.)
        call require(good,'serial common non-MOL two-iteration register close')
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,nested,owners,physical,local,lt,wide, &
          1.e-10_dp,2,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,faces, &
          mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,redistribution='state_redist',mol_mode=.false., &
          state_redist_max_order=2,embedded_wall_temperature=700._dp,embedded_noslip=.true.)
        call require(good,'MPI common non-MOL two-iteration register close')
        else
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.false., &
          [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
          b,w,ch,sc,lim,good,redistribution='state_redist',mol_mode=.false.,state_redist_max_order=2, &
          embedded_wall_temperature=700._dp,embedded_noslip=.true.,state_redist_registers=serial_regs)
        call require(good,'serial common non-MOL two-iteration register close')
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,nested,owners,physical,local,lt,wide, &
          1.e-10_dp,2,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,faces, &
          mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,redistribution='state_redist',mol_mode=.false., &
          state_redist_max_order=2,embedded_wall_temperature=700._dp,embedded_noslip=.true., &
          state_redist_registers=mpi_regs)
        call require(good,'MPI common non-MOL two-iteration register close')
        endif
        if(step==2) then
        do p=1,3
          if(allocated(serial_regs(p)%coarse_delta)) &
            call require(all(serial_regs(p)%coarse_delta==0),'serial common non-MOL coarse consumed')
          if(allocated(serial_regs(p)%fine_delta)) &
            call require(all(serial_regs(p)%fine_delta==0),'serial common non-MOL fine consumed')
          if(owners(p)/=rank) then
            call require(.not.allocated(mpi_regs(p)%coarse_delta).and..not.allocated(mpi_regs(p)%fine_delta), &
              'common non-MOL remote registers absent')
            cycle
          endif
          if(allocated(mpi_regs(p)%coarse_delta)) &
            call require(all(mpi_regs(p)%coarse_delta==0),'MPI common non-MOL coarse consumed')
          if(allocated(mpi_regs(p)%fine_delta)) &
            call require(all(mpi_regs(p)%fine_delta==0),'MPI common non-MOL fine consumed')
        enddo
        endif
      elseif(test_mode=='common_srd') then
        if(step==1) then
          call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,1,.false., &
            [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
            b,w,ch,sc,lim,good,redistribution='state_redist',mol_mode=.true.,state_redist_max_order=2, &
            embedded_wall_temperature=700._dp,embedded_noslip=.true.)
          call require(good,'serial common implicit register route')
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,nested,owners,physical,local,lt,wide, &
            1.e-10_dp,1,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,faces, &
            mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,redistribution='state_redist',mol_mode=.true., &
            state_redist_max_order=2,embedded_wall_temperature=700._dp,embedded_noslip=.true.)
        else
          call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,1,.false., &
            [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
            b,w,ch,sc,lim,good,redistribution='state_redist',mol_mode=.true.,state_redist_max_order=2, &
            embedded_wall_temperature=700._dp,embedded_noslip=.true.,state_redist_registers=serial_regs)
          call require(good,'serial common explicit register output')
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,nested,owners,physical,local,lt,wide, &
            1.e-10_dp,1,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,faces, &
            mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,redistribution='state_redist',mol_mode=.true., &
            state_redist_max_order=2,embedded_wall_temperature=700._dp,embedded_noslip=.true., &
            state_redist_registers=mpi_regs)
        endif
        call require(good,'MPI common MOL accepted')
        if(step==2) then
          do p=1,3
            if(allocated(serial_regs(p)%coarse_delta)) &
              call require(all(serial_regs(p)%coarse_delta==0),'serial common coarse consumed')
            if(allocated(serial_regs(p)%fine_delta)) &
              call require(all(serial_regs(p)%fine_delta==0),'serial common fine consumed')
            if(owners(p)/=rank) cycle
            if(allocated(mpi_regs(p)%coarse_delta)) &
              call require(all(mpi_regs(p)%coarse_delta==0),'MPI common coarse consumed')
            if(allocated(mpi_regs(p)%fine_delta)) &
              call require(all(mpi_regs(p)%fine_delta==0),'MPI common fine consumed')
          enddo
        endif
      elseif(test_mode=='held_common') then
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,1,.false., &
          [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
          b,w,ch,sc,lim,good,redistribution='state_redist',mol_mode=.true.,target_volume_fraction=.8_dp, &
          internal_energy=held_serial,internal_energy_flux=held_serial_flux)
        call require(good,'serial properly nested held-state StateRedist')
        if(layout==1.and.step==1) then
          held_baseline_gas=seed;held_baseline=held_seed
          call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,1,.false., &
            [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp, &
            held_baseline_gas,wide,faces,bb,bw,bch,bsc,blim,good,mol_mode=.true., &
            internal_energy=held_baseline,internal_energy_flux=held_baseline_flux)
          call require(good,'properly nested held-state baseline MOL')
          call require(any(held_serial%fields(3)%values/=held_baseline%fields(3)%values), &
            'properly nested StateRedist changes held field')
        endif
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,nested,owners,physical,local,lt,wide, &
          1.e-10_dp,1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,faces, &
          mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,redistribution='state_redist',mol_mode=.true., &
          target_volume_fraction=.8_dp,internal_energy=held_local,internal_energy_flux=held_local_flux)
        call require(good,'MPI properly nested held-state StateRedist')
      else
      ! This serial trajectory has no owner map or callback state. Reuse the
      ! same two independently computed endpoints for every MPI ownership
      ! layout; all six MPI updates and each layout's rollback still execute.
      if(test_mode=='nonmol_rec'.and.layout>1) then
        candidate=nonmol_reference(step)
        b=reference_b(:,step);w=reference_w(:,step)
        ch=reference_ch(:,step);sc=reference_sc(:,step)
        lim=reference_lim(step);counts=reference_counts(:,step)
      else
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,recursive_iterations,.false., &
        [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
        b,w,ch,sc,lim,good,level_steps=counts,redistribution='state_redist',mol_mode=recursive_mol, &
        state_redist_max_order=2,embedded_wall_temperature=700._dp,embedded_noslip=.true.)
      call require(good,'recursive ledger accepted three-level interval')
      if(test_mode=='nonmol_rec') then
        nonmol_reference(step)=candidate
        reference_b(:,step)=b;reference_w(:,step)=w
        reference_ch(:,step)=ch;reference_sc(:,step)=sc
        reference_lim(step)=lim;reference_counts(:,step)=counts
      endif
      endif
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,nested,owners,physical,local,lt, &
        1.e-10_dp,recursive_iterations,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,wide,faces, &
        mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,level_steps=mpi_counts,redistribution='state_redist',mol_mode=recursive_mol, &
        state_redist_max_order=2,embedded_wall_temperature=700._dp,embedded_noslip=.true.)
      call require(good,'MPI recursive accepted three-level interval')
      call require(all(mpi_counts==counts),'MPI recursive count comparison')
      endif
      call require(maxval(abs(mb-b)/max(1._dp,abs(b)))<2.e-14_dp.and. &
        maxval(abs(mw-w)/max(1._dp,abs(w)))<2.e-14_dp.and.all(mch==ch).and.all(msc==sc).and.mlim==lim, &
        'MPI recursive diagnostics agree with serial')
      do p=1,3
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'recursive owner-local allocation')
        if(owners(p)/=rank) cycle
        call require(maxval(abs(local%fields(p)%values-candidate%nodes(p)%state)/ &
          max(1._dp,abs(candidate%nodes(p)%state)))<2.e-14_dp,'recursive serial/MPI gas comparison')
        call require(maxval(abs(lt%fields(p)%values(1,:,:,:)-candidate%nodes(p)%temperature)/ &
          max(1._dp,abs(candidate%nodes(p)%temperature)))<2.e-14_dp,'recursive serial/MPI temperature comparison')
        if(test_mode=='held_common') call require(all(held_local%fields(p)%values==held_serial%fields(p)%values), &
          'properly nested held-state exact serial/MPI comparison')
        if(test_mode=='held_common') call require(all(held_local_flux(p)%x==held_serial_flux(p)%x).and. &
          all(held_local_flux(p)%y==held_serial_flux(p)%y).and. &
          all(held_local_flux(p)%z==held_serial_flux(p)%z), &
          'properly nested held-state exact serial/MPI flux')
      enddo
      if(test_mode/='common_srd'.and.test_mode/='held_common'.and. &
        test_mode/='spatial_diff'.and.test_mode/='nonmol_diff'.and.test_mode/='nonmol_srd') &
        call require(all(counts==[1_int64,2_int64,4_int64]),'recursive ledger substep counts')
      call integrate_reactive_eb_tree_3d(candidate,final,good)
      call require(good,'recursive ledger final inventory')
      error=maxval(abs(final-initial-b-w-ch-sc)/max(1._dp,abs(initial)))
      call integrate_mpi_eb_tree_3d(nested,owners,physical,local,mpi_final,MPI_COMM_WORLD,good)
      call require(good,'MPI recursive inventory')
      error=max(error,maxval(abs(mpi_final-initial-mb-mw-mch-msc)/max(1._dp,abs(initial))))
      if(rank==0) print *, 'MPI recursive layout, step, composite budget error: ',layout,step,error
      call require(error<2.e-12_dp,'recursive ledger composite gas and body wall budget')
    enddo
    if(test_mode=='held_common') then
      deallocate(local%fields,lt%fields,held_local%fields)
      cycle
    endif
    kept=local;kt=lt;mb=17;mw=19;mch=23;msc=29;aux=31;mlim=37;mpi_counts=41
    source_failure=owners(3)
    if(test_mode=='spatial_diff'.or.test_mode=='recursive_diff'.or.spatial_nonmol) then
      do callback_failure=1,5
      spatial_failure=callback_failure;spatial_calls=0
      if(rank==0.and.spatial_nonmol) print *, 'non-MOL callback failure begin: ',layout,callback_failure
      if(test_mode=='recursive_diff'.or.test_mode=='nonmol_rd') then
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,nested,owners,physical,local,lt, &
          1.e-10_dp,recursive_iterations,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
          1.e-6_dp,1.e-12_dp,wide,faces, &
          mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,level_steps=mpi_counts,redistribution='state_redist',mol_mode=recursive_mol, &
          state_redist_max_order=2,evaluate_spatial_diffusion=spatial_diffusion_source)
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,nested,owners,physical,local,lt,wide, &
        1.e-10_dp,recursive_iterations,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,faces, &
        mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,redistribution='state_redist',mol_mode=recursive_mol, &
        state_redist_max_order=2,evaluate_spatial_diffusion=spatial_diffusion_source)
      endif
      call require(.not.good,'spatial callback failure collectively rejected')
      call require(all(mb==17).and.all(mw==19).and.all(mch==23).and.all(msc==29).and.mlim==37, &
        'spatial callback boundary source and limiter rollback')
      call require(all(mpi_counts==41),'spatial callback substep-count rollback')
      if(rank==0.and.spatial_nonmol) print *, 'non-MOL callback failure rejected: ',layout,callback_failure,spatial_calls
      do p=1,3
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==kept%fields(p)%values).and. &
          all(lt%fields(p)%values==kt%fields(p)%values),'spatial callback outer gas rollback')
      enddo
      enddo
    elseif(test_mode=='nonmol_srd') then
      kept_regs=mpi_regs
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,nested,owners,physical,local,lt,wide, &
        1.e-10_dp,2,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,faces, &
        mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,redistribution='state_redist',mol_mode=.false., &
        state_redist_max_order=2,auxiliary_source=reject_common_register_endpoint,auxiliary_change=aux, &
        state_redist_registers=mpi_regs)
      call require(.not.good,'MPI common non-MOL endpoint rejection')
      do p=1,3
        if(owners(p)/=rank) cycle
        if(allocated(mpi_regs(p)%coarse_delta)) &
          call require(all(mpi_regs(p)%coarse_delta==kept_regs(p)%coarse_delta),'common non-MOL coarse rollback')
        if(allocated(mpi_regs(p)%fine_delta)) &
          call require(all(mpi_regs(p)%fine_delta==kept_regs(p)%fine_delta),'common non-MOL fine rollback')
      enddo
    elseif(test_mode=='common_srd') then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,nested,owners,physical,local,lt,wide, &
        1.e-10_dp,1,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,faces, &
        mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,redistribution='state_redist',mol_mode=.true., &
        state_redist_max_order=2,auxiliary_source=reject_common_register_endpoint,auxiliary_change=aux)
    else
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,nested,owners,physical,local,lt, &
      1.e-10_dp,recursive_iterations,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,wide,faces, &
      mb,mw,mch,msc,mlim,MPI_COMM_WORLD,good,level_steps=mpi_counts,redistribution='state_redist',mol_mode=recursive_mol, &
      state_redist_max_order=2,auxiliary_source=reject_recursive_register_endpoint,auxiliary_change=aux)
    endif
    call require(.not.good,'MPI recursive late descendant rejection')
    call require(all(mb==17).and.all(mw==19).and.all(mch==23).and.all(msc==29).and.all(aux==31).and.mlim==37.and. &
      all(mpi_counts==41),'MPI recursive diagnostic and count rollback')
    do p=1,3
      if(owners(p)/=rank) cycle
      call require(all(local%fields(p)%values==kept%fields(p)%values).and. &
        all(lt%fields(p)%values==kt%fields(p)%values),'MPI recursive outer gas rollback')
    enddo
    deallocate(local%fields,lt%fields)
    enddo
  end subroutine check_recursive_registers

  subroutine reject_common_register_endpoint(snapshot,temperature,interval,stage,depth,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: snapshot,temperature
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p
    good=rank/=source_failure.or.stage==0
    if(.not.good) return
    do p=1,size(rate%fields)
      if(allocated(rate%fields(p)%values)) rate%fields(p)%values=0
    enddo
  end subroutine reject_common_register_endpoint
  subroutine reject_recursive_register_endpoint(snapshot,temperature,interval,stage,depth,time,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: snapshot,temperature
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p
    good=rank/=source_failure.or.stage==0.or.depth/=2
    if(.not.good) return
    do p=1,size(rate%fields)
      if(allocated(rate%fields(p)%values)) rate%fields(p)%values=0
    enddo
  end subroutine reject_recursive_register_endpoint
  subroutine check_rotating_sdc_reaction()
    type(rotating_frame_3d) :: frame
    type(eb_reactive_tree_3d) :: rotated
    type(amr_tree_state_3d) :: rates,reference,reference_t,reference_ir,old,ot,lr,ir,actual,actual_t,kept,kept_t,kept_ir
    type(eb_geometry_3d) :: geometry(count)
    real(dp) :: position(3),potential
    integer :: owners(count),p,x,y,z,layout
    logical :: good
    frame%axis=2;frame%angular_speed=120._dp;frame%axis_point=[.25_dp,.4_dp,.5_dp]
    rotated=saved
    allocate(rates%fields(count),reference_ir%fields(count))
    do p=1,count
      geometry(p)=rotated%nodes(p)%geometry
      allocate(rates%fields(p)%values(nv,4,4,4),source=0._dp)
      allocate(reference_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        position=[geometry(p)%x_lower,geometry(p)%y_lower,geometry(p)%z_lower]+ &
          (real([x,y,z],dp)-.5_dp)*[geometry(p)%dx,geometry(p)%dy,geometry(p)%dz]
        call rotating_potential_3d(frame,position,potential,good)
        call require(good,'rotating reaction potential')
        rotated%nodes(p)%state(5,x,y,z)=rotated%nodes(p)%state(5,x,y,z)- &
          rotated%nodes(p)%state(1,x,y,z)*potential
        rates%fields(p)%values(2,x,y,z)=2.e4_dp
        rates%fields(p)%values(5,x,y,z)=2.e4_dp*rotated%nodes(p)%state(2,x,y,z)/ &
          rotated%nodes(p)%state(1,x,y,z)
      end do
      end do
      end do
    end do
    call react_sdc_eb_tree_3d(species,reactions,rotated,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp, &
      reference,reference_t,reference_ir,good,frame=frame)
    call require(good,'serial rotating SDC reaction')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(old%fields(count),ot%fields(count),lr%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        old%fields(p)%values=rotated%nodes(p)%state
        allocate(ot%fields(p)%values(1,4,4,4))
        ot%fields(p)%values(1,:,:,:)=rotated%nodes(p)%temperature
        lr%fields(p)%values=rates%fields(p)%values
        allocate(ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
      end do
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,actual,actual_t,ir,MPI_COMM_WORLD,good,frame=frame)
      call require(good,'MPI rotating SDC reaction')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(actual%fields(p)%values==reference%fields(p)%values).and. &
          all(actual_t%fields(p)%values==reference_t%fields(p)%values).and. &
          all(ir%fields(p)%values==reference_ir%fields(p)%values),'MPI exact serial rotating reaction')
      end do
      kept=actual;kept_t=actual_t;kept_ir=ir
      if(rank==0) frame%angular_speed=121._dp
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,actual,actual_t,ir,MPI_COMM_WORLD,good,frame=frame)
      frame%angular_speed=120._dp
      call require(.not.good,'MPI rotating frame mismatch rejected')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(actual%fields(p)%values==kept%fields(p)%values).and. &
          all(actual_t%fields(p)%values==kept_t%fields(p)%values).and. &
          all(ir%fields(p)%values==kept_ir%fields(p)%values),'MPI rotating reaction rollback')
      end do
      deallocate(old%fields,ot%fields,lr%fields,ir%fields)
    end do
  end subroutine check_rotating_sdc_reaction
  subroutine reactor_adapter(species,reactions,initial_energy_density,external_species,external_energy, &
      duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
    real(dp), intent(inout) :: species_density(:),temperature
    integer, intent(out) :: completed_steps
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    real(dp) :: timer_start
    reactor_calls=reactor_calls+1
    timer_start=0
    if((test_mode=='held_chem_policy'.or.test_mode=='held_replay'.or.test_mode=='held_prepared').and.rank==0) then
      timer_start=MPI_Wtime()
      if(reactor_calls<=2.or.mod(reactor_calls,32)==0) &
        print *, 'Held coupled native entry: ',reactor_calls,initial_energy_density,external_energy,temperature
    endif
    if(reactor_replay_failure.and.rank==reactor_failure_rank.and. &
        history_recursive_calls==14.and.duration==reactor_replay_duration) then
      reactor_replay_calls=reactor_replay_calls+1
    endif
    if(reactor_failure.and.rank==reactor_failure_rank.and.reactor_calls==2) then
      species_density=-1;temperature=-1;completed_steps=0;ok=.false.
      return
    end if
    call advance_forced_reactor_interval(species,reactions,initial_energy_density,external_species,external_energy, &
      duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps=maximum_steps)
    if(material_refusal_probe.and..not.ok.and.rank==0) &
      print *, 'MATERIAL_REACTOR_STOP',history_recursive_calls,reactor_calls,duration, &
        reactor_replay_calls,reactor_replay_refused
    if((test_mode=='held_chem_policy'.or.test_mode=='held_replay'.or.test_mode=='held_prepared').and.rank==0) then
      if(reactor_calls<=2.or.mod(reactor_calls,32)==0) &
        print *, 'Held coupled native exit: ',reactor_calls,completed_steps,ok,MPI_Wtime()-timer_start
    endif
    if(ok.and.reactor_replay_failure.and.rank==reactor_failure_rank.and. &
        history_recursive_calls==14.and.duration==reactor_replay_duration.and. &
        reactor_replay_calls>=reactor_replay_refusal_call.and..not.reactor_replay_refused) then
      ! Refuse only an input just proven admissible by the native integrator.
      ! The caller must preserve the provider's refusal, not adopt that solve.
      reactor_replay_refused=.true.;reactor_replay_refusal_actual=reactor_replay_calls
      species_density=-1;temperature=-1;completed_steps=0;ok=.false.
    endif
    if(ok.and.prepared_parent_refusal_call>0.and.rank==reactor_failure_rank.and. &
        reactor_calls==prepared_parent_refusal_call) then
      reactor_replay_refused=.true.;reactor_replay_refusal_actual=reactor_calls
      species_density=-1;temperature=-1;completed_steps=0;ok=.false.
    endif
  end subroutine reactor_adapter
  subroutine spatial_diffusion_source(state,temp,geometry,old_halos,interval,depth,rate,spatial_rate,flux, &
      limiter,comm,good,frame)
    type(amr_tree_state_3d), intent(in) :: state,temp,old_halos
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(in) :: interval
    integer, intent(in) :: depth
    type(amr_tree_state_3d), intent(inout) :: rate,spatial_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: limiter
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), allocatable :: temperatures(:,:,:),misplaced(:,:,:,:)
    real(dp) :: q(reactive_nprim(size(species))),recovered,sound,patch_theta
    integer :: p,x,y,z,n(3)
    logical :: faces(6)
    character(len=16) :: labels(6)
    good=.false.;limiter=1;faces=.true.;labels='outflow'
    allocate(spatial_rate%fields(size(state%fields)),flux(size(state%fields)))
    do p=1,size(state%fields)
      if(.not.allocated(old_halos%fields(p)%values)) cycle
      n=[geometry(p)%nx,geometry(p)%ny,geometry(p)%nz]
      allocate(temperatures(n(1),n(2),n(3)),source=300._dp)
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        call reactive_conserved_to_primitive(species,old_halos%fields(p)%values(:,x,y,z),300._dp,q, &
          recovered,sound,good)
        if(.not.good) return
        temperatures(x,y,z)=recovered
      enddo
      enddo
      enddo
      call compute_sdc_eb_transport_halo_3d(species,transport,old_halos%fields(p)%values,temperatures, &
        geometry(p),interval,[.true.,.true.,.true.,.true.],rate%fields(p)%values,flux(p), &
        spatial_rate%fields(p)%values,patch_theta,good,faces,labels)
      if(.not.good) return
      limiter=min(limiter,patch_theta)
      deallocate(temperatures)
    enddo
    good=.true.
    spatial_calls=spatial_calls+1
    if(rank==diffusion_owners(3).and.mod(spatial_calls,2)==0.and.(depth<0.or.depth==2)) then
      if(spatial_failure==1) good=.false.
      if(spatial_failure==2) spatial_rate%fields(3)%values(1,7,7,7)=spatial_rate%fields(3)%values(1,7,7,7)+1
      if(spatial_failure==3) spatial_rate%fields(3)%values(1,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
      if(spatial_failure==4.or.spatial_failure==5) then
        ! Zero the callback contribution so its valid-crop equality cannot
        ! mask the malformed bounds. Both container shapes remain correct.
        rate%fields(3)%values=0;spatial_rate%fields(3)%values=0
        flux(3)%x=0;flux(3)%y=0;flux(3)%z=0
        if(spatial_failure==4) then
          n=shape(spatial_rate%fields(3)%values(1,:,:,:))
          allocate(misplaced(nv,-5:n(1)-6,-5:n(2)-6,-5:n(3)-6),source=0._dp)
          call move_alloc(misplaced,spatial_rate%fields(3)%values)
        else
          n=shape(rate%fields(3)%values(1,:,:,:))
          allocate(misplaced(nv,0:n(1)-1,0:n(2)-1,0:n(3)-1),source=0._dp)
          call move_alloc(misplaced,rate%fields(3)%values)
        endif
      endif
    endif
  end subroutine spatial_diffusion_source

  subroutine diffusion_source(state,temp,interval,depth,rate,flux,limiter,comm,good, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,fractions)
    type(amr_tree_state_3d), intent(in) :: state,temp
    real(dp), intent(in) :: interval
    integer, intent(in) :: depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: limiter
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: fractions(:)
    type(amr_tree_state_3d) :: source_halo,old_halo
    character(len=16) :: faces(6)
    integer :: p
    integer, allocatable :: selected
    faces='outflow'
    if(depth>=0) selected=depth
    ! Supply the same complete transport stencil as the built-in CTU path.
    call compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,topology,diffusion_owners,diffusion_geometry, &
      spatial_geometry,state,temp,interval,[.true.,.true.,.true.,.true.],rate,flux,source_halo,old_halo,limiter, &
      comm,good,boundary_faces=faces,target_depth=selected, &
      ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
      temperature_end=temperature_end,ancestor_fraction=fractions)
    if(good) then
      ! The spatial callback's density entry is deliberately inconsistent.
      ! Native EB postprocessing must replace it with the unchanged species
      ! source sum before stepping, retaining exact serial/MPI gas and ledgers.
      do p=1,size(diffusion_owners)
        if(diffusion_owners(p)/=rank) cycle
        if(depth>=0.and.topology%nodes(p)%level/=depth) cycle
        rate%fields(p)%values(1,:,:,:)=99._dp
      enddo
    endif
    diffusion_calls=diffusion_calls+1
    if(diffusion_failure.and.rank==nranks-1.and.diffusion_calls==2) good=.false.
  end subroutine diffusion_source
  subroutine check_sdc_loop()
    type(rotating_frame_3d) :: frame
    type(eb_reactive_tree_3d) :: initial,reference
    type(eb_geometry_3d) :: geometry(count),extended_finest
    type(eb_geometry_3d), allocatable :: extended(:)
    type(amr_tree_state_3d) :: local,lt,ir,serial_ir,old_ir,old,ot
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),sb(nv),sw(nv),sch(nv),sext(nv),limiter,slimiter
    real(dp) :: initial_inventory(nv),final_inventory(nv),error,force(3),viscosity_coefficient,recovered,sound
    real(dp) :: position(3),potential
    integer :: owners(count),p,x,y,z,layout,failure,expected_reactor_calls,failed_order
    character(len=16) :: faces(6),method
    logical :: good
    method='none'
    if(test_mode=='state_redist') method='state_redist'
    initial=saved;faces='outflow';force=[2.e4_dp,-3.e4_dp,4.e4_dp]
    do p=1,count
      geometry(p)=initial%nodes(p)%geometry
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        call reactive_conserved_to_primitive(species,initial%nodes(p)%state(:,x,y,z), &
          initial%nodes(p)%temperature(x,y,z),primitive,recovered,sound,good)
        call require(good,'SDC fixture primitive')
        primitive(5)=primitive(5)*1000._dp/recovered
        call reactive_primitive_to_conserved(species,primitive,initial%nodes(p)%state(:,x,y,z), &
          initial%nodes(p)%temperature(x,y,z),sound,good)
        call require(good,'SDC fixture EOS')
      end do
      end do
      end do
    end do
    call average_down_reactive_eb_tree_3d(species,initial,good)
    call require(good,'SDC fixture restriction')
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'z',.31_dp,extended_finest,good)
    call require(good,'SDC extended finest')
    call build_eb_tree_extended_geometry_3d(topology,geometry,extended_finest,extended,good)
    call require(good,'SDC extended geometry')
    if(test_mode=='held_diffusion') then
      call check_held_valid_diffusion(initial,geometry,extended,faces)
      return
    endif
    if(test_mode=='redist_history') then
      call check_redist_interval_history(initial,geometry,extended,faces,force)
      return
    endif
    if(test_mode=='routed_replay'.or.test_mode=='routed_predictor'.or. &
     test_mode=='flux_replay'.or.test_mode=='state_replay'.or.test_mode=='state_budget') then
      call check_three_level_routed_replay(initial,geometry,extended,faces,force)
      return
    endif
    if(test_mode=='route_history'.or.test_mode=='route_first') then
      call check_routed_mol_history(initial,geometry,extended,faces,force)
      return
    endif
    if(test_mode=='recursive_history') then
      call check_recursive_interval_history(initial,geometry,extended,faces,force)
      return
    endif
    if(test_mode=='interval_history') then
      call check_interval_history(initial,geometry,extended,faces,force)
      return
    endif
    if(test_mode=='stored_mol') then
      call check_stored_mol(initial,geometry,extended,faces,force)
      return
    endif
    if(test_mode=='held_prepared'.or.test_mode=='held_prepared_recursive'.or. &
    test_mode=='held_recursive_chemistry'.or.test_mode=='held_prepared_serial'.or.test_mode=='held_prepared_parent') then
      call check_held_prepared(initial,geometry,extended,faces)
      return
    endif
    if(test_mode=='held_natural_replay'.or.test_mode=='held_prepared_natural') then
      call check_held_natural_replay(initial,geometry,extended,faces,force)
      return
    endif
    if(test_mode=='held_replay') then
      call check_held_replay(initial,geometry,extended,faces,force)
      return
    endif
    if(test_mode=='held_recursive_policy') then
      call check_held_recursive_policy(initial,geometry,extended,faces)
      return
    endif
    if(test_mode=='held_policy'.or.test_mode=='held_chem_policy') then
      call check_held_common_policy(initial,geometry,extended,faces)
      return
    endif
    if(test_mode=="mol".or.test_mode=="rotating_mol") then
      if(test_mode=="rotating_mol") then
        frame%axis=2;frame%angular_speed=120._dp;frame%axis_point=[.25_dp,.4_dp,.5_dp]
        do p=1,count
        do z=1,4
        do y=1,4
        do x=1,4
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          position=[geometry(p)%x_lower,geometry(p)%y_lower,geometry(p)%z_lower]+ &
            (real([x,y,z],dp)-.5_dp)*[geometry(p)%dx,geometry(p)%dy,geometry(p)%dz]
          call rotating_potential_3d(frame,position,potential,good)
          call require(good,'rotating MOL potential')
          initial%nodes(p)%state(5,x,y,z)=initial%nodes(p)%state(5,x,y,z)- &
            initial%nodes(p)%state(1,x,y,z)*potential
        end do
        end do
        end do
        end do
        call check_mol_parity(initial,geometry,extended,faces,force,frame)
      else
        call check_mol_parity(initial,geometry,extended,faces,force)
      end if
      return
    end if
    if(test_mode=="aux_source".or.test_mode=="seed_aux") then
      call check_auxiliary_sources(initial,geometry,extended,faces,force)
      return
    end if
    if(test_mode=="redist_level") then
      call check_redist_level(initial,geometry,extended,faces,force)
      return
    end if
    if(test_mode=="subcycle".or.test_mode=="redist_recursive".or.test_mode=="state_recursive".or. &
     test_mode=="aux_recursive".or.test_mode=="aux_route") then
      call check_sdc_recursive(initial,geometry,extended,faces,force)
      return
    end if
    allocate(serial_ir%fields(count))
    do p=1,count
      allocate(serial_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
    end do
    reference=initial
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,slimiter,good,serial_ir,redistribution=method, &
      state_redist_max_order=requested_order)
    call require(good,'serial coupled SDC')
    call require_restricted_reactions(reference,serial_ir)
    call integrate_reactive_eb_tree_3d(initial,initial_inventory,good)
    call require(good,'SDC initial inventory')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(local%fields(count),lt%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=initial%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4))
        lt%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
        ir%fields(p)%values=0
      end do
      old=local;ot=lt;old_ir=ir
      reactor_calls=0;expected_reactor_calls=0
      do p=1,count
        if(owners(p)/=rank) cycle
        expected_reactor_calls=expected_reactor_calls+2*sum(merge(1,0,geometry(p)%volume_fraction>0))
      end do
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
        1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,redistribution=method,evaluate_reaction=reactor_adapter, &
      state_redist_max_order=requested_order)
      call require(good,'MPI coupled SDC')
      call require(reactor_calls==expected_reactor_calls,'callback only on owned noncovered cells')
      reactor_calls=0
      if(test_mode=='diffusion') then
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==reference%nodes(p)%state),'built-in serial MPI CTU transport state')
          call require(all(lt%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'built-in CTU temperature')
          call require(all(ir%fields(p)%values==serial_ir%fields(p)%values),'built-in CTU reaction history')
        enddo
        diffusion_owners=owners;diffusion_geometry=geometry;spatial_geometry=extended
        local=old;lt=ot;ir=old_ir
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
          1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces, &
          b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,evaluate_diffusion=diffusion_source, &
      state_redist_max_order=requested_order)
        call require(good,'MPI callback transport SDC')
      end if
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'coupled SDC owner-local state')
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==reference%nodes(p)%state),'exact serial coupled SDC state')
        call require(all(lt%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'exact serial coupled SDC T')
        call require(all(ir%fields(p)%values==serial_ir%fields(p)%values),'exact serial coupled SDC IR')
      end do
      call require(limiter==slimiter,'exact coupled SDC limiter')
      error=max(maxval(abs(b-sb)/max(1._dp,abs(sb))),maxval(abs(w-sw)/max(1._dp,abs(sw))), &
        maxval(abs(ch-sch)/max(1._dp,abs(sch))),maxval(abs(ext-sext)/max(1._dp,abs(sext))))
      call require(error<2.e-12_dp,'serial coupled SDC inventories')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,final_inventory,MPI_COMM_WORLD,good)
      call require(good,'SDC final inventory')
      error=maxval(abs(final_inventory-initial_inventory-b-w-ch-ext)/max(1._dp,abs(initial_inventory)))
      if(rank==0) print *, 'MPI SDC layout, composite residual: ',layout,error
      call require(error<3.e-14_dp,'SDC composite closure')
      if(test_mode/='state_redist') call check_sdc_level(initial,geometry,extended,owners,old,ot)
      do failure=1,merge(5,3,test_mode=="reactor")
        local=old;lt=ot;ir=old_ir;b=17;w=19;ch=23;ext=29;limiter=31
        failed_order=requested_order
        if(failure==merge(5,3,test_mode=='reactor').and.rank==nranks-1) &
          failed_order=merge(4,2,requested_order==2)
        viscosity_coefficient=.1_dp
        if(failure==1) viscosity_coefficient=-1
        if(failure==2.and.rank==nranks-1) force(1)=force(1)+1
        reactor_calls=0;reactor_failure=failure==3
        if(failure==4.and.rank==nranks-1) then
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
          1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
          b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,difmag=viscosity_coefficient,redistribution=method, &
      state_redist_max_order=failed_order)
        else
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
          1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
          b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,difmag=viscosity_coefficient,redistribution=method, &
          evaluate_reaction=reactor_adapter, &
      state_redist_max_order=failed_order)
        end if
        if(failure==4) call require(reactor_calls==0,'reactor mismatch rejected before cells')
        if(failure==merge(5,3,test_mode=='reactor')) &
          call require(reactor_calls==0,'rank-order mismatch rejected before reactor calls')
        reactor_failure=.false.

        force=[2.e4_dp,-3.e4_dp,4.e4_dp]
        call require(.not.good,'collective coupled SDC rejection')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==old%fields(p)%values),'coupled SDC state rollback')
          call require(all(lt%fields(p)%values==ot%fields(p)%values),'coupled SDC T rollback')
          call require(all(ir%fields(p)%values==old_ir%fields(p)%values),'coupled SDC IR rollback')
        end do
        call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.limiter==31, &
          'coupled SDC inventory rollback')
      end do
      if(test_mode=='diffusion') then
        block
          type(amr_tree_state_3d) :: expected,expected_t,expected_ir
          real(dp) :: eb(nv),ew(nv),ech(nv),ee(nv),et
          local=old;lt=ot;ir=old_ir
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
            1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
            eb,ew,ech,ee,et,MPI_COMM_WORLD,good,ir, &
      state_redist_max_order=requested_order)
          call require(good,'native recursive transport reference')
          expected=local;expected_t=lt;expected_ir=ir
          local=old;lt=ot;ir=old_ir
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
            1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
            b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,evaluate_diffusion=diffusion_source, &
      state_redist_max_order=requested_order)
          call require(good,'recursive diffusion callback')
          do p=1,count
            if(owners(p)/=rank) cycle
            call require(all(local%fields(p)%values==expected%fields(p)%values),'recursive callback state exact')
            call require(all(lt%fields(p)%values==expected_t%fields(p)%values),'recursive callback T exact')
            call require(all(ir%fields(p)%values==expected_ir%fields(p)%values),'recursive callback IR exact')
          end do
          call require(all(b==eb).and.all(w==ew).and.all(ch==ech).and.all(ext==ee).and.limiter==et, &
            'recursive callback diagnostics exact')
          diffusion_failure=.true.;diffusion_calls=0
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
            1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
            b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,evaluate_diffusion=diffusion_source, &
      state_redist_max_order=requested_order)
          call require(.not.good,'rank-local callback failure rejects recursive step')
          diffusion_failure=.false.
          do p=1,count
            if(owners(p)/=rank) cycle
            call require(all(local%fields(p)%values==expected%fields(p)%values),'callback state rollback')
            call require(all(lt%fields(p)%values==expected_t%fields(p)%values),'callback T rollback')
            call require(all(ir%fields(p)%values==expected_ir%fields(p)%values),'callback IR rollback')
          end do
          call require(all(b==eb).and.all(w==ew).and.all(ch==ech).and.all(ext==ee).and.limiter==et, &
            'callback diagnostics rollback')
        end block
      end if
      deallocate(local%fields,lt%fields,ir%fields)
    end do
    if(test_mode=="reactor") call check_sdc_recursive(initial,geometry,extended,faces,force)
  end subroutine check_sdc_loop
  subroutine check_three_level_routed_replay(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(eb_reactive_tree_3d) :: nitrogen,reference
    type(amr_tree_state_3d) :: old,ot,local,lt,ir,heat,kept_ir,kept_heat,refusal_initial,refusal_temperature
    type(mpi_eb_mol_interval_3d), allocatable :: history(:),kept(:)
    type(eb_integrated_flux_3d), allocatable :: flux(:),registered(:),kept_flux(:),kept_registered(:)
    real(dp) :: dt
    real(dp) :: q(nv),v(nv),temp,cs,b(nv),w(nv),ch(nv),ext(nv),aux(nv),theta
    real(dp) :: initial_inventory(nv),final_inventory(nv),error,expected_aux(nv),weight,ref_ir(nv-4,4,4,4,5),ref_heat(1,4,4,4,5)
    integer :: owners(count),p,x,y,z,n2,water,layout,r,iterations,normal_calls,total_calls,accepted_calls,failure_case
    integer :: expected_refusal_calls,c
    logical :: good
    character(len=16) :: redistribution
    real(dp) :: transfer_norm,global_transfer_norm,replay_trace
    character(len=64) :: trace_argument
    integer :: peroxide,trace_status
    replay_trace=0;dt=1.e-8_dp
    if(test_mode=='state_replay') then
      replay_trace=1.e-10_dp;dt=1.e-6_dp
      call get_command_argument(3,trace_argument)
      if(len_trim(trace_argument)>0) then
        read(trace_argument,*,iostat=trace_status) replay_trace
        call require(trace_status==0.and.replay_trace>=0.and.replay_trace<1,'physical initial peroxide trace')
      endif
      call get_command_argument(4,trace_argument)
      if(len_trim(trace_argument)>0) then
        read(trace_argument,*,iostat=trace_status) dt
        call require(trace_status==0.and.dt>0,'physical replay interval')
      endif
    endif
    if(rank==0) print *, 'Initial replay peroxide mass fraction and interval: ',replay_trace,dt
    iterations=2
    if(test_mode=='routed_replay'.or.test_mode=='flux_replay'.or. &
      test_mode=='state_replay'.or.test_mode=='state_budget') iterations=1
    redistribution='none'
    if(test_mode=='flux_replay') redistribution='flux_redist'
    if(test_mode=='state_replay'.or.test_mode=='state_budget') redistribution='state_redist'
    n2=0;water=0;peroxide=0
    do p=1,size(species)
      if(trim(species(p)%name)=='N2') n2=p
      if(trim(species(p)%name)=='H2O') water=p
      if(trim(species(p)%name)=='H2O2') peroxide=p
    enddo
    call require(peroxide>0,'replay peroxide species')
    call require(n2>0.and.water>0,'routed replay nitrogen/water species')
    call require(covered(topology,1,[2,2,3]).and..not.covered(topology,2,[4,4,1]).and. &
      .not.covered(topology,2,[1,1,1]).and..not.covered(topology,1,[3,3,3]),'synthetic source predictor/leaf classification')
    weight=merge(1._dp,.5_dp,iterations==1)
    expected_aux=0
    expected_aux(1)=weight*dt*(geometry(4)%volume_fraction(1,1,1)*geometry(4)%dx*geometry(4)%dy*geometry(4)%dz+ &
      geometry(1)%volume_fraction(3,3,3)*geometry(1)%dx*geometry(1)%dy*geometry(1)%dz+ &
      (geometry(2)%volume_fraction(4,4,1)+geometry(2)%volume_fraction(1,1,1))* &
      geometry(2)%dx*geometry(2)%dy*geometry(2)%dz)
    expected_aux(5+water)=expected_aux(1)
    nitrogen=initial
    do p=1,count
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        q=0;q(1)=merge(.93_dp,.31_dp,geometry(p)%volume_fraction(x,y,z)<1)
        q(5)=135000;q(5+n2)=1-replay_trace;q(5+peroxide)=replay_trace
        call reactive_primitive_to_conserved(species,q,v,temp,cs,good)
        call require(good,'routed replay nitrogen state')
        nitrogen%nodes(p)%state(:,x,y,z)=v;nitrogen%nodes(p)%temperature(x,y,z)=temp
      enddo
      enddo
      enddo
    enddo
    call average_down_reactive_eb_tree_3d(species,nitrogen,good)
    call require(good,'routed replay nitrogen restriction')
    reference=nitrogen
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(old%fields(count),ot%fields(count),ir%fields(count),heat%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        old%fields(p)%values=nitrogen%nodes(p)%state
        allocate(ot%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4), &
          heat%fields(p)%values(1,4,4,4))
        ot%fields(p)%values(1,:,:,:)=nitrogen%nodes(p)%temperature
        ir%fields(p)%values=0;heat%fields(p)%values=0
      enddo
      local=old;lt=ot;history_recursive_calls=0;history_post_calls=0;history_failure=.false.
      b=0;w=0;ch=0;ext=0;aux=0;theta=1
      reactor_calls=0;reactor_failure=.false.;normal_calls=0
      do p=1,count
        normal_calls=normal_calls+iterations*2**topology%nodes(p)%level*count_fluid_cells(geometry(p))
      enddo
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,old,initial_inventory,MPI_COMM_WORLD,good)
      call require(good,'routed replay initial inventory')
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
        dt,iterations,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,mol_mode=.true.,auxiliary_source=replay_recursive_source, &
        auxiliary_change=aux,route_auxiliary=.true.,post_level=history_post_level,interval_history=history, &
        integrated_flux=flux,registered_flux=registered,heat_release=heat,evaluate_reaction=reactor_adapter, &
        redistribution=redistribution)
      if(rank==0) print *, 'Three-level routed replay acceptance: ',layout,good,theta,history_post_calls,history_recursive_calls
      call require(good,'three-level routed replay accepts conservative native solve')
      if(rank==0) print *, 'Three-level routed replay records: ',history%replayed
      call require(size(history)==7,'three-level accepted interval schedule')
      call MPI_Allreduce(reactor_calls,total_calls,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
      if(test_mode=='state_budget') then
        call require(.not.any(history%replayed).and.total_calls==normal_calls.and.theta<1, &
          'prepared local spatial budget admits absent peroxide without extra native reacting solves')
      else if(iterations==1) then
        if(test_mode=='flux_replay') then
          call require(all(history%replayed).and.total_calls>normal_calls, &
            'failed redistributed descendant escalates to native root-clock replay')
        else
          call require(all(history%replayed).and.total_calls>normal_calls, &
            'three-level root replay actually reruns native reaction on every interval')
        endif
      else
        call require(.not.any(history%replayed).and.total_calls==normal_calls, &
          'first/later local flux budgets admit the predictor case without interval replay')
      endif
      if(test_mode=='flux_replay') then
        transfer_norm=0
        do r=1,size(history)
          call require(history(r)%redistribution=='flux_redist'.and..not.history(r)%prepared_redistribution, &
            'replayed redistribution records retain actual operator family')
          do p=1,count
            if(owners(p)/=rank) cycle
            transfer_norm=transfer_norm+sum(abs(history(r)%redistribution_transfer%fields(p)%values))
            call require(all(history(r)%synchronization_transfer%fields(p)%values==0), &
              'canonical replay incorporates closure correction without duplicate spill')
          enddo
        enddo
        call MPI_Allreduce(transfer_norm,global_transfer_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
        call require(status==MPI_SUCCESS.and.global_transfer_norm>0,'native replay retains nonzero redistributed impulses')
      endif
      accepted_calls=total_calls
      call require(maxval(abs(aux-expected_aux)/max(1.e-30_dp,abs(expected_aux)))<3.e-14_dp, &
        'physical source inventory excludes virtual covered peroxide and uses accepted endpoint weights')
      if(test_mode=='state_replay'.or.test_mode=='state_budget') then
        call check_prepared_middle_closure(history,owners,geometry,extended,faces)
      else if(owners(2)==rank) then
        v=history(2)%local_state%fields(2)%values(:,1,1,1)+ &
          (history(3)%routed_source%fields(2)%values(:,1,1,1)+history(4)%routed_source%fields(2)%values(:,1,1,1))
        if(test_mode=='flux_replay') v=v+history(3)%redistribution_transfer%fields(2)%values(:,1,1,1)+ &
          history(4)%redistribution_transfer%fields(2)%values(:,1,1,1)
        call require(all(abs(history(5)%initial%fields(2)%values(:,1,1,1)-v)<= &
          8*epsilon(1._dp)*max(abs(v),abs(history(5)%initial%fields(2)%values(:,1,1,1)))), &
          'middle closure applies descendant impulses once before the next middle interval')
      endif
      call require(history_recursive_calls==14.and.history_post_calls==7,'replay retains frozen source and closure schedules')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,final_inventory,MPI_COMM_WORLD,good)
      call require(good,'routed replay final inventory')
      error=maxval(abs(final_inventory-initial_inventory-b-w-ch-ext-aux)/max(1._dp,abs(initial_inventory)))
      if(rank==0) print *, 'Three-level routed replay inventory error: ',error
      call require(error<3.e-14_dp,'routed replay composite gas conservation')
      do p=1,count
        if(layout==1) then
          if(owners(p)==rank) then
            reference%nodes(p)%state=local%fields(p)%values
            reference%nodes(p)%temperature=lt%fields(p)%values(1,:,:,:)
            ref_ir(:,:,:,:,p)=ir%fields(p)%values;ref_heat(:,:,:,:,p)=heat%fields(p)%values
          endif
          call MPI_Bcast(reference%nodes(p)%state,size(reference%nodes(p)%state),MPI_DOUBLE_PRECISION,owners(p), &
            MPI_COMM_WORLD,status)
          call MPI_Bcast(reference%nodes(p)%temperature,size(reference%nodes(p)%temperature),MPI_DOUBLE_PRECISION, &
            owners(p),MPI_COMM_WORLD,status)
          call MPI_Bcast(ref_ir(:,:,:,:,p),size(ref_ir(:,:,:,:,p)),MPI_DOUBLE_PRECISION,owners(p),MPI_COMM_WORLD,status)
          call MPI_Bcast(ref_heat(:,:,:,:,p),size(ref_heat(:,:,:,:,p)),MPI_DOUBLE_PRECISION,owners(p),MPI_COMM_WORLD,status)
        else if(owners(p)==rank) then
          call require(all(local%fields(p)%values==reference%nodes(p)%state),'routed replay gas exact across owners')
          call require(all(lt%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'routed replay T exact across owners')
          call require(all(ir%fields(p)%values==ref_ir(:,:,:,:,p)),'routed replay IR exact across owners')
          call require(all(heat%fields(p)%values==ref_heat(:,:,:,:,p)),'routed replay heat exact across owners')
        endif
      enddo
      ! Repeat the same initial state: final callback rejection follows the
      ! same replay, but must retain the previously accepted caller archive.
      kept=history;kept_flux=flux;kept_registered=registered
      local=old;lt=ot
      do p=1,count
        if(owners(p)/=rank) cycle
        ir%fields(p)%values=0;heat%fields(p)%values=0
      enddo
      kept_ir=ir;kept_heat=heat
      do failure_case=1,merge(3,merge(2,1,test_mode=='flux_replay'.or.test_mode=='state_replay'.or. &
        test_mode=='state_budget'),test_mode=='routed_replay'.or.replay_refusal_only)
      if(replay_refusal_only.and.failure_case/=3) cycle
      local=old;lt=ot;ir=kept_ir;heat=kept_heat;reactor_calls=0
      history_failure_depth=0;history_failure_iteration=0
      if(failure_case==2) then
        history_failure_depth=1;history_failure_iteration=2
      endif
      history_recursive_calls=0;history_post_calls=0;history_failure=.true.
      if(failure_case==3) then
        history_failure=.false.;reactor_replay_failure=.true.
        reactor_failure_rank=owners(1);reactor_replay_duration=dt;reactor_replay_calls=0
        reactor_replay_refused=.false.;reactor_replay_refusal_actual=0
        ! The first high root input withdraws 1.05e-21 kg/m3 of water
        ! from zero stock. Use less water than that withdrawal and seed other
        ! reactive intermediates too; every chosen input must admit natively.
        do p=1,count
          if(owners(p)/=rank) cycle
          do z=1,4
          do y=1,4
          do x=1,4
            if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
            call reactive_conserved_to_primitive(species,local%fields(p)%values(:,x,y,z), &
              lt%fields(p)%values(1,x,y,z),q,temp,cs,good)
            call require(good,'refusal probe initial primitive')
            do c=1,size(species)
              if(c==n2.or.trim(species(c)%name)=='AR') cycle
              q(5+c)=1.e-23_dp
            enddo
            q(5+water)=1.e-22_dp;q(5+n2)=0
            q(5+n2)=1-sum(q(6:))
            call reactive_primitive_to_conserved(species,q,v,temp,cs,good)
            call require(good,'refusal probe reactant stock')
            local%fields(p)%values(:,x,y,z)=v;lt%fields(p)%values(1,x,y,z)=temp
          enddo
          enddo
          enddo
        enddo
        call average_down_mpi_eb_tree_3d(species,topology,owners,geometry,local,lt,MPI_COMM_WORLD,good)
        call require(good,'refusal probe conservative restriction')
        ! Admit the complete low root solve first, then refuse an admissible
        ! cell of a positive-weight trial. Its input must abort the family
        ! search, rather than silently replacing the supplied solver.
        reactor_replay_refusal_call=count_fluid_cells(geometry(1))+1
      endif
      refusal_initial=local;refusal_temperature=lt
      b=17;w=19;ch=23;ext=29;aux=41;theta=31
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
        dt,iterations,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,mol_mode=.true.,auxiliary_source=replay_recursive_source, &
        auxiliary_change=aux,route_auxiliary=.true.,post_level=history_post_level,interval_history=history, &
        integrated_flux=flux,registered_flux=registered,heat_release=heat,evaluate_reaction=reactor_adapter, &
        redistribution=redistribution)
      call MPI_Allreduce(reactor_calls,total_calls,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
      if(failure_case==3) then
        if(rank==reactor_failure_rank) print *, 'Reactor refusal coverage: ',layout,good,history_recursive_calls, &
          history_post_calls,reactor_replay_calls,reactor_replay_refused,reactor_replay_refusal_actual
        call require(.not.good.and.history_recursive_calls==14,'replayed supplied-reactor refusal aborts collectively')
        if(rank==reactor_failure_rank) call require(reactor_replay_refused.and. &
          reactor_replay_calls==reactor_replay_refusal_actual.and. &
          reactor_replay_refusal_actual>=reactor_replay_refusal_call, &
          'native-admissible positive-weight solver refusal aborts without further family trials')
      else
        expected_refusal_calls=accepted_calls
        ! The unredistributed middle closure succeeds before the later root
        ! reflux requests replay. Refusing its post callback stops at the
        ! original seven native intervals, without entering that root replay.
        if(test_mode=='routed_replay'.and.failure_case==2) expected_refusal_calls=normal_calls
        if(rank==0) print *, 'Replay callback refusal schedule: ',failure_case,good, &
          history_post_calls,history_recursive_calls,total_calls,expected_refusal_calls
        call require(.not.good.and.history_post_calls==8-failure_case.and.history_recursive_calls==14.and. &
          total_calls==expected_refusal_calls, &
          'idle-rank root/deferred callback rejection never retries a failed callback or native group')
      endif
      call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.all(aux==41).and.theta==31, &
        'replayed rejection retains caller counters')
      call require(all(history%replayed.eqv.kept%replayed),'replayed rejection retains prior provenance')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==refusal_initial%fields(p)%values).and. &
          all(lt%fields(p)%values==refusal_temperature%fields(p)%values).and. &
          all(ir%fields(p)%values==kept_ir%fields(p)%values).and. &
          all(heat%fields(p)%values==kept_heat%fields(p)%values),'replayed rejection retains gas T IR heat')
        call require(all(flux(p)%x==kept_flux(p)%x).and.all(flux(p)%y==kept_flux(p)%y).and. &
          all(flux(p)%z==kept_flux(p)%z),'replayed rejection retains physical flux')
        call require(all(registered(p)%x==kept_registered(p)%x).and.all(registered(p)%y==kept_registered(p)%y).and. &
          all(registered(p)%z==kept_registered(p)%z),'replayed rejection retains registered flux')
        do r=1,size(history)
          call require(all(history(r)%routed_source%fields(p)%values==kept(r)%routed_source%fields(p)%values), &
            'replayed rejection retains frozen routed impulse')
          if(test_mode=='flux_replay') then
            call require(all(history(r)%redistribution_transfer%fields(p)%values== &
              kept(r)%redistribution_transfer%fields(p)%values).and. &
              all(history(r)%synchronization_transfer%fields(p)%values==kept(r)%synchronization_transfer%fields(p)%values), &
              'replayed rejection retains redistributed source and closure transfers')
          endif
          if(topology%nodes(p)%level/=history(r)%depth) cycle
          call require(all(history(r)%forcing%fields(p)%values==kept(r)%forcing%fields(p)%values), &
            'replayed rejection retains accepted forcing')
        enddo
      enddo
      enddo
      history_failure=.false.;history_failure_depth=0;history_failure_iteration=0
      reactor_replay_failure=.false.;reactor_failure_rank=0
      deallocate(old%fields,ot%fields,ir%fields,heat%fields)
    enddo
  end subroutine check_three_level_routed_replay

  subroutine check_prepared_middle_closure(history,owners,geometry,extended,faces)
    type(mpi_eb_mol_interval_3d), intent(in) :: history(:)
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d) :: gas,temp
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    type(state_redist_register_3d), allocatable :: registers(:)
    type(state_redist_register_3d) :: part
    real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed,error
    integer :: p,r,x,y,z,n(3),level
    logical :: good
    gas=history(2)%clock_initial;temp=history(2)%clock_temperature
    allocate(flux(count),registers(count))
    do p=1,count
      if(owners(p)/=rank) cycle
      n=topology%nodes(p)%shape;level=topology%nodes(p)%level
      call initialize_state_redist_tree_register_3d(topology,p,nv,6,registers(p),good,faces,.5_dp)
      call require(good,'independent middle clock register masks')
      allocate(flux(p)%x(nv,0:n(1),n(2),n(3)),flux(p)%y(nv,n(1),0:n(2),n(3)), &
        flux(p)%z(nv,n(1),n(2),0:n(3)),source=0._dp)
      if(level==1) then
        gas%fields(p)=history(2)%local_state%fields(p);temp%fields(p)=history(2)%local_temperature%fields(p)
        flux(p)=history(2)%spatial_flux(p)
        part=history(2)%redist_registers(p)
        if(allocated(part%fine_delta)) part%fine_delta=0
        call accumulate_state_redist_register_3d(registers(p),part,nv,n+12,good)
        call require(good,'independent middle coarse register contribution')
      else if(level==2) then
        gas%fields(p)=history(4)%local_state%fields(p);temp%fields(p)=history(4)%local_temperature%fields(p)
        do r=3,4
          flux(p)%x=flux(p)%x+history(r)%synchronized_flux(p)%x
          flux(p)%y=flux(p)%y+history(r)%synchronized_flux(p)%y
          flux(p)%z=flux(p)%z+history(r)%synchronized_flux(p)%z
          part=history(r)%redist_registers(p)
          if(allocated(part%coarse_delta)) part%coarse_delta=0
          call accumulate_state_redist_register_3d(registers(p),part,nv,n+12,good)
          call require(good,'independent middle fine register contribution')
        enddo
      endif
    enddo
    call synchronize_mpi_eb_tree_register_flux_3d(species,topology,owners,geometry,extended,gas,temp,flux, &
      registers,1,MPI_COMM_WORLD,good,faces)
    call require(good,'independent native middle StateRedist closure')
    do p=1,count
      if(owners(p)/=rank) cycle
      if(topology%nodes(p)%level>=1) then
        do r=2,4
          gas%fields(p)%values=gas%fields(p)%values+history(r)%routed_source%fields(p)%values
        enddo
      endif
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        call reactive_conserved_to_primitive(species,gas%fields(p)%values(:,x,y,z),temp%fields(p)%values(1,x,y,z), &
          primitive,recovered,speed,good)
        call require(good,'independent middle routed source EOS')
        temp%fields(p)%values(1,x,y,z)=recovered
      enddo
      enddo
      enddo
    enddo
    call average_down_mpi_eb_tree_3d(species,topology,owners,geometry,gas,temp,MPI_COMM_WORLD,good)
    call require(good,'independent middle routed closure restriction')
    do p=1,count
      if(owners(p)/=rank.or.topology%nodes(p)%level/=1) cycle
      error=maxval(abs(gas%fields(p)%values-history(5)%initial%fields(p)%values)/ &
        max(1._dp,abs(gas%fields(p)%values)))
      call require(error<3.e-14_dp,'middle native closure and routed impulses precede next interval')
    enddo
  end subroutine check_prepared_middle_closure

  integer function count_fluid_cells(g) result(n)
    type(eb_geometry_3d), intent(in) :: g
    n=sum(merge(1,0,g%volume_fraction>0))
  end function count_fluid_cells

  subroutine replay_recursive_source(state,temperature,interval,stage,depth,time,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p,water,peroxide
    real(dp) :: amount
    call require(stage==mod(history_recursive_calls,2),'replay source evaluates frozen old/first-new families only')
    call require(interval>0.and.time>=0.and.depth>=0.and.depth<=2,'replay source clocks')
    call require(allocated(temperature%fields).and.comm==MPI_COMM_WORLD,'replay source context')
    history_recursive_calls=history_recursive_calls+1
    water=0;peroxide=0
    do p=1,size(species)
      if(trim(species(p)%name)=='H2O') water=p
      if(trim(species(p)%name)=='H2O2') peroxide=p
    enddo
    call require(water>0.and.peroxide>0,'replay water and peroxide species')
    do p=1,count
      if(.not.allocated(state%fields(p)%values)) cycle
      rate%fields(p)%values=0
      if(depth==0.and.p==1) then
        ! Covered predictor peroxide has no composite inventory. The adjacent
        ! exposed hot cell consumes its provisional inflow before fine reflux
        ! removes it. This synthetic source isolates the native reaction/reflux
        ! interaction without depending on tiny water dissociation intermediates.
        rate%fields(p)%values(1,2,2,3)=1000
        rate%fields(p)%values(5+peroxide,2,2,3)=1000
      endif
      if(depth==1.and.p==2) then
        rate%fields(p)%values(1,4,4,1)=1
        rate%fields(p)%values(5+water,4,4,1)=1
      endif
      amount=1._dp
      if(depth==0.and.p==4) then
        rate%fields(p)%values(1,1,1,1)=amount
        rate%fields(p)%values(5+water,1,1,1)=amount
      endif
      if(depth==1.and.p==1) then
        rate%fields(p)%values(1,3,3,3)=amount
        rate%fields(p)%values(5+water,3,3,3)=amount
      endif
      if(depth==2.and.p==2) then
        rate%fields(p)%values(1,1,1,1)=amount
        rate%fields(p)%values(5+water,1,1,1)=amount
      endif
    enddo
    good=.true.
  end subroutine replay_recursive_source

  subroutine check_routed_mol_history(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(eb_reactive_tree_3d) :: reference
    type(amr_tree_state_3d) :: old,ot,ir,local,lt,profile,serial_ir
    type(mpi_eb_mol_interval_3d), allocatable :: history(:),kept(:)
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp), parameter :: dt=1.e-10_dp
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(nv),theta,sb(nv),sw(nv),sch(nv),sext(nv),saux(nv),stheta
    real(dp) :: expected_aux(nv),initial_inventory(nv),final_inventory(nv),weight,error,expected_impulse(nv),mpi_auxiliary(nv)
    integer :: owners(count),p,r,layout,x,y,z,iterations
    logical :: good
    iterations=2
    if(test_mode=='route_first') iterations=1
    owners=[(mod(p-1,nranks),p=1,count)]
    allocate(profile%fields(count),serial_ir%fields(count))
    do p=1,count
      allocate(serial_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
      if(owners(p)/=rank) cycle
      allocate(profile%fields(p)%values(nv,4,4,4),source=0._dp)
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        profile%fields(p)%values(:,x,y,z)=source_profile(p,x,y,z,0)
      enddo
      enddo
      enddo
    enddo
    call integrate_mpi_eb_tree_3d(topology,owners,geometry,profile,expected_aux,MPI_COMM_WORLD,good)
    call require(good,'independent routed source profile inventory')
    weight=merge(1.5_dp,.5_dp,iterations==1)
    ! Every physical leaf receives the accepted source family from each of
    ! three source depths, each spanning the same outer duration in total.
    expected_aux=3*weight*dt*expected_aux
    reference=initial;history_recursive_calls=0
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,iterations,.true., &
      [.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,stheta,good,serial_ir,mol_mode=.true.,auxiliary_source=history_serial_source, &
      auxiliary_change=saux,route_auxiliary=.true.)
    call require(good.and.history_recursive_calls==14,'serial routed MOL interval schedule')
    error=maxval(abs(saux-expected_aux)/max(1.e-30_dp,abs(expected_aux)))
    if(rank==0) print *, 'Routed MOL accepted source inventory error: ',iterations,error
    call require(error<3.e-14_dp,'routed MOL inventory uses the same first/later accepted source families')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(old%fields(count),ot%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        old%fields(p)%values=initial%nodes(p)%state
        allocate(ot%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4))
        ot%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature;ir%fields(p)%values=0
      enddo
      local=old;lt=ot;history_recursive_calls=0;history_post_calls=0;history_failure=.false.
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,old,initial_inventory,MPI_COMM_WORLD,good)
      call require(good,'routed MOL initial inventory')
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
        dt,iterations,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,mol_mode=.true.,auxiliary_source=history_recursive_source, &
        auxiliary_change=aux,route_auxiliary=.true.,post_level=history_post_level,interval_history=history, &
        integrated_flux=flux)
      call require(good.and.history_recursive_calls==14.and.history_post_calls==7,'MPI routed MOL interval schedule')
      call require(maxval(abs(aux-expected_aux)/max(1.e-30_dp,abs(expected_aux)))<3.e-14_dp, &
        'MPI routed MOL accepted source inventory')
      ! Serial and MPI inventory integrations have different addition groups;
      ! use the existing ledger criterion for that comparison. MPI owner changes
      ! retain the deterministic tree-order integration and must match exactly.
      error=maxval(abs(aux-saux)/max(1.e-30_dp,abs(saux)))
      if(rank==0) print *, 'Routed MOL serial/MPI ledger error: ',layout,error
      call require(error<3.e-14_dp,'routed source inventory matches serial within the existing ledger criterion')
      if(layout==1) then
        mpi_auxiliary=aux
      else
        call require(all(aux==mpi_auxiliary),'routed source inventory remains exact across MPI owner layouts')
      endif
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==reference%nodes(p)%state),'routed MOL serial MPI gas parity')
        call require(all(lt%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'routed MOL serial MPI T parity')
        call require(all(ir%fields(p)%values==serial_ir%fields(p)%values),'routed MOL serial MPI IR parity')
      enddo
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,final_inventory,MPI_COMM_WORLD,good)
      call require(good,'routed MOL final inventory')
      error=maxval(abs(final_inventory-initial_inventory-b-w-ch-ext-aux)/max(1._dp,abs(initial_inventory)))
      call require(error<3.e-14_dp,'routed MOL sources enter gas once with composite conservation')
      do r=1,size(history)
        call require(allocated(history(r)%routed_source%fields),'accepted routed interval retains frozen off-level impulse')
        do p=1,count
          call require(allocated(history(r)%routed_source%fields(p)%values).eqv.(owners(p)==rank), &
            'routed interval impulse remains owner-local with idle ranks')
          if(owners(p)/=rank) cycle
          do z=1,4
          do y=1,4
          do x=1,4
            expected_impulse=0
            if(topology%nodes(p)%level/=history(r)%depth.and.geometry(p)%volume_fraction(x,y,z)>0) then
              if(.not.covered(topology,p,[x,y,z])) &
                expected_impulse=weight*history(r)%interval*source_profile(p,x,y,z,0)
            endif
            call require(all(abs(history(r)%routed_source%fields(p)%values(:,x,y,z)-expected_impulse)<= &
              8*epsilon(1._dp)*max(abs(expected_impulse),abs(history(r)%routed_source%fields(p)%values(:,x,y,z)))), &
              'frozen routed impulse uses accepted weights and excludes source level covered and solid cells')
          enddo
          enddo
          enddo
        enddo
      enddo
      kept=history
      old=local;ot=lt;profile=ir
      history_recursive_calls=0;history_post_calls=0;history_failure=.true.
      b=17;w=19;ch=23;ext=29;aux=41;theta=31
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
        dt,iterations,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,mol_mode=.true.,auxiliary_source=history_recursive_source, &
        auxiliary_change=aux,route_auxiliary=.true.,post_level=history_post_level,interval_history=history)
      call require(.not.good.and.history_post_calls==7,'idle-rank final closure failure rejects routed MOL transaction')
      call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.all(aux==41).and.theta==31, &
        'failed routed MOL transaction retains all inventory counters')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==old%fields(p)%values).and. &
          all(lt%fields(p)%values==ot%fields(p)%values).and.all(ir%fields(p)%values==profile%fields(p)%values), &
          'failed routed MOL transaction retains gas T and IR')
        do r=1,size(history)
          call require(all(history(r)%routed_source%fields(p)%values==kept(r)%routed_source%fields(p)%values), &
            'failed routed MOL transaction retains the prior off-level impulse')
          if(topology%nodes(p)%level/=history(r)%depth) cycle
          call require(all(history(r)%forcing%fields(p)%values==kept(r)%forcing%fields(p)%values), &
            'failed routed MOL transaction retains prior accepted interval forcing')
        enddo
      enddo
      history_failure=.false.
      deallocate(old%fields,ot%fields,ir%fields)
    enddo
  end subroutine check_routed_mol_history

  subroutine history_serial_source(state,interval,stage,depth,time,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p,x,y,z
    call require(stage==mod(history_recursive_calls,2),'serial MOL source family order')
    call require(depth>=0.and.depth<=2.and.interval>0.and.time>=0,'serial MOL source clocks')
    history_recursive_calls=history_recursive_calls+1
    do p=1,size(state%nodes)
      do z=1,4
      do y=1,4
      do x=1,4
        rate%fields(p)%values(:,x,y,z)=source_profile(p,x,y,z,stage)
      enddo
      enddo
      enddo
    enddo
    good=.true.
  end subroutine history_serial_source

  subroutine check_recursive_interval_history(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(amr_tree_state_3d) :: old,ot,zero_ir,local,lt,ir,baseline,baseline_t,baseline_ir
    type(amr_tree_state_3d) :: replay,replay_t,replay_ir,forcing,replay_old,replay_old_t
    type(mpi_eb_mol_interval_3d), allocatable :: history(:),kept(:)
    type(eb_integrated_flux_3d), allocatable :: flux(:),baseline_flux(:)
    integer, parameter :: depths(7)=[0,1,2,2,1,2,2]
    real(dp), parameter :: dt=1.e-10_dp,times(7)=[0._dp,0._dp,0._dp,.25_dp,.5_dp,.5_dp,.75_dp]*dt
    integer(int64) :: steps(3)
    integer :: owners(count),layout,p,r,x,y,z
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(nv),theta
    real(dp) :: bb(nv),bw(nv),bch(nv),be(nv),ba(nv),bt,record_chemical(nv),delta(nv),expected(nv)
    type(eb_integrated_flux_3d), allocatable :: changed_flux(:)
    logical :: good,active
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(old%fields(count),ot%fields(count),zero_ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        old%fields(p)%values=initial%nodes(p)%state
        allocate(ot%fields(p)%values(1,4,4,4),zero_ir%fields(p)%values(nv-4,4,4,4))
        ot%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature;zero_ir%fields(p)%values=0
      enddo
      baseline=old;baseline_t=ot;baseline_ir=zero_ir
      history_failure=.false.;history_post_calls=0;history_recursive_calls=0
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,baseline,baseline_t, &
        dt,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        bb,bw,bch,be,bt,MPI_COMM_WORLD,good,baseline_ir,mol_mode=.true., &
        auxiliary_source=history_recursive_source,auxiliary_change=ba,post_level=history_post_level, &
        integrated_flux=baseline_flux)
      call require(good.and.history_post_calls==7.and.history_recursive_calls==14,'baseline recursive history staging')
      local=old;lt=ot;ir=zero_ir;history_post_calls=0;history_recursive_calls=0
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
        dt,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,mol_mode=.true.,level_steps=steps, &
        auxiliary_source=history_recursive_source,auxiliary_change=aux,post_level=history_post_level, &
        integrated_flux=flux,interval_history=history)
      call require(good.and.size(history)==7,'all recursive history intervals accepted')
      call require(history_post_calls==7.and.history_recursive_calls==14,'history does not rebuild physical sources')
      call require(all(steps==[1_int64,2_int64,4_int64]),'recursive history level interval counts')
      call require(all(b==bb).and.all(w==bw).and.all(ch==bch).and.all(ext==be).and.all(aux==ba).and.theta==bt, &
        'recursive history leaves accepted ledgers exactly unchanged')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==baseline%fields(p)%values).and. &
          all(lt%fields(p)%values==baseline_t%fields(p)%values).and. &
          all(ir%fields(p)%values==baseline_ir%fields(p)%values),'recording leaves accepted gas exactly unchanged')
        call require(all(flux(p)%x==baseline_flux(p)%x).and.all(flux(p)%y==baseline_flux(p)%y).and. &
          all(flux(p)%z==baseline_flux(p)%z),'recording leaves accepted physical flux exactly unchanged')
      enddo
      record_chemical=0
      do r=1,7
        call require(history(r)%depth==depths(r).and.history(r)%time==times(r).and. &
          history(r)%interval==dt/2**depths(r),'recursive history retains parent-first subcycle clocks')
        record_chemical=record_chemical+history(r)%chemical_change
        replay=old;replay_t=ot;replay_ir=zero_ir;forcing=old
        do p=1,count
          active=owners(p)==rank.and.topology%nodes(p)%level==depths(r)
          call require(allocated(history(r)%initial%fields(p)%values).eqv.active,'sparse initial history ownership')
          call require(allocated(history(r)%initial_temperature%fields(p)%values).eqv.active,'sparse T history ownership')
          call require(allocated(history(r)%forcing%fields(p)%values).eqv.active,'sparse forcing history ownership')
          call require(allocated(history(r)%mass_source%fields(p)%values).eqv.active,'sparse mass history ownership')
          call require(allocated(history(r)%local_state%fields(p)%values).eqv.active,'sparse local endpoint ownership')
          call require(allocated(history(r)%local_temperature%fields(p)%values).eqv.active,'sparse local T ownership')
          call require(allocated(history(r)%reaction_source%fields(p)%values).eqv.active,'sparse IR ownership')
          call require(allocated(history(r)%old_flux(p)%x).eqv.active,'sparse old spatial flux ownership')
          call require(allocated(history(r)%spatial_flux(p)%x).eqv.active,'sparse local spatial flux ownership')
          call require(allocated(history(r)%synchronized_flux(p)%x).eqv.active,'sparse closed spatial flux ownership')
          if(owners(p)/=rank) cycle
          forcing%fields(p)%values=0
          if(.not.active) cycle
          replay%fields(p)=history(r)%initial%fields(p)
          replay_t%fields(p)=history(r)%initial_temperature%fields(p)
          forcing%fields(p)=history(r)%forcing%fields(p)
          call require(lbound(history(r)%old_flux(p)%x,2)==0,'history preserves zero-based face bounds')
          call require(all(history(r)%mass_source%fields(p)%values(2:5,:,:,:)==0),'mass history excludes bulk sources')
          if(r==1) call require(all(history(r)%synchronized_flux(p)%x==flux(p)%x).and. &
            all(history(r)%synchronized_flux(p)%y==flux(p)%y).and. &
            all(history(r)%synchronized_flux(p)%z==flux(p)%z),'root history uses the final synchronized flux')
        enddo
        call reconstruct_mpi_eb_mol_interval_3d(topology,owners,geometry,history(r),history(r)%spatial_flux, &
          nv,forcing,MPI_COMM_WORLD,good)
        call require(good,'unchanged interval flux reconstructs the complete reaction forcing')
        do p=1,count
          if(owners(p)/=rank) cycle
          if(topology%nodes(p)%level/=depths(r)) then
            call require(all(forcing%fields(p)%values==0),'reconstructed rates exclude inactive levels')
          else
            call require(all(forcing%fields(p)%values==history(r)%forcing%fields(p)%values), &
              'unchanged flux retains wall and volume forcing exactly')
          endif
        enddo
        replay_old=replay;replay_old_t=replay_t
        call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,replay_old,replay_old_t,forcing, &
          history(r)%interval,2.e-7_dp,1.e-12_dp,replay,replay_t,replay_ir,MPI_COMM_WORLD,good,target_depth=depths(r))
        call require(good,'recorded recursive forcing admits exact local reaction replay')
        do p=1,count
          if(owners(p)/=rank.or.topology%nodes(p)%level/=depths(r)) cycle
          call require(all(replay%fields(p)%values==history(r)%local_state%fields(p)%values).and. &
            all(replay_t%fields(p)%values==history(r)%local_temperature%fields(p)%values).and. &
            all(replay_ir%fields(p)%values==history(r)%reaction_source%fields(p)%values), &
            'every recursive interval replays state temperature and chemistry exactly')
        enddo
      enddo
      ! A changed internal open face contributes equal/opposite integrated
      ! inventories to its two fluid cells without changing any volume source.
      changed_flux=history(3)%spatial_flux
      do p=1,count
        if(owners(p)/=rank.or.topology%nodes(p)%level/=2) cycle
        call require(all(geometry(p)%volume_fraction(1:2,2,4)>0),'test internal face joins two fluid cells')
        changed_flux(p)%x(1,1,2,4)=changed_flux(p)%x(1,1,2,4)+1.e-14_dp
        changed_flux(p)%x(6,1,2,4)=changed_flux(p)%x(6,1,2,4)+1.e-14_dp
      enddo
      call reconstruct_mpi_eb_mol_interval_3d(topology,owners,geometry,history(3),changed_flux,nv, &
        forcing,MPI_COMM_WORLD,good)
      call require(good,'changed face reconstructs interval forcing transactionally')
      do p=1,count
        if(owners(p)/=rank.or.topology%nodes(p)%level/=2) cycle
        delta=geometry(p)%x_face_fraction(1,2,4)* &
          (changed_flux(p)%x(:,1,2,4)-history(3)%spatial_flux(p)%x(:,1,2,4))/geometry(p)%dx/history(3)%interval
        do z=1,4
        do y=1,4
        do x=1,4
          expected=history(3)%forcing%fields(p)%values(:,x,y,z)
          if(y==2.and.z==4.and.x==1) expected=expected-delta/geometry(p)%volume_fraction(x,y,z)
          if(y==2.and.z==4.and.x==2) expected=expected+delta/geometry(p)%volume_fraction(x,y,z)
          call require(maxval(abs(forcing%fields(p)%values(:,x,y,z)-expected)/max(1._dp,abs(expected)))<3.e-14_dp, &
            'face change contributes only signed open-face divergence divided by fluid volume and interval')
        enddo
        enddo
        enddo
      enddo
      replay_old=forcing
      do p=1,count
        if(owners(p)/=rank.or.topology%nodes(p)%level/=2) cycle
        changed_flux(p)%x(1,1,2,4)=ieee_value(0._dp,ieee_quiet_nan)
      enddo
      call reconstruct_mpi_eb_mol_interval_3d(topology,owners,geometry,history(3),changed_flux,nv, &
        forcing,MPI_COMM_WORLD,good)
      call require(.not.good,'rank-local nonfinite face rejects reconstruction on idle ranks too')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(forcing%fields(p)%values==replay_old%fields(p)%values),'failed reconstruction retains prior rates')
      enddo
      call require(all(record_chemical==ch),'preorder reaction inventory matches the accepted outer ledger')
      kept=history
      history_failure=.true.;history_post_calls=0;history_recursive_calls=0
      steps=37;b=17;w=19;ch=23;ext=29;aux=41;theta=31
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
        dt,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,mol_mode=.true.,level_steps=steps, &
        auxiliary_source=history_recursive_source,auxiliary_change=aux,post_level=history_post_level, &
        integrated_flux=flux,interval_history=history)
      call require(.not.good.and.history_post_calls==7,'last root closure failure rejects the entire history transaction')
      call require(all(steps==37).and.all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and. &
        all(aux==41).and.theta==31,'late recursive history failure rolls back counters and ledgers')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==baseline%fields(p)%values).and. &
          all(lt%fields(p)%values==baseline_t%fields(p)%values).and. &
          all(ir%fields(p)%values==baseline_ir%fields(p)%values),'late history failure rolls back gas and reaction fields')
        call require(all(flux(p)%x==baseline_flux(p)%x).and.all(flux(p)%y==baseline_flux(p)%y).and. &
          all(flux(p)%z==baseline_flux(p)%z),'late history failure retains final flux')
        do r=1,7
          if(topology%nodes(p)%level/=depths(r)) cycle
          call require(all(history(r)%initial%fields(p)%values==kept(r)%initial%fields(p)%values).and. &
            all(history(r)%forcing%fields(p)%values==kept(r)%forcing%fields(p)%values).and. &
            all(history(r)%local_state%fields(p)%values==kept(r)%local_state%fields(p)%values), &
            'late failure retains prior accepted interval fields')
          call require(all(history(r)%synchronized_flux(p)%x==kept(r)%synchronized_flux(p)%x).and. &
            all(history(r)%synchronized_flux(p)%y==kept(r)%synchronized_flux(p)%y).and. &
            all(history(r)%synchronized_flux(p)%z==kept(r)%synchronized_flux(p)%z),'late failure retains prior closed flux')
        enddo
      enddo
      history_failure=.false.
      deallocate(old%fields,ot%fields,zero_ir%fields)
    enddo
  end subroutine check_recursive_interval_history

  subroutine history_recursive_source(state,temperature,interval,stage,depth,time,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p,x,y,z
    call require(stage==mod(history_recursive_calls,2),'recursive history source uses only old and first-new families')
    call require(depth>=0.and.depth<=2.and.interval>0.and.time>=0,'recursive history source clocks')
    call require(allocated(temperature%fields).and.comm==MPI_COMM_WORLD,'recursive history source context')
    history_recursive_calls=history_recursive_calls+1
    do p=1,count
      if(.not.allocated(state%fields(p)%values)) cycle
      do z=1,4
      do y=1,4
      do x=1,4
        rate%fields(p)%values(:,x,y,z)=source_profile(p,x,y,z,stage)
      enddo
      enddo
      enddo
    enddo
    good=.true.
  end subroutine history_recursive_source

  subroutine history_post_level(depth,time,interval,iteration,ncycle,comm,good)
    integer, intent(in) :: depth,iteration,ncycle
    real(dp), intent(in) :: time,interval
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    call require(depth>=0.and.depth<=2.and.interval>0.and.time>=0,'history closure clocks')
    call require(iteration>=1.and.iteration<=ncycle.and.comm==MPI_COMM_WORLD,'history closure context')
    history_post_calls=history_post_calls+1
    good=.not.(history_failure.and.depth==history_failure_depth.and.rank==nranks-1.and. &
      (history_failure_iteration==0.or.iteration==history_failure_iteration))
  end subroutine history_post_level

  subroutine check_redist_interval_history(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(amr_tree_state_3d) :: old,ot,zero_ir,local,lt,ir,base,base_t,base_ir,forcing,mass,transfer,reconstructed, &
      base_forcing,base_mass,base_transfer,changed_spatial,fresh_start
    type(mpi_eb_mol_interval_3d) :: history,kept,replayed_history
    type(mpi_eb_mol_interval_3d), allocatable :: optional_history
    type(state_redist_register_3d), allocatable :: registers(:),saved_registers(:),base_registers(:),replay_registers(:)
    type(eb_integrated_flux_3d), allocatable :: flux(:),reference_low_flux(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(nv),theta,expected(nv),error,correction,global_correction,weight
    integer :: owners(count),layout,method_index,pass,p,x,y,z,mpi_status,stage,blend_case
    character(len=16) :: method
    logical :: good,provider_refused
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(old%fields(count),ot%fields(count),zero_ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        old%fields(p)%values=initial%nodes(p)%state
        allocate(ot%fields(p)%values(1,4,4,4),zero_ir%fields(p)%values(nv-4,4,4,4))
        ot%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature;zero_ir%fields(p)%values=0
      enddo
      do method_index=1,2
        method='flux_redist'
        if(method_index==2) method='state_redist'
        do pass=1,2
          local=old;lt=ot;ir=zero_ir;history_calls=0;history_failure=.false.
          ! Optional archive capture must not alter the accepted operator or
          ! its transfer/register family. An unallocated scalar forwards an
          ! absent archive on the first comparison call.
          if(allocated(optional_history)) deallocate(optional_history)
          if(pass==2) allocate(optional_history)
          if(method=='state_redist') then
            call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
              1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
              MPI_COMM_WORLD,good,ir,target_depth=1,ancestor_start=old,ancestor_end=old,temperature_start=ot,temperature_end=ot, &
              fractions_start=[0._dp,0._dp],fractions_end=[1._dp,1._dp],auxiliary_source=history_auxiliary,auxiliary_change=aux, &
              mol_mode=.true.,redistribution=method,target_volume_fraction=.9_dp,state_redist_registers=registers, &
              accepted_forcing=forcing,accepted_mass_source=mass,integrated_flux=flux,interval_history=optional_history)
          else
            call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
              1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
              MPI_COMM_WORLD,good,ir,target_depth=1,ancestor_start=old,ancestor_end=old,temperature_start=ot,temperature_end=ot, &
              fractions_start=[0._dp,0._dp],fractions_end=[1._dp,1._dp],auxiliary_source=history_auxiliary,auxiliary_change=aux, &
              mol_mode=.true.,redistribution=method,target_volume_fraction=.9_dp,redistribution_transfer=transfer, &
              accepted_forcing=forcing,accepted_mass_source=mass,integrated_flux=flux,interval_history=optional_history)
          endif

          if(pass==1) then
            call require(good.and.history_calls==2,'redistributed MOL accepts without archive')
            base=local;base_t=lt;base_ir=ir;base_forcing=forcing;base_mass=mass
            if(method=='flux_redist') base_transfer=transfer
            if(method=='state_redist') base_registers=registers
          else
            if(good) history=optional_history

            call require(good.and.history_calls==2,'redistributed MOL accepts with archive')
            call require(history%redistribution==method.and.history%depth==1.and..not.history%replayed, &
              'archive identifies the actual redistribution family')
            call require(history%prepared_redistribution.eqv.(method=='state_redist'), &
              'archive distinguishes prepared register and unprepared transfer paths')
            correction=0
            do p=1,count
              if(owners(p)/=rank) cycle
              if(method=='flux_redist') then
                correction=correction+sum(abs(history%old_redistribution_transfer%fields(p)%values))+ &
                  sum(abs(history%redistribution_transfer%fields(p)%values))
              elseif(topology%nodes(p)%level==1) then
                correction=correction+sum(abs(history%old_redist_registers(p)%coarse_delta))+ &
                  sum(abs(history%old_redist_registers(p)%fine_delta))+ &
                  sum(abs(history%redist_registers(p)%coarse_delta))+sum(abs(history%redist_registers(p)%fine_delta))
              endif
            enddo
            call MPI_Allreduce(correction,global_correction,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_status)
            call require(mpi_status==MPI_SUCCESS.and.global_correction>0, &
              'archive exercises nonzero cross-level redistribution corrections')
            call reconstruct_mpi_eb_mol_interval_3d(topology,owners,geometry,history,history%spatial_flux,nv, &
              reconstructed,MPI_COMM_WORLD,good)
            call require(good,'redistributed archive reconstructs original complete forcing')
            if(method=='state_redist') then
              call require(allocated(history%prepared_stages),'prepared MOL retains native stage inputs')
              call require(lbound(history%prepared_stages,1)==0.and.ubound(history%prepared_stages,1)==2, &
                'old and every reacting new stage retained separately')
              do stage=0,2
                associate(saved=>history%prepared_stages(stage))
                call require(saved%source_stage==stage.and.saved%depth==1.and.saved%interval==1.e-10_dp.and. &
                  saved%target_volume_fraction==.9_dp.and.saved%register_scale==.5_dp.and.saved%max_order==0, &
                  'prepared stage retains exact operator controls')
                reconstructed=old
                do p=1,count
                  call require(allocated(saved%state_halos%fields(p)%values).eqv. &
                    (owners(p)==rank.and.topology%nodes(p)%level==1),'stage halos keep sparse owners and level')
                  if(owners(p)/=rank) cycle
                  reconstructed%fields(p)%values=0
                  if(topology%nodes(p)%level/=1) cycle
                  call require(all(shape(saved%state_halos%fields(p)%values)==[nv,16,16,16]).and. &
                    all(shape(saved%raw_rate_halos%fields(p)%values)==[nv,16,16,16]), &
                    'stage inputs retain all six halo layers')
                enddo
                if(allocated(replay_registers)) deallocate(replay_registers)
                call compute_state_redist_mpi_level_stage_3d(topology,owners,extended,saved%state_halos, &
                  saved%raw_rate_halos,nv,saved%depth,saved%interval,reconstructed,replay_registers,MPI_COMM_WORLD, &
                  good,faces,saved%target_volume_fraction,saved%max_order,register_scale=saved%register_scale)
                call require(good,'archived input reruns actual StateRedist geometry operator')
                do p=1,count
                  if(owners(p)/=rank.or.topology%nodes(p)%level/=1) cycle
                  call post_eb_spatial_density_rate_3d(geometry(p),reconstructed%fields(p)%values,good)
                  call require(good.and.all(reconstructed%fields(p)%values==saved%spatial_rate%fields(p)%values), &
                    'archived raw stage reproduces exact post-redistribution density and spatial RHS')
                  call require_register_equal(replay_registers(p),saved%registers(p))
                  if(stage==0) call require(all(saved%faces(p)%x==history%old_flux(p)%x).and. &
                    all(saved%faces(p)%y==history%old_flux(p)%y).and.all(saved%faces(p)%z==history%old_flux(p)%z), &
                    'old stage retains actual unweighted spatial faces')
                enddo
                end associate
              enddo
              ! The final reaction forcing uses the old and last new stages;
              ! intermediate coarse register contributions remain separate.
              do p=1,count
                if(owners(p)/=rank.or.topology%nodes(p)%level/=1) cycle
                call require(all(history%old_spatial_rate%fields(p)%values== &
                  history%prepared_stages(0)%spatial_rate%fields(p)%values).and. &
                  all(history%spatial_rate%fields(p)%values==.5_dp*( &
                  history%prepared_stages(0)%spatial_rate%fields(p)%values+ &
                  history%prepared_stages(2)%spatial_rate%fields(p)%values)), &
                  'accepted spatial forcing uses the archived old and final new native operators')
              enddo
              call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,history, &
                history%prepared_stages,nv,reconstructed,replay_registers,MPI_COMM_WORLD,good,faces)
              call require(good,'native stages reconstruct accepted prepared interval operator')
              do p=1,count
                if(owners(p)/=rank) cycle
                call require(all(reconstructed%fields(p)%values==forcing%fields(p)%values), &
                  're-evaluated stage mean reproduces complete accepted forcing exactly')
                if(topology%nodes(p)%level/=1) cycle
                call require_register_equal(replay_registers(p),history%redist_registers(p))
              enddo
              ! Changing a supplied final-stage raw source must invoke the
              ! native redistribution operator and compose its new mean RHS;
              ! an identity-only shortcut cannot satisfy this check.
              kept=history
              if(owners(2)==rank) kept%prepared_stages(2)%raw_rate_halos%fields(2)%values(2,7,7,7)= &
                kept%prepared_stages(2)%raw_rate_halos%fields(2)%values(2,7,7,7)+1._dp
              changed_spatial=old
              do p=1,count
                if(owners(p)==rank) changed_spatial%fields(p)%values=0
              enddo
              if(allocated(replay_registers)) deallocate(replay_registers)
              associate(saved=>kept%prepared_stages(2))
              call compute_state_redist_mpi_level_stage_3d(topology,owners,extended,saved%state_halos, &
                saved%raw_rate_halos,nv,saved%depth,saved%interval,changed_spatial,replay_registers,MPI_COMM_WORLD, &
                good,faces,saved%target_volume_fraction,saved%max_order,register_scale=saved%register_scale)
              end associate
              call require(good,'changed stage has an independently evaluated native spatial operator')
              do p=1,count
                if(owners(p)/=rank.or.topology%nodes(p)%level/=1) cycle
                call post_eb_spatial_density_rate_3d(geometry(p),changed_spatial%fields(p)%values,good)
                call require(good,'changed stage native density consistency')
              enddo
              call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,history, &
                kept%prepared_stages,nv,reconstructed,replay_registers,MPI_COMM_WORLD,good,faces)
              call require(good,'replaced raw stage rebuilds complete prepared interval operator')
              if(owners(2)==rank) then
                call require(any(reconstructed%fields(2)%values/=forcing%fields(2)%values), &
                  'changed native spatial source actually changes interval forcing')
                call require(all(reconstructed%fields(2)%values==history%volume_source%fields(2)%values+.5_dp*( &
                  history%prepared_stages(0)%spatial_rate%fields(2)%values+changed_spatial%fields(2)%values)), &
                  'changed interval uses native mean RHS and original physical volume forcing')
              endif
              ! Restore the accepted outputs before collective rejection.
              call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,history, &
                history%prepared_stages,nv,reconstructed,replay_registers,MPI_COMM_WORLD,good,faces)
              call require(good,'prepared reconstruction restores the exact accepted anchor')
              ! A refusal on the idle last rank in a later stage must leave
              ! the already accepted output forcing and every register intact.
              kept=history
              if(rank==nranks-1) kept%prepared_stages(2)%interval=2.e-10_dp
              call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,history, &
                kept%prepared_stages,nv,reconstructed,replay_registers,MPI_COMM_WORLD,good,faces)
              call require(.not.good,'idle-rank stage control mismatch rejects interval reconstruction collectively')
              do p=1,count
                if(owners(p)/=rank) cycle
                call require(all(reconstructed%fields(p)%values==forcing%fields(p)%values), &
                  'rejected prepared stage reconstruction retains complete output forcing')
                if(topology%nodes(p)%level/=1) cycle
                call require_register_equal(replay_registers(p),history%redist_registers(p))
              enddo
              ! Run the real MOL predictor/native chemistry loop again from
              ! fresh gas/halo states, with frozen volume families. The source
              ! provider would now reject if it were called a second time.
              local=old;lt=ot;ir=zero_ir
              do p=1,count
                if(owners(p)==rank.and.topology%nodes(p)%level==1) ir%fields(p)%values=7
              enddo
              history_calls=0;history_failure=.true.
              call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
                1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
                MPI_COMM_WORLD,good,ir,target_depth=1,ancestor_start=old,ancestor_end=old, &
                temperature_start=ot,temperature_end=ot,fractions_start=[0._dp,0._dp],fractions_end=[1._dp,1._dp], &
                auxiliary_source=history_auxiliary,auxiliary_change=aux,mol_mode=.true.,redistribution=method, &
                target_volume_fraction=.9_dp,state_redist_registers=registers,accepted_forcing=forcing, &
                accepted_mass_source=mass,integrated_flux=flux,interval_history=replayed_history,replay_interval=history)
              call require(good.and.history_calls==0,'prepared native replay never re-executes frozen source provider')
              call require(replayed_history%has_auxiliary_source.and.replayed_history%replayed, &
                'native replay retains frozen source identity and actual rerun provenance')
              do p=1,count
                if(owners(p)/=rank) cycle
                call require(all(local%fields(p)%values==base%fields(p)%values).and. &
                  all(lt%fields(p)%values==base_t%fields(p)%values).and.all(ir%fields(p)%values==base_ir%fields(p)%values), &
                  'fresh native replay restores archived I_R predictor and reproduces gas T IR exactly')
                call require(all(forcing%fields(p)%values==base_forcing%fields(p)%values).and. &
                  all(mass%fields(p)%values==base_mass%fields(p)%values),'native replay forcing families exact')
                if(topology%nodes(p)%level/=1) cycle
                call require_register_equal(registers(p),base_registers(p))
                call require(all(replayed_history%initial_reaction_source%fields(p)%values== &
                  history%initial_reaction_source%fields(p)%values),'native replay records the actual archived initial I_R')
                do stage=0,2
                  call require(all(replayed_history%prepared_stages(stage)%state_halos%fields(p)%values== &
                    history%prepared_stages(stage)%state_halos%fields(p)%values).and. &
                    all(replayed_history%prepared_stages(stage)%raw_rate_halos%fields(p)%values== &
                    history%prepared_stages(stage)%raw_rate_halos%fields(p)%values), &
                    'native replay rebuilds every fresh stage halo and raw RHS exactly')
                enddo
              enddo
              ! An owner-local forced-reactor rejection must abort publication,
              ! even though earlier stages have rebuilt fresh private halos.
              local=old;lt=ot;ir=zero_ir;reactor_calls=0
              reactor_failure=.true.;reactor_failure_rank=owners(2)
              b=17;w=19;ch=23;ext=29;aux=41;theta=31
              call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
                1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
                MPI_COMM_WORLD,good,ir,target_depth=1,ancestor_start=old,ancestor_end=old, &
                temperature_start=ot,temperature_end=ot,fractions_start=[0._dp,0._dp],fractions_end=[1._dp,1._dp], &
                auxiliary_source=history_auxiliary,auxiliary_change=aux,mol_mode=.true.,redistribution=method, &
                target_volume_fraction=.9_dp,state_redist_registers=registers,accepted_forcing=forcing, &
                accepted_mass_source=mass,integrated_flux=flux,interval_history=replayed_history,replay_interval=history, &
                evaluate_reaction=reactor_adapter,provider_refused=provider_refused)
              call require(.not.good.and.history_calls==0,'native replay reactor failure rejects without source callbacks')
              call require(provider_refused,'native-admissible reactor provider refusal is collective and aborts trials')
              call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.all(aux==41).and.theta==31, &
                'native replay failure retains all caller ledgers')
              do p=1,count
                if(owners(p)/=rank) cycle
                call require(all(local%fields(p)%values==old%fields(p)%values).and. &
                  all(lt%fields(p)%values==ot%fields(p)%values).and.all(ir%fields(p)%values==zero_ir%fields(p)%values), &
                  'native replay failure retains caller gas T IR')
                call require(all(forcing%fields(p)%values==base_forcing%fields(p)%values).and. &
                  all(mass%fields(p)%values==base_mass%fields(p)%values),'failed replay retains accepted forcing outputs')
                if(topology%nodes(p)%level/=1) cycle
                call require_register_equal(registers(p),base_registers(p))
                call require(all(replayed_history%forcing%fields(p)%values==history%forcing%fields(p)%values), &
                  'failed native replay retains prior accepted archive')
              enddo
              reactor_failure=.false.;reactor_failure_rank=0;history_failure=.false.
              ! Replay must rebuild halos/operators from the supplied new gas,
              ! rather than silently reuse its archived raw RHS/state halos.
              local=old;lt=ot;ir=zero_ir
              do p=1,count
                if(owners(p)/=rank.or.topology%nodes(p)%level/=1) cycle
                do z=1,4
                do y=1,4
                do x=1,4
                  if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
                  weight=local%fields(p)%values(1,x,y,z)
                  error=local%fields(p)%values(2,x,y,z)
                  local%fields(p)%values(2,x,y,z)=error+.125_dp*weight
                  local%fields(p)%values(5,x,y,z)=local%fields(p)%values(5,x,y,z)+ &
                    .125_dp*error+.5_dp*weight*.125_dp**2
                enddo
                enddo
                enddo
              enddo
              fresh_start=local;history_calls=0
              call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
                1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
                MPI_COMM_WORLD,good,ir,target_depth=1,ancestor_start=old,ancestor_end=old, &
                temperature_start=ot,temperature_end=ot,fractions_start=[0._dp,0._dp],fractions_end=[1._dp,1._dp], &
                auxiliary_source=history_auxiliary,auxiliary_change=aux,mol_mode=.true.,redistribution=method, &
                target_volume_fraction=.9_dp,state_redist_registers=registers,accepted_forcing=forcing, &
                accepted_mass_source=mass,integrated_flux=flux,interval_history=replayed_history,replay_interval=history)
              call require(good.and.history_calls==0,'changed-state native replay retains frozen provider schedule')
              do p=1,count
                if(owners(p)/=rank.or.topology%nodes(p)%level/=1) cycle
                call require(all(replayed_history%initial%fields(p)%values==fresh_start%fields(p)%values).and. &
                  all(replayed_history%prepared_stages(0)%state_halos%fields(p)%values(:,7:10,7:10,7:10)== &
                    fresh_start%fields(p)%values),'native replay uses supplied fresh gas in its initial halo crop')
                call require(any(replayed_history%prepared_stages(0)%raw_rate_halos%fields(p)%values/= &
                  history%prepared_stages(0)%raw_rate_halos%fields(p)%values).and. &
                  any(local%fields(p)%values/=base%fields(p)%values), &
                  'fresh gas drives an actual new native raw RHS and reacting endpoint')
              enddo
              call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,geometry,extended,old,ot, &
                1.e-10_dp,1,changed_spatial,reference_low_flux,MPI_COMM_WORLD,good,faces, &
                target_depth=1,ancestor_start=old,ancestor_end=old,temperature_start=ot,temperature_end=ot, &
                ancestor_fraction=[0._dp,0._dp],lax_friedrichs=.true.,old_state_halos=fresh_start)
              call require(good,'independent first-order LF replay anchor')
              do blend_case=0,2
                weight=.5_dp*blend_case
                local=old;lt=ot;ir=zero_ir;history_calls=0
                call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
                  1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
                  MPI_COMM_WORLD,good,ir,target_depth=1,ancestor_start=old,ancestor_end=old, &
                  temperature_start=ot,temperature_end=ot,fractions_start=[0._dp,0._dp],fractions_end=[1._dp,1._dp], &
                  auxiliary_source=history_auxiliary,auxiliary_change=aux,mol_mode=.true.,redistribution=method, &
                  target_volume_fraction=.9_dp,state_redist_registers=registers,accepted_forcing=forcing, &
                  accepted_mass_source=mass,integrated_flux=flux,interval_history=replayed_history,replay_interval=history, &
                  replay_spatial_weight=weight)
                call require(good.and.history_calls==0,'weighted native replay preserves frozen provider schedule')
                call require(all(replayed_history%wall_change==w),'weighted replay retains actual native interval wall ledger')
                do stage=0,2
                  associate(saved=>replayed_history%prepared_stages(stage))
                  call require(saved%spatial_weight==weight,'each native stage records the same replay weight')
                  if(weight==1) then
                    call require(.not.allocated(saved%low_spatial_rate%fields), &
                      'unit replay weight keeps native stage without constructing a low family')
                    cycle
                  endif
                  changed_spatial=old
                  do p=1,count
                    if(owners(p)==rank) changed_spatial%fields(p)%values=0
                  enddo
                  if(allocated(replay_registers)) deallocate(replay_registers)
                  call compute_state_redist_mpi_level_stage_3d(topology,owners,extended,saved%state_halos, &
                    saved%raw_rate_halos,nv,1,saved%interval,changed_spatial,replay_registers,MPI_COMM_WORLD,good, &
                    faces,saved%target_volume_fraction,saved%max_order,register_scale=saved%register_scale)
                  call require(good,'weighted stage independently recomputes native high family')
                  do p=1,count
                    if(owners(p)/=rank.or.topology%nodes(p)%level/=1) cycle
                    call post_eb_spatial_density_rate_3d(geometry(p),changed_spatial%fields(p)%values,good)
                    call require(good,'weighted stage native high density consistency')
                    call require(all(saved%low_faces(p)%x==reference_low_flux(p)%x).and. &
                      all(saved%low_faces(p)%y==reference_low_flux(p)%y).and. &
                      all(saved%low_faces(p)%z==reference_low_flux(p)%z),'all stages share actual first-order old LF faces')
                    call require(all(saved%faces(p)%x==(1-weight)*saved%low_faces(p)%x+weight*saved%native_faces(p)%x).and. &
                      all(saved%faces(p)%y==(1-weight)*saved%low_faces(p)%y+weight*saved%native_faces(p)%y).and. &
                      all(saved%faces(p)%z==(1-weight)*saved%low_faces(p)%z+weight*saved%native_faces(p)%z), &
                      'all physical face components use the same native/low stage weight')
                    call require(all(saved%spatial_rate%fields(p)%values(2:,:,:,:)== &
                      (1-weight)*saved%low_spatial_rate%fields(p)%values(2:,:,:,:)+ &
                      weight*changed_spatial%fields(p)%values(2:,:,:,:)), &
                      'momentum energy and every species use the same stage weight')
                    call require(all(saved%registers(p)%coarse_delta==weight*replay_registers(p)%coarse_delta).and. &
                      all(saved%registers(p)%fine_delta==weight*replay_registers(p)%fine_delta), &
                      'coarse and fine deltas co-blend with the same conservative stage family')
                    do z=1,4
                    do y=1,4
                    do x=1,4
                      if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
                      call require(saved%spatial_rate%fields(p)%values(1,x,y,z)== &
                        sum(saved%spatial_rate%fields(p)%values(6:,x,y,z)), &
                        'weighted native stage preserves ordered density/species consistency')
                    enddo
                    enddo
                    enddo
                  enddo
                  end associate
                enddo
                call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,replayed_history, &
                  replayed_history%prepared_stages,nv,reconstructed,replay_registers,MPI_COMM_WORLD,good,faces)
                call require(good,'weighted stage archive reconstructs native full interval')
                do p=1,count
                  if(owners(p)/=rank) cycle
                  call require(all(reconstructed%fields(p)%values==forcing%fields(p)%values), &
                    'weighted native archive reproduces full accepted forcing exactly')
                  if(topology%nodes(p)%level/=1) cycle
                  call require_register_equal(replay_registers(p),registers(p))
                  if(weight==1) call require(all(local%fields(p)%values==base%fields(p)%values).and. &
                    all(lt%fields(p)%values==base_t%fields(p)%values).and.all(ir%fields(p)%values==base_ir%fields(p)%values), &
                    'explicit unit replay weight retains original native gas T IR exactly')
                enddo
              enddo
              kept=replayed_history;saved_registers=replay_registers
              if(rank==nranks-1) kept%prepared_stages%spatial_weight=0
              call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,kept, &
                kept%prepared_stages,nv,reconstructed,replay_registers,MPI_COMM_WORLD,good,faces)
              call require(.not.good,'rank-divergent archived weights reject before mixed-family collectives')
              do p=1,count
                if(owners(p)/=rank) cycle
                call require(all(reconstructed%fields(p)%values==forcing%fields(p)%values), &
                  'rejected archived weight preserves output forcing')
                if(topology%nodes(p)%level==1) call require_register_equal(replay_registers(p),saved_registers(p))
              enddo
              local=old;lt=ot;ir=zero_ir;history_calls=0
              weight=0
              if(rank==nranks-1) weight=.5_dp
              b=17;w=19;ch=23;ext=29;aux=41;theta=31
              call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
                1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
                MPI_COMM_WORLD,good,ir,target_depth=1,ancestor_start=old,ancestor_end=old, &
                temperature_start=ot,temperature_end=ot,fractions_start=[0._dp,0._dp],fractions_end=[1._dp,1._dp], &
                auxiliary_source=history_auxiliary,auxiliary_change=aux,mol_mode=.true.,redistribution=method, &
                target_volume_fraction=.9_dp,state_redist_registers=registers,accepted_forcing=forcing, &
                accepted_mass_source=mass,integrated_flux=flux,interval_history=replayed_history,replay_interval=history, &
                replay_spatial_weight=weight)
              call require(.not.good.and.history_calls==0,'rank-divergent spatial weights reject collectively')
              call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.all(aux==41).and.theta==31, &
                'rejected replay weight preserves caller ledgers')
              do p=1,count
                if(owners(p)/=rank) cycle
                call require(all(local%fields(p)%values==old%fields(p)%values).and. &
                  all(lt%fields(p)%values==ot%fields(p)%values).and.all(ir%fields(p)%values==zero_ir%fields(p)%values), &
                  'rejected spatial weight preserves caller state')
              enddo
              local=base;lt=base_t;ir=base_ir;forcing=base_forcing;mass=base_mass;registers=base_registers
            else
              call require(.not.allocated(history%prepared_stages),'unprepared operator has no StateRedist stage archive')
            endif
            call reconstruct_mpi_eb_mol_interval_3d(topology,owners,geometry,history,history%spatial_flux,nv, &
              reconstructed,MPI_COMM_WORLD,good)
            call require(good,'stage rerun leaves caller forcing reconstruction available')
            do p=1,count
              call require(allocated(history%volume_source%fields(p)%values).eqv. &
                (owners(p)==rank.and.topology%nodes(p)%level==1),'redistribution rates retain sparse owner/level layout')
              if(owners(p)/=rank) cycle
              call require(all(local%fields(p)%values==base%fields(p)%values).and. &
                all(lt%fields(p)%values==base_t%fields(p)%values).and.all(ir%fields(p)%values==base_ir%fields(p)%values), &
                'capturing redistribution history leaves gas T IR exact')
              call require(all(reconstructed%fields(p)%values==forcing%fields(p)%values), &
                'complete accepted forcing remains authoritative without residual subtraction')
              call require(all(forcing%fields(p)%values==base_forcing%fields(p)%values).and. &
                all(mass%fields(p)%values==base_mass%fields(p)%values),'archive leaves exported forcing families exact')
              if(method=='flux_redist') then
                call require(all(transfer%fields(p)%values==base_transfer%fields(p)%values), &
                  'archive leaves exported redistribution transfer exact')
                call require(allocated(history%redistribution_transfer%fields(p)%values), &
                  'redistribution impulse retains every owned target level')
                call require(all(history%redistribution_transfer%fields(p)%values==transfer%fields(p)%values), &
                  'archive preserves the exact accepted off-level redistribution impulse')
              endif
              if(topology%nodes(p)%level/=1) cycle
              call require(all(history%volume_source%fields(p)%values(1,:,:,:)==mass%fields(p)%values(1,:,:,:)).and. &
                all(history%volume_source%fields(p)%values(6:,:,:,:)==mass%fields(p)%values(6:,:,:,:)), &
                'volume family retains accepted density/species forcing')
              call require(any(history%volume_source%fields(p)%values(2:4,:,:,:)/=0), &
                'volume family also retains momentum absent from mass-only history')
              do z=1,4
              do y=1,4
              do x=1,4
                if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
                expected=source_profile(p,x,y,z,0)
                call require(history%old_volume_source%fields(p)%values(1,x,y,z)==expected(1).and. &
                  all(history%old_volume_source%fields(p)%values(6:,x,y,z)==expected(6:)), &
                  'old volume family is unweighted actual old source')
                error=maxval(abs(history%forcing%fields(p)%values(:,x,y,z)- &
                  history%volume_source%fields(p)%values(:,x,y,z)-history%spatial_rate%fields(p)%values(:,x,y,z))/ &
                  max(1._dp,abs(history%forcing%fields(p)%values(:,x,y,z))))
                call require(error<3.e-14_dp,'separate volume/spatial families compose the complete accepted forcing')
              enddo
              enddo
              enddo
              if(method=='state_redist') then
                call require_register_equal(registers(p),base_registers(p))
                call require_register_equal(history%redist_registers(p),registers(p))
                call require(history%old_redist_registers(p)%is_valid(nv,[16,16,16]), &
                  'old prepared register retains original halo masks and scale')
              endif
            enddo
          endif
        enddo
        kept=history
        if(allocated(registers)) saved_registers=registers
        local=old;lt=ot;ir=zero_ir;history_calls=0;history_failure=.true.
        b=17;w=19;ch=23;ext=29;aux=41;theta=31
        if(allocated(optional_history)) deallocate(optional_history)
        allocate(optional_history,source=history)
      if(method=='state_redist') then
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
          1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
          MPI_COMM_WORLD,good,ir,target_depth=1,ancestor_start=old,ancestor_end=old,temperature_start=ot,temperature_end=ot, &
          fractions_start=[0._dp,0._dp],fractions_end=[1._dp,1._dp],auxiliary_source=history_auxiliary,auxiliary_change=aux, &
          mol_mode=.true.,redistribution=method,target_volume_fraction=.9_dp,state_redist_registers=registers, &
          accepted_forcing=forcing,accepted_mass_source=mass,integrated_flux=flux,interval_history=optional_history)
      else
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
          1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
          MPI_COMM_WORLD,good,ir,target_depth=1,ancestor_start=old,ancestor_end=old,temperature_start=ot,temperature_end=ot, &
          fractions_start=[0._dp,0._dp],fractions_end=[1._dp,1._dp],auxiliary_source=history_auxiliary,auxiliary_change=aux, &
          mol_mode=.true.,redistribution=method,target_volume_fraction=.9_dp,redistribution_transfer=transfer, &
          accepted_forcing=forcing,accepted_mass_source=mass,integrated_flux=flux,interval_history=optional_history)
      endif

        history=optional_history
        call require(.not.good.and.history_calls==2,'idle-rank source failure rejects redistribution archive transaction')
        call require(all(history%wall_change==kept%wall_change),'failed capture retains prior interval wall ledger')
        call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.all(aux==41).and.theta==31, &
          'failed redistributed capture retains inventory outputs')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==old%fields(p)%values).and.all(lt%fields(p)%values==ot%fields(p)%values).and. &
            all(ir%fields(p)%values==zero_ir%fields(p)%values),'failed redistributed capture retains gas T IR')
          if(method=='flux_redist') call require(all(transfer%fields(p)%values== &
            kept%redistribution_transfer%fields(p)%values),'failed capture retains caller transfer impulse')
          if(method=='flux_redist') call require(all(history%redistribution_transfer%fields(p)%values== &
            kept%redistribution_transfer%fields(p)%values),'failed capture retains prior transfer impulse')
          if(topology%nodes(p)%level/=1) cycle
          call require(all(history%spatial_rate%fields(p)%values==kept%spatial_rate%fields(p)%values).and. &
            all(history%old_spatial_rate%fields(p)%values==kept%old_spatial_rate%fields(p)%values).and. &
            all(history%volume_source%fields(p)%values==kept%volume_source%fields(p)%values).and. &
            all(history%old_volume_source%fields(p)%values==kept%old_volume_source%fields(p)%values), &
            'failed capture retains complete prior source families')
          if(method=='state_redist') then
            call require_register_equal(history%redist_registers(p),kept%redist_registers(p))
            call require_register_equal(history%old_redist_registers(p),kept%old_redist_registers(p))
            call require_register_equal(registers(p),saved_registers(p))
            do stage=0,2
              call require(all(history%prepared_stages(stage)%state_halos%fields(p)%values== &
                kept%prepared_stages(stage)%state_halos%fields(p)%values).and. &
                all(history%prepared_stages(stage)%raw_rate_halos%fields(p)%values== &
                kept%prepared_stages(stage)%raw_rate_halos%fields(p)%values).and. &
                all(history%prepared_stages(stage)%spatial_rate%fields(p)%values== &
                kept%prepared_stages(stage)%spatial_rate%fields(p)%values), &
                'late source failure retains prior complete prepared stage inputs')
              call require_register_equal(history%prepared_stages(stage)%registers(p),kept%prepared_stages(stage)%registers(p))
            enddo
          endif
        enddo
        history_failure=.false.
        if(allocated(registers)) deallocate(registers)
        if(allocated(saved_registers)) deallocate(saved_registers)
        if(allocated(base_registers)) deallocate(base_registers)
      enddo
      deallocate(old%fields,ot%fields,zero_ir%fields)
    enddo
  end subroutine check_redist_interval_history

  subroutine require_register_equal(a,b)
    type(state_redist_register_3d), intent(in) :: a,b
    call require(a%scale==b%scale,'redistribution archive retains exact register scale')
    call require((allocated(a%coarse_flags).eqv.allocated(b%coarse_flags)).and. &
      (allocated(a%fine_ghost).eqv.allocated(b%fine_ghost)).and.(allocated(a%valid_cells).eqv.allocated(b%valid_cells)).and. &
      (allocated(a%coarse_delta).eqv.allocated(b%coarse_delta)).and.(allocated(a%fine_delta).eqv.allocated(b%fine_delta)), &
      'redistribution archive retains register component presence')
    if(allocated(a%coarse_flags)) call require(all(a%coarse_flags==b%coarse_flags),'coarse register masks exact')
    if(allocated(a%fine_ghost)) call require(all(a%fine_ghost.eqv.b%fine_ghost),'fine register masks exact')
    if(allocated(a%valid_cells)) call require(all(a%valid_cells.eqv.b%valid_cells),'valid register masks exact')
    if(allocated(a%coarse_delta)) call require(all(a%coarse_delta==b%coarse_delta),'coarse register deltas exact')
    if(allocated(a%fine_delta)) call require(all(a%fine_delta==b%fine_delta),'fine register deltas exact')
  end subroutine require_register_equal

  subroutine check_interval_history(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(eb_geometry_3d) :: mol_geometry(count)
    type(amr_tree_state_3d) :: old,ot,local,lt,ir,forcing,mass_source,kept_forcing,kept_mass
    type(amr_tree_state_3d) :: replay,replay_t,replay_ir,spatial
    type(eb_integrated_flux_3d), allocatable :: old_flux(:),expected_flux(:),kept_flux(:),flux(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(nv),theta,expected(nv)
    real(dp), parameter :: dt=1.e-10_dp
    integer :: owners(count),layout,p,x,y,z,iterations
    logical :: good
    do p=1,count
      call extract_eb_geometry_3d(extended(p),[4,4,4],[13,13,13],mol_geometry(p),good)
      call require(good,'interval history MOL geometry')
    enddo
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(old%fields(count),ot%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        old%fields(p)%values=initial%nodes(p)%state
        allocate(ot%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4))
        ot%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature;ir%fields(p)%values=0
      enddo
      do iterations=1,2
      do p=1,count
        if(owners(p)/=rank) cycle
        ir%fields(p)%values=0
      enddo
      local=old;lt=ot;history_calls=0;history_failure=.false.
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
        dt,iterations,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,target_depth=2,ancestor_start=old,ancestor_end=old, &
        temperature_start=ot,temperature_end=ot,fractions_start=[.25_dp,.25_dp],fractions_end=[.75_dp,.75_dp], &
        auxiliary_source=history_auxiliary,auxiliary_change=aux,mol_mode=.true., &
        accepted_forcing=forcing,accepted_mass_source=mass_source,old_spatial_flux=old_flux,integrated_flux=flux)
      call require(good.and.history_calls==2,'accepted history uses old and first-new MOL sources once')
      call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,geometry,mol_geometry,old,ot,dt,2, &
        spatial,expected_flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=2, &
        ancestor_start=old,ancestor_end=old,temperature_start=ot,temperature_end=ot,ancestor_fraction=[.25_dp,.25_dp])
      call require(good,'independent unweighted old MOL flux')
      replay=old;replay_t=ot;replay_ir=ir
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,forcing, &
        dt,2.e-7_dp,1.e-12_dp,replay,replay_t,replay_ir,MPI_COMM_WORLD,good,target_depth=2)
      call require(good,'accepted forcing reconstructs the reacted interval')
      do p=1,count
        call require(allocated(forcing%fields(p)%values).eqv.(owners(p)==rank),'history remains owner-local')
        if(owners(p)/=rank) cycle
        call require(all(mass_source%fields(p)%values(2:5,:,:,:)==0),'mass history excludes momentum and energy')
        if(topology%nodes(p)%level/=2) then
          call require(.not.allocated(old_flux(p)%x),'inactive old flux remains absent')
          call require(all(mass_source%fields(p)%values==0),'inactive mass history is zero')
          cycle
        endif
        call require(all(replay%fields(p)%values==local%fields(p)%values).and. &
          all(replay_t%fields(p)%values==lt%fields(p)%values).and.all(replay_ir%fields(p)%values==ir%fields(p)%values), &
          'exported nonchemical forcing reproduces state temperature and IR exactly')
        call require(all(old_flux(p)%x==expected_flux(p)%x).and.all(old_flux(p)%y==expected_flux(p)%y).and. &
          all(old_flux(p)%z==expected_flux(p)%z),'history retains physical f0 without half-stage weighting')
        do z=1,4
        do y=1,4
        do x=1,4
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          ! The first reacting MOL pass includes both source endpoints.
          ! Later Picard passes retain only the old nonspatial family.
          expected=merge(1.5_dp,.5_dp,iterations==1)*source_profile(p,x,y,z,0);expected(2:5)=0
          call require(all(mass_source%fields(p)%values(:,x,y,z)==expected), &
            'accepted mass source follows pinned first/later reacting MOL weights')
        enddo
        enddo
        enddo
      enddo
      kept_forcing=forcing;kept_mass=mass_source;kept_flux=old_flux
      history_calls=0;history_failure=.true.
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
        dt,iterations,.true.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,target_depth=2,ancestor_start=old,ancestor_end=old, &
        temperature_start=ot,temperature_end=ot,fractions_start=[.25_dp,.25_dp],fractions_end=[.75_dp,.75_dp], &
        auxiliary_source=history_auxiliary,auxiliary_change=aux,mol_mode=.true., &
        accepted_forcing=forcing,accepted_mass_source=mass_source,old_spatial_flux=old_flux,integrated_flux=flux)
      call require(.not.good.and.history_calls==2,'late source failure rejects interval history collectively')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(forcing%fields(p)%values==kept_forcing%fields(p)%values).and. &
          all(mass_source%fields(p)%values==kept_mass%fields(p)%values),'failed step retains both history rates')
        if(topology%nodes(p)%level/=2) cycle
        call require(all(old_flux(p)%x==kept_flux(p)%x).and.all(old_flux(p)%y==kept_flux(p)%y).and. &
          all(old_flux(p)%z==kept_flux(p)%z),'failed step retains physical old flux')
      enddo
      history_failure=.false.
      enddo
      deallocate(old%fields,ot%fields,ir%fields)
    enddo
  end subroutine check_interval_history
  subroutine history_auxiliary(state,temperature,interval,stage,depth,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p,x,y,z,local_rank,ierr
    call MPI_Comm_rank(comm,local_rank,ierr)
    call require(ierr==MPI_SUCCESS.and.local_rank==rank,'history callback communicator')
    call require(interval==1.e-10_dp.and.(depth==2.or.(test_mode=='redist_history'.and.depth==1)).and. &
      stage==history_calls,'history callback interval/staging')
    call require(allocated(temperature%fields).and.size(temperature%fields)==count,'history temperature ownership')
    history_calls=history_calls+1
    do p=1,count
      if(.not.allocated(state%fields(p)%values)) cycle
      do z=1,4
      do y=1,4
      do x=1,4
        rate%fields(p)%values(:,x,y,z)=source_profile(p,x,y,z,stage)
      enddo
      enddo
      enddo
    enddo
    good=.not.(history_failure.and.stage==1.and.rank==nranks-1)
  end subroutine history_auxiliary
  subroutine check_auxiliary_sources(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(eb_reactive_tree_3d) :: reference,profile
    type(amr_tree_state_3d) :: local,lt,ir,serial_ir,old_ir,old,ot
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),sb(nv),sw(nv),sch(nv),sext(nv),limiter,slimiter
    real(dp) :: aux(nv),saux(nv),expected_aux(nv),initial_inventory(nv),final_inventory(nv),error
    integer :: owners(count),p,x,y,z,layout,failure
    logical :: good
    serial_source_calls=0;source_failure=0
    allocate(serial_ir%fields(count))
    do p=1,count
      allocate(serial_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
    end do
    reference=initial
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,test_mode/='seed_aux', &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,slimiter,good,serial_ir,auxiliary_source=serial_auxiliary,auxiliary_change=saux, &
      state_redist_max_order=requested_order)
    call require(good,'serial coupled SDC')
    call require_restricted_reactions(reference,serial_ir)
    call require(serial_source_calls==3,'one old and two new serial source calls')
    profile=initial
    do p=1,count
      do z=1,4
      do y=1,4
      do x=1,4
        profile%nodes(p)%state(:,x,y,z)=source_profile(p,x,y,z,0)
        if(profile%nodes(p)%geometry%volume_fraction(x,y,z)<=0) profile%nodes(p)%state(:,x,y,z)=0
      end do
      end do
      end do
    end do
    call integrate_reactive_eb_tree_3d(profile,expected_aux,good)
    call require(good,'source profile inventory')
    expected_aux=2.e-10_dp*expected_aux
    call require(maxval(abs(saux-expected_aux)/max(1.e-30_dp,abs(expected_aux)))<3.e-14_dp, &
      'accepted old/new source inventory excludes intermediate iterations')
    call integrate_reactive_eb_tree_3d(initial,initial_inventory,good)
    call require(good,'SDC initial inventory')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(local%fields(count),lt%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=initial%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4))
        lt%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
        ir%fields(p)%values=0
      end do
      predictor_owners=owners;predictor_geometry=geometry;predictor_metrics=lt
      do p=1,count
        if(owners(p)/=rank) cycle
        deallocate(predictor_metrics%fields(p)%values)
        allocate(predictor_metrics%fields(p)%values(1,10,10,10))
        predictor_metrics%fields(p)%values(1,:,:,:)=extended(p)%volume_fraction(4:13,4:13,4:13)
      end do
      predictor_calls=0
      old=local;ot=lt;old_ir=ir;mpi_source_calls=0;source_failure=0
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
        1.e-10_dp,2,test_mode/='seed_aux',[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,auxiliary_source=mpi_auxiliary,auxiliary_change=aux, &
          predictor_builder=build_predictor_halo, &
      state_redist_max_order=requested_order)
      call require(good,'MPI coupled SDC')
      call require(mpi_source_calls==3,'one old and two new MPI source calls')
      call require(predictor_calls==2,'one predictor build per SDC iteration')
      call require(maxval(abs(aux-saux)/max(1.e-30_dp,abs(saux)))<3.e-14_dp,'separate serial MPI source ledger')
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'coupled SDC owner-local state')
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==reference%nodes(p)%state),'exact serial coupled SDC state')
        call require(all(lt%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'exact serial coupled SDC T')
        call require(all(ir%fields(p)%values==serial_ir%fields(p)%values),'exact serial coupled SDC IR')
      end do
      call require(limiter==slimiter,'exact coupled SDC limiter')
      error=max(maxval(abs(b-sb)/max(1._dp,abs(sb))),maxval(abs(w-sw)/max(1._dp,abs(sw))), &
        maxval(abs(ch-sch)/max(1._dp,abs(sch))),maxval(abs(ext-sext)/max(1._dp,abs(sext))))
      call require(error<2.e-12_dp,'serial coupled SDC inventories')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,final_inventory,MPI_COMM_WORLD,good)
      call require(good,'SDC final inventory')
      error=maxval(abs(final_inventory-initial_inventory-b-w-ch-ext-aux)/max(1._dp,abs(initial_inventory)))
      if(rank==0) print *, 'MPI SDC layout, composite residual: ',layout,error
      call require(error<3.e-14_dp,'SDC composite closure')
      do failure=1,4
        local=old;lt=ot;ir=old_ir;b=17;w=19;ch=23;ext=29;limiter=31
        source_failure=failure;mpi_source_calls=0;predictor_calls=0;aux=37
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
          1.e-10_dp,2,test_mode/='seed_aux',[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
          b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,auxiliary_source=mpi_auxiliary,auxiliary_change=aux, &
          predictor_builder=build_predictor_halo, &
      state_redist_max_order=requested_order)
        call require(.not.good,'collective late auxiliary source rejection')
        call require(mpi_source_calls==merge(2,3,failure==4),'source calls before late failure')
        call require(predictor_calls==2,'late predictor build reached')
        call require(all(aux==37),'auxiliary ledger rollback')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==old%fields(p)%values),'coupled SDC state rollback')
          call require(all(lt%fields(p)%values==ot%fields(p)%values),'coupled SDC T rollback')
          call require(all(ir%fields(p)%values==old_ir%fields(p)%values),'coupled SDC IR rollback')
        end do
        call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.limiter==31, &
          'coupled SDC inventory rollback')
      end do
      deallocate(local%fields,lt%fields,ir%fields)
    end do
  end subroutine check_auxiliary_sources
  subroutine build_predictor_halo(rate,interval,iteration,depth,halo,comm,good)
    type(amr_tree_state_3d), intent(in) :: rate
    real(dp), intent(in) :: interval
    integer, intent(in) :: iteration,depth
    type(amr_tree_state_3d), intent(inout) :: halo
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    call require(interval==1.e-10_dp.and.depth==-1,'predictor interval/depth')
    predictor_calls=predictor_calls+1
    call require(iteration==predictor_calls,'predictor iteration sequence')
    call fill_sdc_mpi_eb_tree_source_halo_3d(topology,predictor_owners,predictor_geometry,rate,nv,3,halo, &
      comm,good,target_volume_fraction=predictor_metrics)
    if(iteration==2.and.source_failure==4.and.rank==nranks-1) good=.false.
  end subroutine
  subroutine serial_auxiliary(state,interval,stage,depth,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p,x,y,z,d(3)
    call require(interval==1.e-10_dp.and.depth==-1,'serial source interval/depth')
    call require(stage==serial_source_calls,'serial source stage sequence')
    serial_source_calls=serial_source_calls+1
    do p=1,size(state%nodes)
      d=state%topology%nodes(p)%shape
      do z=1,d(3)
      do y=1,d(2)
      do x=1,d(1)
        rate%fields(p)%values(:,x,y,z)=source_profile(p,x,y,z,stage)
      end do
      end do
      end do
    end do
    good=.true.
  end subroutine
  subroutine mpi_auxiliary(state,temperature,interval,stage,depth,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p,x,y,z,d(3),local_rank,ierr
    call MPI_Comm_rank(comm,local_rank,ierr)
    call require(ierr==MPI_SUCCESS.and.local_rank==rank,'source communicator')
    call require(interval==1.e-10_dp.and.depth==-1,'MPI source interval/depth')
    call require(stage==mpi_source_calls,'MPI source stage sequence')
    mpi_source_calls=mpi_source_calls+1;good=.true.
    do p=1,size(state%fields)
      if(.not.allocated(state%fields(p)%values)) cycle
      call require(allocated(temperature%fields(p)%values),'source temperature ownership')
      d=shape(temperature%fields(p)%values(1,:,:,:))
      do z=1,d(3)
      do y=1,d(2)
      do x=1,d(1)
        rate%fields(p)%values(:,x,y,z)=source_profile(p,x,y,z,stage)
      end do
      end do
      end do
    end do
    if(stage/=2) return
    if(source_failure==1.and.rank==nranks-1) good=.false.
    if(allocated(state%fields(1)%values)) then
      if(source_failure==2) rate%fields(1)%values(:,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
      if(source_failure==3) deallocate(rate%fields(1)%values)
    end if
  end subroutine
  function source_profile(p,x,y,z,stage) result(value)
    integer, intent(in) :: p,x,y,z,stage
    real(dp) :: value(nv),amount
    amount=1000._dp*real(stage+1,dp)*(p+.01_dp*x+.02_dp*y+.03_dp*z)
    value=0;value(1)=amount;value(2:4)=amount*[2._dp,-3._dp,4._dp]
    value(5)=amount*2.e6_dp;value(6)=amount
    if(tree%nodes(p)%geometry%volume_fraction(x,y,z)<=0) value=ieee_value(0._dp,ieee_quiet_nan)
  end function
  subroutine check_redist_level(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(eb_reactive_tree_3d) :: reference,ledger
    type(amr_tree_state_3d) :: local,lt,ir,serial_ir,transfer,serial_transfer,old,ot,old_ir,old_transfer
    type(eb_integrated_flux_3d), allocatable :: flux(:),serial_flux(:)
    integer :: owners(count),p,layout,x,y,z,c
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),sb(nv),sw(nv),sch(nv),sext(nv),limiter,slimiter
    real(dp) :: flux_inventory(nv),change(nv),error
    logical :: good
    reference=initial
    allocate(serial_ir%fields(count))
    do p=1,count
      allocate(serial_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
    end do
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,slimiter,good,serial_ir,target_depth=1,integrated_flux=serial_flux, &
      redistribution='flux_redist',redistribution_transfer=serial_transfer, &
      state_redist_max_order=requested_order)
    call require(good,'serial redistributed selected SDC')
    ledger=initial
    do p=1,count
      ledger%nodes(p)%state=0
      if(topology%nodes(p)%level/=1) cycle
      associate(g=>geometry(p),f=>serial_flux(p))
      do z=1,4
      do y=1,4
      do x=1,4
        if(g%volume_fraction(x,y,z)<=0) cycle
        do c=1,nv
          ledger%nodes(p)%state(c,x,y,z)=( &
            g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(c,x-1,y,z)-g%x_face_fraction(x,y,z)*f%x(c,x,y,z))+ &
            g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(c,x,y-1,z)-g%y_face_fraction(x,y,z)*f%y(c,x,y,z))+ &
            g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(c,x,y,z-1)-g%z_face_fraction(x,y,z)*f%z(c,x,y,z)))/ &
            (g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz)
        end do
      end do
      end do
      end do
      end associate
    end do
    call integrate_reactive_eb_tree_3d(ledger,flux_inventory,good)
    call require(good,'selected exported flux divergence')
    do p=1,count
      ledger%nodes(p)%state=reference%nodes(p)%state-initial%nodes(p)%state+serial_transfer%fields(p)%values
      if(topology%nodes(p)%level==1) then
        call require(all(serial_transfer%fields(p)%values==0),'selected level has no deferred transfer')
      else
        call require(all(reference%nodes(p)%state==initial%nodes(p)%state),'other gas waits for transfer')
      end if
    end do
    call integrate_reactive_eb_tree_3d(ledger,change,good)
    call require(good,'selected state plus transfer inventory')
    error=maxval(abs(change-flux_inventory-sw-sch-sext)/max(1._dp,abs(flux_inventory)))
    print *, 'Selected redistribution flux/ledger residual: ',error
    call require(error<2.e-9_dp,'selected state plus transfer closes exported fluxes')
    call require(any(serial_transfer%fields(1)%values/=0),'coarse transfer is active')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(local%fields(count),lt%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=initial%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4))
        lt%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
        ir%fields(p)%values=0
      end do
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
        1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,target_depth=1,integrated_flux=flux, &
        redistribution='flux_redist',redistribution_transfer=transfer, &
      state_redist_max_order=requested_order)
      call require(good,'MPI redistributed selected SDC')
      do p=1,count
        call require(allocated(transfer%fields(p)%values).eqv.(owners(p)==rank),'owner-local stage transfer')
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==reference%nodes(p)%state),'exact selected redistributed gas')
        call require(all(lt%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'exact selected redistributed T')
        call require(all(ir%fields(p)%values==serial_ir%fields(p)%values),'exact selected redistributed IR')
        call require(all(transfer%fields(p)%values==serial_transfer%fields(p)%values),'exact stage transfer')
      end do
      error=max(maxval(abs(b-sb)/max(1._dp,abs(sb))),maxval(abs(w-sw)/max(1._dp,abs(sw))), &
        maxval(abs(ch-sch)/max(1._dp,abs(sch))),maxval(abs(ext-sext)/max(1._dp,abs(sext))))
      call require(error<2.e-12_dp.and.limiter==slimiter,'selected redistributed diagnostics')
      old=local;ot=lt;old_ir=ir;old_transfer=transfer;b=17;w=19;ch=23;ext=29;limiter=31
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
        1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,difmag=-1._dp,target_depth=1,integrated_flux=flux, &
        redistribution='flux_redist',redistribution_transfer=transfer, &
      state_redist_max_order=requested_order)
      call require(.not.good,'selected redistribution failed stage')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==old%fields(p)%values),'selected redistributed gas rollback')
        call require(all(lt%fields(p)%values==ot%fields(p)%values),'selected redistributed T rollback')
        call require(all(ir%fields(p)%values==old_ir%fields(p)%values),'selected redistributed IR rollback')
        call require(all(transfer%fields(p)%values==old_transfer%fields(p)%values),'stage transfer rollback')
      end do
      call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.limiter==31, &
        'selected redistributed diagnostics rollback')
      deallocate(local%fields,lt%fields,ir%fields)
    end do
  end subroutine check_redist_level
  subroutine recursive_source_schedule(interval,stage,depth,time,calls)
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    integer, intent(inout) :: calls
    integer, parameter :: depths(7)=[0,1,2,2,1,2,2],ticks(7)=[0,0,0,1,2,2,3]
    integer :: n
    n=calls/3+1
    call require(n<=7,'recursive source interval count')
    call require(stage==mod(calls,3).and.depth==depths(n),'recursive source stage/depth order')
    call require(abs(time-2.5e-11_dp*ticks(n))<1e-25_dp,'recursive source start time')
    call require(abs(interval-1e-10_dp/2**depth)<1e-25_dp,'recursive source interval')
    calls=calls+1
  end subroutine
  subroutine serial_recursive_aux(state,interval,stage,depth,time,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p,x,y,z
    call recursive_source_schedule(interval,stage,depth,time,serial_source_calls)
    do p=1,size(state%nodes)
      do z=1,4
      do y=1,4
      do x=1,4
        rate%fields(p)%values(:,x,y,z)=(1+time/1e-10_dp)*source_profile(p,x,y,z,stage)
      end do
      end do
      end do
    end do
    good=.true.
  end subroutine
  subroutine mpi_recursive_aux(state,temperature,interval,stage,depth,time,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p,x,y,z
    call recursive_source_schedule(interval,stage,depth,time,mpi_source_calls)
    do p=1,size(state%fields)
      if(.not.allocated(state%fields(p)%values)) cycle
      do z=1,4
      do y=1,4
      do x=1,4
        rate%fields(p)%values(:,x,y,z)=(1+time/1e-10_dp)*source_profile(p,x,y,z,stage)
      end do
      end do
      end do
    end do
    good=.not.(source_failure==1.and.mpi_source_calls==21.and.rank==nranks-1)
  end subroutine
  subroutine check_mol_parity(initial,geometry,extended,faces,force,frame)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(rotating_frame_3d) :: bad_frame
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(eb_reactive_tree_3d) :: reference
    type(amr_tree_state_3d) :: local,lt,snapshot,snapshot_t
    type(eb_integrated_flux_3d), allocatable :: serial_registered(:),mpi_registered(:), &
      serial_physical(:),mpi_physical(:)
    integer :: owners(count),p,variant
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),sb(nv),sw(nv),sch(nv),sext(nv),limiter,slimiter,error
    real(dp) :: initial_inventory(nv),final_inventory(nv)
    logical :: good,reacting,lf,flags(4)
    character(len=16) :: method
    call integrate_reactive_eb_tree_3d(initial,initial_inventory,good)
    call require(good,'MOL initial inventory')
    do variant=1,merge(6,4,present(frame))
      reacting=mod(variant,2)==0;lf=variant==2.or.variant==3
      flags=variant>1;method='none'
      if(variant==3) method='flux_redist'
      if(variant==4) method='state_redist'
      if(variant==6) then
        reacting=.false.
        flags=.false.
        method='state_redist'
      end if
      reference=initial
      if(variant<=2) then
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,reacting,flags,force, &
          2.e-7_dp,1.e-12_dp,reference,extended,faces,sb,sw,sch,sext,slimiter,good, &
          redistribution=method,mol_mode=.true.,mol_lax_friedrichs=lf,frame=frame, &
      state_redist_max_order=requested_order)
      else
        if(variant==5) then
          call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,2,reacting,flags,force, &
            2.e-7_dp,1.e-12_dp,reference,extended,faces,sb,sw,sch,sext,slimiter,good, &
            redistribution=method,mol_mode=.true.,mol_lax_friedrichs=lf,frame=frame, &
            embedded_wall_temperature=1100._dp,embedded_noslip=.true., &
      state_redist_max_order=requested_order)
        else
          call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,2,reacting,flags,force, &
            2.e-7_dp,1.e-12_dp,reference,extended,faces,sb,sw,sch,sext,slimiter,good, &
            redistribution=method,mol_mode=.true.,mol_lax_friedrichs=lf,frame=frame, &
      state_redist_max_order=requested_order,registered_flux=serial_registered,integrated_flux=serial_physical, &
      target_volume_fraction=merge(.8_dp,.5_dp,variant==6))
        end if
      end if
      call require(good,'serial MOL evolution')
      owners=[(mod(p-1,nranks),p=1,count)]
      allocate(local%fields(count),lt%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=initial%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4))
        lt%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
      end do
      if(variant<=2) then
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
          1.e-10_dp,2,reacting,flags,force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,limiter,MPI_COMM_WORLD,good, &
          redistribution=method,mol_mode=.true.,mol_lax_friedrichs=lf,frame=frame, &
      state_redist_max_order=requested_order)
      else
        if(variant==5) then
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
            1.e-10_dp,2,reacting,flags,force,2.e-7_dp,1.e-12_dp,extended,faces,b,w,ch,ext,limiter,MPI_COMM_WORLD,good, &
            redistribution=method,mol_mode=.true.,mol_lax_friedrichs=lf,frame=frame, &
            embedded_wall_temperature=1100._dp,embedded_noslip=.true., &
      state_redist_max_order=requested_order)
        else
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
            1.e-10_dp,2,reacting,flags,force,2.e-7_dp,1.e-12_dp,extended,faces,b,w,ch,ext,limiter,MPI_COMM_WORLD,good, &
            redistribution=method,mol_mode=.true.,mol_lax_friedrichs=lf,frame=frame, &
      state_redist_max_order=requested_order,registered_flux=mpi_registered,integrated_flux=mpi_physical, &
      target_volume_fraction=merge(.8_dp,.5_dp,variant==6))
        end if
      end if
      call require(good,'MPI MOL evolution')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==reference%nodes(p)%state),'exact serial MPI MOL state')
        call require(all(lt%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'exact serial MPI MOL temperature')
        if(variant==3.or.variant==6) then
          call require(all(mpi_registered(p)%x==serial_registered(p)%x).and. &
            all(mpi_registered(p)%y==serial_registered(p)%y).and. &
            all(mpi_registered(p)%z==serial_registered(p)%z),'exact serial MPI MOL registered flux')
          call require(all(mpi_physical(p)%x==serial_physical(p)%x).and. &
            all(mpi_physical(p)%y==serial_physical(p)%y).and. &
            all(mpi_physical(p)%z==serial_physical(p)%z),'exact serial MPI MOL physical flux')
        endif
      end do
      if((variant==3.or.variant==6).and.rank==owners(1)) call require( &
        any(mpi_registered(1)%x/=mpi_physical(1)%x).or. &
        any(mpi_registered(1)%y/=mpi_physical(1)%y).or. &
        any(mpi_registered(1)%z/=mpi_physical(1)%z), &
        'MOL old-half registration remains distinct from physical flux')
      call require(limiter==slimiter,'exact MOL limiter')
      error=max(maxval(abs(b-sb)/max(1._dp,abs(sb))),maxval(abs(w-sw)/max(1._dp,abs(sw))), &
        maxval(abs(ch-sch)/max(1._dp,abs(sch))),maxval(abs(ext-sext)/max(1._dp,abs(sext))))
      call require(error<2.e-12_dp,'serial MPI MOL ledgers')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,final_inventory,MPI_COMM_WORLD,good)
      call require(good,'MOL final inventory')
      error=maxval(abs(final_inventory-initial_inventory-b-w-ch-ext)/max(1._dp,abs(initial_inventory)))
      if(rank==0) print *, 'MOL variant, composite residual: ',variant,error
      call require(error<3.e-14_dp,'MOL composite closure')
      deallocate(local%fields,lt%fields)
    end do
    if(present(frame)) then
      bad_frame=frame
      if(rank==0) bad_frame%angular_speed=bad_frame%angular_speed+1._dp
      owners=[(mod(p-1,nranks),p=1,count)]
      allocate(local%fields(count),lt%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=initial%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4))
        lt%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
      end do
      snapshot=local;snapshot_t=lt
      b=17;w=19;ch=23;ext=29;limiter=.7_dp
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
        1.e-10_dp,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,mol_mode=.true.,frame=bad_frame, &
      state_redist_max_order=requested_order)
      call require(.not.good,'MPI SDC frame mismatch rejected')
      call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.limiter==.7_dp, &
        'MPI SDC frame mismatch ledger rollback')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==snapshot%fields(p)%values).and. &
          all(lt%fields(p)%values==snapshot_t%fields(p)%values),'MPI SDC frame mismatch state rollback')
      end do
    end if
  end subroutine check_mol_parity

  subroutine serial_held_auxiliary(state,interval,stage,depth,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p,x,y,z,n(3)
    real(dp), parameter :: forcing=1.e6_dp
    auxiliary_family_calls=auxiliary_family_calls+1
    good=.true.
    do p=1,size(state%nodes)
      call require(size(state%nodes(p)%state,1)==nv+1,'serial source receives actual held input row')
      rate%fields(p)%values=0
      if(depth>=0.and.state%topology%nodes(p)%level/=depth) cycle
      n=state%topology%nodes(p)%shape
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        if(state%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
        call require((state%nodes(p)%state(nv+1,x,y,z)<0).eqv.auxiliary_family_water, &
          'serial actual held input retains signed formation energy')
        rate%fields(p)%values(5,x,y,z)=merge(forcing,-forcing,auxiliary_family_reacting)
        rate%fields(p)%values(nv+1,x,y,z)=merge(forcing,2*forcing,auxiliary_family_reacting)
      enddo
      enddo
      enddo
    enddo
  end subroutine serial_held_auxiliary

  subroutine mpi_held_auxiliary(state,temp,interval,stage,depth,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temp
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p,x,y,z,n(3)
    real(dp), allocatable :: wrong(:,:,:,:)
    real(dp), parameter :: forcing=1.e6_dp
    auxiliary_family_calls=auxiliary_family_calls+1;good=.true.
    do p=1,size(state%fields)
      if(auxiliary_family_owners(p)/=rank) cycle
      call require(size(state%fields(p)%values,1)==nv+1,'MPI source receives actual held input row')
      rate%fields(p)%values=0
      if(depth>=0.and.topology%nodes(p)%level/=depth) cycle
      n=topology%nodes(p)%shape
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        if(auxiliary_family_geometry(p)%volume_fraction(x,y,z)<=0) cycle
        call require((state%fields(p)%values(nv+1,x,y,z)<0).eqv.auxiliary_family_water, &
          'MPI actual held input retains signed formation energy')
        rate%fields(p)%values(5,x,y,z)=merge(forcing,-forcing,auxiliary_family_reacting)
        rate%fields(p)%values(nv+1,x,y,z)=merge(forcing,2*forcing,auxiliary_family_reacting)
      enddo
      enddo
      enddo
    enddo
    if(auxiliary_family_calls==2) then
      if(auxiliary_family_failure==1.and.rank==nranks-1) good=.false.
      if(auxiliary_family_failure==2.and.rank==auxiliary_family_owners(1)) &
        rate%fields(1)%values(nv+1,1,1,2)=ieee_value(0._dp,ieee_quiet_nan)
      if(auxiliary_family_failure==3.and.rank==auxiliary_family_owners(1)) then
        n=topology%nodes(1)%shape
        allocate(wrong(nv+1,0:n(1)-1,0:n(2)-1,0:n(3)-1),source=0._dp)
        call move_alloc(wrong,rate%fields(1)%values)
      endif
    endif
  end subroutine mpi_held_auxiliary

  subroutine serial_routed_held_auxiliary(state,interval,stage,depth,time,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p,x,y,z,n(3),target
    target=merge(2,1,depth==0);auxiliary_family_calls=auxiliary_family_calls+1;good=.true.
    do p=1,size(state%nodes)
      call require(size(state%nodes(p)%state,1)==nv+1,'routed serial source actual held row')
      rate%fields(p)%values=0;if(p/=target) cycle
      n=state%topology%nodes(p)%shape
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        if(state%nodes(p)%geometry%volume_fraction(x,y,z)<=0.or.covered(state%topology,p,[x,y,z])) cycle
        rate%fields(p)%values(5,x,y,z)=merge(1.e6_dp,-1.e6_dp,auxiliary_family_reacting)
        rate%fields(p)%values(nv+1,x,y,z)=merge(1.e6_dp,2.e6_dp,auxiliary_family_reacting)
      enddo
      enddo
      enddo
    enddo
  end subroutine serial_routed_held_auxiliary

  subroutine mpi_routed_held_auxiliary(state,temp,interval,stage,depth,time,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temp
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p,x,y,z,n(3),target
    target=merge(2,1,depth==0);auxiliary_family_calls=auxiliary_family_calls+1;good=.true.
    do p=1,size(state%fields)
      if(auxiliary_family_owners(p)/=rank) cycle
      call require(size(state%fields(p)%values,1)==nv+1,'routed MPI source actual held row')
      rate%fields(p)%values=0;if(p/=target) cycle
      n=topology%nodes(p)%shape
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        if(auxiliary_family_geometry(p)%volume_fraction(x,y,z)<=0.or.covered(topology,p,[x,y,z])) cycle
        rate%fields(p)%values(5,x,y,z)=merge(1.e6_dp,-1.e6_dp,auxiliary_family_reacting)
        rate%fields(p)%values(nv+1,x,y,z)=merge(1.e6_dp,2.e6_dp,auxiliary_family_reacting)
      enddo
      enddo
      enddo
    enddo
    if(auxiliary_family_failure==1.and.auxiliary_family_calls==6.and.rank==nranks-1) good=.false.
  end subroutine mpi_routed_held_auxiliary

  subroutine mpi_halo_routed_held_auxiliary(state,temp,interval,stage,depth,time,rate,comm,good, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,fractions)
    type(amr_tree_state_3d), intent(in) :: state,temp,ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), intent(in) :: interval,time,fractions(:)
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p,x,y,z,n(3)
    call require(size(fractions)==1.and.all(fractions>=0).and.all(fractions<=1), &
      'paired auxiliary halo callback retains ancestor time fractions')
    if(depth==1) call require(abs(fractions(1)-(time+merge(interval,0._dp,stage>0))/1.e-8_dp)<1.e-14_dp, &
      'paired auxiliary halo callback samples actual old/predictor source clock')
    do p=1,count
      if(auxiliary_family_owners(p)/=rank) cycle
      call require(size(ancestor_start%fields(p)%values,1)==nv+1.and. &
        size(ancestor_end%fields(p)%values,1)==nv+1.and.size(temperature_start%fields(p)%values,1)==1.and. &
        size(temperature_end%fields(p)%values,1)==1,'paired auxiliary halo callback receives actual held ancestor rows')
      if(depth/=1.or.topology%nodes(p)%level/=0) cycle
      n=topology%nodes(p)%shape
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        if(auxiliary_family_geometry(p)%volume_fraction(x,y,z)<=0) cycle
        call require(ancestor_start%fields(p)%values(nv+1,x,y,z)-ancestor_start%fields(p)%values(5,x,y,z)>1.and. &
          ancestor_end%fields(p)%values(nv+1,x,y,z)-ancestor_end%fields(p)%values(5,x,y,z)>1, &
          'auxiliary ancestors carry independent UEINT instead of inferring it from UEDEN')
      enddo
      enddo
      enddo
    enddo
    call mpi_routed_held_auxiliary(state,temp,interval,stage,depth,time,rate,comm,good)
  end subroutine mpi_halo_routed_held_auxiliary

  subroutine check_held_auxiliary_family()
    type(amr_patch_3d) :: box
    type(eb_geometry_3d) :: fine,extended_fine
    type(eb_geometry_3d), allocatable :: extended(:)
    type(eb_reactive_tree_3d) :: seed,reference,routed_reference
    type(amr_tree_state_3d) :: seed_held,reference_held,routed_held,gas,temp,held,old,ot,oh
    type(eb_integrated_flux_3d), allocatable :: faces(:),scalar(:),reference_faces(:),reference_scalar(:)
    type(eb_integrated_flux_3d), allocatable :: kept_faces(:),kept_scalar(:)
    type(state_redist_register_3d), allocatable :: registers(:),held_registers(:)
    type(mpi_eb_mol_interval_3d) :: record
    type(mpi_eb_mol_interval_3d), allocatable :: routed_records(:)
    type(eb_mol_held_interval_3d), allocatable :: serial_records(:),replayed_records(:)
    type(eb_reactive_tree_3d) :: replayed
    type(amr_tree_state_3d) :: replayed_held
    integer(int64) :: source_steps(2)
    real(dp) :: replay_prefix,weight
    type(dual_energy_policy) :: policy
    real(dp) :: q(reactive_nprim(size(species))),v(nv),u0(nv,2,2,2),t0(2,2,2),sound,recovered
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(nv+1),limiter,correction,sb(nv),sw(nv),sch(nv),sext(nv)
    real(dp) :: routed_saux(nv+1)
    real(dp) :: saux(nv+1),slimiter,sprefix,before(nv),after(nv),volume,error,expected_source
    real(dp), parameter :: interval=1.e-8_dp,force(3)=0._dp
    character(len=16) :: boundaries(6),method
    integer :: p,x,y,z,n(3),mode,layout,index,failure,trial,r
    logical :: good,reacting,prepared
    boundaries='outflow';count=2
    call initialize_amr_patch_3d(2,2,2,1,1,1,1,1,1,2,box,good)
    call require(good,'paired source partial child')
    call initialize_amr_patch_tree_3d([2,2,2],[1],[box],topology,good)
    call require(good,'paired source two-level tree')
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.31_dp,fine,good)
    call require(good,'paired source cut geometry')
    call build_axis_plane_eb_geometry_3d(28,28,28,-3._dp,4._dp,-3._dp,4._dp,-3._dp,4._dp, &
      'z',.31_dp,extended_fine,good)
    call require(good,'paired source extended geometry')
    do mode=0,2
      auxiliary_family_water=mode==1;reacting=mode==2;prepared=reacting
      auxiliary_family_reacting=reacting;method='none';if(prepared) method='state_redist'
      index=0
      do p=1,size(species)
        if((.not.auxiliary_family_water.and.trim(species(p)%name)=='N2').or. &
          (auxiliary_family_water.and.trim(species(p)%name)=='H2O')) index=p
      enddo
      call require(index>0,'paired source species')
      q=0;q(1)=.31_dp;q(5)=60000;q(5+index)=1
      call reactive_primitive_to_conserved(species,q,v,recovered,sound,good)
      call require(good,'paired source independent primitive seed')
      do z=1,2
      do y=1,2
      do x=1,2
        u0(:,x,y,z)=v;t0(x,y,z)=recovered
      enddo
      enddo
      enddo
      call initialize_reactive_eb_tree_3d(species,topology,fine,u0,t0,seed,good)
      call require(good,'paired source native seed')
      if(allocated(seed_held%fields)) deallocate(seed_held%fields)
      allocate(seed_held%fields(count),auxiliary_family_geometry(count))
      do p=1,count
        auxiliary_family_geometry(p)=seed%nodes(p)%geometry;n=topology%nodes(p)%shape
        allocate(seed_held%fields(p)%values(1,n(1),n(2),n(3)),source=0._dp)
        do z=1,n(3)
        do y=1,n(2)
        do x=1,n(1)
          if(seed%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          seed_held%fields(p)%values(1,x,y,z)=seed%nodes(p)%state(5,x,y,z)
          seed%nodes(p)%state(5,x,y,z)=seed%nodes(p)%state(5,x,y,z)-11
        enddo
        enddo
        enddo
      enddo
      call build_eb_tree_extended_geometry_3d(topology,auxiliary_family_geometry,extended_fine,extended,good)
      call require(good,'paired source grow6 geometry')
      call integrate_reactive_eb_tree_3d(seed,before,good)
      call require(good,'paired source initial inventory')
      volume=.69_dp;policy=dual_energy_policy();policy%eta2=2
      ! The negative-energy policy intentionally projects UEINT from total.
      ! Positive N2 uses the pinned held-energy repair branch instead, so
      ! distinct gas/UEINT sources and the depleted initial total stay visible.
      if(.not.auxiliary_family_water) policy%allow_negative_energy=.false.
      reference=seed;reference_held=seed_held;sprefix=7;auxiliary_family_calls=0
      if(.not.prepared) then
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,2,.false., &
          [.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,reference,extended,boundaries, &
          sb,sw,sch,sext,slimiter,good,mol_mode=.true.,internal_energy=reference_held, &
          internal_energy_flux=reference_scalar,integrated_flux=reference_faces,energy_policy=policy, &
          energy_correction=sprefix,auxiliary_source=serial_held_auxiliary,auxiliary_change=saux, &
          auxiliary_has_internal_energy=.true.)
        call require(good.and.auxiliary_family_calls==2,'serial packed gas/UEINT source shared evaluation')
        call require(abs(saux(nv+1)-2.e6_dp*interval*volume)<1.e-14_dp,'serial independent held volume ledger')
      endif
      do layout=1,2
        auxiliary_family_owners=[(mod(p-1,nranks),p=1,count)];if(layout==2) auxiliary_family_owners=0
        allocate(old%fields(count),ot%fields(count),oh%fields(count))
        do p=1,count
          if(auxiliary_family_owners(p)/=rank) cycle
          old%fields(p)%values=seed%nodes(p)%state;n=topology%nodes(p)%shape
          allocate(ot%fields(p)%values(1,n(1),n(2),n(3)))
          ot%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
          oh%fields(p)=seed_held%fields(p)
        enddo
        gas=old;temp=ot;held=oh;correction=7;auxiliary_family_calls=0;auxiliary_family_failure=0
        b=0;w=0;ch=0;ext=0;aux=0;limiter=1
        if(prepared) then
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,auxiliary_family_owners, &
            auxiliary_family_geometry,gas,temp,extended,interval,2,.true.,[.false.,.false.,.false.,.false.], &
            force,2.e-7_dp,1.e-12_dp,boundaries,b,w,ch,ext,limiter,MPI_COMM_WORLD,good,mol_mode=.true., &
            redistribution=method,target_depth=0,integrated_flux=faces,state_redist_registers=registers, &
            held_state_redist_registers=held_registers,internal_energy=held,internal_energy_flux=scalar, &
            energy_policy=policy,energy_correction=correction,auxiliary_source=mpi_held_auxiliary, &
            auxiliary_change=aux,auxiliary_has_internal_energy=.true.,interval_history=record, &
            evaluate_reaction=held_diffusion_inert_reactor)
        else
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,auxiliary_family_owners, &
            auxiliary_family_geometry,gas,temp,extended,interval,2,.false.,[.false.,.false.,.false.,.false.], &
            force,2.e-7_dp,1.e-12_dp,boundaries,b,w,ch,ext,limiter,MPI_COMM_WORLD,good,mol_mode=.true., &
            internal_energy=held,internal_energy_flux=scalar,integrated_flux=faces,energy_policy=policy, &
            energy_correction=correction,auxiliary_source=mpi_held_auxiliary,auxiliary_change=aux, &
            auxiliary_has_internal_energy=.true.)
        endif
        call require(good.and.auxiliary_family_calls==2,'MPI packed gas/UEINT source shared old/predictor evaluation')
        do p=1,count
          if(auxiliary_family_owners(p)/=rank) cycle
          if(.not.prepared) then
            call require(all(gas%fields(p)%values==reference%nodes(p)%state).and. &
              all(temp%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature).and. &
              all(held%fields(p)%values==reference_held%fields(p)%values),'exact serial/MPI packed source endpoint')
          endif
        enddo
        if(.not.prepared) then
          call require(all(aux==saux),'exact serial/MPI gas and held volume ledgers')
          call integrate_mpi_eb_tree_3d(topology,auxiliary_family_owners,auxiliary_family_geometry, &
            gas,after,MPI_COMM_WORLD,good)
          call require(good,'paired source final inventory')
          after=after-before-b-w-ch-ext-aux(:nv);after(5)=after(5)-(correction-7)
          error=maxval(abs(after)/max(1._dp,abs(before)))
          call require(error<3.e-14_dp,'paired source gas physical/numerical independent inventory')
        else
          call require(record%auxiliary_has_internal_energy.and.record%has_auxiliary_source, &
            'prepared archive retains explicit held source layout')
          do p=1,count
            if(auxiliary_family_owners(p)/=rank.or.topology%nodes(p)%level/=0) cycle
            call require(size(record%prepared_stages(0)%auxiliary_rate%fields(p)%values,1)==nv+1.and. &
              all(record%prepared_stages(2)%auxiliary_rate%fields(p)%values==0), &
              'prepared archive retains packed old source and skips later new source')
          enddo
          gas=old;temp=ot;held=oh;correction=7+record%initial_energy_correction
          do p=1,count
            if(auxiliary_family_owners(p)/=rank.or.topology%nodes(p)%level/=0) cycle
            gas%fields(p)=record%initial%fields(p)
            temp%fields(p)=record%initial_temperature%fields(p)
            held%fields(p)=record%initial_internal_energy%fields(p)
          enddo
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,auxiliary_family_owners, &
            auxiliary_family_geometry,gas,temp,extended,interval,2,.true.,[.false.,.false.,.false.,.false.], &
            force,2.e-7_dp,1.e-12_dp,boundaries,b,w,ch,ext,limiter,MPI_COMM_WORLD,good,mol_mode=.true., &
            redistribution=method,target_depth=0,integrated_flux=faces,state_redist_registers=registers, &
            held_state_redist_registers=held_registers,internal_energy=held,internal_energy_flux=scalar, &
            energy_policy=policy,energy_correction=correction,auxiliary_source=mpi_held_auxiliary, &
            auxiliary_change=aux,auxiliary_has_internal_energy=.true.,replay_interval=record,replay_spatial_weight=0._dp, &
            evaluate_reaction=held_diffusion_inert_reactor)
          call require(good.and.auxiliary_family_calls==2,'zero spatial weight replays frozen held source without callback')
          expected_source=.5e6_dp*interval
          do p=1,count
            if(auxiliary_family_owners(p)/=rank.or.topology%nodes(p)%level/=0) cycle
            n=topology%nodes(p)%shape
            do z=1,n(3)
            do y=1,n(2)
            do x=1,n(1)
              if(auxiliary_family_geometry(p)%volume_fraction(x,y,z)<=0) cycle
              if(abs(held%fields(p)%values(1,x,y,z)-oh%fields(p)%values(1,x,y,z)-expected_source)>= &
                3.e-14_dp*max(1._dp,abs(oh%fields(p)%values(1,x,y,z)))) then
                print *, 'held source replay observed/expected: ',p,x,y,z, &
                  held%fields(p)%values(1,x,y,z)-oh%fields(p)%values(1,x,y,z),expected_source, &
                  aux(nv+1),correction
              endif
              call require(abs(held%fields(p)%values(1,x,y,z)-oh%fields(p)%values(1,x,y,z)-expected_source)< &
                3.e-14_dp*max(1._dp,abs(oh%fields(p)%values(1,x,y,z))), &
                'frozen held volume source survives zero spatial blend')
            enddo
            enddo
            enddo
          enddo
        endif
        if(mode==0) then
          kept_faces=faces;kept_scalar=scalar
          do failure=1,3
            gas=old;temp=ot;held=oh;correction=7;b=17;w=19;ch=23;ext=29;aux=31;limiter=.7_dp
            auxiliary_family_calls=0;auxiliary_family_failure=failure
            call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,auxiliary_family_owners, &
              auxiliary_family_geometry,gas,temp,extended,interval,2,.false.,[.false.,.false.,.false.,.false.], &
              force,2.e-7_dp,1.e-12_dp,boundaries,b,w,ch,ext,limiter,MPI_COMM_WORLD,good,mol_mode=.true., &
              internal_energy=held,internal_energy_flux=scalar,integrated_flux=faces,energy_policy=policy, &
              energy_correction=correction,auxiliary_source=mpi_held_auxiliary,auxiliary_change=aux, &
              auxiliary_has_internal_energy=.true.)
            call require(.not.good.and.correction==7.and.all(b==17).and.all(w==19).and.all(ch==23).and. &
              all(ext==29).and.all(aux==31).and.limiter==.7_dp,'late held source refusal preserves every ledger')
            do p=1,count
              if(auxiliary_family_owners(p)/=rank) cycle
              call require(all(gas%fields(p)%values==old%fields(p)%values).and. &
                all(temp%fields(p)%values==ot%fields(p)%values).and.all(held%fields(p)%values==oh%fields(p)%values), &
                'late held source refusal preserves exact paired fields')
              call require(all(faces(p)%x==kept_faces(p)%x).and.all(faces(p)%y==kept_faces(p)%y).and. &
                all(faces(p)%z==kept_faces(p)%z).and.all(scalar(p)%x==kept_scalar(p)%x).and. &
                all(scalar(p)%y==kept_scalar(p)%y).and.all(scalar(p)%z==kept_scalar(p)%z), &
                'late held source refusal preserves gas and scalar faces')
            enddo
          enddo
        endif
        if(mode==0) then
          gas=old;temp=ot;held=oh;correction=7;auxiliary_family_calls=0;auxiliary_family_failure=0
          routed_reference=seed;routed_held=seed_held;sprefix=7
          call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,interval,2,.false., &
            [.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,routed_reference,extended,boundaries, &
            sb,sw,sch,sext,slimiter,good,mol_mode=.true.,internal_energy=routed_held, &
            internal_energy_flux=reference_scalar,integrated_flux=reference_faces,energy_policy=policy, &
            energy_correction=sprefix,auxiliary_source=serial_routed_held_auxiliary,auxiliary_change=routed_saux, &
            auxiliary_has_internal_energy=.true.,route_auxiliary=.true.)
          call require(good.and.auxiliary_family_calls==6,'serial routed held source common child clocks')
          call require(abs(routed_saux(nv+1)-2.e6_dp*interval*volume)<1.e-14_dp,'routed held physical inventory')
          auxiliary_family_calls=0
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,auxiliary_family_owners, &
            auxiliary_family_geometry,gas,temp,interval,2,.false.,[.false.,.false.,.false.,.false.], &
            force,2.e-7_dp,1.e-12_dp,extended,boundaries,b,w,ch,ext,limiter,MPI_COMM_WORLD,good, &
            mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar,integrated_flux=faces, &
            energy_policy=policy,energy_correction=correction,auxiliary_source=mpi_routed_held_auxiliary, &
            auxiliary_change=aux,auxiliary_has_internal_energy=.true.,route_auxiliary=.true.)
          call require(good.and.auxiliary_family_calls==6,'MPI routed held source common child clocks')
          do p=1,count
            if(auxiliary_family_owners(p)/=rank) cycle
            call require(all(gas%fields(p)%values==routed_reference%nodes(p)%state).and. &
              all(temp%fields(p)%values(1,:,:,:)==routed_reference%nodes(p)%temperature).and. &
              all(held%fields(p)%values==routed_held%fields(p)%values), &
              'exact serial/MPI routed source gas temperature and held scalar')
          enddo
          call require(all(aux==routed_saux),'exact routed physical gas/UEINT ledgers')
          call integrate_mpi_eb_tree_3d(topology,auxiliary_family_owners,auxiliary_family_geometry, &
            gas,after,MPI_COMM_WORLD,good)
          call require(good,'routed source final inventory')
          after=after-before-b-w-ch-ext-aux(:nv);after(5)=after(5)-(correction-7)
          error=maxval(abs(after)/max(1._dp,abs(before)))
          call require(error<3.e-14_dp,'routed source gas and numerical energy composite inventory')
          gas=old;temp=ot;held=oh;correction=7;b=17;w=19;ch=23;ext=29;aux=31;limiter=.7_dp
          auxiliary_family_calls=0;auxiliary_family_failure=1
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,auxiliary_family_owners, &
            auxiliary_family_geometry,gas,temp,interval,2,.false.,[.false.,.false.,.false.,.false.], &
            force,2.e-7_dp,1.e-12_dp,extended,boundaries,b,w,ch,ext,limiter,MPI_COMM_WORLD,good, &
            mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar,integrated_flux=faces, &
            energy_policy=policy,energy_correction=correction,auxiliary_source=mpi_routed_held_auxiliary, &
            auxiliary_change=aux,auxiliary_has_internal_energy=.true.,route_auxiliary=.true.)
          call require(.not.good.and.correction==7.and.all(aux==31).and.all(b==17).and.all(w==19).and. &
            all(ch==23).and.all(ext==29).and.limiter==.7_dp,'late descendant routed source refusal keeps every ledger')
          do p=1,count
            if(auxiliary_family_owners(p)/=rank) cycle
            call require(all(gas%fields(p)%values==old%fields(p)%values).and. &
              all(temp%fields(p)%values==ot%fields(p)%values).and.all(held%fields(p)%values==oh%fields(p)%values), &
              'late descendant routed source refusal keeps paired fields')
          enddo
          auxiliary_family_failure=0;auxiliary_family_calls=0
          routed_reference=seed;routed_held=seed_held
          call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,interval,2,.false., &
            [.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,routed_reference,extended,boundaries, &
            sb,sw,sch,sext,slimiter,good,mol_mode=.true.,internal_energy=routed_held, &
            internal_energy_flux=reference_scalar,auxiliary_source=serial_routed_held_auxiliary, &
            auxiliary_change=routed_saux,auxiliary_has_internal_energy=.true.,route_auxiliary=.true.)
          call require(good.and.auxiliary_family_calls==6,'serial routed held sources permit an absent numerical ledger')
          gas=old;temp=ot;held=oh;auxiliary_family_calls=0
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,auxiliary_family_owners, &
            auxiliary_family_geometry,gas,temp,interval,2,.false.,[.false.,.false.,.false.,.false.], &
            force,2.e-7_dp,1.e-12_dp,extended,boundaries,b,w,ch,ext,limiter,MPI_COMM_WORLD,good, &
            mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar, &
            halo_source=mpi_halo_routed_held_auxiliary,auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
            route_auxiliary=.true.)
          call require(good.and.auxiliary_family_calls==6,'MPI routed held sources omit numerical reset without a ledger')
          do p=1,count
            if(auxiliary_family_owners(p)/=rank) cycle
            call require(all(gas%fields(p)%values==routed_reference%nodes(p)%state).and. &
              all(temp%fields(p)%values(1,:,:,:)==routed_reference%nodes(p)%temperature).and. &
              all(held%fields(p)%values==routed_held%fields(p)%values), &
              'routed gas temperature and UEINT agree exactly without a numerical ledger')
            n=topology%nodes(p)%shape
            do z=1,n(3)
            do y=1,n(2)
            do x=1,n(1)
              if(auxiliary_family_geometry(p)%volume_fraction(x,y,z)<=0) cycle
              call require(abs(held%fields(p)%values(1,x,y,z)-oh%fields(p)%values(1,x,y,z)-2.e6_dp*interval)< &
                3.e-14_dp*max(1._dp,abs(oh%fields(p)%values(1,x,y,z))), &
                'independent routed UEINT volume survives without total-energy reset')
              call temperature_from_internal_energy(species,gas%fields(p)%values(6:,x,y,z)/ &
                gas%fields(p)%values(1,x,y,z),held%fields(p)%values(1,x,y,z)/gas%fields(p)%values(1,x,y,z), &
                ot%fields(p)%values(1,x,y,z),recovered,good)
              call require(good.and.abs(temp%fields(p)%values(1,x,y,z)-recovered)< &
                3.e-14_dp*max(1._dp,abs(recovered)),'routed no-ledger temperature uses actual held UEINT')
            enddo
            enddo
            enddo
          enddo
          call integrate_mpi_eb_tree_3d(topology,auxiliary_family_owners,auxiliary_family_geometry, &
            gas,after,MPI_COMM_WORLD,good)
          call require(good.and.all(aux==routed_saux),'routed source no-ledger physical inventory and budgets')
          after=after-before-b-w-ch-ext-aux(:nv)
          call require(maxval(abs(after)/max(1._dp,abs(before)))<3.e-14_dp, &
            'routed source absent-ledger gas conservation has no numerical correction')
        endif
        if(mode==2) then
          routed_reference=seed;routed_held=seed_held;sprefix=7
          auxiliary_family_calls=0;auxiliary_family_failure=0
          call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,interval,2,.true., &
            [.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,routed_reference,extended,boundaries, &
            sb,sw,sch,sext,slimiter,good,mol_mode=.true.,redistribution='state_redist', &
            internal_energy=routed_held,internal_energy_flux=reference_scalar,integrated_flux=reference_faces, &
            energy_policy=policy,energy_correction=sprefix,auxiliary_source=serial_routed_held_auxiliary, &
            auxiliary_change=routed_saux,auxiliary_has_internal_energy=.true.,route_auxiliary=.true., &
            held_interval_history=serial_records,level_steps=source_steps,evaluate_reaction=held_diffusion_inert_reactor)
          call require(good.and.auxiliary_family_calls==6.and.all(source_steps==[1_int64,2_int64]), &
            'prepared serial routed source captures actual old/predictor physical clocks')
          call require(size(serial_records)==3.and.all(serial_records%routed_auxiliary).and. &
            all(serial_records%auxiliary_has_internal_energy),'prepared routed archive retains paired volume layout')
          call require(abs(routed_saux(nv+1)-.5e6_dp*interval*volume)<1.e-14_dp, &
            'reacting routed old-only source budget excludes later MOL auxiliary stages')
          gas=old;temp=ot;held=oh;correction=7;auxiliary_family_calls=0
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,auxiliary_family_owners, &
            auxiliary_family_geometry,gas,temp,interval,2,.true.,[.false.,.false.,.false.,.false.], &
            force,2.e-7_dp,1.e-12_dp,extended,boundaries,b,w,ch,ext,limiter,MPI_COMM_WORLD,good, &
            mol_mode=.true.,redistribution='state_redist',internal_energy=held,internal_energy_flux=scalar, &
            integrated_flux=faces,energy_policy=policy,energy_correction=correction, &
            auxiliary_source=mpi_routed_held_auxiliary,auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
            route_auxiliary=.true.,interval_history=routed_records,level_steps=source_steps, &
            evaluate_reaction=held_diffusion_inert_reactor)
          call require(good.and.auxiliary_family_calls==6.and.all(source_steps==[1_int64,2_int64]), &
            'prepared MPI routed source retains actual source clocks')
          call require(size(routed_records)==3.and.all(routed_records%routed_auxiliary).and. &
            all(routed_records%auxiliary_has_internal_energy),'MPI paired routed archives retain complete source flags')
          do p=1,count
            if(auxiliary_family_owners(p)/=rank) cycle
            call require(all(gas%fields(p)%values==routed_reference%nodes(p)%state).and. &
              all(temp%fields(p)%values(1,:,:,:)==routed_reference%nodes(p)%temperature).and. &
              all(held%fields(p)%values==routed_held%fields(p)%values), &
              'prepared routed serial/MPI gas temperature and UEINT agree exactly')
            do r=1,3
              call require(size(routed_records(r)%routed_source%fields(p)%values,1)==nv+1.and. &
                all(routed_records(r)%routed_source%fields(p)%values==serial_records(r)%routed_source%fields(p)%values), &
                'prepared serial/MPI frozen routed impulse agrees exactly')
            enddo
          enddo
          call require(all(aux==routed_saux),'prepared serial/MPI packed physical source budget agrees exactly')
          do trial=1,2
            weight=merge(1._dp,0._dp,trial==1)
            replayed=routed_reference;replayed_held=routed_held;replay_prefix=29
            auxiliary_family_calls=0
            call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,interval,2,.true., &
              [.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,replayed,extended,boundaries, &
              b,w,ch,ext,limiter,good,mol_mode=.true.,redistribution='state_redist', &
              internal_energy=replayed_held,internal_energy_flux=scalar,integrated_flux=faces, &
              energy_policy=policy,energy_correction=replay_prefix,auxiliary_source=serial_routed_held_auxiliary, &
              auxiliary_change=aux,auxiliary_has_internal_energy=.true.,route_auxiliary=.true., &
              held_interval_history=replayed_records,replay_interval_history=serial_records,replay_spatial_weight=weight, &
              level_steps=source_steps,evaluate_reaction=held_diffusion_inert_reactor)
            call require(good.and.auxiliary_family_calls==0.and.all(replayed_records%replayed).and. &
              all(source_steps==[1_int64,2_int64]),'complete prepared routed replay freezes callbacks and source clocks')
            call require(all(aux==routed_saux),'frozen routed volume budget survives zero spatial weight')
            call integrate_reactive_eb_tree_3d(replayed,after,good)
            call require(good,'prepared routed replay physical inventory')
            after=after-before-b-w-ch-ext-aux(:nv);after(5)=after(5)-(replay_prefix-7)
            error=maxval(abs(after)/max(1._dp,abs(before)))
            call require(error<3.e-14_dp,'complete paired routed replay physical/numerical inventory')
            if(trial/=1) cycle
            do p=1,count
              call require(all(replayed%nodes(p)%state==routed_reference%nodes(p)%state).and. &
                all(replayed%nodes(p)%temperature==routed_reference%nodes(p)%temperature).and. &
                all(replayed_held%fields(p)%values==routed_held%fields(p)%values), &
                'unit weight paired routed complete replay exactly reproduces native fields')
            enddo
          enddo
        endif
        deallocate(old%fields,ot%fields,oh%fields)
      enddo
      deallocate(auxiliary_family_geometry)
    enddo
    auxiliary_family_failure=0
  end subroutine check_held_auxiliary_family

  subroutine held_diffusion_inert_reactor(species,reactions,initial_energy_density,external_species,external_energy, &
      duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
    real(dp), intent(inout) :: species_density(:),temperature
    integer, intent(out) :: completed_steps
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    real(dp) :: candidate(size(species)),density,energy,next_temperature
    integer :: p
    ! This structural callback test uses the exact forced ODE solution for
    ! pure N2, whose chemical rates vanish. Native reaction integration is
    ! qualified separately by the application lifecycle checks.
    completed_steps=0;ok=.false.
    if(size(species_density)/=size(species).or.size(external_species)/=size(species)) return
    if(.not.all(ieee_is_finite([initial_energy_density,external_energy,duration,temperature]))) return
    if(duration<=0) return
    do p=1,size(species)
      if(trim(species(p)%name)=='N2') cycle
      if(species_density(p)/=0.or.external_species(p)/=0) return
    enddo
    candidate=species_density+duration*external_species
    if(.not.all(ieee_is_finite(candidate)).or.any(candidate<0)) return
    density=sum(candidate);if(density<=0) return
    energy=(initial_energy_density+duration*external_energy)/density
    call temperature_from_internal_energy(species,candidate/density,energy,temperature,next_temperature,ok)
    if(.not.ok) return
    species_density=candidate;temperature=next_temperature;completed_steps=1
  end subroutine held_diffusion_inert_reactor

  subroutine held_valid_diffusion(state,temp,interval,depth,rate,flux,limiter,comm,good, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,fractions)
    type(amr_tree_state_3d), intent(in) :: state,temp
    real(dp), intent(in) :: interval
    integer, intent(in) :: depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: limiter
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: fractions(:)
    integer :: p,n(3)
    ! Constant x-directed total-energy flux has zero valid divergence for
    ! this plane geometry. UEINT has no extra diffusion source in this API.
    call require(.not.allocated(flux),'additional diffusion receives independent empty face storage')
    allocate(flux(size(state%fields)))
    do p=1,size(state%fields)
      if(.not.allocated(state%fields(p)%values)) cycle
      rate%fields(p)%values=0
      if(depth>=0.and.topology%nodes(p)%level/=depth) cycle
      n=topology%nodes(p)%shape
      allocate(flux(p)%x(nv,0:n(1),n(2),n(3)),flux(p)%y(nv,n(1),0:n(2),n(3)), &
        flux(p)%z(nv,n(1),n(2),0:n(3)),source=0._dp)
      flux(p)%x(5,:,:,:)=.03125_dp*interval
    enddo
    limiter=1;good=.true.;diffusion_calls=diffusion_calls+1
    if(held_diffusion_failure.and.rank==nranks-1.and.diffusion_calls==2) good=.false.
  end subroutine held_valid_diffusion

  subroutine check_held_valid_diffusion(initial,geometry,extended,faces)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    type(eb_reactive_tree_3d) :: seed
    type(amr_tree_state_3d) :: seed_held,gas,temp,held,old,ot,oh,baseline_gas,baseline_temp,baseline_held
    type(eb_integrated_flux_3d), allocatable :: physical(:),scalar(:),kept_faces(:),kept_scalar(:)
    type(eb_integrated_flux_3d), allocatable :: baseline_faces(:),baseline_scalar(:)
    type(dual_energy_policy) :: policy
    real(dp) :: q(reactive_nprim(size(species))),v(nv),t0,sound,b(nv),w(nv),ch(nv),ext(nv),limiter,correction
    real(dp) :: before(nv),after(nv),error,flux_norm,global_norm
    real(dp), parameter :: interval=1.e-10_dp,force(3)=0._dp
    integer :: owners(count),p,x,y,z,n2,layout,mode,expected
    logical :: good,recursive,classified,flags(4)
    n2=0
    do p=1,size(species)
      if(trim(species(p)%name)=='N2') n2=p
    enddo
    call require(n2>0,'held valid diffusion inert species')
    q=0;q(1)=.31_dp;q(5)=60000;q(5+n2)=1
    call reactive_primitive_to_conserved(species,q,v,t0,sound,good)
    call require(good,'held valid diffusion independent inert seed')
    seed=initial
    allocate(seed_held%fields(count))
    do p=1,count
      allocate(seed_held%fields(p)%values(1,4,4,4),source=0._dp)
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        seed%nodes(p)%state(:,x,y,z)=v;seed%nodes(p)%temperature(x,y,z)=t0
        seed%nodes(p)%state(5,x,y,z)=v(5)-11
        seed_held%fields(p)%values(1,x,y,z)=v(5)
      enddo
      enddo
      enddo
    enddo
    policy%eta2=2;policy%allow_negative_energy=.false.
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)];if(layout==2) owners=0
      allocate(old%fields(count),ot%fields(count),oh%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        old%fields(p)%values=seed%nodes(p)%state
        allocate(ot%fields(p)%values(1,4,4,4))
        ot%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
        oh%fields(p)=seed_held%fields(p)
      enddo
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,old,before,MPI_COMM_WORLD,good)
      call require(good,'held valid diffusion initial inventory')
      do mode=0,3
        recursive=btest(mode,0);classified=btest(mode,1);flags=btest(mode,0).eqv.btest(mode,1)
        baseline_gas=old;baseline_temp=ot;baseline_held=oh;correction=7
        b=0;w=0;ch=0;ext=0;limiter=1
        if(recursive) then
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry, &
            baseline_gas,baseline_temp,interval,2,.true.,flags,force,2.e-7_dp,1.e-12_dp,extended,faces, &
            b,w,ch,ext,limiter,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=baseline_held, &
            internal_energy_flux=baseline_scalar,integrated_flux=baseline_faces, &
            energy_policy=policy,energy_correction=correction,evaluate_reaction=held_diffusion_inert_reactor)
        else
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry, &
            baseline_gas,baseline_temp,extended,interval,2,.true.,flags,force,2.e-7_dp,1.e-12_dp,faces, &
            b,w,ch,ext,limiter,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=baseline_held, &
            internal_energy_flux=baseline_scalar,integrated_flux=baseline_faces, &
            energy_policy=policy,energy_correction=correction,evaluate_reaction=held_diffusion_inert_reactor)
        endif
        call require(good,'held diffusion independent native no-callback baseline')
        gas=old;temp=ot;held=oh;correction=7;diffusion_calls=0;held_diffusion_failure=.false.
        b=0;w=0;ch=0;ext=0;limiter=1
        if(recursive) then
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
            interval,2,.true.,flags,force,2.e-7_dp,1.e-12_dp,extended,faces, &
            b,w,ch,ext,limiter,MPI_COMM_WORLD,good,mol_mode=.true.,evaluate_diffusion=held_valid_diffusion, &
            internal_energy=held,internal_energy_flux=scalar,integrated_flux=physical, &
            energy_policy=policy,energy_correction=correction,spatial_diffusion_is_source=classified, &
            evaluate_reaction=held_diffusion_inert_reactor)
        else
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
            interval,2,.true.,flags,force,2.e-7_dp,1.e-12_dp,faces, &
            b,w,ch,ext,limiter,MPI_COMM_WORLD,good,mol_mode=.true.,evaluate_diffusion=held_valid_diffusion, &
            internal_energy=held,internal_energy_flux=scalar,integrated_flux=physical, &
            energy_policy=policy,energy_correction=correction,spatial_diffusion_is_source=classified, &
            evaluate_reaction=held_diffusion_inert_reactor)
        endif
        call require(good,'held valid-only diffusion accepted with optional molecular transport')
        expected=merge(2,3,classified)*merge(7,1,recursive)
        call require(diffusion_calls==expected,'held valid-only diffusion stage/recursive forwarding policy')
        flux_norm=0
        do p=1,count
          call require(allocated(held%fields(p)%values).eqv.(owners(p)==rank),'held callback retains owner-local scalar')
          if(owners(p)/=rank) cycle
          flux_norm=flux_norm+sum(abs(physical(p)%x(5,:,:,:)-baseline_faces(p)%x(5,:,:,:)))
          call require(all(held%fields(p)%values==baseline_held%fields(p)%values).and. &
            all(scalar(p)%x==baseline_scalar(p)%x).and.all(scalar(p)%y==baseline_scalar(p)%y).and. &
            all(scalar(p)%z==baseline_scalar(p)%z), &
            'total-energy-only callback adds no UEINT source or face beyond native paired baseline')
        enddo
        call MPI_Allreduce(flux_norm,global_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
        call require(global_norm>0,'nonzero callback physical faces survive MOL assembly')
        call integrate_mpi_eb_tree_3d(topology,owners,geometry,gas,after,MPI_COMM_WORLD,good)
        call require(good,'held valid diffusion final inventory')
        after=after-before-b-w-ch-ext;after(5)=after(5)-(correction-7)
        error=maxval(abs(after)/max(1._dp,abs(before)))
        call require(error<3.e-14_dp,'held valid-only diffusion independent composite balance')
        ! Refuse on an idle rank after successful stage 0; every payload and
        ! previously published face family must survive the outer transaction.
        kept_faces=physical;kept_scalar=scalar
        gas=old;temp=ot;held=oh;correction=7;b=17;w=19;ch=23;ext=29;limiter=.7_dp
        held_diffusion_failure=.true.;diffusion_calls=0
        if(recursive) then
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
            interval,2,.true.,flags,force,2.e-7_dp,1.e-12_dp,extended,faces, &
            b,w,ch,ext,limiter,MPI_COMM_WORLD,good,mol_mode=.true.,evaluate_diffusion=held_valid_diffusion, &
            internal_energy=held,internal_energy_flux=scalar,integrated_flux=physical, &
            energy_policy=policy,energy_correction=correction,spatial_diffusion_is_source=classified, &
            evaluate_reaction=held_diffusion_inert_reactor)
        else
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
            interval,2,.true.,flags,force,2.e-7_dp,1.e-12_dp,faces, &
            b,w,ch,ext,limiter,MPI_COMM_WORLD,good,mol_mode=.true.,evaluate_diffusion=held_valid_diffusion, &
            internal_energy=held,internal_energy_flux=scalar,integrated_flux=physical, &
            energy_policy=policy,energy_correction=correction,spatial_diffusion_is_source=classified, &
            evaluate_reaction=held_diffusion_inert_reactor)
        endif
        call require(.not.good.and.correction==7.and.all(b==17).and.all(w==19).and.all(ch==23).and. &
          all(ext==29).and.limiter==.7_dp,'late callback refusal retains numerical prefix and physical ledgers')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(gas%fields(p)%values==old%fields(p)%values).and. &
            all(temp%fields(p)%values==ot%fields(p)%values).and.all(held%fields(p)%values==oh%fields(p)%values), &
            'late callback refusal retains exact gas temperature UEINT')
          call require(all(physical(p)%x==kept_faces(p)%x).and.all(physical(p)%y==kept_faces(p)%y).and. &
            all(physical(p)%z==kept_faces(p)%z).and.all(scalar(p)%x==kept_scalar(p)%x).and. &
            all(scalar(p)%y==kept_scalar(p)%y).and.all(scalar(p)%z==kept_scalar(p)%z), &
            'late callback refusal retains previously published gas and UEINT faces')
        enddo
      enddo
      deallocate(old%fields,ot%fields,oh%fields)
    enddo
    held_diffusion_failure=.false.
    if(rank==0) print *, 'Held valid-only diffusion common/recursive stage policy and rollback PASS'
  end subroutine check_held_valid_diffusion

  subroutine check_held_common_policy(initial,geometry,extended,faces)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    type(dual_energy_policy) :: policy,rank_policy
    type(eb_reactive_tree_3d) :: seed,reference,failed
    type(amr_tree_state_3d) :: seed_held,reference_held,gas,temp,held,inventory_gas,molecular_rate
    type(eb_integrated_flux_3d), allocatable :: reference_flux(:),scalar_flux(:),molecular_faces(:)
    type(eb_geometry_3d), allocatable :: invalid_geometry(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),theta,sb(nv),sw(nv),sch(nv),sext(nv),stheta
    real(dp) :: before(nv),after(nv),error,correction,reference_correction,uniform(nv),uniform_t,energy,volume
    real(dp), parameter :: interval=1.e-7_dp,force(3)=[2.e9_dp,-3.e9_dp,0._dp]
    integer :: owners(count),p,x,y,z,layout
    logical :: good,reacting,flags(4),converted
    real(dp) :: uniform_pv(nv),local_pv(nv),local_sound,target_t
    reacting=test_mode=='held_chem_policy';flags=reacting
    if(reacting.and.rank==0) print *, 'Held coupled seed start: ',MPI_Wtime()
    seed=initial;uniform=initial%nodes(1)%state(:,1,1,4);uniform_t=initial%nodes(1)%temperature(1,1,4)
    call require(geometry(1)%volume_fraction(1,1,4)>0,'held policy uniform fluid seed')
    energy=uniform(5)-.5_dp*sum(uniform(2:4)**2)/uniform(1)
    if(reacting) then
      call reactive_conserved_to_primitive(species,uniform,uniform_t,uniform_pv,target_t,local_sound,good)
      call require(good,'held coupled fixture primitive')
    endif
    allocate(seed_held%fields(count))
    volume=0
    do p=1,count
      allocate(seed_held%fields(p)%values(1,4,4,4),source=0._dp)
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        seed%nodes(p)%state(:,x,y,z)=uniform;seed%nodes(p)%temperature(x,y,z)=uniform_t
        if(reacting) then
          local_pv=uniform_pv
          target_t=uniform_t*(1+.01_dp*(geometry(p)%x_lower+(x-.5_dp)*geometry(p)%dx)+ &
            .02_dp*(geometry(p)%z_lower+(z-.5_dp)*geometry(p)%dz))
          local_pv(5)=uniform_pv(5)*target_t/uniform_t
          local_pv(2)=5*(geometry(p)%x_lower+(x-.5_dp)*geometry(p)%dx)
          local_pv(3)=-3*(geometry(p)%x_lower+(x-.5_dp)*geometry(p)%dx)
          call reactive_primitive_to_conserved(species,local_pv,seed%nodes(p)%state(:,x,y,z), &
            seed%nodes(p)%temperature(x,y,z),local_sound,good)
          call require(good,'held coupled nonuniform molecular fixture')
          seed_held%fields(p)%values(1,x,y,z)=seed%nodes(p)%state(5,x,y,z)- &
            .5_dp*sum(seed%nodes(p)%state(2:4,x,y,z)**2)/seed%nodes(p)%state(1,x,y,z)
        else
          seed_held%fields(p)%values(1,x,y,z)=energy
        endif
        seed%nodes(p)%state(5,x,y,z)=seed%nodes(p)%state(5,x,y,z)-11
        if(.not.covered(topology,p,[x,y,z])) volume=volume+geometry(p)%volume_fraction(x,y,z)* &
          geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
      enddo
      enddo
      enddo
    enddo
    call integrate_reactive_eb_tree_3d(seed,before,good)
    call require(good,'held policy initial depleted gas inventory')
    if(reacting) then
      allocate(inventory_gas%fields(count))
      do p=1,count
        inventory_gas%fields(p)%values=seed%nodes(p)%state
      enddo
      call held_policy_compensated_inventory(inventory_gas,geometry,before)
    endif
    if(reacting) then
      if(rank==0) print *, 'Held coupled molecular source start: ',MPI_Wtime()
      call compute_sdc_eb_tree_transport_3d(species,transport,interval,.true.,.true.,.true.,.true., &
        seed,molecular_rate,molecular_faces,theta,good,boundary_faces=faces)
      call require(good,'held coupled fixture molecular spatial source')
      converted=.false.
      do p=1,count
        converted=converted.or.any(molecular_faces(p)%x(5,:,:,:)/=0).or. &
          any(molecular_faces(p)%y(5,:,:,:)/=0).or.any(molecular_faces(p)%z(5,:,:,:)/=0)
      enddo
      call require(converted,'held coupled fixture has actual molecular energy flux')
      if(rank==0) print *, 'Held coupled molecular source complete: ',MPI_Wtime()
    endif
    policy%allow_negative_energy=.false.;policy%eta2=2
    reference=seed;reference_held=seed_held;reference_correction=7
    if(reacting) then
      if(rank==0) print *, 'Held coupled serial MOL start: ',MPI_Wtime()
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,2,reacting, &
      flags,force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,stheta,good,mol_mode=.true.,internal_energy=reference_held, &
      internal_energy_flux=reference_flux,energy_policy=policy,energy_correction=reference_correction, &
      evaluate_reaction=reactor_adapter)
    else
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,2,reacting, &
      flags,force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,stheta,good,mol_mode=.true.,internal_energy=reference_held, &
      internal_energy_flux=reference_flux,energy_policy=policy,energy_correction=reference_correction)
    endif
    call require(good,'serial held policy common MOL step')
    call integrate_reactive_eb_tree_3d(reference,after,good)
    call require(good,'serial held policy final inventory')
    if(reacting) then
      do p=1,count
        inventory_gas%fields(p)%values=reference%nodes(p)%state
      enddo
      call held_policy_compensated_inventory(inventory_gas,geometry,after)
    endif
    after=after-before-sb-sw-sch-sext;after(5)=after(5)-(reference_correction-7)
    error=maxval(abs(after)/max(1._dp,abs(before)))
    print *, 'Held common policy serial closure and correction: ',error,reference_correction
    call require(error<3.e-14_dp,'serial held policy independent numerical-energy budget')
    if(reacting) then
      converted=.false.
      do p=1,count
        converted=converted.or.any(reference%nodes(p)%state(6:,:,:,:)/=seed%nodes(p)%state(6:,:,:,:))
      enddo
      call require(converted.and.any(sch(6:)/=0),'held coupled common step performs chemical conversion')
      call require(abs(sch(5))/max(1._dp,abs(before(5)))<3.e-14_dp, &
        'held coupled repair is separate from chemical total-energy exchange')
    else
      call require(abs(reference_correction-7-11*volume)<3.e-14_dp*max(1._dp,abs(before(5))), &
        'predictor and final forced kinetic reset corrections cancel with pinned weights')
      call require(maxval(abs(sch)/max(1._dp,abs(before)))<3.e-14_dp, &
        'held numerical-energy repair is not chemical exchange')
    endif
    invalid_geometry=extended;invalid_geometry(count)%nx=invalid_geometry(count)%nx-1
    failed=seed;held=seed_held;correction=7;b=17;w=19;ch=23;ext=29;theta=.7_dp
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,2,reacting, &
      flags,force,2.e-7_dp,1.e-12_dp,failed,invalid_geometry,faces, &
      b,w,ch,ext,theta,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux, &
      energy_policy=policy,energy_correction=correction)
    call require(.not.good.and.correction==7,'serial failure after initial reset retains numerical ledger')
    call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.theta==.7_dp.and. &
      .not.allocated(scalar_flux),'serial held policy failure retains diagnostics and unallocated flux')
    do p=1,count
      call require(all(failed%nodes(p)%state==seed%nodes(p)%state).and. &
        all(failed%nodes(p)%temperature==seed%nodes(p)%temperature).and. &
        all(held%fields(p)%values==seed_held%fields(p)%values),'serial initial reset rejection keeps gas and held scalar')
    enddo
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)];if(layout==2) owners=0
      allocate(gas%fields(count),temp%fields(count))
      if(allocated(held%fields)) deallocate(held%fields)
      allocate(held%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        gas%fields(p)%values=seed%nodes(p)%state
        allocate(temp%fields(p)%values(1,4,4,4));temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
        held%fields(p)=seed_held%fields(p)
      enddo
      correction=7
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
        interval,2,reacting,flags,force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux, &
        energy_policy=policy,energy_correction=correction)
      call require(good,'MPI held policy common MOL step')
      call require(abs(correction-reference_correction)<3.e-14_dp*max(1._dp,abs(before(5))), &
        'MPI held policy composite correction matches serial')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(gas%fields(p)%values==reference%nodes(p)%state).and. &
          all(temp%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature).and. &
          all(held%fields(p)%values==reference_held%fields(p)%values),'MPI held policy exact owner/idle state')
        call require(all(scalar_flux(p)%x==reference_flux(p)%x).and.all(scalar_flux(p)%y==reference_flux(p)%y).and. &
          all(scalar_flux(p)%z==reference_flux(p)%z),'MPI held policy exact scalar face family')
      enddo
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,gas,after,MPI_COMM_WORLD,good)
      call require(good,'MPI held policy final inventory')
      if(reacting) call held_policy_compensated_inventory(gas,geometry,after,MPI_COMM_WORLD)
      after=after-before-b-w-ch-ext;after(5)=after(5)-(correction-7)
      error=maxval(abs(after)/max(1._dp,abs(before)))
      call require(error<3.e-14_dp,'MPI held policy independent numerical-energy budget')
      do p=1,count
        if(owners(p)/=rank) cycle
        gas%fields(p)%values=seed%nodes(p)%state
        temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
        held%fields(p)=seed_held%fields(p)
      enddo
      correction=7;b=17;w=19;ch=23;ext=29;theta=.7_dp
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,invalid_geometry, &
        interval,2,reacting,flags,force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux, &
        energy_policy=policy,energy_correction=correction)
      call require(.not.good.and.correction==7,'MPI initial reset geometry rejection retains numerical ledger')
      if(reacting) then
        reactor_calls=0;reactor_failure=.true.;reactor_failure_rank=0
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
          interval,2,reacting,flags,force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta,MPI_COMM_WORLD,good, &
          mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux, &
          energy_policy=policy,energy_correction=correction,evaluate_reaction=reactor_adapter)
        call require(.not.good.and.correction==7,'held coupled native reactor refusal retains numerical ledger')
        if(rank==0) call require(reactor_calls==2,'held coupled refusal is not retried')
        reactor_failure=.false.
      endif
      rank_policy=policy;if(rank==nranks-1) rank_policy%eta2=3
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
        interval,2,reacting,flags,force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux, &
        energy_policy=rank_policy,energy_correction=correction)
      call require(.not.good.and.correction==7,'MPI policy mismatch including idle rank rejects collectively')
      call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.theta==.7_dp, &
        'MPI held policy failure retains all diagnostics')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(gas%fields(p)%values==seed%nodes(p)%state).and. &
          all(temp%fields(p)%values(1,:,:,:)==seed%nodes(p)%temperature).and. &
          all(held%fields(p)%values==seed_held%fields(p)%values),'MPI held policy failure keeps all owner data')
        call require(all(scalar_flux(p)%x==reference_flux(p)%x).and.all(scalar_flux(p)%y==reference_flux(p)%y).and. &
          all(scalar_flux(p)%z==reference_flux(p)%z),'MPI held policy failure keeps published scalar faces')
      enddo
      deallocate(gas%fields,temp%fields,held%fields)
    enddo
  end subroutine check_held_common_policy

  subroutine check_held_prepared(initial,geometry,extended,faces)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    type(eb_reactive_tree_3d) :: seed,reference,kept_reference
    type(amr_tree_state_3d) :: seed_held,gas,temp,held,ir,reconstructed,scalar_rate,kept_gas,kept_t,kept_held,kept_ir
    type(mpi_eb_mol_interval_3d) :: record,replayed,kept,bad_record
    type(eb_mol_held_interval_3d) :: serial_replayed
    type(eb_mol_held_interval_3d), allocatable :: serial_frozen
    type(mpi_eb_mol_interval_3d), allocatable :: recursive_records(:),kept_records(:)
    type(eb_mol_held_interval_3d), allocatable :: serial_records(:),kept_serial_records(:)
    type(amr_tree_state_3d) :: heat,kept_heat,reference_ir,reference_heat
    type(eb_integrated_flux_3d), allocatable :: registered(:),kept_faces(:)
    type(eb_integrated_flux_3d), allocatable :: reference_faces(:),reference_scalar_faces(:),reference_registered(:)
    type(amr_tree_state_3d) :: reference_held,discarded,spill,discarded_increment
    type(state_redist_register_3d), allocatable :: discarded_registers(:)
    type(eb_integrated_flux_3d), allocatable :: zero_faces(:),discarded_faces(:)
    real(dp) :: discarded_inventory(nv),discarded_total(nv)
    integer :: diagnostic_depth,first_depth,funding_fault
    logical :: has_spill
    integer(int64) :: steps(3),reference_steps(3)
    real(dp) :: v0(nv),v1(nv),error,reference_correction
    real(dp) :: rb(nv),rw(nv),rch(nv),rext(nv),rtheta
    character(len=16) :: method
    type(state_redist_register_3d), allocatable :: registers(:),held_registers(:),rebuilt(:),held_rebuilt(:)
    type(state_redist_register_3d), allocatable :: replay_regs(:),replay_held_regs(:),kept_held_regs(:)
    type(state_redist_register_3d), allocatable :: reference_regs(:),reference_held_regs(:)
    type(state_redist_register_3d), allocatable :: kept_replay_regs(:),kept_rebuilt(:),kept_held_rebuilt(:)
    type(eb_integrated_flux_3d), allocatable :: kept_replay_faces(:),kept_replay_scalar_faces(:)
    type(state_redist_register_3d), allocatable :: registered_regs(:),registered_held_regs(:)
    type(state_redist_register_3d), allocatable :: serial_registered(:),serial_held_registered(:)
    type(eb_mol_held_interval_3d) :: bad_serial
    type(state_redist_register_3d) :: partial_registration
    real(dp), allocatable :: malformed_delta(:,:,:,:)
    integer :: malformed_dims(3)
    type(amr_tree_state_3d) :: kept_rate,kept_scalar_rate
    type(eb_integrated_flux_3d), allocatable :: flux(:),held_faces(:),replay_faces(:),replay_scalar_faces(:)
    type(dual_energy_policy) :: policy
    real(dp) :: pv(nv),state(nv),t,sound,point(3),b(nv),w(nv),ch(nv),ext(nv),theta,correction
    real(dp) :: local_norm,global_norm,weight,kept_correction,kept_chemical(nv),case_force(3)
    integer :: owners(count),p,x,y,z,sp,layout,stage,trial,recursive_iterations
    character(len=32) :: recursive_diagnostic
    logical :: good,refused
    real(dp), parameter :: dt=1.e-10_dp,force(3)=[2.e4_dp,-3.e4_dp,4.e4_dp]
    method='state_redist';case_force=force
    ! Serial and MPI replacement retain the same actual volume stages even
    ! when the weighted spatial update changes force-work at the new endpoint.
    if(test_mode=='held_recursive_chemistry') method='none'
    recursive_iterations=2
    call get_command_argument(3,recursive_diagnostic)
    if(trim(recursive_diagnostic)=='first') recursive_iterations=1
    seed=initial
    allocate(seed_held%fields(count))
    do p=1,count
      allocate(seed_held%fields(p)%values(1,4,4,4),source=0._dp)
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        point=[geometry(p)%x_lower,geometry(p)%y_lower,geometry(p)%z_lower]+ &
          (real([x,y,z],dp)-.5_dp)*[geometry(p)%dx,geometry(p)%dy,geometry(p)%dz]
        pv=0;pv(1)=.31_dp;pv(2:4)=[2+5*point(1),-1._dp,.5_dp]
        do sp=1,size(species)
          if(trim(species(sp)%name)=='N2') pv(5+sp)=1
        enddo
        pv(5)=pv(1)*1000*(1+.01_dp*point(1)+.02_dp*point(3))* &
          sum(pv(6:)*8.31446261815324e3_dp/species%molecular_weight)
        call reactive_primitive_to_conserved(species,pv,state,t,sound,good)
        call require(good,'held prepared inert gradient EOS')
        seed%nodes(p)%state(:,x,y,z)=state;seed%nodes(p)%temperature(x,y,z)=t
        seed_held%fields(p)%values(1,x,y,z)=state(5)-.5_dp*sum(state(2:4)**2)/state(1)
        seed%nodes(p)%state(5,x,y,z)=state(5)-11
      enddo
      enddo
      enddo
    enddo
    policy%eta2=2;policy%allow_negative_energy=.false.
    if(test_mode=='held_prepared_parent') then
      call check_prepared_parent_replay(seed,seed_held,extended,faces,policy,dt,force)
      return
    endif
    if(test_mode=='held_prepared_serial') &
      call check_held_prepared_simultaneous(seed,seed_held,geometry,extended,faces,policy,dt,force)
    if(test_mode=='held_recursive_chemistry'.or.test_mode=='held_prepared_serial') then
      reference=seed;reference_held=seed_held;reference_heat=seed_held
      allocate(reference_ir%fields(count))
      do p=1,count
        allocate(reference_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
        reference_heat%fields(p)%values=0
      enddo
      reference_correction=7
      call integrate_reactive_eb_tree_3d(seed,v0,good)
      call require(good,'serial held reacting initial independent inventory')
      if(method=='none') then
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,recursive_iterations,.true., &
        [.true.,.true.,.true.,.true.],case_force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
        rb,rw,rch,rext,rtheta,good,reference_ir,mol_mode=.true.,redistribution=method, &
        target_volume_fraction=.9_dp,state_redist_max_order=0, &
        internal_energy=reference_held,internal_energy_flux=reference_scalar_faces,energy_policy=policy, &
        energy_correction=reference_correction,level_steps=reference_steps,integrated_flux=reference_faces, &
        registered_flux=reference_registered,heat_release=reference_heat,held_interval_history=serial_records)
      else
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,recursive_iterations,.true., &
        [.true.,.true.,.true.,.true.],case_force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
        rb,rw,rch,rext,rtheta,good,reference_ir,mol_mode=.true.,redistribution=method, &
        target_volume_fraction=.9_dp,state_redist_max_order=0, &
        internal_energy=reference_held,internal_energy_flux=reference_scalar_faces,energy_policy=policy, &
        energy_correction=reference_correction,level_steps=reference_steps,integrated_flux=reference_faces, &
        registered_flux=reference_registered,heat_release=reference_heat,held_interval_history=serial_records)
      endif
      call require(good,'serial held native recursive transport/chemistry accepts')
      call require(all(reference_steps==[1_int64,2_int64,4_int64]),'serial held native retains all physical clocks')
      call integrate_reactive_eb_tree_3d(reference,v1,good)
      call require(good,'serial held native accepted composite inventory')
      v1=v1-v0-rb-rw-rch-rext;v1(5)=v1(5)-(reference_correction-7)
      error=maxval(abs(v1)/max(1._dp,abs(v0)))
      if(rank==0) print *, 'Serial held recursive native inventory/correction: ',error,reference_correction
      call require(error<3.e-14_dp,'serial held native independent chemical/numeric inventory')
      kept_reference=reference;kept_held=reference_held;kept_ir=reference_ir;kept_heat=reference_heat
      kept_correction=reference_correction;kept_faces=reference_scalar_faces
      kept_serial_records=serial_records
      reactor_failure=.true.;reactor_failure_rank=rank;reactor_calls=0
      if(method=='none') then
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,recursive_iterations,.true., &
        [.true.,.true.,.true.,.true.],case_force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
        rb,rw,rch,rext,rtheta,good,reference_ir,mol_mode=.true.,redistribution=method, &
        target_volume_fraction=.9_dp,state_redist_max_order=0, &
        internal_energy=reference_held,internal_energy_flux=reference_scalar_faces,energy_policy=policy, &
        energy_correction=reference_correction,level_steps=reference_steps,integrated_flux=reference_faces, &
        registered_flux=reference_registered,heat_release=reference_heat,evaluate_reaction=reactor_adapter, &
        held_interval_history=serial_records)
      else
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,recursive_iterations,.true., &
        [.true.,.true.,.true.,.true.],case_force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
        rb,rw,rch,rext,rtheta,good,reference_ir,mol_mode=.true.,redistribution=method, &
        target_volume_fraction=.9_dp,state_redist_max_order=0, &
        internal_energy=reference_held,internal_energy_flux=reference_scalar_faces,energy_policy=policy, &
        energy_correction=reference_correction,level_steps=reference_steps,integrated_flux=reference_faces, &
        registered_flux=reference_registered,heat_release=reference_heat,evaluate_reaction=reactor_adapter, &
        held_interval_history=serial_records)
      endif
      reactor_failure=.false.
      call require(.not.good.and.reactor_calls==2.and.reference_correction==kept_correction, &
        'serial held recursive native refusal is not retried and retains accepted numeric prefix')
      do p=1,count
        call require(all(reference%nodes(p)%state==kept_reference%nodes(p)%state).and. &
          all(reference%nodes(p)%temperature==kept_reference%nodes(p)%temperature).and. &
          all(reference_held%fields(p)%values==kept_held%fields(p)%values).and. &
          all(reference_ir%fields(p)%values==kept_ir%fields(p)%values).and. &
          all(reference_heat%fields(p)%values==kept_heat%fields(p)%values), &
          'serial held recursive refusal rolls back gas/T/held/IR/heat together')
        call require(all(reference_scalar_faces(p)%x==kept_faces(p)%x).and. &
          all(reference_scalar_faces(p)%y==kept_faces(p)%y).and.all(reference_scalar_faces(p)%z==kept_faces(p)%z), &
          'serial held recursive refusal retains accepted scalar face family')
      enddo
    endif
    if(test_mode=='held_recursive_chemistry'.or.test_mode=='held_prepared_serial') then
      call require(size(serial_records)==7.and.all(serial_records%depth==kept_serial_records%depth).and. &
        all(serial_records%time==kept_serial_records%time),'serial held refusal retains all native clock metadata')
      do stage=1,size(serial_records)
        do p=1,count
          call require(all(serial_records(stage)%clock_internal_energy%fields(p)%values== &
            kept_serial_records(stage)%clock_internal_energy%fields(p)%values), &
            'serial held refusal retains the complete accepted scalar clock')
        enddo
      enddo
      if(method=='none') call check_serial_held_group(seed,seed_held,serial_records,faces,rw,rext)
    endif
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(gas%fields(count),temp%fields(count),held%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        gas%fields(p)%values=seed%nodes(p)%state
        allocate(temp%fields(p)%values(1,4,4,4));temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
        held%fields(p)=seed_held%fields(p)
        allocate(ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
      enddo
      if(test_mode=='held_prepared_recursive'.or. &
    test_mode=='held_recursive_chemistry'.or.test_mode=='held_prepared_serial') then
        heat=held
        do p=1,count
          if(owners(p)==rank) heat%fields(p)%values=0
        enddo
        call integrate_mpi_eb_tree_3d(topology,owners,geometry,gas,v0,MPI_COMM_WORLD,good)
        call require(good,'held prepared recursive initial independent composite inventory')
        correction=7
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
          dt,recursive_iterations,.true.,[.true.,.true.,.true.,.true.],case_force,2.e-7_dp,1.e-12_dp,extended,faces, &
          b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,mol_mode=.true.,redistribution=method, &
          target_volume_fraction=.9_dp,state_redist_max_order=0,internal_energy=held,internal_energy_flux=held_faces, &
          energy_policy=policy,energy_correction=correction,interval_history=recursive_records, &
          level_steps=steps,integrated_flux=flux,registered_flux=registered,heat_release=heat)
        call require(good,'native held prepared three-level recursive transport/chemistry')
        call require(all(recursive_records%accepted_coarse_stages.eqv.(method=='state_redist')), &
          'paired recursive ledgers explicitly identify accepted mean stages')
        call require(size(recursive_records)==7.and.all(steps==[1_int64,2_int64,4_int64]), &
          'held prepared recursive retains seven native clock intervals')
        call integrate_mpi_eb_tree_3d(topology,owners,geometry,gas,v1,MPI_COMM_WORLD,good)
        call require(good,'held prepared recursive accepted composite inventory')
        v1=v1-v0-b-w-ch-ext;v1(5)=v1(5)-(correction-7)
        error=maxval(abs(v1)/max(1._dp,abs(v0)))
        if(rank==0) print *, 'Held prepared recursive inventory/correction/replayed: ',layout,error,correction, &
          recursive_records%replayed
        if(rank==0) then
          print *, 'Held prepared recursive inventory residual: ',v1
          print *, 'Held prepared recursive initial inventory: ',v0
          print *, 'Held prepared recursive boundary/wall/chemical/external: ',b,w,ch,ext
        endif
        if(recursive_iterations>1.and.method=='state_redist') then
          zero_faces=flux
          do p=1,count
            if(owners(p)/=rank) cycle
            zero_faces(p)%x=0;zero_faces(p)%y=0;zero_faces(p)%z=0
          enddo
          allocate(discarded_registers(count));discarded_total=0
          do diagnostic_depth=0,1
            do p=1,count
              if(owners(p)/=rank) cycle
              call initialize_state_redist_tree_register_3d(topology,p,nv,6,discarded_registers(p),good,faces,scale=.5_dp)
              call require(good,'discarded prepared coarse native register geometry')
            enddo
            do stage=1,size(recursive_records)
              if(recursive_records(stage)%depth/=diagnostic_depth) cycle
              do p=1,count
                if(owners(p)/=rank.or.topology%nodes(p)%level/=diagnostic_depth) cycle
                call accumulate_state_redist_coarse_register_3d(discarded_registers(p), &
                  recursive_records(stage)%prepared_stages(1)%registers(p),nv,[16,16,16],good)
                call require(good,'isolate actual intermediate native coarse contribution')
              enddo
            enddo
            call assemble_state_redist_mpi_interface_rereflux_3d(topology,owners,extended,discarded_registers,nv, &
              diagnostic_depth,discarded,MPI_COMM_WORLD,good,faces,spill,has_spill)
            call require(good,'assemble isolated intermediate coarse debt with native masks')
            first_depth=diagnostic_depth
            if(has_spill) then
              first_depth=diagnostic_depth-1
              do p=1,count
                if(owners(p)/=rank) cycle
                if(topology%nodes(p)%level==first_depth) discarded%fields(p)=spill%fields(p)
              enddo
            endif
            call assemble_mpi_eb_tree_reflux_3d(species,topology,owners,geometry,gas,temp,zero_faces,discarded_increment, &
              discarded_faces,MPI_COMM_WORLD,good,diagnostic_depth,diagnostic_depth+1,faces, &
              rereflux_delta=discarded,rereflux_first_level=first_depth)
            call require(good,'route isolated intermediate coarse debt to exposed descendant leaves')
            call integrate_mpi_eb_tree_3d(topology,owners,geometry,discarded_increment,discarded_inventory,MPI_COMM_WORLD,good)
            call require(good,'integrate independent discarded coarse debt')
            discarded_total=discarded_total+discarded_inventory
          enddo
          if(rank==0) print *, 'Held prepared isolated intermediate coarse defect: ',discarded_total
          deallocate(discarded_registers)
        endif
        call require(error<3.e-14_dp,'held prepared recursive independent chemical/numeric composite inventory')
        local_norm=0
        do stage=1,size(recursive_records)
          call require(recursive_records(stage)%has_internal_energy.and. &
            allocated(recursive_records(stage)%clock_internal_energy%fields), &
            'prepared recursive archive retains complete held clock')
          if(recursive_records(stage)%depth>0.and.method=='state_redist') call require( &
            allocated(recursive_records(stage)%open_ancestor_held_registers), &
            'prepared recursive archive retains still-open parent scalar registers')
          do p=1,count
            call require(allocated(recursive_records(stage)%clock_internal_energy%fields(p)%values).eqv. &
              (owners(p)==rank),'prepared recursive full held clock is owner sparse')
            if(test_mode=='held_prepared_serial'.and.owners(p)==rank) then
              call require(all(serial_records(stage)%clock_initial%nodes(p)%state== &
                recursive_records(stage)%clock_initial%fields(p)%values).and. &
                all(serial_records(stage)%clock_internal_energy%fields(p)%values== &
                recursive_records(stage)%clock_internal_energy%fields(p)%values), &
                'serial/MPI prepared archives retain identical complete gas and scalar clock anchors')
              call require(serial_records(stage)%prepared_redistribution.and. &
                serial_records(stage)%accepted_coarse_stages.and. &
                size(serial_records(stage)%prepared_stages)==recursive_iterations+1, &
                'serial prepared physical clocks retain every actual operator and coarse-stage contract')
              call require(allocated(serial_records(stage)%ancestor_start%nodes).and. &
                allocated(serial_records(stage)%ancestor_end%nodes).and. &
                allocated(serial_records(stage)%ancestor_internal_energy_start%fields).and. &
                allocated(serial_records(stage)%ancestor_internal_energy_end%fields).and. &
                allocated(serial_records(stage)%fractions_start).and.allocated(serial_records(stage)%fractions_end), &
                'serial prepared clock keeps actual ancestor gas/scalar trajectories and interpolation fractions')
              call require(size(serial_records(stage)%fractions_start)==serial_records(stage)%depth.and. &
                size(serial_records(stage)%fractions_end)==serial_records(stage)%depth.and. &
                all(serial_records(stage)%fractions_end>serial_records(stage)%fractions_start), &
                'serial prepared descendant archives retain positive actual ancestor clock exposure')
              if(recursive_records(stage)%depth>0) then
                call require_register_equal(serial_records(stage)%open_ancestor_registers(p), &
                  recursive_records(stage)%open_ancestor_registers(p))
                call require_register_equal(serial_records(stage)%open_ancestor_held_registers(p), &
                  recursive_records(stage)%open_ancestor_held_registers(p))
              endif
            endif
            if(owners(p)/=rank.or.recursive_records(stage)%depth==0.or.method=='none') cycle
            if(allocated(recursive_records(stage)%open_ancestor_held_registers(p)%coarse_delta)) &
              local_norm=local_norm+sum(abs(recursive_records(stage)%open_ancestor_held_registers(p)%coarse_delta))
            if(allocated(recursive_records(stage)%open_ancestor_held_registers(p)%fine_delta)) &
              local_norm=local_norm+sum(abs(recursive_records(stage)%open_ancestor_held_registers(p)%fine_delta))
          enddo
        enddo
        call MPI_Allreduce(local_norm,global_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
        if(method=='state_redist') call require(global_norm>0, &
          'prepared recursive fixture records actual nonzero open scalar register debt')
        if(method=='state_redist'.and.test_mode/='held_prepared_serial') then
        if(layout==1) then
          reference=seed;reference_held=seed_held;reference_correction=correction
          do p=1,count
            if(owners(p)==rank) then
              reference%nodes(p)%state=gas%fields(p)%values
              reference%nodes(p)%temperature=temp%fields(p)%values(1,:,:,:)
              reference_held%fields(p)=held%fields(p)
            endif
            call MPI_Bcast(reference%nodes(p)%state,size(reference%nodes(p)%state),MPI_DOUBLE_PRECISION, &
              owners(p),MPI_COMM_WORLD,status)
            call MPI_Bcast(reference%nodes(p)%temperature,size(reference%nodes(p)%temperature),MPI_DOUBLE_PRECISION, &
              owners(p),MPI_COMM_WORLD,status)
            call MPI_Bcast(reference_held%fields(p)%values,size(reference_held%fields(p)%values),MPI_DOUBLE_PRECISION, &
              owners(p),MPI_COMM_WORLD,status)
          enddo
        else
          call require(abs(correction-reference_correction)<3.e-14_dp*max(1._dp,abs(reference_correction)), &
            'prepared recursive idle layout retains accepted numeric ledger')
          do p=1,count
            if(owners(p)/=rank) cycle
            call require(all(gas%fields(p)%values==reference%nodes(p)%state).and. &
              all(temp%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature).and. &
              all(held%fields(p)%values==reference_held%fields(p)%values), &
              'prepared recursive exact distributed/owner-zero gas/T/held endpoints')
          enddo
        endif
        else
          call require(abs(correction-reference_correction)<3.e-14_dp*max(1._dp,abs(reference_correction)), &
            'serial/MPI held reacting numeric ledger agrees across owner layouts')
          call require(all(steps==reference_steps),'serial/MPI held reacting clocks agree')
          do p=1,count
            if(owners(p)/=rank) cycle
            call require(all(gas%fields(p)%values==reference%nodes(p)%state).and. &
              all(temp%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature).and. &
              all(held%fields(p)%values==reference_held%fields(p)%values).and. &
              all(ir%fields(p)%values==reference_ir%fields(p)%values).and. &
              all(heat%fields(p)%values==reference_heat%fields(p)%values), &
              'serial/MPI native recursive gas/T/held/IR/heat endpoints agree exactly')
            call require(all(held_faces(p)%x==reference_scalar_faces(p)%x).and. &
              all(held_faces(p)%y==reference_scalar_faces(p)%y).and.all(held_faces(p)%z==reference_scalar_faces(p)%z), &
              'serial/MPI native recursive scalar faces agree exactly')
            call require(all(flux(p)%x==reference_faces(p)%x).and.all(flux(p)%y==reference_faces(p)%y).and. &
              all(flux(p)%z==reference_faces(p)%z),'serial/MPI native recursive physical gas faces agree exactly')
            call require(all(registered(p)%x==reference_registered(p)%x).and. &
              all(registered(p)%y==reference_registered(p)%y).and.all(registered(p)%z==reference_registered(p)%z), &
              'serial/MPI native recursive registered gas faces agree exactly')
          enddo
        endif
        if(test_mode=='held_prepared_serial') then
          call reconstruct_reactive_eb_mol_held_registration_3d(reference,serial_records(2), &
            serial_registered,serial_held_registered,good)
          call require(good,'serial native clock reconstructs the pinned admitted registration')
          call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,recursive_records(2), &
            recursive_records(2)%prepared_stages,nv,reconstructed,registered_regs,MPI_COMM_WORLD,good,faces, &
            internal_energy_rate=scalar_rate,internal_energy_registers=registered_held_regs, &
            register_provisional_coarse=.true.)
          call require(good,'MPI middle clock reconstructs chronological registered gas/held deltas')
          do p=1,count
            if(owners(p)/=rank.or.topology%nodes(p)%level/=1) cycle
            call require_register_equal(registered_regs(p),serial_registered(p))
            call require_register_equal(registered_held_regs(p),serial_held_registered(p))
          enddo
          reference_regs=serial_registered;reference_held_regs=serial_held_registered
          bad_serial=serial_records(2)
          do p=1,count
            if(topology%nodes(p)%level/=1) cycle
            if(.not.allocated(bad_serial%prepared_stages(1)%registers(p)%coarse_delta)) cycle
            malformed_dims=topology%nodes(p)%shape+12
            partial_registration=serial_registered(p)
            deallocate(partial_registration%coarse_delta)
            call accumulate_state_redist_mol_registration_stage_3d(partial_registration, &
              bad_serial%prepared_stages(1)%registers(p),1,2,nv,malformed_dims,good)
            call require(.not.good.and.allocated(partial_registration%coarse_flags).and. &
              .not.allocated(partial_registration%coarse_delta), &
              'partly allocated destination registration refuses without evaluating an absent delta bound')
            allocate(malformed_delta(0:nv-1,malformed_dims(1),malformed_dims(2),malformed_dims(3)))
            malformed_delta=bad_serial%prepared_stages(1)%registers(p)%coarse_delta
            call move_alloc(malformed_delta,bad_serial%prepared_stages(1)%registers(p)%coarse_delta)
            exit
          enddo
          call require(p<=count,'serial registration refusal probe reaches an actual provisional coarse ledger')
          call reconstruct_reactive_eb_mol_held_registration_3d(reference,bad_serial, &
            serial_registered,serial_held_registered,good)
          call require(.not.good,'late zero-based source registration refuses before indexing its delta')
          do p=1,count
            call require_register_equal(serial_registered(p),reference_regs(p))
            call require_register_equal(serial_held_registered(p),reference_held_regs(p))
          enddo
        endif
        kept_gas=gas;kept_t=temp;kept_held=held;kept_ir=ir;kept_heat=heat
        kept_records=recursive_records;kept_faces=held_faces
        kept_correction=correction;kept_chemical=ch
        reactor_failure=.true.;reactor_failure_rank=1;reactor_calls=0
        if(layout==2) reactor_failure_rank=0
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
          dt,recursive_iterations,.true.,[.true.,.true.,.true.,.true.],case_force,2.e-7_dp,1.e-12_dp,extended,faces, &
          b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,mol_mode=.true.,redistribution=method, &
          target_volume_fraction=.9_dp,state_redist_max_order=0,internal_energy=held,internal_energy_flux=held_faces, &
          energy_policy=policy,energy_correction=correction,interval_history=recursive_records, &
          level_steps=steps,integrated_flux=flux,registered_flux=registered,heat_release=heat,evaluate_reaction=reactor_adapter)
        reactor_failure=.false.
        call require(.not.good.and.correction==kept_correction.and.all(ch==kept_chemical), &
          'prepared recursive native provider refusal retains both ledgers')
        if(rank==reactor_failure_rank) call require(reactor_calls==2,'prepared recursive native refusal is not retried')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
            all(temp%fields(p)%values==kept_t%fields(p)%values).and. &
            all(held%fields(p)%values==kept_held%fields(p)%values).and. &
            all(ir%fields(p)%values==kept_ir%fields(p)%values).and. &
            all(heat%fields(p)%values==kept_heat%fields(p)%values), &
            'prepared recursive provider refusal discards all gas/T/held/IR/heat clock work')
          call require(all(held_faces(p)%x==kept_faces(p)%x).and.all(held_faces(p)%y==kept_faces(p)%y).and. &
            all(held_faces(p)%z==kept_faces(p)%z),'prepared recursive failure retains accepted scalar face family')
          do stage=1,size(recursive_records)
            call require(all(recursive_records(stage)%clock_internal_energy%fields(p)%values== &
              kept_records(stage)%clock_internal_energy%fields(p)%values), &
              'prepared recursive failure retains accepted complete held clock archives')
          enddo
        enddo
        deallocate(gas%fields,temp%fields,held%fields,ir%fields)
        cycle
      endif
      correction=7
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
        dt,2,.true.,[.true.,.true.,.true.,.true.],case_force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
        MPI_COMM_WORLD,good,ir,target_depth=1,mol_mode=.true.,redistribution='state_redist', &
        target_volume_fraction=.9_dp,state_redist_max_order=0,state_redist_registers=registers, &
        held_state_redist_registers=held_registers,integrated_flux=flux,internal_energy=held, &
        internal_energy_flux=held_faces,energy_policy=policy,energy_correction=correction,interval_history=record)
      call require(good,'native held prepared transport/chemistry interval')
      call require(record%accepted_coarse_stages,'held prepared archive records its accepted coarse-stage contract')
      call require(record%has_internal_energy.and.record%prepared_redistribution.and. &
        size(record%prepared_stages)==3,'held prepared chemistry archives every native stage')
      call require(correction==7+record%numerical_energy_correction, &
        'held prepared chemistry retains separate accepted numeric reset ledger')
      local_norm=0
      do p=1,count
        do stage=0,2
          call require(record%prepared_stages(stage)%has_internal_energy,'every prepared stage carries held family')
          call require(allocated(record%prepared_stages(stage)%held_raw_rate_halos%fields(p)%values).eqv. &
            (owners(p)==rank.and.topology%nodes(p)%level==1),'held prepared raw halos are sparse by owner/level')
        enddo
        if(owners(p)/=rank.or.topology%nodes(p)%level/=1) cycle
        if(allocated(record%prepared_stages(1)%held_registers(p)%coarse_delta)) local_norm=local_norm+ &
          sum(abs(record%prepared_stages(1)%held_registers(p)%coarse_delta))
      enddo
      call MPI_Allreduce(local_norm,global_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
      print *, 'Held prepared intermediate coarse ledger norm: ',layout,global_norm
      call require(global_norm>0,'fixture exercises actual intermediate scalar coarse-register accumulation')
      call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,record, &
        record%prepared_stages,nv,reconstructed,rebuilt,MPI_COMM_WORLD,good,faces, &
        internal_energy_rate=scalar_rate,internal_energy_registers=held_rebuilt)
      call require(good,'prepared held native stages reconstruct both source/register families')
      do p=1,count
        if(owners(p)/=rank.or.topology%nodes(p)%level/=1) cycle
        call require(all(reconstructed%fields(p)%values==record%forcing%fields(p)%values), &
          'prepared paired identity retains actual gas forcing exactly')
        call require(all(scalar_rate%fields(p)%values==record%internal_energy_spatial_rate%fields(p)%values), &
          'prepared paired identity retains independent held mean exactly')
        call require_register_equal(rebuilt(p),registers(p))
        call require_register_equal(held_rebuilt(p),held_registers(p))
        if(allocated(registers(p)%coarse_delta)) call require(all(registers(p)%coarse_delta== &
          record%prepared_stages(0)%registers(p)%coarse_delta+record%prepared_stages(2)%registers(p)%coarse_delta), &
          'accepted gas coarse ledger equals the actual old/final forcing family')
        if(allocated(held_registers(p)%coarse_delta)) call require(all(held_registers(p)%coarse_delta== &
          record%prepared_stages(0)%held_registers(p)%coarse_delta+record%prepared_stages(2)%held_registers(p)%coarse_delta), &
          'accepted scalar coarse ledger equals the actual old/final forcing family')
      enddo
      call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,record, &
        record%prepared_stages,nv,reconstructed,registered_regs,MPI_COMM_WORLD,good,faces, &
        internal_energy_rate=scalar_rate,internal_energy_registers=registered_held_regs,register_provisional_coarse=.true.)
      call require(good,'held prepared registration reconstructs the pinned all-evaluation coarse family')
      local_norm=0
      do p=1,count
        if(owners(p)/=rank.or.topology%nodes(p)%level/=1) cycle
        call require(all(reconstructed%fields(p)%values==record%forcing%fields(p)%values).and. &
          all(scalar_rate%fields(p)%values==record%internal_energy_spatial_rate%fields(p)%values), &
          'registration selection preserves gas forcing and independent held physical mean exactly')
        if(allocated(registered_regs(p)%coarse_delta)) call require(all(registered_regs(p)%coarse_delta== &
          (record%prepared_stages(0)%registers(p)%coarse_delta+record%prepared_stages(1)%registers(p)%coarse_delta)+ &
          record%prepared_stages(2)%registers(p)%coarse_delta),'pinned gas coarse registration includes provisional stage')
        if(allocated(registered_held_regs(p)%coarse_delta)) then
          call require(all(registered_held_regs(p)%coarse_delta== &
            (record%prepared_stages(0)%held_registers(p)%coarse_delta+ &
            record%prepared_stages(1)%held_registers(p)%coarse_delta)+ &
            record%prepared_stages(2)%held_registers(p)%coarse_delta), &
            'pinned held coarse registration includes provisional stage in evaluation order')
          local_norm=local_norm+sum(abs(registered_held_regs(p)%coarse_delta-held_registers(p)%coarse_delta))
        endif
        if(allocated(registered_regs(p)%fine_delta)) call require(all(registered_regs(p)%fine_delta== &
          registers(p)%fine_delta),'pinned fine gas registration excludes zero-flux provisional evaluations')
        if(allocated(registered_held_regs(p)%fine_delta)) call require(all(registered_held_regs(p)%fine_delta== &
          held_registers(p)%fine_delta),'pinned fine held registration excludes zero-flux provisional evaluations')
      enddo
      call MPI_Allreduce(local_norm,global_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
      call require(global_norm>0,'actual nonuniform native fixture distinguishes pinned and accepted scalar ledgers')
      kept_rate=reconstructed;kept_scalar_rate=scalar_rate
      replay_regs=registered_regs;replay_held_regs=registered_held_regs
      call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,record, &
        record%prepared_stages,nv,reconstructed,registered_regs,MPI_COMM_WORLD,good,faces, &
        internal_energy_rate=scalar_rate,internal_energy_registers=registered_held_regs, &
        register_provisional_coarse=rank/=nranks-1)
      call require(.not.good,'idle rank registration-policy disagreement refuses collectively')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(reconstructed%fields(p)%values==kept_rate%fields(p)%values).and. &
          all(scalar_rate%fields(p)%values==kept_scalar_rate%fields(p)%values), &
          'registration-policy disagreement preserves both reconstructed RHS outputs')
        call require_register_equal(registered_regs(p),replay_regs(p))
        call require_register_equal(registered_held_regs(p),replay_held_regs(p))
      enddo
      call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,record, &
        record%prepared_stages,nv,reconstructed,rebuilt,MPI_COMM_WORLD,good,faces)
      call require(.not.good,'prepared held archive refuses gas-only reconstruction')
      do trial=1,3
        weight=1
        if(trial==2) weight=.5_dp
        if(trial==3) weight=0
        do p=1,count
          if(owners(p)/=rank) cycle
          gas%fields(p)%values=seed%nodes(p)%state;temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
          held%fields(p)=seed_held%fields(p);ir%fields(p)%values=0
          if(topology%nodes(p)%level/=1) ir%fields(p)%values=17+trial
        enddo
        correction=7
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
          dt,2,.true.,[.true.,.true.,.true.,.true.],case_force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
          MPI_COMM_WORLD,good,ir,target_depth=1,mol_mode=.true.,redistribution='state_redist', &
          target_volume_fraction=.9_dp,state_redist_max_order=0,state_redist_registers=replay_regs, &
          held_state_redist_registers=replay_held_regs,integrated_flux=replay_faces,internal_energy=held, &
          internal_energy_flux=replay_scalar_faces,energy_policy=policy,energy_correction=correction, &
          interval_history=replayed,replay_interval=record,replay_spatial_weight=weight,provider_refused=refused)
        call require(good.and..not.refused.and.replayed%replayed,'native prepared held replacement replay accepts')
        if(rank==0) print *, 'Held prepared accepted replay layout/weight: ',layout,weight
        if(layout==1) then
          reference=seed;reference_held=seed_held;reference_heat=seed_held
          if(allocated(reference_ir%fields)) deallocate(reference_ir%fields)
          allocate(reference_ir%fields(count))
          do p=1,count
            allocate(reference_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
            reference_heat%fields(p)%values=0
            if(topology%nodes(p)%level/=1) reference_ir%fields(p)%values=17+trial
          enddo
          reference_correction=7
          call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,2,.true., &
            [.true.,.true.,.true.,.true.],case_force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
            rb,rw,rch,rext,rtheta,good,reference_ir,target_depth=1,mol_mode=.true.,redistribution='state_redist', &
            target_volume_fraction=.9_dp,state_redist_max_order=0,state_redist_registers=rebuilt, &
            held_state_redist_registers=held_rebuilt,integrated_flux=reference_faces,internal_energy=reference_held, &
            internal_energy_flux=reference_scalar_faces,energy_policy=policy,energy_correction=reference_correction, &
            heat_release=reference_heat,prepared_spatial_weight=weight,held_interval_history=serial_replayed, &
            replay_interval=serial_frozen,provider_refused=refused)
          call require(good.and..not.refused,'serial prepared native common transport/chemistry family accepts')
          reference_regs=rebuilt;reference_held_regs=held_rebuilt
          if(trial==1) allocate(serial_frozen,source=serial_replayed)
          call require(serial_replayed%prepared_redistribution.and.serial_replayed%accepted_coarse_stages.and. &
            size(serial_replayed%prepared_stages)==3,'serial prepared capture retains every actual native stage')
          call require(serial_replayed%replayed.eqv.(trial>1),'serial prepared replacement identifies actual replay')
          call require(all(serial_replayed%prepared_stages%spatial_weight==weight), &
            'serial prepared replacement uses one weight for all gas/scalar stages')
          call require(abs(reference_correction-correction)<3.e-14_dp*max(1._dp,abs(correction)), &
            'serial/MPI prepared weighted common numerical reset ledger agrees')
          do p=1,count
            if(topology%nodes(p)%level/=1) then
              call require(all(reference_ir%fields(p)%values==17+trial), &
                'serial prepared replay retains current inactive reaction histories instead of frozen tails')
              if(owners(p)==rank) call require(all(ir%fields(p)%values==17+trial), &
                'MPI prepared replay retains the same current inactive reaction histories')
              cycle
            endif
            if(owners(p)/=rank) cycle
            call require(all(reference%nodes(p)%state==gas%fields(p)%values).and. &
              all(reference%nodes(p)%temperature==temp%fields(p)%values(1,:,:,:)).and. &
              all(reference_held%fields(p)%values==held%fields(p)%values).and. &
              all(reference_ir%fields(p)%values==ir%fields(p)%values), &
              'serial/MPI prepared weights preserve exact gas/T/held/IR endpoints')
            do stage=0,2
              call require(all(serial_replayed%prepared_stages(stage)%external_rate%fields(p)%values== &
                replayed%prepared_stages(stage)%external_rate%fields(p)%values), &
                'serial/MPI weighted replacement retains identical actual force-work stages')
              call require(all(serial_replayed%prepared_stages(stage)%raw%fields(p)%values== &
                replayed%prepared_stages(stage)%raw_rate_halos%fields(p)%values).and. &
                all(serial_replayed%prepared_stages(stage)%scalar_raw%fields(p)%values== &
                replayed%prepared_stages(stage)%held_raw_rate_halos%fields(p)%values), &
                'serial/MPI prepared replacement rebuilds identical native gas and scalar raw halos')
            enddo
            call require_register_equal(rebuilt(p),replay_regs(p))
            call require_register_equal(held_rebuilt(p),replay_held_regs(p))
            call require(all(reference_faces(p)%x==replay_faces(p)%x).and. &
              all(reference_faces(p)%y==replay_faces(p)%y).and.all(reference_faces(p)%z==replay_faces(p)%z), &
              'serial/MPI prepared weighted gas faces agree exactly')
            call require(all(reference_scalar_faces(p)%x==replay_scalar_faces(p)%x).and. &
              all(reference_scalar_faces(p)%y==replay_scalar_faces(p)%y).and. &
              all(reference_scalar_faces(p)%z==replay_scalar_faces(p)%z), &
              'serial/MPI prepared weighted scalar faces agree exactly')
          enddo
        endif
        call require(all(replayed%prepared_stages%spatial_weight==weight), &
          'native prepared held replay uses one weight for every gas/scalar stage')
        call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,replayed, &
          replayed%prepared_stages,nv,reconstructed,rebuilt,MPI_COMM_WORLD,good,faces, &
          internal_energy_rate=scalar_rate,internal_energy_registers=held_rebuilt)
        call require(good,'accepted blended held stages reconstruct native paired operator')
        do p=1,count
          if(owners(p)/=rank.or.topology%nodes(p)%level/=1) cycle
          call require(all(reconstructed%fields(p)%values==replayed%forcing%fields(p)%values).and. &
            all(scalar_rate%fields(p)%values==replayed%internal_energy_spatial_rate%fields(p)%values), &
            'blended held archive identity preserves admitted gas and independent scalar means')
          call require_register_equal(rebuilt(p),replay_regs(p))
          call require_register_equal(held_rebuilt(p),replay_held_regs(p))
          if(weight==0) then
            if(allocated(replay_held_regs(p)%coarse_delta)) &
              call require(all(replay_held_regs(p)%coarse_delta==0),'zero family leaves no unfunded scalar coarse ledger')
            if(allocated(replay_held_regs(p)%fine_delta)) &
              call require(all(replay_held_regs(p)%fine_delta==0),'zero family leaves no unfunded scalar fine ledger')
          endif
        enddo
      enddo
      if(layout==1) call check_serial_prepared_replay_refusal(reference,reference_held,reference_ir,reference_heat, &
        serial_replayed,serial_frozen,reference_faces,reference_scalar_faces,reference_regs,reference_held_regs, &
        reference_correction,extended,faces,policy,dt,case_force)
      kept=replayed;kept_gas=gas;kept_t=temp;kept_held=held;kept_ir=ir
      kept_correction=correction;kept_chemical=ch;bad_record=record
      if(rank==nranks-1) bad_record%prepared_stages(2)%has_internal_energy=.false.
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
        dt,2,.true.,[.true.,.true.,.true.,.true.],case_force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
        MPI_COMM_WORLD,good,ir,target_depth=1,mol_mode=.true.,redistribution='state_redist', &
        target_volume_fraction=.9_dp,state_redist_max_order=0,state_redist_registers=replay_regs, &
        held_state_redist_registers=replay_held_regs,integrated_flux=replay_faces,internal_energy=held, &
        internal_energy_flux=replay_scalar_faces,energy_policy=policy,energy_correction=correction, &
        interval_history=replayed,replay_interval=bad_record,replay_spatial_weight=.5_dp)
      call require(.not.good.and.correction==kept_correction.and.all(ch==kept_chemical), &
        'prepared owner/idle stage mismatch retains accepted chemical/numeric ledger')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
          all(temp%fields(p)%values==kept_t%fields(p)%values).and. &
          all(held%fields(p)%values==kept_held%fields(p)%values).and. &
          all(ir%fields(p)%values==kept_ir%fields(p)%values),'prepared held replay mismatch rolls back all fields')
        if(topology%nodes(p)%level/=1) cycle
        call require(all(replayed%prepared_stages(2)%held_raw_rate_halos%fields(p)%values== &
          kept%prepared_stages(2)%held_raw_rate_halos%fields(p)%values), &
          'prepared held replay mismatch keeps the accepted scalar stage payload')
      enddo
      kept_rate=reconstructed;kept_scalar_rate=scalar_rate;bad_record=kept
      if(rank==nranks-1) bad_record%accepted_coarse_stages=.false.
      call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,bad_record, &
        bad_record%prepared_stages,nv,reconstructed,rebuilt,MPI_COMM_WORLD,good,faces, &
        internal_energy_rate=scalar_rate,internal_energy_registers=held_rebuilt)
      call require(.not.good,'paired reconstruction collectively refuses mixed coarse-stage contracts')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(reconstructed%fields(p)%values==kept_rate%fields(p)%values).and. &
          all(scalar_rate%fields(p)%values==kept_scalar_rate%fields(p)%values), &
          'coarse-stage contract refusal retains both admitted source families')
      enddo
      kept_replay_regs=replay_regs;kept_held_regs=replay_held_regs
      kept_replay_faces=replay_faces;kept_replay_scalar_faces=replay_scalar_faces
      kept_rebuilt=rebuilt;kept_held_rebuilt=held_rebuilt
      do funding_fault=1,7
        bad_record=kept
        if(rank==nranks-1) then
          select case(funding_fault)
          case(1);bad_record%has_source_funding=.true.
          case(2);allocate(bad_record%source_funding_stock%fields(0))
          case(3);allocate(bad_record%source_funding_impulse%fields(0))
          case(4);bad_record%native_iterations=1
          case(5);allocate(bad_record%auxiliary_stages(0:0))
          case(6);allocate(bad_record%provider_stages(2,0:0))
          case(7);allocate(bad_record%native_provider_geometry(0))
          end select
        endif
        call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,bad_record, &
          bad_record%prepared_stages,nv,reconstructed,rebuilt,MPI_COMM_WORLD,good,faces, &
          internal_energy_rate=scalar_rate,internal_energy_registers=held_rebuilt)
        call require(.not.good,'prepared reconstruction refuses idle-only funding flag or stale packet')
        reactor_calls=0
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
          dt,2,.true.,[.true.,.true.,.true.,.true.],case_force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
          MPI_COMM_WORLD,good,ir,target_depth=1,mol_mode=.true.,redistribution='state_redist', &
          target_volume_fraction=.9_dp,state_redist_max_order=0,state_redist_registers=replay_regs, &
          held_state_redist_registers=replay_held_regs,integrated_flux=replay_faces,internal_energy=held, &
          internal_energy_flux=replay_scalar_faces,energy_policy=policy,energy_correction=correction, &
          interval_history=replayed,replay_interval=bad_record,replay_spatial_weight=.5_dp, &
          evaluate_reaction=reactor_adapter)
        call require(.not.good.and.reactor_calls==0.and.correction==kept_correction.and.all(ch==kept_chemical), &
          'prepared stale funding metadata refuses before callbacks and ledger publication')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
            all(temp%fields(p)%values==kept_t%fields(p)%values).and. &
            all(held%fields(p)%values==kept_held%fields(p)%values).and. &
            all(ir%fields(p)%values==kept_ir%fields(p)%values), &
            'prepared stale funding refusal preserves all gas/T/H/IR fields')
          if(topology%nodes(p)%level/=1) cycle
          call require(all(reconstructed%fields(p)%values==kept_rate%fields(p)%values).and. &
            all(scalar_rate%fields(p)%values==kept_scalar_rate%fields(p)%values), &
            'prepared stale funding reconstruction preserves accepted gas/H forcing')
          call require_register_equal(rebuilt(p),kept_rebuilt(p))
          call require_register_equal(held_rebuilt(p),kept_held_rebuilt(p))
          call require_register_equal(replay_regs(p),kept_replay_regs(p))
          call require_register_equal(replay_held_regs(p),kept_held_regs(p))
          call require(all(replay_faces(p)%x==kept_replay_faces(p)%x).and. &
            all(replay_faces(p)%y==kept_replay_faces(p)%y).and.all(replay_faces(p)%z==kept_replay_faces(p)%z).and. &
            all(replay_scalar_faces(p)%x==kept_replay_scalar_faces(p)%x).and. &
            all(replay_scalar_faces(p)%y==kept_replay_scalar_faces(p)%y).and. &
            all(replay_scalar_faces(p)%z==kept_replay_scalar_faces(p)%z), &
            'prepared stale funding refusal preserves complete paired faces')
          call require(all(replayed%prepared_stages(2)%held_raw_rate_halos%fields(p)%values== &
            kept%prepared_stages(2)%held_raw_rate_halos%fields(p)%values), &
            'prepared stale funding refusal preserves accepted scalar stage archive')
        enddo
      enddo
      bad_record=kept
      if(owners(2)==rank) bad_record%prepared_stages(2)%held_raw_rate_halos%fields(2)%values(1,7,7,7)= &
        ieee_value(0._dp,ieee_quiet_nan)
      call reconstruct_mpi_eb_mol_state_redist_interval_3d(topology,owners,geometry,extended,bad_record, &
        bad_record%prepared_stages,nv,reconstructed,rebuilt,MPI_COMM_WORLD,good,faces, &
        internal_energy_rate=scalar_rate,internal_energy_registers=held_rebuilt)
      call require(.not.good,'prepared paired reconstruction collectively refuses remote scalar NaN')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(reconstructed%fields(p)%values==kept_rate%fields(p)%values).and. &
          all(scalar_rate%fields(p)%values==kept_scalar_rate%fields(p)%values), &
          'failed paired reconstruction preserves both accepted owner-local source families')
      enddo
      kept_held_regs=replay_held_regs
      reactor_calls=0;reactor_failure=.true.;reactor_failure_rank=1
      if(layout==2) reactor_failure_rank=0
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
        dt,2,.true.,[.true.,.true.,.true.,.true.],case_force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
        MPI_COMM_WORLD,good,ir,target_depth=1,mol_mode=.true.,redistribution='state_redist', &
        target_volume_fraction=.9_dp,state_redist_max_order=0,state_redist_registers=replay_regs, &
        held_state_redist_registers=replay_held_regs,integrated_flux=replay_faces,internal_energy=held, &
        internal_energy_flux=replay_scalar_faces,energy_policy=policy,energy_correction=correction, &
        interval_history=replayed,replay_interval=record,replay_spatial_weight=.5_dp, &
        evaluate_reaction=reactor_adapter,provider_refused=refused)
      reactor_failure=.false.
      call require(.not.good.and.refused.and.correction==kept_correction.and.all(ch==kept_chemical), &
        'prepared held provider refusal retains chemical/numeric ledger and refuses all retries')
      if(rank==reactor_failure_rank) call require(reactor_calls==2,'prepared held provider refusal is not retried')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
          all(temp%fields(p)%values==kept_t%fields(p)%values).and. &
          all(held%fields(p)%values==kept_held%fields(p)%values).and. &
          all(ir%fields(p)%values==kept_ir%fields(p)%values),'prepared provider refusal rolls back all coupled fields')
        if(topology%nodes(p)%level/=1) cycle
        call require_register_equal(replay_held_regs(p),kept_held_regs(p))
        call require(all(replayed%prepared_stages(2)%held_raw_rate_halos%fields(p)%values== &
          kept%prepared_stages(2)%held_raw_rate_halos%fields(p)%values), &
          'prepared provider refusal retains accepted raw scalar stages and registers')
      enddo
      deallocate(gas%fields,temp%fields,held%fields,ir%fields)
    enddo
    if(rank==0) print *, 'PASS: prepared held chemistry/stages, scalar registers, native intervals and selected mode'

  end subroutine check_held_prepared

  subroutine check_prepared_parent_replay(seed,seed_held,extended,faces,policy,dt,force)
    type(eb_reactive_tree_3d), intent(in) :: seed
    type(amr_tree_state_3d), intent(in) :: seed_held
    type(eb_geometry_3d), intent(in) :: extended(:)
    character(len=*), intent(in) :: faces(6)
    type(dual_energy_policy), intent(in) :: policy
    real(dp), intent(in) :: dt,force(3)
    type(eb_reactive_tree_3d) :: reference,work,kept,next_clock
    type(amr_tree_state_3d) :: reference_held,reference_ir,reference_heat,held,ir,heat,kept_held,kept_ir,kept_heat
    type(eb_mol_held_interval_3d), allocatable :: native(:),records(:),bad(:),kept_records(:)
    type(eb_integrated_flux_3d), allocatable :: reference_flux(:),reference_scalar(:),reference_registered(:)
    type(eb_integrated_flux_3d), allocatable :: flux(:),scalar(:),registered(:),kept_flux(:),kept_scalar(:)
    type(state_redist_register_3d), allocatable :: parent(:),parent_scalar(:),kept_parent(:),kept_parent_scalar(:)
    type(amr_tree_state_3d) :: next_held
    real(dp) :: reference_prefix,prefix,kept_prefix,b(nv),w(nv),ch(nv),ext(nv),theta,v0(nv),v1(nv),error,weight,failure_dt
    integer(int64) :: steps(3)
    integer :: p,r,trial,failure,x,y,z
    real(dp) :: next_prefix,cell_change,pv(reactive_nprim(size(species))),temp,sound
    logical :: good,refused
    reference=seed;reference_held=seed_held;reference_heat=seed_held
    allocate(reference_ir%fields(count))
    do p=1,count
      allocate(reference_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
      reference_heat%fields(p)%values=0
    enddo
    reference_prefix=7
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,2,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces,b,w,ch,ext,theta,good, &
      reference_ir,mol_mode=.true.,redistribution='state_redist',target_volume_fraction=.9_dp, &
      state_redist_max_order=0,internal_energy=reference_held,internal_energy_flux=reference_scalar, &
      energy_policy=policy,energy_correction=reference_prefix,heat_release=reference_heat, &
      held_interval_history=native,integrated_flux=reference_flux,registered_flux=reference_registered,level_steps=steps)
    call require(good,'prepared serial parent fixture admits its original native interval')
    call require(size(native)==7.and.all(steps==[1_int64,2_int64,4_int64]), &
      'prepared serial parent fixture captures all actual physical intervals')
    call integrate_reactive_eb_tree_3d(seed,v0,good)
    call require(good,'prepared whole-clock independent initial inventory')
    do trial=1,3
      weight=1;if(trial==2) weight=.5_dp;if(trial==3) weight=0
      work=reference;held=reference_held;ir=reference_ir;heat=reference_heat;records=native
      do p=1,count
        work%nodes(p)%state(5,:,:,:)=work%nodes(p)%state(5,:,:,:)+13
        held%fields(p)%values=held%fields(p)%values+17;ir%fields(p)%values=19;heat%fields(p)%values=23
      enddo
      prefix=29;b=17;w=19;ch=23;ext=29;theta=.7_dp
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,2,.true., &
        [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,work,extended,faces,b,w,ch,ext,theta,good,ir, &
        mol_mode=.true.,redistribution='state_redist',target_volume_fraction=.9_dp,state_redist_max_order=0, &
        internal_energy=held,internal_energy_flux=scalar,energy_policy=policy,energy_correction=prefix,heat_release=heat, &
        held_interval_history=records,replay_interval_history=native,replay_spatial_weight=weight, &
        integrated_flux=flux,registered_flux=registered,level_steps=steps,provider_refused=refused)
      call require(good.and..not.refused.and.all(records%replayed).and.all(steps==[1_int64,2_int64,4_int64]), &
        'prepared native complete-window replay restores every actual clock without changing physical counts')
      call require(all(records%time==native%time).and.all(records%interval==native%interval), &
        'prepared complete-window replay retains each actual interval start and duration')
      do r=1,size(records)
        call require(all(records(r)%prepared_stages%spatial_weight==weight), &
          'prepared complete-window gas/scalar operators use one weight at every actual stage')
      enddo
      call integrate_reactive_eb_tree_3d(work,v1,good)
      call require(good,'prepared complete-window accepted independent inventory')
      v1=v1-v0-b-w-ch-ext;v1(5)=v1(5)-(prefix-7)
      error=maxval(abs(v1)/max(1._dp,abs(v0)))
      if(rank==0) print *, 'Serial prepared parent replay weight/inventory/prefix: ',weight,error,prefix
      call require(error<3.e-14_dp,'prepared complete-window separates physical, chemical and numerical inventory')
      if(trial/=1) cycle
      do p=1,count
        call require(all(work%nodes(p)%state==reference%nodes(p)%state).and. &
          all(work%nodes(p)%temperature==reference%nodes(p)%temperature).and. &
          all(held%fields(p)%values==reference_held%fields(p)%values).and. &
          all(ir%fields(p)%values==reference_ir%fields(p)%values).and. &
          all(heat%fields(p)%values==reference_heat%fields(p)%values), &
          'prepared unit-weight complete-window gas/T/UEINT/IR/heat equals the original native run exactly')
        call require(all(flux(p)%x==reference_flux(p)%x).and.all(flux(p)%y==reference_flux(p)%y).and. &
          all(flux(p)%z==reference_flux(p)%z).and.all(scalar(p)%x==reference_scalar(p)%x).and. &
          all(scalar(p)%y==reference_scalar(p)%y).and.all(scalar(p)%z==reference_scalar(p)%z), &
          'prepared unit-weight complete-window physical gas/scalar faces equal native exactly')
        call require(all(registered(p)%x==reference_registered(p)%x).and. &
          all(registered(p)%y==reference_registered(p)%y).and.all(registered(p)%z==reference_registered(p)%z), &
          'prepared unit-weight complete-window registered gas faces equal native exactly')
      enddo
    enddo
    ! Replay only the first completed level-one window. Preserve the still-open
    ! level-zero trajectory and its original contribution, including siblings.
    prefix=29;records=native(2:4)
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,native(2)%interval,2,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,work,extended,faces,b,w,ch,ext,theta,good,ir, &
      mol_mode=.true.,redistribution='state_redist',target_volume_fraction=.9_dp,state_redist_max_order=0, &
      internal_energy=held,internal_energy_flux=scalar,energy_policy=policy,energy_correction=prefix,heat_release=heat, &
      held_interval_history=records,replay_interval_history=native(2:4),replay_spatial_weight=1._dp, &
      replay_parent_registers=parent,replay_parent_held_registers=parent_scalar, &
      integrated_flux=flux,registered_flux=registered,level_steps=steps,provider_refused=refused)
    call require(good.and..not.refused.and.all(records%replayed).and.all(steps==[0_int64,1_int64,2_int64]), &
      'prepared sub-parent replay retains its own complete physical window and enclosing clock')
    ! The next sibling archive follows its canonical level-one reset, which
    ! repairs covered cells too, without charging their composite inventory.
    next_clock=work;next_held=held;next_prefix=prefix
    do p=1,count
      if(topology%nodes(p)%level/=1) cycle
      associate(g=>next_clock%nodes(p)%geometry)
      do z=1,g%nz
      do y=1,g%ny
      do x=1,g%nx
        if(g%volume_fraction(x,y,z)<=0) cycle
        cell_change=0
        call reset_nasa7_dual_energy_cell(species,next_clock%nodes(p)%state(:,x,y,z), &
          next_held%fields(p)%values(1,x,y,z),next_clock%nodes(p)%temperature(x,y,z),policy,cell_change,good)
        call require(good,'independent next sibling canonical cell reset admits')
        if(.not.covered(topology,p,[x,y,z])) next_prefix=next_prefix+g%dx*g%dy*g%dz*g%volume_fraction(x,y,z)*cell_change
        call recover_reactive_eb_cell_3d(species,g,[x,y,z],next_clock%nodes(p)%state(:,x,y,z), &
          next_clock%nodes(p)%temperature(x,y,z),pv,temp,sound,good, &
          internal_energy_density=next_held%fields(p)%values(1,x,y,z))
        call require(good,'independent next sibling canonical EOS recovery admits')
        next_clock%nodes(p)%temperature(x,y,z)=temp
      enddo
      enddo
      enddo
      end associate
    enddo
    call require(abs(next_prefix-native(5)%clock_energy_correction)<3.e-14_dp*max(1._dp,abs(next_prefix)), &
      'prepared sub-parent replay and next canonical reset restore the actual numerical prefix')
    do p=1,count
      if(rank==0) print *, 'Prepared sub-parent pre/post-next-reset state/T difference: ',p, &
        maxval(abs(work%nodes(p)%state-native(5)%clock_initial%nodes(p)%state)), &
        maxval(abs(next_clock%nodes(p)%state-native(5)%clock_initial%nodes(p)%state)), &
        maxval(abs(work%nodes(p)%temperature-native(5)%clock_initial%nodes(p)%temperature)), &
        maxval(abs(next_clock%nodes(p)%temperature-native(5)%clock_initial%nodes(p)%temperature))
      call require(all(next_clock%nodes(p)%state==native(5)%clock_initial%nodes(p)%state).and. &
        all(next_clock%nodes(p)%temperature==native(5)%clock_initial%nodes(p)%temperature).and. &
        all(next_held%fields(p)%values==native(5)%clock_internal_energy%fields(p)%values).and. &
        all(ir%fields(p)%values==native(5)%clock_reaction_source%fields(p)%values).and. &
        all(heat%fields(p)%values==native(5)%clock_heat%fields(p)%values), &
        'prepared sub-parent replay with independent canonical reset exactly reaches the next native sibling clock')
      call require_register_equal(parent(p),native(5)%open_ancestor_registers(p))
      call require_register_equal(parent_scalar(p),native(5)%open_ancestor_held_registers(p))
    enddo
    kept=work;kept_held=held;kept_ir=ir;kept_heat=heat;kept_records=records;kept_prefix=prefix
    kept_flux=flux;kept_scalar=scalar;kept_parent=parent;kept_parent_scalar=parent_scalar
    do failure=1,3
      bad=native;failure_dt=dt
      if(failure==1) bad(7)%time=bad(7)%time+bad(7)%interval/10
      if(failure==2) then
        bad=native(2:4);failure_dt=bad(1)%interval
        bad(1)%open_ancestor_held_registers(1)%scale=.25_dp
      endif
      reactor_calls=0;reactor_failure_rank=rank;reactor_replay_refused=.false.;reactor_replay_refusal_actual=0
      prepared_parent_refusal_call=0
      if(failure==3) prepared_parent_refusal_call=2*sum(merge(1,0,seed%nodes(1)%geometry%volume_fraction>0))+2
      b=17;w=19;ch=23;ext=29;theta=.7_dp
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,failure_dt,2,.true., &
        [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,work,extended,faces,b,w,ch,ext,theta,good,ir, &
        mol_mode=.true.,redistribution='state_redist',target_volume_fraction=.9_dp,state_redist_max_order=0, &
        internal_energy=held,internal_energy_flux=scalar,energy_policy=policy,energy_correction=prefix,heat_release=heat, &
        held_interval_history=records,replay_interval_history=bad,replay_spatial_weight=.5_dp, &
        replay_parent_registers=parent,replay_parent_held_registers=parent_scalar, &
        integrated_flux=flux,evaluate_reaction=reactor_adapter,provider_refused=refused)
      call require(.not.good.and.prefix==kept_prefix.and.all(b==17).and.all(w==19).and.all(ch==23).and. &
        all(ext==29).and.theta==.7_dp,'prepared failed whole-window replay retains every accepted ledger')
      if(failure<3) then
        call require(.not.refused.and.reactor_calls==0,'prepared incomplete clock preflight precedes every reactor call')
      else
        call require(refused.and.reactor_replay_refused.and.reactor_calls==prepared_parent_refusal_call.and. &
          reactor_replay_refusal_actual==prepared_parent_refusal_call, &
          'prepared native-admitted child refusal after completed root physics is not retried')
      endif
      prepared_parent_refusal_call=0
      do p=1,count
        call require(all(work%nodes(p)%state==kept%nodes(p)%state).and. &
          all(work%nodes(p)%temperature==kept%nodes(p)%temperature).and. &
          all(held%fields(p)%values==kept_held%fields(p)%values).and.all(ir%fields(p)%values==kept_ir%fields(p)%values).and. &
          all(heat%fields(p)%values==kept_heat%fields(p)%values), &
          'prepared failed whole-window replay preserves gas/T/UEINT/IR/heat together')
        call require(all(flux(p)%x==kept_flux(p)%x).and.all(flux(p)%y==kept_flux(p)%y).and. &
          all(flux(p)%z==kept_flux(p)%z).and.all(scalar(p)%x==kept_scalar(p)%x).and. &
          all(scalar(p)%y==kept_scalar(p)%y).and.all(scalar(p)%z==kept_scalar(p)%z), &
          'prepared failed whole-window replay retains both accepted face families')
        call require_register_equal(parent(p),kept_parent(p))
        call require_register_equal(parent_scalar(p),kept_parent_scalar(p))
        do r=1,size(records)
          call require(all(records(r)%clock_initial%nodes(p)%state==kept_records(r)%clock_initial%nodes(p)%state).and. &
            all(records(r)%clock_internal_energy%fields(p)%values==kept_records(r)%clock_internal_energy%fields(p)%values), &
            'prepared failed whole-window replay retains all previously accepted clock records')
        enddo
      enddo
    enddo
  end subroutine check_prepared_parent_replay

  subroutine check_serial_prepared_replay_refusal(accepted,accepted_held,accepted_ir,accepted_heat, &
      accepted_record,frozen,accepted_faces,accepted_scalar,accepted_regs,accepted_held_regs,prefix, &
      extended,faces,policy,dt,force)
    type(eb_reactive_tree_3d), intent(in) :: accepted
    type(amr_tree_state_3d), intent(in) :: accepted_held,accepted_ir,accepted_heat
    type(eb_mol_held_interval_3d), intent(in) :: accepted_record,frozen
    type(eb_integrated_flux_3d), intent(in) :: accepted_faces(:),accepted_scalar(:)
    type(state_redist_register_3d), intent(in) :: accepted_regs(:),accepted_held_regs(:)
    type(eb_geometry_3d), intent(in) :: extended(:)
    character(len=*), intent(in) :: faces(6)
    type(dual_energy_policy), intent(in) :: policy
    real(dp), intent(in) :: prefix,dt,force(3)
    type(eb_reactive_tree_3d) :: work
    type(amr_tree_state_3d) :: held,ir,heat
    type(eb_mol_held_interval_3d) :: record,bad
    type(eb_integrated_flux_3d), allocatable :: flux(:),scalar(:)
    type(state_redist_register_3d), allocatable :: registers(:),scalar_registers(:)
    real(dp) :: correction,b(nv),w(nv),ch(nv),ext(nv),theta
    integer :: failure,p,stage
    logical :: good,refused
    do failure=1,2
      if(rank==0) print *, 'Serial prepared refusal case entry: ',failure,MPI_Wtime()
      work=accepted;held=accepted_held;ir=accepted_ir;heat=accepted_heat;record=accepted_record;bad=frozen
      flux=accepted_faces;scalar=accepted_scalar;registers=accepted_regs;scalar_registers=accepted_held_regs
      correction=prefix;b=17;w=19;ch=23;ext=29;theta=.7_dp
      if(failure==1) bad%prepared_stages(2)%scalar_raw%fields(2)%values(1,7,7,7)= &
        ieee_value(0._dp,ieee_quiet_nan)
      reactor_calls=0;reactor_failure=failure==2;reactor_failure_rank=rank
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,2,.true., &
        [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,work,extended,faces, &
        b,w,ch,ext,theta,good,ir,target_depth=1,mol_mode=.true.,redistribution='state_redist', &
        target_volume_fraction=.9_dp,state_redist_max_order=0,state_redist_registers=registers, &
        held_state_redist_registers=scalar_registers,integrated_flux=flux,internal_energy=held, &
        internal_energy_flux=scalar,energy_policy=policy,energy_correction=correction,heat_release=heat, &
        prepared_spatial_weight=.5_dp,held_interval_history=record,replay_interval=bad, &
        evaluate_reaction=reactor_adapter,provider_refused=refused)
      reactor_failure=.false.
      if(rank==0) print *, 'Serial prepared refusal case result/calls: ',failure,good,refused,reactor_calls,MPI_Wtime()
      call require(.not.good.and.correction==prefix.and.all(b==17).and.all(w==19).and.all(ch==23).and. &
        all(ext==29).and.theta==.7_dp,'serial prepared replay refusal retains every accepted diagnostic')
      if(failure==1) then
        call require(.not.refused.and.reactor_calls==0,'serial malformed scalar stage refuses before reaction')
      else
        call require(refused.and.reactor_calls==2,'serial prepared native-admitted refusal is not retried')
      endif
      do p=1,count
        call require(all(work%nodes(p)%state==accepted%nodes(p)%state).and. &
          all(work%nodes(p)%temperature==accepted%nodes(p)%temperature).and. &
          all(held%fields(p)%values==accepted_held%fields(p)%values).and. &
          all(ir%fields(p)%values==accepted_ir%fields(p)%values).and. &
          all(heat%fields(p)%values==accepted_heat%fields(p)%values), &
          'serial prepared failed replay retains gas/T/UEINT/IR/heat')
        ! Selected-level face and register payloads are intentionally sparse.
        ! Inactive nodes have no face arrays to compare.
        if(topology%nodes(p)%level/=1) cycle
        call require(all(flux(p)%x==accepted_faces(p)%x).and.all(flux(p)%y==accepted_faces(p)%y).and. &
          all(flux(p)%z==accepted_faces(p)%z).and.all(scalar(p)%x==accepted_scalar(p)%x).and. &
          all(scalar(p)%y==accepted_scalar(p)%y).and.all(scalar(p)%z==accepted_scalar(p)%z), &
          'serial prepared failed replay retains paired accepted faces')
        call require_register_equal(registers(p),accepted_regs(p))
        call require_register_equal(scalar_registers(p),accepted_held_regs(p))
        do stage=0,2
          call require(all(record%prepared_stages(stage)%scalar_raw%fields(p)%values== &
            accepted_record%prepared_stages(stage)%scalar_raw%fields(p)%values).and. &
            all(record%prepared_stages(stage)%raw%fields(p)%values== &
            accepted_record%prepared_stages(stage)%raw%fields(p)%values), &
            'serial prepared failed replay retains every accepted native stage payload')
        enddo
      enddo
    enddo
  end subroutine check_serial_prepared_replay_refusal

  subroutine check_serial_held_group(seed,seed_held,native_records,faces,wall,external)
    type(eb_reactive_tree_3d), intent(in) :: seed
    type(amr_tree_state_3d), intent(in) :: seed_held
    type(eb_mol_held_interval_3d), intent(in) :: native_records(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: wall(nv),external(nv)
    type(eb_reactive_tree_3d) :: work,kept
    type(amr_tree_state_3d) :: held,ir,heat,kept_held,kept_ir,kept_heat
    type(amr_tree_state_3d) :: expected_heat,expected_ir,expected_held,reaction_u,reaction_t
    type(eb_mol_held_interval_3d), allocatable :: records(:),kept_records(:),bad(:)
    type(eb_integrated_flux_3d), allocatable :: flux(:),scalar(:),kept_flux(:),kept_scalar(:)
    real(dp) :: weight,prefix,kept_prefix,chemical(nv),kept_chemical(nv),v0(nv),v1(nv),boundary(nv),error
    integer :: trial,p,c,n(3),r
    logical :: good,refused
    call integrate_reactive_eb_tree_3d(seed,v0,good)
    call require(good,'serial paired group independent initial inventory')
    do trial=1,3
      weight=1;if(trial==2) weight=.5_dp;if(trial==3) weight=0
      work=seed;held=seed_held;heat=seed_held;ir=native_records(1)%clock_reaction_source;records=native_records
      ! A completed descendant tail must be discarded in every field and in
      ! the numerical prefix, rather than becoming a replay initial condition.
      do p=1,count
        work%nodes(p)%state(5,:,:,:)=work%nodes(p)%state(5,:,:,:)+13
        held%fields(p)%values=held%fields(p)%values+17
        ir%fields(p)%values=19;heat%fields(p)%values=23
      enddo
      prefix=29;chemical=31
      call replay_reactive_eb_mol_held_group_3d(species,reactions,work,held,ir,prefix,records,weight, &
        flux,scalar,chemical,2.e-7_dp,1.e-12_dp,good,faces,heat_release=heat,provider_refused=refused)
      call require(good.and..not.refused.and.all(records%replayed),'serial native paired whole-clock group accepts')
      call require(all(records%time==native_records%time).and.all(records%interval==native_records%interval).and. &
        sum(merge(1,0,records%depth==0))==1.and.sum(merge(1,0,records%depth==1))==2.and. &
        sum(merge(1,0,records%depth==2))==4, &
        'serial paired group retains every physical clock duration and ordering')
      boundary=0;n=work%topology%nodes(1)%shape
      associate(g=>work%nodes(1)%geometry,f=>flux(1))
      do c=1,nv
        boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
          g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
          g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)-g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
          g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)-g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
      enddo
      end associate
      call integrate_reactive_eb_tree_3d(work,v1,good)
      call require(good,'serial paired group independent accepted inventory')
      v1=v1-v0-boundary-wall-chemical-external;v1(5)=v1(5)-(prefix-7)
      error=maxval(abs(v1)/max(1._dp,abs(v0)))
      if(rank==0) print *, 'Serial paired full-clock weight/inventory/correction: ',weight,error,prefix
      call require(error<3.e-14_dp,'serial paired group retains independent mass/species/momentum/energy ledger')
    enddo
    ! Compare heat against actual native solves on the accepted forcing and
    ! old UEINT, rather than assume the inert diagnostic rounds to exact zero.
    expected_heat=native_records(1)%clock_heat
    do r=1,size(records)
      reaction_u=records(r)%forcing;reaction_t=records(r)%initial_internal_energy
      expected_ir=records(r)%initial_reaction_source;expected_held=records(r)%initial_internal_energy
      do p=1,count
        reaction_u%fields(p)%values=records(r)%initial%nodes(p)%state
        reaction_t%fields(p)%values(1,:,:,:)=records(r)%initial%nodes(p)%temperature
      enddo
      call react_sdc_eb_tree_3d(species,reactions,records(r)%initial,records(r)%forcing,records(r)%interval, &
        2.e-7_dp,1.e-12_dp,reaction_u,reaction_t,expected_ir,good,target_depth=records(r)%depth, &
        heat_release=expected_heat,internal_energy=expected_held,old_internal_energy=records(r)%initial_internal_energy)
      call require(good,'accepted serial paired group has independent native heat reference')
    enddo
    call average_down_reactive_eb_tree_field_3d(work,expected_heat,good)
    call require(good,'independent native heat reference restricts the same physical hierarchy')
    do p=1,count
      if(rank==0) print *, 'Serial paired replay/native heat maximum: ',p,maxval(abs(heat%fields(p)%values)), &
        maxval(abs(expected_heat%fields(p)%values))
      call require(all(heat%fields(p)%values==expected_heat%fields(p)%values), &
        'serial paired replay discards stale heat and equals independently evaluated native heat exactly')
    enddo
    kept=work;kept_held=held;kept_ir=ir;kept_heat=heat;kept_records=records
    kept_flux=flux;kept_scalar=scalar;kept_prefix=prefix;kept_chemical=chemical
    reactor_calls=0;reactor_failure=.true.;reactor_failure_rank=rank
    call replay_reactive_eb_mol_held_group_3d(species,reactions,work,held,ir,prefix,records,.5_dp, &
      flux,scalar,chemical,2.e-7_dp,1.e-12_dp,good,faces,evaluate_reaction=reactor_adapter, &
      heat_release=heat,provider_refused=refused)
    reactor_failure=.false.
    call require(.not.good.and.refused.and.reactor_calls==2.and.prefix==kept_prefix.and.all(chemical==kept_chemical), &
      'serial paired group native-admitted provider refusal is not retried and retains both ledgers')
    bad=records;bad(size(bad))%time=bad(size(bad))%time+bad(size(bad))%interval/10
    reactor_calls=0
    call replay_reactive_eb_mol_held_group_3d(species,reactions,work,held,ir,prefix,bad,0._dp, &
      flux,scalar,chemical,2.e-7_dp,1.e-12_dp,good,faces,evaluate_reaction=reactor_adapter,heat_release=heat)
    call require(.not.good.and.reactor_calls==0.and.prefix==kept_prefix.and.all(chemical==kept_chemical), &
      'serial paired group refuses an incomplete physical clock before native reaction')
    do p=1,count
      call require(all(work%nodes(p)%state==kept%nodes(p)%state).and. &
        all(work%nodes(p)%temperature==kept%nodes(p)%temperature).and. &
        all(held%fields(p)%values==kept_held%fields(p)%values).and.all(ir%fields(p)%values==kept_ir%fields(p)%values).and. &
        all(heat%fields(p)%values==kept_heat%fields(p)%values),'serial paired group rejection retains every clock field')
      call require(all(flux(p)%x==kept_flux(p)%x).and.all(flux(p)%y==kept_flux(p)%y).and. &
        all(flux(p)%z==kept_flux(p)%z).and.all(scalar(p)%x==kept_scalar(p)%x).and. &
        all(scalar(p)%y==kept_scalar(p)%y).and.all(scalar(p)%z==kept_scalar(p)%z), &
        'serial paired group rejection retains both accepted face families')
      do r=1,size(records)
        call require(all(records(r)%clock_internal_energy%fields(p)%values== &
          kept_records(r)%clock_internal_energy%fields(p)%values).and. &
          all(records(r)%clock_initial%nodes(p)%state==kept_records(r)%clock_initial%nodes(p)%state), &
          'serial paired group provider refusal keeps all previously accepted full-clock records')
      enddo
    enddo
  end subroutine check_serial_held_group

  subroutine check_held_prepared_simultaneous(seed,seed_held,geometry,extended,faces,policy,dt,force)
    type(eb_reactive_tree_3d), intent(in) :: seed
    type(amr_tree_state_3d), intent(in) :: seed_held
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    type(dual_energy_policy), intent(in) :: policy
    real(dp), intent(in) :: dt,force(3)
    type(eb_reactive_tree_3d) :: reference
    type(amr_tree_state_3d) :: reference_held,reference_ir,reference_heat,gas,temp,held,ir,heat
    type(eb_integrated_flux_3d), allocatable :: reference_flux(:),reference_scalar(:),flux(:),scalar(:)
    real(dp) :: reference_correction,correction,b(nv),w(nv),ch(nv),ext(nv),theta,v0(nv),v1(nv),error
    integer :: owners(count),p,layout
    logical :: good
    reference=seed;reference_held=seed_held;reference_heat=seed_held;reference_correction=7
    allocate(reference_ir%fields(count))
    do p=1,count
      allocate(reference_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
      reference_heat%fields(p)%values=0
    enddo
    call integrate_reactive_eb_tree_3d(seed,v0,good)
    call require(good,'held common paired initial inventory')
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,2,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      b,w,ch,ext,theta,good,reference_ir,mol_mode=.true.,redistribution='state_redist', &
      target_volume_fraction=.9_dp,state_redist_max_order=0,internal_energy=reference_held, &
      internal_energy_flux=reference_scalar,energy_policy=policy,energy_correction=reference_correction, &
      integrated_flux=reference_flux,heat_release=reference_heat)
    call require(good,'serial simultaneous prepared held chemistry accepts paired closure')
    call integrate_reactive_eb_tree_3d(reference,v1,good)
    call require(good,'serial simultaneous held inventory')
    v1=v1-v0-b-w-ch-ext;v1(5)=v1(5)-(reference_correction-7)
    error=maxval(abs(v1)/max(1._dp,abs(v0)))
    if(rank==0) print *, 'Serial held common prepared inventory/correction: ',error,reference_correction
    call require(error<3.e-14_dp,'serial simultaneous held independent chemical/numeric inventory')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)];if(layout==2) owners=0
      allocate(gas%fields(count),temp%fields(count),held%fields(count),ir%fields(count),heat%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        gas%fields(p)%values=seed%nodes(p)%state
        allocate(temp%fields(p)%values(1,4,4,4))
        temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
        held%fields(p)=seed_held%fields(p);heat%fields(p)=seed_held%fields(p)
        heat%fields(p)%values=0
        allocate(ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
      enddo
      correction=7
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
        dt,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta, &
        MPI_COMM_WORLD,good,ir,mol_mode=.true.,redistribution='state_redist', &
        target_volume_fraction=.9_dp,state_redist_max_order=0,internal_energy=held,internal_energy_flux=scalar, &
        energy_policy=policy,energy_correction=correction,integrated_flux=flux,heat_release=heat)
      call require(good,'MPI simultaneous prepared held chemistry accepts paired closure')
      call require(abs(correction-reference_correction)<3.e-14_dp*max(1._dp,abs(reference_correction)), &
        'serial/MPI simultaneous prepared held numerical ledger agrees')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,gas,v1,MPI_COMM_WORLD,good)
      call require(good,'MPI simultaneous held inventory')
      v1=v1-v0-b-w-ch-ext;v1(5)=v1(5)-(correction-7)
      error=maxval(abs(v1)/max(1._dp,abs(v0)))
      if(rank==0) print *, 'MPI held common prepared layout/inventory/correction: ',layout,error,correction
      call require(error<3.e-14_dp,'MPI simultaneous held independent chemical/numeric inventory')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(gas%fields(p)%values==reference%nodes(p)%state).and. &
          all(temp%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature).and. &
          all(held%fields(p)%values==reference_held%fields(p)%values).and. &
          all(ir%fields(p)%values==reference_ir%fields(p)%values).and. &
          all(heat%fields(p)%values==reference_heat%fields(p)%values), &
          'serial/MPI simultaneous held prepared gas/T/UEINT/IR/heat agree exactly')
        call require(all(flux(p)%x==reference_flux(p)%x).and.all(flux(p)%y==reference_flux(p)%y).and. &
          all(flux(p)%z==reference_flux(p)%z),'serial/MPI simultaneous prepared physical gas faces agree exactly')
        call require(all(scalar(p)%x==reference_scalar(p)%x).and.all(scalar(p)%y==reference_scalar(p)%y).and. &
          all(scalar(p)%z==reference_scalar(p)%z),'serial/MPI simultaneous prepared scalar faces agree exactly')
      enddo
      deallocate(gas%fields,temp%fields,held%fields,ir%fields,heat%fields)
    enddo
  end subroutine check_held_prepared_simultaneous

  subroutine seed_material_root(seed,stored,geometry,nitrogen,argon,source_temperature)
    type(eb_reactive_tree_3d), intent(inout) :: seed
    type(amr_tree_state_3d), intent(inout) :: stored
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: nitrogen,argon
    real(dp), intent(out) :: source_temperature
    real(dp) :: q(nv),v(nv),temperature,sound,gas_constant,density
    integer :: p,x,y,z
    logical :: good
    q=0;q(1)=.93_dp;q(5)=135000;q(5+nitrogen)=.5_dp;q(5+argon)=.5_dp
    call reactive_primitive_to_conserved(species,q,v,source_temperature,sound,good)
    call require(good,'material source mixed inert target temperature')
    call require(geometry(1)%volume_fraction(material_target(1),material_target(2),material_target(3))>0.and. &
      .not.covered(topology,1,material_target),'material withdrawal owns exposed coarse fluid')
    call require(covered(topology,2,[3,3,2]),'material fixture covered fine source has actual grandchild consumers')
    do p=1,count
    do z=1,4
    do y=1,4
    do x=1,4
      if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
      q=0;density=merge(.93_dp,.31_dp,geometry(p)%volume_fraction(x,y,z)<1)
      if(p==1) then
        q(5+argon)=1
        if(all([x,y,z]==material_target)) then
          q(5+argon)=.5_dp;q(5+nitrogen)=.5_dp
        endif
      else
        density=4;q(5+nitrogen)=1
      endif
      q(1)=density
      gas_constant=mixture_specific_gas_constant(species,q(6:),good)
      call require(good,'material donor inert gas constant')
      q(5)=density*gas_constant*source_temperature
      call reactive_primitive_to_conserved(species,q,v,temperature,sound,good)
      call require(good,'material actual current-stock donor gas/UEINT')
      seed%nodes(p)%state(:,x,y,z)=v;seed%nodes(p)%temperature(x,y,z)=temperature
      stored%fields(p)%values(1,x,y,z)=v(5)
      seed%nodes(p)%state(5,x,y,z)=v(5)-11
    enddo
    enddo
    enddo
    enddo
  end subroutine seed_material_root

  subroutine check_held_natural_replay(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(eb_reactive_tree_3d) :: seed,reference
    type(amr_tree_state_3d) :: seed_held,reference_held,reference_ir,reference_heat,held,gas,temp,ir,heat
    type(eb_mol_held_interval_3d), allocatable :: serial_records(:),trial_records(:)
    type(eb_reactive_tree_3d) :: scalar_ledger,replay_work
    type(amr_tree_state_3d) :: replay_held,replay_ir,mpi_scalar_ledger
    type(eb_integrated_flux_3d), allocatable :: replay_gas_faces(:),replay_held_faces(:)
    real(dp) :: before_h(1),after_h(1),advection(1),volume_rate(1),cell_impulse(1),budget_error
    real(dp) :: duration,cell_volume,saved_start,replay_prefix,replay_chemical(nv)
    integer :: clock,repeat
    type(mpi_eb_mol_interval_3d), allocatable :: records(:)
    type(eb_integrated_flux_3d), allocatable :: reference_flux(:),reference_scalar(:),reference_registered(:)
    type(eb_integrated_flux_3d), allocatable :: flux(:),scalar(:),registered(:)
    type(eb_held_funding_context_3d), allocatable :: serial_context,mpi_context
    type(dual_energy_policy) :: policy
    real(dp) :: dt,q(nv),v(nv),t,cs,v0(nv),v1(nv),error,b(nv),w(nv),ch(nv),ext(nv),theta
    real(dp), allocatable :: aux(:),expected_aux(:)
    real(dp) :: source_temperature,source_r,loan_max,pending_max,nitrogen_loan,held_loan
    real(dp) :: reference_correction,correction
    integer :: owners(count),layout,p,x,y,z,n2,iterations,ios,r,last,argon,bucket
    real(dp) :: projection_total,reset_total,scalar_difference
    integer(int64) :: steps(3),reference_steps(3)
    character(len=64) :: argument
    character(len=16) :: method
    logical :: good
    dt=1.e-8_dp;iterations=1
    call get_command_argument(3,argument)
    if(len_trim(argument)>0) then
      read(argument,*,iostat=ios) dt
      call require(ios==0.and.dt>0,'held natural replay diagnostic physical interval')
    endif
    allocate(aux(nv+merge(1,0,native_material_replay_test)),expected_aux(nv+merge(1,0,native_material_replay_test)))
    n2=0;argon=0
    do p=1,size(species)
      if(trim(species(p)%name)=='N2') n2=p
      if(trim(species(p)%name)=='AR') argon=p
    enddo
    call require(n2>0,'held natural replay initial inert species')
    if(native_material_replay_test) call require(argon>0,'material root fixture identifies inert argon')
    method='none';if(test_mode=='held_prepared_natural') method='state_redist'
    seed=initial;allocate(seed_held%fields(count),reference_ir%fields(count))
    do p=1,count
      allocate(seed_held%fields(p)%values(1,4,4,4),reference_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        q=0;q(1)=merge(.93_dp,.31_dp,geometry(p)%volume_fraction(x,y,z)<1)
        q(5)=135000;q(5+n2)=1
        call reactive_primitive_to_conserved(species,q,v,t,cs,good)
        call require(good,'held natural replay native inert pressure/density state')
        seed%nodes(p)%state(:,x,y,z)=v;seed%nodes(p)%temperature(x,y,z)=t
        seed_held%fields(p)%values(1,x,y,z)=v(5)-.5_dp*sum(v(2:4)**2)/v(1)
        seed%nodes(p)%state(5,x,y,z)=v(5)-11
      enddo
      enddo
      enddo
    enddo
    if(native_material_replay_test) then
      call seed_material_root(seed,seed_held,geometry,n2,argon,source_temperature)
      q=0;q(1)=1;q(5+n2)=1
      source_r=mixture_specific_gas_constant(species,q(6:),good)
      call require(good,'material source pure nitrogen gas constant')
      q(5)=source_r*source_temperature
      call reactive_primitive_to_conserved(species,q,v,t,cs,good)
      call require(good,'material source gas/actual UEINT at source temperature')
      allocate(material_packet(nv+1));material_packet(:nv)=v;material_packet(nv+1)=v(5)
      cell_volume=geometry(1)%volume_fraction(material_target(1),material_target(2),material_target(3))* &
        geometry(1)%dx*geometry(1)%dy*geometry(1)%dz
      expected_aux=-.5_dp*dt*material_rate*cell_volume*material_packet
    endif
    if(native_material_replay_test) then
      call average_down_reactive_eb_tree_3d(species,seed,good,internal_energy=seed_held)
    else
      call average_down_reactive_eb_tree_3d(species,seed,good)
    endif
    call require(good,'held natural source fixture composite initialization')
    policy%eta2=2;policy%allow_negative_energy=.false.
    if(scalar_replay_budget_test) policy%eta2=0
    if(material_refusal_probe) then
      allocate(serial_records(0),records(0),reference_scalar(count),scalar(count))
      do p=1,count
        allocate(reference_scalar(p)%x(1,0:4,4,4),reference_scalar(p)%y(1,4,0:4,4), &
          reference_scalar(p)%z(1,4,4,0:4),source=0._dp)
        if(rank==0) scalar(p)=reference_scalar(p)
      enddo
      call check_held_natural_failures(seed,seed_held,geometry,extended,faces,force,policy,dt, &
        serial_records,records,reference_scalar,scalar)
      if(rank==0) print *, 'PASS_NATIVE_MOL_AUTOMATIC_MATERIAL_REFUSAL'
      deallocate(material_packet)
      return
    endif
    call integrate_reactive_eb_tree_3d(seed,v0,good)
    call require(good,'held natural replay independent initial inventory')
    reference=seed;reference_held=seed_held;reference_heat=seed_held;reference_correction=7
    do p=1,count
      reference_heat%fields(p)%values=0
    enddo
    if(native_funding_replay_test) allocate(serial_context)
    serial_source_calls=0;source_failure=0
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,iterations,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces,b,w,ch,ext,theta,good, &
      reference_ir,mol_mode=.true.,redistribution=method,target_volume_fraction=.9_dp,state_redist_max_order=0, &
      internal_energy=reference_held,internal_energy_flux=reference_scalar, &
      energy_policy=policy,energy_correction=reference_correction,held_interval_history=serial_records, &
      auxiliary_source=held_natural_serial_source,auxiliary_change=aux,level_steps=reference_steps, &
      integrated_flux=reference_flux,registered_flux=reference_registered,heat_release=reference_heat, &
      source_funding_context=serial_context,auxiliary_has_internal_energy=native_material_replay_test)
    if(rank==0) print *, 'Serial held natural acceptance/source calls/weight/correction: ', &
      good,serial_source_calls,theta,reference_correction
    call require(good,'serial held auxiliary native closure accepts complete clock replay')
    if(rank==0) print *, 'Serial held natural replayed physical clocks: ',serial_records%replayed
    if(native_funding_replay_test) then
      call require(all(serial_records%has_source_funding),'automatic native replacement retains funding history')
      call require(theta>=0.and.theta<1,'automatic native parent selects a reduced common spatial weight')
      call require(all(serial_context%open%clock<0),'automatic native root closes all funding face clocks')
      do r=lbound(serial_context%pending,1),ubound(serial_context%pending,1)
        do p=1,count
          call require(all(serial_context%pending(r)%fields(p)%values==0), &
            'automatic native root pays every pending receipt exactly once')
        enddo
      enddo
    endif
    if(native_material_replay_test) then
      call require(.not.serial_records(1)%replayed.and.all(serial_records(2:)%replayed), &
        'material fixture reaches actual failed native sub-parent windows under an open material-funded root')
    else
      call require(all(serial_records%replayed),'serial held fixture reaches actual failed native root-window replay')
    endif
    if(scalar_replay_budget_test) then
    if(rank==0) print *, 'ACTUAL_MOL_REPLAY_START_RESET',serial_records(1)%held_clock_start_change, &
      serial_records(1)%held_clock_reset_after,serial_records(1)%held_clock_projection_after
    call require(abs(serial_records(1)%held_clock_start_change)>3.e-14_dp*max(1._dp,abs(v0(5))), &
      'actual accepted root replay retains a resolved historical initial H repair')
    call require(all(serial_records%held_clock_events_closed).and.all(serial_records%held_clock_prefixes_valid), &
      'accepted actual serial MOL replacement remeasures every scalar entry/exit')
    endif
    call require(serial_source_calls==14.and.all(reference_steps==[1_int64,2_int64,4_int64]), &
      'serial held natural replay freezes source callbacks and keeps physical clocks')
    if(native_material_replay_test) then
      call require(maxval(abs(aux-expected_aux)/max(1._dp,abs(expected_aux)))<3.e-14_dp, &
        'actual nitrogen gas/UEINT source is integrated exactly once')
      loan_max=0;pending_max=0
      do r=1,size(serial_records)
        do p=2,count
          loan_max=max(loan_max,maxval(abs(serial_records(r)%source_funding_impulse%fields(p)%values)))
          do bucket=lbound(serial_context%pending,1),ubound(serial_context%pending,1)
            pending_max=max(pending_max,maxval(abs(serial_records(r)%clock_source_funding_context% &
              pending(bucket)%fields(p)%values)))
          enddo
        enddo
      enddo
      if(rank==0) print *, 'AUTOMATIC_MATERIAL_LOAN_PENDING',loan_max,pending_max
      call require(loan_max>1.e-6_dp.and.pending_max>1.e-6_dp, &
        'automatic accepted replay retains resolved future fine gas/UEINT debt')
      nitrogen_loan=0;held_loan=0
      do p=2,count
        nitrogen_loan=max(nitrogen_loan,maxval(abs(serial_records(1)%source_funding_impulse% &
          fields(p)%values(5+n2,:,:,:))))
        held_loan=max(held_loan,maxval(abs(serial_records(1)%source_funding_impulse%fields(p)%values(nv+1,:,:,:))))
        call require(all(serial_records(2)%clock_source_funding_context%pending(0)%fields(p)%values== &
          serial_records(1)%source_funding_impulse%fields(p)%values), &
          'automatic sub-parent anchor retains the actual coarse material receipt')
        if(topology%nodes(p)%level/=1) cycle
        call require(all(serial_records(2)%initial%nodes(p)%state(6:,:,:,:)== &
          serial_records(2)%clock_initial%nodes(p)%state(6:,:,:,:)+ &
          serial_records(1)%source_funding_impulse%fields(p)%values(6:nv,:,:,:)), &
          'automatic sub-parent entry pays actual species receipt exactly once')
        call require(all(serial_records(2)%initial_internal_energy%fields(p)%values== &
          serial_records(2)%clock_internal_energy%fields(p)%values+ &
          serial_records(1)%source_funding_impulse%fields(p)%values(nv+1:nv+1,:,:,:)), &
          'automatic sub-parent entry pays actual UEINT receipt exactly once')
      enddo
      call require(nitrogen_loan>1.e-6_dp.and.held_loan>1.e-6_dp, &
        'automatic root withdraws resolved nitrogen and actual UEINT from future fine material')
    else
      call require(all(aux==0),'covered predictor injection has no composite physical source inventory')
    endif
    call integrate_reactive_eb_tree_3d(reference,v1,good)
    call require(good,'serial held natural replay final inventory')
    do p=1,size(species)
      if(p==n2.or.(native_material_replay_test.and.p==argon)) cycle
      call require(v1(5+p)==0,'serial covered-source replay leaves exactly no composite chemical products')
    enddo
    v1=v1-v0-b-w-ch-ext-aux(:nv);v1(5)=v1(5)-(reference_correction-7)
    error=maxval(abs(v1)/max(1._dp,abs(v0)))
    if(rank==0) print *, 'Serial held natural replay independent inventory: ',error
    call require(error<3.e-14_dp,'serial held actual failed closure preserves species/energy inventory')
    if(scalar_replay_budget_test.or.native_funding_replay_test) then
    do p=1,count
      if(topology%nodes(p)%level/=2) cycle
      call require(any(reference_scalar(p)%x/=0).or.any(reference_scalar(p)%y/=0).or. &
        any(reference_scalar(p)%z/=0), 'ordinary/prepared MOL replay publishes actual grandchild scalar faces')
      call require(any(reference_flux(p)%x/=0).or.any(reference_flux(p)%y/=0).or. &
        any(reference_flux(p)%z/=0), 'ordinary/prepared MOL replay publishes actual grandchild gas faces')
    enddo
    ! Independent H inventory uses only the original and accepted endpoints,
    ! the complete accepted scalar face family, and each accepted local
    ! pressure/volume rate. No final H is used to derive projection/reset.
    scalar_ledger=seed
    do p=1,count
      scalar_ledger%nodes(p)%state=seed_held%fields(p)%values
    enddo
    call integrate_reactive_eb_tree_3d(scalar_ledger,before_h,good)
    call require(good,'actual MOL replay original independent H inventory')
    do p=1,count
      scalar_ledger%nodes(p)%state=reference_held%fields(p)%values
    enddo
    call integrate_reactive_eb_tree_3d(scalar_ledger,after_h,good)
    call require(good,'actual MOL replay accepted independent H inventory')
    do p=1,count
      scalar_ledger%nodes(p)%state=0
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        cell_volume=geometry(p)%volume_fraction(x,y,z)*geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
        cell_impulse=scalar_face_inventory(geometry(p),reference_scalar(p),x,y,z)
        scalar_ledger%nodes(p)%state(1,x,y,z)=cell_impulse(1)/cell_volume
      enddo
      enddo
      enddo
    enddo
    call integrate_reactive_eb_tree_3d(scalar_ledger,advection,good)
    call require(good,'actual MOL replay complete accepted scalar face inventory')
    do p=1,count
      scalar_ledger%nodes(p)%state=0
    enddo
    do clock=1,size(serial_records)
      duration=serial_records(clock)%interval
      do p=1,count
        if(topology%nodes(p)%level/=serial_records(clock)%depth) cycle
        do z=1,4
        do y=1,4
        do x=1,4
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          cell_volume=geometry(p)%volume_fraction(x,y,z)*geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
          cell_impulse=scalar_face_inventory(geometry(p),serial_records(clock)%internal_energy_flux(p),x,y,z)
          scalar_ledger%nodes(p)%state(1,x,y,z)=scalar_ledger%nodes(p)%state(1,x,y,z)+ &
            duration*serial_records(clock)%internal_energy_spatial_rate%fields(p)%values(1,x,y,z)- &
            cell_impulse(1)/cell_volume
        enddo
        enddo
        enddo
      enddo
    enddo
    call integrate_reactive_eb_tree_3d(scalar_ledger,volume_rate,good)
    call require(good,'actual MOL replay accepted pressure/volume inventory')
    budget_error=abs(after_h(1)-before_h(1)-advection(1)-volume_rate(1)- &
      serial_records(1)%held_clock_projection_after-serial_records(1)%held_clock_reset_after)/ &
      max(1._dp,abs(before_h(1)))
    if(rank==0) print *, 'ACTUAL_MOL_REPLAY_H_BUDGET',budget_error,before_h,after_h,advection,volume_rate, &
      serial_records(1)%held_clock_projection_after,serial_records(1)%held_clock_reset_after
    call require(budget_error<3.e-14_dp,'actual nonzero-entry MOL replacement closes independent H inventory')
    if(method=='none'.and..not.native_funding_replay_test) then
      ! Replaying the canonical root twice retains its historical start once.
      ! Public local groups leave the actual root EOF to their caller.
      trial_records=serial_records;replay_work=reference;replay_held=reference_held;replay_ir=reference_ir
      replay_prefix=reference_correction;replay_chemical=0;saved_start=serial_records(1)%held_clock_start_change
      do repeat=1,2
        call replay_reactive_eb_mol_held_group_3d(species,reactions,replay_work,replay_held,replay_ir, &
          replay_prefix,trial_records,0._dp,replay_gas_faces,replay_held_faces,replay_chemical, &
          2.e-7_dp,1.e-12_dp,good,boundary_faces=faces)
        call require(good,'ordinary actual MOL canonical repeated root replay accepts')
        call require(trial_records(1)%held_clock_start_change==saved_start, &
          'ordinary repeated canonical replay retains historical H start exactly once')
        call require(.not.trial_records(1)%held_clock_events_closed.and. &
          .not.any(trial_records%held_clock_prefixes_valid), &
          'ordinary local group does not publish a root EOF before caller closure')
      enddo
      if(rank==0) print *, 'ACTUAL_MOL_REPEAT_START_RESET',saved_start,trial_records(1)%held_clock_start_change
    endif

    endif
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)];if(layout==2) owners=0
      allocate(gas%fields(count),temp%fields(count),held%fields(count),ir%fields(count),heat%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        gas%fields(p)%values=seed%nodes(p)%state
        allocate(temp%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
        temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
        held%fields(p)=seed_held%fields(p);heat%fields(p)=seed_held%fields(p);heat%fields(p)%values=0
      enddo
      if(native_funding_replay_test) allocate(mpi_context)
      mpi_source_calls=0;correction=7
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
        dt,iterations,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,mol_mode=.true.,redistribution=method,target_volume_fraction=.9_dp, &
        state_redist_max_order=0, &
      internal_energy=held,internal_energy_flux=scalar, &
        energy_policy=policy,energy_correction=correction,interval_history=records, &
        auxiliary_source=held_natural_mpi_source,auxiliary_change=aux,level_steps=steps, &
        integrated_flux=flux,registered_flux=registered,heat_release=heat,source_funding_context=mpi_context, &
        auxiliary_has_internal_energy=native_material_replay_test)
      if(rank==0) print *, 'MPI held natural layout/acceptance/source calls/weight/correction: ', &
        layout,good,mpi_source_calls,theta,correction
      call require(good,'MPI held fixture accepts native failed-closure replay')
      if(native_material_replay_test) then
        call require(.not.records(1)%replayed.and.all(records(2:)%replayed), &
          'MPI owner/idle material fixture reaches automatic sub-parent recovery under an open funded root')
      else
        call require(all(records%replayed),'MPI held fixture reaches native failed-closure complete root replay')
      endif
      call require(all(records%held_clock_events_closed).and.all(records%held_clock_prefixes_valid), &
        'accepted actual owner/idle MOL replacement remeasures every scalar entry/exit')
      do r=1,size(records)
        last=r
        do p=r+1,size(records)
          if(records(p)%depth<=records(r)%depth) exit
          last=p
        enddo
        projection_total=records(r)%held_clock_projection_before+ &
          sum(records(r:last)%internal_energy_projection_change)-sum(records(r:last)%initial_internal_energy_change)
        reset_total=records(r)%held_clock_reset_before+sum(records(r:last)%internal_energy_closure_change)+ &
          sum(records(r:last)%held_clock_end_change)
        if(r<last) reset_total=reset_total+sum(records(r+1:last)%held_clock_start_change)
        call require(max(abs(projection_total-records(r)%held_clock_projection_after), &
          abs(reset_total-records(r)%held_clock_reset_after))/max(1._dp,abs(v0(5)))<3.e-14_dp, &
          'actual MOL replacement subtree ledger uses accepted local and delayed EOF events')
        scalar_difference=maxval(abs([records(r)%held_clock_start_change,records(r)%held_clock_end_change, &
          records(r)%held_clock_projection_before,records(r)%held_clock_reset_before, &
          records(r)%held_clock_projection_after,records(r)%held_clock_reset_after]- &
          [serial_records(r)%held_clock_start_change,serial_records(r)%held_clock_end_change, &
          serial_records(r)%held_clock_projection_before,serial_records(r)%held_clock_reset_before, &
          serial_records(r)%held_clock_projection_after,serial_records(r)%held_clock_reset_after]))
        call require(scalar_difference/max(1._dp,abs(v0(5)))<3.e-14_dp, &
          'actual serial/MPI MOL replacement scalar ledger agrees at original inventory tolerance')
      enddo
      call require(mpi_source_calls==14.and.all(steps==reference_steps), &
        'MPI held natural replay retains source schedule, physical clocks and covered-source exclusion')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,gas,v1,MPI_COMM_WORLD,good)
      call require(good,'MPI held natural replay final inventory')
      do p=1,size(species)
        if(p==n2.or.(native_material_replay_test.and.p==argon)) cycle
        call require(v1(5+p)==0,'MPI covered-source replay leaves exactly no composite chemical products')
      enddo
      if(native_material_replay_test) call require( &
        maxval(abs(aux-expected_aux)/max(1._dp,abs(expected_aux)))<3.e-14_dp, &
        'MPI owner/idle gas/UEINT source is integrated exactly once')
      if(.not.native_material_replay_test) &
        call require(all(aux==0),'MPI covered predictor injection has no physical source inventory')
      v1=v1-v0-b-w-ch-ext-aux(:nv);v1(5)=v1(5)-(correction-7)
      error=maxval(abs(v1)/max(1._dp,abs(v0)))
      if(rank==0) print *, 'MPI held natural layout/independent inventory: ',layout,error
      call require(error<3.e-14_dp,'MPI actual held failed closure independent inventory')
      if(scalar_replay_budget_test.or.native_funding_replay_test) then
      ! Integrate owner-local H and accepted transport directly; a matching
      ! serial endpoint alone does not establish the distributed ledger.
      mpi_scalar_ledger=held
      do p=1,count
        if(owners(p)/=rank) cycle
        mpi_scalar_ledger%fields(p)=seed_held%fields(p)
      enddo
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,mpi_scalar_ledger,before_h,MPI_COMM_WORLD,good)
      call require(good,'actual MPI MOL replay original independent H inventory')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,held,after_h,MPI_COMM_WORLD,good)
      call require(good,'actual MPI MOL replay accepted independent H inventory')
      do p=1,count
        if(owners(p)/=rank) cycle
        mpi_scalar_ledger%fields(p)%values=0
        do z=1,4
        do y=1,4
        do x=1,4
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          cell_volume=geometry(p)%volume_fraction(x,y,z)*geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
          cell_impulse=scalar_face_inventory(geometry(p),scalar(p),x,y,z)
          mpi_scalar_ledger%fields(p)%values(1,x,y,z)=cell_impulse(1)/cell_volume
        enddo
        enddo
        enddo
      enddo
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,mpi_scalar_ledger,advection,MPI_COMM_WORLD,good)
      call require(good,'actual MPI MOL replay complete accepted scalar face inventory')
      do p=1,count
        if(owners(p)==rank) mpi_scalar_ledger%fields(p)%values=0
      enddo
      do clock=1,size(records)
        duration=records(clock)%interval
        do p=1,count
          if(owners(p)/=rank.or.topology%nodes(p)%level/=records(clock)%depth) cycle
          do z=1,4
          do y=1,4
          do x=1,4
            if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
            cell_volume=geometry(p)%volume_fraction(x,y,z)*geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
            cell_impulse=scalar_face_inventory(geometry(p),records(clock)%internal_energy_flux(p),x,y,z)
            mpi_scalar_ledger%fields(p)%values(1,x,y,z)=mpi_scalar_ledger%fields(p)%values(1,x,y,z)+ &
              duration*records(clock)%internal_energy_spatial_rate%fields(p)%values(1,x,y,z)- &
              cell_impulse(1)/cell_volume
          enddo
          enddo
          enddo
        enddo
      enddo
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,mpi_scalar_ledger,volume_rate,MPI_COMM_WORLD,good)
      call require(good,'actual MPI MOL replay accepted pressure/volume inventory')
      budget_error=abs(after_h(1)-before_h(1)-advection(1)-volume_rate(1)- &
        records(1)%held_clock_projection_after-records(1)%held_clock_reset_after)/max(1._dp,abs(before_h(1)))
      if(rank==0) print *, 'ACTUAL_MPI_MOL_REPLAY_H_BUDGET',layout,budget_error
      call require(budget_error<3.e-14_dp,'actual MPI nonzero-entry MOL replacement closes independent H inventory')
      endif

      call require(abs(correction-reference_correction)<3.e-14_dp*max(1._dp,abs(reference_correction)), &
        'serial/MPI natural replay accepted numerical prefix agrees')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(gas%fields(p)%values==reference%nodes(p)%state).and. &
          all(temp%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature).and. &
          all(held%fields(p)%values==reference_held%fields(p)%values).and. &
          all(ir%fields(p)%values==reference_ir%fields(p)%values).and. &
          all(heat%fields(p)%values==reference_heat%fields(p)%values), &
          'serial/MPI native natural replay gas/T/held/IR/heat agree exactly')
        if(any(flux(p)%x/=reference_flux(p)%x).or.any(flux(p)%y/=reference_flux(p)%y).or. &
          any(flux(p)%z/=reference_flux(p)%z)) print *, 'Natural serial/MPI physical flux difference: ',layout,p, &
          maxval(abs(flux(p)%x-reference_flux(p)%x)),maxval(abs(flux(p)%y-reference_flux(p)%y)), &
          maxval(abs(flux(p)%z-reference_flux(p)%z))
        call require(all(flux(p)%x==reference_flux(p)%x).and.all(flux(p)%y==reference_flux(p)%y).and. &
          all(flux(p)%z==reference_flux(p)%z),'serial/MPI natural replay physical gas faces agree exactly')
        call require(all(scalar(p)%x==reference_scalar(p)%x).and.all(scalar(p)%y==reference_scalar(p)%y).and. &
          all(scalar(p)%z==reference_scalar(p)%z),'serial/MPI natural replay scalar faces agree exactly')
        call require(all(registered(p)%x==reference_registered(p)%x).and. &
          all(registered(p)%y==reference_registered(p)%y).and.all(registered(p)%z==reference_registered(p)%z), &
          'serial/MPI natural replay registered gas faces agree exactly')
      enddo
      if(allocated(mpi_context)) then
        call require(all(mpi_context%open%clock<0),'MPI automatic native root closes owner/idle face clocks')
        do r=lbound(mpi_context%pending,1),ubound(mpi_context%pending,1)
          do p=1,count
            if(owners(p)/=rank) cycle
            call require(all(mpi_context%pending(r)%fields(p)%values==0), &
              'MPI automatic native root pays every owner receipt exactly once')
          enddo
        enddo
        if(rank==0) then
          if(native_material_replay_test) then
            print *, 'PASS_NATIVE_MOL_AUTOMATIC_MATERIAL',layout,theta,error
          else
            print *, 'PASS_NATIVE_MOL_AUTOMATIC_PARENT',layout,theta,error
          endif
        endif
        deallocate(mpi_context)
      endif
      deallocate(gas%fields,temp%fields,held%fields,ir%fields,heat%fields)
    enddo
    call check_held_natural_failures(seed,seed_held,geometry,extended,faces,force,policy,dt, &
      serial_records,records,reference_scalar,scalar)
    if(allocated(material_packet)) deallocate(material_packet)
  end subroutine check_held_natural_replay

  function scalar_face_inventory(g,f,x,y,z) result(value)
    type(eb_geometry_3d), intent(in) :: g
    type(eb_integrated_flux_3d), intent(in) :: f
    integer, intent(in) :: x,y,z
    real(dp) :: value(1)
    value=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(:,x-1,y,z)-g%x_face_fraction(x,y,z)*f%x(:,x,y,z))+ &
      g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(:,x,y-1,z)-g%y_face_fraction(x,y,z)*f%y(:,x,y,z))+ &
      g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(:,x,y,z-1)-g%z_face_fraction(x,y,z)*f%z(:,x,y,z))
  end function scalar_face_inventory

  subroutine check_held_natural_failures(seed,seed_held,geometry,extended,faces,force,policy,dt, &
      accepted_serial,accepted_mpi,accepted_scalar,accepted_mpi_scalar)
    type(eb_reactive_tree_3d), intent(in) :: seed
    type(amr_tree_state_3d), intent(in) :: seed_held
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3),dt
    type(dual_energy_policy), intent(in) :: policy
    type(eb_mol_held_interval_3d), intent(in) :: accepted_serial(:)
    type(mpi_eb_mol_interval_3d), intent(in) :: accepted_mpi(:)
    type(eb_integrated_flux_3d), intent(in) :: accepted_scalar(:),accepted_mpi_scalar(:)
    type(eb_held_funding_context_3d), allocatable :: serial_context,mpi_context
    type(eb_reactive_tree_3d) :: serial
    type(amr_tree_state_3d) :: held,ir,heat,gas,temp
    type(eb_mol_held_interval_3d), allocatable :: serial_records(:)
    type(mpi_eb_mol_interval_3d), allocatable :: records(:)
    type(eb_integrated_flux_3d), allocatable :: scalar(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),theta,correction
    real(dp), allocatable :: aux(:)
    integer :: failure,p,r,owners(count)
    character(len=16) :: method
    logical :: good
    procedure(forced_reactor_integrator), pointer :: refusal_reactor
    allocate(aux(nv+merge(1,0,native_material_replay_test)))
    method='none';if(test_mode=='held_prepared_natural') method='state_redist'
    owners=0
    do failure=merge(2,1,material_refusal_probe),2
      nullify(refusal_reactor)
      if(failure==2.or..not.native_material_replay_test) refusal_reactor=>reactor_adapter
      if(native_funding_replay_test) allocate(serial_context,mpi_context)
      serial=seed;held=seed_held;heat=seed_held;serial_records=accepted_serial;scalar=accepted_scalar
      allocate(ir%fields(count))
      do p=1,count
        allocate(ir%fields(p)%values(nv-4,4,4,4),source=0._dp);heat%fields(p)%values=0
      enddo
      serial_source_calls=0;source_failure=0;history_recursive_calls=0
      reactor_replay_failure=failure==2;reactor_failure_rank=rank
      reactor_replay_duration=merge(.5_dp*dt,dt,native_material_replay_test)
      reactor_replay_calls=0;reactor_replay_refusal_call=2
      reactor_replay_refusal_actual=0;reactor_replay_refused=.false.
      if(native_material_replay_test.and.failure==2) then
        ! A supplied reactor owns its signed-species solve. The native path
        ! does not preview it to decide whether composition needs borrowing.
        ! This separate refusal fixture also has a density deficit, so native
        ! current-stock admission must fund it before invoking that provider.
        call require(.5_dp*dt*material_refusal_rate>seed%nodes(1)%state(1, &
          material_target(1),material_target(2),material_target(3)), &
          'material provider refusal fixture requires native density funding')
      endif
      if(failure==1) source_failure=7
      correction=7;b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,1,.true., &
        [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,serial,extended,faces,b,w,ch,ext,theta,good,ir, &
        mol_mode=.true.,redistribution=method,target_volume_fraction=.9_dp,state_redist_max_order=0, &
      internal_energy=held,internal_energy_flux=scalar,energy_policy=policy, &
        energy_correction=correction,held_interval_history=serial_records,heat_release=heat, &
        auxiliary_source=held_natural_serial_source,auxiliary_change=aux,evaluate_reaction=refusal_reactor, &
        source_funding_context=serial_context,auxiliary_has_internal_energy=native_material_replay_test)
      call require(.not.good.and.correction==7.and.all(b==17).and.all(w==19).and.all(ch==23).and. &
        all(ext==29).and.all(aux==31).and.theta==.7_dp,'serial native source/replay refusal retains all accepted ledgers')
      call require(size(serial_records)==size(accepted_serial), 'serial refusal retains accepted clock count')
      call require(all(transfer([serial_records%internal_energy_projection_change, &
          serial_records%internal_energy_closure_change, &
          serial_records%initial_internal_energy_change, &
          serial_records%held_clock_start_change, &
          serial_records%held_clock_end_change, &
          serial_records%held_clock_projection_before, &
          serial_records%held_clock_reset_before, &
          serial_records%held_clock_projection_after, &
          serial_records%held_clock_reset_after], [0_int64],9*size(serial_records))== &
        transfer([accepted_serial%internal_energy_projection_change, &
          accepted_serial%internal_energy_closure_change, &
          accepted_serial%initial_internal_energy_change, &
          accepted_serial%held_clock_start_change, &
          accepted_serial%held_clock_end_change, &
          accepted_serial%held_clock_projection_before, &
          accepted_serial%held_clock_reset_before, &
          accepted_serial%held_clock_projection_after, &
          accepted_serial%held_clock_reset_after], [0_int64],9*size(accepted_serial))).and. &
        all(serial_records%held_clock_events_closed.eqv.accepted_serial%held_clock_events_closed).and. &
        all(serial_records%held_clock_prefixes_valid.eqv.accepted_serial%held_clock_prefixes_valid), &
        'serial refusal retains exact accepted scalar entry/exit and local-repair metadata')
      if(failure==1) then
        if(rank==0) print *, 'MATERIAL_SERIAL_REFUSAL_STAGE',failure,serial_source_calls,reactor_calls
        call require(serial_source_calls==7,'serial held source refusal is not retried')
      else
        if(rank==0) print *, 'MATERIAL_SERIAL_REACTOR_REFUSAL_STAGE',serial_source_calls,reactor_calls, &
          reactor_replay_calls,reactor_replay_refused,reactor_replay_refusal_actual
        call require(serial_source_calls==14.and.reactor_replay_refused.and.reactor_replay_refusal_actual==2, &
          'serial actual full-clock replay provider refusal is native-admitted and not retried')
      endif
      do p=1,count
        call require(all(serial%nodes(p)%state==seed%nodes(p)%state).and. &
          all(serial%nodes(p)%temperature==seed%nodes(p)%temperature).and. &
          all(held%fields(p)%values==seed_held%fields(p)%values).and.all(ir%fields(p)%values==0).and. &
          all(heat%fields(p)%values==0),'serial source/replay refusal discards all completed physical clock fields')
        call require(all(scalar(p)%x==accepted_scalar(p)%x).and.all(scalar(p)%y==accepted_scalar(p)%y).and. &
          all(scalar(p)%z==accepted_scalar(p)%z),'serial failed native group retains accepted scalar faces')
        do r=1,size(serial_records)
          call require(all(serial_records(r)%clock_internal_energy%fields(p)%values== &
            accepted_serial(r)%clock_internal_energy%fields(p)%values), &
            'serial source/replay refusal retains previously accepted complete scalar clocks')
        enddo
      enddo
      if(allocated(serial_context)) call require(.not.allocated(serial_context%pending).and. &
        .not.allocated(serial_context%open),'serial native provider refusal leaves caller funding context unchanged')
      deallocate(ir%fields)
      allocate(gas%fields(count),temp%fields(count),ir%fields(count))
      held=seed_held;heat=seed_held
      do p=1,count
        if(rank/=0) then
          deallocate(held%fields(p)%values,heat%fields(p)%values);cycle
        endif
        gas%fields(p)%values=seed%nodes(p)%state
        allocate(temp%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
        temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature;heat%fields(p)%values=0
      enddo
      records=accepted_mpi;scalar=accepted_mpi_scalar;mpi_source_calls=0;history_recursive_calls=0
      source_failure=0;if(failure==1) source_failure=7
      reactor_failure_rank=0;reactor_replay_calls=0;reactor_replay_refusal_actual=0;reactor_replay_refused=.false.
      correction=7;b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
        dt,1,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,mol_mode=.true.,redistribution=method,target_volume_fraction=.9_dp, &
        state_redist_max_order=0, &
      internal_energy=held,internal_energy_flux=scalar, &
        energy_policy=policy,energy_correction=correction,interval_history=records,heat_release=heat, &
        auxiliary_source=held_natural_mpi_source,auxiliary_change=aux,evaluate_reaction=refusal_reactor, &
        source_funding_context=mpi_context,auxiliary_has_internal_energy=native_material_replay_test)
      call require(.not.good.and.correction==7.and.all(b==17).and.all(w==19).and.all(ch==23).and. &
        all(ext==29).and.all(aux==31).and.theta==.7_dp,'MPI native source/replay refusal retains all accepted ledgers')
      call require(size(records)==size(accepted_mpi), 'MPI owner/idle refusal retains accepted clock count')
      call require(all(transfer([records%internal_energy_projection_change, &
          records%internal_energy_closure_change, &
          records%initial_internal_energy_change, &
          records%held_clock_start_change, &
          records%held_clock_end_change, &
          records%held_clock_projection_before, &
          records%held_clock_reset_before, &
          records%held_clock_projection_after, &
          records%held_clock_reset_after], [0_int64],9*size(records))== &
        transfer([accepted_mpi%internal_energy_projection_change, &
          accepted_mpi%internal_energy_closure_change, &
          accepted_mpi%initial_internal_energy_change, &
          accepted_mpi%held_clock_start_change, &
          accepted_mpi%held_clock_end_change, &
          accepted_mpi%held_clock_projection_before, &
          accepted_mpi%held_clock_reset_before, &
          accepted_mpi%held_clock_projection_after, &
          accepted_mpi%held_clock_reset_after], [0_int64],9*size(accepted_mpi))).and. &
        all(records%held_clock_events_closed.eqv.accepted_mpi%held_clock_events_closed).and. &
        all(records%held_clock_prefixes_valid.eqv.accepted_mpi%held_clock_prefixes_valid), &
        'MPI owner/idle refusal retains exact accepted scalar entry/exit and local-repair metadata')
      if(failure==1) then
        call require(mpi_source_calls==7,'MPI idle source refusal is not retried')
      else
        call require(mpi_source_calls==14,'MPI actual failed-closure replay retains its frozen source schedule')
        if(rank==0) then
          print *, 'MATERIAL_MPI_REACTOR_REFUSAL_STAGE',mpi_source_calls,reactor_replay_calls, &
            reactor_replay_refused,reactor_replay_refusal_actual
          call require(reactor_replay_refused.and.reactor_replay_refusal_actual==2, &
            'MPI actual replay preserves native-admitted provider refusal without another weight')
        endif
      endif
      do p=1,count
        if(rank/=0) cycle
        call require(all(gas%fields(p)%values==seed%nodes(p)%state).and. &
          all(temp%fields(p)%values(1,:,:,:)==seed%nodes(p)%temperature).and. &
          all(held%fields(p)%values==seed_held%fields(p)%values).and.all(ir%fields(p)%values==0).and. &
          all(heat%fields(p)%values==0),'MPI source/replay refusal discards all owned physical clock work')
        call require(all(scalar(p)%x==accepted_mpi_scalar(p)%x).and.all(scalar(p)%y==accepted_mpi_scalar(p)%y).and. &
          all(scalar(p)%z==accepted_mpi_scalar(p)%z),'MPI failed native group retains accepted scalar faces')
        do r=1,size(records)
          call require(all(records(r)%clock_internal_energy%fields(p)%values== &
            accepted_mpi(r)%clock_internal_energy%fields(p)%values), &
            'MPI source/replay refusal retains previously accepted scalar clocks')
        enddo
      enddo
      if(allocated(mpi_context)) then
        call require(.not.allocated(mpi_context%pending).and..not.allocated(mpi_context%open), &
          'MPI native provider refusal leaves owner/idle caller funding context unchanged')
        deallocate(serial_context,mpi_context)
      endif
      deallocate(gas%fields,temp%fields,ir%fields)
    enddo
    reactor_replay_failure=.false.;source_failure=0
  end subroutine check_held_natural_failures

  ! Controlled covered-coarse injection. All source calls report success;
  ! any replay is caused by the native chemistry/reflux closure itself.
  subroutine held_natural_serial_source(state,interval,stage,depth,time,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p,peroxide
    serial_source_calls=serial_source_calls+1;history_recursive_calls=serial_source_calls;peroxide=0
    do p=1,size(species)
      if(trim(species(p)%name)=='H2O2') peroxide=p
    enddo
    call require(peroxide>0.and.interval>0.and.time>=0.and.(stage==0.or.stage==1), &
      'serial held source has the actual old/first-new physical stage')
    do p=1,count
      rate%fields(p)%values=0
    enddo
    if(native_material_replay_test) then
      if(depth==1) then
        rate%fields(2)%values(1,3,3,2)=1000
        rate%fields(2)%values(5+peroxide,3,3,2)=1000
      endif
    else if(depth==0) then
      rate%fields(1)%values(1,2,2,3)=1000
      rate%fields(1)%values(5+peroxide,2,2,3)=1000
    endif
    if(native_material_replay_test.and.depth==0.and.stage==1) then
      if(allocated(rate%fields(1)%values)) rate%fields(1)%values(:, &
        material_target(1),material_target(2),material_target(3))= &
        -merge(material_refusal_rate,material_rate,reactor_replay_failure)*material_packet
    endif
    good=.true.
    if(source_failure>0.and.serial_source_calls==source_failure) good=.false.
  end subroutine held_natural_serial_source

  subroutine held_natural_mpi_source(state,temperature,interval,stage,depth,time,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p,peroxide
    mpi_source_calls=mpi_source_calls+1;history_recursive_calls=mpi_source_calls;peroxide=0
    do p=1,size(species)
      if(trim(species(p)%name)=='H2O2') peroxide=p
    enddo
    call require(peroxide>0.and.interval>0.and.time>=0.and.comm==MPI_COMM_WORLD.and.(stage==0.or.stage==1), &
      'MPI held source has the actual old/first-new physical stage')
    do p=1,count
      if(.not.allocated(state%fields(p)%values)) cycle
      rate%fields(p)%values=0
      if(native_material_replay_test) then
        if(depth==1.and.p==2) then
          rate%fields(p)%values(1,3,3,2)=1000
          rate%fields(p)%values(5+peroxide,3,3,2)=1000
        endif
      else if(depth==0.and.p==1) then
        rate%fields(p)%values(1,2,2,3)=1000
        rate%fields(p)%values(5+peroxide,2,2,3)=1000
      endif
    enddo
    if(native_material_replay_test.and.depth==0.and.stage==1) then
      if(allocated(rate%fields(1)%values)) rate%fields(1)%values(:, &
        material_target(1),material_target(2),material_target(3))= &
        -merge(material_refusal_rate,material_rate,reactor_replay_failure)*material_packet
    endif
    good=.true.
    if(source_failure>0.and.mpi_source_calls==source_failure.and.rank==nranks-1) good=.false.
  end subroutine held_natural_mpi_source

  subroutine check_held_replay(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(amr_tree_state_3d) :: gas,temp,held,ir,heat
    type(mpi_eb_mol_interval_3d), allocatable :: records(:)
    type(eb_integrated_flux_3d), allocatable :: scalar_faces(:),gas_faces(:),registered(:)
    type(dual_energy_policy) :: policy
    real(dp) :: interval,correction,b(nv),w(nv),ch(nv),ext(nv),theta,v0(nv),v1(nv),error
    character(len=64) :: argument
    integer :: owners(count),p,x,y,z,iterations,ios
    integer(int64) :: steps(3)
    logical :: good
    interval=1.e-8_dp;iterations=2
    call get_command_argument(3,argument)
    if(len_trim(argument)>0) then
      read(argument,*,iostat=ios) interval
      call require(ios==0.and.interval>0,'held replay physical clock interval')
    endif
    call get_command_argument(4,argument)
    if(len_trim(argument)>0) then
      read(argument,*,iostat=ios) iterations
      call require(ios==0.and.iterations>0,'held replay physical correction count')
    endif
    owners=[(mod(p-1,nranks),p=1,count)]
    allocate(gas%fields(count),temp%fields(count),held%fields(count),ir%fields(count),heat%fields(count))
    do p=1,count
      if(owners(p)/=rank) cycle
      gas%fields(p)%values=initial%nodes(p)%state
      allocate(temp%fields(p)%values(1,4,4,4),held%fields(p)%values(1,4,4,4), &
        ir%fields(p)%values(nv-4,4,4,4),heat%fields(p)%values(1,4,4,4),source=0._dp)
      temp%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        held%fields(p)%values(1,x,y,z)=gas%fields(p)%values(5,x,y,z)- &
          .5_dp*sum(gas%fields(p)%values(2:4,x,y,z)**2)/gas%fields(p)%values(1,x,y,z)
        gas%fields(p)%values(5,x,y,z)=gas%fields(p)%values(5,x,y,z)-11
      enddo
      enddo
      enddo
    enddo
    policy%allow_negative_energy=.false.;policy%eta2=2;correction=7
    call integrate_mpi_eb_tree_3d(topology,owners,geometry,gas,v0,MPI_COMM_WORLD,good)
    call require(good,'held replay initial composite inventory')
    if(rank==0) print *, 'Held native replay fixture clock/iterations: ',interval,iterations
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
      interval,iterations,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
      b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,mol_mode=.true.,internal_energy=held, &
      internal_energy_flux=scalar_faces,energy_policy=policy,energy_correction=correction, &
      interval_history=records,level_steps=steps,integrated_flux=gas_faces,registered_flux=registered,heat_release=heat, &
      evaluate_reaction=reactor_adapter)
    if(rank==0) print *, 'Held native replay acceptance and correction: ',good,theta,correction
    call require(good,'held physical native recursive transport/chemistry step')
    if(rank==0) print *, 'Held native replay interval records: ',records%replayed
    call integrate_mpi_eb_tree_3d(topology,owners,geometry,gas,v1,MPI_COMM_WORLD,good)
    call require(good,'held replay final composite inventory')
    v1=v1-v0-b-w-ch-ext;v1(5)=v1(5)-(correction-7)
    error=maxval(abs(v1)/max(1._dp,abs(v0)))
    if(rank==0) print *, 'Held native replay independent inventory: ',error
    call require(error<3.e-14_dp,'held native replay preserves independent composite inventory')
    call require(all(steps==[1_int64,2_int64,4_int64]),'held native replay retains actual level clock counts')
  end subroutine check_held_replay

  subroutine check_held_recursive_policy(initial,geometry,extended,faces)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    type(dual_energy_policy) :: policy,rank_policy
    type(eb_reactive_tree_3d) :: seed,reference,failed
    type(amr_tree_state_3d) :: seed_held,reference_held,gas,temp,held,inventory_gas
    type(eb_integrated_flux_3d), allocatable :: reference_flux(:),scalar_flux(:)
    type(eb_geometry_3d), allocatable :: invalid_geometry(:)
    type(eb_integrated_flux_3d), allocatable :: reference_registered(:),reference_physical(:)
    type(eb_integrated_flux_3d), allocatable :: registered(:),physical(:)
    integer(int64) :: steps(3),reference_steps(3)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),theta,sb(nv),sw(nv),sch(nv),sext(nv),stheta
    real(dp) :: before(nv),after(nv),error,correction,reference_correction,uniform(nv),uniform_t,energy,volume
    real(dp), parameter :: interval=1.e-7_dp,force(3)=[2.e9_dp,-3.e9_dp,0._dp]
    integer :: owners(count),p,x,y,z,layout
    logical :: good
    seed=initial;uniform=initial%nodes(1)%state(:,1,1,4);uniform_t=initial%nodes(1)%temperature(1,1,4)
    call require(geometry(1)%volume_fraction(1,1,4)>0,'held recursive policy uniform fluid seed')
    energy=uniform(5)-.5_dp*sum(uniform(2:4)**2)/uniform(1)
    allocate(seed_held%fields(count))
    volume=0
    do p=1,count
      allocate(seed_held%fields(p)%values(1,4,4,4),source=0._dp)
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        seed%nodes(p)%state(:,x,y,z)=uniform;seed%nodes(p)%temperature(x,y,z)=uniform_t
        seed_held%fields(p)%values(1,x,y,z)=energy
        seed%nodes(p)%state(5,x,y,z)=uniform(5)-11
        if(.not.covered(topology,p,[x,y,z])) volume=volume+geometry(p)%volume_fraction(x,y,z)* &
          geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
      enddo
      enddo
      enddo
    enddo
    call integrate_reactive_eb_tree_3d(seed,before,good)
    call require(good,'held recursive policy initial depleted gas inventory')
    allocate(inventory_gas%fields(count))
    do p=1,count
      inventory_gas%fields(p)%values=seed%nodes(p)%state
    enddo
    call held_policy_compensated_inventory(inventory_gas,geometry,before)
    policy%allow_negative_energy=.false.;policy%eta2=2
    reference=seed;reference_held=seed_held;reference_correction=7
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,interval,2,.false., &
      [.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,stheta,good,mol_mode=.true.,internal_energy=reference_held, &
      internal_energy_flux=reference_flux,energy_policy=policy,energy_correction=reference_correction, &
      level_steps=reference_steps,registered_flux=reference_registered,integrated_flux=reference_physical)
    call require(good,'serial held recursive policy recursive MOL step')
    call require(all(reference_steps==[1_int64,2_int64,4_int64]),'held recursive policy level clocks')
    call integrate_reactive_eb_tree_3d(reference,after,good)
    call require(good,'serial held recursive policy final inventory')
    do p=1,count
      inventory_gas%fields(p)%values=reference%nodes(p)%state
    enddo
    call held_policy_compensated_inventory(inventory_gas,geometry,after)
    after=after-before-sb-sw-sch-sext;after(5)=after(5)-(reference_correction-7)
    error=maxval(abs(after)/max(1._dp,abs(before)))
    print *, 'Held recursive policy serial closure and correction: ',error,reference_correction
    print *, 'Held recursive policy residual: ',after
    print *, 'Held recursive policy initial inventory: ',before
    print *, 'Held recursive policy source inventory: ',sext
    call require(error<3.e-14_dp,'serial held recursive policy independent numerical-energy budget')
    call require(reference_correction>7,'held recursive policy charges an actual numerical repair')
    print *, 'Held recursive correction beyond initial reset: ',reference_correction-7-11*volume
    call require(maxval(abs(sch)/max(1._dp,abs(before)))<3.e-14_dp, &
      'held numerical-energy repair is not chemical exchange')
    invalid_geometry=extended;invalid_geometry(count)%nx=invalid_geometry(count)%nx-1
    failed=seed;held=seed_held;correction=7;b=17;w=19;ch=23;ext=29;theta=.7_dp;steps=13
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,interval,2,.false., &
      [.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,failed,invalid_geometry,faces, &
      b,w,ch,ext,theta,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux, &
      energy_policy=policy,energy_correction=correction, &
      level_steps=steps,registered_flux=registered,integrated_flux=physical)
    call require(.not.good.and.correction==7,'serial failure after initial reset retains numerical ledger')
    call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.theta==.7_dp.and. &
      .not.allocated(scalar_flux),'serial held recursive policy failure retains diagnostics and unallocated flux')
    call require(all(steps==13).and..not.allocated(registered).and..not.allocated(physical), &
      'serial held recursive policy failure retains counters and flux families')
    do p=1,count
      call require(all(failed%nodes(p)%state==seed%nodes(p)%state).and. &
        all(failed%nodes(p)%temperature==seed%nodes(p)%temperature).and. &
        all(held%fields(p)%values==seed_held%fields(p)%values),'serial initial reset rejection keeps gas and held scalar')
    enddo
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)];if(layout==2) owners=0
      allocate(gas%fields(count),temp%fields(count))
      if(allocated(held%fields)) deallocate(held%fields)
      allocate(held%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        gas%fields(p)%values=seed%nodes(p)%state
        allocate(temp%fields(p)%values(1,4,4,4));temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
        held%fields(p)=seed_held%fields(p)
      enddo
      correction=7
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
        interval,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux, &
        energy_policy=policy,energy_correction=correction, &
        level_steps=steps,registered_flux=registered,integrated_flux=physical)
      call require(good,'MPI held recursive policy recursive MOL step')
      call require(all(steps==reference_steps),'MPI held recursive policy level clocks')
      call require(abs(correction-reference_correction)<3.e-14_dp*max(1._dp,abs(before(5))), &
        'MPI held recursive policy composite correction matches serial')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(gas%fields(p)%values==reference%nodes(p)%state).and. &
          all(temp%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature).and. &
          all(held%fields(p)%values==reference_held%fields(p)%values),'MPI held recursive policy exact owner/idle state')
        call require(all(scalar_flux(p)%x==reference_flux(p)%x).and.all(scalar_flux(p)%y==reference_flux(p)%y).and. &
          all(scalar_flux(p)%z==reference_flux(p)%z),'MPI held recursive policy exact scalar face family')
        call require(all(registered(p)%x==reference_registered(p)%x).and. &
          all(registered(p)%y==reference_registered(p)%y).and.all(registered(p)%z==reference_registered(p)%z), &
          'MPI held recursive policy exact registered faces')
        call require(all(physical(p)%x==reference_physical(p)%x).and. &
          all(physical(p)%y==reference_physical(p)%y).and.all(physical(p)%z==reference_physical(p)%z), &
          'MPI held recursive policy exact physical faces')
      enddo
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,gas,after,MPI_COMM_WORLD,good)
      call require(good,'MPI held recursive policy final inventory')
      call held_policy_compensated_inventory(gas,geometry,after,MPI_COMM_WORLD)
      after=after-before-b-w-ch-ext;after(5)=after(5)-(correction-7)
      error=maxval(abs(after)/max(1._dp,abs(before)))
      call require(error<3.e-14_dp,'MPI held recursive policy independent numerical-energy budget')
      do p=1,count
        if(owners(p)/=rank) cycle
        gas%fields(p)%values=seed%nodes(p)%state
        temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
        held%fields(p)=seed_held%fields(p)
      enddo
      correction=7;b=17;w=19;ch=23;ext=29;theta=.7_dp;steps=13
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
        interval,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,invalid_geometry,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux, &
        energy_policy=policy,energy_correction=correction, &
        level_steps=steps,registered_flux=registered,integrated_flux=physical)
      call require(.not.good.and.correction==7,'MPI initial reset geometry rejection retains numerical ledger')
      rank_policy=policy;if(rank==nranks-1) rank_policy%eta2=3
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
        interval,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux, &
        energy_policy=rank_policy,energy_correction=correction, &
        level_steps=steps,registered_flux=registered,integrated_flux=physical)
      call require(.not.good.and.correction==7,'MPI policy mismatch including idle rank rejects collectively')
      if(rank==nranks-1) correction=8
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
        interval,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux, &
        energy_policy=policy,energy_correction=correction, &
        level_steps=steps,registered_flux=registered,integrated_flux=physical)
      call require(.not.good.and.correction==merge(8._dp,7._dp,rank==nranks-1), &
        'MPI held recursive policy mismatched correction prefix rejects collectively')
      correction=7
      if(rank==nranks-1) then
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
          interval,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
          b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux, &
          energy_correction=correction,level_steps=steps,registered_flux=registered,integrated_flux=physical)
      else
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
          interval,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
          b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux, &
          energy_policy=policy,energy_correction=correction, &
          level_steps=steps,registered_flux=registered,integrated_flux=physical)
      endif
      call require(.not.good.and.correction==7,'MPI held recursive policy presence mismatch rejects collectively')
      call require(all(steps==13),'MPI held recursive policy failure retains level counters')
      call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.theta==.7_dp, &
        'MPI held recursive policy failure retains all diagnostics')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(gas%fields(p)%values==seed%nodes(p)%state).and. &
          all(temp%fields(p)%values(1,:,:,:)==seed%nodes(p)%temperature).and. &
          all(held%fields(p)%values==seed_held%fields(p)%values),'MPI held recursive policy failure keeps all owner data')
        call require(all(scalar_flux(p)%x==reference_flux(p)%x).and.all(scalar_flux(p)%y==reference_flux(p)%y).and. &
          all(scalar_flux(p)%z==reference_flux(p)%z),'MPI held recursive policy failure keeps published scalar faces')
        call require(all(registered(p)%x==reference_registered(p)%x).and. &
          all(registered(p)%y==reference_registered(p)%y).and.all(registered(p)%z==reference_registered(p)%z), &
          'MPI held recursive policy failure keeps registered faces')
        call require(all(physical(p)%x==reference_physical(p)%x).and. &
          all(physical(p)%y==reference_physical(p)%y).and.all(physical(p)%z==reference_physical(p)%z), &
          'MPI held recursive policy failure keeps physical faces')
      enddo
      deallocate(gas%fields,temp%fields,held%fields)
    enddo
  end subroutine check_held_recursive_policy


  ! Independent inventory oracle: compensate leaf-cell summation so a large
  ! forced momentum increment does not turn summation drift into a residual.
  subroutine held_policy_compensated_inventory(gas,geometry,values,comm)
    type(amr_tree_state_3d), intent(in) :: gas
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(out) :: values(nv)
    type(MPI_Comm), optional, intent(in) :: comm
    real(dp) :: local(nv),carry(nv),term(nv),next(nv),volume
    integer :: p,x,y,z,status
    local=0;carry=0
    do p=1,count
      if(.not.allocated(gas%fields(p)%values)) cycle
      volume=geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0.or.covered(topology,p,[x,y,z])) cycle
        term=volume*geometry(p)%volume_fraction(x,y,z)*gas%fields(p)%values(:,x,y,z)-carry
        next=local+term;carry=(next-local)-term;local=next
      enddo
      enddo
      enddo
    enddo
    values=local
    if(present(comm)) then
      call MPI_Allreduce(local,values,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
      call require(status==MPI_SUCCESS,'held policy compensated inventory reduction')
    endif
  end subroutine held_policy_compensated_inventory

  subroutine check_stored_mol(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(eb_reactive_tree_3d) :: serial_gas
    type(amr_tree_state_3d) :: gas,temp,held,initial_held,reference_held,saved_gas,saved_held,serial_held,serial_initial
    type(eb_integrated_flux_3d), allocatable :: scalar_flux(:),reference_flux(:),saved_flux(:),serial_flux(:)
    type(eb_integrated_flux_3d), allocatable :: serial_registered(:),serial_physical(:),mpi_registered(:),mpi_physical(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),theta,kinetic,sb(nv),sw(nv),sch(nv),sext(nv),stheta
    real(dp) :: initial_inventory(nv),final_inventory(nv),error
    integer :: owners(count),p,x,y,z,layout,mpi_status
    logical :: good,changed,global_changed
    allocate(reference_held%fields(count),reference_flux(count),serial_held%fields(count))
    serial_gas=initial
    do p=1,count
      allocate(serial_held%fields(p)%values(1,4,4,4),source=0._dp)
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        kinetic=.5_dp*sum(initial%nodes(p)%state(2:4,x,y,z)**2)/initial%nodes(p)%state(1,x,y,z)
        serial_held%fields(p)%values(1,x,y,z)=(initial%nodes(p)%state(5,x,y,z)-kinetic)* &
          (1._dp+.1_dp*real(x-1,dp)+.2_dp*real(z-1,dp))
      enddo
      enddo
      enddo
    enddo
    serial_initial=serial_held
    sb=0;sw=0;sch=0;sext=0;stheta=1
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.false., &
      [.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,serial_gas,extended,faces, &
      sb,sw,sch,sext,stheta,good,mol_mode=.true.,internal_energy=serial_held,internal_energy_flux=serial_flux)
    call require(good,'serial stored internal-energy MOL common step')
    changed=.false.
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(gas%fields(count),temp%fields(count),held%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        gas%fields(p)%values=initial%nodes(p)%state
        allocate(temp%fields(p)%values(1,4,4,4))
        temp%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
        held%fields(p)=serial_initial%fields(p)
      enddo
      initial_held=held
      b=0;w=0;ch=0;ext=0;theta=1
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
        1.e-10_dp,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held, &
        internal_energy_flux=scalar_flux)
      call require(good,'MPI stored internal-energy MOL common step')
      do p=1,count
        call require(allocated(held%fields(p)%values).eqv.(owners(p)==rank), &
          'stored internal-energy owner-local output')
        call require(allocated(scalar_flux(p)%x).eqv.(owners(p)==rank), &
          'stored internal-energy owner-local flux')
        if(owners(p)==rank) then
          call require(all(held%fields(p)%values==serial_held%fields(p)%values), &
            'stored internal-energy exact serial/MPI state')
          call require(all(scalar_flux(p)%x==serial_flux(p)%x).and. &
            all(scalar_flux(p)%y==serial_flux(p)%y).and. &
            all(scalar_flux(p)%z==serial_flux(p)%z),'stored internal-energy exact serial/MPI flux')
        endif
        if(layout==1) then
          allocate(reference_held%fields(p)%values(1,4,4,4))
          allocate(reference_flux(p)%x(1,0:4,4,4),reference_flux(p)%y(1,4,0:4,4), &
            reference_flux(p)%z(1,4,4,0:4))
          if(owners(p)==rank) then
            reference_held%fields(p)%values=held%fields(p)%values
            reference_flux(p)=scalar_flux(p)
            if(p==count) changed=changed.or.any(held%fields(p)%values/=initial_held%fields(p)%values)
          endif
          call MPI_Bcast(reference_held%fields(p)%values,64,MPI_DOUBLE_PRECISION,owners(p),MPI_COMM_WORLD,mpi_status)
          call require(mpi_status==MPI_SUCCESS,'stored internal-energy reference field broadcast')
          call MPI_Bcast(reference_flux(p)%x,size(reference_flux(p)%x),MPI_DOUBLE_PRECISION, &
            owners(p),MPI_COMM_WORLD,mpi_status)
          call require(mpi_status==MPI_SUCCESS,'stored internal-energy x-flux broadcast')
          call MPI_Bcast(reference_flux(p)%y,size(reference_flux(p)%y),MPI_DOUBLE_PRECISION, &
            owners(p),MPI_COMM_WORLD,mpi_status)
          call require(mpi_status==MPI_SUCCESS,'stored internal-energy y-flux broadcast')
          call MPI_Bcast(reference_flux(p)%z,size(reference_flux(p)%z),MPI_DOUBLE_PRECISION, &
            owners(p),MPI_COMM_WORLD,mpi_status)
          call require(mpi_status==MPI_SUCCESS,'stored internal-energy z-flux broadcast')
        else
          if(owners(p)/=rank) cycle
          call require(all(held%fields(p)%values==reference_held%fields(p)%values), &
            'stored internal-energy rank-layout field equality')
          call require(all(scalar_flux(p)%x==reference_flux(p)%x).and. &
            all(scalar_flux(p)%y==reference_flux(p)%y).and. &
            all(scalar_flux(p)%z==reference_flux(p)%z),'stored internal-energy rank-layout flux equality')
        endif
      enddo
      if(layout==1.and.nranks>1) then
        saved_gas=gas;saved_held=held;saved_flux=scalar_flux
        b=17;w=19;ch=23;ext=29;theta=.7_dp
        if(rank==0) then
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
            1.e-10_dp,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces, &
            b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held, &
            internal_energy_flux=scalar_flux)
        else
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
            1.e-10_dp,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces, &
            b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.)
        endif
        call require(.not.good,'stored internal-energy rank-disagreeing option rejected')
        call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.theta==.7_dp, &
          'stored internal-energy rejected ledger rollback')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(gas%fields(p)%values==saved_gas%fields(p)%values).and. &
            all(held%fields(p)%values==saved_held%fields(p)%values), &
            'stored internal-energy rejected state rollback')
          call require(all(scalar_flux(p)%x==saved_flux(p)%x).and. &
            all(scalar_flux(p)%y==saved_flux(p)%y).and. &
            all(scalar_flux(p)%z==saved_flux(p)%z),'stored internal-energy rejected flux rollback')
        enddo
      endif
      deallocate(gas%fields,temp%fields,held%fields,scalar_flux)
    enddo
    call MPI_Allreduce(changed,global_changed,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_status)
    call require(mpi_status==MPI_SUCCESS.and.global_changed,'stored internal-energy step is nontrivial')
    saved_flux=serial_flux
    serial_gas=initial;serial_held=serial_initial
    sb=0;sw=0;sch=0;sext=0;stheta=1
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.false., &
      [.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,serial_gas,extended,faces, &
      sb,sw,sch,sext,stheta,good,mol_mode=.true.,redistribution='state_redist',target_volume_fraction=.8_dp, &
      internal_energy=serial_held,internal_energy_flux=serial_flux)
    call require(good,'boundary-touching serial StateRedist accepted')
    call integrate_reactive_eb_tree_3d(initial,initial_inventory,good)
    call require(good,'boundary-touching initial inventory')
    call integrate_reactive_eb_tree_3d(serial_gas,final_inventory,good)
    call require(good,'boundary-touching serial inventory')
    error=maxval(abs(final_inventory-initial_inventory-sb-sw-sch-sext)/max(1._dp,abs(initial_inventory)))
    call require(error<3.e-14_dp,'boundary-touching serial composite closure')
    owners=[(mod(p-1,nranks),p=1,count)]
    allocate(gas%fields(count),temp%fields(count),held%fields(count))
    do p=1,count
      if(owners(p)/=rank) cycle
      gas%fields(p)%values=initial%nodes(p)%state
      allocate(temp%fields(p)%values(1,4,4,4))
      temp%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
      held%fields(p)=serial_initial%fields(p)
    enddo
    b=0;w=0;ch=0;ext=0;theta=1
    call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
      1.e-10_dp,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,faces, &
      b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,redistribution='state_redist', &
      target_volume_fraction=.8_dp, &
      internal_energy=held,internal_energy_flux=scalar_flux)
    call require(good,'boundary-touching MPI StateRedist accepted')
    do p=1,count
      if(owners(p)/=rank) cycle
      call require(all(gas%fields(p)%values==serial_gas%nodes(p)%state).and. &
        all(held%fields(p)%values==serial_held%fields(p)%values), &
        'boundary-touching serial MPI gas and held parity')
      call require(all(scalar_flux(p)%x==serial_flux(p)%x).and. &
        all(scalar_flux(p)%y==serial_flux(p)%y).and.all(scalar_flux(p)%z==serial_flux(p)%z), &
        'boundary-touching serial MPI scalar flux parity')
    enddo
    call integrate_mpi_eb_tree_3d(topology,owners,geometry,gas,final_inventory,MPI_COMM_WORLD,good)
    call require(good,'boundary-touching MPI inventory')
    error=maxval(abs(final_inventory-initial_inventory-b-w-ch-ext)/max(1._dp,abs(initial_inventory)))
    call require(error<3.e-14_dp,'boundary-touching MPI composite closure')
    deallocate(gas%fields,temp%fields,held%fields)
    serial_gas=initial;serial_held=serial_initial
    sb=0;sw=0;sch=0;sext=0;stheta=1
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,2,.false., &
      [.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,serial_gas,extended,faces, &
      sb,sw,sch,sext,stheta,good,mol_mode=.true.,internal_energy=serial_held,internal_energy_flux=serial_flux)
    call require(good,'serial stored internal-energy recursive MOL step')
    call integrate_reactive_eb_tree_3d(serial_gas,final_inventory,good)
    call require(good,'stored internal-energy recursive baseline inventory')
    owners=[(mod(p-1,nranks),p=1,count)]
    allocate(gas%fields(count),temp%fields(count),held%fields(count))
    do p=1,count
      if(owners(p)/=rank) cycle
      gas%fields(p)%values=initial%nodes(p)%state
      allocate(temp%fields(p)%values(1,4,4,4))
      temp%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
      held%fields(p)=serial_initial%fields(p)
    enddo
    b=0;w=0;ch=0;ext=0;theta=1
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
      1.e-10_dp,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
      b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux)
    call require(good,'MPI stored internal-energy recursive MOL step')
    do p=1,count
      if(owners(p)/=rank) cycle
      call require(all(held%fields(p)%values==serial_held%fields(p)%values), &
        'stored internal-energy recursive exact serial/MPI state')
      call require(all(scalar_flux(p)%x==serial_flux(p)%x).and. &
        all(scalar_flux(p)%y==serial_flux(p)%y).and. &
        all(scalar_flux(p)%z==serial_flux(p)%z),'stored internal-energy recursive exact serial/MPI flux')
    enddo
    if(nranks>1) then
      saved_gas=gas;saved_held=held;saved_flux=scalar_flux
      b=17;w=19;ch=23;ext=29;theta=.7_dp
      if(rank==0) then
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
          1.e-10_dp,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
          b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,internal_energy=held,internal_energy_flux=scalar_flux)
      else
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
          1.e-10_dp,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
          b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.)
      endif
      call require(.not.good,'recursive stored internal-energy rank-disagreeing option rejected')
      call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.theta==.7_dp, &
        'recursive stored internal-energy rejected ledger rollback')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(gas%fields(p)%values==saved_gas%fields(p)%values).and. &
          all(held%fields(p)%values==saved_held%fields(p)%values), &
          'recursive stored internal-energy rejected state rollback')
        call require(all(scalar_flux(p)%x==saved_flux(p)%x).and. &
          all(scalar_flux(p)%y==saved_flux(p)%y).and. &
          all(scalar_flux(p)%z==saved_flux(p)%z),'recursive stored internal-energy rejected flux rollback')
      enddo
    endif
    deallocate(gas%fields,temp%fields,held%fields,scalar_flux)
    serial_gas=initial;serial_held=serial_initial
    sb=0;sw=0;sch=0;sext=0;stheta=1
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,2,.false., &
      [.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,serial_gas,extended,faces, &
      sb,sw,sch,sext,stheta,good,mol_mode=.true.,redistribution='state_redist',target_volume_fraction=.8_dp, &
      internal_energy=serial_held,internal_energy_flux=serial_flux,registered_flux=serial_registered, &
      integrated_flux=serial_physical)
    call require(good,'serial boundary-touching recursive held StateRedist')
    call require(any(serial_registered(1)%x/=serial_physical(1)%x), &
      'recursive held StateRedist preserves distinct registered and physical gas fluxes')
    call integrate_reactive_eb_tree_3d(serial_gas,final_inventory,good)
    call require(good,'recursive boundary-touching serial inventory')
    error=maxval(abs(final_inventory-initial_inventory-sb-sw-sch-sext)/max(1._dp,abs(initial_inventory)))
    call require(error<3.e-14_dp,'recursive boundary-touching serial composite closure')
    owners=[(mod(p-1,nranks),p=1,count)]
    allocate(gas%fields(count),temp%fields(count),held%fields(count))
    do p=1,count
      if(owners(p)/=rank) cycle
      gas%fields(p)%values=initial%nodes(p)%state
      allocate(temp%fields(p)%values(1,4,4,4))
      temp%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
      held%fields(p)=serial_initial%fields(p)
    enddo
    b=0;w=0;ch=0;ext=0;theta=1
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
      1.e-10_dp,2,.false.,[.false.,.false.,.false.,.false.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
      b,w,ch,ext,theta,MPI_COMM_WORLD,good,mol_mode=.true.,redistribution='state_redist', &
      target_volume_fraction=.8_dp,internal_energy=held,internal_energy_flux=scalar_flux, &
      registered_flux=mpi_registered,integrated_flux=mpi_physical)
    call require(good,'MPI boundary-touching recursive held StateRedist')
    do p=1,count
      if(owners(p)/=rank) cycle
      error=maxval(abs(gas%fields(p)%values-serial_gas%nodes(p)%state)/ &
        max(1._dp,abs(serial_gas%nodes(p)%state)))
      call require(error<2.e-14_dp,'recursive boundary-touching serial MPI gas parity')
      call require(all(held%fields(p)%values==serial_held%fields(p)%values), &
        'recursive boundary-touching serial MPI held parity')
      call require(all(scalar_flux(p)%x==serial_flux(p)%x).and. &
        all(scalar_flux(p)%y==serial_flux(p)%y).and.all(scalar_flux(p)%z==serial_flux(p)%z), &
        'recursive boundary-touching serial MPI scalar flux parity')
      call require(all(mpi_registered(p)%x==serial_registered(p)%x).and. &
        all(mpi_registered(p)%y==serial_registered(p)%y).and.all(mpi_registered(p)%z==serial_registered(p)%z), &
        'recursive boundary-touching serial MPI registered gas flux parity')
      call require(all(mpi_physical(p)%x==serial_physical(p)%x).and. &
        all(mpi_physical(p)%y==serial_physical(p)%y).and.all(mpi_physical(p)%z==serial_physical(p)%z), &
        'recursive boundary-touching serial MPI physical gas flux parity')
    enddo
    call integrate_mpi_eb_tree_3d(topology,owners,geometry,gas,final_inventory,MPI_COMM_WORLD,good)
    call require(good,'recursive boundary-touching MPI inventory')
    error=maxval(abs(final_inventory-initial_inventory-b-w-ch-ext)/max(1._dp,abs(initial_inventory)))
    call require(error<3.e-14_dp,'recursive boundary-touching MPI composite closure')
  end subroutine check_stored_mol

  subroutine check_sdc_recursive(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    procedure(eb_sdc_recursive_source), pointer :: serial_provider=>null()
    procedure(mpi_eb_sdc_recursive_source), pointer :: mpi_provider=>null()
    real(dp), allocatable :: aux(:),saux(:)
    real(dp) :: coefficient,expected_aux(nv),source_factor
    integer :: x,y,z,failed_order
    type(eb_reactive_tree_3d) :: reference,profile
    type(amr_tree_state_3d) :: local,lt,ir,serial_ir,old,ot,old_ir,heat,serial_heat,old_heat
    integer(int64) :: counts(3),serial_counts(3)
    integer :: owners(count),p,layout
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),sb(nv),sw(nv),sch(nv),sext(nv),limiter,slimiter,error
    real(dp) :: initial_inventory(nv),final_inventory(nv)
    character(len=16) :: method
    logical :: good
    method="none"
    if(test_mode=="redist_recursive") method="flux_redist"
    if(test_mode=="state_recursive") method="state_redist"
    nullify(serial_provider,mpi_provider)
    if(test_mode=='aux_recursive'.or.test_mode=='aux_route') then
      serial_provider=>serial_recursive_aux;mpi_provider=>mpi_recursive_aux
      allocate(aux(nv),saux(nv))
    end if
    serial_source_calls=0;source_failure=0
    reference=initial
    allocate(serial_ir%fields(count),serial_heat%fields(count))
    do p=1,count
      allocate(serial_heat%fields(p)%values(1,4,4,4),source=0._dp)
      allocate(serial_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
    end do
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,slimiter,good,serial_ir,level_steps=serial_counts,redistribution=method, &
      auxiliary_source=serial_provider,auxiliary_change=saux,route_auxiliary=test_mode=='aux_route', &
      evaluate_reaction=reactor_adapter, &
      state_redist_max_order=requested_order,heat_release=serial_heat)
    call require(good.and.reactor_calls>0,'serial recursive SDC callback')
    call require_restricted_reactions(reference,serial_ir)
    call require_restricted_reactions(reference,serial_heat)
    reactor_calls=0
    if(allocated(saux)) then
      call require(serial_source_calls==21,'serial recursive source calls')
      profile=initial
      do p=1,count
        source_factor=3._dp-1._dp/2**topology%nodes(p)%level
        if(test_mode=='aux_route') source_factor=7.25_dp
        do z=1,4
        do y=1,4
        do x=1,4
          profile%nodes(p)%state(:,x,y,z)=1e-10_dp*source_factor* &
            source_profile(p,x,y,z,0)
          if(profile%nodes(p)%geometry%volume_fraction(x,y,z)<=0) profile%nodes(p)%state(:,x,y,z)=0
        end do
        end do
        end do
      end do
      call integrate_reactive_eb_tree_3d(profile,expected_aux,good)
      call require(good,'recursive expected inventory')
      call require(maxval(abs(saux-expected_aux)/max(1e-30_dp,abs(expected_aux)))<3e-14_dp, &
        'recursive sources integrate accepted stages and leaf time intervals')
    end if
    call require(all(serial_counts==[1_int64,2_int64,4_int64]),'serial level counts')
    call integrate_reactive_eb_tree_3d(initial,initial_inventory,good)
    call require(good,'recursive initial inventory')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(local%fields(count),lt%fields(count),ir%fields(count),heat%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        allocate(heat%fields(p)%values(1,4,4,4),source=0._dp)
        local%fields(p)%values=initial%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4))
        lt%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
        ir%fields(p)%values=0
      end do
      mpi_source_calls=0;source_failure=0
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
        1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,level_steps=counts,redistribution=method, &
        auxiliary_source=mpi_provider,auxiliary_change=aux,route_auxiliary=test_mode=='aux_route', &
      state_redist_max_order=requested_order,heat_release=heat)
      call require(good,'MPI recursive SDC')
      if(allocated(aux)) then
        call require(mpi_source_calls==21,'MPI recursive source calls')
        call require(maxval(abs(aux-saux)/max(1e-30_dp,abs(saux)))<3e-14_dp,'recursive auxiliary serial MPI')
      end if
      call require(all(counts==serial_counts),'MPI recursive level counts')
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'recursive owner-local fields')
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==reference%nodes(p)%state),'exact recursive state')
        call require(all(lt%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'exact recursive T')
        call require(all(ir%fields(p)%values==serial_ir%fields(p)%values),'exact recursive IR')
        call require(all(heat%fields(p)%values==serial_heat%fields(p)%values),'exact recursive heat')
      end do
      call require(limiter==slimiter,'recursive limiter')
      error=max(maxval(abs(b-sb)/max(1._dp,abs(sb))),maxval(abs(w-sw)/max(1._dp,abs(sw))), &
        maxval(abs(ch-sch)/max(1._dp,abs(sch))),maxval(abs(ext-sext)/max(1._dp,abs(sext))))
      call require(error<2.e-12_dp,'recursive serial ledgers')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,final_inventory,MPI_COMM_WORLD,good)
      call require(good,'recursive final inventory')
      if(allocated(aux)) final_inventory=final_inventory-aux
      error=maxval(abs(final_inventory-initial_inventory-b-w-ch-ext)/max(1._dp,abs(initial_inventory)))
      if(rank==0) print *, 'Recursive SDC layout, residual: ',layout,error
      call require(error<3.e-14_dp,'recursive closure')
      coefficient=-1;mpi_source_calls=0
      if(allocated(aux)) then
        coefficient=.1_dp;source_failure=1;aux=41
      end if
      reactor_calls=0;reactor_failure=test_mode=="reactor"
      if(reactor_failure) coefficient=.1_dp
      old_heat=heat
      old=local;ot=lt;old_ir=ir;b=17;w=19;ch=23;ext=29;limiter=31;counts=37
      failed_order=requested_order
      if(rank==nranks-1) failed_order=merge(4,2,requested_order==2)
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
        1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,level_steps=counts,redistribution=method, &
        state_redist_max_order=failed_order,heat_release=heat)
      call require(.not.good.and.reactor_calls==0,'recursive rank-order mismatch rejects before stages')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==old%fields(p)%values),'rank-order recursive state rollback')
        call require(all(lt%fields(p)%values==ot%fields(p)%values),'rank-order recursive temperature rollback')
        call require(all(ir%fields(p)%values==old_ir%fields(p)%values),'rank-order recursive reaction rollback')
      enddo
      call require(all(counts==37).and.all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.limiter==31, &
        'rank-order recursive counters and ledger rollback')
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
        1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,difmag=coefficient,level_steps=counts,redistribution=method, &
        auxiliary_source=mpi_provider,auxiliary_change=aux,route_auxiliary=test_mode=='aux_route', &
        evaluate_reaction=reactor_adapter, &
      state_redist_max_order=requested_order,heat_release=heat)
      reactor_failure=.false.
      call require(.not.good,'recursive rejected stage')
      if(allocated(aux)) then
        call require(mpi_source_calls==21.and.all(aux==41),'last descendant source failure and ledger rollback')
      end if
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==old%fields(p)%values),'recursive state rollback')
        call require(all(lt%fields(p)%values==ot%fields(p)%values),'recursive T rollback')
        call require(all(ir%fields(p)%values==old_ir%fields(p)%values),'recursive IR rollback')
        call require(all(heat%fields(p)%values==old_heat%fields(p)%values),'recursive heat rollback')
      end do
      call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.limiter==31.and. &
        all(counts==37),'recursive diagnostics rollback')
      deallocate(local%fields,lt%fields,ir%fields,heat%fields)
    end do
  end subroutine check_sdc_recursive
  subroutine check_sdc_level(initial,geometry,extended,owners,old,ot)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(in) :: old,ot
    type(eb_reactive_tree_3d) :: reference,ending
    type(amr_tree_state_3d) :: state,temp,ir,serial_ir,finish,source0,source1,local0,local1,op,np,sop,snp
    type(eb_integrated_flux_3d), allocatable :: flux(:),serial_flux(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),sb(nv),sw(nv),sch(nv),sext(nv),theta,stheta,error
    real(dp) :: force(3)
    integer :: p,attempt
    character(len=16) :: faces(6)
    logical :: good
    reference=initial;ending=initial;state=old;temp=ot;finish=old
    faces='outflow';force=[2.e4_dp,-3.e4_dp,4.e4_dp]
    allocate(serial_ir%fields(count),ir%fields(count),source0%fields(count),source1%fields(count), &
      local0%fields(count),local1%fields(count))
    do p=1,count
      ending%nodes(p)%state=1.01_dp*initial%nodes(p)%state
      allocate(serial_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
      if(topology%nodes(p)%level/=2) serial_ir%fields(p)%values=7
      allocate(source0%fields(p)%values(nv,4,4,4),source=0._dp)
      allocate(source1%fields(p)%values(nv,4,4,4),source=0._dp)
      source0%fields(p)%values(2,:,:,:)=1.e4_dp;source1%fields(p)%values(2,:,:,:)=2.e4_dp
      if(owners(p)/=rank) cycle
      finish%fields(p)%values=ending%nodes(p)%state
      ir%fields(p)%values=serial_ir%fields(p)%values
      local0%fields(p)%values=source0%fields(p)%values
      local1%fields(p)%values=source1%fields(p)%values
    end do
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,stheta,good,serial_ir,target_depth=2,ancestor_start=initial,ancestor_end=ending, &
      fractions_start=[.25_dp,.25_dp],fractions_end=[.75_dp,.75_dp],source_start=source0,source_end=source1, &
      integrated_flux=serial_flux,old_predictor=sop,new_predictor=snp, &
      state_redist_max_order=requested_order)
    call require(good,'serial selected coupled SDC')
    do attempt=1,2
      if(attempt==1) then
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temp,extended, &
          1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
          b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,target_depth=2,ancestor_start=old,ancestor_end=finish, &
          temperature_start=ot,temperature_end=ot,fractions_start=[.25_dp,.25_dp],fractions_end=[.75_dp,.75_dp], &
          source_start=local0,source_end=local1,integrated_flux=flux,old_predictor=op,new_predictor=np, &
      state_redist_max_order=requested_order)
        call require(good,'MPI selected coupled SDC')
      else
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temp,extended, &
          1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
          b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,difmag=-1._dp,target_depth=2, &
          integrated_flux=flux,old_predictor=op,new_predictor=np, &
      state_redist_max_order=requested_order)
        call require(.not.good,'selected loop rejects after old transport')
      end if
      call require(theta==stheta,'selected loop limiter retained')
      error=max(maxval(abs(b-sb)/max(1._dp,abs(sb))),maxval(abs(w-sw)/max(1._dp,abs(sw))), &
        maxval(abs(ch-sch)/max(1._dp,abs(sch))),maxval(abs(ext-sext)/max(1._dp,abs(sext))))
      call require(error<2.e-12_dp,'selected loop inventories match or roll back')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(state%fields(p)%values==reference%nodes(p)%state),'selected loop exact state')
        call require(all(temp%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'selected loop exact T')
        call require(all(ir%fields(p)%values==serial_ir%fields(p)%values),'selected loop exact IR and inactive retention')
        call require(all(op%fields(p)%values==sop%fields(p)%values),'selected loop exact old predictor')
        call require(all(np%fields(p)%values==snp%fields(p)%values),'selected loop exact new predictor')
        if(topology%nodes(p)%level/=2) then
          call require(.not.allocated(flux(p)%x),'selected loop inactive flux absent')
          cycle
        end if
        call require(all(flux(p)%x==serial_flux(p)%x).and.all(flux(p)%y==serial_flux(p)%y).and. &
          all(flux(p)%z==serial_flux(p)%z),'selected loop exact flux')
      end do
    end do
  end subroutine check_sdc_level
  subroutine check_sdc_reaction()
    type(amr_tree_state_3d) :: rates,ref_state,ref_temp,ref_ir,local,lt,lr,ir,old,ot,kept,kt,ki,ref_held,held, &
      policy_state,policy_temp,policy_ir,policy_held
    type(amr_tree_state_3d) :: old_held,owner_old_held,input_state,input_temp,input_ir,input_held, &
      trial_state,trial_temp,trial_ir,trial_held
    type(eb_geometry_3d) :: geometry(count)
    real(dp) :: tolerance,ref_correction,mpi_correction
    integer :: p,layout,owners(count)
    logical :: good
    allocate(rates%fields(count),ref_ir%fields(count))
    do p=1,count
      geometry(p)=saved%nodes(p)%geometry
      allocate(rates%fields(p)%values(nv,4,4,4),source=0._dp)
      allocate(ref_ir%fields(p)%values(nv-4,4,4,4),source=7._dp)
      rates%fields(p)%values(2,:,:,:)=2.e4_dp
      rates%fields(p)%values(5,:,:,:)=2.e4_dp*saved%nodes(p)%state(2,:,:,:)/saved%nodes(p)%state(1,:,:,:)
    end do
    call react_sdc_eb_tree_3d(species,reactions,saved,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp, &
      ref_state,ref_temp,ref_ir,good)
    call require(good,'serial SDC reaction reference')
    allocate(ref_held%fields(count))
    do p=1,count
      allocate(ref_held%fields(p)%values(1,4,4,4),source=23._dp)
    enddo
    call react_sdc_eb_tree_3d(species,reactions,saved,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp, &
      ref_state,ref_temp,ref_ir,good,internal_energy=ref_held)
    call require(good,'serial SDC retained reactor oracle')
    policy_state=ref_state;policy_temp=ref_temp;policy_ir=ref_ir;policy_held=ref_held;ref_correction=7
    call react_sdc_eb_tree_3d(species,reactions,saved,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp, &
      policy_state,policy_temp,policy_ir,good,internal_energy=policy_held,energy_correction=ref_correction)
    call require(good,'serial retained EOS policy oracle')
    old_held=ref_held
    do p=1,count
      old_held%fields(p)%values(1,:,:,:)=saved%nodes(p)%state(5,:,:,:)- &
        .5_dp*sum(saved%nodes(p)%state(2:4,:,:,:)**2,dim=1)/ &
        max(saved%nodes(p)%state(1,:,:,:),tiny(1._dp))+16384._dp
    enddo
    input_ir=ref_ir;input_held=ref_held
    call react_sdc_eb_tree_3d(species,reactions,saved,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp, &
      input_state,input_temp,input_ir,good,internal_energy=input_held,old_internal_energy=old_held)
    call require(good,'serial packed old UEINT reaction oracle')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(old%fields(count),ot%fields(count),lr%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        old%fields(p)%values=saved%nodes(p)%state
        allocate(ot%fields(p)%values(1,4,4,4))
        ot%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
        lr%fields(p)%values=rates%fields(p)%values
        allocate(ir%fields(p)%values(nv-4,4,4,4),source=7._dp)
      end do
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good)
      call require(good,'MPI SDC reaction')
      allocate(held%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        allocate(held%fields(p)%values(1,4,4,4),source=23._dp)
      enddo
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good,internal_energy=held)
      call require(good,'MPI SDC retained reactor endpoint')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(held%fields(p)%values==ref_held%fields(p)%values), &
          'MPI retained reactor exact serial endpoint')
      enddo
      if(rank==0) then
        call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
          1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good,internal_energy=held)
      else
        call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
          1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good)
      endif
      call require(.not.good,'MPI retained presence disagreement rejected')
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'SDC owner-local allocation')
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==ref_state%fields(p)%values).and. &
          all(lt%fields(p)%values==ref_temp%fields(p)%values).and.all(ir%fields(p)%values==ref_ir%fields(p)%values), &
          'MPI SDC exact serial state temperature and reaction')
      end do
      mpi_correction=7
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good, &
        internal_energy=held,energy_correction=mpi_correction)
      call require(good.and.mpi_correction==ref_correction,'MPI retained EOS correction ledger')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==policy_state%fields(p)%values).and. &
          all(lt%fields(p)%values==policy_temp%fields(p)%values).and. &
          all(ir%fields(p)%values==policy_ir%fields(p)%values).and. &
          all(held%fields(p)%values==policy_held%fields(p)%values), &
          'MPI retained EOS exact serial endpoint')
      enddo
      if(rank==0) then
        call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
          1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good, &
          internal_energy=held,energy_correction=mpi_correction)
      else
        call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
          1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good,internal_energy=held)
      endif
      call require(.not.good.and.mpi_correction==ref_correction,'MPI retained policy presence disagreement rollback')
      allocate(owner_old_held%fields(count))
      do p=1,count
        if(owners(p)==rank) owner_old_held%fields(p)=old_held%fields(p)
      enddo
      trial_state=local;trial_temp=lt;trial_ir=ir;trial_held=held
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,trial_state,trial_temp,trial_ir,MPI_COMM_WORLD,good, &
        internal_energy=trial_held,old_internal_energy=owner_old_held)
      call require(good,'MPI packed old UEINT native reaction')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(trial_state%fields(p)%values==input_state%fields(p)%values).and. &
          all(trial_temp%fields(p)%values==input_temp%fields(p)%values).and. &
          all(trial_ir%fields(p)%values==input_ir%fields(p)%values).and. &
          all(trial_held%fields(p)%values==input_held%fields(p)%values), &
          'MPI packed old UEINT exact serial state and endpoint')
        call require(all(owner_old_held%fields(p)%values==old_held%fields(p)%values), &
          'MPI packed old UEINT input unchanged')
      enddo
      if(rank==0) then
        call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
          1.e-9_dp,2.e-7_dp,1.e-12_dp,trial_state,trial_temp,trial_ir,MPI_COMM_WORLD,good, &
          internal_energy=trial_held,old_internal_energy=owner_old_held)
      else
        call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
          1.e-9_dp,2.e-7_dp,1.e-12_dp,trial_state,trial_temp,trial_ir,MPI_COMM_WORLD,good,internal_energy=trial_held)
      endif
      call require(.not.good,'MPI packed old UEINT presence mismatch aborts collectively')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(trial_state%fields(p)%values==input_state%fields(p)%values).and. &
          all(trial_temp%fields(p)%values==input_temp%fields(p)%values).and. &
          all(trial_ir%fields(p)%values==input_ir%fields(p)%values).and. &
          all(trial_held%fields(p)%values==input_held%fields(p)%values), &
          'MPI packed old UEINT presence mismatch keeps all outputs')
      enddo
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=-17;lt%fields(p)%values=-19;ir%fields(p)%values=7
      end do
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good,target_depth=2)
      call require(good,'MPI selected level SDC reaction')
      do p=1,count
        if(owners(p)/=rank) cycle
        if(topology%nodes(p)%level==2) then
          call require(all(local%fields(p)%values==ref_state%fields(p)%values).and. &
            all(ir%fields(p)%values==ref_ir%fields(p)%values),'MPI selected SDC exact reaction')
        else
          call require(all(local%fields(p)%values==-17).and.all(lt%fields(p)%values==-19).and. &
            all(ir%fields(p)%values==7),'MPI SDC inactive level retained')
        end if
      end do
      kept=local;kt=lt;ki=ir
      if(owners(count)==rank) old%fields(count)%values(1,4,4,4)=-1
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good)
      call require(.not.good,'MPI failing owner rejects globally')
      if(owners(count)==rank) old%fields(count)%values=saved%nodes(count)%state
      tolerance=2.e-7_dp
      if(rank==nranks-1) tolerance=-1
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,tolerance,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good)
      call require(.not.good,'MPI invalid control including idle rank rejects')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==kept%fields(p)%values).and. &
          all(lt%fields(p)%values==kt%fields(p)%values).and.all(ir%fields(p)%values==ki%fields(p)%values), &
          'MPI SDC collective rollback')
      end do
      deallocate(old%fields,ot%fields,lr%fields,ir%fields,held%fields,owner_old_held%fields)
    end do
  end subroutine check_sdc_reaction

  subroutine check_mpi_full()
    type(amr_tree_state_3d) :: local,lt,initial,it
    type(eb_geometry_3d) :: geometry(count)
    real(dp) :: b(nv),w(nv),chem(nv),limiter,inventory(nv),error
    integer :: owners(count),p
    logical :: good,enabled
    owners=[(mod(p-1,nranks),p=1,count)]
    if (integrator_index==2) owners=nranks-1
    allocate(local%fields(count),lt%fields(count))
    do p=1,count
      geometry(p)=saved%nodes(p)%geometry
      if (owners(p)/=rank) cycle
      local%fields(p)%values=saved%nodes(p)%state
      allocate(lt%fields(p)%values(1,4,4,4))
      lt%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
    end do
    initial=local;it=lt
    call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,'hllc',1.e-8_dp, &
      chemistry_on,2.e-7_dp,1.e-12_dp,transport_on,.true.,.true.,.true.,.true., &
      local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,integrators(integrator_index))
    call require(good,'MPI full step')
    do p=1,count
      call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'owner-local full U')
      if (owners(p)/=rank) cycle
      call require(all(local%fields(p)%values==tree%nodes(p)%state),'exact serial full U')
      call require(all(lt%fields(p)%values(1,:,:,:)==tree%nodes(p)%temperature),'exact serial full T')
    end do
    call require(all(b==boundary).and.limiter==theta,'exact boundary and limiter')
    call require(maxval(abs(w-wall)/max(1._dp,abs(before)))<3.e-14_dp,'wall parity')
    call require(maxval(abs(chem-chemical)/max(1._dp,abs(before)))<3.e-14_dp,'chemical ledger parity')
    call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,inventory,MPI_COMM_WORLD,good)
    call require(good,'MPI inventory')
    error=maxval(abs(inventory-before-b-w-chem)/max(1._dp,abs(before)))
    call require(error<3.e-14_dp,'MPI full conservation')
    if (rank==0) print *, 'MPI integrator, mode, residual: ',integrators(integrator_index),mode,error
    if (integrator_index==1.and.mode==3) then
      local=initial;lt=it;b=17._dp;w=19._dp;chem=23._dp;limiter=29._dp
      call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,'invalid',1.e-8_dp, &
        .true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
        local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit')
      call require(.not.good,'late hydro rejection after reaction and transport')
      do p=1,count
        if (owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==initial%fields(p)%values),'full U rollback')
        call require(all(lt%fields(p)%values==it%fields(p)%values),'full T rollback')
      end do
      call require(all(b==17._dp).and.all(w==19._dp).and.all(chem==23._dp).and.limiter==29._dp,'full ledger rollback')
      enabled=rank==0
      call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,'hllc',1.e-8_dp, &
        enabled,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
        local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit')
      call require(.not.good,'different chemistry flags rejected')
    end if
  end subroutine check_mpi_full
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine require
  subroutine require_restricted_reactions(layout,field)
    type(eb_reactive_tree_3d), intent(in) :: layout
    type(amr_tree_state_3d), intent(in) :: field
    integer :: child,parent,r,lo(3),hi(3),i,j,k,first(3),last(3),c
    real(dp) :: expected,volume,scale
    do child=2,size(layout%nodes)
      parent=layout%topology%nodes(child)%parent
      associate(box=>layout%topology%nodes(child)%patch)
        r=box%refinement_ratio
        lo=[box%coarse_i_lower,box%coarse_j_lower,box%coarse_k_lower]
        hi=[box%coarse_i_upper,box%coarse_j_upper,box%coarse_k_upper]
      end associate
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        first=([i,j,k]-lo)*r+1;last=first+r-1
        associate(v=>layout%nodes(child)%geometry%volume_fraction(first(1):last(1),first(2):last(2),first(3):last(3)))
          volume=sum(v)
          do c=1,size(field%fields(child)%values,1)
            expected=0
            if(volume>0) expected=sum(v*field%fields(child)%values(c, &
              first(1):last(1),first(2):last(2),first(3):last(3)))/volume
            scale=max(tiny(1._dp),abs(expected))
            call require(abs(field%fields(parent)%values(c,i,j,k)-expected)<=32*epsilon(1._dp)*scale, &
              'synchronized reaction source must equal physical fine-volume average')
          enddo
        end associate
      enddo
      enddo
      enddo
    enddo
  end subroutine require_restricted_reactions
end program test_mpi_eb_tree_full_3d
