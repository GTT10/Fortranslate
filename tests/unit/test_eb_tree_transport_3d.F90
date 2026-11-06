program test_eb_tree_transport_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use eb_reactive_transport_3d_mod, only: compute_sdc_eb_transport_halo_3d
  use eb_wall_heat_source_3d_mod, only: compute_eb_wall_transport_source_3d
  use eb_reactive_redistribution_3d_mod, only: state_redist_register_3d,compute_state_redist_level_stage_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    reactive_eb_tree_timestep_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_transport_3d, &
    compute_sdc_eb_tree_transport_halos_3d,fill_reactive_eb_tree_halo_3d,build_eb_tree_extended_geometry_3d, &
    compute_sdc_eb_tree_transport_3d,synchronize_reactive_eb_tree_flux_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use eb_reactive_patch_tree_3d_mod, only: advance_reactive_eb_tree_sdc_3d,compute_mol_eb_tree_hydro_3d, &
    advance_reactive_eb_tree_sdc_subcycle_3d,compute_sdc_eb_tree_hydro_3d, &
    compute_godunov_state_redist_eb_tree_hydro_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use elementary_kinetics_mod, only: elementary_reaction
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:),boundary(:),primitive(:)
  real(dp) :: residual,theta,speed,sdc_dt,split_dt,common_dt
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k
  character(len=16) :: periodic_faces(6)
  real(dp) :: wall_t(6),base_dt,thermal_dt,plane_cut
  integer :: bound_mode
  integer :: observed_endpoint=0
  logical :: tightened=.false.
  logical :: ok
  character(len=32) :: test_mode
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),primitive(reactive_nprim(size(species))))
  call get_command_argument(1,test_mode)
  plane_cut=.31_dp
  if((test_mode=='multi_sdc'.or.test_mode=='inert_iters'.or.test_mode=='godunov_halo')) plane_cut=.39_dp
  do axis=1,3
  do count=4,5
  if((test_mode=='godunov_halo'.or.test_mode=='spatial_halo'.or. &
    test_mode=='prepared_sdc'.or. &
    (test_mode=='multi_sdc'.or.test_mode=='inert_iters').or. &
    test_mode=='recursive_srd'.or.test_mode=='common_srd').and.(axis/=3.or.count/=5)) cycle
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,2,3,2,boxes(1),ok)
  call require(ok,'child left')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,2,3,2,boxes(2),ok)
  call require(ok,'child right')
  call initialize_amr_patch_3d(4,4,4,3,4,3,4,2,3,2,boxes(3),ok)
  call require(ok,'grandchild left')
  call initialize_amr_patch_3d(4,4,4,1,2,3,4,2,3,2,boxes(4),ok)
  call require(ok,'grandchild right')
  if((test_mode=='multi_sdc'.or.test_mode=='inert_iters'.or.test_mode=='godunov_halo')) then
    ! Put the coarse cut donor immediately below the refined regular cells.
    ! The original boxes covered both members of its vertical neighborhood.
    call initialize_amr_patch_3d(4,4,4,1,2,1,2,3,4,2,boxes(1),ok)
    call require(ok,'coarse-interface child left')
    call initialize_amr_patch_3d(4,4,4,3,4,1,2,3,4,2,boxes(2),ok)
    call require(ok,'coarse-interface child right')
  endif
  if (count==5) then
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3],boxes,topology,ok)
  else
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],boxes(:3),topology,ok)
  end if
  call require(ok,'topology')
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,axes(axis),plane_cut,finest,ok)
  call require(ok,'finest geometry')
  call coarsen_eb_geometry_3d(finest,2,middle,ok)
  call require(ok,'middle geometry')
  call coarsen_eb_geometry_3d(middle,2,coarse,ok)
  call require(ok,'coarse geometry')
  allocate(u(nv,4,4,4),t(4,4,4))
  config%nx=4;config%ny=4;config%nz=4
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
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
  if((test_mode=='godunov_halo'.or.test_mode=='spatial_halo'.or. &
    test_mode=='prepared_sdc'.or. &
    (test_mode=='multi_sdc'.or.test_mode=='inert_iters').or.test_mode=='recursive_srd'.or.test_mode=='common_srd')) then
    call check_prepared_halos()
    stop
  endif
  call reactive_eb_tree_timestep_3d(species,transport,tree,.5_dp,.35_dp,.true.,.true.,.true.,.true., &
    common_dt,ok,subcycling=.false.)
  call require(ok,'common timestep bound')
  call reactive_eb_tree_timestep_3d(species,transport,tree,.5_dp,.35_dp,.true.,.true.,.true.,.true., &
    split_dt,ok,subcycling=.true.)
  call require(ok,'split recursive timestep bound')
  call reactive_eb_tree_timestep_3d(species,transport,tree,.5_dp,.35_dp,.true.,.true.,.true.,.true., &
    sdc_dt,ok,subcycling=.true.,coupled_subcycling=.true.)
  call require(ok.and.sdc_dt>=common_dt.and.sdc_dt<=split_dt,'SDC uses coupled transport substep scaling')

  boundary=17._dp;theta=19._dp
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok)
  call require(ok,'tree transport step')
  call require(all(boundary==0._dp),'closed physical diffusion faces')
  call require(theta>=0._dp.and.theta<=1._dp,'limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'final inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'axis, nodes, inventory residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'transport conservation')
  call require(any(tree%nodes(count)%state/=saved%nodes(count)%state),'nontrivial deep transport update')
  if(axis==3.and.count==5) call check_sdc_transport()
  periodic_faces='no_slip_wall'
  tree=saved
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok,boundary_faces=periodic_faces)
  call require(ok,'wall transport update')
  call require(boundary(1)==0.and.boundary(5)==0.and.all(boundary(6:)==0),'impermeable adiabatic wall flux')
  call require(theta>=0.and.theta<=1,'periodic limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'wall inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'wall axis, nodes, residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'wall composite transport conservation')
  periodic_faces='isothermal_wall'
  wall_t=0
  do node=1,count
    wall_t=max(wall_t,maxval(saved%nodes(node)%temperature)+100)
  end do
  tree=saved
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok,boundary_faces=periodic_faces,wall_temperature=wall_t)
  call require(ok,'isothermal transport update')
  call require(boundary(1)==0.and.boundary(5)>0.and.all(boundary(6:)==0),'impermeable heated wall flux')
  call require(theta>=0.and.theta<=1,'periodic limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'isothermal inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'isothermal axis, nodes, residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'isothermal composite transport conservation')
  do bound_mode=0,1
    call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,1.e-12_dp,.true., &
      .false.,.true.,.false.,base_dt,ok,subcycling=bound_mode==1)
    call require(ok,'baseline tree timestep')
    call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,1.e-12_dp,.true., &
      .false.,.true.,.false.,thermal_dt,ok,subcycling=bound_mode==1, &
      boundary_faces=periodic_faces,wall_temperature=wall_t)
    call require(ok.and.thermal_dt>0.and.thermal_dt<=base_dt,'thermal timestep restriction')
    tightened=tightened.or.thermal_dt<base_dt
    print *, 'thermal timestep axis, nodes, mode, baseline, wall: ',axes(axis),count,bound_mode,base_dt,thermal_dt
    thermal_dt=17
    call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,1.e-12_dp,.true., &
      .false.,.true.,.false.,thermal_dt,ok,subcycling=bound_mode==1,boundary_faces=periodic_faces)
    call require(.not.ok.and.thermal_dt==17,'missing wall temperature preserves timestep')
  end do

  periodic_faces='periodic';periodic_faces(2*axis-1:2*axis)='outflow'
  tree=saved
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok,boundary_faces=periodic_faces)
  call require(ok,'periodic transport update')
  call require(theta>=0.and.theta<=1,'periodic limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'periodic inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'periodic axis, nodes, residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'periodic composite transport conservation')
  tree=saved;boundary=17;theta=19
  periodic_faces(2*mod(axis,3)+1)='outflow'
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok,boundary_faces=periodic_faces)
  call require(.not.ok.and.all(boundary==17).and.theta==19,'unpaired periodic transport rollback')
  call unchanged_tree()
  tree=saved;boundary=17._dp;theta=19._dp
  call advance_reactive_eb_tree_transport_3d(species,transport,-1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok)
  call require(.not.ok,'invalid interval rejected')
  call unchanged_tree()
  call require(all(boundary==17._dp).and.theta==19._dp,'diagnostic rollback')
  tree%nodes(count)%state(1,4,4,4)=-1._dp;saved=tree
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok)
  call require(.not.ok,'invalid deep state rejected')
  call unchanged_tree()
  call require(all(boundary==17._dp).and.theta==19._dp,'deep failure diagnostic rollback')
  deallocate(u,t)
  end do
  end do
  call require(tightened,'wall heat restricts an accepted tree timestep')
  print *, 'test_eb_tree_transport_3d: PASS'
contains
  subroutine check_recursive_registers()
    type(amr_patch_tree_3d) :: nested
    type(amr_patch_3d) :: patches(2)
    type(eb_geometry_3d) :: fg,mg,cg,padded,physical(3)
    type(eb_geometry_3d), allocatable :: wide(:)
    type(eb_reactive_tree_3d) :: candidate,kept
    type(elementary_reaction) :: reactions(0)
    type(state_redist_register_3d), allocatable :: registers(:)
    real(dp) :: root_state(nv,4,4,4),root_temp(4,4,4),q(reactive_nprim(size(species))),pos(3),tt,cs, &
      b(nv),w(nv),ch(nv),sc(nv),aux(nv),lim,initial(nv),final(nv),error
    integer(int64) :: counts(3)
    integer :: p,x,y,z,c,step,n(3)
    character(len=32) :: faces(6)
    logical :: good
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
          if((test_mode=='multi_sdc'.or.test_mode=='inert_iters'.or.test_mode=='godunov_halo')) then
            q(1)=1+.01_dp*pos(1)**2+.02_dp*pos(3)**2
            q(5)=250000+1000*pos(2)**2+3000*pos(3)**2
          endif
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
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'z',.43_dp,padded,good)
    call require(good,'recursive ledger padded geometry')
    call build_eb_tree_extended_geometry_3d(nested,physical,padded,wide,good)
    call require(good,'recursive ledger extended metrics')
    faces='outflow'
    do step=1,2
      call integrate_reactive_eb_tree_3d(candidate,initial,good)
      call require(good,'recursive ledger initial inventory')
      if(test_mode=='common_srd') then
        if(step==1) then
          call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,1,.false., &
            [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
            b,w,ch,sc,lim,good,redistribution='state_redist',mol_mode=.true.,state_redist_max_order=2, &
            embedded_wall_temperature=700._dp,embedded_noslip=.true.)
        else
          call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,1,.false., &
            [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
            b,w,ch,sc,lim,good,redistribution='state_redist',mol_mode=.true.,state_redist_max_order=2, &
            embedded_wall_temperature=700._dp,embedded_noslip=.true.,state_redist_registers=registers)
        endif
        call require(good,'common MOL interface interval accepted')
        if(step==2) then
          do p=1,3
            if(allocated(registers(p)%coarse_delta)) call require(all(registers(p)%coarse_delta==0),'common coarse consumed')
            if(allocated(registers(p)%fine_delta)) call require(all(registers(p)%fine_delta==0),'common fine consumed')
          enddo
        endif
      else
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,1,.false., &
        [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
        b,w,ch,sc,lim,good,level_steps=counts,redistribution='state_redist',mol_mode=.true., &
        state_redist_max_order=2,embedded_wall_temperature=700._dp,embedded_noslip=.true.)
      call require(good,'recursive ledger accepted three-level MOL interval')
      call require(all(counts==[1_int64,2_int64,4_int64]),'recursive ledger substep counts')
      endif
      call integrate_reactive_eb_tree_3d(candidate,final,good)
      call require(good,'recursive ledger final inventory')
      error=maxval(abs(final-initial-b-w-ch-sc)/max(1._dp,abs(initial)))
      print *, 'recursive register step, composite budget error: ',step,error
      call require(error<2.e-12_dp,'recursive ledger composite gas and body wall budget')
    enddo
    kept=candidate;b=17;w=19;ch=23;sc=29;aux=31;lim=37;counts=41
    if(test_mode=='common_srd') then
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,1,.false., &
        [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
        b,w,ch,sc,lim,good,redistribution='state_redist',mol_mode=.true.,state_redist_max_order=2, &
        auxiliary_source=reject_common_endpoint,auxiliary_change=aux)
    else
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,1,.false., &
      [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces, &
      b,w,ch,sc,lim,good,level_steps=counts,redistribution='state_redist',mol_mode=.true.,state_redist_max_order=2, &
      auxiliary_source=reject_recursive_endpoint,auxiliary_change=aux)
    endif
    call require(.not.good,'recursive ledger late descendant rejection')
    call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(sc==29).and.all(aux==31).and.lim==37.and. &
      all(counts==41),'recursive ledger full diagnostic and count rollback')
    do p=1,3
      call require(all(candidate%nodes(p)%state==kept%nodes(p)%state).and. &
        all(candidate%nodes(p)%temperature==kept%nodes(p)%temperature),'recursive ledger outer gas rollback')
    enddo
  end subroutine check_recursive_registers
  subroutine reject_common_endpoint(snapshot,interval,stage,depth,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: snapshot
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p
    good=stage==0
    if(.not.good) return
    do p=1,size(rate%fields)
      if(allocated(rate%fields(p)%values)) rate%fields(p)%values=0
    enddo
  end subroutine reject_common_endpoint
  subroutine reject_recursive_endpoint(snapshot,interval,stage,depth,time,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: snapshot
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p
    good=stage==0.or.depth/=2
    if(.not.good) return
    do p=1,size(snapshot%nodes)
      rate%fields(p)%values=0
    enddo
  end subroutine reject_recursive_endpoint
  subroutine native_post_density_control(snapshot,rate,depth)
    type(eb_reactive_tree_3d), intent(in) :: snapshot
    type(amr_tree_state_3d), intent(inout) :: rate
    integer, intent(in) :: depth
    integer :: p,i,j,k,n
    ! EB.cpp posts density after constructing the redistribution registers.
    do p=1,size(snapshot%nodes)
      if(topology%nodes(p)%level/=depth) cycle
      do k=1,size(rate%fields(p)%values,4)
      do j=1,size(rate%fields(p)%values,3)
      do i=1,size(rate%fields(p)%values,2)
        if(snapshot%nodes(p)%geometry%volume_fraction(i,j,k)<=0) then
          rate%fields(p)%values(:,i,j,k)=0
        else
          rate%fields(p)%values(1,i,j,k)=0
          do n=6,size(rate%fields(p)%values,1)
            rate%fields(p)%values(1,i,j,k)=rate%fields(p)%values(1,i,j,k)+rate%fields(p)%values(n,i,j,k)
          enddo
        endif
      enddo
      enddo
      enddo
    enddo
  end subroutine native_post_density_control
  subroutine check_nonmol_sdc(fixture,narrow,wide)
    type(eb_reactive_tree_3d), intent(in) :: fixture
    type(eb_geometry_3d), intent(in) :: narrow(:),wide(:)
    type(eb_reactive_tree_3d) :: candidate,predicted,kept
    type(elementary_reaction) :: reactions(0)
    type(amr_tree_state_3d) :: d0,d1,d2,h1,h2,forcing,halo,old
    type(eb_integrated_flux_3d), allocatable :: fd0(:),fd1(:),fd2(:),fh1(:),fh2(:),actual_flux(:),kept_flux(:)
    type(state_redist_register_3d), allocatable :: rd0(:),rd1(:),rd2(:),rh1(:),rh2(:),actual_regs(:),kept_regs(:)
    real(dp), parameter :: dt=1.e-10_dp
    real(dp) :: limiter,bc(nv),wc(nv),ch(nv),sc(nv),q(reactive_nprim(size(species))),tt,cs
    real(dp) :: expected(nv,4,4,4)
    integer :: p,i,j,k
    logical :: good
    character(len=32) :: faces(6)
    faces='outflow'
    call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,narrow,dt, &
      [.true.,.true.,.true.,.true.],d0,fd0,halo,old,limiter,good,boundary_faces=faces,target_depth=0, &
      registers=rd0,max_order=2)
    call require(good,'Non-MOL independent diffusion stage')
    call native_post_density_control(fixture,d0,0)
    forcing=d0
    call compute_godunov_state_redist_eb_tree_hydro_3d(species,fixture,forcing,wide,dt,h1,fh1,good, &
      boundary_faces=faces,target_depth=0,registers=rh1,max_order=2)
    call require(good,'Non-MOL independent hydro stage')
    predicted=fixture
    predicted%nodes(1)%state=fixture%nodes(1)%state+dt*(d0%fields(1)%values+h1%fields(1)%values)
    do k=1,4
    do j=1,4
    do i=1,4
      if(fixture%nodes(1)%geometry%volume_fraction(i,j,k)<=0) cycle
      call reactive_conserved_to_primitive(species,predicted%nodes(1)%state(:,i,j,k), &
        predicted%nodes(1)%temperature(i,j,k),q,tt,cs,good)
      call require(good,'Non-MOL independent predictor EOS')
      predicted%nodes(1)%temperature(i,j,k)=tt
    enddo
    enddo
    enddo
    call compute_sdc_eb_tree_transport_halos_3d(species,transport,predicted,narrow,dt, &
      [.true.,.true.,.true.,.true.],d1,fd1,halo,old,limiter,good,boundary_faces=faces,target_depth=0, &
      registers=rd1,max_order=2)
    call require(good,'Non-MOL independent diffusion stage')
    call native_post_density_control(fixture,d1,0)
    forcing=d0
    forcing%fields(1)%values=.5_dp*(d0%fields(1)%values+d1%fields(1)%values)
    call compute_godunov_state_redist_eb_tree_hydro_3d(species,fixture,forcing,wide,dt,h2,fh2,good, &
      boundary_faces=faces,target_depth=0,registers=rh2,max_order=2)
    call require(good,'Non-MOL independent hydro stage')
    predicted=fixture
    predicted%nodes(1)%state=fixture%nodes(1)%state+dt*(forcing%fields(1)%values+h2%fields(1)%values)
    do k=1,4
    do j=1,4
    do i=1,4
      if(fixture%nodes(1)%geometry%volume_fraction(i,j,k)<=0) cycle
      call reactive_conserved_to_primitive(species,predicted%nodes(1)%state(:,i,j,k), &
        predicted%nodes(1)%temperature(i,j,k),q,tt,cs,good)
      call require(good,'Non-MOL independent predictor EOS')
      predicted%nodes(1)%temperature(i,j,k)=tt
    enddo
    enddo
    enddo
    call compute_sdc_eb_tree_transport_halos_3d(species,transport,predicted,narrow,dt, &
      [.true.,.true.,.true.,.true.],d2,fd2,halo,old,limiter,good,boundary_faces=faces,target_depth=0, &
      registers=rd2,max_order=2)
    call require(good,'Non-MOL independent diffusion stage')
    call native_post_density_control(fixture,d2,0)
    expected=fixture%nodes(1)%state+dt*(h2%fields(1)%values+.5_dp*(d0%fields(1)%values+d2%fields(1)%values))
    candidate=fixture
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,2,.false.,[.true.,.true.,.true.,.true.], &
      [0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces,bc,wc,ch,sc,limiter,good, &
      target_depth=0,integrated_flux=actual_flux,redistribution='state_redist',mol_mode=.false., &
      state_redist_max_order=2,state_redist_registers=actual_regs)
    call require(good,'Prepared non-MOL two-iteration SDC')
    call require(maxval(abs(candidate%nodes(1)%state-expected)/max(1._dp,abs(expected)))<1.e-13_dp, &
      'Non-MOL independently constructed two-iteration gas')
    call require(maxval(abs(actual_flux(1)%x-(fh2(1)%x+.5_dp*(fd0(1)%x+fd2(1)%x))))+ &
      maxval(abs(actual_flux(1)%y-(fh2(1)%y+.5_dp*(fd0(1)%y+fd2(1)%y))))+ &
      maxval(abs(actual_flux(1)%z-(fh2(1)%z+.5_dp*(fd0(1)%z+fd2(1)%z))))<1.e-13_dp, &
      'Non-MOL final hydro and old/latest diffusion faces')
    call require(maxval(abs(rd1(1)%coarse_delta))+maxval(abs(rh1(1)%coarse_delta))>0, &
      'Non-MOL provisional coarse contributions nonzero')
    call require(all(actual_regs(1)%coarse_delta== &
      (rd0(1)%coarse_delta+(rd2(1)%coarse_delta+rd1(1)%coarse_delta))+ &
      (rh2(1)%coarse_delta+rh1(1)%coarse_delta)),'Non-MOL all diffusion/hydro coarse evaluations retained')
    do p=2,count
      call require(all(candidate%nodes(p)%state==fixture%nodes(p)%state),'Non-MOL inactive state preserved')
    enddo
    kept=candidate;kept_flux=actual_flux;kept_regs=actual_regs
    observed_endpoint=0
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,2,.false.,[.true.,.true.,.true.,.true.], &
      [0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,candidate,wide,faces,bc,wc,ch,sc,limiter,good, &
      target_depth=0,integrated_flux=actual_flux,redistribution='state_redist',mol_mode=.false., &
      state_redist_max_order=2,state_redist_registers=actual_regs,auxiliary_source=reject_endpoint, &
      auxiliary_change=sc)
    call require(.not.good.and.observed_endpoint==1,'Non-MOL first endpoint failure observed')
    call require(all(candidate%nodes(1)%state==kept%nodes(1)%state).and. &
      all(candidate%nodes(1)%temperature==kept%nodes(1)%temperature),'Non-MOL joint gas/temperature rollback')
    call require(all(actual_flux(1)%x==kept_flux(1)%x).and.all(actual_flux(1)%y==kept_flux(1)%y).and. &
      all(actual_flux(1)%z==kept_flux(1)%z),'Non-MOL flux rollback')
    call require(all(actual_regs(1)%coarse_delta==kept_regs(1)%coarse_delta),'Non-MOL register rollback')
  end subroutine check_nonmol_sdc
  subroutine check_godunov_halos(fixture,physical,narrow)
    type(eb_reactive_tree_3d), intent(in) :: fixture
    type(eb_geometry_3d), intent(in) :: physical(:),narrow(:)
    type(state_redist_register_3d), allocatable :: expected_regs(:),actual_regs(:),kept_regs(:)
    type(eb_geometry_3d) :: padded
    type(eb_geometry_3d), allocatable :: wide(:),insufficient(:)
    type(amr_tree_state_3d) :: predictor,rate,control,halo,old,kept,kept_halo,kept_old
    type(eb_integrated_flux_3d), allocatable :: flux(:),reference_flux(:),baseline(:),reference_baseline(:),kept_flux(:)
    real(dp) :: dt,scale,limiter
    real(dp) :: prepared(nv,-5:10,-5:10,-5:10),prepared_t(-5:10,-5:10,-5:10),prepared_vf(-5:10,-5:10,-5:10)
    logical :: good
    integer :: p
    character(len=32) :: faces(6)
    dt=1.e-10_dp;faces='outflow'
    call build_axis_plane_eb_geometry_3d(96,96,96,-2.5_dp,3.5_dp,-2.5_dp,3.5_dp, &
      -2.5_dp,3.5_dp,'z',plane_cut,padded,good)
    call require(good,'Godunov grow-10 padded geometry')
    call build_eb_tree_extended_geometry_3d(topology,physical,padded,wide,good,halo_width=10)
    call require(good,'Godunov grow-10 patch geometry')
    allocate(predictor%fields(count))
    do p=1,count
      allocate(predictor%fields(p)%values(nv,4,4,4),source=0._dp)
    enddo
    call compute_sdc_eb_tree_hydro_3d(species,fixture,predictor,narrow,dt,control,reference_flux,good, &
      boundary_faces=faces,target_depth=0,baseline_flux=reference_baseline)
    call require(good,'Godunov original valid-only control')
    call compute_sdc_eb_tree_hydro_3d(species,fixture,predictor,wide,dt,rate,flux,good, &
      boundary_faces=faces,target_depth=0,baseline_flux=baseline,spatial_rate_halos=halo,old_state_halos=old)
    call require(good,'Godunov spatial source export')
    call require(all(shape(halo%fields(1)%values)==[nv,16,16,16]),'Godunov grow-6 exported layout')
    call require(all(halo%fields(1)%values(:,7:10,7:10,7:10)==rate%fields(1)%values), &
      'Godunov exact valid RHS crop')
    scale=max(1._dp,maxval(abs(control%fields(1)%values)))
    call require(maxval(abs(rate%fields(1)%values-control%fields(1)%values))/scale<1.e-13_dp, &
      'Godunov wider evaluation valid RHS control')
    call require(maxval(abs(flux(1)%x-reference_flux(1)%x))+ &
      maxval(abs(flux(1)%y-reference_flux(1)%y))+maxval(abs(flux(1)%z-reference_flux(1)%z))<1.e-13_dp, &
      'Godunov wider evaluation accepted face control')
    call require(all(baseline(1)%x==reference_baseline(1)%x).and. &
      all(baseline(1)%y==reference_baseline(1)%y).and.all(baseline(1)%z==reference_baseline(1)%z), &
      'Godunov baseline face crop control')
    call require(any(halo%fields(1)%values(:,4:6,7:10,7:10)/=0), &
      'Godunov nonzero locally evaluated spatial ghost RHS')
    call fill_reactive_eb_tree_halo_3d(species,fixture,1,prepared,prepared_t,good, &
      boundary_faces=faces,ghost_volume_fraction=prepared_vf,halo_width=6, &
      target_volume_fraction=narrow(1)%volume_fraction,restrict_snapshot=.false.)
    call require(good,'Godunov independent prepared old state')
    call require(all(old%fields(1)%values==prepared),'Godunov exact prepared old halo crop')
    do p=2,count
      call require(all(rate%fields(p)%values==0).and..not.allocated(halo%fields(p)%values).and. &
        .not.allocated(old%fields(p)%values),'Godunov inactive selected levels')
    enddo
    control=rate
    call compute_state_redist_level_stage_3d(topology,narrow,old,halo,nv,0,dt,control,expected_regs,good, &
      boundary_faces=faces,max_order=2,register_scale=1._dp)
    call require(good,'Godunov independent ML stage control')
    call require(maxval(abs(expected_regs(1)%coarse_delta))>0,'Godunov nonzero coarse stage ledger')
    call compute_godunov_state_redist_eb_tree_hydro_3d(species,fixture,predictor,wide,dt,rate,flux,good, &
      boundary_faces=faces,target_depth=0,registers=actual_regs,max_order=2,correct_forward_flux=.true.)
    call require(good,'Godunov ML source/register composition')
    call require(all(rate%fields(1)%values==control%fields(1)%values),'Godunov redistributed valid RHS control')
    call require(actual_regs(1)%scale==1.and. &
      all(actual_regs(1)%coarse_delta==expected_regs(1)%coarse_delta),'Godunov exact factor1 coarse ledger')
    do p=2,count
      if(allocated(actual_regs(p)%fine_delta)) call require(all(actual_regs(p)%fine_delta==expected_regs(p)%fine_delta), &
        'Godunov fine ledger independent control')
    enddo
    ! A grow9 CTU input cannot supply the grow4 allowances needed to
    ! correct the outer grow3 RHS used by StateRedistribute.
    allocate(insufficient(count))
    do p=1,count
      call extract_eb_geometry_3d(wide(p),[2,2,2],topology%nodes(p)%shape+19,insufficient(p),good)
      call require(good,'Godunov grow9 correction rejection fixture')
    enddo
    kept=rate;kept_flux=flux;kept_regs=actual_regs
    call compute_godunov_state_redist_eb_tree_hydro_3d(species,fixture,predictor,insufficient,dt,rate,flux,good, &
      boundary_faces=faces,target_depth=0,registers=actual_regs,max_order=2,correct_forward_flux=.true.)
    call require(.not.good,'Godunov insufficient correction geometry rejected')
    call require(all(rate%fields(1)%values==kept%fields(1)%values),'Godunov insufficient geometry RHS rollback')
    call require(all(flux(1)%x==kept_flux(1)%x).and.all(flux(1)%y==kept_flux(1)%y).and. &
      all(flux(1)%z==kept_flux(1)%z),'Godunov insufficient geometry flux rollback')
    call require(all(actual_regs(1)%coarse_delta==kept_regs(1)%coarse_delta), &
      'Godunov insufficient geometry register rollback')
    call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,narrow,dt, &
      [.true.,.true.,.true.,.true.],control,reference_flux,halo,old,limiter,good,boundary_faces=faces,target_depth=0)
    call require(good,'Non-MOL diffusion raw spatial control')
    call compute_state_redist_level_stage_3d(topology,narrow,old,halo,nv,0,dt,control,expected_regs,good, &
      boundary_faces=faces,max_order=2,register_scale=1._dp)
    call require(good,'Non-MOL diffusion independent factor1 ML control')
    call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,narrow,dt, &
      [.true.,.true.,.true.,.true.],rate,flux,halo,old,limiter,good,boundary_faces=faces,target_depth=0, &
      registers=actual_regs,max_order=2)
    call require(good,'Non-MOL diffusion ML register export')
    call require(all(rate%fields(1)%values==control%fields(1)%values),'Non-MOL diffusion redistributed RHS control')
    call require(actual_regs(1)%scale==1.and.maxval(abs(actual_regs(1)%coarse_delta))>0, &
      'Non-MOL diffusion nonzero factor1 ledger')
    call require(all(actual_regs(1)%coarse_delta==expected_regs(1)%coarse_delta), &
      'Non-MOL diffusion exact coarse ledger control')
    call check_nonmol_sdc(fixture,narrow,wide)
    kept=rate;kept_halo=halo;kept_old=old;kept_flux=flux
    predictor%fields(1)%values(1,1,1,4)=ieee_value(0._dp,ieee_quiet_nan)
    call compute_sdc_eb_tree_hydro_3d(species,fixture,predictor,wide,dt,rate,flux,good, &
      boundary_faces=faces,target_depth=0,spatial_rate_halos=halo,old_state_halos=old)
    call require(.not.good,'Godunov invalid source rejected')
    call require(all(rate%fields(1)%values==kept%fields(1)%values).and. &
      all(halo%fields(1)%values==kept_halo%fields(1)%values).and. &
      all(old%fields(1)%values==kept_old%fields(1)%values),'Godunov joint field rollback')
    call require(all(flux(1)%x==kept_flux(1)%x).and.all(flux(1)%y==kept_flux(1)%y).and. &
      all(flux(1)%z==kept_flux(1)%z),'Godunov joint flux rollback')
  end subroutine check_godunov_halos
  subroutine check_selected_sdc(fixture,wide)
    type(eb_reactive_tree_3d), intent(in) :: fixture
    type(eb_geometry_3d), intent(in) :: wide(:)
    type(eb_reactive_tree_3d) :: candidate,predicted,kept
    type(elementary_reaction) :: reactions(1)
    type(elementary_reaction), allocatable :: loaded(:)
    type(amr_tree_state_3d) :: r0,r1,r2,halo,old_halo,extra,mol_halo,mol_old
    type(eb_integrated_flux_3d), allocatable :: f0(:),f1(:),f2(:),fx(:),actual_flux(:),kept_flux(:)
    type(state_redist_register_3d), allocatable :: reg0(:),reg1(:),reg2(:),actual_regs(:),kept_regs(:)
    real(dp), parameter :: dt=1.e-10_dp
    real(dp) :: bc(nv),wc(nv),ch(nv),sc(nv),aux(nv),lim,tt,cs,expected(nv,4,4,4), &
      q(reactive_nprim(size(species))),exchange(nv),frequency,provisional_coarse,endpoint_factor,fine_endpoint_factor
    logical :: good
    integer :: depth,p,x,y,z,flags_id,requested_iterations
    logical :: flags(4)
    character(len=32) :: faces(6)
    ! Keep a valid stoichiometric reaction while disabling its rate so the
    ! reacting iteration lifecycle is isolated from chemical evolution.
    call load_h2o2_elementary_mechanism(loaded,good)
    call require(good,'load zero-rate iteration mechanism')
    reactions(1)=loaded(1)
    reactions(1)%forward_rate%pre_exponential=0
    reactions(1)%reversible=.false.
    faces='outflow'
    provisional_coarse=0
    requested_iterations=1;endpoint_factor=1;fine_endpoint_factor=1
    if(test_mode=='inert_iters') then
      ! The nonreacting state still accepts both trapezoidal endpoints.
      ! The state/physical flux still use both independently evaluated RHSs.
      ! Diffusion.cpp adds coarse drho on every call; fine dm is inserted only
      ! when the nonzero reflux_factor registers the corresponding face flux.
      requested_iterations=2;endpoint_factor=1;fine_endpoint_factor=0
    endif
    do flags_id=0,1
    flags=flags_id==1
    do depth=0,2
      ! Independent two-stage construction with unweighted RHS and .5 ledgers.
      candidate=fixture
      call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,wide,dt,flags,r0,f0,halo,old_halo,lim, &
        good,boundary_faces=faces,target_depth=depth)
      call require(good,'prepared SDC old transport control')
      call compute_mol_eb_tree_hydro_3d(species,fixture,wide,dt,2,extra,fx,good,boundary_faces=faces,target_depth=depth, &
        spatial_rate_halos=mol_halo,old_state_halos=mol_old)
      call require(good,'prepared SDC old hydro control')
      do p=1,count
        if(topology%nodes(p)%level/=depth) cycle
        r0%fields(p)%values=r0%fields(p)%values+extra%fields(p)%values
        halo%fields(p)%values=halo%fields(p)%values+mol_halo%fields(p)%values
        f0(p)%x=f0(p)%x+fx(p)%x;f0(p)%y=f0(p)%y+fx(p)%y;f0(p)%z=f0(p)%z+fx(p)%z
      enddo
      call compute_state_redist_level_stage_3d(topology,wide,old_halo,halo,nv,depth,dt,r0,reg0,good, &
        boundary_faces=faces,max_order=2,register_scale=.5_dp)
      call require(good,'prepared SDC old StateRedist control')
      call native_post_density_control(fixture,r0,depth)
      predicted=fixture
      do p=1,count
        if(topology%nodes(p)%level/=depth) cycle
        predicted%nodes(p)%state=fixture%nodes(p)%state+dt*r0%fields(p)%values
        do z=1,4
        do y=1,4
        do x=1,4
          if(fixture%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          call reactive_conserved_to_primitive(species,predicted%nodes(p)%state(:,x,y,z), &
            predicted%nodes(p)%temperature(x,y,z),q,tt,cs,good)
          call require(good,'prepared SDC control predictor EOS')
          predicted%nodes(p)%temperature(x,y,z)=tt
        enddo
        enddo
        enddo
      enddo
      call compute_sdc_eb_tree_transport_halos_3d(species,transport,predicted,wide,dt,flags,r1,f1,halo,old_halo,lim, &
        good,boundary_faces=faces,target_depth=depth)
      call require(good,'prepared SDC new transport control')
      call compute_mol_eb_tree_hydro_3d(species,predicted,wide,dt,2,extra,fx,good, &
        boundary_faces=faces,target_depth=depth,spatial_rate_halos=mol_halo,old_state_halos=mol_old)
      call require(good,'prepared SDC new hydro control')
      do p=1,count
        if(topology%nodes(p)%level/=depth) cycle
        r1%fields(p)%values=r1%fields(p)%values+extra%fields(p)%values
        halo%fields(p)%values=halo%fields(p)%values+mol_halo%fields(p)%values
        f1(p)%x=f1(p)%x+fx(p)%x;f1(p)%y=f1(p)%y+fx(p)%y;f1(p)%z=f1(p)%z+fx(p)%z
      enddo
      call compute_state_redist_level_stage_3d(topology,wide,old_halo,halo,nv,depth,dt,r1,reg1,good, &
        boundary_faces=faces,max_order=2,register_scale=.5_dp)
      call require(good,'prepared SDC new StateRedist control')
      call native_post_density_control(fixture,r1,depth)
      do p=1,count
        if(allocated(reg1(p)%coarse_delta)) provisional_coarse=max(provisional_coarse, &
          maxval(abs(reg1(p)%coarse_delta)))
      enddo
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,requested_iterations, &
        test_mode=='multi_sdc',flags,[0._dp,0._dp,0._dp], &
        1.e-6_dp,1.e-12_dp,candidate,wide,faces,bc,wc,ch,sc,lim,good,target_depth=depth, &
        integrated_flux=actual_flux,redistribution='state_redist',mol_mode=.true.,state_redist_max_order=2, &
        state_redist_registers=actual_regs)
      call require(good,'prepared selected-level SDC step')
      do p=1,count
        if(topology%nodes(p)%level/=depth) then
          call require(all(candidate%nodes(p)%state==fixture%nodes(p)%state),'prepared SDC inactive level preserved')
          cycle
        endif
        expected=fixture%nodes(p)%state+.5_dp*dt*(r0%fields(p)%values+r1%fields(p)%values)
        call require(maxval(abs(candidate%nodes(p)%state-expected)/max(1._dp,abs(expected)))<1.e-13_dp, &
          'prepared SDC actual two-stage gas control')
        call require(maxval(abs(actual_flux(p)%x-.5_dp*(f0(p)%x+endpoint_factor*f1(p)%x)))<1.e-13_dp.and. &
          maxval(abs(actual_flux(p)%y-.5_dp*(f0(p)%y+endpoint_factor*f1(p)%y)))<1.e-13_dp.and. &
          maxval(abs(actual_flux(p)%z-.5_dp*(f0(p)%z+endpoint_factor*f1(p)%z)))<1.e-13_dp,'prepared SDC actual stage flux')
        if(allocated(actual_regs(p)%coarse_delta)) then
          call require(all(actual_regs(p)%coarse_delta==reg0(p)%coarse_delta+reg1(p)%coarse_delta), &
            'prepared SDC coarse stage register sum')
        endif
        if(allocated(actual_regs(p)%fine_delta)) then
          call require(all(actual_regs(p)%fine_delta==reg0(p)%fine_delta+fine_endpoint_factor*reg1(p)%fine_delta), &
            'prepared SDC fine stage register sum')
        endif
      enddo
      if(test_mode=='multi_sdc') then
        ! Evaluate the second endpoint at the independently checked first
        ! corrector, then restart the actual two-iteration interval from old.
        call compute_sdc_eb_tree_transport_halos_3d(species,transport,candidate,wide,dt,flags,r2,f2,halo,old_halo,lim, &
          good,boundary_faces=faces,target_depth=depth)
        call require(good,'second endpoint transport control')
        call compute_mol_eb_tree_hydro_3d(species,candidate,wide,dt,2,extra,fx,good, &
          boundary_faces=faces,target_depth=depth,spatial_rate_halos=mol_halo,old_state_halos=mol_old)
        call require(good,'second endpoint hydro control')
        do p=1,count
          if(topology%nodes(p)%level/=depth) cycle
          r2%fields(p)%values=r2%fields(p)%values+extra%fields(p)%values
          halo%fields(p)%values=halo%fields(p)%values+mol_halo%fields(p)%values
          f2(p)%x=f2(p)%x+fx(p)%x;f2(p)%y=f2(p)%y+fx(p)%y;f2(p)%z=f2(p)%z+fx(p)%z
        enddo
        call compute_state_redist_level_stage_3d(topology,wide,old_halo,halo,nv,depth,dt,r2,reg2,good, &
          boundary_faces=faces,max_order=2,register_scale=.5_dp)
        call require(good,'second endpoint StateRedist control')
      call native_post_density_control(fixture,r2,depth)
        candidate=fixture
        observed_endpoint=0
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,2,.true.,flags,[0._dp,0._dp,0._dp], &
          1.e-6_dp,1.e-12_dp,candidate,wide,faces,bc,wc,ch,sc,lim,good,target_depth=depth, &
          integrated_flux=actual_flux,redistribution='state_redist',mol_mode=.true.,state_redist_max_order=2, &
          state_redist_registers=actual_regs,auxiliary_source=observe_endpoint,auxiliary_change=aux)
        call require(good,'selected two-iteration prepared step')
        ! Later MOL passes rebuild spatial terms; the auxiliary new-time
        ! callback belongs only to the first predictor, as in Advance.cpp.
        call require(observed_endpoint==1,'MOL auxiliary callback stops after the first endpoint')
        do p=1,count
          if(topology%nodes(p)%level/=depth) cycle
          expected=fixture%nodes(p)%state+.5_dp*dt*(r0%fields(p)%values+r2%fields(p)%values)
          call require(maxval(abs(candidate%nodes(p)%state-expected)/max(1._dp,abs(expected)))<1.e-13_dp, &
            'two-iteration gas uses latest endpoint RHS')
          call require(maxval(abs(actual_flux(p)%x-.5_dp*(f0(p)%x+f2(p)%x)))<1.e-13_dp.and. &
            maxval(abs(actual_flux(p)%y-.5_dp*(f0(p)%y+f2(p)%y)))<1.e-13_dp.and. &
            maxval(abs(actual_flux(p)%z-.5_dp*(f0(p)%z+f2(p)%z)))<1.e-13_dp,'two-iteration latest accepted flux')
          if(allocated(actual_regs(p)%coarse_delta)) then
            call require(all(actual_regs(p)%coarse_delta==reg0(p)%coarse_delta+ &
              (reg2(p)%coarse_delta+reg1(p)%coarse_delta)),'all coarse evaluations remain in register')
          endif
          if(allocated(actual_regs(p)%fine_delta)) then
            call require(all(actual_regs(p)%fine_delta==reg0(p)%fine_delta+reg2(p)%fine_delta), &
              'provisional fine delta is excluded')
          endif
        enddo
      endif
    enddo
    enddo
    if(test_mode=='multi_sdc'.or.test_mode=='inert_iters') &
      call require(provisional_coarse>0,'nonzero provisional coarse register control')
    kept=candidate;kept_regs=actual_regs;kept_flux=actual_flux
    bc=17;wc=19;ch=23;sc=29;aux=31;lim=37
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,1,.false.,flags,[0._dp,0._dp,0._dp], &
      1.e-6_dp,1.e-12_dp,candidate,wide,faces,bc,wc,ch,sc,lim,good,target_depth=2, &
      integrated_flux=actual_flux,redistribution='state_redist',mol_mode=.true.,state_redist_max_order=2, &
      state_redist_registers=actual_regs,auxiliary_source=reject_endpoint,auxiliary_change=aux)
    call require(.not.good,'prepared SDC late endpoint rejection')
    call require(all(bc==17).and.all(wc==19).and.all(ch==23).and.all(sc==29).and.all(aux==31).and.lim==37, &
      'prepared SDC diagnostics joint rollback')
    do p=1,count
      call require(all(candidate%nodes(p)%state==kept%nodes(p)%state),'prepared SDC gas joint rollback')
      call require(all(candidate%nodes(p)%temperature==kept%nodes(p)%temperature),'prepared SDC EOS joint rollback')
      if(allocated(actual_regs(p)%coarse_delta)) then
        call require(all(actual_regs(p)%coarse_delta==kept_regs(p)%coarse_delta),'prepared SDC coarse register rollback')
      endif
      if(allocated(actual_regs(p)%fine_delta)) then
        call require(all(actual_regs(p)%fine_delta==kept_regs(p)%fine_delta),'prepared SDC fine register rollback')
      endif
      if(allocated(actual_flux(p)%x)) then
        call require(all(actual_flux(p)%x==kept_flux(p)%x).and.all(actual_flux(p)%y==kept_flux(p)%y).and. &
          all(actual_flux(p)%z==kept_flux(p)%z),'prepared SDC flux rollback')
      endif
    enddo
    print *, 'PASS: selected-level prepared MOL SDC, two-stage gas/flux/register controls and joint rollback'
  end subroutine check_selected_sdc
  subroutine observe_endpoint(snapshot,interval,stage,depth,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: snapshot
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p
    observed_endpoint=max(observed_endpoint,stage)
    do p=1,size(snapshot%nodes)
      if(allocated(rate%fields(p)%values)) rate%fields(p)%values=0
    enddo
    good=.true.
  end subroutine observe_endpoint
  subroutine reject_endpoint(snapshot,interval,stage,depth,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: snapshot
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p
    observed_endpoint=stage
    good=stage==0
    if(.not.good) return
    do p=1,size(snapshot%nodes)
      rate%fields(p)%values=0
    enddo
  end subroutine reject_endpoint
  subroutine check_total_transport(fixture,wide)
    type(eb_reactive_tree_3d), intent(in) :: fixture
    type(eb_geometry_3d), intent(in) :: wide(:)
    type(amr_tree_state_3d) :: held,reference,actual,legacy,halo,old
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp) :: theta,delta,scale
    integer :: p,x,y,z
    logical :: good
    allocate(held%fields(size(fixture%nodes)))
    do p=1,size(fixture%nodes)
      allocate(held%fields(p)%values(1,4,4,4),source=0._dp)
      do z=1,4
      do y=1,4
      do x=1,4
        if(fixture%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
        held%fields(p)%values(1,x,y,z)=(fixture%nodes(p)%state(5,x,y,z)- &
          sum(fixture%nodes(p)%state(2:4,x,y,z)**2)/(2*fixture%nodes(p)%state(1,x,y,z)))*(1-.001_dp*x*x)
      enddo
      enddo
      enddo
    enddo
    call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,wide,1.e-9_dp, &
      [.false.,.true.,.false.,.false.],reference,flux,halo,old,theta,good,constant_coefficients=[0._dp,1._dp])
    call require(good,'total transport reference')
    call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,wide,1.e-9_dp, &
      [.false.,.true.,.false.,.false.],legacy,flux,halo,old,theta,good, &
      constant_coefficients=[0._dp,1._dp],internal_energy=held)
    call require(good,'held transport retained')
    call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,wide,1.e-9_dp, &
      [.false.,.true.,.false.,.false.],actual,flux,halo,old,theta,good, &
      constant_coefficients=[0._dp,1._dp],internal_energy=held,primitive_from_total_energy=.true.)
    call require(good,'total transport ignores inconsistent held thermal state')
    delta=0
    do p=1,size(fixture%nodes)
      scale=max(1._dp,maxval(abs(reference%fields(p)%values)))
      call require(maxval(abs(actual%fields(p)%values-reference%fields(p)%values))<1.e-9_dp*scale, &
        'total transport matches independent no-held reference')
      delta=max(delta,maxval(abs(legacy%fields(p)%values-reference%fields(p)%values))/scale)
    enddo
    call require(delta>1.e-4_dp,'held discrepancy regression is discriminating')
  end subroutine check_total_transport
  subroutine check_prepared_halos()
    type(eb_reactive_tree_3d) :: fixture,ledger
    type(eb_geometry_3d), allocatable :: physical(:),wide(:)
    type(eb_geometry_3d) :: padded
    type(amr_tree_state_3d) :: rates,old,halos,plain,kept_rate,kept_old,kept_halo,redistributed
    type(eb_integrated_flux_3d), allocatable :: flux(:),plain_flux(:),kept_flux(:)
    type(eb_integrated_flux_3d) :: direct_flux
    real(dp), allocatable :: direct_state(:,:,:,:),direct_temp(:,:,:),direct_vf(:,:,:),ordinary(:,:,:,:),direct_halo(:,:,:,:)
    real(dp) :: direct_rate(nv,4,4,4),direct_body(nv,4,4,4),direct_exchange(nv),direct_theta,direct_frequency
    type(state_redist_register_3d), allocatable :: regs(:)
    type(rotating_frame_3d), allocatable :: spin
    real(dp), allocatable :: phi
    real(dp) :: q(reactive_nprim(size(species))),pos(3),tt,cs,exchange(nv),frequency,limiter, &
      kept_exchange(nv),kept_frequency,kept_limiter,inventory(nv)
    integer :: p,x,y,z,c,rotation_id,depth
    logical :: good
    character(len=32) :: faces(6)
    allocate(physical(count))
    do p=1,count
      physical(p)=saved%nodes(p)%geometry
    enddo
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'z',plane_cut,padded,good)
    call require(good,'tree transport halo padded geometry')
    call build_eb_tree_extended_geometry_3d(topology,physical,padded,wide,good)
    call require(good,'tree transport halo geometry')
    do rotation_id=0,1
      fixture=saved
      if(rotation_id==1) then
        allocate(spin,phi)
        spin%axis=3;spin%angular_speed=100;spin%axis_point=[.5_dp,.5_dp,.5_dp]
      endif
      do p=1,count
        do z=1,4
        do y=1,4
        do x=1,4
          if(physical(p)%volume_fraction(x,y,z)<=0) cycle
          pos=[physical(p)%x_lower+(x-.5_dp)*physical(p)%dx, &
            physical(p)%y_lower+(y-.5_dp)*physical(p)%dy,physical(p)%z_lower+(z-.5_dp)*physical(p)%dz]
          q=0;q(1)=1+.01_dp*pos(1);q(2:4)=.1_dp*pos;q(5)=250000+1000*pos(2)
          if((test_mode=='multi_sdc'.or.test_mode=='inert_iters'.or.test_mode=='godunov_halo')) then
            q(1)=1+.01_dp*pos(1)**2+.02_dp*pos(3)**2
            q(5)=250000+1000*pos(2)**2+3000*pos(3)**2
          endif
          do c=1,size(species)
            if(species(c)%name=='N2') q(c+5)=.95_dp
            if(species(c)%name=='H2') q(c+5)=.05_dp
          enddo
          if(allocated(spin)) then
            call rotating_potential_3d(spin,pos,phi,good)
            call require(good,'tree transport halo potential')
          endif
          call reactive_primitive_to_conserved(species,q,fixture%nodes(p)%state(:,x,y,z),tt,cs,good)
          call require(good,'tree transport halo initial EOS')
          fixture%nodes(p)%temperature(x,y,z)=tt
        enddo
        enddo
        enddo
      enddo
      ! Match the ordinary gas after restriction before changing energy
      ! representation. Restricting rho*Phi first would change this control.
      call average_down_reactive_eb_tree_3d(species,fixture,good)
      call require(good,'tree transport halo initial restriction')
      if(allocated(spin)) then
        do p=1,count
        do z=1,4
        do y=1,4
        do x=1,4
          if(physical(p)%volume_fraction(x,y,z)<=0) cycle
          pos=[physical(p)%x_lower+(x-.5_dp)*physical(p)%dx, &
            physical(p)%y_lower+(y-.5_dp)*physical(p)%dy,physical(p)%z_lower+(z-.5_dp)*physical(p)%dz]
          call rotating_potential_3d(spin,pos,phi,good)
          call require(good,'tree transport matched potential')
          fixture%nodes(p)%state(5,x,y,z)=fixture%nodes(p)%state(5,x,y,z)-fixture%nodes(p)%state(1,x,y,z)*phi
        enddo
        enddo
        enddo
        enddo
      endif
      if(.not.allocated(spin).and.test_mode=='spatial_halo') call check_total_transport(fixture,wide)
      faces='outflow'
      if(test_mode=='godunov_halo') then
        call check_godunov_halos(fixture,physical,wide)
        return
      endif
      if(test_mode=='recursive_srd'.or.test_mode=='common_srd') then
        call check_recursive_registers()
        return
      endif
      if(test_mode=='prepared_sdc'.or. &
    (test_mode=='multi_sdc'.or.test_mode=='inert_iters')) then
        call check_selected_sdc(fixture,wide)
        return
      endif
      do depth=-1,2
        if(depth==-1) then
          call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,wide,1.e-8_dp, &
            [.true.,.true.,.true.,.true.],rates,flux,halos,old,limiter,good,boundary_faces=faces,frame=spin, &
            embedded_wall_temperature=700._dp,embedded_noslip=.true., &
            embedded_exchange=exchange,embedded_frequency=frequency)
        else
          call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,wide,1.e-8_dp, &
            [.true.,.true.,.true.,.true.],rates,flux,halos,old,limiter,good,boundary_faces=faces,frame=spin, &
            target_depth=depth,embedded_wall_temperature=700._dp,embedded_noslip=.true., &
            embedded_exchange=exchange,embedded_frequency=frequency)
        endif
        call require(good,'tree transport prepared source and body wall')
        do p=1,count
          if(depth>=0.and.topology%nodes(p)%level/=depth) then
            call require(.not.allocated(halos%fields(p)%values),'tree transport halo inactive source absent')
            cycle
          endif
          call require(all(lbound(old%fields(p)%values)==1),'tree transport halo old state indexing')
          call require(maxval(abs(old%fields(p)%values(:,7:10,7:10,7:10)-fixture%nodes(p)%state), &
            mask=spread(physical(p)%volume_fraction>0,1,nv))==0, &
            'tree transport halo published energy retained')
          call require(all(halos%fields(p)%values(:,7:10,7:10,7:10)==rates%fields(p)%values), &
            'tree transport halo valid source equals accepted rate')
        enddo
        if(depth>=0) then
          redistributed=rates
          call compute_state_redist_level_stage_3d(topology,wide,old,halos,nv,depth,1.e-8_dp, &
            redistributed,regs,good,boundary_faces=faces,max_order=2,register_scale=.5_dp)
          call require(good,'tree transport source feeds StateRedist with MOL register factor')
        else
          call require(frequency>0.and.maxval(abs(exchange(2:4)))>0,'tree transport no-slip impulse active')
          call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,wide,1.e-8_dp, &
            [.true.,.true.,.true.,.true.],plain,plain_flux,kept_halo,kept_old,kept_limiter,good, &
            boundary_faces=faces,frame=spin)
          call require(good,'tree transport no-body control')
          ledger=fixture
          do p=1,count
            ledger%nodes(p)%state=rates%fields(p)%values-plain%fields(p)%values
          enddo
          call integrate_reactive_eb_tree_3d(ledger,inventory,good)
          call require(good,'tree transport body composite inventory')
          call require(maxval(abs(1.e-8_dp*inventory-exchange))<1.e-9_dp,'tree transport valid body ledger')
          if(.not.allocated(direct_state)) allocate(direct_state(nv,-5:10,-5:10,-5:10), &
            direct_temp(-5:10,-5:10,-5:10),direct_vf(-5:10,-5:10,-5:10),ordinary(nv,16,16,16))
          call fill_reactive_eb_tree_halo_3d(species,fixture,1,direct_state,direct_temp,good,boundary_faces=faces, &
            ghost_volume_fraction=direct_vf,halo_width=6,target_volume_fraction=wide(1)%volume_fraction, &
            restrict_snapshot=.false.,frame=spin)
          call require(good,'tree transport independent prepared root halo')
          ordinary=direct_state
          if(allocated(spin)) then
            do z=1,16
            do y=1,16
            do x=1,16
              if(wide(1)%volume_fraction(x,y,z)<=0) cycle
              pos=[wide(1)%x_lower+(x-.5_dp)*wide(1)%dx,wide(1)%y_lower+(y-.5_dp)*wide(1)%dy, &
                wide(1)%z_lower+(z-.5_dp)*wide(1)%dz]
              call rotating_potential_3d(spin,pos,phi,good)
              call require(good,'tree transport independent ghost potential')
              ordinary(5,x,y,z)=ordinary(5,x,y,z)+ordinary(1,x,y,z)*phi
            enddo
            enddo
            enddo
          endif
          call compute_sdc_eb_transport_halo_3d(species,transport,ordinary,direct_temp,wide(1),1.e-8_dp, &
            [.true.,.true.,.true.,.true.],direct_rate,direct_flux,direct_halo,direct_theta,good, &
            faces=[.true.,.true.,.true.,.true.,.true.,.true.],physical_faces=faces)
          call require(good,'tree transport independent ordinary-energy control')
          call compute_eb_wall_transport_source_3d(species,transport,ordinary,direct_temp,wide(1),[7,7,7],[10,10,10], &
            .true.,.true.,700._dp,direct_body,direct_exchange,direct_frequency,good)
          call require(good,'tree transport independent body wall control')
          direct_rate=direct_rate+direct_body
          call require(maxval(abs(rates%fields(1)%values-direct_rate)/max(1._dp,abs(direct_rate)))<1.e-12_dp, &
            'tree transport private energy matches prepared ghost control')
        endif
      enddo
      kept_rate=rates;kept_halo=halos;kept_old=old;kept_flux=flux
      kept_exchange=exchange;kept_frequency=frequency;kept_limiter=limiter
      fixture%nodes(count)%state(1,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
      call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,wide,1.e-8_dp, &
        [.true.,.true.,.true.,.true.],rates,flux,halos,old,limiter,good,boundary_faces=faces,frame=spin, &
        embedded_wall_temperature=700._dp,embedded_noslip=.true.,embedded_exchange=exchange,embedded_frequency=frequency)
      call require(.not.good,'tree transport late patch EOS rollback')
      do p=1,count
        call require(all(rates%fields(p)%values==kept_rate%fields(p)%values),'tree transport all rates rollback')
        if(allocated(halos%fields(p)%values)) then
          call require(all(halos%fields(p)%values==kept_halo%fields(p)%values),'tree transport RHS rollback')
          call require(all(old%fields(p)%values==kept_old%fields(p)%values),'tree transport old state rollback')
        endif
        if(allocated(flux(p)%x)) then
          call require(all(flux(p)%x==kept_flux(p)%x).and.all(flux(p)%y==kept_flux(p)%y).and. &
            all(flux(p)%z==kept_flux(p)%z),'tree transport flux rollback')
        endif
      enddo
      call require(all(exchange==kept_exchange).and.frequency==kept_frequency.and.limiter==kept_limiter, &
        'tree transport wall ledger and limiter rollback')
    enddo
    print *, 'PASS: prepared tree transport halos, rotating energy, body walls and StateRedist stage'
  end subroutine check_prepared_halos

  subroutine check_sdc_transport()
    type(amr_tree_state_3d) :: rates,old_rates
    type(eb_integrated_flux_3d), allocatable :: flux(:),old_flux(:)
    type(eb_reactive_tree_3d) :: result,bad
    real(dp) :: limiter
    type(amr_tree_state_3d) :: level_rates,reference_rates
    type(eb_integrated_flux_3d), allocatable :: level_flux(:),reference_flux(:)
    type(eb_reactive_tree_3d) :: endpoint,interpolated
    integer :: p
    call compute_sdc_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
      saved,rates,flux,limiter,ok)
    call require(ok,'SDC provisional transport')
    call require(limiter==theta,'SDC transport limiter matches Euler')
    result=saved
    do p=1,count
      result%nodes(p)%state=result%nodes(p)%state+1.e-8_dp*rates%fields(p)%values
    end do
    call synchronize_reactive_eb_tree_flux_3d(species,result,flux,ok)
    call require(ok,'SDC transport reconstructed synchronization')
    do p=1,count
      call require(all(abs(result%nodes(p)%state-tree%nodes(p)%state)<= &
        2.e-12_dp*max(1._dp,abs(tree%nodes(p)%state))),'SDC transport reconstructs Euler state')
    end do
    endpoint=saved;interpolated=saved
    do p=1,count
      if(saved%topology%nodes(p)%level>=2) cycle
      endpoint%nodes(p)%state=1.01_dp*saved%nodes(p)%state
      interpolated%nodes(p)%state=.75_dp*saved%nodes(p)%state+.25_dp*endpoint%nodes(p)%state
      interpolated%nodes(p)%temperature=.75_dp*saved%nodes(p)%temperature+.25_dp*endpoint%nodes(p)%temperature
    end do
    call compute_sdc_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
      interpolated,reference_rates,reference_flux,limiter,ok,target_depth=2)
    call require(ok,'explicit ancestor transport reference')
    call compute_sdc_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
      saved,level_rates,level_flux,limiter,ok,target_depth=2, &
      ancestor_start=saved,ancestor_end=endpoint,ancestor_fraction=[.25_dp,.25_dp])
    call require(ok,'selected level transport with ancestor time interpolation')
    do p=1,count
      if(saved%topology%nodes(p)%level==2) then
        call require(all(level_rates%fields(p)%values==reference_rates%fields(p)%values).and. &
          all(level_flux(p)%x==reference_flux(p)%x).and.all(level_flux(p)%y==reference_flux(p)%y).and. &
          all(level_flux(p)%z==reference_flux(p)%z),'exact interpolated level transport')
      else
        call require(all(level_rates%fields(p)%values==0).and..not.allocated(level_flux(p)%x), &
          'inactive level has zero rates and absent flux')
      end if
    end do
    old_rates=rates;old_flux=flux;bad=saved
    bad%nodes(count)%state(1,4,4,4)=-1
    limiter=17
    call compute_sdc_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
      bad,rates,flux,limiter,ok)
    call require(.not.ok.and.limiter==17,'SDC transport invalid fluid density')
    do p=1,count
      call require(all(rates%fields(p)%values==old_rates%fields(p)%values).and. &
        all(flux(p)%x==old_flux(p)%x).and.all(flux(p)%y==old_flux(p)%y).and. &
        all(flux(p)%z==old_flux(p)%z),'SDC transport complete output rollback')
    end do
  end subroutine check_sdc_transport

  subroutine unchanged_tree()
    integer :: n
    do n=1,size(tree%nodes)
      call require(all(tree%nodes(n)%state==saved%nodes(n)%state).and. &
        all(tree%nodes(n)%temperature==saved%nodes(n)%temperature),'all state and temperature rollback')
    end do
  end subroutine unchanged_tree
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine require
end program test_eb_tree_transport_3d
