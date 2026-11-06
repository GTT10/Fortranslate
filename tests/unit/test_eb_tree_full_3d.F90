program test_eb_tree_full_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_full_3d, &
    advance_reactive_eb_tree_transport_3d,advance_reactive_eb_tree_hydro_3d,react_reactive_eb_tree_3d,react_sdc_eb_tree_3d, &
    compute_mol_eb_tree_hydro_3d,fill_reactive_eb_tree_halo_3d,build_eb_tree_extended_geometry_3d, &
    advance_reactive_eb_tree_sdc_3d,advance_reactive_eb_tree_sdc_subcycle_3d
  use eb_mol_hydro_3d_mod, only: compute_mol_eb_hydro_3d
  use eb_reactive_redistribution_3d_mod, only: compute_state_redist_level_stage_3d,state_redist_register_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d,react_sdc_eb_state_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved,reference
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:),boundary(:),primitive(:),wall(:),chemical(:),exchange(:),ref_b(:),ref_w(:),ref_c(:)
  real(dp) :: residual,theta,speed,ref_theta,stage_theta
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k,mode,integrator_index,half
  character(len=16) :: test_mode,order_arg
  integer :: redist_order=0
  logical :: ok,chemistry_on,transport_on
  integer :: callback_mode=0,callback_calls=0
  character(len=8), parameter :: integrators(2)=[character(len=8)::'explicit','implicit']
  call get_command_argument(2,order_arg)
  if(len_trim(order_arg)>0) read(order_arg,*) redist_order
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  call get_command_argument(1,test_mode)
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),primitive(reactive_nprim(size(species))))
  allocate(wall(nv),chemical(nv),exchange(nv),ref_b(nv),ref_w(nv),ref_c(nv))
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
  if(test_mode=='spatial_callback') then
    call check_spatial_callback()
    stop
  endif
  if(test_mode=='spatial_halo') then
    call check_spatial_halos()
    stop
  endif
  if(test_mode/='redist'.and.test_mode/='redist_recursive'.and.test_mode/='state_redist'.and. &
    test_mode/='mol_redist') call check_sdc_reaction()
  if(test_mode=='stored_reactor') then
    call check_sdc_reaction()
    deallocate(u,t)
    cycle
  endif
  call check_sdc_loop()
  if(test_mode=='redist'.or.test_mode=='redist_recursive'.or.test_mode=='state_redist'.or.test_mode=='mol_redist') then
    deallocate(u,t)
    cycle
  end if
  do integrator_index=1,2
  do mode=0,3
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
    ! Check split order against separately invoked public operators once.
    if (mode==3.and.integrator_index==1) then
      reference=saved;ref_c=0._dp;ref_b=0._dp;ref_theta=1._dp
      do half=1,2
        if (half==1) call reaction()
        call advance_reactive_eb_tree_transport_3d(species,transport,0.5e-8_dp,.true.,.true.,.true.,.true., &
          reference,exchange,stage_theta,ok)
        call require(ok,'reference transport')
        ref_b=ref_b+exchange;ref_theta=min(ref_theta,stage_theta)
        if (half==1) then
          call advance_reactive_eb_tree_hydro_3d(species,'hllc',1.e-8_dp,reference,exchange,ref_w,ok)
          call require(ok,'reference hydro')
          ref_b=ref_b+exchange
        else
          call reaction()
        end if
      end do
      do node=1,count
        call require(all(tree%nodes(node)%state==reference%nodes(node)%state).and. &
          all(tree%nodes(node)%temperature==reference%nodes(node)%temperature),'exact independent split state')
      end do
      call require(all(boundary==ref_b).and.all(wall==ref_w).and.all(chemical==ref_c).and. &
        theta==ref_theta,'exact independent split diagnostics')
    end if
  end do
  end do
  ! Hydro rejection follows successful private reaction and transport halves.
  tree=saved;boundary=17._dp;wall=19._dp;chemical=23._dp;theta=29._dp
  call advance_reactive_eb_tree_full_3d(species,reactions,transport,'invalid',1.e-8_dp,.true., &
    2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true.,tree,boundary,wall,chemical,theta,ok,'explicit')
  call require(.not.ok,'failure after earlier stages')
  call unchanged_tree()
  call require(all(boundary==17._dp).and.all(wall==19._dp).and.all(chemical==23._dp).and. &
    theta==29._dp,'complete diagnostic rollback')
  deallocate(u,t)
  end do
  end do
  print *, 'test_eb_tree_full_3d: PASS'
contains
  subroutine check_spatial_callback()
    type(eb_reactive_tree_3d) :: result,baseline
    type(eb_geometry_3d) :: padded
    type(eb_geometry_3d), allocatable :: physical(:),extended(:)
    real(dp) :: b(nv),w(nv),ch(nv),forcing(nv),limit,v0(nv),v1(nv)
    integer :: p,recursive,method,failure,width
    logical :: mol
    allocate(physical(count))
    do p=1,count
      physical(p)=saved%nodes(p)%geometry
    enddo
    call build_axis_plane_eb_geometry_3d(96,96,96,-2.5_dp,3.5_dp,-2.5_dp,3.5_dp,-2.5_dp,3.5_dp, &
      axes(axis),.31_dp,padded,ok)
    call require(ok,'serial spatial callback padded geometry')
    call integrate_reactive_eb_tree_3d(saved,v0,ok)
    call require(ok,'serial spatial callback initial inventory')
    do method=0,1
      mol=method==1;width=10
      if(mol) then
        width=6
        call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
          axes(axis),.31_dp,padded,ok)
        call require(ok,'serial MOL callback padded geometry')
      endif
      call build_eb_tree_extended_geometry_3d(topology,physical,padded,extended,ok,halo_width=width)
      call require(ok,'serial spatial callback extended geometry')
      do recursive=0,1
        result=saved;callback_mode=3;callback_calls=0
        call callback_advance(result,extended,recursive,mol,b,w,ch,forcing,limit,ok)
        call require(ok.and.callback_calls>=2,'serial spatial callback old/new execution')
        baseline=result
        result=saved;callback_mode=0;callback_calls=0
        call callback_advance(result,extended,recursive,mol,b,w,ch,forcing,limit,ok)
        call require(ok.and.callback_calls>=2,'serial spatial callback nonzero operator')
        call require(any(result%nodes(1)%state/=baseline%nodes(1)%state),'serial spatial callback changes gas')
        call integrate_reactive_eb_tree_3d(result,v1,ok)
        call require(ok,'serial spatial callback final inventory')
        call require(maxval(abs(v1-v0-b-w-ch-forcing)/max(1._dp,abs(v0),abs(v1)))<3.e-14_dp, &
          'serial spatial callback composite flux budget')
        do failure=1,2
          result=saved;callback_mode=failure;callback_calls=0
          b=17;w=19;ch=23;forcing=29;limit=37
          call callback_advance(result,extended,recursive,mol,b,w,ch,forcing,limit,ok)
          call require(.not.ok,'serial spatial callback malformed output rejected')
          do p=1,count
            call require(all(result%nodes(p)%state==saved%nodes(p)%state).and. &
              all(result%nodes(p)%temperature==saved%nodes(p)%temperature),'serial spatial callback state rollback')
          enddo
          call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(forcing==29).and.limit==37, &
            'serial spatial callback diagnostic rollback')
        enddo
      enddo
    enddo
  end subroutine check_spatial_callback
  subroutine callback_advance(result,extended,recursive,mol,b,w,ch,forcing,limit,good)
      type(eb_reactive_tree_3d), intent(inout) :: result
      type(eb_geometry_3d), intent(in) :: extended(:)
      integer, intent(in) :: recursive
      logical, intent(in) :: mol
      real(dp), intent(inout) :: b(:),w(:),ch(:),forcing(:),limit
      logical, intent(out) :: good
      if(recursive==1) then
        call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,2,.false., &
          [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp, &
          result,extended,[character(len=7)::'outflow','outflow','outflow','outflow','outflow','outflow'], &
          b,w,ch,forcing,limit,good,redistribution='state_redist',mol_mode=mol, &
          evaluate_spatial_diffusion=spatial_callback)
      else
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.false., &
          [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp, &
          result,extended,[character(len=7)::'outflow','outflow','outflow','outflow','outflow','outflow'], &
          b,w,ch,forcing,limit,good,redistribution='state_redist',mol_mode=mol, &
          evaluate_spatial_diffusion=spatial_callback)
      endif
  end subroutine callback_advance
  subroutine spatial_callback(state,geometry,old_halos,interval,depth,rate,spatial_rate,flux,limiter,good,frame)
    type(eb_reactive_tree_3d), intent(in) :: state
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: old_halos
    real(dp), intent(in) :: interval
    integer, intent(in) :: depth
    type(amr_tree_state_3d), intent(inout) :: rate,spatial_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: limiter
    logical, intent(out) :: good
    type(rotating_frame_3d), optional, intent(in) :: frame
    integer :: p,d(3),x,y,z
    real(dp) :: lo,hi,amplitude,pi
    real(dp), allocatable :: misplaced(:,:,:,:)
    callback_calls=callback_calls+1;good=.false.
    if(callback_mode==2.and.modulo(callback_calls,2)==0) return
    amplitude=1.e8_dp;pi=acos(-1._dp)
    if(callback_mode/=0) amplitude=0
    allocate(spatial_rate%fields(count),flux(count))
    do p=1,count
      rate%fields(p)%values=0
      if(depth>=0.and.state%topology%nodes(p)%level/=depth) cycle
      d=state%topology%nodes(p)%shape
      call require(allocated(old_halos%fields(p)%values),'serial spatial callback actual old donor')
      allocate(spatial_rate%fields(p)%values(nv,d(1)+12,d(2)+12,d(3)+12),source=0._dp)
      do z=4,d(3)+9
      do y=4,d(2)+9
      do x=4,d(1)+9
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        lo=geometry(p)%y_lower+real(y-1,dp)*geometry(p)%dy;hi=lo+geometry(p)%dy
        spatial_rate%fields(p)%values(5,x,y,z)=amplitude*(sin(2*pi*lo)-sin(2*pi*hi))/geometry(p)%dy
      enddo
      enddo
      enddo
      rate%fields(p)%values=spatial_rate%fields(p)%values(:,7:d(1)+6,7:d(2)+6,7:d(3)+6)
      allocate(flux(p)%x(nv,0:d(1),d(2),d(3)),flux(p)%y(nv,d(1),0:d(2),d(3)), &
        flux(p)%z(nv,d(1),d(2),0:d(3)))
      flux(p)%x=0;flux(p)%y=0;flux(p)%z=0
      do y=0,d(2)
        lo=state%nodes(p)%geometry%y_lower+real(y,dp)*geometry(p)%dy
        flux(p)%y(5,:,y,:)=interval*amplitude*sin(2*pi*lo)
      enddo
      if(callback_mode==1.and.modulo(callback_calls,2)==0) then
        allocate(misplaced(nv,0:d(1)+11,0:d(2)+11,0:d(3)+11),source=0._dp)
        call move_alloc(misplaced,spatial_rate%fields(p)%values)
      endif
    enddo
    limiter=1;good=.true.
  end subroutine spatial_callback
  subroutine check_spatial_halos()
    type(eb_geometry_3d) :: padded
    type(eb_geometry_3d), allocatable :: physical(:),wide(:)
    type(amr_tree_state_3d) :: rate,halos,old,kept_rate,kept_halos,kept_old,redistributed, &
      aux_rate,baseline_rate,kept_aux_rate,aux_halos,kept_aux_halos,held_state,held_end,aux_old,kept_aux_old, &
      aux_redistributed
    type(eb_integrated_flux_3d), allocatable :: fluxes(:),kept_fluxes(:),aux_flux(:),baseline_flux(:),kept_aux_flux(:)
    type(state_redist_register_3d), allocatable :: regs(:)
    type(state_redist_register_3d), allocatable :: aux_regs(:)
    type(eb_integrated_flux_3d) :: direct_flux
    real(dp), allocatable :: direct_halo(:,:,:,:),direct_state(:,:,:,:),direct_temperature(:,:,:),direct_vf(:,:,:), &
      direct_held(:,:,:)
    real(dp) :: direct_rate(nv,4,4,4),b(nv),w(nv),q(reactive_nprim(size(species))),tt,cs
    integer :: p,i,j,k,c,depth,solver
    real(dp) :: position(3),aux_inventory,aux_boundary,coarse_old
    logical :: good
    allocate(physical(count))
    ! Smooth physical field isolates halo/source assembly from the unrelated
    ! large pressure jumps used by the full chemistry/reflux test fixture.
    do p=1,count
      physical(p)=tree%nodes(p)%geometry
      do k=1,4
      do j=1,4
      do i=1,4
        if(physical(p)%volume_fraction(i,j,k)<=0) cycle
        position=[physical(p)%x_lower+(i-.5_dp)*physical(p)%dx, &
          physical(p)%y_lower+(j-.5_dp)*physical(p)%dy,physical(p)%z_lower+(k-.5_dp)*physical(p)%dz]
        q=0;q(1)=1+.01_dp*position(1);q(2:4)=.1_dp*position;q(5)=135000
        do c=1,size(species)
          if(species(c)%name=='N2') q(c+5)=1
        enddo
        call reactive_primitive_to_conserved(species,q,tree%nodes(p)%state(:,i,j,k),tree%nodes(p)%temperature(i,j,k),cs,good)
        call require(good,'hydro halo smooth state')
      enddo
      enddo
      enddo
    enddo
    call average_down_reactive_eb_tree_3d(species,tree,good)
    call require(good,'hydro halo smooth restriction')
    allocate(held_state%fields(count))
    do p=1,count
      allocate(held_state%fields(p)%values(1,4,4,4))
      held_state%fields(p)%values=0
      do k=1,4
      do j=1,4
      do i=1,4
        if(physical(p)%volume_fraction(i,j,k)<=0) cycle
        held_state%fields(p)%values(1,i,j,k)=tree%nodes(p)%state(5,i,j,k)- &
          .5_dp*sum(tree%nodes(p)%state(2:4,i,j,k)**2)/tree%nodes(p)%state(1,i,j,k)+real(p,dp)
      enddo
      enddo
      enddo
    enddo
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'z',.31_dp,padded,good)
    call require(good,'hydro halo padded geometry')
    call build_eb_tree_extended_geometry_3d(topology,physical,padded,wide,good)
    call require(good,'hydro halo patch geometry')
    do solver=0,1
      call compute_mol_eb_tree_hydro_3d(species,tree,wide,1.e-9_dp,2,rate,fluxes,good, &
        lax_friedrichs=solver==1,spatial_rate_halos=halos,old_state_halos=old)
      call require(good,'hydro halo whole tree spatial evaluation')
      baseline_rate=rate;baseline_flux=fluxes
      ! Observer-only outputs must not change gas evolution. Supplying a held
      ! state below also selects spatial coarse interpolation, so it is a
      ! different input contract and is not an observer-only comparison.
      call compute_mol_eb_tree_hydro_3d(species,tree,wide,1.e-9_dp,2,rate,fluxes,good, &
        lax_friedrichs=solver==1,spatial_rate_halos=halos,old_state_halos=old, &
        internal_energy_rate=aux_rate,internal_energy_flux=aux_flux, &
        internal_energy_spatial_rate_halos=aux_halos)
      call require(good,'hydro halo observer-only spatial evaluation')
      do p=1,count
        call require(all(rate%fields(p)%values==baseline_rate%fields(p)%values), &
          'retained energy output preserves tree gas rate')
        call require(all(fluxes(p)%x==baseline_flux(p)%x).and.all(fluxes(p)%y==baseline_flux(p)%y).and. &
          all(fluxes(p)%z==baseline_flux(p)%z),'retained energy output preserves tree gas flux')
      enddo
      call compute_mol_eb_tree_hydro_3d(species,tree,wide,1.e-9_dp,2,rate,fluxes,good, &
        lax_friedrichs=solver==1,spatial_rate_halos=halos,old_state_halos=old, &
        internal_energy_rate=aux_rate,internal_energy_flux=aux_flux, &
        internal_energy_spatial_rate_halos=aux_halos,internal_energy_state=held_state, &
        internal_energy_old_state_halos=aux_old)
      call require(good,'hydro halo retained energy spatial evaluation')
      do p=1,count
        call require(all(shape(aux_rate%fields(p)%values)==[1,4,4,4]),'hydro retained energy rate shape')
        call require(all(shape(aux_halos%fields(p)%values)==[1,16,16,16]), &
          'hydro retained energy halo shape')
        call require(all(shape(aux_old%fields(p)%values)==[1,16,16,16]), &
          'hydro retained energy old halo shape')
        call require(all(aux_halos%fields(p)%values(1,7:10,7:10,7:10)== &
          aux_rate%fields(p)%values(1,:,:,:)),'hydro retained energy valid halo agrees with corrected rate')
        aux_inventory=1.e-9_dp*physical(p)%dx*physical(p)%dy*physical(p)%dz* &
          sum(physical(p)%volume_fraction*aux_rate%fields(p)%values(1,:,:,:))
        aux_boundary=physical(p)%dy*physical(p)%dz*sum( &
          physical(p)%x_face_fraction(0,:,:)*aux_flux(p)%x(1,0,:,:)- &
          physical(p)%x_face_fraction(4,:,:)*aux_flux(p)%x(1,4,:,:))+ &
          physical(p)%dx*physical(p)%dz*sum( &
          physical(p)%y_face_fraction(:,0,:)*aux_flux(p)%y(1,:,0,:)- &
          physical(p)%y_face_fraction(:,4,:)*aux_flux(p)%y(1,:,4,:))+ &
          physical(p)%dx*physical(p)%dy*sum( &
          physical(p)%z_face_fraction(:,:,0)*aux_flux(p)%z(1,:,:,0)- &
          physical(p)%z_face_fraction(:,:,4)*aux_flux(p)%z(1,:,:,4))
        call require(abs(aux_inventory-aux_boundary)<1.e-14_dp,'retained energy tree spatial conservation')
        call require(all(shape(halos%fields(p)%values)==[nv,16,16,16]),'hydro halo whole tree source shape')
        call require(all(lbound(old%fields(p)%values)==1),'hydro halo state uses one-based prepared coordinates')
        call require(all(rate%fields(p)%values==halos%fields(p)%values(:,7:10,7:10,7:10)), &
          'hydro halo valid source equals published rate')
      enddo
      if(solver==0) then
        if(.not.allocated(direct_state)) allocate(direct_state(nv,-5:10,-5:10,-5:10), &
          direct_temperature(-5:10,-5:10,-5:10),direct_vf(-5:10,-5:10,-5:10),direct_held(-5:10,-5:10,-5:10))
        call fill_reactive_eb_tree_halo_3d(species,tree,2,direct_state,direct_temperature,good, &
          ghost_volume_fraction=direct_vf,halo_width=6,target_volume_fraction=wide(2)%volume_fraction, &
          restrict_snapshot=.false.,internal_energy=held_state,ghost_internal_energy=direct_held, &
          coarse_geometry=wide,preserve_temperature_seed=.true.)
        call require(good,'hydro halo direct fill')
        call require(all(direct_state==old%fields(2)%values),'hydro halo old state equals prepared snapshot')
        call require(all(direct_held(1:4,1:4,1:4)==held_state%fields(2)%values(1,:,:,:)), &
          'hydro halo retained energy valid cells use stored density')
        call require(all(direct_held==aux_old%fields(2)%values(1,:,:,:)), &
          'hydro halo retained energy tree old halo equals direct fill')
        call require(abs(direct_held(5,2,2)-held_state%fields(3)%values(1,1,2,2))<1.e-10_dp, &
          'hydro halo retained energy peer donor')
        coarse_old=direct_held(2,2,5)
        call require(coarse_old>0,'hydro halo retained energy coarse donor')
        held_end=held_state
        held_end%fields(1)%values=held_end%fields(1)%values+1._dp
        call fill_reactive_eb_tree_halo_3d(species,tree,2,direct_state,direct_temperature,good, &
          ghost_volume_fraction=direct_vf,halo_width=6,target_volume_fraction=wide(2)%volume_fraction, &
          restrict_snapshot=.false.,internal_energy=held_state,ghost_internal_energy=direct_held, &
          ancestor_start=tree,ancestor_end=tree,ancestor_fraction=[.5_dp], &
          ancestor_internal_energy_start=held_state,ancestor_internal_energy_end=held_end, &
          coarse_geometry=wide,preserve_temperature_seed=.true.)
        call require(good,'hydro halo retained energy ancestor interpolation')
        call require(abs(direct_held(2,2,5)-coarse_old-.5_dp)<1.e-10_dp, &
          'hydro halo retained energy time-interpolated coarse donor')
        call compute_mol_eb_hydro_3d(species,wide(2),direct_state,direct_temperature,1.e-9_dp,2, &
          direct_rate,direct_flux,b,w,good,spatial_rate_halo=direct_halo,primitive_from_total_energy=.true.)
        call require(good,'hydro halo direct local spatial calculation')
        call require(maxval(abs(halos%fields(2)%values(:,7,7,11)-direct_halo(:,7,7,11))/ &
          max(1._dp,abs(direct_halo(:,7,7,11))))<1.e-12_dp,'hydro coarse ghost RHS evaluated locally')
      endif
      do depth=0,2
        redistributed=rate
        call compute_state_redist_level_stage_3d(topology,wide,old,halos,nv,depth,1.e-9_dp,redistributed,regs,good,max_order=2)
        call require(good,'hydro halo feeds per-level StateRedist')
        aux_redistributed=aux_rate
        call compute_state_redist_level_stage_3d(topology,wide,aux_old,aux_halos,1,depth,1.e-9_dp, &
          aux_redistributed,aux_regs,good,max_order=2)
        call require(good,'retained energy halo feeds per-level StateRedist')
      enddo
    enddo
    kept_rate=rate;kept_halos=halos;kept_old=old;kept_fluxes=fluxes
    kept_aux_rate=aux_rate;kept_aux_flux=aux_flux;kept_aux_halos=aux_halos;kept_aux_old=aux_old
    tree%nodes(count)%state(1,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
    call compute_mol_eb_tree_hydro_3d(species,tree,wide,1.e-9_dp,2,rate,fluxes,good, &
      spatial_rate_halos=halos,old_state_halos=old,internal_energy_rate=aux_rate,internal_energy_flux=aux_flux, &
      internal_energy_spatial_rate_halos=aux_halos,internal_energy_state=held_state, &
      internal_energy_old_state_halos=aux_old)
    call require(.not.good,'hydro halo late patch failure')
    do p=1,count
      call require(all(rate%fields(p)%values==kept_rate%fields(p)%values),'hydro halo tree rates rollback')
      call require(all(halos%fields(p)%values==kept_halos%fields(p)%values),'hydro halo tree source rollback')
      call require(all(old%fields(p)%values==kept_old%fields(p)%values),'hydro halo tree old state rollback')
      call require(all(fluxes(p)%x==kept_fluxes(p)%x).and.all(fluxes(p)%y==kept_fluxes(p)%y).and. &
        all(fluxes(p)%z==kept_fluxes(p)%z),'hydro halo tree flux rollback')
      call require(all(aux_rate%fields(p)%values==kept_aux_rate%fields(p)%values), &
        'hydro halo tree retained energy rate rollback')
      call require(all(aux_halos%fields(p)%values==kept_aux_halos%fields(p)%values), &
        'hydro halo tree retained energy halo rollback')
      call require(all(aux_old%fields(p)%values==kept_aux_old%fields(p)%values), &
        'hydro halo tree retained energy old halo rollback')
      call require(all(aux_flux(p)%x==kept_aux_flux(p)%x).and.all(aux_flux(p)%y==kept_aux_flux(p)%y).and. &
        all(aux_flux(p)%z==kept_aux_flux(p)%z),'hydro halo tree retained energy flux rollback')
    enddo
    tree=saved
    print *, 'PASS: tree spatial halos, local coarse ghost source, StateRedist stage and rollback'
  end subroutine check_spatial_halos
  subroutine check_sdc_redistribution(initial,extended)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: extended(:)
    type(eb_reactive_tree_3d) :: result,baseline,previous
    type(amr_tree_state_3d) :: ir,saved_ir
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),v0(nv),v1(nv),limiter,error,change
    integer(int64) :: counts(3)
    integer :: p,step,c
    logical :: good,mol
    character(len=16) :: faces(6),method
    faces='outflow'
    mol=test_mode=='mol_redist'
    allocate(ir%fields(count))
    do p=1,count
      allocate(ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
    end do
    saved_ir=ir
    call integrate_reactive_eb_tree_3d(initial,v0,good)
    call require(good,'redistribution initial inventory')
    do step=1,2
      result=initial;ir=saved_ir
      method='none'
      if(step==2) method='flux_redist'
      if(step==2.and.test_mode=='state_redist') method='state_redist'
      if(test_mode=='redist_recursive') then
        call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,2,.true., &
          [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended,faces, &
          b,w,ch,ext,limiter,good,ir,level_steps=counts,redistribution=method,state_redist_max_order=redist_order)
        call require(all(counts==[1_int64,2_int64,4_int64]),'redistributed recursive counts')
      else
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
        [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended,faces, &
        b,w,ch,ext,limiter,good,ir,redistribution=method,state_redist_max_order=redist_order,mol_mode=mol)
      end if
      call require(good,'SDC redistributed coupled step')
      call integrate_reactive_eb_tree_3d(result,v1,good)
      call require(good,'redistributed final inventory')
      error=maxval(abs(v1-v0-b-w-ch-ext)/max(1._dp,abs(v0)))
      print *, 'SDC redistribution method, closure: ',method,error
      call require(error<2.e-11_dp,'redistributed coupled closure')
      call require(any(ch(6:)/=0).and.any(ext(2:4)/=0),'redistributed active reaction and forcing')
      do p=1,count
        do c=6,nv
          call require(all(pack(result%nodes(p)%state(c,:,:,:),result%nodes(p)%geometry%volume_fraction>0)>=0), &
            'redistributed published species nonnegative')
        enddo
      enddo
      if(step==1) baseline=result
    end do
    change=0
    do p=1,count
      change=max(change,maxval(abs(result%nodes(p)%state-baseline%nodes(p)%state)))
    end do
    call require(change>0,'redistribution changes nonuniform stage')
    previous=result;saved_ir=ir;b=17;w=19;ch=23;ext=29;limiter=31
    if(mol) then
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
        [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended,faces, &
        b,w,ch,ext,limiter,good,ir,redistribution=method,mol_mode=.true.,evaluate_reaction=refuse_reactor)
    else if(test_mode=='redist_recursive') then
      counts=37
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,2,.true., &
        [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended,faces, &
        b,w,ch,ext,limiter,good,ir,difmag=-1._dp,level_steps=counts,redistribution=method,state_redist_max_order=redist_order)
      call require(all(counts==37),'redistributed recursive counts rollback')
    else
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended,faces, &
      b,w,ch,ext,limiter,good,ir,difmag=-1._dp,redistribution=method,state_redist_max_order=redist_order)
    end if
    call require(.not.good,'redistributed rejected stage')
    do p=1,count
      call require(all(result%nodes(p)%state==previous%nodes(p)%state),'redistributed state rollback')
      call require(all(result%nodes(p)%temperature==previous%nodes(p)%temperature),'redistributed T rollback')
      call require(all(ir%fields(p)%values==saved_ir%fields(p)%values),'redistributed IR rollback')
    end do
    call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.limiter==31, &
      'redistributed diagnostic rollback')
  end subroutine check_sdc_redistribution

  subroutine refuse_reactor(species,reactions,initial_energy_density,external_species, &
      external_energy,duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
    real(dp), intent(inout) :: species_density(:),temperature
    integer, intent(out) :: completed_steps
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    species_density=-1;temperature=-1;completed_steps=0;ok=.false.
  end subroutine refuse_reactor
  subroutine check_sdc_loop()
    type(eb_reactive_tree_3d) :: initial,result
    type(eb_geometry_3d) :: padded
    type(eb_geometry_3d), allocatable :: physical(:),extended(:)
    type(amr_tree_state_3d) :: ir,previous_ir,heat,saved_heat
    real(dp) :: b(nv),w(nv),chem(nv),forcing(nv),v0(nv),v1(nv),limit,q(size(primitive)),temp,sound
    type(eb_reactive_tree_3d) :: endpoint,interpolated,reference_level
    type(amr_tree_state_3d) :: op,np,ref_op,ref_np,ref_ir
    type(eb_integrated_flux_3d), allocatable :: flux(:),ref_flux(:)
    real(dp) :: ref_diag(nv,4),ref_limit
    integer(int64) :: counts(3)
    integer :: p,x,y,z
    initial=saved
    do p=1,count
      do z=1,4
      do y=1,4
      do x=1,4
        if(initial%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
        call reactive_conserved_to_primitive(species,initial%nodes(p)%state(:,x,y,z), &
          initial%nodes(p)%temperature(x,y,z),q,temp,sound,ok)
        call require(ok,'SDC loop initial EOS')
        q(5)=q(5)*1000._dp/temp
        call reactive_primitive_to_conserved(species,q,initial%nodes(p)%state(:,x,y,z), &
          initial%nodes(p)%temperature(x,y,z),sound,ok)
        call require(ok,'SDC loop warm fixture')
      end do
      end do
      end do
    end do
    call average_down_reactive_eb_tree_3d(species,initial,ok)
    call require(ok,'SDC loop initial restriction')
    allocate(physical(count),ir%fields(count),heat%fields(count))
    do p=1,count
      allocate(heat%fields(p)%values(1,4,4,4),source=0._dp)
      physical(p)=initial%nodes(p)%geometry
      allocate(ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
    end do
    previous_ir=ir
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'z',.31_dp,padded,ok)
    call require(ok,'SDC loop padded geometry')
    call build_eb_tree_extended_geometry_3d(topology,physical,padded,extended,ok)
    call require(ok,'SDC loop extended patches')
    if(test_mode=='redist'.or.test_mode=='redist_recursive'.or.test_mode=='state_redist'.or.test_mode=='mol_redist') then
      call check_sdc_redistribution(initial,extended)
      return
    end if
    result=initial
    call integrate_reactive_eb_tree_3d(initial,v0,ok)
    call require(ok,'SDC loop initial inventory')
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended, &
      [character(len=7)::'outflow','outflow','outflow','outflow','outflow','outflow'], &
      b,w,chem,forcing,limit,ok,ir,heat_release=heat)
    call require(ok,'coupled EB tree SDC step')
    call require_restricted_reactions(result,ir)
    call require_restricted_reactions(result,heat)
    call integrate_reactive_eb_tree_3d(result,v1,ok)
    call require(ok,'SDC loop final inventory')
    print *, 'SDC tree inventory residual:',maxval(abs(v1-v0-b-w-chem-forcing)/max(1._dp,abs(v0)))
    call require(all(abs(v1-v0-b-w-chem-forcing)<=2.e-11_dp*max(1._dp,abs(v0))),'SDC composite ledger closure')
    call require(any(chem(6:)/=0).and.any(forcing(2:4)/=0),'SDC active chemistry and forcing')
    call require(limit>=0.and.limit<=1,'SDC loop limiter')
    result=initial;ir=previous_ir;counts=-1
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended, &
      [character(len=7)::'outflow','outflow','outflow','outflow','outflow','outflow'], &
      b,w,chem,forcing,limit,ok,ir,level_steps=counts,heat_release=heat)
    call require(ok,'recursive EB tree SDC')
    call require_restricted_reactions(result,ir)
    call require_restricted_reactions(result,heat)
    call require(all(counts==[1_int64,2_int64,4_int64]),'recursive SDC level counts')
    call integrate_reactive_eb_tree_3d(result,v1,ok)
    call require(ok,'recursive SDC inventory')
    print *, 'Recursive SDC inventory residual:',maxval(abs(v1-v0-b-w-chem-forcing)/max(1._dp,abs(v0)))
    call require(all(abs(v1-v0-b-w-chem-forcing)<=2.e-11_dp*max(1._dp,abs(v0))), &
      'recursive SDC composite conservation')
    saved_heat=heat
    result=initial;ir=previous_ir;counts=17;b=19;w=23;chem=29;forcing=31;limit=37
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended, &
      [character(len=7)::'outflow','outflow','outflow','outflow','outflow','outflow'], &
      b,w,chem,forcing,limit,ok,ir,difmag=-1._dp,level_steps=counts,heat_release=heat)
    call require(.not.ok,'recursive SDC failed stage')
    do p=1,count
      call require(all(heat%fields(p)%values==saved_heat%fields(p)%values),'recursive heat rollback')
      call require(all(result%nodes(p)%state==initial%nodes(p)%state).and. &
        all(result%nodes(p)%temperature==initial%nodes(p)%temperature).and. &
        all(ir%fields(p)%values==previous_ir%fields(p)%values),'recursive SDC state rollback')
    end do
    call require(all(counts==17).and.all(b==19).and.all(w==23).and.all(chem==29).and. &
      all(forcing==31).and.limit==37,'recursive SDC counts and ledgers rollback')
    saved_heat=heat
    endpoint=initial;interpolated=initial
    do p=1,count
      if(initial%topology%nodes(p)%level>=2) cycle
      endpoint%nodes(p)%state=1.01_dp*initial%nodes(p)%state
      interpolated%nodes(p)%state=.75_dp*initial%nodes(p)%state+.25_dp*endpoint%nodes(p)%state
      interpolated%nodes(p)%temperature=.75_dp*initial%nodes(p)%temperature+.25_dp*endpoint%nodes(p)%temperature
    end do
    reference_level=interpolated;ref_ir=previous_ir
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,reference_level,extended, &
      [character(len=7)::'outflow','outflow','outflow','outflow','outflow','outflow'], &
      ref_diag(:,1),ref_diag(:,2),ref_diag(:,3),ref_diag(:,4),ref_limit,ok,ref_ir, &
      target_depth=2,integrated_flux=ref_flux,old_predictor=ref_op,new_predictor=ref_np,heat_release=heat)
    call require(ok,'SDC explicit ancestor level reference')
    result=initial;ir=previous_ir
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended, &
      [character(len=7)::'outflow','outflow','outflow','outflow','outflow','outflow'], &
      b,w,chem,forcing,limit,ok,ir,target_depth=2,ancestor_start=initial,ancestor_end=endpoint, &
      fractions_start=[.25_dp,.25_dp],fractions_end=[.25_dp,.25_dp],integrated_flux=flux, &
      old_predictor=op,new_predictor=np,heat_release=heat)
    call require(ok,'SDC coupled level with ancestor interpolation')
    call require(all(b==ref_diag(:,1)).and.all(w==ref_diag(:,2)).and.all(chem==ref_diag(:,3)).and. &
      all(forcing==ref_diag(:,4)).and.limit==ref_limit,'SDC exact level diagnostics')
    do p=1,count
      if(initial%topology%nodes(p)%level/=2) &
        call require(all(heat%fields(p)%values==saved_heat%fields(p)%values),'inactive heat retained')
      if(initial%topology%nodes(p)%level==2) then
        call require(all(result%nodes(p)%state==reference_level%nodes(p)%state).and. &
          all(result%nodes(p)%temperature==reference_level%nodes(p)%temperature).and. &
          all(ir%fields(p)%values==ref_ir%fields(p)%values),'SDC exact level gas and reaction')
        call require(all(flux(p)%x==ref_flux(p)%x).and.all(flux(p)%y==ref_flux(p)%y).and. &
          all(flux(p)%z==ref_flux(p)%z),'SDC exact level flux')
        call require(all(op%fields(p)%values==ref_op%fields(p)%values).and. &
          all(np%fields(p)%values==ref_np%fields(p)%values),'SDC exact predictors')
      else
        call require(all(result%nodes(p)%state==initial%nodes(p)%state).and. &
          all(result%nodes(p)%temperature==initial%nodes(p)%temperature).and. &
          all(ir%fields(p)%values==previous_ir%fields(p)%values),'SDC inactive level remains unchanged')
        call require(.not.allocated(flux(p)%x).and.all(op%fields(p)%values==0).and. &
          all(np%fields(p)%values==0),'SDC inactive flux and predictors')
      end if
    end do
    result=initial;ir=previous_ir;b=17;w=19;chem=23;forcing=29;limit=31
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended, &
      [character(len=7)::'outflow','outflow','outflow','outflow','outflow','outflow'], &
      b,w,chem,forcing,limit,ok,ir,difmag=-1._dp,heat_release=heat)
    call require(.not.ok,'SDC hydro failure after old transport')
    do p=1,count
      call require(all(result%nodes(p)%state==initial%nodes(p)%state).and. &
        all(result%nodes(p)%temperature==initial%nodes(p)%temperature).and. &
        all(ir%fields(p)%values==previous_ir%fields(p)%values),'SDC complete tree rollback')
    end do
    call require(all(b==17).and.all(w==19).and.all(chem==23).and.all(forcing==29).and.limit==31, &
      'SDC diagnostic rollback')
  end subroutine check_sdc_loop

  subroutine check_sdc_reaction()
    type(amr_tree_state_3d) :: rates,states,temps,ir,rs,rt,ri,heat,rh,held
    type(eb_reactive_tree_3d) :: bad
    real(dp), allocatable :: expected_ir(:,:,:,:),eu(:,:,:,:),et(:,:,:),expected_heat(:,:,:)
    real(dp) :: correction
    integer :: p,c,depth
    allocate(rates%fields(count),ir%fields(count),heat%fields(count))
    do p=1,count
      allocate(heat%fields(p)%values(1,4,4,4),source=31._dp)
      allocate(rates%fields(p)%values(nv,4,4,4),ir%fields(p)%values(nv-4,4,4,4))
      rates%fields(p)%values=0;ir%fields(p)%values=7
      rates%fields(p)%values(2,:,:,:)=2.e4_dp
      rates%fields(p)%values(5,:,:,:)=2.e4_dp*saved%nodes(p)%state(2,:,:,:)/saved%nodes(p)%state(1,:,:,:)
      ! Rates in solid cells are deliberately invalid and must not participate.
      do c=1,nv
        where(saved%nodes(p)%geometry%volume_fraction<=0) &
          rates%fields(p)%values(c,:,:,:)=ieee_value(0._dp,ieee_quiet_nan)
      end do
    end do
    call react_sdc_eb_tree_3d(species,reactions,saved,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp,states,temps,ir,ok,heat_release=heat)
    call require(ok,'SDC tree reaction')
    allocate(held%fields(count))
    do p=1,count
      allocate(held%fields(p)%values(1,4,4,4),source=23._dp)
    enddo
    call react_sdc_eb_tree_3d(species,reactions,saved,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp, &
      states,temps,ir,ok,heat_release=heat,internal_energy=held)
    call require(ok,'SDC tree retained reactor endpoint')
    do p=1,count
      call require(all(pack(held%fields(p)%values(1,:,:,:),saved%nodes(p)%geometry%volume_fraction<=0)==23), &
        'tree covered retained scalar unchanged')
    enddo
    do p=1,count
      eu=saved%nodes(p)%state;et=saved%nodes(p)%temperature
      allocate(expected_heat(4,4,4),source=31._dp)
      allocate(expected_ir(nv-4,4,4,4),source=7._dp)
      call react_sdc_eb_state_3d(species,reactions,saved%nodes(p)%state,saved%nodes(p)%temperature, &
        saved%nodes(p)%geometry,rates%fields(p)%values,1.e-9_dp,2.e-7_dp,1.e-12_dp,eu,et,expected_ir,ok,heat_release=expected_heat)
      call require(ok,'SDC patch reference')
      call require(all(states%fields(p)%values==eu).and.all(temps%fields(p)%values(1,:,:,:)==et).and. &
        all(ir%fields(p)%values==expected_ir),'SDC tree exact patch reference including refined parents')
      call require(all(heat%fields(p)%values(1,:,:,:)==expected_heat),'tree heat matches accepted patch stage')
      call require(all(pack(expected_heat,saved%nodes(p)%geometry%volume_fraction<=0)==0),'solid heat zero')
      deallocate(expected_ir,expected_heat)
    end do
    rs=states;rt=temps;ri=ir;correction=7
    call react_sdc_eb_tree_3d(species,reactions,saved,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp, &
      states,temps,ir,ok,internal_energy=held,energy_correction=correction)
    call require(ok.and.correction==7,'SDC tree retained EOS policy and correction ledger')
    do p=1,count
      call require(all(states%fields(p)%values==rs%fields(p)%values), &
        'retained EOS recovery does not change accepted gas without energy repair')
      call require(all(pack(held%fields(p)%values(1,:,:,:),saved%nodes(p)%geometry%volume_fraction<=0)==23), &
        'retained EOS policy preserves covered scalar')
    enddo
    states=rs;temps=rt;ir=ri
    do p=1,count
      held%fields(p)%values=29
    enddo
    call react_sdc_eb_tree_3d(species,reactions,saved,rates,-1.e-9_dp,2.e-7_dp,1.e-12_dp, &
      states,temps,ir,ok,internal_energy=held)
    call require(.not.ok,'SDC tree retained invalid interval rejects')
    do p=1,count
      call require(all(held%fields(p)%values==29),'SDC tree retained rollback')
    enddo
    rs=states;rt=temps;ri=ir;rh=heat
    do depth=0,2
      do p=1,count
        states%fields(p)%values=-17;temps%fields(p)%values=-19;ir%fields(p)%values=7
        heat%fields(p)%values=31
      end do
      call react_sdc_eb_tree_3d(species,reactions,saved,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp, &
        states,temps,ir,ok,target_depth=depth,heat_release=heat)
      call require(ok,'SDC selected level reaction')
      do p=1,count
        if(saved%topology%nodes(p)%level==depth) then
          call require(all(heat%fields(p)%values==rh%fields(p)%values),'selected heat exact')
          call require(all(states%fields(p)%values==rs%fields(p)%values).and. &
            all(temps%fields(p)%values==rt%fields(p)%values).and. &
            all(ir%fields(p)%values==ri%fields(p)%values),'SDC selected exact state and sources')
        else
          call require(all(heat%fields(p)%values==31),'inactive stage heat retained')
          call require(all(states%fields(p)%values==-17).and.all(temps%fields(p)%values==-19).and. &
            all(ir%fields(p)%values==7),'SDC inactive level retained')
        end if
      end do
    end do
    states=rs;temps=rt;ir=ri;heat=rh;bad=saved
    bad%nodes(count)%state(1,4,4,4)=-1
    call require(bad%nodes(count)%geometry%volume_fraction(4,4,4)>0,'late failing cell is fluid')
    call react_sdc_eb_tree_3d(species,reactions,bad,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp,states,temps,ir,ok,heat_release=heat)
    do p=1,count
      call require(all(heat%fields(p)%values==rh%fields(p)%values),'late reaction heat rollback')
    enddo
    call require(.not.ok,'SDC late node failure')
    do p=1,count
      call require(all(states%fields(p)%values==rs%fields(p)%values).and. &
        all(temps%fields(p)%values==rt%fields(p)%values).and. &
        all(ir%fields(p)%values==ri%fields(p)%values),'SDC complete tree rollback')
    end do
  end subroutine check_sdc_reaction

  subroutine reaction()
    call react_reactive_eb_tree_3d(species,reactions,0.5e-8_dp,2.e-7_dp,1.e-12_dp, &
      reference,exchange,ok,'explicit')
    call require(ok,'reference reaction')
    ref_c=ref_c+exchange
  end subroutine reaction
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
end program test_eb_tree_full_3d
