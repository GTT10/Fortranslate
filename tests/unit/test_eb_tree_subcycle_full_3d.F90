program test_eb_tree_subcycle_full_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_full_3d, &
    advance_reactive_eb_tree_transport_3d,advance_reactive_eb_tree_hydro_3d, &
    advance_reactive_eb_tree_subcycle_transport_3d,advance_reactive_eb_tree_subcycle_hydro_3d,react_reactive_eb_tree_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
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
  real(dp), allocatable :: ref_h(:)
  real(dp) :: residual,theta,speed,ref_theta,stage_theta
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k,mode,integrator_index,half
  logical :: ok,chemistry_on,transport_on,recursive_steps
  real(dp) :: mach(6),tin(6),pin(6)
  character(len=16) :: test_mode
  character(len=8), parameter :: integrators(2)=[character(len=8)::'explicit','implicit']
  call get_command_argument(1,test_mode)
  recursive_steps=trim(test_mode)/='mach-common'
  mach=-1;tin=400;pin=0
  if (index(test_mode,'mach')==1) then
    mach(1)=.2_dp;mach(4)=.1_dp;pin(4)=120000
  end if
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),primitive(reactive_nprim(size(species))))
  allocate(ref_h(nv))
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
  do integrator_index=1,2
  do mode=0,3
    chemistry_on=btest(mode,0);transport_on=btest(mode,1)
    tree=saved
    call advance_reactive_eb_tree_full_3d(species,reactions,transport,'hllc',1.e-8_dp,chemistry_on, &
      2.e-7_dp,1.e-12_dp,transport_on,.true.,.true.,.true.,.true.,tree,boundary,wall,chemical,theta,ok, &
      integrators(integrator_index),subcycling=recursive_steps, &
      inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
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
        if (recursive_steps) then
          call advance_reactive_eb_tree_subcycle_transport_3d(species,transport,0.5e-8_dp,.true.,.true.,.true.,.true., &
            reference,exchange,stage_theta,ok,mach,tin,pin)
        else
          call advance_reactive_eb_tree_transport_3d(species,transport,0.5e-8_dp,.true.,.true.,.true.,.true., &
            reference,exchange,stage_theta,ok,mach,tin,pin)
        end if
        call require(ok,'reference transport')
        if (any(mach>=0)) call require(abs(exchange(5))>1e-12_dp,'active inlet transport energy flux')
        ref_b=ref_b+exchange;ref_theta=min(ref_theta,stage_theta)
        if (half==1) then
          if (recursive_steps) then
            call advance_reactive_eb_tree_subcycle_hydro_3d(species,'hllc',1.e-8_dp,reference,exchange,ref_w,ok,mach,tin,pin)
          else
            call advance_reactive_eb_tree_hydro_3d(species,'hllc',1.e-8_dp,reference,exchange,ref_w,ok,mach,tin,pin)
          end if
          call require(ok,'reference hydro')
          ref_h=exchange
        else
          call reaction()
        end if
      end do
      ! Match the public ledger definition H + (T/2 + T/2); addition is not associative.
      ref_b=ref_h+ref_b
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
    2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true.,tree,boundary,wall,chemical,theta,ok, &
      'explicit',subcycling=recursive_steps, &
      inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
  call require(.not.ok,'failure after earlier stages')
  call unchanged_tree()
  call require(all(boundary==17._dp).and.all(wall==19._dp).and.all(chemical==23._dp).and. &
    theta==29._dp,'complete diagnostic rollback')
  if (any(mach>=0)) then
    tin(1)=ieee_value(0._dp,ieee_quiet_nan)
    call advance_reactive_eb_tree_full_3d(species,reactions,transport,'hllc',1.e-8_dp,.true., &
      2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true.,tree,boundary,wall,chemical,theta,ok, &
      'explicit',subcycling=recursive_steps,inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
    call require(.not.ok,'invalid inlet after private reaction')
    call unchanged_tree()
    call require(all(boundary==17).and.all(wall==19).and.all(chemical==23).and.theta==29, &
      'inlet failure diagnostic rollback')
  end if
  deallocate(u,t)
  end do
  end do
  print *, 'test_eb_tree_subcycle_full_3d: PASS'
contains
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
end program test_eb_tree_subcycle_full_3d
