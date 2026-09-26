program test_eb_tree_full_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
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
    build_eb_tree_extended_geometry_3d,advance_reactive_eb_tree_sdc_3d,advance_reactive_eb_tree_sdc_subcycle_3d
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
  logical :: ok,chemistry_on,transport_on
  character(len=8), parameter :: integrators(2)=[character(len=8)::'explicit','implicit']
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
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
  call check_sdc_reaction()
  call check_sdc_loop()
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
  subroutine check_sdc_loop()
    type(eb_reactive_tree_3d) :: initial,result
    type(eb_geometry_3d) :: padded
    type(eb_geometry_3d), allocatable :: physical(:),extended(:)
    type(amr_tree_state_3d) :: ir,previous_ir
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
    allocate(physical(count),ir%fields(count))
    do p=1,count
      physical(p)=initial%nodes(p)%geometry
      allocate(ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
    end do
    previous_ir=ir
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'z',.31_dp,padded,ok)
    call require(ok,'SDC loop padded geometry')
    call build_eb_tree_extended_geometry_3d(topology,physical,padded,extended,ok)
    call require(ok,'SDC loop extended patches')
    result=initial
    call integrate_reactive_eb_tree_3d(initial,v0,ok)
    call require(ok,'SDC loop initial inventory')
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended, &
      [character(len=7)::'outflow','outflow','outflow','outflow','outflow','outflow'], &
      b,w,chem,forcing,limit,ok,ir)
    call require(ok,'coupled EB tree SDC step')
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
      b,w,chem,forcing,limit,ok,ir,level_steps=counts)
    call require(ok,'recursive EB tree SDC')
    call require(all(counts==[1_int64,2_int64,4_int64]),'recursive SDC level counts')
    call integrate_reactive_eb_tree_3d(result,v1,ok)
    call require(ok,'recursive SDC inventory')
    print *, 'Recursive SDC inventory residual:',maxval(abs(v1-v0-b-w-chem-forcing)/max(1._dp,abs(v0)))
    call require(all(abs(v1-v0-b-w-chem-forcing)<=2.e-11_dp*max(1._dp,abs(v0))), &
      'recursive SDC composite conservation')
    result=initial;ir=previous_ir;counts=17;b=19;w=23;chem=29;forcing=31;limit=37
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended, &
      [character(len=7)::'outflow','outflow','outflow','outflow','outflow','outflow'], &
      b,w,chem,forcing,limit,ok,ir,difmag=-1._dp,level_steps=counts)
    call require(.not.ok,'recursive SDC failed stage')
    do p=1,count
      call require(all(result%nodes(p)%state==initial%nodes(p)%state).and. &
        all(result%nodes(p)%temperature==initial%nodes(p)%temperature).and. &
        all(ir%fields(p)%values==previous_ir%fields(p)%values),'recursive SDC state rollback')
    end do
    call require(all(counts==17).and.all(b==19).and.all(w==23).and.all(chem==29).and. &
      all(forcing==31).and.limit==37,'recursive SDC counts and ledgers rollback')
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
      target_depth=2,integrated_flux=ref_flux,old_predictor=ref_op,new_predictor=ref_np)
    call require(ok,'SDC explicit ancestor level reference')
    result=initial;ir=previous_ir
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],[2.e4_dp,-3.e4_dp,4.e4_dp],2.e-7_dp,1.e-12_dp,result,extended, &
      [character(len=7)::'outflow','outflow','outflow','outflow','outflow','outflow'], &
      b,w,chem,forcing,limit,ok,ir,target_depth=2,ancestor_start=initial,ancestor_end=endpoint, &
      fractions_start=[.25_dp,.25_dp],fractions_end=[.25_dp,.25_dp],integrated_flux=flux, &
      old_predictor=op,new_predictor=np)
    call require(ok,'SDC coupled level with ancestor interpolation')
    call require(all(b==ref_diag(:,1)).and.all(w==ref_diag(:,2)).and.all(chem==ref_diag(:,3)).and. &
      all(forcing==ref_diag(:,4)).and.limit==ref_limit,'SDC exact level diagnostics')
    do p=1,count
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
      b,w,chem,forcing,limit,ok,ir,difmag=-1._dp)
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
    type(amr_tree_state_3d) :: rates,states,temps,ir,rs,rt,ri
    type(eb_reactive_tree_3d) :: bad
    real(dp), allocatable :: expected_ir(:,:,:,:),eu(:,:,:,:),et(:,:,:)
    integer :: p,c,depth
    allocate(rates%fields(count),ir%fields(count))
    do p=1,count
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
    call react_sdc_eb_tree_3d(species,reactions,saved,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp,states,temps,ir,ok)
    call require(ok,'SDC tree reaction')
    do p=1,count
      eu=saved%nodes(p)%state;et=saved%nodes(p)%temperature
      allocate(expected_ir(nv-4,4,4,4),source=7._dp)
      call react_sdc_eb_state_3d(species,reactions,saved%nodes(p)%state,saved%nodes(p)%temperature, &
        saved%nodes(p)%geometry,rates%fields(p)%values,1.e-9_dp,2.e-7_dp,1.e-12_dp,eu,et,expected_ir,ok)
      call require(ok,'SDC patch reference')
      call require(all(states%fields(p)%values==eu).and.all(temps%fields(p)%values(1,:,:,:)==et).and. &
        all(ir%fields(p)%values==expected_ir),'SDC tree exact patch reference including refined parents')
      deallocate(expected_ir)
    end do
    rs=states;rt=temps;ri=ir
    do depth=0,2
      do p=1,count
        states%fields(p)%values=-17;temps%fields(p)%values=-19;ir%fields(p)%values=7
      end do
      call react_sdc_eb_tree_3d(species,reactions,saved,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp, &
        states,temps,ir,ok,target_depth=depth)
      call require(ok,'SDC selected level reaction')
      do p=1,count
        if(saved%topology%nodes(p)%level==depth) then
          call require(all(states%fields(p)%values==rs%fields(p)%values).and. &
            all(temps%fields(p)%values==rt%fields(p)%values).and. &
            all(ir%fields(p)%values==ri%fields(p)%values),'SDC selected exact state and sources')
        else
          call require(all(states%fields(p)%values==-17).and.all(temps%fields(p)%values==-19).and. &
            all(ir%fields(p)%values==7),'SDC inactive level retained')
        end if
      end do
    end do
    states=rs;temps=rt;ir=ri;bad=saved
    bad%nodes(count)%state(1,4,4,4)=-1
    call require(bad%nodes(count)%geometry%volume_fraction(4,4,4)>0,'late failing cell is fluid')
    call react_sdc_eb_tree_3d(species,reactions,bad,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp,states,temps,ir,ok)
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
end program test_eb_tree_full_3d
