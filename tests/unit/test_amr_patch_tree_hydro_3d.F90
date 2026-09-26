program test_amr_patch_tree_hydro_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: advance_reactive_boundary_hydro_3d,compute_reactive_cfl_timestep_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: advance_amr_patch_tree_hydro_3d,compute_amr_patch_tree_cfl_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_sdc_amr_patch_tree_hydro_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,reflux_amr_patch_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,temperature,saved,saved_t
  type(amr_patch_3d) :: patches(5)
  real(dp) :: root(12,4,4,4),rt(1,4,4,4),serial(12,4,4,4),st(4,4,4),q(12),sound
  real(dp) :: spacing(3),dt,root_dt,expected_dt,local_dt,h(3),budget(12),serial_budget(12),largest,initial(12),final(12)
  real(dp) :: outlet(6),mach(6),tin(6),pin(6)
  integer :: i,j,k,mode,geometry,node,n(3)
  logical :: ok,ppm
  character(len=32) :: faces(6),reconstruction,arg
  call get_command_argument(1,arg)
  ppm=trim(arg)=='ppm'
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  spacing=[.001_dp,.0015_dp,.002_dp]
  do k=1,4
    do j=1,4
      do i=1,4
        q=0;q(1)=.3_dp*(1+.02_dp*(i+j+k));q(5)=101325*(1+.02_dp*k)
        q(2:4)=[10.0_dp+i,20.0_dp+j,-4.0_dp+.5_dp*k]
        q(6)=.025_dp+.003_dp*i;q(9)=.22_dp;q(12)=1-q(6)-q(9)
        call reactive_primitive_to_conserved(species,q,root(:,i,j,k),rt(1,i,j,k),sound,ok)
        call require(ok,'initial state')
      end do
    end do
  end do
  call initialize_amr_patch_3d(4,4,4,1,2,2,3,1,4,2,patches(1),ok)
  call initialize_amr_patch_3d(4,4,4,3,4,2,3,1,4,2,patches(2),ok)
  call initialize_amr_patch_3d(4,4,8,3,4,1,2,1,2,2,patches(3),ok)
  call initialize_amr_patch_3d(4,4,8,1,2,1,2,1,2,2,patches(4),ok)
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,1,2,2,patches(5),ok)
  do geometry=1,2
    if (geometry==1) then
      call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches(:0),tree,ok)
    else
      call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3,4],patches,tree,ok)
    end if
    call require(ok,'topology')
    do mode=1,8
      if (ppm.and.(mode==2.or.mode==4.or.mode==6)) cycle
      outlet=0;mach=-1;tin=700;pin=0
      if (mode>4) then
        outlet(2)=102325;outlet(6)=102325
      end if
      faces='periodic'
      if (mode>2) then
        faces='outflow';faces(3:4)='slip_wall'
      end if
      reconstruction='pcm'
      if (mod(mode,2)==0) reconstruction='characteristic_plm'
      if(mode>6) then
        outlet=0;faces(1)='inflow';mach(1)=.2_dp
        if(mode==8) pin(1)=120000
        reconstruction='pcm'
        if(mode==8) reconstruction='characteristic_plm'
      end if
      if (ppm) reconstruction='characteristic_ppm'
      call prolong_amr_patch_tree_3d(tree,root,state,ok)
      call require(ok,'state initialization')
      call prolong_amr_patch_tree_3d(tree,rt,temperature,ok)
      call require(ok,'temperature initialization')
      call compute_amr_patch_tree_cfl_3d(tree,species,state,temperature,spacing,.15_dp,dt,ok)
      call require(ok.and.dt>0,'CFL')
      if(ppm.and.mode==1) call check_sdc_hydro_stage()
      if (geometry==1) root_dt=dt
      expected_dt=huge(1.0_dp)
      do node=1,size(tree%nodes)
        n=tree%nodes(node)%shape;h=spacing/2**tree%nodes(node)%level
        call compute_reactive_cfl_timestep_3d(species,state%fields(node)%values, &
          temperature%fields(node)%values(1,:,:,:),n(1),n(2),n(3),h(1),h(2),h(3),.15_dp,local_dt,ok)
        call require(ok,'local CFL oracle')
        expected_dt=min(expected_dt,local_dt)
      end do
      call require(dt==expected_dt,'hierarchy CFL minimum')
      if (geometry==2) call require(dt<root_dt/4,'finest grid CFL controls step')
      call composite_integrals_amr_patch_tree_3d(tree,state,spacing,initial,ok)
      call require(ok,'initial inventory')
      call advance_amr_patch_tree_hydro_3d(tree,species,state,temperature,spacing,dt,faces, &
        reconstruction,'mc','hllc',budget,largest,ok,outlet_pressure=outlet, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
      call require(ok,'hydro step')
      if (geometry==1) then
        serial=root;st=rt(1,:,:,:)
        call advance_reactive_boundary_hydro_3d(species,serial,st,spacing,dt,faces, &
          reconstruction,'mc','hllc',serial_budget,ok,outlet_pressure=outlet, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
        call require(ok,'serial oracle')
        call require(maxval(abs(state%fields(1)%values-serial)/max(1.0_dp,abs(serial)))<2e-13_dp,'serial state')
        call require(maxval(abs(temperature%fields(1)%values(1,:,:,:)-st))<1e-9_dp,'serial temperature')
        call require(maxval(abs(budget-serial_budget))<1e-14_dp,'serial boundary ledger')
      end if
      call composite_integrals_amr_patch_tree_3d(tree,state,spacing,final,ok)
      call require(ok,'final inventory')
      call require(maxval(abs(final-initial-budget)/max(1e-12_dp,abs(initial)))<2e-12_dp,'composite conservation')
      do node=1,size(tree%nodes)
        call require(all(temperature%fields(node)%values>0),'positive temperature')
        call require(maxval(abs(sum(state%fields(node)%values(6:,:,:,:),dim=1)- &
          state%fields(node)%values(1,:,:,:)))<1e-12_dp,'species closure')
      end do
      saved=state;saved_t=temperature;budget=-17;largest=-19
      call advance_amr_patch_tree_hydro_3d(tree,species,state,temperature,spacing,1.0_dp,faces, &
        reconstruction,'mc','hllc',budget,largest,ok,outlet_pressure=outlet, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
      call require(.not.ok,'nonphysical candidate rejection')
      do node=1,size(tree%nodes)
        call require(all(state%fields(node)%values==saved%fields(node)%values),'state rollback')
        call require(all(temperature%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
      end do
      call require(all(budget==-17).and.largest==-19,'diagnostics rollback')
    end do
  end do
contains
  subroutine check_sdc_hydro_stage()
    type(amr_tree_state_3d) :: source,rates,candidate,expected,et,forced,level_rates,packed,future,source_future,retained
    type(amr_tree_flux_3d), allocatable :: flux(:),forced_flux(:),level_flux(:)
    real(dp) :: boundary(12),correction(12),maximum,change
    integer :: peer,depth,nv,dims(3)
    real(dp), allocatable :: fractions(:)
    source=state
    do peer=1,size(tree%nodes)
      source%fields(peer)%values=0
    end do
    call compute_sdc_amr_patch_tree_hydro_3d(tree,species,state,temperature,source,spacing,dt,faces, &
      rates,flux,ok)
    call require(ok,'SDC hierarchy CTU stage')
    candidate=state
    do peer=1,size(tree%nodes)
      candidate%fields(peer)%values=state%fields(peer)%values+dt*rates%fields(peer)%values
    end do
    call reflux_amr_patch_tree_3d(tree,candidate,flux,spacing,faces,correction,maximum,ok)
    call require(ok,'SDC stage integrated flux reflux')
    expected=state;et=temperature
    call advance_amr_patch_tree_hydro_3d(tree,species,expected,et,spacing,dt,faces, &
      'characteristic_ppm','mc','pelec',boundary,maximum,ok)
    call require(ok,'SDC hierarchy CTU reference')
    do peer=1,size(tree%nodes)
      call require(maxval(abs(candidate%fields(peer)%values-expected%fields(peer)%values)/ &
        max(1._dp,abs(expected%fields(peer)%values)))<1e-12_dp,'SDC hydro stage plus reflux matches existing tree step')
      source%fields(peer)%values(2,:,:,:)=20000
    end do
    call compute_sdc_amr_patch_tree_hydro_3d(tree,species,state,temperature,source,spacing,dt,faces, &
      forced,forced_flux,ok)
    call require(ok,'SDC hierarchy signed source interpolation')
    change=0
    do peer=1,size(tree%nodes)
      change=max(change,maxval(abs(forced%fields(peer)%values-rates%fields(peer)%values)))
    end do
    call require(change>1e-10_dp,'SDC predictor source changes CTU rate')
    do depth=0,tree%level_count()-1
      call compute_sdc_amr_patch_tree_hydro_3d(tree,species,state,temperature,source,spacing,dt,faces, &
        level_rates,level_flux,ok,target_depth=depth)
      call require(ok,'level-selective SDC hydro')
      do peer=1,size(tree%nodes)
        if(tree%nodes(peer)%level==depth) then
          call require(all(level_rates%fields(peer)%values==forced%fields(peer)%values),'selected level rates exact')
          call require(all(level_flux(peer)%x==forced_flux(peer)%x),'selected level flux exact')
        else
          call require(all(level_rates%fields(peer)%values==0),'inactive rate zero')
          call require(.not.allocated(level_flux(peer)%x),'inactive flux absent')
        end if
      end do
    end do
    depth=tree%level_count()-1
    if(depth==0) return
    allocate(packed%fields(size(tree%nodes)),fractions(depth))
    nv=size(state%fields(1)%values,1)
    do peer=1,size(tree%nodes)
      dims=tree%nodes(peer)%shape
      allocate(packed%fields(peer)%values(nv+1,dims(1),dims(2),dims(3)))
      packed%fields(peer)%values(:nv,:,:,:)=state%fields(peer)%values
      packed%fields(peer)%values(nv+1,:,:,:)=temperature%fields(peer)%values(1,:,:,:)
    end do
    future=packed;source_future=source
    do peer=1,size(tree%nodes)
      if(tree%nodes(peer)%level>=depth) cycle
      future%fields(peer)%values(5,:,:,:)=1.0001_dp*packed%fields(peer)%values(5,:,:,:)
      source_future%fields(peer)%values=2*source%fields(peer)%values
    end do
    fractions=0
    call compute_sdc_amr_patch_tree_hydro_3d(tree,species,state,temperature,source,spacing,dt,faces, &
      level_rates,level_flux,ok,target_depth=depth,coarse_start=packed,coarse_end=future,coarse_fractions=fractions, &
      source_start=source,source_end=source_future,source_fractions=fractions)
    call require(ok,'coarse old endpoint sampling')
    do peer=1,size(tree%nodes)
      if(tree%nodes(peer)%level/=depth) cycle
      call require(all(level_rates%fields(peer)%values==forced%fields(peer)%values),'old endpoint exact')
    end do
    fractions=.5_dp
    retained=level_rates
    call compute_sdc_amr_patch_tree_hydro_3d(tree,species,state,temperature,source,spacing,dt,faces, &
      level_rates,level_flux,ok,target_depth=depth,coarse_start=packed,coarse_end=future,coarse_fractions=fractions, &
      source_start=source,source_end=source_future,source_fractions=fractions)
    call require(ok,'coarse intermediate state and signed source sampling')
    change=0
    do peer=1,size(tree%nodes)
      change=max(change,maxval(abs(level_rates%fields(peer)%values-retained%fields(peer)%values)))
    end do
    call require(change>1e-10_dp,'coarse temporal donors change fine rate')
    retained=level_rates;fractions=1.1_dp
    call compute_sdc_amr_patch_tree_hydro_3d(tree,species,state,temperature,source,spacing,dt,faces, &
      level_rates,level_flux,ok,target_depth=depth,coarse_start=packed,coarse_end=future,coarse_fractions=fractions)
    call require(.not.ok,'invalid coarse time rejected')
    do peer=1,size(tree%nodes)
      call require(all(level_rates%fields(peer)%values==retained%fields(peer)%values),'invalid coarse time rollback')
    end do
  end subroutine check_sdc_hydro_stage
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
