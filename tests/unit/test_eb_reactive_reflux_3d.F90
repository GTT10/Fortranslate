program test_eb_reactive_reflux_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_nan
  use precision_mod, only: dp
  use constants_mod, only: density_floor
  use nasa7_thermo_mod, only: nasa7_species,configure_nasa7_eos_policy
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar
  use mixture_thermo_mod, only: mixture_mass_properties
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d, &
    extract_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,restrict_eb_integrated_flux_3d, &
    reflux_reactive_eb_state_3d,reconcile_reactive_eb_peer_flux_3d,compute_reactive_eb_reflux_increment_3d, &
    reconcile_reactive_eb_interface_budget_3d,supply_reactive_eb_interface_donors_3d, &
    supply_reactive_eb_endpoint_donor_3d,supply_reactive_eb_gas_eos_donor_3d,supply_reactive_eb_held_eos_donor_3d, &
    reactive_eb_exchange_interval_3d,assemble_eb_rereflux_3d,add_eb_rereflux_to_fine_3d,apply_reactive_eb_rereflux_pair_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    assemble_reactive_eb_tree_reflux_3d,synchronize_reactive_eb_tree_flux_3d,integrate_reactive_eb_tree_3d
  use eb_reactive_patch_tree_3d_mod, only: build_eb_tree_extended_geometry_3d, &
    synchronize_reactive_eb_tree_register_flux_3d,synchronize_reactive_eb_tree_scalar_flux_3d, &
    synchronize_reactive_eb_tree_scalar_register_flux_3d
  use eb_reactive_redistribution_3d_mod, only: state_redist_register_3d,initialize_state_redist_tree_register_3d, &
    assemble_state_redist_interface_rereflux_3d,assemble_state_redist_patch_rereflux_3d, &
    state_redist_box_3d,initialize_state_redist_register_3d,srd_coarse_boundary
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: fine,coarse
  type(eb_integrated_flux_3d) :: ff,cf,restricted,saved_flux
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:),saved(:,:,:,:),saved_t(:,:,:),reference(:),delta(:)
  real(dp), allocatable :: increment(:,:,:,:),saved_increment(:,:,:,:)
  logical, allocatable :: touched(:,:,:),saved_touched(:,:,:)
  logical :: unrefined(3,3,3)
  real(dp) :: fine_change,coarse_change,volume,scale
  logical :: refined(3,3,3),ok
  integer :: ratio,n,nvar,i,j,k,c,trace_species
  call check_rereflux()
  call check_physical_domain_rereflux()
  call check_redistribution_increment()
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo load')
  call check_trace_energy_donors()
  call check_reactor_energy_donors()
  call check_extrapolated_energy_donors()
  call check_redistribution_eos()
  call check_rereflux_pair()
  call check_tree_rereflux()
  nvar=reactive_nvar(size(species))
  allocate(state(nvar,3,3,3),temperature(3,3,3),saved(nvar,3,3,3),saved_t(3,3,3),reference(nvar),delta(nvar))
  allocate(cf%x(nvar,0:3,3,3),cf%y(nvar,3,0:3,3),cf%z(nvar,3,3,0:3))
  config%nx=3;config%ny=3;config%nz=3
  refined=.false.;refined(2,2,2)=.true.
  do ratio=2,3
    n=3*ratio
    call build_axis_plane_eb_geometry_3d(n,n,n,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.37_dp,fine,ok)
    call require(ok,'fine cut geometry')
    call coarsen_eb_geometry_3d(fine,ratio,coarse,ok)
    call require(ok,'coarse cut geometry')
    call initialize_reactive_eb_density_sheet_3d(species,config,coarse,state,temperature,ok)
    call require(ok,'initial coarse state')
    saved=state;saved_t=temperature;reference=state(:,2,2,2)
    allocate(ff%x(nvar,0:n,n,n),ff%y(nvar,n,0:n,n),ff%z(nvar,n,n,0:n))
    do k=1,n
      do j=1,n
        do i=0,n
          ff%x(:,i,j,k)=1.e-6_dp*(i+2*j+3*k)*reference
        end do
      end do
    end do
    do k=1,n
      do j=0,n
        do i=1,n
          ff%y(:,i,j,k)=2.e-6_dp*(i+2*j+3*k)*reference
        end do
      end do
    end do
    do k=0,n
      do j=1,n
        do i=1,n
          ff%z(:,i,j,k)=3.e-6_dp*(i+2*j+3*k)*reference
        end do
      end do
    end do
    call restrict_eb_integrated_flux_3d(fine,ratio,ff,restricted,ok)
    call require(ok,'aperture weighted flux restriction')
    call require(all(restricted%z(:,:, :,0)==0._dp),'closed face flux zero')
    cf%x=0._dp;cf%y=0._dp;cf%z=0._dp
    call compute_reactive_eb_reflux_increment_3d(coarse,refined,cf,restricted,increment,ok,touched=touched)
    call require(ok,'signed reflux increment assembly')
    call reflux_reactive_eb_state_3d(species,coarse,refined,cf,restricted,state,temperature,ok)
    call require(ok,'cut cell reflux')
    call require(all(state==saved+increment),'applied reflux equals the assembled signed increment exactly')
    call require(maxval(abs(state-saved))>1.e-10_dp,'reflux is nontrivial')
    call require(all(state(:,2,2,2)==saved(:,2,2,2)).and.temperature(2,2,2)==saved_t(2,2,2), &
      'refined coarse storage unchanged')
    call require(all(state(:,1,1,3)==saved(:,1,1,3)).and.temperature(1,1,3)==saved_t(1,1,3), &
      'noninterface storage unchanged')
    call require(all(state(:,:,:,1)==saved(:,:,:,1)),'covered storage unchanged')
    ! Independent fine-volume divergence over the refined cube. Its integrated
    ! inventory change must cancel the unrefined coarse-cell corrections.
    delta=0._dp
    do k=ratio+1,2*ratio
      do j=ratio+1,2*ratio
        do i=ratio+1,2*ratio
          delta=delta-fine%dy*fine%dz*(fine%x_face_fraction(i,j,k)*ff%x(:,i,j,k)- &
            fine%x_face_fraction(i-1,j,k)*ff%x(:,i-1,j,k)) &
            -fine%dx*fine%dz*(fine%y_face_fraction(i,j,k)*ff%y(:,i,j,k)- &
            fine%y_face_fraction(i,j-1,k)*ff%y(:,i,j-1,k)) &
            -fine%dx*fine%dy*(fine%z_face_fraction(i,j,k)*ff%z(:,i,j,k)- &
            fine%z_face_fraction(i,j,k-1)*ff%z(:,i,j,k-1))
        end do
      end do
    end do
    volume=coarse%dx*coarse%dy*coarse%dz
    do c=1,nvar
      coarse_change=sum((state(c,:,:,:)-saved(c,:,:,:))*coarse%volume_fraction)*volume
      fine_change=delta(c)
      scale=max(1._dp,sum(abs(saved(c,:,:,:))*coarse%volume_fraction)*volume)
      call require(abs(coarse_change+fine_change)<=2.e-12_dp*scale,'composite flux conservation')
    end do
    state=saved;temperature=saved_t
    call reflux_reactive_eb_state_3d(species,coarse,refined,restricted,restricted,state,temperature,ok)
    call require(ok.and.all(state==saved),'matching flux no state correction')
    state=saved;temperature=saved_t
    saved_flux=restricted
    trace_species=5+minloc(saved(6:,1,2,2),dim=1)
    call require(saved(trace_species,1,2,2)==0,'zero tracer fixture')
    restricted%x=0;restricted%y=0;restricted%z=0
    restricted%x(trace_species,1,2,2)=1e-19_dp
    call reflux_reactive_eb_state_3d(species,coarse,refined,cf,restricted,state,temperature,ok)
    call require(.not.ok.and.all(state==saved).and.all(temperature==saved_t), &
      'trace-negative reflux rejected transactionally without clipping')
    call compute_reactive_eb_reflux_increment_3d(coarse,refined,cf,restricted,increment,ok,touched=touched)
    call require(ok.and.increment(trace_species,1,2,2)<0,'raw reflux retains debt before donor reconciliation')
    unrefined=.false.
    call compute_reactive_eb_reflux_increment_3d(coarse,unrefined,cf,restricted,increment,ok, &
      all_faces=.true.,touched=touched)
    call require(ok.and.count(touched)==2,'all-face interior transfer touches exactly two exposed cells')
    volume=product([coarse%dx,coarse%dy,coarse%dz])
    call require(increment(trace_species,1,2,2)<0.and.increment(trace_species,2,2,2)>0, &
      'signed all-face transfer retains both donor and recipient increments')
    call require(abs(sum(increment(trace_species,:,:,:)*coarse%volume_fraction)*volume)<1e-34_dp, &
      'raw interior reflux transfer conserves its species inventory')
    saved_increment=increment;saved_touched=touched
    call compute_reactive_eb_reflux_increment_3d(coarse,refined(:2,:,:),cf,restricted,increment,ok,touched=touched)
    call require(.not.ok.and.all(increment==saved_increment).and.all(touched.eqv.saved_touched), &
      'raw reflux invalid mask retains increment and touched storage')
    restricted%x(trace_species,1,2,2)=ieee_value(0._dp,ieee_quiet_nan)
    call compute_reactive_eb_reflux_increment_3d(coarse,unrefined,cf,restricted,increment,ok,touched=touched)
    call require(.not.ok.and.all(increment==saved_increment).and.all(touched.eqv.saved_touched), &
      'raw reflux nonfinite flux retains increment and touched storage')

    restricted=saved_flux
    restricted%x=restricted%x*1.e12_dp
    restricted%y=restricted%y*1.e12_dp
    restricted%z=restricted%z*1.e12_dp
    call reflux_reactive_eb_state_3d(species,coarse,refined,cf,restricted,state,temperature,ok)
    call require(.not.ok.and.all(state==saved).and.all(temperature==saved_t),'nonphysical reflux rollback')
    restricted=saved_flux
    call restrict_eb_integrated_flux_3d(fine,0,ff,restricted,ok)
    call require(.not.ok.and.all(restricted%x==saved_flux%x).and.all(restricted%y==saved_flux%y).and. &
      all(restricted%z==saved_flux%z),'invalid restriction rollback')
    call reflux_reactive_eb_state_3d(species,coarse,refined(:2,:,:),cf,restricted,state,temperature,ok)
    call require(.not.ok.and.all(state==saved).and.all(temperature==saved_t),'bad mask rollback')
    temperature(1,2,2)=ieee_value(0._dp,ieee_quiet_nan)
    call reflux_reactive_eb_state_3d(species,coarse,refined,cf,restricted,state,temperature,ok)
    call require(.not.ok.and.all(state==saved).and.ieee_is_nan(temperature(1,2,2)), &
      'nonfinite interface temperature rejected without floating point trap')
    temperature=saved_t
    ff%x(1,0,1,1)=ieee_value(0._dp,ieee_quiet_nan)
    call restrict_eb_integrated_flux_3d(fine,ratio,ff,restricted,ok)
    call require(.not.ok.and.all(restricted%x==saved_flux%x),'nonfinite flux rollback')
    deallocate(ff%x,ff%y,ff%z)
  end do
  call check_exchange_interval()
  call check_interior_donors()
  call check_interface_budget()
  call check_endpoint_donor()
  call check_uniform_inventory_packet()
  call check_passive_budgets()
  call check_peer_flux()
  write(*,'(a)') 'test_eb_reactive_reflux_3d: PASS'
contains

  ! React.cpp packs old UEINT but forces it with the projected gas-energy
  ! change. An independently healthy gas/H endpoint can therefore still be
  ! inadmissible to the reactor. These are caloric/material admission checks,
  ! not a qualification of reaction rates or a full reacting native step.
  subroutine check_trace_energy_donors()
    type(eb_geometry_3d) :: metrics
    integer :: ng,ns,c,trace,inert,mode,attempt
    real(dp) :: y(size(species)),mw,rg,cp,cv,gamma,h,e,s,eh
    real(dp) :: target(reactive_nvar(size(species))+1),donor(size(target)),debt(size(target)),zero(size(target))
    real(dp) :: funded(size(target)),remaining(size(target)),packet(size(target)),post(size(target)),residual(size(target))
    real(dp) :: stock(size(target)),pv(5+size(species)),temp,sound,bias(2)
    real(dp), parameter :: target_volume=.125_dp,donor_volume=2._dp
    logical :: good,complete
    ns=size(species);ng=reactive_nvar(ns);trace=0;inert=0
    do c=1,ns
      if(species(c)%name=='H2') trace=c
      if(species(c)%name=='N2') inert=c
    enddo
    call require(trace>0.and.inert>0,'trace donor species identities')
    call build_axis_plane_eb_geometry_3d(1,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'z',-1._dp,metrics,good)
    call require(good,'trace donor regular metrics')
    y=0;y(inert)=1
    call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,h,e,s,good)
    call require(good,'trace target gas caloric state')
    call mixture_mass_properties(species,y,700._dp,mw,rg,cp,cv,gamma,h,eh,s,good)
    call require(good,'trace target independent held state')
    target=0;target(1)=1;target(5)=e;target(6:ng)=y;target(ng+1)=eh
    y(trace)=1e-6_dp;y(inert)=1-y(trace)
    call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,h,e,s,good)
    call require(good,'trace donor gas caloric state')
    call mixture_mass_properties(species,y,700._dp,mw,rg,cp,cv,gamma,h,eh,s,good)
    call require(good,'trace donor independent held state')
    donor=0;donor(1)=1;donor(5)=e;donor(6:ng)=y;donor(ng+1)=eh
    debt=0;debt(1)=-1e-21_dp;debt(5+trace)=-1e-21_dp;zero=0
    post=target+debt
    call recover_reactive_eb_cell_3d(species,metrics,[1,1,1],post(:ng),1000._dp,pv,temp,sound,good)
    call require(good.and.post(5+trace)<0,'EOS tolerance alone does not settle a signed trace debt')
    do mode=1,3
    do attempt=1,2
      stock=donor
      if(attempt==2) stock=target
      bias=[target(ng+1)-target(5),stock(ng+1)-stock(5)]
      funded=-17;remaining=-19;packet=-23
      select case(mode)
      case(1)
        call supply_reactive_eb_gas_eos_donor_3d(species,metrics,[1,1,1],1000._dp,metrics,[1,1,1],1000._dp, &
          target(:ng),debt(:ng),target_volume,stock(:ng),zero(:ng),donor_volume, &
          funded(:ng),remaining(:ng),packet(:ng),complete,good)
        if(good) then
          funded(ng+1)=target(ng+1);remaining(ng+1)=stock(ng+1);packet(ng+1)=0
        endif
      case(2)
        call supply_reactive_eb_gas_eos_donor_3d(species,metrics,[1,1,1],1000._dp,metrics,[1,1,1],1000._dp, &
          target,debt,target_volume,stock,zero,donor_volume,funded,remaining,packet,complete,good, &
          reactor_energy_bias=bias)
      case(3)
        call supply_reactive_eb_held_eos_donor_3d(species,metrics,[1,1,1],1000._dp,metrics,[1,1,1],1000._dp, &
          target,debt,target_volume,stock,zero,donor_volume,funded,remaining,packet,complete,good)
      end select
      if(attempt==2) then
        call require(.not.good.and..not.complete.and.all(funded==-17).and.all(remaining==-19).and.all(packet==-23), &
          'missing trace inventory refuses without clipping or publishing outputs')
        cycle
      endif
      call require(good.and.complete,'actual material donor settles an EOS-admissible negative trace')
      post=funded+debt
      call require(all(post(6:ng)>=0).and.all(remaining(6:ng)>=0).and.packet(5+trace)>0, &
        'trace settlement publishes a nonnegative state and a real material packet')
      residual=target_volume*(funded-target)+donor_volume*(remaining-donor)
      call require(maxval(abs(residual)/max(1._dp,abs(target_volume*target)+abs(donor_volume*donor)))< &
        64*epsilon(1._dp),'trace donor conserves every transported component')
      call require(maxval(abs(target_volume*(funded-target)-packet)/max(1._dp,abs(target_volume*target)))< &
        64*epsilon(1._dp),'trace target increment agrees with the actual face packet')
      call require(maxval(abs(donor_volume*(remaining-donor)+packet)/max(1._dp,abs(donor_volume*donor)))< &
        64*epsilon(1._dp),'trace donor debit agrees with the opposite face packet')
    enddo
    enddo
    print *, 'PASS: gas, held and reactor-energy donors settle signed trace debt conservatively'
  end subroutine check_trace_energy_donors

  subroutine check_extrapolated_energy_donors()
    type(nasa7_species),allocatable::original_eos_species(:),eos_species(:)
    type(eb_geometry_3d)::g
    real(dp),allocatable::gas(:),donor(:),debt(:),zero(:),a(:),b(:),packet(:),y(:),pv(:),post(:)
    real(dp)::requested,mw,rg,cp,cv,gamma,h,e,s,t,sound,residual
    integer::ng,sp,case_id,mode,passed
    logical::ok,complete
    original_eos_species=species
    eos_species=original_eos_species
    call configure_nasa7_eos_policy(eos_species,.true.,ok,200._dp,3500._dp)
    if(.not.ok) error stop 'policy'
    ng=reactive_nvar(size(eos_species))
    allocate(gas(ng+1),donor(ng+1),debt(ng+1),zero(ng+1),a(ng+1),b(ng+1),packet(ng+1), &
      y(size(eos_species)),pv(5+size(eos_species)),post(ng+1))
    y=0
    do sp=1,size(eos_species)
      if(eos_species(sp)%name=='N2') y(sp)=1
    enddo
    call build_axis_plane_eb_geometry_3d(1,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',-1._dp,g,ok)
    if(.not.ok) error stop 'geometry'
    passed=0
    do case_id=1,6
      select case(case_id)
      case(1);requested=150
      case(2);requested=190
      case(3);requested=3600
      case(4);requested=5000
      case(5);requested=300
      case(6);requested=2200
      end select
      if(case_id==5) then
        call configure_nasa7_eos_policy(eos_species,.true.,ok,500._dp,1500._dp)
        call require(ok,'explicit donor EOS inverse limits')
      endif
      call mixture_mass_properties(eos_species,y,requested,mw,rg,cp,cv,gamma,h,e,s,ok)
      if(.not.ok) error stop 'active initial properties'
      gas=0;gas(1)=1;gas(5)=e;gas(6:ng)=y;gas(ng+1)=e
      donor=gas;debt=-1.5_dp*gas;zero=0
      do mode=1,3
        a=-17;b=-19;packet=-23
        select case(mode)
        case(1)
          call supply_reactive_eb_gas_eos_donor_3d(eos_species,g,[1,1,1],requested,g,[1,1,1],requested, &
            gas(:ng),debt(:ng),1._dp,donor(:ng),zero(:ng),2._dp,a(:ng),b(:ng),packet(:ng),complete,ok)
        case(2)
          call supply_reactive_eb_gas_eos_donor_3d(eos_species,g,[1,1,1],requested,g,[1,1,1],requested, &
            gas,debt,1._dp,donor,zero,2._dp,a,b,packet,complete,ok,reactor_energy_bias=[0._dp,0._dp])
        case(3)
          debt=0;debt(ng+1)=-2*abs(gas(ng+1))-1e6_dp
          call supply_reactive_eb_held_eos_donor_3d(eos_species,g,[1,1,1],requested,g,[1,1,1],requested, &
            gas,debt,1._dp,donor,zero,20._dp,a,b,packet,complete,ok)
        end select
        print *, 'case/mode/T/ok/complete=',case_id,mode,requested,ok,complete
        if(.not.ok.or..not.complete) error stop 'valid extrapolated donor refused'
        if(mode==1) then
          residual=maxval(abs((a(:ng)-gas(:ng))+2*(b(:ng)-donor(:ng)))/max(1._dp,abs(gas(:ng))))
        else
          residual=maxval(abs((a-gas)+merge(20._dp,2._dp,mode==3)*(b-donor))/max(1._dp,abs(gas)))
        endif
        if(residual>3e-12_dp) error stop 'packet conservation'
        post=a+debt
        call recover_reactive_eb_cell_3d(eos_species,g,[1,1,1],post(:ng),requested,pv,t,sound,ok)
        if(.not.ok) error stop 'post gas EOS'
        if(mode/=1) then
          call recover_reactive_eb_cell_3d(eos_species,g,[1,1,1],post(:ng),requested,pv,t,sound,ok, &
            internal_energy_density=post(ng+1))
          if(.not.ok) error stop 'post actual H EOS'
        endif
        passed=passed+1
        debt=-1.5_dp*gas
      enddo
    enddo
    call require(passed==18,'extrapolated gas/held/reactor donor cases')
    ! Removing an artificial fitted ceiling must not admit invalid heat capacity
    ! or enable extrapolation on the original strict table. Refused requests
    ! retain every published target, donor and packet component.
    call configure_nasa7_eos_policy(eos_species,.true.,ok,200._dp,3500._dp)
    call require(ok,'restore active donor limits')
    do case_id=1,2
      requested=150
      call mixture_mass_properties(eos_species,y,requested,mw,rg,cp,cv,gamma,h,e,s,ok)
      call require(ok,'refusal reference active properties')
      gas=0;gas(1)=1;gas(5)=e;gas(6:ng)=y;gas(ng+1)=e
      donor=gas;debt=-1.5_dp*gas;zero=0
      if(case_id==1) then
        call mixture_mass_properties(eos_species,y,3500._dp,mw,rg,cp,cv,gamma,h,e,s,ok)
        call require(ok,'heat capacity refusal upper edge')
        donor(5)=e+cv*(20000._dp-3500._dp)
        donor(ng+1)=donor(5)
        call recover_reactive_eb_cell_3d(eos_species,g,[1,1,1],donor(:ng),requested,pv,t,sound,ok)
        call require(.not.ok,'positive extrapolated temperature with invalid Cv refused')
      else
        eos_species=original_eos_species
      endif
      do mode=1,3
        a=-17;b=-19;packet=-23
        select case(mode)
        case(1)
          call supply_reactive_eb_gas_eos_donor_3d(eos_species,g,[1,1,1],requested,g,[1,1,1],requested, &
            gas(:ng),debt(:ng),1._dp,donor(:ng),zero(:ng),2._dp,a(:ng),b(:ng),packet(:ng),complete,ok)
        case(2)
          call supply_reactive_eb_gas_eos_donor_3d(eos_species,g,[1,1,1],requested,g,[1,1,1],requested, &
            gas,debt,1._dp,donor,zero,2._dp,a,b,packet,complete,ok,reactor_energy_bias=[0._dp,0._dp])
        case(3)
          debt=0;debt(ng+1)=-2*abs(gas(ng+1))-1e6_dp
          call supply_reactive_eb_held_eos_donor_3d(eos_species,g,[1,1,1],requested,g,[1,1,1],requested, &
            gas,debt,1._dp,donor,zero,20._dp,a,b,packet,complete,ok)
        end select
        call require(.not.ok.and..not.complete,'invalid heat capacity or strict donor refused')
        call require(all(a==-17).and.all(b==-19).and.all(packet==-23),'donor refusal preserves all outputs')
      enddo
    enddo

  end subroutine check_extrapolated_energy_donors

  subroutine check_reactor_energy_donors()
    type(eb_geometry_3d) :: metrics
    integer :: ng,ns,c,case_id,loans
    real(dp) :: composition(size(species)),mw,rg,cp,cv,gamma,h,e,entropy,eg,eh,ec
    real(dp) :: original(reactive_nvar(size(species))+1),donor(reactive_nvar(size(species))+1)
    real(dp) :: target(size(original)),debt(size(original)),donor_debt(size(original)),post(size(original))
    real(dp) :: next_target(size(original)),next_donor(size(original)),packet(size(original)),residual(size(original))
    real(dp) :: saved_target(size(original)),saved_donor(size(original)),saved_packet(size(original)),bias(2)
    real(dp) :: pv(5+size(species)),temp,sound
    logical :: good,complete
    ns=size(species);ng=reactive_nvar(ns);composition=0
    do c=1,ns
      if(species(c)%name=='N2') composition(c)=1
    enddo
    call require(sum(composition)==1,'reactor donor inert composition')
    call build_axis_plane_eb_geometry_3d(1,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'z',-1._dp,metrics,good)
    call require(good,'reactor donor regular metrics')
    call mixture_mass_properties(species,composition,1000._dp,mw,rg,cp,cv,gamma,h,eg,entropy,good)
    call require(good,'reactor donor gas caloric state')
    call mixture_mass_properties(species,composition,700._dp,mw,rg,cp,cv,gamma,h,eh,entropy,good)
    call require(good,'reactor donor independent held state')
    call mixture_mass_properties(species,composition,400._dp,mw,rg,cp,cv,gamma,h,ec,entropy,good)
    call require(good,'reactor donor forced caloric state')
    original=0;original(1)=1;original(5)=eg;original(6:ng)=composition;original(ng+1)=eh
    bias=eh-eg;donor_debt=0
    ! All three predicates and actual complete-packet conservation must hold
    ! after a loan; merely repairing density is insufficient.
    do case_id=1,3
      target=original;donor=original;debt=0
      if(case_id==1) then
        debt(5)=ec-eg
        post=target+debt
        call recover_reactive_eb_cell_3d(species,metrics,[1,1,1],post(:ng),1000._dp,pv,temp,sound,good)
        call require(good,'forced gas alone admits')
        call recover_reactive_eb_cell_3d(species,metrics,[1,1,1],post(:ng),700._dp,pv,temp,sound,good, &
          internal_energy_density=post(ng+1))
        call require(good,'forced actual UEINT alone admits')
        post(5)=post(5)+bias(1)
        call recover_reactive_eb_cell_3d(species,metrics,[1,1,1],post(:ng),700._dp,pv,temp,sound,good)
        call require(.not.good,'forced reactor budget independently refuses')
      else
        debt=-(real(case_id,dp)-.5_dp)*original
      endif
      residual=0;complete=.false.
      do loans=1,6
        next_target=-17;next_donor=-19;packet=-23
        call supply_reactive_eb_gas_eos_donor_3d(species,metrics,[1,1,1],1000._dp,metrics,[1,1,1],1000._dp, &
          target,debt,1._dp,donor,donor_debt,1._dp,next_target,next_donor,packet,complete,good, &
          reactor_energy_bias=bias)
        call require(good,'joint reactor/gas/held donor proposal')
        residual=residual+(next_target-target)+(next_donor-donor)
        call require(maxval(abs(residual)/max(1._dp,abs(original)))<3.e-12_dp, &
          'complete gas and independent UEINT packet conservation')
        call require(maxval(abs(packet-packet(1)*original)/max(1._dp,abs(original)))<3.e-12_dp, &
          'reactor offset never enters the material packet')
        call check_joint_endpoint(next_donor,bias(2))
        call check_joint_endpoint(next_target,bias(1))
        if(case_id==3.and.loans==1) then
          call require(.not.complete.and.next_target(1)+debt(1)<0, &
            'partial loan progresses while pending density remains negative')
        endif
        if(case_id==2.and.loans==1) call require(.not.complete,'reactor budget requires multiple healthy donors')
        target=next_target
        if(complete) exit
      enddo
      call require(complete,'pooled donors admit the complete pending endpoint')
      call check_joint_endpoint(target+debt,bias(1))
    enddo
    ! A donor's own pending gas can remain positive and gas/H-admissible but
    ! exhaust its reactor energy. That provider must not publish any packet.
    target=original;debt=0;debt(5)=ec-eg;donor=original;donor_debt=0;donor_debt(5)=eg-ec
    next_target=-17;next_donor=-19;packet=-23
    saved_target=next_target;saved_donor=next_donor;saved_packet=packet
    call supply_reactive_eb_gas_eos_donor_3d(species,metrics,[1,1,1],1000._dp,metrics,[1,1,1],1000._dp, &
      target,debt,1._dp,donor,donor_debt,1._dp,next_target,next_donor,packet,complete,good, &
      reactor_energy_bias=bias)
    call require(.not.good.and..not.complete,'pending donor reactor budget refusal')
    call require(all(next_target==saved_target).and.all(next_donor==saved_donor).and.all(packet==saved_packet), &
      'pending donor refusal leaves every output unchanged')
    donor_debt=0;bias(2)=ieee_value(0._dp,ieee_quiet_nan)
    call supply_reactive_eb_gas_eos_donor_3d(species,metrics,[1,1,1],1000._dp,metrics,[1,1,1],1000._dp, &
      target,debt,1._dp,donor,donor_debt,1._dp,next_target,next_donor,packet,complete,good, &
      reactor_energy_bias=bias)
    call require(.not.good.and..not.complete,'nonfinite reactor bias refusal')
    call require(all(next_target==saved_target).and.all(next_donor==saved_donor).and.all(packet==saved_packet), &
      'nonfinite reactor bias leaves every output unchanged')
  end subroutine check_reactor_energy_donors

  subroutine check_joint_endpoint(value,bias)
    real(dp), intent(in) :: value(:),bias
    type(eb_geometry_3d) :: metrics
    real(dp) :: gas(reactive_nvar(size(species))),pv(5+size(species)),temp,sound
    logical :: good
    integer :: ng
    ng=reactive_nvar(size(species))
    call build_axis_plane_eb_geometry_3d(1,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'z',-1._dp,metrics,good)
    call require(good,'joint endpoint geometry')
    gas=value(:ng)
    call recover_reactive_eb_cell_3d(species,metrics,[1,1,1],gas,1000._dp,pv,temp,sound,good)
    call require(good,'joint endpoint projected gas EOS')
    call recover_reactive_eb_cell_3d(species,metrics,[1,1,1],gas,700._dp,pv,temp,sound,good, &
      internal_energy_density=value(ng+1))
    call require(good,'joint endpoint independent actual UEINT EOS')
    gas(5)=gas(5)+bias
    call recover_reactive_eb_cell_3d(species,metrics,[1,1,1],gas,700._dp,pv,temp,sound,good)
    call require(good,'joint endpoint forced-reactor energy EOS')
  end subroutine check_joint_endpoint
  subroutine check_tree_rereflux()
    type(amr_patch_3d) :: patches(2)
    type(amr_patch_tree_3d) :: topology
    type(eb_geometry_3d) :: fg,cg,mg
    type(reactive_eb_3d_config) :: cc
    type(eb_reactive_tree_3d) :: tree,before_tree,gas_scalar_reference
    type(amr_tree_state_3d) :: correction,result,plain,old_result,scalar,scalar_keep,scalar_correction
    type(eb_integrated_flux_3d), allocatable :: flux(:),replacement(:),plain_flux(:),old_flux(:), &
      scalar_flux(:),scalar_flux_keep(:),gas_scalar_flux(:)
    real(dp), allocatable :: cs(:,:,:,:),ct(:,:,:),density(:)
    real(dp) :: before(reactive_nvar(size(species))),after(size(before)),expected(size(before)), &
      cv,fv,err
    type(eb_geometry_3d) :: padded,physical(3)
    type(eb_geometry_3d), allocatable :: wide(:)
    type(state_redist_register_3d) :: registers(3),kept_registers(3),scalar_registers(3)
    integer :: nv,p,c
    logical :: good
    nv=reactive_nvar(size(species));cc%nx=4;cc%ny=4;cc%nz=4
    call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.37_dp,fg,good)
    call require(good,'tree rereflux finest geometry')
    call coarsen_eb_geometry_3d(fg,2,mg,good)
    call require(good,'tree rereflux intermediate geometry')
    call coarsen_eb_geometry_3d(mg,2,cg,good)
    call require(good,'tree rereflux coarse geometry')
    call initialize_amr_patch_3d(4,4,4,2,3,2,3,2,3,2,patches(1),good)
    call require(good,'tree rereflux child bounds')
    call initialize_amr_patch_3d(4,4,4,3,4,1,2,1,2,2,patches(2),good)
    call require(good,'tree rereflux grandchild bounds')
    call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,topology,good)
    call require(good,'tree rereflux topology')
    allocate(cs(nv,4,4,4),ct(4,4,4))
    call initialize_reactive_eb_density_sheet_3d(species,cc,cg,cs,ct,good)
    call require(good,'tree rereflux root state')
    call initialize_reactive_eb_tree_3d(species,topology,fg,cs,ct,tree,good)
    call require(good,'tree rereflux initialized tree')
    allocate(flux(3),correction%fields(3))
    do p=1,3
      allocate(flux(p)%x(nv,0:4,4,4),flux(p)%y(nv,4,0:4,4),flux(p)%z(nv,4,4,0:4))
      flux(p)%x=0;flux(p)%y=0;flux(p)%z=0
      allocate(correction%fields(p)%values(nv,4,4,4),source=0._dp)
    enddo
    allocate(scalar%fields(3),scalar_flux(3))
    do p=1,3
      allocate(scalar%fields(p)%values(1,4,4,4))
      scalar%fields(p)%values(1,:,:,:)=tree%nodes(p)%state(5,:,:,:)
      allocate(scalar_flux(p)%x(1,0:4,4,4),scalar_flux(p)%y(1,4,0:4,4),scalar_flux(p)%z(1,4,4,0:4))
      scalar_flux(p)%x=0;scalar_flux(p)%y=0;scalar_flux(p)%z=0
    enddo
    gas_scalar_reference=tree;gas_scalar_flux=flux
    scalar_flux(3)%x(1,0,1,1)=1.e-12_dp*max(1._dp,abs(scalar%fields(3)%values(1,1,1,1)))
    gas_scalar_flux(3)%x(5,0,1,1)=scalar_flux(3)%x(1,0,1,1)
    call synchronize_reactive_eb_tree_flux_3d(species,gas_scalar_reference,gas_scalar_flux,good)
    call require(good,'scalar reflux gas energy control')
    call synchronize_reactive_eb_tree_scalar_flux_3d(species,tree,scalar,scalar_flux,good)
    call require(good,'scalar reflux accepted')
    call require(any(scalar%fields(1)%values(1,:,:,:)/=tree%nodes(1)%state(5,:,:,:)), &
      'scalar reflux changes the exposed coarse field')
    do p=1,3
      call require(maxval(abs(scalar%fields(p)%values(1,:,:,:)-gas_scalar_reference%nodes(p)%state(5,:,:,:)))/ &
        max(1._dp,maxval(abs(gas_scalar_reference%nodes(p)%state(5,:,:,:))))<1.e-12_dp, &
        'scalar reflux matches gas energy face synchronization')
      call require(all(scalar_flux(p)%x(1,:,:,:)==gas_scalar_flux(p)%x(5,:,:,:)).and. &
        all(scalar_flux(p)%y(1,:,:,:)==gas_scalar_flux(p)%y(5,:,:,:)).and. &
        all(scalar_flux(p)%z(1,:,:,:)==gas_scalar_flux(p)%z(5,:,:,:)), &
        'scalar reflux effective faces match gas energy')
    enddo
    scalar_keep=scalar;scalar_flux_keep=scalar_flux
    scalar_flux(3)%x(1,0,1,1)=ieee_value(0._dp,ieee_quiet_nan)
    call synchronize_reactive_eb_tree_scalar_flux_3d(species,tree,scalar,scalar_flux,good)
    call require(.not.good,'scalar reflux invalid face rejected')
    do p=1,3
      call require(all(scalar%fields(p)%values==scalar_keep%fields(p)%values), &
        'scalar reflux state rollback')
    enddo
    scalar_flux=scalar_flux_keep
    density=1.e-7_dp*tree%nodes(1)%state(:,3,2,2)
    correction%fields(1)%values(:,2,1,1)=density
    correction%fields(1)%values(:,2,2,2)=2*density
    correction%fields(1)%values(:,3,2,2)=3*density
    correction%fields(2)%values(:,3,1,1)=4*density
    before_tree=tree
    call assemble_reactive_eb_tree_reflux_3d(species,tree,flux,result,replacement,good,rereflux_delta=correction)
    call require(good,'tree rereflux assembled before EOS')
    call require(all(result%fields(1)%values(:,2,1,1)==density),'exposed root receives its density correction')
    call require(all(result%fields(1)%values(:,2:3,2:3,2:3)==0),'covered root receives no duplicate correction')
    call require(all(result%fields(2)%values(:,1,1,1)==2*density),'covered cut root forwards density to fine')
    call require(all(result%fields(2)%values(:,3:4,1:2,1:2)==0),'covered child forwards to descendant')
    err=maxval(abs(result%fields(3)%values(:,1,1,1)-7*density)/max(1._dp,abs(density)))
    call require(err<16*epsilon(1._dp)*maxval(abs(density)),'nested corrections sum once at descendant')
    call require(all(result%fields(3)%values(:,3,3,3)==3*density),'ancestor correction reaches entire grandchild')
    do p=1,3
      call require(all(tree%nodes(p)%state==before_tree%nodes(p)%state),'assembly leaves input tree unchanged')
    enddo
    call integrate_reactive_eb_tree_3d(tree,before,good)
    call require(good,'tree rereflux initial inventory')
    cv=cg%dx*cg%dy*cg%dz;fv=tree%nodes(2)%geometry%dx**3
    do c=1,nv
      expected(c)=cv*sum(cg%volume_fraction*correction%fields(1)%values(c,:,:,:))+ &
        fv*sum(tree%nodes(2)%geometry%volume_fraction*correction%fields(2)%values(c,:,:,:))
    enddo
    call synchronize_reactive_eb_tree_flux_3d(species,tree,flux,good,rereflux_delta=correction)
    call require(good,'tree rereflux joint state and flux acceptance')
    call integrate_reactive_eb_tree_3d(tree,after,good)
    call require(good,'tree rereflux final inventory')
    call require(maxval(abs(after-before-expected)/max(abs(before),1._dp))<2.e-14_dp, &
      'nested cut-cell rereflux conserves every composite component')
    ! A nonzero fine interface flux and the register correction share the same
    ! signed assembly. Compare the increment with/without register input.
    tree=before_tree;flux(3)%x(:,0,1,1)=1.e-9_dp*tree%nodes(3)%state(:,1,1,1)
    call assemble_reactive_eb_tree_reflux_3d(species,tree,flux,plain,plain_flux,good)
    call require(good,'tree rereflux nonzero interface baseline')
    old_result=result
    call assemble_reactive_eb_tree_reflux_3d(species,tree,flux,result,replacement,good,rereflux_delta=correction)
    call require(good,'tree rereflux interface plus register assembly')
    do p=1,3
      err=maxval(abs(result%fields(p)%values-plain%fields(p)%values-old_result%fields(p)%values))
      call require(err<1.e-15_dp,'tree rereflux adds register to signed flux mismatch once')
    enddo
    old_flux=flux
    correction%fields(1)%values(:,3,2,2)=-2*tree%nodes(3)%state(:,1,1,1)
    correction%fields(2)%values=0
    call synchronize_reactive_eb_tree_flux_3d(species,tree,flux,good,rereflux_delta=correction)
    call require(.not.good,'descendant EOS rejection after provisional root acceptance')
    do p=1,3
      call require(all(tree%nodes(p)%state==before_tree%nodes(p)%state).and. &
        all(tree%nodes(p)%temperature==before_tree%nodes(p)%temperature),'tree rereflux whole-state rollback')
      call require(all(flux(p)%x==old_flux(p)%x).and.all(flux(p)%y==old_flux(p)%y).and. &
        all(flux(p)%z==old_flux(p)%z),'tree rereflux effective-flux rollback')
    enddo
    old_result=result
    correction%fields(1)%values=ieee_value(0._dp,ieee_quiet_nan)
    call assemble_reactive_eb_tree_reflux_3d(species,tree,flux,result,replacement,good,rereflux_delta=correction)
    call require(.not.good,'tree rereflux invalid register rejected')
    do p=1,3
      call require(all(result%fields(p)%values==old_result%fields(p)%values),'invalid register assembly rollback')
    enddo
    ! Actual interval ledger assembly -> EB re-reflux -> fine injection ->
    ! accepted gas/flux and one-time ledger consumption on a nested tree.
    tree=before_tree
    do p=1,3
      flux(p)%x=0;flux(p)%y=0;flux(p)%z=0
      physical(p)=tree%nodes(p)%geometry
      call initialize_state_redist_tree_register_3d(topology,p,nv,6,registers(p),good,scale=.5_dp)
      call require(good,'interval close tree register')
    enddo
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'x',.37_dp,padded,good)
    call require(good,'interval close padded geometry')
    call build_eb_tree_extended_geometry_3d(topology,physical,padded,wide,good)
    call require(good,'interval close prepared geometry')
    ! The physical-domain cut-cell correction now redistributes its complete
    ! signed source among finite-domain recipients before closing the interval.
    density=1.e-7_dp*tree%nodes(1)%state(:,2,1,2)
    registers(1)%coarse_delta(:,8,7,8)=density
    do p=1,3
      correction%fields(p)%values=0._dp
    enddo
    kept_registers=registers
    call assemble_state_redist_interface_rereflux_3d(topology,wide,registers,nv,0,correction,good)
    call require(good,'interval close retains physical-boundary register inventory')
    do c=1,nv
      expected(c)=cv*sum(cg%volume_fraction*correction%fields(1)%values(c,:,:,:))
    enddo
    call require(maxval(abs(expected-cv*cg%volume_fraction(2,1,2)*density)/ &
      max(tiny(1._dp),abs(cv*cg%volume_fraction(2,1,2)*density)))<128*epsilon(1._dp), &
      'physical-boundary register conserves every source component')
    do p=1,3
      call require(all(tree%nodes(p)%state==before_tree%nodes(p)%state), &
        'physical-boundary assembly retains input gas')
      if(allocated(registers(p)%coarse_delta)) call require( &
        all(registers(p)%coarse_delta==kept_registers(p)%coarse_delta),'physical-boundary assembly retains coarse debt')
      if(allocated(registers(p)%fine_delta)) call require( &
        all(registers(p)%fine_delta==kept_registers(p)%fine_delta),'physical-boundary assembly retains fine debt')
    enddo
    registers(1)%coarse_delta(:,8,7,8)=0._dp
    ! Exercise a separate open x-side interval with fine-register injection.
    density=1.e-7_dp*tree%nodes(1)%state(:,4,2,2)
    registers(1)%coarse_delta(:,10,8,8)=density
    registers(2)%fine_delta(:,11,7,7)=.2_dp*density
    registers(2)%coarse_delta(:,9,7,7)=.3_dp*density
    registers(3)%fine_delta(:,6,7,7)=.4_dp*density
    do p=1,3
      scalar%fields(p)%values(1,:,:,:)=tree%nodes(p)%state(5,:,:,:)
      scalar_flux(p)%x=0;scalar_flux(p)%y=0;scalar_flux(p)%z=0
      call initialize_state_redist_tree_register_3d(topology,p,1,6,scalar_registers(p),good,scale=.5_dp)
      call require(good,'scalar interval register initialization')
      if(allocated(registers(p)%coarse_delta)) &
        scalar_registers(p)%coarse_delta=registers(p)%coarse_delta(5:5,:,:,:)
      if(allocated(registers(p)%fine_delta)) &
        scalar_registers(p)%fine_delta=registers(p)%fine_delta(5:5,:,:,:)
    enddo
    kept_registers=registers
    call assemble_state_redist_interface_rereflux_3d(topology,wide,registers,nv,0,correction,good)
    call require(good,'interval close assembled EB correction')
    call assemble_state_redist_interface_rereflux_3d(topology,wide,scalar_registers,1,0,scalar_correction,good)
    call require(good,'interval close assembled scalar correction')
    do p=1,3
      call require(all(scalar_correction%fields(p)%values(1,:,:,:)==correction%fields(p)%values(5,:,:,:)), &
        'scalar and gas energy register assembly agree')
    enddo
    do c=1,nv
      expected(c)=cg%dx*cg%dy*cg%dz*sum(cg%volume_fraction*correction%fields(1)%values(c,:,:,:))
    enddo
    call integrate_reactive_eb_tree_3d(tree,before,good)
    call require(good,'interval close initial inventory')
    call synchronize_reactive_eb_tree_register_flux_3d(species,tree,flux,wide,registers,0,good)
    call require(good,'interval close gas flux and registers accepted')
    call synchronize_reactive_eb_tree_scalar_register_flux_3d(species,before_tree,scalar,scalar_flux, &
      wide,scalar_registers,0,good)
    call require(good,'interval close scalar flux and registers accepted')
    do p=1,3
      call require(maxval(abs(scalar%fields(p)%values(1,:,:,:)-tree%nodes(p)%state(5,:,:,:)))/ &
        max(1._dp,maxval(abs(tree%nodes(p)%state(5,:,:,:))))<1.e-12_dp, &
        'scalar interval register matches gas energy transfer')
    enddo
    call require(all(scalar_registers(1)%coarse_delta==0).and.all(scalar_registers(2)%fine_delta==0), &
      'scalar interval consumes accepted interface exactly once')
    scalar_keep=scalar;scalar_flux_keep=scalar_flux
    scalar_registers(1)%coarse_delta(1,10,8,8)=ieee_value(0._dp,ieee_quiet_nan)
    call synchronize_reactive_eb_tree_scalar_register_flux_3d(species,before_tree,scalar,scalar_flux, &
      wide,scalar_registers,0,good)
    call require(.not.good,'scalar interval invalid register rejected')
    do p=1,3
      call require(all(scalar%fields(p)%values==scalar_keep%fields(p)%values), &
        'scalar interval rejected register preserves state')
      call require(all(scalar_flux(p)%x==scalar_flux_keep(p)%x).and. &
        all(scalar_flux(p)%y==scalar_flux_keep(p)%y).and. &
        all(scalar_flux(p)%z==scalar_flux_keep(p)%z), &
        'scalar interval rejected register preserves flux')
    enddo
    call require(ieee_is_nan(scalar_registers(1)%coarse_delta(1,10,8,8)), &
      'scalar interval rejected register remains unconsumed')
    call integrate_reactive_eb_tree_3d(tree,after,good)
    call require(good.and.maxval(abs(after-before-expected)/max(abs(before),1._dp))<2.e-14_dp, &
      'interval close composite ledger conservation after nested injection')
    call require(all(registers(1)%coarse_delta==0).and.all(registers(2)%fine_delta==0), &
      'interval close consumes accepted interface exactly once')
    call require(all(registers(2)%coarse_delta==kept_registers(2)%coarse_delta).and. &
      all(registers(3)%fine_delta==kept_registers(3)%fine_delta),'interval close preserves next interface ledgers')
    before_tree=tree;old_flux=flux
    registers(1)%coarse_delta(:,9,7,8)=-100._dp*tree%nodes(1)%state(:,3,1,2)
    kept_registers=registers
    call synchronize_reactive_eb_tree_register_flux_3d(species,tree,flux,wide,registers,0,good)
    call require(.not.good,'interval close rejects nonphysical gas correction')
    do p=1,3
      call require(all(tree%nodes(p)%state==before_tree%nodes(p)%state).and. &
        all(tree%nodes(p)%temperature==before_tree%nodes(p)%temperature),'interval close gas EOS rollback')
      call require(all(flux(p)%x==old_flux(p)%x).and.all(flux(p)%y==old_flux(p)%y).and. &
        all(flux(p)%z==old_flux(p)%z),'interval close flux rollback')
      if(allocated(registers(p)%coarse_delta)) then
        call require(all(registers(p)%coarse_delta==kept_registers(p)%coarse_delta),'interval close coarse ledger rollback')
      endif
      if(allocated(registers(p)%fine_delta)) then
        call require(all(registers(p)%fine_delta==kept_registers(p)%fine_delta),'interval close fine ledger rollback')
      endif
    enddo
  end subroutine check_tree_rereflux
  subroutine check_rereflux_pair()
    type(eb_geometry_3d) :: cg,fg
    type(reactive_eb_3d_config) :: cc,fc
    real(dp), allocatable :: cs(:,:,:,:),ct(:,:,:),fs(:,:,:,:),ft(:,:,:),increment(:,:,:,:), &
      old_cs(:,:,:,:),old_ct(:,:,:),old_fs(:,:,:,:),old_ft(:,:,:)
    real(dp) :: before(reactive_nvar(size(species))),after(reactive_nvar(size(species))),expected(size(before)),cv,fv
    logical :: inside(4,3,3),good
    integer :: nv,c
    nv=reactive_nvar(size(species));cc%nx=4;cc%ny=3;cc%nz=3;fc%nx=2;fc%ny=3;fc%nz=2
    call build_axis_plane_eb_geometry_3d(4,3,3,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.37_dp,cg,good)
    call require(good,'rereflux pair coarse geometry')
    call build_axis_plane_eb_geometry_3d(2,3,2,.5_dp,.75_dp,1._dp/3,2._dp/3,1._dp/3,2._dp/3,'x',.37_dp,fg,good)
    call require(good,'rereflux pair fine geometry')
    allocate(cs(nv,4,3,3),ct(4,3,3),fs(nv,2,3,2),ft(2,3,2),increment(nv,4,3,3))
    call initialize_reactive_eb_density_sheet_3d(species,cc,cg,cs,ct,good)
    call require(good,'rereflux pair coarse state')
    call initialize_reactive_eb_density_sheet_3d(species,fc,fg,fs,ft,good)
    call require(good,'rereflux pair fine state')
    inside=.false.;inside(3,2,2)=.true.;increment=0
    increment(:,2,2,2)=1.e-7_dp*cs(:,2,2,2);increment(:,3,2,2)=2.e-7_dp*cs(:,3,2,2)
    cv=cg%dx*cg%dy*cg%dz;fv=fg%dx*fg%dy*fg%dz
    do c=1,nv
      before(c)=cv*sum(cg%volume_fraction*cs(c,:,:,:),mask=.not.inside)+fv*sum(fg%volume_fraction*fs(c,:,:,:))
      expected(c)=cv*sum(cg%volume_fraction*increment(c,:,:,:))
    enddo
    call apply_reactive_eb_rereflux_pair_3d(species,cg,fg,increment,inside,[2,3,2],[0,0,0],[4,3,2], &
      cs,ct,fs,ft,good)
    call require(good,'joint coarse/fine rereflux EOS acceptance')
    do c=1,nv
      after(c)=cv*sum(cg%volume_fraction*cs(c,:,:,:),mask=.not.inside)+fv*sum(fg%volume_fraction*fs(c,:,:,:))
    enddo
    call require(all(abs(after-before-expected)<2.e-13_dp*max(1._dp,abs(before))), &
      'rereflux pair composite component inventory')
    old_cs=cs;old_ct=ct;old_fs=fs;old_ft=ft;increment=0
    increment(:,2,2,2)=1.e-7_dp*cs(:,2,2,2)
    increment(:,3,2,2)=-2*fs(:,1,1,1)
    call apply_reactive_eb_rereflux_pair_3d(species,cg,fg,increment,inside,[2,3,2],[0,0,0],[4,3,2], &
      cs,ct,fs,ft,good)
    call require(.not.good.and.all(cs==old_cs).and.all(ct==old_ct).and.all(fs==old_fs).and.all(ft==old_ft), &
      'fine EOS rejection rolls back accepted provisional coarse update')
  end subroutine check_rereflux_pair
  subroutine check_physical_domain_rereflux()
    type(eb_geometry_3d) :: padded
    type(state_redist_register_3d) :: reg
    type(state_redist_box_3d) :: box,bad_box
    real(dp) :: source(3,16,16,16),vector(3),mass(3),kappa
    real(dp), allocatable :: increment(:,:,:,:),reference(:,:,:,:),previous(:,:,:,:),halo(:,:,:,:),kept_halo(:,:,:,:)
    integer :: flags(16,16,16),cell(3),scenario,c
    logical :: owned(16,16,16),good
    call build_axis_plane_eb_geometry_3d(16,16,16,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'x',.37_dp,padded,good)
    call require(good,'physical register padded geometry')
    owned=.false.;owned(7:10,7:10,7:10)=.true.
    box%valid_lower=7;box%valid_upper=10;box%domain_lower=7;box%domain_upper=10;box%periodic=.false.
    vector=[2._dp,-3._dp,4._dp]
    do scenario=0,2
      cell=[8,8,8]
      if(scenario==1) cell(2)=7
      if(scenario==2) cell(2:3)=7
      flags=0;flags(cell(1),cell(2),cell(3))=srd_coarse_boundary
      call initialize_state_redist_register_3d(3,[16,16,16],reg,good,coarse_flags=flags,valid_cells=owned)
      call require(good,'physical register flags')
      source=0;source(:,cell(1),cell(2),cell(3))=vector
      kappa=padded%volume_fraction(cell(1),cell(2),cell(3))
      call assemble_state_redist_patch_rereflux_3d(padded,reg,[4,4,4],source,reference,good)
      call require(good,'physical register unbounded kernel')
      call assemble_state_redist_patch_rereflux_3d(padded,reg,[4,4,4],source,increment,good,box=box)
      call require(good,'physical register bounded kernel')
      do c=1,3
        mass(c)=sum(padded%volume_fraction(7:10,7:10,7:10)*increment(c,:,:,:))
      enddo
      call require(maxval(abs(mass-kappa*vector))<=64*epsilon(1._dp)*maxval(abs(kappa*vector)), &
        'physical register preserves complete signed component inventory')
      if(scenario==0) call require(all(increment==reference),'interior register retains exact kernel output')
      if(scenario==1) then
        ! The pinned AMReX geometry/kernel probe independently gives this ratio.
        mass(1)=sum(padded%volume_fraction(7:10,7:10,7:10)*reference(1,:,:,:))/(kappa*vector(1))
        call require(abs(mass(1)-.83367781155015186_dp)<=64*epsilon(1._dp), &
          'physical register unbounded kernel reproduces pinned boundary loss')
        box%periodic(2)=.true.
        call assemble_state_redist_patch_rereflux_3d(padded,reg,[4,4,4],source,increment,good,box=box)
        call require(good.and.all(increment==reference),'periodic register images retain original kernel output')
        box%periodic=.false.
      endif
    enddo
    previous=increment
    allocate(halo(3,16,16,16),source=-99._dp);kept_halo=halo
    bad_box=box;bad_box%valid_upper(2)=9
    call assemble_state_redist_patch_rereflux_3d(padded,reg,[4,4,4],source,increment,good,halo,bad_box)
    call require(.not.good.and.all(increment==previous).and.all(halo==kept_halo), &
      'physical register invalid domain box rolls back both outputs')
  end subroutine check_physical_domain_rereflux

  subroutine check_rereflux()
    type(eb_geometry_3d) :: g
    real(dp) :: source(1,4,3,3),kappa,wtot,expected,mass,fine_delta(1,2,3,2),saved_fine(1,2,3,2)
    real(dp), allocatable :: increment(:,:,:,:),reference(:,:,:,:),previous(:,:,:,:)
    logical :: boundary(4,3,3),connected(3,3,3,4,3,3),inside(4,3,3),good
    integer :: i,j,k
    call build_axis_plane_eb_geometry_3d(4,3,3,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.37_dp,g,good)
    call require(good,'cut-cell rereflux geometry')
    boundary=.false.;boundary(2,2,2)=.true.;source=0;source(1,2,2,2)=2;connected=.false.;wtot=0
    do k=-1,1
    do j=-1,1
    do i=-1,1
      if(i==0.and.j==0.and.k==0) cycle
      if(g%volume_fraction(2+i,2+j,2+k)<=0) cycle
      connected(i+2,j+2,k+2,2,2,2)=.true.
      wtot=wtot+g%volume_fraction(2+i,2+j,2+k)
    enddo
    enddo
    enddo
    kappa=g%volume_fraction(2,2,2);expected=2*kappa*(1-kappa)/wtot
    call assemble_eb_rereflux_3d(g,source,boundary,connected,increment,good)
    call require(good,'cut-cell rereflux assembly')
    call require(abs(increment(1,2,2,2)-2*kappa)<1.e-14_dp,'cut-cell rereflux self volume fraction')
    call require(abs(increment(1,3,2,2)-expected)<1.e-14_dp,'cut-cell rereflux connected-neighbor density')
    mass=sum(g%volume_fraction*increment(1,:,:,:))
    call require(abs(mass-2*kappa)<2.e-14_dp,'cut-cell rereflux volume-weighted conservation')
    reference=increment
    call assemble_eb_rereflux_3d(g,source,boundary,connected,increment,good,[3,1,1],[4,3,3])
    call require(good.and.all(increment(:,3:4,:,:)==reference(:,3:4,:,:)).and.all(increment(:,1:2,:,:)==0), &
      'rereflux valid-box assembly receives halo roots without clipping denominator')
    boundary(3,2,2)=.true.;source(1,3,2,2)=.5_dp
    call assemble_eb_rereflux_3d(g,source,boundary,connected,increment,good)
    call require(good.and.abs(increment(1,3,2,2)-expected-.5_dp)<1.e-14_dp,'regular rereflux adds density directly')
    previous=increment;connected=.false.
    call assemble_eb_rereflux_3d(g,source,boundary,connected,increment,good)
    call require(.not.good.and.all(increment==previous),'isolated cut-cell rereflux rollback')
    inside=.false.;inside(3,2,2)=.true.;fine_delta=.125_dp
    call add_eb_rereflux_to_fine_3d(reference,inside,[2,3,2],[-4,-3,-3],[-4,-6,-4],fine_delta,good)
    call require(good.and.maxval(abs(fine_delta-.125_dp-expected))<1.e-14_dp, &
      'anisotropic rereflux injects coarse density with negative global coordinates')
    call require(abs(sum(fine_delta-.125_dp)/12-expected)<1.e-14_dp,'fine child mass equals fine-covered coarse transfer')
    saved_fine=fine_delta;inside=.false.
    call add_eb_rereflux_to_fine_3d(reference,inside,[2,3,2],[-4,-3,-3],[-4,-6,-4],fine_delta,good)
    call require(good.and.all(fine_delta==saved_fine),'rereflux fine inside mask')
    call add_eb_rereflux_to_fine_3d(reference,inside,[2,3,2],[-4,-3,-3],[-100,-6,-4],fine_delta,good)
    call require(.not.good.and.all(fine_delta==saved_fine),'rereflux missing coarse donor rollback')
  end subroutine check_rereflux
  subroutine check_redistribution_eos()
    type(eb_geometry_3d) :: g
    type(eb_integrated_flux_3d) :: zero_flux
    type(reactive_eb_3d_config) :: cfg
    real(dp), allocatable :: u(:,:,:,:),t(:,:,:),old(:,:,:,:),ot(:,:,:),correction(:,:,:,:)
    logical :: mask(4,2,2),good
    integer :: nv
    nv=reactive_nvar(size(species));cfg%nx=4;cfg%ny=2;cfg%nz=2
    call build_axis_plane_eb_geometry_3d(4,2,2,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.37_dp,g,good)
    call require(good,'StateRedist reflux EOS geometry')
    allocate(u(nv,4,2,2),t(4,2,2),correction(nv,4,2,2))
    call initialize_reactive_eb_density_sheet_3d(species,cfg,g,u,t,good)
    call require(good,'StateRedist reflux EOS state')
    old=u;ot=t;mask=.false.;correction=0;correction(:,2,1,1)=1.e-7_dp*u(:,2,1,1)
    allocate(zero_flux%x(nv,0:4,2,2),zero_flux%y(nv,4,0:2,2),zero_flux%z(nv,4,2,0:2))
    zero_flux%x=0;zero_flux%y=0;zero_flux%z=0
    call reflux_reactive_eb_state_3d(species,g,mask,zero_flux,zero_flux,u,t,good,redistribution_delta=correction)
    call require(good.and.all(u==old+correction),'StateRedist register reaches accepted reflux state')
    old=u;ot=t;correction=0;correction(:,2,1,1)=-2*u(:,2,1,1)
    call reflux_reactive_eb_state_3d(species,g,mask,zero_flux,zero_flux,u,t,good,redistribution_delta=correction)
    call require(.not.good.and.all(u==old).and.all(t==ot),'StateRedist reflux EOS failure rollback')
  end subroutine check_redistribution_eos
  subroutine check_redistribution_increment()
    type(eb_geometry_3d) :: g
    type(eb_integrated_flux_3d) :: zero_flux
    real(dp) :: correction(1,4,2,2)
    real(dp), allocatable :: result(:,:,:,:),previous(:,:,:,:)
    logical :: mask(4,2,2),good
    logical, allocatable :: updated(:,:,:)
    call build_axis_plane_eb_geometry_3d(4,2,2,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.37_dp,g,good)
    call require(good,'redistribution reflux fixture')
    allocate(zero_flux%x(1,0:4,2,2),zero_flux%y(1,4,0:2,2),zero_flux%z(1,4,2,0:2))
    zero_flux%x=0;zero_flux%y=0;zero_flux%z=0;mask=.false.;correction=0
    correction(1,2,1,1)=.125_dp
    call compute_reactive_eb_reflux_increment_3d(g,mask,zero_flux,zero_flux,result,good,touched=updated, &
      redistribution_delta=correction)
    call require(good.and.all(result==correction).and.updated(2,1,1),'reflux assembles StateRedist correction once')
    previous=result;mask(2,1,1)=.true.
    call compute_reactive_eb_reflux_increment_3d(g,mask,zero_flux,zero_flux,result,good, &
      redistribution_delta=correction)
    call require(.not.good.and.all(result==previous),'unexported fine-covered correction rejection rollback')
    mask=.false.;correction(1,1,1,1)=.25_dp
    call compute_reactive_eb_reflux_increment_3d(g,mask,zero_flux,zero_flux,result,good, &
      redistribution_delta=correction)
    call require(.not.good.and.all(result==previous),'covered-solid correction rejection rollback')
  end subroutine check_redistribution_increment
  subroutine check_exchange_interval()
    real(dp) :: baseline(8,2),gradient(8,2),lo,hi,kept_lo,kept_hi,endpoint(8,2)
    logical :: good
    baseline=0;gradient=0;baseline(1,:)=1
    baseline(6,:)=[-.2_dp,.4_dp];baseline(7,:)=1-baseline(6,:)
    gradient(6,:)=[.5_dp,-.5_dp];gradient(7,:)=-gradient(6,:)
    lo=19;hi=20
    call reactive_eb_exchange_interval_3d(baseline,gradient,lo,hi,good)
    call require(good.and.abs(lo-.4_dp)<1e-15_dp.and.abs(hi-.8_dp)<1e-15_dp, &
      'multiple physical endpoints impose both a funding floor and a withdrawal ceiling')
    endpoint=baseline+lo*gradient
    call require(all(endpoint(1,:)>0).and.all(endpoint(6:,:)>=0),'lower exchange endpoint is represented admissible')
    endpoint=baseline+hi*gradient
    call require(all(endpoint(1,:)>0).and.all(endpoint(6:,:)>=0),'upper exchange endpoint is represented admissible')
    baseline=0;gradient=0;baseline(1,:)=[-.1_dp,1._dp];gradient(1,:)=[.5_dp,-1._dp]
    call reactive_eb_exchange_interval_3d(baseline,gradient,lo,hi,good)
    call require(good.and.lo>.2_dp.and.hi<1,'strict density opens both analytic interval boundaries')
    call require(all(baseline(1,:)+lo*gradient(1,:)>0).and.all(baseline(1,:)+hi*gradient(1,:)>0), &
      'represented density endpoints stay strictly positive')
    kept_lo=lo;kept_hi=hi
    baseline(6,1)=-.9_dp;gradient(6,1)=.5_dp
    call reactive_eb_exchange_interval_3d(baseline,gradient,lo,hi,good)
    call require(.not.good.and.lo==kept_lo.and.hi==kept_hi,'empty exchange interval preserves published bounds')
    baseline(6,1)=ieee_value(0._dp,ieee_quiet_nan)
    call reactive_eb_exchange_interval_3d(baseline,gradient,lo,hi,good)
    call require(.not.good.and.lo==kept_lo.and.hi==kept_hi,'nonfinite exchange predictor rejects transactionally')
    baseline=0;gradient=0;baseline(1,:)=1
    baseline(6,1)=-1;gradient(6,1)=1;baseline(6,2)=1;gradient(6,2)=-1
    call reactive_eb_exchange_interval_3d(baseline,gradient,lo,hi,good)
    call require(good.and.lo==1.and.hi==1,'species permit an exact zero singleton exchange interval')
  end subroutine check_exchange_interval

  subroutine check_interior_donors()
    integer, parameter :: nv=8,nf=3
    real(dp) :: coarse(nv),fine(nv,nf),donor(nv,nf),old(nv),raw(nv,nf),weights(nf),volumes(nf),dv(nf),theta
    real(dp) :: balance(nv),scale(nv),work(nv,nf),original_donor(nv,nf)
    real(dp) :: mixed_coarse(nv),mixed_fine(nv,nf),mixed_raw(nv,nf),mixed_donor(nv,nf)
    real(dp), allocatable :: funded(:,:),remaining(:,:),exchange(:,:),faces(:,:),dc(:),df(:,:)
    real(dp), allocatable :: kept_funded(:,:),kept_remaining(:,:),kept_exchange(:,:)
    logical :: good
    integer :: q
    coarse=0;coarse(1)=1;coarse(6)=.1_dp;coarse(7)=.9_dp
    fine=0;fine(1,:)=1;fine(7,:)=1;fine(2,:)=2;fine(5,:)=30
    donor=fine;donor(6,:)=[.5_dp,.8_dp,.8_dp];donor(7,:)=1-donor(6,:);donor(2,:)=3;donor(5,:)=40
    weights=[.2_dp,.3_dp,.5_dp];volumes=[.25_dp,.5_dp,.25_dp];dv=[.25_dp,.25_dp,.5_dp]
    old=0;raw=0;raw(6,:)=-.2_dp*weights;raw(7,:)=-raw(6,:)
    theta=1
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,fine,volumes,raw,weights,faces,dc,df,theta,good)
    call require(.not.good,'exhausted fine boundary requires interior donors')
    call supply_reactive_eb_interface_donors_3d(coarse,1._dp,old,fine,volumes,raw,weights,donor,dv, &
      funded,remaining,exchange,good)
    call require(good,'interior donors fund the jointly admissible interface')
    call require(all(funded(6:,:)>=0).and.all(remaining(6:,:)>=0),'convex donor exchange preserves all species')
    balance=0;scale=1
    do q=1,nf
      balance=balance+volumes(q)*(funded(:,q)-fine(:,q))+dv(q)*(remaining(:,q)-donor(:,q))
      scale=max(scale,abs(volumes(q)*fine(:,q))+abs(dv(q)*donor(:,q)))
      call require(maxval(abs(volumes(q)*(funded(:,q)-fine(:,q))-exchange(:,q))/scale)<64*epsilon(1._dp), &
        'published donor inventory equals the target state increment')
      call require(maxval(abs(dv(q)*(remaining(:,q)-donor(:,q))+exchange(:,q))/scale)<64*epsilon(1._dp), &
        'published donor inventory equals the opposite interior state increment')
    enddo
    call require(maxval(abs(balance)/scale)<64*epsilon(1._dp),'interior exchange conserves every component')
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,funded,volumes,raw,weights,faces,dc,df,theta,good)
    call require(good.and.all(coarse(6:)+dc(6:)>=0),'funded interface now admits strict coarse species')
    do q=1,nf
      call require(all(funded(6:,q)+df(6:,q)>=0),'funded interface retains strict fine species')
    enddo
    work=fine;original_donor=donor
    do q=1,nf
      donor=fine;donor(:,q)=original_donor(:,q)
      call supply_reactive_eb_interface_donors_3d(coarse,1._dp,old,work,volumes,raw,weights,donor,dv, &
        funded,remaining,exchange,good,target_subface=q)
      call require(good,'selected subface funding does not require other donor inventories')
      call require(all(funded(:,:q-1)==work(:,:q-1)).and.all(funded(:,q+1:)==work(:,q+1:)), &
        'selected subface funding preserves every other target')
      call require(all(exchange(:,:q-1)==0).and.all(exchange(:,q+1:)==0), &
        'selected subface funding publishes only its own extensive exchange')
      work=funded
    enddo
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,work,volumes,raw,weights,faces,dc,df,theta,good)
    call require(good,'sequential selected subfaces admit a joint conservative interface')
    donor=fine;donor(6,1)=.1_dp;donor(7,1)=.9_dp;donor(2,1)=3;donor(5,1)=40
    call supply_reactive_eb_interface_donors_3d(coarse,1._dp,old,fine,volumes,raw,weights,donor,dv, &
      funded,remaining,exchange,good,target_subface=1)
    call require(.not.good,'a weak donor cannot finish the exhausted subface alone')
    call supply_reactive_eb_interface_donors_3d(coarse,1._dp,old,fine,volumes,raw,weights,donor,dv, &
      funded,remaining,exchange,good,target_subface=1,allow_partial=.true.)
    call require(good.and.funded(6,1)>0.and.funded(6,1)<.08_dp, &
      'partial convex funding improves the deficit without declaring interface acceptance')
    balance=volumes(1)*(funded(:,1)-fine(:,1))+dv(1)*(remaining(:,1)-donor(:,1))
    call require(maxval(abs(balance))<64*epsilon(1._dp), 'partial funding conserves all components')
    work=funded;donor(6,1)=.2_dp;donor(7,1)=.8_dp
    call supply_reactive_eb_interface_donors_3d(coarse,1._dp,old,work,volumes,raw,weights,donor,dv, &
      funded,remaining,exchange,good,target_subface=1)
    call require(good.and.funded(6,1)>=.08_dp,'a second donor completes the partially funded subface')
    mixed_coarse=0;mixed_coarse(1)=1;mixed_coarse(6:8)=[.1_dp,.1_dp,.8_dp]
    mixed_fine=0;mixed_fine(1,:)=1;mixed_fine(8,:)=1
    mixed_raw=0;mixed_raw(6,:)=-.2_dp*weights;mixed_raw(7,:)=-.2_dp*weights
    mixed_raw(8,:)=.4_dp*weights
    mixed_donor=mixed_fine;mixed_donor(6,1)=.4_dp;mixed_donor(8,1)=.6_dp
    call supply_reactive_eb_interface_donors_3d(mixed_coarse,1._dp,old,mixed_fine,volumes,mixed_raw,weights, &
      mixed_donor,dv,funded,remaining,exchange,good,target_subface=1,allow_partial=.true.)
    call require(good.and.funded(6,1)>0.and.funded(7,1)==0, &
      'partial donor can improve one chemical deficit while leaving another unchanged')
    balance=volumes(1)*(funded(:,1)-mixed_fine(:,1))+dv(1)*(remaining(:,1)-mixed_donor(:,1))
    call require(maxval(abs(balance))<64*epsilon(1._dp),'mixed deficit partial exchange conserves every species')
    work=funded;mixed_donor=mixed_fine;mixed_donor(7,1)=.4_dp;mixed_donor(8,1)=.6_dp
    call supply_reactive_eb_interface_donors_3d(mixed_coarse,1._dp,old,work,volumes,mixed_raw,weights, &
      mixed_donor,dv,funded,remaining,exchange,good,target_subface=1)
    call require(good.and.funded(6,1)>=.08_dp.and.funded(7,1)>=.08_dp, &
      'a complementary donor completes both chemical deficits without creating a new deficit')
    kept_funded=funded;kept_remaining=remaining;kept_exchange=exchange
    call supply_reactive_eb_interface_donors_3d(coarse,1._dp,old,fine,volumes,raw,weights,donor,dv, &
      funded,remaining,exchange,good,target_subface=0)
    call require(.not.good.and.all(exchange==kept_exchange),'invalid selected subface rejects without mutation')
    donor(6,:)=0;donor(7,:)=1
    call supply_reactive_eb_interface_donors_3d(coarse,1._dp,old,fine,volumes,raw,weights,donor,dv, &
      funded,remaining,exchange,good)
    call require(.not.good.and.all(funded==kept_funded).and.all(remaining==kept_remaining).and. &
      all(exchange==kept_exchange),'unfunded interior donors preserve all prior outputs')
    dv(1)=0
    call supply_reactive_eb_interface_donors_3d(coarse,1._dp,old,fine,volumes,raw,weights,donor,dv, &
      funded,remaining,exchange,good)
    call require(.not.good.and.all(exchange==kept_exchange),'invalid interior volume rejects transactionally')
  end subroutine check_interior_donors

  subroutine check_interface_budget()
    integer, parameter :: nv=8,nf=3
    real(dp) :: coarse(nv),fine(nv,nf),old(nv),raw(nv,nf),weights(nf),volumes(nf),theta,kept_theta
    real(dp) :: actual(nv),expected(nv),endpoint(nv),basis(nv),amount,open_sum
    real(dp), allocatable :: faces(:,:),dc(:),df(:,:),kept_faces(:,:),kept_dc(:),kept_df(:,:)
    logical :: good
    integer :: q,c
    coarse=0;coarse(1)=1;coarse(6)=.1_dp;coarse(7)=.9_dp
    fine=0;fine(1,:)=1;fine(6,:)=.4_dp;fine(7,:)=.6_dp
    weights=[.2_dp,.3_dp,.5_dp];volumes=[.25_dp,.5_dp,.25_dp];old=0
    raw=0;raw(6,:)=[-.03_dp,-.06_dp,-.11_dp];raw(7,:)=-raw(6,:)
    raw(2,:)=2*raw(6,:);raw(3,:)=-3*raw(6,:);raw(4,:)=4*raw(6,:);raw(5,:)=7*raw(6,:)
    theta=19
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,fine,volumes,raw,weights,faces,dc,df,theta,good)
    call require(good.and.theta>0.and.theta<1,'joint coarse/fine interface budget has an admissible blend')
    call require(abs(theta-.5_dp)<1e-14_dp,'coarse tracer binds the common interface correction')
    call require(all(coarse(6:)+dc(6:)>=0),'corrected coarse species remain nonnegative')
    do q=1,nf
      call require(all(fine(6:,q)+df(6:,q)>=0),'corrected fine subface donor species remain nonnegative')
      call require(maxval(abs(faces(:,q)-theta*raw(:,q)))<1e-15_dp, &
        'one correction factor applies to density, species, momentum and energy')
    enddo
    actual=dc
    do q=1,nf
      actual=actual+volumes(q)*df(:,q)
    enddo
    expected=sum(raw,dim=2)-old
    call require(maxval(abs(actual-expected))<1e-15_dp,'joint correction preserves every reflux inventory')
    kept_faces=faces;kept_dc=dc;kept_df=df;kept_theta=theta
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,fine,volumes,raw,weights, &
      faces,dc,df,theta,good,minimum_theta=.6_dp)
    call require(.not.good,'descendant theta floor rejects an empty coupled interval')
    call require(all(faces==kept_faces).and.all(dc==kept_dc).and.all(df==kept_df).and.theta==kept_theta, &
      'descendant theta floor rejection preserves every kernel output')
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,fine,volumes,raw,weights, &
      faces,dc,df,theta,good,minimum_theta=-.1_dp)
    call require(.not.good.and.theta==kept_theta,'invalid descendant theta floor rejects without mutation')
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,fine,volumes,raw,weights, &
      faces,dc,df,theta,good,minimum_theta=.4_dp)
    call require(good.and.theta>=.4_dp,'compatible descendant theta floor preserves an admissible interval')
    coarse(6)=1;coarse(7)=0
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,fine,volumes,raw,weights,faces,dc,df,theta,good)
    call require(good.and.theta==1.and.all(faces==raw).and.all(df==0), &
      'admissible original fine transfers retain exact identity')
    kept_faces=faces;kept_dc=dc;kept_df=df;kept_theta=theta
    ! The earlier stressed MOL pair has no admissible shared-face interval:
    ! coarse fuel 5.15e-5 cannot pay 6.18e-5, while fine fuel is exhausted.
    coarse(6)=5.15e-5_dp;coarse(7)=1-coarse(6)
    fine(6,:)=0;fine(7,:)=1
    raw=0;raw(6,:)=-6.18e-5_dp*weights;raw(7,:)=-raw(6,:)
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,fine,volumes,raw,weights,faces,dc,df,theta,good)
    call require(.not.good,'incompatible coarse/fine bounds require a wider donor network')
    call require(all(faces==kept_faces).and.all(dc==kept_dc).and.all(df==kept_df).and.theta==kept_theta, &
      'infeasible joint correction retains all prior outputs')
    weights(1)=ieee_value(0._dp,ieee_quiet_nan)
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,fine,volumes,raw,weights,faces,dc,df,theta,good)
    call require(.not.good.and.all(faces==kept_faces).and.theta==kept_theta,'invalid interface weights rollback')
    weights=[.2_dp,.3_dp,.5_dp]
    ! A single-point feasible interval must preserve exact zero species rather
    ! than become an artificial failure through an inward safety margin.
    coarse=0;coarse(1)=1;coarse(6)=.125_dp;coarse(7)=.875_dp
    fine=0;fine(1,:)=1;fine(6,:)=.125_dp;fine(7,:)=.875_dp
    volumes=1;weights=1._dp/3;raw=0;raw(6,1)=-.25_dp;raw(7,1)=.25_dp
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,fine,volumes,raw,weights,faces,dc,df,theta,good)
    call require(good.and.theta==.5_dp.and.coarse(6)+dc(6)==0.and.fine(6,1)+df(6,1)==0, &
      'single-point species interval is accepted without clipping')
    ! Density uses a strict positive endpoint even when its analytic bound
    ! permits zero. Other faces/species still share the same correction.
    coarse=0;coarse(1)=.125_dp;coarse(6)=coarse(1)
    fine=0;fine(1,:)=1;fine(6,:)=1;raw=0;raw(1,1)=-.25_dp;raw(6,1)=raw(1,1)
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,fine,volumes,raw,weights,faces,dc,df,theta,good)
    call require(good.and.coarse(1)+dc(1)>0.and.coarse(6)+dc(6)>=0,'joint density bound remains strictly positive')
  end subroutine check_interface_budget

  subroutine check_endpoint_donor()
    real(dp) :: target(7),donor(7),second(7),correction(7),funded(7),remaining(7),packet(7)
    real(dp) :: prior(7),sentinel(7),budget(7)
    logical :: complete,valid
    target=0;target(1)=1;target(5)=1.e6_dp
    target(6)=5.15e-5_dp;target(7)=1-target(6)
    correction=0;correction(6)=-6.18e-5_dp;correction(7)=-correction(6)
    donor=target;donor(6)=6.e-5_dp;donor(7)=1-donor(6)
    second=target;second(6)=2.e-4_dp;second(7)=1-second(6)
    sentinel=-17;funded=sentinel;remaining=sentinel;packet=sentinel
    prior=target+donor
    call supply_reactive_eb_endpoint_donor_3d(target,correction,1._dp,donor,1._dp, &
      funded,remaining,packet,complete,valid)
    call require(valid.and..not.complete,'weak coarse donor gives a partial conservative packet')
    call require(funded(6)>target(6).and.funded(6)+correction(6)<0, &
      'partial packet advances the actual H2 deficit without clipping')
    call require(maxval(abs(funded+remaining-prior))<1.e-15_dp, &
      'partial coarse donor exchange conserves every component')
    prior=funded+second;donor=remaining
    call supply_reactive_eb_endpoint_donor_3d(funded,correction,1._dp,second,1._dp, &
      target,remaining,packet,complete,valid)
    call require(valid.and.complete.and.target(6)+correction(6)>=0, &
      'second coarse donor funds the remaining H2 reflux deficit')
    call require(maxval(abs(target+remaining-prior))<1.e-15_dp, &
      'completed coarse donor exchange conserves every component')
    budget=target;funded=sentinel;remaining=sentinel;packet=sentinel
    donor(6)=0;donor(7)=1
    call supply_reactive_eb_endpoint_donor_3d(budget,correction,1._dp,donor,1._dp, &
      funded,remaining,packet,complete,valid)
    call require(.not.valid.and..not.complete.and.all(funded==sentinel).and. &
      all(remaining==sentinel).and.all(packet==sentinel), &
      'nonfunding donor leaves all proposed outputs untouched')
  end subroutine check_endpoint_donor

  subroutine check_uniform_inventory_packet()
    real(dp) :: target(8),donor(8),debt(8),funded(8),remaining(8),packet(8),initial(8)
    logical :: complete,valid
    target=0;target(1)=1;target(5)=1e6_dp;target(7)=1;target(8)=-200
    donor=target;debt=-1.01_dp*target
    funded=-17;remaining=-17;packet=-17
    call supply_reactive_eb_endpoint_donor_3d(target,debt,.125_dp,donor,2._dp, &
      funded,remaining,packet,complete,valid,passive_components=1)
    call require(.not.valid.and..not.complete.and.all(funded==-17).and. &
      all(remaining==-17).and.all(packet==-17),'uniform inventory borrowing is opt-in and atomic')
    initial=.125_dp*target+2*donor
    call supply_reactive_eb_endpoint_donor_3d(target,debt,.125_dp,donor,2._dp, &
      funded,remaining,packet,complete,valid,passive_components=1,allow_inventory_packet=.true.)
    call require(valid.and.complete.and.funded(1)+debt(1)>density_floor.and.all(funded(6:7)+debt(6:7)>=0), &
      'uniform registered endpoint is funded from positive physical inventory')
    call require(remaining(1)>0.and.all(remaining(6:7)>=0).and.packet(8)<0, &
      'complete physical packet retains donor bounds and signed independent energy')
    call require(maxval(abs(.125_dp*funded+2*remaining-initial)/max(1._dp,abs(initial)))<64*epsilon(1._dp), &
      'unequal-volume registered endpoint funding conserves every component')
    debt=-target
    call supply_reactive_eb_endpoint_donor_3d(target,debt,1._dp,donor,1._dp, &
      funded,remaining,packet,complete,valid,passive_components=1,allow_inventory_packet=.true.)
    call require(valid.and.complete.and.funded(1)+debt(1)>density_floor, &
      'exactly zero post-register density receives a representable positive margin')
    funded=-17;remaining=-17;packet=-17
    debt=-2*target
    call supply_reactive_eb_endpoint_donor_3d(target,debt,1._dp,donor,.001_dp, &
      funded,remaining,packet,complete,valid,passive_components=1,allow_inventory_packet=.true.)
    call require(.not.valid.and..not.complete.and.all(funded==-17).and. &
      all(remaining==-17).and.all(packet==-17),'insufficient complete inventory refuses without partial spending')
    debt=0;debt(6)=-1e-4_dp
    call supply_reactive_eb_endpoint_donor_3d(target,debt,1._dp,donor,1._dp, &
      funded,remaining,packet,complete,valid,passive_components=1,allow_inventory_packet=.true.)
    call require(.not.valid.and..not.complete.and.all(funded==-17).and. &
      all(remaining==-17).and.all(packet==-17),'missing donor species refuses without invented inventory')
  end subroutine check_uniform_inventory_packet

  subroutine check_passive_budgets()
    real(dp) :: target(8),donor(8),correction(8),funded(8),remaining(8),packet(8),kept(8)
    real(dp) :: gas_target(7),gas_donor(7),gas_packet(7),lower,upper,gas_lower,gas_upper,fraction,theta,gas_theta
    real(dp) :: coarse(8),fine(8,2),old(8),raw(8,2),weights(2),volumes(2),total(8)
    real(dp), allocatable :: faces(:,:),dc(:),df(:,:),gas_faces(:,:),gas_dc(:),gas_df(:,:)
    real(dp), allocatable :: supplied(:,:),supplied_donor(:,:),exchange(:,:)
    real(dp), allocatable :: gas_supplied(:,:),gas_supplied_donor(:,:),gas_exchange(:,:)
    logical :: complete,gas_complete,good
    ! A signed independent energy has no species lower bound. It still shares
    ! the exact gas allowance and every component of the conservative transfer.
    target=0;target(1)=1;target(5)=1e6_dp;target(6)=.01_dp;target(7)=.99_dp;target(8)=-100
    correction=0;correction(6)=-.02_dp;correction(7)=.02_dp;correction(8)=200
    lower=-17;upper=-17;gas_lower=-17;gas_upper=-17
    call reactive_eb_exchange_interval_3d(reshape(target(:7),[7,1]),reshape(correction(:7),[7,1]), &
      gas_lower,gas_upper,good)
    call require(good,'gas exchange interval reference')
    call reactive_eb_exchange_interval_3d(reshape(target,[8,1]),reshape(correction,[8,1]), &
      lower,upper,good,passive_components=1)
    call require(good.and.lower==gas_lower.and.upper==gas_upper.and.lower==0, &
      'independent signed energy shares gas interval without a species bound')
    target(6)=0;target(7)=1;donor=target;donor(6)=.2_dp;donor(7)=.8_dp;donor(8)=-200
    correction=0;correction(6)=-.01_dp;correction(7)=.01_dp
    call supply_reactive_eb_endpoint_donor_3d(target(:7),correction(:7),1._dp,donor(:7),1._dp, &
      gas_target,gas_donor,gas_packet,gas_complete,good)
    call require(good.and.gas_complete,'gas donor reference')
    call supply_reactive_eb_endpoint_donor_3d(target,correction,1._dp,donor,1._dp, &
      funded,remaining,packet,complete,good,passive_components=1)
    call require(good.and.complete.and.all(funded(:7)==gas_target).and.all(remaining(:7)==gas_donor).and. &
      all(packet(:7)==gas_packet),'paired scalar donor retains exact gas result')
    fraction=packet(6)/(donor(6)-target(6))
    call require(abs(packet(8)-fraction*(donor(8)-target(8)))<64*epsilon(1._dp)*abs(packet(8)), &
      'held donor uses the gas equalization fraction even at equal densities')
    call require(funded(8)<0.and.remaining(8)<0.and. &
      abs(funded(8)+remaining(8)-target(8)-donor(8))<64*epsilon(1._dp)*300, &
      'signed held energy is transported conservatively without clipping')
    kept=-17;funded=kept;remaining=kept;packet=kept
    call supply_reactive_eb_endpoint_donor_3d(target,correction,1._dp,donor,1._dp, &
      funded,remaining,packet,complete,good)
    call require(.not.good.and.all(funded==kept).and.all(remaining==kept).and.all(packet==kept), &
      'ordinary donor admission remains strict and transactional')
    call supply_reactive_eb_endpoint_donor_3d(target,correction,1._dp,donor,1._dp, &
      funded,remaining,packet,complete,good,passive_components=3)
    call require(.not.good.and.all(packet==kept),'invalid passive boundary rejects before transfer')
    donor(8)=ieee_value(1._dp,ieee_quiet_nan)
    call supply_reactive_eb_endpoint_donor_3d(target,correction,1._dp,donor,1._dp, &
      funded,remaining,packet,complete,good,passive_components=1)
    call require(.not.good.and.all(packet==kept),'nonfinite passive donor rollback')
    coarse=0;coarse(1)=1;coarse(5)=1e6_dp;coarse(6)=.125_dp;coarse(7)=.875_dp;coarse(8)=-50
    fine=0;fine(1,:)=1;fine(5,:)=1e6_dp;fine(6,:)=.5_dp;fine(7,:)=.5_dp;fine(8,:)=[-100._dp,-200._dp]
    old=0;old(8)=11;raw=0;raw(6,1)=-.25_dp;raw(7,1)=.25_dp;raw(8,:)=[123._dp,-31._dp]
    weights=.5_dp;volumes=1;theta=-17;gas_theta=-17
    call reconcile_reactive_eb_interface_budget_3d(coarse(:7),1._dp,old(:7),fine(:7,:),volumes,raw(:7,:), &
      weights,gas_faces,gas_dc,gas_df,gas_theta,good)
    call require(good,'gas face allowance reference')
    call reconcile_reactive_eb_interface_budget_3d(coarse,1._dp,old,fine,volumes,raw,weights, &
      faces,dc,df,theta,good,passive_components=1)
    call require(good.and.theta==gas_theta.and.all(faces(:7,:)==gas_faces).and. &
      all(dc(:7)==gas_dc).and.all(df(:7,:)==gas_df),'independent scalar shares exact gas face allowance')
    call require(maxval(abs(faces(8,:)-(weights*old(8)+theta*(raw(8,:)-weights*old(8)))))< &
      64*epsilon(1._dp)*maxval(abs(raw(8,:))),'independent face uses the same blend')
    total=dc+sum(df,dim=2)-(sum(raw,dim=2)-old)
    call require(maxval(abs(total))<64*epsilon(1._dp)*maxval(abs(raw)), &
      'gas and scalar face reconciliation conserve the imposed transfer')
    coarse=target;fine(:,1)=target;fine(6,1)=.01_dp;fine(7,1)=.99_dp;fine(8,1)=-10
    donor=fine(:,1);donor(6)=.5_dp;donor(7)=.5_dp;donor(8)=-50
    old=0;raw=0;raw(6,1)=-.1_dp;raw(7,1)=.1_dp;weights(1)=1
    call supply_reactive_eb_interface_donors_3d(coarse(:7),1._dp,old(:7),fine(:7,:1),volumes(:1), &
      raw(:7,:1),weights(:1),reshape(donor(:7),[7,1]),volumes(:1), &
      gas_supplied,gas_supplied_donor,gas_exchange,good)
    call require(good,'gas interface donor reference')
    call supply_reactive_eb_interface_donors_3d(coarse,1._dp,old,fine(:,:1),volumes(:1), &
      raw(:,:1),weights(:1),reshape(donor,[8,1]),volumes(:1),supplied,supplied_donor,exchange,good, &
      passive_components=1)
    call require(good.and.all(supplied(:7,:)==gas_supplied).and.all(supplied_donor(:7,:)==gas_supplied_donor).and. &
      all(exchange(:7,:)==gas_exchange),'paired interface funding retains exact gas packets')
    fraction=exchange(6,1)/(donor(6)-fine(6,1))
    call require(abs(exchange(8,1)-fraction*(donor(8)-fine(8,1)))< &
      64*epsilon(1._dp)*abs(exchange(8,1)),'interface funding transports independent scalar with common fraction')
    call require(abs(supplied(8,1)+supplied_donor(8,1)-fine(8,1)-donor(8))<64*epsilon(1._dp)*60, &
      'interface donor held inventory is conserved')
  end subroutine check_passive_budgets

  subroutine check_peer_flux()
    type(eb_geometry_3d) :: domain,owner,receiver
    type(eb_integrated_flux_3d) :: owned,received,kept_flux
    real(dp), allocatable :: u(:,:,:,:),t(:,:,:),kept(:,:,:,:),kept_t(:,:,:)
    real(dp), allocatable :: period(:)
    real(dp) :: vector(nvar),expected(nvar),actual(nvar),area,h(3),rvolume,shift
    integer :: trial,axis,side,ol(3),oh(3),rl(3),rh(3),offset(3),rn(3),first(3),last(3),component
    logical :: valid
    call build_axis_plane_eb_geometry_3d(6,6,6,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.37_dp,domain,valid)
    call require(valid,'peer flux domain geometry')
    allocate(owned%x(nvar,0:3,3,3),owned%y(nvar,3,0:3,3),owned%z(nvar,3,3,0:3))
    do trial=1,2
    if(trial==2) allocate(period(3))
    do axis=1,3
      do side=0,1
        ol=1;oh=3;rl=2;rh=3
        if (side==1) then
          rl(axis)=4;rh(axis)=6
        else
          ol(axis)=4;oh(axis)=6;rl(axis)=1;rh(axis)=3
        end if
        offset=rl-ol;rn=rh-rl+1
        call extract_eb_geometry_3d(domain,ol,oh,owner,valid)
        call require(valid,'owner geometry')
        call extract_eb_geometry_3d(domain,rl,rh,receiver,valid)
        call require(valid,'receiver partial face geometry')
        allocate(u(nvar,rn(1),rn(2),rn(3)),t(rn(1),rn(2),rn(3)))
        allocate(received%x(nvar,0:rn(1),rn(2),rn(3)),received%y(nvar,rn(1),0:rn(2),rn(3)), &
          received%z(nvar,rn(1),rn(2),0:rn(3)))
        config%nx=rn(1);config%ny=rn(2);config%nz=rn(3)
        call initialize_reactive_eb_density_sheet_3d(species,config,receiver,u,t,valid)
        call require(valid,'receiver physical state')
        if(allocated(period)) then
          period=0;period(axis)=1
          shift=-real(2*side-1,dp)
          select case(axis)
          case(1);receiver%x_lower=receiver%x_lower+shift;receiver%x_upper=receiver%x_upper+shift
            where(receiver%boundary_area>0) receiver%boundary_centroid_x=receiver%boundary_centroid_x+shift
          case(2);receiver%y_lower=receiver%y_lower+shift;receiver%y_upper=receiver%y_upper+shift
            where(receiver%boundary_area>0) receiver%boundary_centroid_y=receiver%boundary_centroid_y+shift
          case(3);receiver%z_lower=receiver%z_lower+shift;receiver%z_upper=receiver%z_upper+shift
            where(receiver%boundary_area>0) receiver%boundary_centroid_z=receiver%boundary_centroid_z+shift
          end select
        end if
        vector=1.e-7_dp*u(:,rn(1),rn(2),rn(3))
        do component=1,nvar
          owned%x(component,:,:,:)=vector(component)
          owned%y(component,:,:,:)=vector(component)
          owned%z(component,:,:,:)=vector(component)
        end do
        received%x=0._dp;received%y=0._dp;received%z=0._dp
        kept=u;kept_t=t
        call reconcile_reactive_eb_peer_flux_3d(species,owner,receiver,offset,owned,received,u,t,valid,periodic_length=period)
        call require(valid,'peer flux synchronization')
        first=offset+1;last=offset+rn;first(axis)=side*3;last(axis)=first(axis)
        select case(axis)
        case(1);area=sum(owner%x_face_fraction(first(1):last(1),first(2):last(2),first(3):last(3)))
        case(2);area=sum(owner%y_face_fraction(first(1):last(1),first(2):last(2),first(3):last(3)))
        case(3);area=sum(owner%z_face_fraction(first(1):last(1),first(2):last(2),first(3):last(3)))
        end select
        h=[owner%dx,owner%dy,owner%dz];area=area*product(h)/h(axis)
        expected=real(2*side-1,dp)*area*vector
        rvolume=receiver%dx*receiver%dy*receiver%dz
        do component=1,nvar
          actual(component)=rvolume*sum((u(component,:,:,:)-kept(component,:,:,:))*receiver%volume_fraction)
        end do
        call require(maxval(abs(expected))>1.e-10_dp,'peer correction is nonzero')
        call require(all(abs(actual-expected)<3.e-12_dp*max(1._dp,abs(expected))),'receiver inventory equals face difference')
        call require(all(t>0._dp),'corrected peer temperatures positive')
        kept=u;kept_t=t;kept_flux=received
        ! Matching fluxes are idempotent (up to EOS temperature recovery).
        call reconcile_reactive_eb_peer_flux_3d(species,owner,receiver,offset,owned,received,u,t,valid,periodic_length=period)
        call require(valid.and.all(u==kept),'synchronized conserved state is idempotent')
        kept_t=t
        if(allocated(period)) then
          period(axis)=-1
          call reconcile_reactive_eb_peer_flux_3d(species,owner,receiver,offset,owned,received,u,t,valid, &
            periodic_length=period)
          call require(.not.valid.and.all(u==kept).and.all(t==kept_t),'invalid period rolls back state')
          call require(all(received%x==kept_flux%x).and.all(received%y==kept_flux%y).and. &
            all(received%z==kept_flux%z),'invalid period rolls back flux')
          period(axis)=1
        end if
        do component=1,nvar
          owned%x(component,:,:,:)=-real(2*side-1,dp)*1.e13_dp*vector(component)
          owned%y(component,:,:,:)=-real(2*side-1,dp)*1.e13_dp*vector(component)
          owned%z(component,:,:,:)=-real(2*side-1,dp)*1.e13_dp*vector(component)
        end do
        call reconcile_reactive_eb_peer_flux_3d(species,owner,receiver,offset,owned,received,u,t,valid,periodic_length=period)
        call require(.not.valid.and.all(u==kept).and.all(t==kept_t),'nonphysical correction rolls back state')
        call require(all(received%x==kept_flux%x).and.all(received%y==kept_flux%y).and. &
          all(received%z==kept_flux%z),'failed correction rolls back fluxes')
        call reconcile_reactive_eb_peer_flux_3d(species,owner,receiver,[0,0,0],owned,received,u,t,valid,periodic_length=period)
        call require(.not.valid.and.all(u==kept).and.all(t==kept_t),'overlapping peer flux boxes rejected')
        deallocate(u,t,received%x,received%y,received%z)
      end do
    end do
    end do
  end subroutine check_peer_flux

  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not.condition) error stop message
  end subroutine require
end program test_eb_reactive_reflux_3d
