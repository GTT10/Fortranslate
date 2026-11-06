module sdc_test_source_mod
  use precision_mod, only: dp
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d
  implicit none
  logical :: sdc_failure=.false.
  integer :: sdc_calls=0
  integer :: halo_calls=0,source_halo_calls=0
contains
  subroutine sdc_test_halos(state,temp,time,width,viscous,ghost,gt,valid, &
      internal_energy_density,internal_energy_ghost)
    real(dp), optional, intent(in) :: internal_energy_density(:,:,:)
    real(dp), optional, intent(out) :: internal_energy_ghost(:,:,:)
    real(dp), intent(in) :: state(:,:,:,:),temp(:,:,:),time
    integer, intent(in) :: width
    logical, intent(in) :: viscous
    real(dp), intent(out) :: ghost(:,:,:,:),gt(:,:,:)
    logical, intent(out) :: valid
    character(len=8) :: faces(6)
    halo_calls=halo_calls+1;valid=.false.;ghost=0;gt=0
    if(sdc_failure.and.time>0) return
    faces='periodic'
    call fill_reactive_boundary_ghosts_3d(state,temp,faces,ghost,gt,valid,ghost_width=width,viscous_wall=viscous, &
      internal_energy_density=internal_energy_density,internal_energy_ghost=internal_energy_ghost)
  end subroutine sdc_test_halos
  subroutine sdc_test_source_halos(source,time,ghost,valid)
    real(dp), intent(in) :: source(:,:,:,:),time
    real(dp), intent(out) :: ghost(:,:,:,:)
    logical, intent(out) :: valid
    integer :: i,j,k,n(3),index(3)
    source_halo_calls=source_halo_calls+1
    n=shape(source(1,:,:,:))
    do k=1,n(3)+8
    do j=1,n(2)+8
    do i=1,n(1)+8
      index=modulo([i,j,k]-5,n)+1
      ghost(:,i,j,k)=source(:,index(1),index(2),index(3))
    end do
    end do
    end do
    valid=time>=0
  end subroutine sdc_test_source_halos
  subroutine sdc_test_sources(state,temp,time,sources,valid)
    real(dp), intent(in) :: state(:,:,:,:),temp(:,:,:),time
    real(dp), intent(out) :: sources(:,:,:,:,:)
    logical, intent(out) :: valid
    sdc_calls=sdc_calls+1;sources=0;valid=.false.
    if(sdc_failure.and.time>0) return
    if(any(temp<=0)) return
    sources(:,:,:,:,1)=12._dp*state;sources(:,:,:,:,2)=-2._dp*state
    valid=.true.
  end subroutine sdc_test_sources
  subroutine sdc_transport_sources(state,temp,time,sources,valid)
    use reactive_3d_mod, only: compute_reactive_external_source_3d
    real(dp), intent(in) :: state(:,:,:,:),temp(:,:,:),time
    real(dp), intent(out) :: sources(:,:,:,:,:)
    logical, intent(out) :: valid
    sdc_calls=sdc_calls+1;sources=0;valid=.false.
    if(sdc_failure.and.time>0) return
    if(size(sources,5)/=2.or.any(temp<=0)) return
    call compute_reactive_external_source_3d(state,[10._dp,-20._dp,30._dp],sources(:,:,:,:,1),valid)
    if(.not.valid) return
    sources(:,:,:,:,2)=.75_dp*sources(:,:,:,:,1)
    sources(:,:,:,:,1)=.25_dp*sources(:,:,:,:,1)
  end subroutine sdc_transport_sources
end module sdc_test_source_mod

program test_reactive_ppm_step_3d
  use sdc_test_source_mod, only: sdc_test_sources,sdc_failure,sdc_calls
  use sdc_test_source_mod, only: sdc_transport_sources
  use sdc_test_source_mod, only: sdc_test_halos,sdc_test_source_halos,halo_calls,source_halo_calls
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use reactive_3d_mod, only: advance_reactive_euler_ppm_ctu_3d,advance_reactive_full_3d,advance_reactive_boundary_hydro_3d
  use reactive_3d_mod, only: compute_reactive_external_source_3d,construct_reactive_sdc_state_3d,advance_reactive_sdc_periodic_3d
  use reactive_3d_mod, only: advance_reactive_sdc_3d,react_sdc_state_3d
  use reactive_transport_3d_mod, only: advance_reactive_transport_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),saved(:,:,:,:),st(:,:,:),reference(:,:,:,:),rt(:,:,:),q(:),before(:),change(:)
  real(dp) :: theta,stage_theta,speed,phase,dt,residual,after
  integer :: nv,i,j,k,m,step,policy
  character(len=5) :: solver
  logical :: ok
  character(len=8) :: faces(6)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  allocate(reactions(0));nv=reactive_nvar(size(species))
  allocate(u(nv,4,4,4),t(4,4,4),q(nv),before(nv),change(nv));dt=1e-5_dp
  do k=1,4
  do j=1,4
  do i=1,4
    phase=2*acos(-1._dp)*real(i+j+k,dp)/4
    q=0;q(1)=1+.05_dp*sin(phase);q(2:4)=[20._dp,10._dp,-5._dp];q(5)=500000
    q(6)=.05_dp+.01_dp*cos(phase);q(9)=.1_dp+.007_dp*sin(phase);q(nv)=1-q(6)-q(9)
    call reactive_primitive_to_conserved(species,q,u(:,i,j,k),t(i,j,k),speed,ok)
    call require(ok,'initial state')
  end do
  end do
  end do
  call check_sdc_iteration()
  call check_sdc_reactions()
  call check_sdc_outflow()
  call check_sdc_flux_export()
  call check_sdc_halos()
  call check_sdc_callback_transport()
  call check_sdc_assembly()
  call check_source_step()
  saved=u;st=t
  do m=1,nv
    before(m)=sum(u(m,:,:,:))
  end do
  do policy=1,2
  solver="hllc"
  if (policy==2) solver="pelec"
  u=saved;t=st
  do step=1,4
    call advance_reactive_euler_ppm_ctu_3d(species,u,t,[.25_dp,.25_dp,.25_dp],dt, &
      .false.,.true.,solver,theta,ok,difmag=merge(0._dp,.1_dp,mod(step,2)==0))
    call require(ok.and.theta>=0.and.theta<=1.and.all(t>0),'accepted CTU step')
    do m=1,nv
      after=sum(u(m,:,:,:));residual=abs(after-before(m))/max(1._dp,abs(before(m)))
      call require(residual<3e-14_dp,'multistep conserved inventory')
    end do
    call require(maxval(abs(sum(u(6:,:,:,:),dim=1)-u(1,:,:,:)))<3e-14_dp,'species mass closure')
  end do
  call require(any(u/=saved),'nontrivial time evolution')
  u=saved;t=st;theta=17
  call advance_reactive_euler_ppm_ctu_3d(species,u,t,[.25_dp,.25_dp,.25_dp],1._dp, &
    .false.,.true.,solver,theta,ok)
  call require(.not.ok.and.all(u==saved).and.all(t==st).and.theta==17,'excessive dt complete rollback')
  call advance_reactive_euler_ppm_ctu_3d(species,u,t,[.25_dp,.25_dp,.25_dp],0._dp, &
    .false.,.true.,solver,theta,ok)
  call require(ok.and.all(u==saved).and.all(t==st).and.theta==1,'zero step exact identity')
  reference=saved;rt=st
  call advance_reactive_transport_3d(species,transport,reference,rt,4,4,4,.25_dp,.25_dp,.25_dp,.5_dp*dt, &
    .true.,.true.,.true.,.true.,stage_theta,ok)
  call require(ok,'first transport half')
  call advance_reactive_euler_ppm_ctu_3d(species,reference,rt,[.25_dp,.25_dp,.25_dp],dt, &
    .false.,.true.,solver,theta,ok)
  call require(ok,'reference PPM hydro')
  call advance_reactive_transport_3d(species,transport,reference,rt,4,4,4,.25_dp,.25_dp,.25_dp,.5_dp*dt, &
    .true.,.true.,.true.,.true.,stage_theta,ok)
  call require(ok,'second transport half')
  call advance_reactive_full_3d(species,reactions,transport,u,t,4,4,4,.25_dp,.25_dp,.25_dp,dt,solver, &
    .false.,1e-7_dp,1e-12_dp,.true.,.true.,.true.,.true.,.true.,theta,ok,reconstruction='characteristic_ppm')
  call require(ok.and.all(u==reference).and.all(t==rt),'exact transport PPM split dispatch')
  u=saved;t=st;faces='outflow'
  reference=saved;rt=st
  call advance_reactive_transport_3d(species,transport,reference,rt,4,4,4,.25_dp,.25_dp,.25_dp,.5_dp*dt, &
    .true.,.true.,.true.,.true.,stage_theta,ok,faces)
  call require(ok,'physical first transport half')
  call advance_reactive_boundary_hydro_3d(species,reference,rt,[.25_dp,.25_dp,.25_dp],dt, &
    faces,'characteristic_ppm','mc',solver,change,ok)
  call require(ok,'physical PPM hydro')
  call advance_reactive_transport_3d(species,transport,reference,rt,4,4,4,.25_dp,.25_dp,.25_dp,.5_dp*dt, &
    .true.,.true.,.true.,.true.,stage_theta,ok,faces)
  call require(ok,'physical second transport half')
  call advance_reactive_full_3d(species,reactions,transport,u,t,4,4,4,.25_dp,.25_dp,.25_dp,dt,solver, &
    .false.,1e-7_dp,1e-12_dp,.true.,.true.,.true.,.true.,.true.,theta,ok, &
    reconstruction='characteristic_ppm',boundary_faces=faces)
  call require(ok.and.all(u==reference).and.all(t==rt),'exact physical transport PPM split dispatch')
  end do
  print *, 'test_reactive_ppm_step_3d: PASS'
contains
  subroutine check_sdc_callback_transport()
    real(dp) :: state(nv,4,4,4),temp(4,4,4),expected(nv,4,4,4),et(4,4,4)
    real(dp) :: change(nv),reference_change(nv),boundary(nv),reference_boundary(nv)
    real(dp) :: fx(nv,0:4,4,4),fy(nv,4,0:4,4),fz(nv,4,4,0:4)
    real(dp) :: rx(nv,0:4,4,4),ry(nv,4,0:4,4),rz(nv,4,4,0:4)
    real(dp) :: estimate(nv,4,4,4,3)
    character(len=8) :: bc(6)
    integer :: mode
    do mode=1,2
      bc='periodic'
      if(mode==2) bc(1:2)='outflow'
      expected=u;et=t
      call advance_reactive_sdc_3d(species,expected,et,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
        2,2,ok=ok,external_force=[10._dp,-20._dp,30._dp],transport=transport,boundary_faces=bc, &
        source_change=reference_change,boundary_change=reference_boundary,flux_x=rx,flux_y=ry,flux_z=rz)
      call require(ok,'builtin force plus transport reference')
      state=u;temp=t;sdc_calls=0
      call advance_reactive_sdc_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
        2,3,sdc_transport_sources,ok,transport=transport,boundary_faces=bc, &
        source_change=change,boundary_change=boundary,flux_x=fx,flux_y=fy,flux_z=fz)
      call require(ok.and.sdc_calls==3,'callback receives two families alongside final transport family')
      call require(maxval(abs(state-expected)/max(1._dp,abs(expected)))<3e-14_dp.and. &
        maxval(abs(temp-et)/et)<3e-14_dp,'callback transport state matches builtin force')
      call require(all(abs(change-reference_change)<=1e-12_dp*max(1._dp,abs(reference_change))).and. &
        all(abs(boundary-reference_boundary)<=1e-12_dp*max(1._dp,abs(reference_boundary))), &
        'callback transport source and boundary inventories')
      call require(all(abs(fx-rx)<=1e-12_dp*max(1._dp,abs(rx))).and. &
        all(abs(fy-ry)<=1e-12_dp*max(1._dp,abs(ry))).and. &
        all(abs(fz-rz)<=1e-12_dp*max(1._dp,abs(rz))),'callback retains full hydro transport flux export')
    end do
    state=u;temp=t;estimate=0;sdc_failure=.true.
    call advance_reactive_sdc_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
      2,3,sdc_transport_sources,ok,transport=transport,new_source_estimate=estimate)
    call require(.not.ok.and.all(state==u).and.all(temp==t).and.all(estimate==0),'callback transport late rollback')
    sdc_failure=.false.
  end subroutine check_sdc_callback_transport

  subroutine check_sdc_halos()
    real(dp) :: state(nv,4,4,4),temp(4,4,4),expected(nv,4,4,4),et(4,4,4)
    expected=u;et=t
    call advance_reactive_sdc_3d(species,expected,et,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
      2,2,ok=ok,external_force=[0._dp,0._dp,0._dp],transport=transport)
    call require(ok,'SDC native periodic halo baseline')
    state=u;temp=t;halo_calls=0;source_halo_calls=0
    call advance_reactive_sdc_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
      2,2,ok=ok,external_force=[0._dp,0._dp,0._dp],transport=transport, &
      fill_halos=sdc_test_halos,fill_source_halos=sdc_test_source_halos)
    call require(ok.and.halo_calls==4.and.source_halo_calls==2,'SDC staged state and source halo callbacks')
    call require(maxval(abs(state-expected)/max(1._dp,abs(expected)))<1e-12_dp.and. &
      maxval(abs(temp-et))<1e-7_dp,'SDC caller halos reproduce native periodic update')
    state=u;temp=t;sdc_failure=.true.
    call advance_reactive_sdc_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
      2,2,ok=ok,external_force=[0._dp,0._dp,0._dp],transport=transport, &
      fill_halos=sdc_test_halos,fill_source_halos=sdc_test_source_halos)
    call require(.not.ok.and.all(state==u).and.all(temp==t),'SDC new-time halo failure rolls back')
    sdc_failure=.false.
  end subroutine check_sdc_halos

  subroutine check_sdc_flux_export()
    real(dp) :: state(nv,4,4,4),temp(4,4,4),expected(nv,4,4,4)
    real(dp) :: fx(nv,0:4,4,4),fy(nv,4,0:4,4),fz(nv,4,4,0:4)
    character(len=8) :: bc(6)
    integer :: mode,ii,jj,kk
    do mode=1,2
      bc='periodic'
      if(mode==2) bc(1:2)='outflow'
      state=u;temp=t
      call advance_reactive_sdc_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
        2,2,ok=ok,external_force=[0._dp,0._dp,0._dp],boundary_faces=bc,transport=transport, &
        flux_x=fx,flux_y=fy,flux_z=fz)
      call require(ok,'SDC exported hydro plus transport fluxes')
      do kk=1,4
      do jj=1,4
      do ii=1,4
        expected(:,ii,jj,kk)=u(:,ii,jj,kk)-dt/.25_dp*(fx(:,ii,jj,kk)-fx(:,ii-1,jj,kk)+ &
          fy(:,ii,jj,kk)-fy(:,ii,jj-1,kk)+fz(:,ii,jj,kk)-fz(:,ii,jj,kk-1))
      end do
      end do
      end do
      call require(maxval(abs(state-expected)/max(1._dp,abs(expected)))<1e-12_dp, &
        'SDC full-face divergence reconstructs every cell update')
    end do
    state=u;temp=t;fx=1;fy=1;fz=1
    call advance_reactive_sdc_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
      2,1,ok=ok,external_force=[0._dp,0._dp,0._dp],flux_x=fx(:,1:,:,:),flux_y=fy,flux_z=fz)
    call require(.not.ok.and.all(state==u).and.all(temp==t).and.all(fx(:,1:,:,:)==0).and. &
      all(fy==0).and.all(fz==0),'SDC malformed flux output rejects without state mutation')
  end subroutine check_sdc_flux_export

  subroutine check_sdc_outflow()
    real(dp) :: state(nv,4,4,4),temp(4,4,4),expected(nv,4,4,4),et(4,4,4)
    real(dp) :: boundary(nv),baseline_boundary(nv),source(nv),difference
    character(len=8) :: bc(6)
    integer :: component
    bc='periodic';bc(1:2)='outflow'
    state=u;temp=t;expected=u;et=t
    call advance_reactive_boundary_hydro_3d(species,expected,et,[.25_dp,.25_dp,.25_dp],dt, &
      bc,'characteristic_ppm','mc','pelec',baseline_boundary,ok)
    call require(ok,'SDC outflow reference hydro')
    call advance_reactive_sdc_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
      2,1,ok=ok,external_force=[0._dp,0._dp,0._dp],boundary_faces=bc,boundary_change=boundary)
    call require(ok,'SDC outflow step')
    call require(maxval(abs(state-expected)/max(1._dp,abs(expected)))<1e-13_dp, &
      'SDC zero-source outflow matches boundary CTU')
    call require(maxval(abs(boundary-baseline_boundary))<1e-10_dp,'SDC outflow flux inventory')
    state=u;temp=t
    call advance_reactive_sdc_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
      2,1,ok=ok,external_force=[20000._dp,-30000._dp,40000._dp],boundary_faces=bc, &
      boundary_change=boundary,source_change=source)
    call require(ok,'SDC forced outflow')
    do component=1,nv
      difference=sum(state(component,:,:,:)-u(component,:,:,:))*.25_dp**3-boundary(component)-source(component)
      call require(abs(difference)/max(1._dp,abs(sum(u(component,:,:,:))*.25_dp**3))<1e-12_dp, &
        'SDC outflow boundary and volume source inventories')
    end do
  end subroutine check_sdc_outflow

  subroutine check_sdc_reactions()
    type(elementary_reaction), allocatable :: active_reactions(:)
    real(dp) :: state(nv,4,4,4),temp(4,4,4),old(nv,4,4,4),ot(4,4,4),ir(nv-4,4,4,4)
    real(dp) :: cell(nv-5),zero_rates(nv-5),cell_t,rhoe,h,prior_ir(nv-4,4,4,4)
    real(dp) :: primitive(nv),initial(nv),initial_t,sound
    real(dp) :: stage_state(nv,1,1,1),stage_t(1,1,1),stage_ir(nv-4,1,1,1),stage_rates(nv,1,1,1)
    real(dp) :: heat(4,4,4),stage_heat(1,1,1),prior_heat(4,4,4)
    integer :: count
    call load_h2o2_elementary_mechanism(active_reactions,ok)
    call require(ok,'SDC reaction mechanism')
    primitive=0;primitive(1)=1;primitive(5)=500000
    primitive(6)=.05_dp;primitive(7:8)=1e-5_dp;primitive(9)=.1_dp;primitive(10)=1e-5_dp
    primitive(nv)=1-sum(primitive(6:nv-1))
    call reactive_primitive_to_conserved(species,primitive,initial,initial_t,sound,ok)
    call require(ok,'SDC reactive mixture with radical seeds')
    do k=1,4
    do j=1,4
    do i=1,4
      old(:,i,j,k)=initial;ot(i,j,k)=initial_t
    end do
    end do
    end do
    cell=old(6:,1,1,1);cell_t=ot(1,1,1);h=1e-8_dp;zero_rates=0
    rhoe=old(5,1,1,1)-.5_dp*sum(old(2:4,1,1,1)**2)/old(1,1,1,1)
    call advance_forced_reactor_interval(species,active_reactions,rhoe,zero_rates,0._dp,h,1e-6_dp, &
      1e-12_dp*old(1,1,1,1),cell,cell_t,count,ok)
    call require(ok,'SDC independent cell integration')
    stage_state=old(:,1:1,1:1,1:1);stage_t=ot(1:1,1:1,1:1);stage_ir=0;stage_rates=0
    call react_sdc_state_3d(species,active_reactions,old(:,1:1,1:1,1:1),ot(1:1,1:1,1:1), &
      stage_rates,h,1e-6_dp,1e-12_dp,stage_state,stage_t,stage_ir,ok,heat_release=stage_heat)
    call require(ok,'SDC standalone reaction stage')
    call require(maxval(abs(stage_state(6:,1,1,1)-cell))<1e-12_dp.and.abs(stage_t(1,1,1)-cell_t)<1e-6_dp, &
      'SDC standalone reaction matches full-interval cell integration')
    stage_state=7;stage_t=8;stage_ir=9
    call react_sdc_state_3d(species,active_reactions,old(:,1:1,1:1,1:1),ot(1:1,1:1,1:1), &
      stage_rates,-h,1e-6_dp,1e-12_dp,stage_state,stage_t,stage_ir,ok)
    call require(.not.ok.and.all(stage_state==7).and.all(stage_t==8).and.all(stage_ir==9), &
      'SDC reaction stage failure preserves all outputs')
    state=old;temp=ot;ir=0
    call advance_reactive_sdc_periodic_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],h,0._dp, &
      2,1,ok=ok,external_force=[0._dp,0._dp,0._dp],reactions=active_reactions,reaction_source=ir,heat_release=heat)
    call require(ok,'SDC reacting step')
    call require(all(abs(heat-stage_heat(1,1,1))<=32*epsilon(1._dp)*max(1._dp,abs(stage_heat(1,1,1)))), &
      'full SDC reports accepted reaction-stage heat release')
    call require(maxval(abs(state(6:,1,1,1)-cell))<1e-12_dp.and.abs(temp(1,1,1)-cell_t)<1e-6_dp, &
      'SDC uniform chemistry equals single interval without repeated reaction advancement')
    call require(maxval(abs(cell-old(6:,1,1,1)))>1e-14_dp,'SDC reaction is active')
    call require(maxval(abs(ir(:nv-5,1,1,1)-(cell-old(6:,1,1,1))/h))<1e-7_dp,'SDC reaction source inventory')
    state=old;temp=ot;prior_ir=ir;prior_heat=heat;sdc_failure=.true.
    call advance_reactive_sdc_periodic_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],h,0._dp, &
      2,2,sdc_test_sources,ok,reactions=active_reactions,reaction_source=ir,heat_release=heat)
    call require(.not.ok.and.all(state==old).and.all(temp==ot).and.all(ir==prior_ir).and.all(heat==prior_heat), &
      'SDC failure preserves gas, source history and heat release')
    sdc_failure=.false.
  end subroutine check_sdc_reactions

  subroutine check_sdc_iteration()
    real(dp) :: state(nv,4,4,4),temp(4,4,4),old(nv,4,4,4),ot(4,4,4),factor,expected_factor
    real(dp) :: estimate(nv,4,4,4,2),prior_estimate(nv,4,4,4,2),expected(nv,4,4,4)
    integer :: iterations,iteration
    do k=1,4
    do j=1,4
    do i=1,4
      old(:,i,j,k)=u(:,1,1,1);ot(i,j,k)=t(1,1,1)
    end do
    end do
    end do
    do iterations=1,3,2
      state=old;temp=ot;sdc_calls=0
      call advance_reactive_sdc_periodic_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
        iterations,2,sdc_test_sources,ok)
      factor=1
      do iteration=1,iterations
        factor=1+5*dt*(1+factor)
      end do
      expected_factor=1+5*dt*(1+factor)
      call require(ok.and.sdc_calls==iterations+1,'SDC old source once and new source each iteration')
      call require(maxval(abs(state-expected_factor*old)/max(1._dp,abs(old)))<3e-14_dp, &
        'SDC rebuilds from old state instead of accumulating steps')
    end do
    expected=state
    call sdc_test_sources(old,ot,0._dp,estimate,ok)
    call require(ok,'SDC initialize external iteration estimate')
    do iteration=1,3
      state=old;temp=ot
      call advance_reactive_sdc_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
        1,2,sdc_test_sources,ok,new_source_estimate=estimate)
      call require(ok,'SDC externally coordinated iteration')
    end do
    call require(maxval(abs(state-expected)/max(1._dp,abs(expected)))<1e-14_dp, &
      'SDC three coordinated iterations match three internal iterations')
    prior_estimate=estimate
    state=old;temp=ot;sdc_failure=.true.
    call advance_reactive_sdc_periodic_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,0._dp, &
      2,2,sdc_test_sources,ok,new_source_estimate=estimate)
    call require(.not.ok.and.all(state==old).and.all(temp==ot),'SDC new source failure rolls back entire step')
    call require(all(estimate==prior_estimate),'SDC failure preserves external source estimate')
    sdc_failure=.false.
  end subroutine check_sdc_iteration
  subroutine check_sdc_assembly()
    real(dp) :: old(nv,2,2,2),state(nv,2,2,2),os(nv,2,2,2,2),ns(nv,2,2,2,2)
    real(dp) :: hydro(nv,2,2,2),reaction(nv-4,2,2,2),expected(nv,2,2,2)
    integer :: component
    do component=1,nv
      old(component,:,:,:)=component
    end do
    os(:,:,:,:,1)=2;ns(:,:,:,:,1)=6
    os(:,:,:,:,2)=-4;ns(:,:,:,:,2)=8
    hydro=3;reaction=0;reaction(1,:,:,:)=2;reaction(nv-5,:,:,:)=-2;reaction(nv-4,:,:,:)=5
    state=-1;expected=old+1.125_dp
    expected(6,:,:,:)=expected(6,:,:,:)+.25_dp
    expected(nv,:,:,:)=expected(nv,:,:,:)-.25_dp
    expected(5,:,:,:)=expected(5,:,:,:)+.625_dp
    call construct_reactive_sdc_state_3d(old,.125_dp,os,ns,state,ok,hydro,reaction)
    call require(ok.and.all(state==expected),'SDC multiple sources hydro and reaction assembly')
    call construct_reactive_sdc_state_3d(old,.125_dp,os(:,:,:,:,:0),ns(:,:,:,:,:0),state,ok)
    call require(ok.and.all(state==old),'SDC empty source list identity')
    call construct_reactive_sdc_state_3d(old,.125_dp,os,ns,state,ok,reaction_source=reaction)
    call require(ok.and.all(state(1:4,:,:,:)==old(1:4,:,:,:)+.75_dp),'SDC reactions exclude density and momentum')
    expected=state
    call construct_reactive_sdc_state_3d(old,-.125_dp,os,ns,state,ok)
    call require(.not.ok.and.all(state==expected),'SDC invalid interval preserves candidate')
    call construct_reactive_sdc_state_3d(old,.125_dp,os,ns(:,:,:,:,:1),state,ok)
    call require(.not.ok.and.all(state==expected),'SDC inconsistent source lists preserve candidate')
  end subroutine check_sdc_assembly
  subroutine check_source_step()
    real(dp) :: state(nv,4,4,4),temp(4,4,4),source(nv,4,4,4),old(nv,4,4,4),ot(4,4,4),lim
    real(dp) :: halo_source(nv,-3:8,-3:8,-3:8),boundary(nv)
    character(len=8) :: walls(6)
    real(dp) :: composed(nv,4,4,4),composed_t(4,4,4),boundary_ref(nv),transport_limit
    real(dp) :: force(3),work,expected_energy,source_total(nv)
    integer :: component,bc
    state=u;temp=t
    call compute_reactive_external_source_3d(state,[2._dp,-3._dp,4._dp],source,ok)
    call require(ok,'external force source generation')
    call require(all(source(1,:,:,:)==0).and.all(source(6:,:,:,:)==0),'force leaves mass and species unchanged')
    call require(all(source(2,:,:,:)==2).and.all(source(3,:,:,:)==-3).and.all(source(4,:,:,:)==4), &
      'force is per volume independent of density')
    call require(maxval(abs(source(5,:,:,:)+10._dp))<1e-13_dp,'external source mechanical work')
    state(1,4,4,4)=-1
    call compute_reactive_external_source_3d(state,[2._dp,-3._dp,4._dp],source,ok)
    call require(.not.ok.and.all(source==0),'external source rejects invalid density transactionally')
    ! Nonuniform fixture: periodic flux divergence must leave only source inventory.
    state=u;temp=t;source=0
    source(2,:,:,:)=.2_dp;source(5,:,:,:)=100._dp
    source(6,:,:,:)=.01_dp;source(nv,:,:,:)=-.01_dp
    call advance_reactive_euler_ppm_ctu_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt, &
      .false.,.true.,'pelec',lim,ok,conserved_source=source)
    call require(ok,'source CTU step')
    do component=1,nv
      call require(abs(sum(state(component,:,:,:))-sum(u(component,:,:,:))- &
        dt*sum(source(component,:,:,:)))<3e-14_dp*max(1._dp,abs(sum(u(component,:,:,:)))), &
        'periodic source inventory')
    end do
    call require(maxval(abs(sum(state(6:,:,:,:),dim=1)-state(1,:,:,:)))<3e-14_dp,'source species closure')
    state=u;temp=t;walls='outflow'
    do k=-3,8
    do j=-3,8
    do i=-3,8
      halo_source(:,i,j,k)=source(:,1,1,1)
    end do
    end do
    end do
    call advance_reactive_boundary_hydro_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,walls, &
      'characteristic_ppm','mc','pelec',boundary,ok,conserved_source=halo_source)
    call require(ok,'boundary source CTU step')
    do component=1,nv
      call require(abs((sum(state(component,:,:,:))-sum(u(component,:,:,:)))*.25_dp**3- &
        boundary(component)-dt*sum(source(component,:,:,:))*.25_dp**3)< &
        3e-14_dp*max(1._dp,abs(sum(u(component,:,:,:)))*.25_dp**3),'boundary and source inventory separation')
    end do
    do bc=1,2
      walls='periodic'
      if(bc==2) walls='outflow'
      state=u;temp=t
      call advance_reactive_transport_3d(species,transport,state,temp,4,4,4,.25_dp,.25_dp,.25_dp,.5_dp*dt, &
        .true.,.true.,.true.,.true.,transport_limit,ok,walls)
      call require(ok,'source split first transport')
      if(bc==1) then
        call advance_reactive_euler_ppm_ctu_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt, &
          .false.,.true.,'pelec',lim,ok,conserved_source=source)
      else
        call advance_reactive_boundary_hydro_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,walls, &
          'characteristic_ppm','mc','pelec',boundary_ref,ok,conserved_source=halo_source)
      end if
      call require(ok,'source split reference hydro')
      call advance_reactive_transport_3d(species,transport,state,temp,4,4,4,.25_dp,.25_dp,.25_dp,.5_dp*dt, &
        .true.,.true.,.true.,.true.,transport_limit,ok,walls)
      call require(ok,'source split second transport')
      composed=u;composed_t=t
      call advance_reactive_full_3d(species,reactions,transport,composed,composed_t,4,4,4, &
        .25_dp,.25_dp,.25_dp,dt,'pelec',.false.,1e-7_dp,1e-12_dp,.true.,.true.,.true.,.true.,.true.,lim,ok, &
        reconstruction='characteristic_ppm',boundary_faces=walls,conserved_source=halo_source)
      call require(ok.and.all(composed==state).and.all(composed_t==temp),'exact source transport split')
    end do
    state=u;temp=t;halo_source=0;halo_source(5,:,:,:)=-1e14_dp
    call advance_reactive_boundary_hydro_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,walls, &
      'characteristic_ppm','mc','pelec',boundary,ok,conserved_source=halo_source)
    call require(.not.ok.and.all(state==u).and.all(temp==t).and.all(boundary==0),'boundary source rollback')
    call advance_reactive_full_3d(species,reactions,transport,state,temp,4,4,4, &
      .25_dp,.25_dp,.25_dp,dt,'pelec',.false.,1e-7_dp,1e-12_dp,.true.,.true.,.true.,.true.,.true.,lim,ok, &
      reconstruction='characteristic_ppm',boundary_faces=walls,conserved_source=halo_source,boundary_change=boundary)
    call require(.not.ok.and.all(state==u).and.all(temp==t).and.all(boundary==0),'source split complete rollback')
    ! Uniform fields isolate the source increment, including density injection.
    do k=1,4
    do j=1,4
    do i=1,4
      old(:,i,j,k)=u(:,1,1,1);ot(i,j,k)=t(1,1,1)
      source(:,i,j,k)=.1_dp*old(:,i,j,k)
    end do
    end do
    end do
    state=old;temp=ot
    call advance_reactive_euler_ppm_ctu_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt, &
      .false.,.true.,'pelec',lim,ok,conserved_source=source)
    call require(ok.and.maxval(abs(state-old-dt*source)/max(1._dp,abs(old)))<3e-14_dp, &
      'uniform source applied exactly once')
    force=[2e4_dp,-3e4_dp,4e4_dp]
    state=old;temp=ot
    call advance_reactive_full_3d(species,reactions,transport,state,temp,4,4,4, &
      .25_dp,.25_dp,.25_dp,dt,'pelec',.false.,1e-7_dp,1e-12_dp,.false.,.true.,.true.,.true.,.true.,lim,ok, &
      reconstruction='characteristic_ppm',external_force=force,source_change=source_total)
    call require(ok,'time-centered external force update')
    work=sum(old(2:4,1,1,1)/old(1,1,1,1)*force)
    expected_energy=dt*work+.5_dp*dt**2*sum(force**2)/old(1,1,1,1)
    call require(maxval(abs(state(2:4,1,1,1)-old(2:4,1,1,1)-dt*force))<1e-13_dp, &
      'exact external force momentum increment')
    call require(abs(state(5,1,1,1)-old(5,1,1,1)-expected_energy)< &
      3e-14_dp*abs(old(5,1,1,1)),'external work includes quadratic kinetic increment')
    call require(maxval(abs(temp-ot))<3e-10_dp,'external force preserves uniform internal temperature')
    do component=1,nv
      call require(abs(sum(state(component,:,:,:)-old(component,:,:,:))*.25_dp**3-source_total(component))< &
        3e-14_dp*max(1._dp,abs(old(component,1,1,1))),'external force source ledger')
    end do
    state=u;temp=t;walls='outflow'
    call advance_reactive_full_3d(species,reactions,transport,state,temp,4,4,4, &
      .25_dp,.25_dp,.25_dp,dt,'pelec',.false.,1e-7_dp,1e-12_dp,.false.,.true.,.true.,.true.,.true.,lim,ok, &
      reconstruction='characteristic_ppm',boundary_faces=walls,external_force=force,source_change=source_total, &
      boundary_change=boundary)
    call require(ok,'outflow external force update')
    do component=1,nv
      call require(abs(sum(state(component,:,:,:)-u(component,:,:,:))*.25_dp**3- &
        source_total(component)-boundary(component))<3e-14_dp*max(1._dp,abs(sum(u(component,:,:,:)))*.25_dp**3), &
        'outflow external force inventory')
    end do
    state=old;temp=ot;lim=17
    source=0;source(5,:,:,:)=-1e14_dp
    call advance_reactive_euler_ppm_ctu_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt, &
      .false.,.true.,'pelec',lim,ok,conserved_source=source)
    call require(.not.ok.and.all(state==old).and.all(temp==ot).and.lim==17,'source failure complete rollback')
  end subroutine check_source_step
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
