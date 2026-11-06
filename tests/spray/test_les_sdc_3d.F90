module les_sdc_halo_test_mod
  use precision_mod, only: dp
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d
  implicit none
  logical :: sdc_failure=.false.
  integer :: sdc_calls=0
  integer :: halo_calls=0,source_halo_calls=0
  integer :: auxiliary_calls=0,auxiliary_stages(3)=0
  logical :: auxiliary_failure=.false.
contains
  subroutine sdc_test_auxiliary(state,temp,time,interval,stage,rate,valid)
    real(dp), intent(in) :: state(:,:,:,:),temp(:,:,:),time,interval
    integer, intent(in) :: stage
    real(dp), intent(out) :: rate(:,:,:,:)
    logical, intent(out) :: valid
    auxiliary_calls=auxiliary_calls+1
    if(auxiliary_calls<=3) auxiliary_stages(auxiliary_calls)=stage
    rate=0;rate(2,:,:,:)=merge(100._dp,50._dp,stage==2)
    valid=time>=0.and.interval>0.and..not.(auxiliary_failure.and.stage==2)
  end subroutine sdc_test_auxiliary
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
end module les_sdc_halo_test_mod
program test_les_sdc_3d
  use les_sdc_halo_test_mod, only: sdc_test_halos,sdc_test_source_halos,sdc_failure,halo_calls, &
    sdc_test_auxiliary,auxiliary_calls,auxiliary_stages,auxiliary_failure
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: advance_reactive_sdc_3d
  use les_sdc_3d_mod, only: advance_les_sdc_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp) :: initial(12,4,4,4),it(4,4,4),u(12,4,4,4),t(4,4,4),reference(12,4,4,4),rt(4,4,4)
  real(dp) :: q(12),sound,h(3),dt,boundary(12),source(12),diffusion,before(12),after(12),scale(12)
  real(dp) :: ref_boundary(12),ref_source(12),wall_t(6),maximum,auxiliary(12)
  real(dp) :: fx(12,0:4,4,4),fy(12,4,0:4,4),fz(12,4,4,0:4),expected(12,4,4,4)
  real(dp) :: estimate(12,4,4,4,3)
  integer :: i,j,k,c,mode,b,iteration
  character(len=16) :: faces(6),models(4)
  logical :: ok
  models=[character(len=16)::'smagorinsky','wale','vreman','dynamic']
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok);call require(ok,'transport')
  do k=1,4
  do j=1,4
  do i=1,4
    q=0;q(1)=.3_dp;q(2)=2*j;q(3)=3*i;q(4)=k;q(5)=101325+10*j;q(12)=1
    call reactive_primitive_to_conserved(species,q,initial(:,i,j,k),it(i,j,k),sound,ok)
    call require(ok,'initial state')
  end do
  end do
  end do
  h=.001_dp;dt=1e-9_dp;wall_t=1100;maximum=0
  do c=1,12
    before(c)=sum(initial(c,:,:,:))*product(h)
  end do
  do b=1,2
    faces='periodic'
    if(b==2) faces(1:2)='isothermal_wall'
    reference=initial;rt=it
    call advance_reactive_sdc_3d(species,reference,rt,h,dt,0._dp,2,2,ok=ok, &
      external_force=[10._dp,-20._dp,30._dp],transport=transport,boundary_faces=faces, &
      wall_temperature=wall_t,boundary_change=ref_boundary,source_change=ref_source)
    call require(ok,'zero LES reference')
    u=initial;t=it
    call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,0._dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
      external_force=[10._dp,-20._dp,30._dp],transport=transport,boundary_faces=faces,wall_temperature=wall_t)
    call require(ok.and.diffusion==0,'zero LES accepted')
    call require(maxval(abs(u-reference)/max(1._dp,abs(reference)))<3e-14_dp,'zero LES preserves gas SDC')
    do mode=1,4
      u=initial;t=it;boundary=42;source=43;diffusion=44
      call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,.17_dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
        model=trim(models(mode)),external_force=[10._dp,-20._dp,30._dp],transport=transport, &
        boundary_faces=faces,wall_temperature=wall_t)
      call require(ok.and.diffusion>0,'LES closure participates in SDC')
      call require(any(u/=reference),'nonzero LES changes coupled state')
      do c=1,12
        after(c)=sum(u(c,:,:,:))*product(h)
      end do
      scale=max(abs(before),abs(after),1e-20_dp)
      scale(2:4)=max(scale(2:4),sqrt(before(1)*before(5)))
      maximum=max(maximum,maxval(abs(after-before-boundary-source)/scale))
      call require(all(abs(after-before-boundary-source)<=3e-12_dp*scale),'gas LES boundary force inventory')
      u=initial;t=it
      call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,.17_dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
        model=trim(models(mode)),transport=transport,boundary_faces=faces,wall_temperature=wall_t, &
        flux_x=fx,flux_y=fy,flux_z=fz)
      call require(ok,'combined LES SDC face flux export')
      do k=1,4
      do j=1,4
      do i=1,4
        expected(:,i,j,k)=initial(:,i,j,k)-dt*((fx(:,i,j,k)-fx(:,i-1,j,k))/h(1)+ &
          (fy(:,i,j,k)-fy(:,i,j-1,k))/h(2)+(fz(:,i,j,k)-fz(:,i,j,k-1))/h(3))
      end do
      end do
      end do
      call require(maxval(abs(u-expected)/max(1._dp,abs(expected)))<1e-12_dp, &
        'combined face divergence reconstructs every LES SDC cell')
      if(b==1) then
        expected=u;u=initial;t=it;halo_calls=0
        call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,.17_dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
          model=trim(models(mode)),transport=transport,flux_x=fx,flux_y=fy,flux_z=fz, &
          fill_halos=sdc_test_halos,fill_source_halos=sdc_test_source_halos)
        call require(ok.and.halo_calls>3,'LES supplied stage halos used')
        call require(maxval(abs(u-expected)/max(1._dp,abs(expected)))<1e-12_dp,'supplied halo LES matches periodic path')
      end if
    end do
  end do
  u=initial;t=it;estimate=0
  call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,.17_dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
    transport=transport,new_source_estimate=estimate,fill_halos=sdc_test_halos,fill_source_halos=sdc_test_source_halos)
  call require(ok,'externally seeded reference')
  expected=u;estimate=0
  do iteration=1,2
    u=initial;t=it
    call advance_les_sdc_3d(species,u,t,h,dt,0._dp,1,.17_dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
      transport=transport,new_source_estimate=estimate,fill_halos=sdc_test_halos,fill_source_halos=sdc_test_source_halos)
    call require(ok,'externally coordinated LES iteration')
  end do
  call require(maxval(abs(u-expected)/max(1._dp,abs(expected)))<1e-14_dp,'external and internal LES iterations match')
  u=initial;t=it;estimate=0;sdc_failure=.true.;boundary=42;source=43;diffusion=44;fx=45;fy=46;fz=47
  call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,.17_dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
    transport=transport,new_source_estimate=estimate,flux_x=fx,flux_y=fy,flux_z=fz, &
    fill_halos=sdc_test_halos,fill_source_halos=sdc_test_source_halos)
  call require(.not.ok.and.all(u==initial).and.all(t==it).and.all(estimate==0),'new-time halo failure rolls back')
  call require(all(boundary==42).and.all(source==43).and.diffusion==44.and.all(fx==45),'halo failure diagnostics rollback')
  sdc_failure=.false.
  u=initial;t=it;boundary=42;source=43;diffusion=44;fx=45;fy=46;fz=47
  call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,.17_dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
    model='dynamic',test_ratio=-1,transport=transport,flux_x=fx,flux_y=fy,flux_z=fz)
  call require(.not.ok.and.all(u==initial).and.all(t==it),'invalid dynamic filter rolls back state')
  call require(all(boundary==42).and.all(source==43).and.diffusion==44,'failure preserves diagnostics')
  call require(all(fx==45).and.all(fy==46).and.all(fz==47),'failed LES evaluation preserves flux outputs')
  call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,.17_dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
    flux_x=fx(:,1:,:,:),flux_y=fy,flux_z=fz)
  call require(.not.ok.and.all(u==initial).and.all(fx==45),'malformed LES face shape rejects atomically')
  u=initial;t=it;auxiliary=0;auxiliary_calls=0;auxiliary_stages=-1
  call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,0._dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
    external_force=[10._dp,-20._dp,30._dp],transport=transport,transport_flags=[.false.,.false.,.false.,.false.], &
    auxiliary_source=sdc_test_auxiliary,auxiliary_change=auxiliary)
  call require(ok,'serial LES SDC auxiliary source')
  call require(auxiliary_calls==3.and.all(auxiliary_stages==[0,1,2]),'auxiliary old and iteration stage indices')
  call require(abs(auxiliary(2)-75._dp*dt*64*product(h))<1e-28_dp,'auxiliary half-kick inventory')
  call require(abs(source(2)-10._dp*dt*64*product(h))<1e-28_dp,'force inventory excludes auxiliary source')
  do c=1,12
    after(c)=sum(u(c,:,:,:))*product(h)
  enddo
  scale=max(abs(before),abs(after),1e-30_dp)
  call require(maxval(abs(after-before-source-boundary-auxiliary)/scale)<3e-12_dp, &
    'separate serial gas external and auxiliary inventories')
  reference=u;rt=t;boundary=41;source=42;auxiliary=43;diffusion=44
  auxiliary_failure=.true.;auxiliary_calls=0
  call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,0._dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
    external_force=[10._dp,-20._dp,30._dp],auxiliary_source=sdc_test_auxiliary,auxiliary_change=auxiliary)
  call require(.not.ok.and.auxiliary_calls==3,'late serial auxiliary failure')
  call require(all(u==reference).and.all(t==rt),'auxiliary failure gas temperature rollback')
  call require(all(boundary==41).and.all(source==42).and.all(auxiliary==43).and.diffusion==44, &
    'auxiliary failure ledgers rollback')
  auxiliary_failure=.false.
  print *, 'PASS: LES SDC four closures, transport, boundaries; residual ',maximum
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine require
end program test_les_sdc_3d
