program test_les_sdc_3d
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
  real(dp) :: ref_boundary(12),ref_source(12),wall_t(6),maximum
  real(dp) :: fx(12,0:4,4,4),fy(12,4,0:4,4),fz(12,4,4,0:4),expected(12,4,4,4)
  integer :: i,j,k,c,mode,b
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
    end do
  end do
  u=initial;t=it;boundary=42;source=43;diffusion=44;fx=45;fy=46;fz=47
  call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,.17_dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
    model='dynamic',test_ratio=-1,transport=transport,flux_x=fx,flux_y=fy,flux_z=fz)
  call require(.not.ok.and.all(u==initial).and.all(t==it),'invalid dynamic filter rolls back state')
  call require(all(boundary==42).and.all(source==43).and.diffusion==44,'failure preserves diagnostics')
  call require(all(fx==45).and.all(fy==46).and.all(fz==47),'failed LES evaluation preserves flux outputs')
  call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,.17_dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
    flux_x=fx(:,1:,:,:),flux_y=fy,flux_z=fz)
  call require(.not.ok.and.all(u==initial).and.all(fx==45),'malformed LES face shape rejects atomically')
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
