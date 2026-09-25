program test_reactive_pcm_plm_viscosity_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: advance_reactive_boundary_hydro_3d,advance_reactive_full_3d
  implicit none
  integer, parameter :: n=4,nv=12
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: initial(nv,n,n,n),initial_t(n,n,n),u(nv,n,n,n),t(n,n,n)
  real(dp) :: baseline(nv,n,n,n),bt(n,n,n),reference(nv,n,n,n),rt(n,n,n)
  real(dp) :: fx(nv,0:n,n,n),fy(nv,n,0:n,n),fz(nv,n,n,0:n)
  real(dp) :: q(nv),sound,change(nv),other_change(nv),prediction(nv),inventory(nv),theta
  real(dp), parameter :: spacing(3)=[.01_dp,.02_dp,.03_dp],dt=1e-7_dp
  character(len=32) :: faces(6),scheme
  logical :: ok
  integer :: i,j,k,m,mode,method
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  allocate(reactions(0))
  do k=1,n
  do j=1,n
  do i=1,n
    q=0;q(1)=1+.03_dp*cos(real(i+j+k,dp));q(5)=500000
    q(2:4)=[10*sin(real(i,dp)),7*sin(real(j,dp)),4*sin(real(k,dp))]
    q(6)=.05_dp;q(9)=.1_dp;q(nv)=.85_dp
    call reactive_primitive_to_conserved(species,q,initial(:,i,j,k),initial_t(i,j,k),sound,ok)
    call require(ok,'initial state')
  end do
  end do
  end do
  do method=1,2
    scheme='pcm'
    if(method==2) scheme='characteristic_plm'
    do mode=1,2
      faces='periodic'
      if(mode==2) then
        faces(1)='no_slip_wall';faces(2)='outflow'
      end if
      baseline=initial;bt=initial_t
      call advance_reactive_boundary_hydro_3d(species,baseline,bt,spacing,dt,faces,scheme,'mc','pelec',change,ok)
      call require(ok,'omitted coefficient')
      u=initial;t=initial_t
      call advance_reactive_boundary_hydro_3d(species,u,t,spacing,dt,faces,scheme,'mc','pelec',change,ok,difmag=0._dp)
      call require(ok.and.any(u/=baseline),'zero disables default correction')
      u=initial;t=initial_t
      call advance_reactive_boundary_hydro_3d(species,u,t,spacing,dt,faces,scheme,'mc','pelec',change,ok,difmag=.1_dp)
      call require(ok.and.all(u==baseline).and.all(t==bt),'explicit reference equals omitted policy')
      u=initial;t=initial_t
      call advance_reactive_full_3d(species,reactions,transport,u,t,n,n,n,spacing(1),spacing(2),spacing(3),dt, &
        'pelec',.false.,0._dp,0._dp,.false.,.false.,.false.,.false.,.false.,theta,ok, &
        reconstruction=scheme,limiter='mc',boundary_faces=faces,boundary_change=other_change)
      call require(ok.and.all(u==baseline).and.all(t==bt),'omitted full composition reference parity')
      u=initial;t=initial_t
      call advance_reactive_boundary_hydro_3d(species,u,t,spacing,dt,faces,scheme,'mc','pelec',change,ok, &
        flux_x=fx,flux_y=fy,flux_z=fz,difmag=.037_dp)
      call require(ok.and.any(u/=baseline),'positive coefficient changes compressive flow')
      do k=1,n
      do j=1,n
      do i=1,n
        prediction=initial(:,i,j,k)-dt*((fx(:,i,j,k)-fx(:,i-1,j,k))/spacing(1)+ &
          (fy(:,i,j,k)-fy(:,i,j-1,k))/spacing(2)+(fz(:,i,j,k)-fz(:,i,j,k-1))/spacing(3))
        call require(maxval(abs(u(:,i,j,k)-prediction)/max(1._dp,abs(prediction)))<2e-13_dp, &
          'published mean flux reconstructs both-stage conservative update')
      end do
      end do
      end do
      do m=1,nv
        inventory(m)=sum(u(m,:,:,:)-initial(m,:,:,:))*product(spacing)
      end do
      call require(maxval(abs(inventory-change)/max(1._dp,abs(change)))<2e-13_dp,'boundary ledger closure')
      if(mode==1) call require(maxval(abs(change))<2e-13_dp,'periodic cancellation')
      call require(maxval(abs(sum(u(6:,:,:,:),dim=1)-u(1,:,:,:)))<2e-13_dp,'species closure')
      reference=u;rt=t;u=initial;t=initial_t
      call advance_reactive_full_3d(species,reactions,transport,u,t,n,n,n,spacing(1),spacing(2),spacing(3),dt, &
        'pelec',.false.,1e-6_dp,1e-12_dp,.false.,.false.,.false.,.false.,.false.,theta,ok,scheme,'mc', &
        boundary_faces=faces,boundary_change=other_change,difmag=.037_dp)
      call require(ok.and.all(u==reference).and.all(t==rt).and.all(change==other_change),'coupled dispatch parity')
      u=initial;t=initial_t
      call advance_reactive_boundary_hydro_3d(species,u,t,spacing,dt,faces,scheme,'mc','pelec',change,ok,difmag=-1._dp)
      call require(.not.ok.and.all(u==initial).and.all(t==initial_t).and.all(change==0),'negative input rollback')
    end do
  end do
  print *, 'PASS: PCM/PLM viscosity, mean flux, conservation, coupled dispatch and rollback'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(condition) return
    print *, 'FAIL: ',label
    error stop 1
  end subroutine
end program
