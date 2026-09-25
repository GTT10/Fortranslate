program test_spray_mpi_ppm
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: advance_reactive_boundary_hydro_3d
  use spray_mpi_mod, only: advance_mpi_spray_hydro
  implicit none
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: global(12,8,4,4),gt(8,4,4),reference(12,8,4,4),rt(8,4,4),q(12),sound
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),saved(:,:,:,:),st(:,:,:)
  real(dp) :: budget(12),expected(12),mach(6),tin(6),pin(6),coefficient
  real(dp), parameter :: spacing(3)=[.125_dp,.25_dp,.25_dp],dt=1e-5_dp
  character(len=32) :: faces(6),reconstruction,solver
  integer :: rank,nranks,ierr,nx,first,i,j,k,mode,policy,face
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==4,'four short slabs')
  nx=8/nranks;first=rank*nx+1
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  do k=1,4
  do j=1,4
  do i=1,8
    q=0;q(1)=1+.02_dp*sin(real(i+j+k,dp));q(2:4)=[20._dp,10._dp,-5._dp]
    q(5)=500000;q(6)=.05_dp;q(9)=.1_dp;q(12)=.85_dp
    call reactive_primitive_to_conserved(species,q,global(:,i,j,k),gt(i,j,k),sound,ok)
    call require(ok,'initial EOS')
  end do
  end do
  end do
  do policy=1,3
    solver="hllc"
    if(policy>=2) solver="pelec"
    coefficient=0
    if(policy==2) coefficient=.1_dp
    if(policy==3) coefficient=.037_dp
  do mode=1,8
    faces='periodic';mach=-1;tin=0;pin=0
    if(mode>=2) then
      faces='outflow';faces(3:4)='slip_wall'
    end if
    if(mode>=3) then
      face=mode-2
      faces(face)='inflow';mach(face)=.1_dp;tin(face)=700;pin(face)=500000
    end if
    reference=global;rt=gt
    call advance_reactive_boundary_hydro_3d(species,reference,rt,spacing,dt,faces, &
      'characteristic_ppm','mc',solver,expected,ok,inflow_temperature=tin,inflow_mach=mach,inflow_pressure=pin,difmag=coefficient)
    call require(ok,'serial reference')
    u=global(:,first:first+nx-1,:,:);t=gt(first:first+nx-1,:,:)
    saved=u;st=t
    if(mode==1) then
      call advance_mpi_spray_hydro(MPI_COMM_WORLD,species,u,t,spacing,[8,4,4],dt,solver,'mc',ok, &
        reconstruction='characteristic_ppm',difmag=coefficient)
    else
      call advance_mpi_spray_hydro(MPI_COMM_WORLD,species,u,t,spacing,[8,4,4],dt,solver,'mc',ok, &
        faces,budget,inflow_temperature=tin,inflow_mach=mach,inflow_pressure=pin, &
        reconstruction='characteristic_ppm',difmag=coefficient)
      call require(maxval(abs(budget-expected)/max(1._dp,abs(expected)))<2e-13_dp,'surface budget parity')
    end if
    call require(ok,'MPI PPM')
    call require(all(u==reference(:,first:first+nx-1,:,:)),'exact serial state')
    call require(all(t==rt(first:first+nx-1,:,:)),'exact serial temperature')
    u=saved;t=st
    call advance_mpi_spray_hydro(MPI_COMM_WORLD,species,u,t,spacing,[8,4,4],1._dp,solver,'mc',ok, &
      faces,budget,inflow_temperature=tin,inflow_mach=mach,inflow_pressure=pin, &
        reconstruction='characteristic_ppm',difmag=coefficient)
    call require(.not.ok.and.all(u==saved).and.all(t==st).and.all(budget==0),'collective CFL rollback')
    reconstruction='characteristic_ppm'
    if(rank==nranks-1) reconstruction='characteristic_plm'
    call advance_mpi_spray_hydro(MPI_COMM_WORLD,species,u,t,spacing,[8,4,4],dt,solver,'mc',ok, &
      faces,budget,inflow_temperature=tin,inflow_mach=mach,inflow_pressure=pin,reconstruction=reconstruction)
    call require(.not.ok.and.all(u==saved).and.all(t==st),'collective selector mismatch')
  end do
  end do
  coefficient=.037_dp
  if(rank==nranks-1) coefficient=.05_dp
  call advance_mpi_spray_hydro(MPI_COMM_WORLD,species,u,t,spacing,[8,4,4],dt,solver,'mc',ok, &
    reconstruction='characteristic_ppm',difmag=coefficient)
  call require(.not.ok.and.all(u==saved).and.all(t==st),'collective coefficient mismatch rollback')
  coefficient=.037_dp
  if(rank==nranks-1) coefficient=-1._dp
  call advance_mpi_spray_hydro(MPI_COMM_WORLD,species,u,t,spacing,[8,4,4],dt,solver,'mc',ok, &
    reconstruction='characteristic_ppm',difmag=coefficient)
  call require(.not.ok.and.all(u==saved).and.all(t==st),'collective invalid coefficient rollback')
  if(rank==0) print *, 'PASS: short-slab PPM, serial parity, boundary budget and collective rollback'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: status
    if(condition) return
    print *, 'FAIL rank ',rank,': ',label
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
end program
