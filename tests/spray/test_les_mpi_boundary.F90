program test_les_mpi_boundary
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use spray_partition_mod, only: spray_cartesian_block
  use les_mpi_dynamic_mod, only: compute_mpi_dynamic_les_boundary_rhs_3d,advance_mpi_dynamic_les_ssprk2_3d
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_boundary_rhs_3d,advance_dynamic_les_ssprk2_3d
  use les_mpi_boundary_mod, only: fill_mpi_les_boundary_ghosts_3d
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d,fill_pressure_outlet_ghosts_3d
  implicit none
  integer, parameter :: cells(3)=[5,3,2],width=7,nv=12
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: full(nv,5,3,2),temperature(5,3,2),q(nv),sound,walls(3,6),pressure(6),inlet(nv,6),inlet_t(6)
  real(dp) :: reference(nv,19,17,16),reference_t(19,17,16)
  real(dp), allocatable :: local(:,:,:,:),local_t(:,:,:),ghost(:,:,:,:),gt(:,:,:),saved(:,:,:,:),saved_t(:,:,:)
  real(dp) :: rhs_full(nv,5,3,2),serial_rate(nv),rate(nv),diffusion,serial_diffusion,h(3),wall_t(6)
  real(dp) :: serial_u(nv,5,3,2),serial_t(5,3,2),serial_budget(nv),budget(nv),saved_budget(nv),scale
  real(dp), allocatable :: rhs(:,:,:,:),u(:,:,:,:),t(:,:,:)
  character(len=32) :: faces(6)
  integer :: rank,nranks,ierr,first(3),extent(3),i,j,k,mode,lo,hi,requested_width
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,ok)
  call require(ok,'partition')
  do k=1,2
    do j=1,3
      do i=1,5
        q=0;q(1)=.3_dp;q(2:4)=[real(i,dp),real(j,dp),real(k,dp)];q(5)=101325+100*i;q(12)=1
        call reactive_primitive_to_conserved(species,q,full(:,i,j,k),temperature(i,j,k),sound,ok)
        call require(ok,'initial state')
      end do
    end do
  end do
  lo=first(1)+1;hi=first(1)+extent(1)
  local=full(:,lo:hi,:,:);local_t=temperature(lo:hi,:,:)
  inlet=spread(full(:,1,1,1),2,6);inlet_t=temperature(1,1,1)
  allocate(rhs,mold=local)
  h=[.001_dp,.0015_dp,.002_dp];wall_t=1100
  do mode=1,6
    walls=0;pressure=0;faces='periodic'
    select case(mode)
    case(2)
      faces='outflow';faces(1)='inflow'
    case(3)
      faces='slip_wall'
    case(4)
      faces='no_slip_wall';walls(2,1)=1;walls(2,2)=3;walls(1,3)=2
    case(5)
      faces='isothermal_wall';walls(2,1)=1;walls(2,2)=3;walls(1,3)=2
    case(6)
      faces='outflow';faces(3:4)='no_slip_wall';walls(1,3)=2;pressure(2)=101000
    end select
    call fill_reactive_boundary_ghosts_3d(full,temperature,faces,reference,reference_t,ok,viscous_wall=.true., &
      inflow_state=inlet,inflow_temperature=inlet_t,wall_velocity=walls,ghost_width=width)
    call require(ok,'serial ordinary donors')
    call fill_pressure_outlet_ghosts_3d(species,cells,faces,pressure,reference,reference_t,ok,ghost_width=width)
    call require(ok,'serial pressure donors')
    call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,local,local_t,cells,faces,width,ghost,gt,ok, &
      inlet,inlet_t,walls,pressure)
    call require(ok,'distributed boundary donors')
    call require(all(ghost==reference(:,lo:hi+2*width,:,:)),'exact full-width serial states')
    call require(all(gt==reference_t(lo:hi+2*width,:,:)),'exact full-width serial temperatures')
    call compute_dynamic_les_boundary_rhs_3d(species,full,temperature,h,'box',2,rhs_full,serial_diffusion,ok, &
      boundary_faces=faces,boundary_rate=serial_rate,inflow_state=inlet,inflow_temperature=inlet_t, &
      wall_temperature=wall_t,wall_velocity=walls,outlet_pressure=pressure)
    call require(ok,'serial physical RHS')
    call compute_mpi_dynamic_les_boundary_rhs_3d(MPI_COMM_WORLD,species,local,local_t,h,cells,'box',2, &
      rhs,diffusion,ok,boundary_faces=faces,boundary_rate=rate,inflow_state=inlet,inflow_temperature=inlet_t, &
      wall_temperature=wall_t,wall_velocity=walls,outlet_pressure=pressure)
    call require(ok,'distributed physical RHS')
    call require(all(rhs==rhs_full(:,lo:hi,:,:)),'exact serial RHS')
    scale=max(maxval(abs(serial_rate)),1e-30_dp)
    call require(maxval(abs(rate-serial_rate))<=1e-12_dp*scale,'global boundary rate')
    if(mode==2.or.mode>=5) then
      serial_u=full;serial_t=temperature;u=local;t=local_t
      call advance_dynamic_les_ssprk2_3d(species,serial_u,serial_t,h,'box',2,1e-9_dp,ok, &
        boundary_faces=faces,boundary_change=serial_budget,inflow_state=inlet,inflow_temperature=inlet_t, &
        wall_temperature=wall_t,wall_velocity=walls,outlet_pressure=pressure)
      call require(ok,'serial physical step')
      call advance_mpi_dynamic_les_ssprk2_3d(MPI_COMM_WORLD,species,u,t,h,cells,'box',2,1e-9_dp,ok, &
        boundary_faces=faces,boundary_change=budget,inflow_state=inlet,inflow_temperature=inlet_t, &
        wall_temperature=wall_t,wall_velocity=walls,outlet_pressure=pressure)
      call require(ok,'distributed physical step')
      call require(all(u==serial_u(:,lo:hi,:,:)).and.all(t==serial_t(lo:hi,:,:)),'exact serial SSPRK2')
      scale=max(maxval(abs(serial_budget)),1e-30_dp)
      call require(maxval(abs(budget-serial_budget))<=1e-12_dp*scale,'integrated boundary budget')
      saved_budget=budget
      call advance_mpi_dynamic_les_ssprk2_3d(MPI_COMM_WORLD,species,u,t,h,cells,'box',2,1._dp,ok, &
        boundary_faces=faces,boundary_change=budget,inflow_state=inlet,inflow_temperature=inlet_t, &
        wall_temperature=wall_t,wall_velocity=walls,outlet_pressure=pressure)
      call require(.not.ok,'excessive dt rejects')
      call require(all(u==serial_u(:,lo:hi,:,:)).and.all(t==serial_t(lo:hi,:,:)).and. &
        all(budget==saved_budget),'step rollback')
    end if
  end do
  saved=ghost;saved_t=gt;requested_width=width
  if(rank==nranks-1) requested_width=width-1
  call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,local,local_t,cells,faces,requested_width,ghost,gt,ok, &
    inlet,inlet_t,walls,pressure)
  call require(.not.ok.and.all(ghost==saved).and.all(gt==saved_t),'width disagreement rollback')
  if(rank==nranks-1) walls(1,3)=3
  call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,local,local_t,cells,faces,width,ghost,gt,ok, &
    inlet,inlet_t,walls,pressure)
  call require(.not.ok.and.all(ghost==saved).and.all(gt==saved_t),'wall disagreement rollback')
  walls(1,3)=2
  if(rank==nranks-1) pressure(2)=102000
  call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,local,local_t,cells,faces,width,ghost,gt,ok, &
    inlet,inlet_t,walls,pressure)
  call require(.not.ok.and.all(ghost==saved).and.all(gt==saved_t),'pressure disagreement rollback')
  rhs=123;diffusion=456;rate=789
  pressure=0
  if(rank==nranks-1) wall_t(3)=1200
  call compute_mpi_dynamic_les_boundary_rhs_3d(MPI_COMM_WORLD,species,local,local_t,h,cells,'box',2, &
    rhs,diffusion,ok,boundary_faces=faces,boundary_rate=rate,wall_temperature=wall_t,wall_velocity=walls)
  call require(.not.ok.and.all(rhs==123).and.diffusion==456.and.all(rate==789),'RHS disagreement rollback')
  if(rank==0) print *, 'PASS: distributed physical donors, Dynamic fluxes and SSPRK2 match serial with rollback'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer :: status
    if(condition) return
    print *,rank,'FAIL: '//message
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
end program
