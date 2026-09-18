program test_les_mpi_dynamic
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use spray_partition_mod, only: spray_cartesian_block
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_coefficients_3d
  use les_mpi_dynamic_mod, only: compute_mpi_dynamic_les_coefficients_3d
  implicit none
  integer, parameter :: cells(3)=[5,3,2]
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: full(15,5,3,2),temps(5,3,2),reference(3,5,3,2),q(15),sound,spacing(3)
  real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:),coeff(:,:,:,:),saved(:,:,:,:)
  integer :: rank,nranks,ierr,i,j,k,first(3),extent(3),lo,hi,m,ratio,grid
  character(len=8) :: model
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermo')
  call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,ok)
  call require(ok,'partition')
  lo=first(1)+1;hi=lo+extent(1)-1;spacing=[.001_dp,.002_dp,.003_dp]
  do k=1,2
    do j=1,3
      do i=1,5
        q=0;q(1)=.3_dp*(1+.1_dp*sin(real(i,dp)));q(5)=101325*(1+.05_dp*cos(real(j,dp)));q(15)=1
        q(2)=20*sin(real(j,dp));q(3)=15*cos(real(i,dp));q(4)=5*sin(real(i+k,dp))
        call reactive_primitive_to_conserved(species,q,full(:,i,j,k),temps(i,j,k),sound,ok)
        call require(ok,'state')
      end do
    end do
  end do
  state=full(:,lo:hi,:,:);temperature=temps(lo:hi,:,:)
  allocate(coeff(3,extent(1),3,2))
  do m=1,2
    model='box';ratio=2;grid=1
    if (m==2) then
      model='gaussian';ratio=4;grid=2
    end if
    call compute_dynamic_les_coefficients_3d(species,full,temps,spacing,model,ratio,reference,ok,grid)
    call require(ok,'serial dynamic coefficient field')
    call compute_mpi_dynamic_les_coefficients_3d(MPI_COMM_WORLD,species,state,temperature,spacing,cells, &
      model,ratio,coeff,ok,grid)
    call require(ok.and.all(coeff==reference(:,lo:hi,:,:)),'exact serial dynamic coefficient parity')
    call require(maxval(reference)>1e-6_dp,'nontrivial coefficients')
  end do
  call require(all(state==full(:,lo:hi,:,:)).and.all(temperature==temps(lo:hi,:,:)),'inputs preserved')
  saved=coeff
  grid=1
  if (rank==nranks-1) grid=2
  call compute_mpi_dynamic_les_coefficients_3d(MPI_COMM_WORLD,species,state,temperature,spacing,cells, &
    'box',2,coeff,ok,grid)
  call require(.not.ok.and.all(coeff==saved),'different grid widths collective rollback')
  model='box'
  if (rank==nranks-1) model='invalid'
  call compute_mpi_dynamic_les_coefficients_3d(MPI_COMM_WORLD,species,state,temperature,spacing,cells, &
    model,2,coeff,ok)
  call require(.not.ok.and.all(coeff==saved),'invalid remote model collective rollback')
  if (rank==nranks-1) temperature(1,1,1)=ieee_value(1._dp,ieee_quiet_nan)
  call compute_mpi_dynamic_les_coefficients_3d(MPI_COMM_WORLD,species,state,temperature,spacing,cells, &
    'box',2,coeff,ok)
  call require(.not.ok.and.all(coeff==saved),'remote EOS failure collective rollback')
  if (rank==0) print *, 'PASS: distributed dynamic LES coefficient fields'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer :: abort_error
    if (condition) return
    print *,rank,'FAIL: ',message
    call MPI_Abort(MPI_COMM_WORLD,1,abort_error)
  end subroutine
end program
