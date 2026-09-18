program test_les_mpi_models
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use les_smagorinsky_3d_mod, only: advance_les_ssprk2_3d
  use spray_mpi_mod, only: advance_mpi_spray_les
  use spray_partition_mod, only: spray_cartesian_block
  implicit none
  type(nasa7_species), allocatable :: species(:)
  integer, parameter :: cells(3)=[8,4,4]
  real(dp) :: initial(15,8,4,4),initial_t(8,4,4),reference(15,8,4,4),reference_t(8,4,4)
  real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:),saved(:,:,:,:),saved_t(:,:,:)
  real(dp) :: q(15),sound,change(15),serial_change(15),wall_t(6),wall_v(3,6),ci,cs,trial_ci
  real(dp), parameter :: spacing(3)=[.001_dp,.001_dp,.001_dp],dt=1e-7_dp
  character(len=16), parameter :: models(3)=[character(len=16)::'smagorinsky','wale','vreman']
  character(len=16) :: model,faces(6)
  integer :: rank,nranks,ierr,i,j,k,m,b,first(3),extent(3),lo,hi
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermo')
  call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,ok)
  call require(ok,'partition')
  lo=first(1)+1;hi=lo+extent(1)-1
  do k=1,4
    do j=1,4
      do i=1,8
        q=0;q(1)=.3_dp;q(5)=101325;q(15)=1
        q(2)=10*sin(2*acos(-1._dp)*(j-.5_dp)/4)
        q(3)=10*sin(2*acos(-1._dp)*(i-.5_dp)/8)
        call reactive_primitive_to_conserved(species,q,initial(:,i,j,k),initial_t(i,j,k),sound,ok)
        call require(ok,'initial state')
      end do
    end do
  end do
  wall_t=900;wall_v=0;wall_v(1,3:4)=2
  do m=1,9
    model=models(mod(m-1,3)+1)
    ci=0;cs=.17_dp
    if (m>3) ci=.1_dp
    if (m>6) cs=0
    do b=1,2
      reference=initial;reference_t=initial_t
      state=initial(:,lo:hi,:,:);temperature=initial_t(lo:hi,:,:)
      change=0;serial_change=0
      if (b==1) then
        call advance_les_ssprk2_3d(species,reference,reference_t,spacing,cs,.7_dp,.7_dp,dt,ok,model=model,ci=ci)
        call require(ok,'serial periodic')
        call advance_mpi_spray_les(MPI_COMM_WORLD,species,state,temperature,spacing,cells,dt,cs,.7_dp,.7_dp, &
          ok,model=model,ci=ci)
      else
        faces='outflow';faces(3:4)='isothermal_wall'
        call advance_les_ssprk2_3d(species,reference,reference_t,spacing,cs,.7_dp,.7_dp,dt,ok,faces,serial_change, &
          wall_temperature=wall_t,wall_velocity=wall_v,model=model,ci=ci)
        call require(ok,'serial walls')
        call advance_mpi_spray_les(MPI_COMM_WORLD,species,state,temperature,spacing,cells,dt,cs,.7_dp,.7_dp, &
          ok,faces,change,wall_temperature=wall_t,wall_velocity=wall_v,model=model,ci=ci)
      end if
      call require(ok,'MPI LES')
      call require(maxval(abs(state-reference(:,lo:hi,:,:))/max(1._dp,abs(reference(:,lo:hi,:,:))))<2e-12_dp, &
        'state parity')
      call require(maxval(abs(temperature-reference_t(lo:hi,:,:)))<1e-8_dp,'temperature parity')
      call require(all(abs(change-serial_change)<1e-10_dp*max(1e-12_dp,abs(serial_change))),'budget parity')
      saved=state;saved_t=temperature
      trial_ci=ci
      if (rank==nranks-1) trial_ci=ci+.01_dp
      call advance_mpi_spray_les(MPI_COMM_WORLD,species,state,temperature,spacing,cells,0._dp,0._dp,.7_dp,.7_dp, &
        ok,model=model,ci=trial_ci)
      call require(.not.ok.and.all(state==saved).and.all(temperature==saved_t),'CI mismatch at disabled step')
      trial_ci=ci
      if (rank==nranks-1) trial_ci=ieee_value(0._dp,ieee_quiet_nan)
      call advance_mpi_spray_les(MPI_COMM_WORLD,species,state,temperature,spacing,cells,dt,cs,.7_dp,.7_dp, &
        ok,model=model,ci=trial_ci)
      call require(.not.ok.and.all(state==saved).and.all(temperature==saved_t),'CI NaN collective rollback')
      if (ci>0) then
        call require(any(state/=initial(:,lo:hi,:,:)),'isotropic evolution active')
        call advance_mpi_spray_les(MPI_COMM_WORLD,species,state,temperature,spacing,cells,1._dp,cs,.7_dp,.7_dp, &
          ok,model=model,ci=ci)
        call require(.not.ok.and.all(state==saved).and.all(temperature==saved_t),'CI excessive dt rollback')
      end if
      if (rank==nranks-1) model='invalid'
      call advance_mpi_spray_les(MPI_COMM_WORLD,species,state,temperature,spacing,cells,0._dp,0._dp,.7_dp,.7_dp, &
        ok,model=model,ci=ci)
      call require(.not.ok.and.all(state==saved).and.all(temperature==saved_t),'invalid model collective rollback')
      model=models(mod(m-1,3)+1)
      if (rank==nranks-1) model=models(mod(m,3)+1)
      call advance_mpi_spray_les(MPI_COMM_WORLD,species,state,temperature,spacing,cells,dt,cs,.7_dp,.7_dp, &
        ok,model=model,ci=ci)
      call require(.not.ok.and.all(state==saved).and.all(temperature==saved_t),'mismatched model collective rollback')
      model=models(mod(m-1,3)+1)
    end do
  end do
  if (rank==0) print *, 'PASS: three LES models, CI variants, periodic/wall parity, collective rollback'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: status
    if (condition) return
    print *, 'FAIL rank ',rank,': ',label
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
end program
