program test_mpi_amr_patch_tree_csv_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use mpi_amr_patch_tree_csv_3d_mod, only: write_mpi_amr_tree_csv_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(0)
  type(amr_tree_state_3d) :: state,temp
  real(dp) :: q(12),u(12),t,sound,time
  character(len=128) :: path,other
  character(len=16) :: line
  integer :: rank,nranks,ierr,unit,owners(1),mode
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call initialize_amr_patch_tree_3d([2,2,2],[integer::],patches,tree,ok)
  call require(ok,'tree')
  owners=nranks-1
  allocate(state%fields(1),temp%fields(1))
  if (rank==owners(1)) then
    q=0;q(1)=.3_dp;q(5)=101325;q(6)=.03_dp;q(9)=.22_dp;q(12)=.75_dp
    call reactive_primitive_to_conserved(species,q,u,t,sound,ok)
    call require(ok,'state')
    allocate(state%fields(1)%values(12,2,2,2),temp%fields(1)%values(1,2,2,2))
    state%fields(1)%values=spread(spread(spread(u,2,2),3,2),4,2)
    temp%fields(1)%values=t
  end if
  write(path,'(a,i0,a)') 'mpi_tree_csv_preflight_',nranks,'.csv'
  if (rank==0) then
    open(newunit=unit,file=trim(path),status='replace')
    write(unit,'(a)') 'preserve'
    close(unit)
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  do mode=1,3
    time=0;other=path
    if (mode==1.and.rank==owners(1)) state%fields(1)%values(1,2,2,2)=-1
    if (mode==2.and.rank==owners(1)) time=1
    if (mode==3.and.rank==owners(1)) other='different.csv'
    call write_mpi_amr_tree_csv_3d(trim(other),tree,owners,species,state,temp,[0._dp,0._dp,0._dp], &
      [.001_dp,.001_dp,.001_dp],time,MPI_COMM_WORLD,ok)
    call require(.not.ok,'collective preflight rejection')
    if (rank==owners(1)) state%fields(1)%values(1,2,2,2)=u(1)
    if (rank==0) then
      open(newunit=unit,file=trim(path),status='old')
      read(unit,'(a)') line
      call require(trim(line)=='preserve','existing output preserved')
      read(unit,'(a)',iostat=ierr) line
      call require(ierr<0,'existing output length preserved')
      close(unit)
    end if
  end do
  if (rank==0) print *, 'PASS: CSV invalid owner state, time/path mismatch preserve existing output'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: status
    if (condition) return
    print *, 'FAIL rank ',rank,': '//label
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
end program
