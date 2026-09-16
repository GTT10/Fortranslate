program test_mpi_amr_tagging_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_reactive_3d_mod, only: temperature_tag_patch_3d
  use mpi_amr_sparse_reactive_3d_mod, only: mpi_amr_sparse_distribution_3d, &
    initialize_mpi_amr_sparse_distribution_3d,temperature_tag_mpi_amr_sparse_patch_3d
  implicit none
  type(amr_patch_3d) :: patch,expected,actual
  type(mpi_amr_sparse_distribution_3d) :: distribution
  real(dp) :: global_t(8,8,8),threshold
  real(dp), allocatable :: local_t(:,:,:)
  character(len=32), allocatable :: faces(:)
  integer :: ierr,rank,mode,case_index,buffer,first,last
  logical :: ok,serial_ok,changed,serial_changed

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call initialize_amr_patch_3d(8,8,8,2,3,2,3,2,3,2,patch,ok)
  call require(ok,'initial patch')
  call initialize_mpi_amr_sparse_distribution_3d(patch,MPI_COMM_WORLD,distribution,ok)
  call require(ok,'distribution')
  first=distribution%coarse_first;last=distribution%coarse_last
  global_t=300
  global_t(1,2,3)=800;global_t(8,6,4)=900
  global_t(3,4,2)=750;global_t(4,1,8)=850
  local_t=global_t(first:last,:,:)
  do mode=1,3
    if (allocated(faces)) deallocate(faces)
    if (mode>1) then
      allocate(faces(6));faces='periodic';faces(1:2)='outflow'
      if (mode==3) faces='outflow'
    end if
    do case_index=1,4
      threshold=750;buffer=case_index-1
      if (case_index==3) buffer=8
      if (case_index==4) threshold=2000
      call temperature_tag_patch_3d(patch,global_t,threshold,buffer,expected,serial_changed,serial_ok,faces)
      call temperature_tag_mpi_amr_sparse_patch_3d(distribution,patch,local_t,threshold,buffer, &
        actual,changed,ok,faces)
      call require(ok.and.serial_ok,'valid tag calls')
      call require(all(bounds(actual)==bounds(expected)).and.(changed.eqv.serial_changed),'serial tag equality')
      if (case_index==4) call require(.not.changed,'empty tags preserve patch')
    end do
  end do
  call temperature_tag_mpi_amr_sparse_patch_3d(distribution,expected,local_t,2000.0_dp,0,actual,changed,ok,faces)
  call require(ok.and..not.changed,'empty repeated tag')
  local_t=300
  do mode=max(first,2),min(last,3)
    local_t(mode-first+1,2:3,2:3)=750
  end do
  call temperature_tag_mpi_amr_sparse_patch_3d(distribution,patch,local_t,750.0_dp,0,actual,changed,ok,faces)
  call require(ok.and..not.changed.and.all(bounds(actual)==bounds(patch)),'unchanged nonempty tags')
  local_t=global_t(first:last,:,:)
  threshold=750
  if (rank==0) threshold=751
  call reject(threshold,1)
  buffer=1
  if (rank==0) buffer=2
  call reject(750.0_dp,buffer)
  if (rank==0) local_t(1,1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
  call reject(750.0_dp,1)
  local_t=global_t(first:last,:,:)
  if (rank==0) faces(1:2)='periodic'
  call reject(750.0_dp,1)
  faces='outflow'
  if (rank==0) faces(1)='invalid'
  call reject(750.0_dp,1)
  if (rank==0) print '(a)','test_mpi_amr_tagging_3d: PASS'
  call MPI_Finalize(ierr)
contains
  subroutine reject(value,padding)
    real(dp), intent(in) :: value
    integer, intent(in) :: padding
    call temperature_tag_mpi_amr_sparse_patch_3d(distribution,patch,local_t,value,padding,actual,changed,ok,faces)
    call require(.not.ok.and..not.changed.and.all(bounds(actual)==bounds(patch)),'collective rejection preserves patch')
  end subroutine reject
  function bounds(p) result(values)
    type(amr_patch_3d), intent(in) :: p
    integer :: values(10)
    values=[p%coarse_nx,p%coarse_ny,p%coarse_nz,p%coarse_i_lower,p%coarse_i_upper, &
      p%coarse_j_lower,p%coarse_j_upper,p%coarse_k_lower,p%coarse_k_upper,p%refinement_ratio]
  end function bounds
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer :: status
    if (condition) return
    print '(a,i0,2a)','rank ',rank,': ',message
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine require
end program test_mpi_amr_tagging_3d
