program test_mpi_amr_regrid_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d,prolong_pcm_3d
  use amr_reactive_3d_mod, only: regrid_amr_reactive_3d
  use mpi_amr_sparse_reactive_3d_mod, only: mpi_amr_sparse_distribution_3d,mpi_amr_sparse_hierarchy_3d, &
    initialize_mpi_amr_sparse_distribution_3d,scatter_mpi_amr_sparse_hierarchy_3d, &
    gather_mpi_amr_sparse_hierarchy_3d,regrid_mpi_amr_sparse_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: patch,serial_patch,next_patch,bad,saved_patch
  type(mpi_amr_sparse_distribution_3d) :: distribution,saved_distribution
  type(mpi_amr_sparse_hierarchy_3d) :: hierarchy,saved
  real(dp) :: c(12,8,8,8),t(8,8,8),q(12),sound,inlet(12,6),it(6)
  real(dp), allocatable :: f(:,:,:,:),ft(:,:,:),gc(:,:,:,:),gt(:,:,:),gf(:,:,:,:),gft(:,:,:)
  character(len=32), allocatable :: faces(:)
  character(len=32) :: method
  integer :: ierr,rank,i,j,k,r,mode,lo,hi
  logical :: ok

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermodynamics')
  do r=2,3
    call initialize_amr_patch_3d(8,8,8,2,3,2,3,2,3,r,patch,ok)
    call require(ok,'initial geometry');serial_patch=patch
    call initialize_mpi_amr_sparse_distribution_3d(patch,MPI_COMM_WORLD,distribution,ok)
    call require(ok,'initial distribution')
    do k=1,8
      do j=1,8
        do i=1,8
          q=0;q(1)=.3_dp+.01_dp*i+.005_dp*j+.002_dp*k;q(5)=101325*q(1)/.3_dp
          q(2)=20;q(3)=-10;q(6)=.03_dp;q(9)=.22_dp;q(12)=.75_dp
          call reactive_primitive_to_conserved(species,q,c(:,i,j,k),t(i,j,k),sound,ok)
          call require(ok,'initial state')
        end do
      end do
    end do
    allocate(f(12,2*r,2*r,2*r),ft(2*r,2*r,2*r))
    call prolong_pcm_3d(c,patch,f,ok);call require(ok,'initial fine state')
    ft=t(1,1,1)
    do i=1,2*r
      f(:,i,:,:)=f(:,i,:,:)*(1.0_dp+.01_dp*(-1.0_dp)**i)
    end do
    call scatter_mpi_amr_sparse_hierarchy_3d(distribution,species,patch,0,c,t,f,ft,hierarchy,ok)
    call require(ok,'scatter')
    do mode=1,5
      method='characteristic_plm';lo=3;hi=5
      if (mode==2) then
        allocate(faces(6));faces='outflow';faces(1)='slip_wall';faces(5)='no_slip_wall'
        lo=1;hi=3
      else if (mode==3) then
        lo=2;hi=3;method='pcm'
      else if (mode==4) then
        lo=6;hi=8
      else if (mode==5) then
        lo=1;hi=8
      end if
      call initialize_amr_patch_3d(8,8,8,lo,hi,lo,hi,lo,hi,r,next_patch,ok)
      call require(ok,'target geometry')
      if (rank==0) then
        call regrid_amr_reactive_3d(species,serial_patch,c,t,f,ft,next_patch,ok,method,'mc',faces)
        call require(ok,'serial remap')
      end if
      call regrid_mpi_amr_sparse_3d(distribution,species,patch,hierarchy,next_patch,ok,method,'mc',faces)
      call require(ok,'distributed remap')
      call gather_mpi_amr_sparse_hierarchy_3d(distribution,patch,0,hierarchy,gc,gt,gf,gft,ok)
      call require(ok,'gather remap')
      if (rank==0) call require(all(c==gc).and.all(t==gt).and.all(f==gf).and.all(ft==gft), &
        'states and temperatures exactly match serial')
      call require(distribution%is_valid(patch),'remapped ownership')
    end do
    saved=hierarchy;saved_distribution=distribution;saved_patch=patch
    bad=patch
    if (rank==0) bad%coarse_i_lower=0
    call regrid_mpi_amr_sparse_3d(distribution,species,patch,hierarchy,bad,ok,'pcm','mc',faces)
    call rollback()
    bad=patch
    if (rank==0) bad%coarse_i_lower=2
    call regrid_mpi_amr_sparse_3d(distribution,species,patch,hierarchy,bad,ok,'pcm','mc',faces)
    call rollback()
    faces(1)='inflow';inlet=spread(hierarchy%coarse_state(:,1,1,1),2,6);it=hierarchy%coarse_temperature(1,1,1)
    ! Make reservoir metadata identical, then invalidate just one rank.
    call MPI_Bcast(inlet,size(inlet),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(it,size(it),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (rank==0) inlet(1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
    bad=patch
    call regrid_mpi_amr_sparse_3d(distribution,species,patch,hierarchy,bad,ok,'pcm','mc',faces,inlet,it)
    call rollback()
    deallocate(f,ft,faces)
  end do
  if (rank==0) print '(a)','test_mpi_amr_regrid_3d: PASS'
  call MPI_Finalize(ierr)
contains
  subroutine rollback()
    call require(.not.ok,'invalid remap rejected')
    call require(all(hierarchy%coarse_state==saved%coarse_state).and. &
      all(hierarchy%coarse_temperature==saved%coarse_temperature).and. &
      all(hierarchy%fine_state==saved%fine_state).and.all(hierarchy%fine_temperature==saved%fine_temperature), &
      'failure preserves hierarchy')
    call require(all(distribution%fine_counts==saved_distribution%fine_counts).and. &
      all(distribution%fine_displacements==saved_distribution%fine_displacements).and. &
      patch%coarse_i_lower==saved_patch%coarse_i_lower.and.patch%coarse_i_upper==saved_patch%coarse_i_upper, &
      'failure preserves geometry and ownership')
  end subroutine rollback
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer :: status
    if (condition) return
    print '(a,i0,2a)','rank ',rank,': ',message
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine require
end program test_mpi_amr_regrid_3d
