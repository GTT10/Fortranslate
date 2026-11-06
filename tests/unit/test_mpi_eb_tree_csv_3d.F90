program test_mpi_eb_tree_csv_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use mpi_eb_tree_csv_3d_mod, only: write_mpi_eb_tree_csv_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(0)
  type(amr_tree_state_3d) :: state,temp,saved_heat,saved_held
  type(amr_tree_state_3d), allocatable :: heat,held
  type(eb_geometry_3d) :: geometry(1),original
  real(dp) :: q(12),u(12),t,sound,value
  character(len=128) :: path,other
  character(len=2048) :: line
  integer :: rank,nranks,ierr,unit,owners(1),mode,lines
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call initialize_amr_patch_tree_3d([2,2,2],[integer::],patches,tree,ok)
  call require(ok,'tree')
  call build_axis_plane_eb_geometry_3d(2,2,2,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.61_dp,geometry(1),ok)
  call require(ok,'geometry')
  original=geometry(1)
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
  allocate(heat)
  allocate(heat%fields(1))
  if(rank==owners(1)) allocate(heat%fields(1)%values(1,2,2,2),source=-1.e9_dp)
  saved_heat=heat
  allocate(held)
  held=heat
  if(rank==owners(1)) held%fields(1)%values=-314159._dp
  saved_held=held
  write(path,'(a,i0,a)') 'mpi_eb_tree_csv_preflight_',nranks,'.csv'
  if (rank==0) then
    open(newunit=unit,file=trim(path),status='replace')
    write(unit,'(a)') 'preserve'
    close(unit)
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  do mode=1,8
    other=path
    if (mode==1.and.rank==owners(1)) state%fields(1)%values(1,2,2,2)=ieee_value(0._dp,ieee_quiet_nan)
    if (mode==2.and.rank==owners(1)) species(1)%name='mismatch'
    if (mode==3.and.rank==owners(1)) other='different.csv'
    if (mode==4.and.rank==owners(1)) geometry(1)%x_lower=.1_dp
    if(mode==5.and.rank==owners(1)) heat%fields(1)%values(1,2,2,2)=ieee_value(0._dp,ieee_quiet_nan)
    if(mode==6.and.rank==0) deallocate(heat)
    if(mode==7.and.rank==owners(1)) held%fields(1)%values(1,2,2,2)=ieee_value(0._dp,ieee_quiet_nan)
    if(mode==8.and.rank==0) deallocate(held)
    call write_mpi_eb_tree_csv_3d(trim(other),tree,owners,species,state,temp,geometry,MPI_COMM_WORLD,ok,heat,held)
    if(.not.allocated(held)) allocate(held)
    held=saved_held
    if(.not.allocated(heat)) allocate(heat)
    heat=saved_heat
    call require(.not.ok,'collective preflight rejection')
    geometry(1)=original
    call load_h2o2_elementary_thermo(species,ok)
    call require(ok,'restore species')
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
  ! Solid payloads are intentionally ignored, including nonfinite temperatures.
  if (rank==owners(1)) then
    state%fields(1)%values(:,:,:,1)=ieee_value(0._dp,ieee_quiet_nan)
    temp%fields(1)%values(:,:,:,1)=ieee_value(0._dp,ieee_quiet_nan)
    heat%fields(1)%values(:,:,:,1)=ieee_value(0._dp,ieee_quiet_nan)
    held%fields(1)%values(:,:,:,1)=ieee_value(0._dp,ieee_quiet_nan)
  end if
  call write_mpi_eb_tree_csv_3d(trim(path),tree,owners,species,state,temp,geometry,MPI_COMM_WORLD,ok)
  call require(ok,'owner-local output with empty root and solid NaNs')
  if (rank==0) then
    open(newunit=unit,file=trim(path),status='old')
    lines=0
    do
      read(unit,'(a)',iostat=ierr) line
      if (ierr/=0) exit
      lines=lines+1
    end do
    close(unit)
    call require(lines==5,'header and four fluid cells only')
    print *, 'PASS: EB CSV preflight preserves output; empty ranks and solid payloads supported'
  end if
  call write_mpi_eb_tree_csv_3d(trim(path),tree,owners,species,state,temp,geometry,MPI_COMM_WORLD,ok,heat)
  call require(ok,'signed heat CSV ignores solid NaNs')
  if(rank==0) then
    open(newunit=unit,file=trim(path),status='old')
    read(unit,'(a)') line
    call require(index(line,',heat_release_W_m3')>0,'heat column includes units')
    do mode=1,4
      read(unit,'(a)',iostat=ierr) line
      call require(ierr==0,'heat row present')
      read(line(index(trim(line),',',back=.true.)+1:),*,iostat=ierr) value
      call require(ierr==0.and.value==-1.e9_dp,'signed stage heat retained')
    end do
    close(unit)
  end if
  call write_mpi_eb_tree_csv_3d(trim(path),tree,owners,species,state,temp,geometry,MPI_COMM_WORLD,ok,heat,held)
  call require(ok,'held and heat CSV ignores solid NaNs')
  if(rank==0) then
    open(newunit=unit,file=trim(path),status='old')
    read(unit,'(a)') line
    call require(index(line,',heat_release_W_m3,rho_internal_energy_J_m3')>0,'held column identity')
    do mode=1,4
      read(unit,'(a)',iostat=ierr) line
      call require(ierr==0,'held row present')
      read(line(index(trim(line),',',back=.true.)+1:),*,iostat=ierr) value
      call require(ierr==0.and.value==-314159._dp,'independent signed held scalar retained')
    enddo
    close(unit)
  endif
  call write_mpi_eb_tree_csv_3d(trim(path),tree,owners,species,state,temp,geometry,MPI_COMM_WORLD,ok,internal_energy=held)
  call require(ok,'held CSV without heat')
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
