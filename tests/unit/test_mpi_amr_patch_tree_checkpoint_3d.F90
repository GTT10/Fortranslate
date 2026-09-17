program test_mpi_amr_patch_tree_checkpoint_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_checkpoint_3d_mod, only: &
    amr_tree_checkpoint_3d,read_amr_tree_checkpoint_3d,read_amr_tree_checkpoint_header_3d
  use mpi_amr_patch_tree_checkpoint_3d_mod, only: &
    write_mpi_amr_tree_checkpoint_3d,read_mpi_amr_tree_checkpoint_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_tree_checkpoint_3d) :: metadata,restored,saved,header_only
  type(amr_patch_3d) :: patches(1)
  type(amr_tree_state_3d) :: state,temp
  real(dp) :: q(12),u(12),t,sound
  character(len=128) :: path
  character(len=8) :: context
  character(len=16) :: line
  character(len=28*13-1) :: corrupt_row
  integer :: rank,nranks,ierr,unit,owners(2),mode,node
  integer, allocatable :: read_owners(:),saved_owners(:)
  integer(MPI_OFFSET_KIND) :: cell_position,before_size,after_size
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call initialize_amr_patch_3d(2,2,2,1,1,1,1,1,1,2,patches(1),ok)
  call require(ok,'patch')
  call initialize_amr_patch_tree_3d([2,2,2],[1],patches,metadata%tree,ok)
  call require(ok,'tree')
  allocate(metadata%initial(12),metadata%boundary(12));metadata%initial=0;metadata%boundary=0
  owners=[0,nranks-1]
  allocate(state%fields(2),temp%fields(2))
  q=0;q(1)=.3_dp;q(5)=101325;q(6)=.03_dp;q(9)=.22_dp;q(12)=.75_dp
  call reactive_primitive_to_conserved(species,q,u,t,sound,ok)
  call require(ok,'state')
  do node=1,2
    if (rank/=owners(node)) cycle
    allocate(state%fields(node)%values(12,2,2,2),temp%fields(node)%values(1,2,2,2))
    state%fields(node)%values=spread(spread(spread(u,2,2),3,2),4,2)
    temp%fields(node)%values=t
  end do
  write(path,'(a,i0,a)') 'mpi_tree_checkpoint_preflight_',nranks,'.chk'
  if (rank==0) then
    open(newunit=unit,file=trim(path),status='replace')
    write(unit,'(a)') 'preserve'
    close(unit)
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  do mode=1,4
    context='minimal'//achar(10)
    if (mode==1.and.rank==owners(2)) temp%fields(2)%values=1.1_dp*t
    if (mode==2.and.rank==0) then
      ! Keep EOS-consistent U/T, but make one covered parent cell disagree with its children.
      state%fields(1)%values(:,1,1,1)=1.01_dp*u
    end if
    if (mode==3) metadata%steps=1
    if (mode==4.and.rank==owners(2)) context='changed'//achar(10)
    call write_mpi_amr_tree_checkpoint_3d(trim(path),species,context,metadata,owners,state,temp,MPI_COMM_WORLD,ok)
    call require(.not.ok,'collective snapshot rejection')
    metadata%steps=0
    if (rank==owners(2)) temp%fields(2)%values=t
    if (rank==0) state%fields(1)%values(:,1,1,1)=u
    if (rank==0) then
      open(newunit=unit,file=trim(path),status='old')
      read(unit,'(a)') line
      call require(trim(line)=='preserve','existing output preserved')
      read(unit,'(a)',iostat=ierr) line
      call require(ierr<0,'existing output length preserved')
      close(unit)
    end if
  end do
  context='minimal'//achar(10)
  call write_mpi_amr_tree_checkpoint_3d(trim(path),species,context,metadata,owners,state,temp,MPI_COMM_WORLD,ok)
  call require(ok,'owner-local checkpoint write')
  if (rank==0) then
    call read_amr_tree_checkpoint_3d(trim(path),species,context,[2,2,2],restored,ok)
    call require(ok,'serial reader accepts distributed checkpoint')
    do node=1,2
      call require(all(restored%state%fields(node)%values==spread(spread(spread(u,2,2),3,2),4,2)),'exact state')
      call require(all(restored%temperature%fields(node)%values==t),'exact temperature')
    end do
    print *, 'PASS: MPI checkpoint preflight rollback, EOS/averages/context, serial interchange and empty owners'
  end if
  call read_mpi_amr_tree_checkpoint_3d(trim(path),species,context,[2,2,2],'weighted',2, &
    restored,read_owners,MPI_COMM_WORLD,ok)
  call require(ok,'distributed indexed read')
  call require(all(read_owners==[1,0]),'changed read ownership')
  do node=1,2
    call require(allocated(restored%state%fields(node)%values).eqv.(read_owners(node)==rank),'owned read payload only')
    if (read_owners(node)/=rank) cycle
    call require(all(restored%state%fields(node)%values==spread(spread(spread(u,2,2),3,2),4,2)),'read state')
    call require(all(restored%temperature%fields(node)%values==t),'read temperature')
  end do
  saved=restored;saved_owners=read_owners
  if (rank==0) then
    inquire(file=trim(path),size=before_size)
    open(newunit=unit,file=trim(path),status='old',access='stream',form='formatted',action='read')
    call read_amr_tree_checkpoint_header_3d(unit,species,context,[2,2,2],header_only,ok)
    call require(ok,'corruption header position')
    inquire(unit=unit,pos=cell_position)
    close(unit)
    write(corrupt_row,'(*(es27.18e3,1x))') u,1.1_dp*t
    open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',action='write')
    write(unit,pos=cell_position) corrupt_row//achar(10)
    close(unit)
    inquire(file=trim(path),size=after_size)
    call require(before_size==after_size,'corrupt only temperature without truncation')
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call read_mpi_amr_tree_checkpoint_3d(trim(path),species,context,[2,2,2],'round_robin',0, &
    restored,read_owners,MPI_COMM_WORLD,ok)
  call require(.not.ok,'late EOS failure rejects read candidate')
  call require(all(read_owners==saved_owners),'read owner rollback')
  do node=1,2
    if (read_owners(node)/=rank) cycle
    call require(all(restored%state%fields(node)%values==saved%state%fields(node)%values),'read state rollback')
    call require(all(restored%temperature%fields(node)%values==saved%temperature%fields(node)%values),'read T rollback')
  end do
  if (rank==0) print *, 'PASS: distributed checkpoint read, changed ownership and late validation rollback'
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
