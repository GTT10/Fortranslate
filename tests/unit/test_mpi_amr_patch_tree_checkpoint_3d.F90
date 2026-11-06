program test_mpi_amr_patch_tree_checkpoint_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_checkpoint_3d_mod, only: &
    amr_tree_checkpoint_3d,read_amr_tree_checkpoint_3d,read_amr_tree_checkpoint_header_3d,write_amr_tree_checkpoint_3d
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
  logical :: ok,exists
  character(len=16) :: test_mode
  call get_command_argument(1,test_mode)
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
  do mode=1,5
    context='minimal'//achar(10)
    if (mode==1.and.rank==owners(2)) temp%fields(2)%values=1.1_dp*t
    if (mode==2.and.rank==0) then
      ! Keep EOS-consistent U/T, but make one covered parent cell disagree with its children.
      state%fields(1)%values(:,1,1,1)=1.01_dp*u
    end if
    if (mode==3) metadata%steps=1
    if (mode==4.and.rank==owners(2)) context='changed'//achar(10)
    call write_mpi_amr_tree_checkpoint_3d(trim(path),species,context,metadata,owners,state,temp,MPI_COMM_WORLD,ok, &
      mol_energy_reset=mode==5.and.rank==0)
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
  if(test_mode=='publication') then
    if(rank==0) then
      open(newunit=unit,file=trim(path)//'.tmp',status='new');write(unit,'(a)') 'temporary';close(unit)
    endif
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call write_mpi_amr_tree_checkpoint_3d(trim(path),species,context,metadata,owners,state,temp,MPI_COMM_WORLD,ok)
    call require(.not.ok,'temporary collision rejected collectively')
    if(rank==0) then
      open(newunit=unit,file=trim(path),status='old');read(unit,'(a)') line;close(unit)
      call require(trim(line)=='preserve','collision preserves final')
      inquire(file=trim(path),size=ierr)
      call require(ierr==9,'collision preserves final size')
      open(newunit=unit,file=trim(path)//'.tmp',status='old');read(unit,'(a)') line
      call require(trim(line)=='temporary','collision preserves other temporary')
      inquire(unit=unit,size=ierr)
      call require(ierr==10,'collision preserves other temporary size')
      close(unit,status='delete')
      call execute_command_line('mkdir publication-directory',exitstat=ierr)
      call require(ierr==0,'create rename-rejection fixture')
    endif
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call write_mpi_amr_tree_checkpoint_3d('publication-directory',species,context,metadata,owners,state,temp, &
      MPI_COMM_WORLD,ok)
    call require(.not.ok,'rename failure rejected collectively')
    if(rank==0) then
      inquire(file='publication-directory.tmp',exist=exists)
      call require(.not.exists,'rename failure removes own temporary')
      call execute_command_line('rmdir publication-directory',exitstat=ierr)
      call require(ierr==0,'rename failure preserves directory')
      print *, 'PASS: MPI publication collision and rename rollback'
    endif
    call MPI_Finalize(ierr)
    stop
  endif
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
  metadata%time=1e-9_dp;metadata%steps=3
  metadata%previous_timestep=2e-10_dp;metadata%post_regrid=.true.
  if(rank==nranks-1) metadata%previous_timestep=3e-10_dp
  call write_mpi_amr_tree_checkpoint_3d(trim(path),species,context,metadata,owners,state,temp,MPI_COMM_WORLD,ok)
  call require(.not.ok,'rank disagreement in interval rejects publication')
  metadata%previous_timestep=2e-10_dp
  call write_mpi_amr_tree_checkpoint_3d(trim(path),species,context,metadata,owners,state,temp,MPI_COMM_WORLD,ok,indexed=.true.)
  call require(ok,'indexed collective history write')
  call read_mpi_amr_tree_checkpoint_3d(trim(path),species,context,[2,2,2],'round_robin',0, &
    restored,read_owners,MPI_COMM_SELF,ok)
  call require(ok.and.restored%previous_timestep==2e-10_dp.and.restored%post_regrid,'changed-rank history restore')
  restored%previous_timestep=53
  call read_mpi_amr_tree_checkpoint_3d(trim(path),species,'changed'//achar(10),[2,2,2],'round_robin',0, &
    restored,read_owners,MPI_COMM_WORLD,ok)
  call require(.not.ok.and.restored%previous_timestep==53.and.restored%post_regrid,'collective history rollback')
  metadata%time=0;metadata%steps=0;metadata%post_regrid=.false.;deallocate(metadata%previous_timestep)
  call check_internal_energy()
  call check_sdc()
  call MPI_Finalize(ierr)
contains
  subroutine check_internal_energy()
    type(amr_tree_checkpoint_3d) :: payload,actual,prior,aux_header
    type(amr_tree_state_3d) :: retained_t
    real(dp) :: recovered(2,2,2)
    integer :: variant,peer
    integer(MPI_OFFSET_KIND) :: payload_position
    integer, allocatable :: prior_owners(:)
    character(len=27) :: bad_record
    payload=metadata;payload%internal_energy=temp;payload%energy_correction=-7._dp
    retained_t=temp;payload%previous_timestep=0._dp
    do peer=1,size(owners)
      if(owners(peer)/=rank) cycle
      payload%internal_energy%fields(peer)%values(1,:,:,:)=.8_dp*state%fields(peer)%values(5,:,:,:)
      call recover_reactive_temperatures_3d(species,state%fields(peer)%values,temp%fields(peer)%values(1,:,:,:), &
        2,2,2,recovered,ok,internal_energy_density=payload%internal_energy%fields(peer)%values(1,:,:,:))
      call require(ok,'MPI retained EOS fixture')
      retained_t%fields(peer)%values(1,:,:,:)=recovered
    enddo
    do variant=0,7
      if(variant==4) payload%pressure_work_change=-.25_dp
      if(variant==2) then
        payload%source=metadata%initial*0._dp
        payload%heat_release=temp
        do peer=1,size(owners)
          if(owners(peer)/=rank) cycle
          payload%heat_release%fields(peer)%values=-1._dp
        enddo
      endif
      call write_mpi_amr_tree_checkpoint_3d(trim(path),species,context,payload,owners,state,retained_t, &
        MPI_COMM_WORLD,ok,indexed=mod(variant,2)==1)
      call require(ok,'MPI retained checkpoint write')
      if(rank==0) then
        call read_amr_tree_checkpoint_3d(trim(path),species,context,[2,2,2],actual,ok,allow_internal_energy=.true.)
        call require(ok,'serial reads MPI retained checkpoint')
        call require(actual%energy_correction==-7._dp,'serial retained correction exact')
        call require(allocated(actual%pressure_work_change).eqv.(variant>=4),'serial pressure-work presence exact')
        if(variant>=4) call require(actual%pressure_work_change==-.25_dp,'serial signed pressure-work exact')
        call write_amr_tree_checkpoint_3d(trim(path)//'.retained_serial',species,context,actual,ok,indexed=mod(variant,2)==1)
        call require(ok,'retained serial interchange write')
      endif
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call read_mpi_amr_tree_checkpoint_3d(trim(path),species,context,[2,2,2],'weighted',2, &
        actual,read_owners,MPI_COMM_WORLD,ok)
      call require(.not.ok,'unsupported MPI retained reader rejected')
      call read_mpi_amr_tree_checkpoint_3d(trim(path),species,context,[2,2,2],'weighted',2, &
        actual,read_owners,MPI_COMM_WORLD,ok,allow_internal_energy=.true.)
      call require(ok,'MPI retained ownership-change read')
      call require(actual%energy_correction==-7._dp,'MPI retained correction exact')
      call require(allocated(actual%pressure_work_change).eqv.(variant>=4),'MPI pressure-work presence exact')
      if(variant>=4) call require(actual%pressure_work_change==-.25_dp,'MPI signed pressure-work exact')
      do peer=1,size(owners)
        call require(allocated(actual%internal_energy%fields(peer)%values).eqv.(read_owners(peer)==rank), &
          'retained checkpoint owner-local payload')
        if(read_owners(peer)/=rank) cycle
        call require(all(actual%internal_energy%fields(peer)%values==.8_dp*u(5)),'MPI retained scalar exact')
      enddo
      call read_mpi_amr_tree_checkpoint_3d(trim(path)//'.retained_serial',species,context,[2,2,2],'round_robin',0, &
        actual,read_owners,MPI_COMM_SELF,ok,allow_internal_energy=.true.)
      call require(ok.and.all(read_owners==0),'retained serial checkpoint changed-rank read')
      if(variant>=4) call require(actual%pressure_work_change==-.25_dp,'pressure-work changed-rank exact')
      prior=actual
      if(rank==nranks-1) payload%energy_correction=-8._dp
      call write_mpi_amr_tree_checkpoint_3d(trim(path),species,context,payload,owners,state,retained_t,MPI_COMM_WORLD,ok)
      call require(.not.ok,'rank-disagreeing retained correction rejected before replace')
      payload%energy_correction=-7._dp
      if(variant>=4) then
        if(rank==nranks-1) payload%pressure_work_change=-.5_dp
        call write_mpi_amr_tree_checkpoint_3d(trim(path),species,context,payload,owners,state,retained_t,MPI_COMM_WORLD,ok)
        call require(.not.ok,'rank-disagreeing pressure work refuses before replace')
        if(rank==nranks-1) payload%pressure_work_change=ieee_value(0._dp,ieee_quiet_nan)
        call write_mpi_amr_tree_checkpoint_3d(trim(path),species,context,payload,owners,state,retained_t,MPI_COMM_WORLD,ok)
        call require(.not.ok,'rank-local nonfinite pressure work refuses collectively')
        payload%pressure_work_change=-.25_dp
      endif
      call read_mpi_amr_tree_checkpoint_3d(trim(path),species,context,[2,2,2],'weighted',2, &
        actual,read_owners,MPI_COMM_WORLD,ok,allow_internal_energy=.true.)
      call require(ok,'retained file survives rejected collective write')
      if(mod(variant,2)==1) then
        prior=actual;prior_owners=read_owners
        if(rank==0) then
          open(newunit=unit,file=trim(path),status='old',access='stream',form='formatted')
          call read_amr_tree_checkpoint_header_3d(unit,species,context,[2,2,2],aux_header,ok, &
            allow_sdc=.true.,allow_internal_energy=.true.)
          call require(ok,'retained corruption header')
          inquire(unit=unit,pos=payload_position)
          close(unit)
          payload_position=payload_position+aux_header%record_offsets(size(owners)+1)-28
          write(bad_record,'(es27.18e3)') ieee_value(0._dp,ieee_quiet_nan)
          open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',action='write')
          write(unit,pos=payload_position) bad_record
          close(unit)
        endif
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call read_mpi_amr_tree_checkpoint_3d(trim(path),species,context,[2,2,2],'round_robin',0, &
          actual,read_owners,MPI_COMM_WORLD,ok,allow_internal_energy=.true.)
        call require(.not.ok,'late retained scalar corruption rejected')
        if(variant>=4) call require(actual%pressure_work_change==prior%pressure_work_change, &
          'late corruption preserves pressure-work history')
        call require(all(read_owners==prior_owners).and.actual%energy_correction==prior%energy_correction, &
          'retained corrupt read metadata rollback')
        do peer=1,size(owners)
          if(prior_owners(peer)/=rank) cycle
          call require(all(actual%state%fields(peer)%values==prior%state%fields(peer)%values).and. &
            all(actual%temperature%fields(peer)%values==prior%temperature%fields(peer)%values).and. &
            all(actual%internal_energy%fields(peer)%values==prior%internal_energy%fields(peer)%values), &
            'retained corrupt read field rollback')
        enddo
      endif
    enddo
  end subroutine check_internal_energy

  subroutine check_sdc()
    type(MPI_Comm) :: smaller
    integer :: c,variant,color,small_rank
    real(dp) :: corrupt_source(8)
    character(len=28*8-1) :: source_record
    allocate(metadata%source(12),metadata%reaction_source%fields(2))
    metadata%source=[(real(-c,dp)*.01_dp,c=1,12)]
    do node=1,2
      if(owners(node)/=rank) cycle
      allocate(metadata%reaction_source%fields(node)%values(8,2,2,2))
      do c=1,8
        metadata%reaction_source%fields(node)%values(c,:,:,:)=real((-1)**c*c*node,dp)
      end do
    end do
    do variant=1,3
      if(variant==2) then
        allocate(metadata%heat_release%fields(size(owners)))
        do node=1,size(owners)
          if(owners(node)/=rank) cycle
          metadata%heat_release%fields(node)%values=temp%fields(node)%values
          metadata%heat_release%fields(node)%values=real((-1)**node*node,dp)*1.e9_dp
        end do
      end if
      if(variant==3) deallocate(metadata%reaction_source%fields)
      call write_mpi_amr_tree_checkpoint_3d(trim(path),species,context,metadata,owners,state,temp,MPI_COMM_WORLD,ok, &
        indexed=variant==2)
      call require(ok,'MPI SDC checkpoint write')
      if(rank==0) then
        call read_amr_tree_checkpoint_3d(trim(path),species,context,[2,2,2],restored,ok)
        call require(ok,'serial accepts MPI SDC payload')
        call require(all(restored%source==metadata%source),'serial source ledger exact')
        call require(allocated(restored%reaction_source%fields).eqv.(variant/=3),'serial source presence')
        if(variant/=3) then
          do node=1,2
          do c=1,8
            call require(all(restored%reaction_source%fields(node)%values(c,:,:,:)==real((-1)**c*c*node,dp)), &
              'serial signed reaction source exact')
          end do
          end do
        end if
        call write_amr_tree_checkpoint_3d(trim(path)//'.serial',species,context,restored,ok,indexed=variant==2)
        call require(ok,'serial SDC interchange write')
      end if
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      ! File was written by the full communicator; read with two ranks.
      color=MPI_UNDEFINED
      if(rank<2) color=0
      call MPI_Comm_split(MPI_COMM_WORLD,color,rank,smaller,ierr)
      if(rank<2) then
        call MPI_Comm_rank(smaller,small_rank,ierr)
        call read_mpi_amr_tree_checkpoint_3d(trim(path),species,context,[2,2,2],'weighted',2, &
          restored,read_owners,smaller,ok)
        call require(ok,'MPI SDC changed-rank read')
        call verify_sdc(variant,small_rank)
        call read_mpi_amr_tree_checkpoint_3d(trim(path)//'.serial',species,context,[2,2,2],'round_robin',0, &
          restored,read_owners,smaller,ok)
        call require(ok,'MPI reads serial SDC checkpoint')
        call verify_sdc(variant,small_rank)
        call MPI_Comm_free(smaller,ierr)
      end if
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      if(variant==2) then
        saved=restored;saved_owners=read_owners
        if(rank==0) then
          inquire(file=trim(path),size=before_size)
          cell_position=before_size-len('END_AMR_TREE_CHECKPOINT')-1-28+1
          corrupt_source=0;corrupt_source(1)=ieee_value(0._dp,ieee_quiet_nan)
          write(source_record,'(es27.18e3)') corrupt_source(1)
          open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',action='write')
          write(unit,pos=cell_position) source_record(:27)//achar(10)
          close(unit)
        end if
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call read_mpi_amr_tree_checkpoint_3d(trim(path),species,context,[2,2,2],'weighted',2, &
          restored,read_owners,MPI_COMM_WORLD,ok)
        call require(.not.ok,'nonfinite late SDC payload rejected')
        call require(all(read_owners==saved_owners),'SDC read owner rollback')
        call require(allocated(restored%reaction_source%fields).eqv.allocated(saved%reaction_source%fields), &
          'SDC source allocation rollback')
        if(allocated(saved%reaction_source%fields)) then
          do node=1,2
            if(.not.allocated(saved%reaction_source%fields(node)%values)) cycle
            call require(all(restored%reaction_source%fields(node)%values==saved%reaction_source%fields(node)%values), &
              'SDC source payload rollback')
            call require(all(restored%heat_release%fields(node)%values==saved%heat_release%fields(node)%values), &
              'late heat payload rollback')
          end do
        end if
      end if
    end do
    if(rank==0) print *, 'PASS: SDC V3/V4, signed sources, serial interchange and changed-rank read'
  end subroutine
  subroutine verify_sdc(variant,reader_rank)
    integer, intent(in) :: variant,reader_rank
    integer :: peer,c
    call require(all(restored%source==metadata%source),'MPI source ledger exact')
    call require(allocated(restored%reaction_source%fields).eqv.(variant/=3),'MPI source presence')
    call require(allocated(restored%heat_release%fields).eqv.(variant>=2),'MPI heat presence')
    if(variant>=2) then
      do peer=1,2
        call require(allocated(restored%heat_release%fields(peer)%values).eqv.(read_owners(peer)==reader_rank), &
          'owner-only heat payload')
        if(read_owners(peer)/=reader_rank) cycle
        call require(all(restored%heat_release%fields(peer)%values==real((-1)**peer*peer,dp)*1.e9_dp), &
          'signed heat exact')
      end do
    end if
    if(variant==3) return
    do peer=1,2
      call require(allocated(restored%reaction_source%fields(peer)%values).eqv.(read_owners(peer)==reader_rank), &
        'owner-only SDC checkpoint payload')
      if(read_owners(peer)/=reader_rank) cycle
      do c=1,8
        call require(all(restored%reaction_source%fields(peer)%values(c,:,:,:)==real((-1)**c*c*peer,dp)), &
          'MPI signed reaction source exact')
      end do
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: status
    if (condition) return
    print *, 'FAIL rank ',rank,': '//label
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
end program
