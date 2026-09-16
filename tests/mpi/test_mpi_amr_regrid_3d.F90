program test_mpi_amr_regrid_3d
  use, intrinsic :: iso_fortran_env, only: int64
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
    gather_mpi_amr_sparse_hierarchy_3d,regrid_mpi_amr_sparse_3d,repartition_mpi_amr_sparse_3d, &
    balance_mpi_amr_sparse_3d,write_mpi_amr_sparse_csv_3d
  use reactive_csv_io_3d_mod, only: write_reactive_3d_csv
  use mpi_amr_checkpoint_3d_mod, only: write_mpi_amr_checkpoint_3d,read_mpi_amr_checkpoint_3d
  use amr_reactive_3d_checkpoint_mod, only: write_amr_reactive_3d_checkpoint,read_amr_reactive_3d_checkpoint, &
    read_amr_checkpoint_metadata
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use simulation_config_amr_reactive_3d_mod, only: amr_reactive_3d_config
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
  integer, allocatable :: counts(:)
  logical :: ok,changed

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
    call check_csv(.true.)
    allocate(counts(distribution%nranks))
    do mode=1,5
      counts=1
      if (modulo(mode,2)==1) then
        counts(1)=9-distribution%nranks
      else
        counts(distribution%nranks)=9-distribution%nranks
      end if
      call repartition_mpi_amr_sparse_3d(distribution,patch,hierarchy,counts,ok)
      call require(ok,'direct repartition')
      call require(all(distribution%coarse_counts==counts),'requested coarse ownership')
      call gather_mpi_amr_sparse_hierarchy_3d(distribution,patch,0,hierarchy,gc,gt,gf,gft,ok)
      call require(ok,'gather repartition')
      if (rank==0) call require(all(c==gc).and.all(t==gt).and.all(f==gf).and.all(ft==gft), &
        'repartition preserves all state and temperature values')
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
      call require(all(distribution%coarse_counts==counts),'remap preserves chosen coarse partition')
      call check_balance()
      call check_csv(.false.)
      call check_checkpoint()
    end do
    saved=hierarchy;saved_distribution=distribution;saved_patch=patch
    counts=distribution%coarse_counts
    if (rank==0) counts(1)=0
    call repartition_mpi_amr_sparse_3d(distribution,patch,hierarchy,counts,ok)
    call rollback()
    if (distribution%nranks<8) then
      counts=distribution%coarse_counts
      if (rank==0) then
        counts(1)=counts(1)-1;counts(2)=counts(2)+1
      end if
      call repartition_mpi_amr_sparse_3d(distribution,patch,hierarchy,counts,ok)
      call rollback()
    end if
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
    deallocate(f,ft,faces,counts)
  end do
  if (rank==0) print '(a)','test_mpi_amr_regrid_3d: PASS'
  call MPI_Finalize(ierr)
contains
  subroutine check_checkpoint()
    type(reactive_3d_config) :: flow
    type(amr_reactive_3d_config) :: amr
    character(len=128) :: path
    character(len=256) :: message
    real(dp) :: integrals(12),clock,saved_t,saved_cell(12)
    integer :: parent,loaded_steps,rows,unit,output_unit,io_status
    character(len=32768) :: line
    integer(int64) :: coarse_offset,fine_offset
    real(dp) :: loaded_time,loaded_reflux,loaded_integrals(12),ledger(12)
    real(dp), allocatable :: ref_c(:,:,:,:),ref_t(:,:,:),ref_f(:,:,:,:),ref_ft(:,:,:)
    type(mpi_amr_sparse_hierarchy_3d) :: loaded,before
    type(mpi_amr_sparse_distribution_3d) :: loaded_distribution
    type(amr_patch_3d) :: loaded_patch
    character :: first_byte
    flow%nx=8;flow%ny=8;flow%nz=8;flow%thermo_model='elementary'
    flow%boundary_condition='outflow';flow%chemistry_enabled=.false.;flow%transport_enabled=.false.
    flow%final_time=1.0_dp
    amr%coarse_i_lower=patch%coarse_i_lower;amr%coarse_i_upper=patch%coarse_i_upper
    amr%coarse_j_lower=patch%coarse_j_lower;amr%coarse_j_upper=patch%coarse_j_upper
    amr%coarse_k_lower=patch%coarse_k_lower;amr%coarse_k_upper=patch%coarse_k_upper
    amr%refinement_ratio=r
    integrals=0;clock=.1_dp
    write(path,'(a,i0,a,i0,a)') 'distributed_checkpoint_np',distribution%nranks,'_r',r,'.chk'
    call write_mpi_amr_checkpoint_3d(distribution,hierarchy,path,species,flow,amr,patch, &
      clock,1,integrals,0.0_dp,ok,message,boundary_change=integrals)
    call require(ok,'distributed checkpoint write: '//trim(message))
    if (rank==0) then
      call write_amr_reactive_3d_checkpoint(trim(path)//'.serial',species,flow,amr,patch,c,t,f,ft, &
        clock,1,integrals,0.0_dp,ok,message,boundary_change=integrals)
      call require(ok,'serial checkpoint write: '//trim(message))
      call same_file(path)
    end if
    loaded_distribution=distribution;loaded_patch=patch
    loaded_time=-1;loaded_steps=-1;loaded_reflux=-1;loaded_integrals=-1;ledger=-1
    call read_mpi_amr_checkpoint_3d(loaded_distribution,loaded,path,species,flow,amr,loaded_patch, &
      loaded_time,loaded_steps,loaded_integrals,loaded_reflux,ok,message,boundary_change=ledger)
    call require(ok,'distributed checkpoint read: '//trim(message))
    call require(loaded_time==clock.and.loaded_steps==1.and.loaded_reflux==0.and. &
      all(loaded_integrals==0).and.all(ledger==0),'distributed restart diagnostics')
    call gather_mpi_amr_sparse_hierarchy_3d(loaded_distribution,loaded_patch,0,loaded,gc,gt,gf,gft,ok)
    call require(ok,'gather checkpoint candidate')
    if (rank==0) then
      ref_c=c;ref_t=t;ref_f=f;ref_ft=ft
      call read_amr_reactive_3d_checkpoint(path,species,flow,amr,patch,ref_c,ref_t,ref_f,ref_ft, &
        loaded_time,loaded_steps,loaded_integrals,loaded_reflux,ok,message,boundary_change=ledger)
      call require(ok,'independent serial restart')
      call require(all(ref_c==gc).and.all(ref_t==gt).and.all(ref_f==gf).and.all(ref_ft==gft), &
        'distributed restart exactly matches serial state and recovered temperatures')
    end if
    if (mode==1) then
      before=loaded
      ! Corrupt a coarse cell owned by the final rank, keeping metadata and length intact.
      if (rank==0) then
        call read_amr_checkpoint_metadata(path,species,flow,amr,loaded_patch,coarse_offset,fine_offset,rows, &
          loaded_time,loaded_steps,loaded_integrals,loaded_reflux,ok,message,boundary_change=ledger)
        call require(ok,'metadata byte offsets')
        coarse_offset=coarse_offset+int(distribution%coarse_displacements(distribution%nranks),int64)*rows+1
        open(newunit=unit,file=path,status='old',access='stream',form='unformatted',action='readwrite')
        read(unit,pos=coarse_offset) first_byte
        write(unit,pos=coarse_offset) '!'
        close(unit)
      end if
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call read_mpi_amr_checkpoint_3d(loaded_distribution,loaded,path,species,flow,amr,loaded_patch, &
        loaded_time,loaded_steps,loaded_integrals,loaded_reflux,ok,message,boundary_change=ledger)
      call require(.not.ok,'corrupt non-root-owned cell rejects collectively')
      call require(all(loaded%coarse_state==before%coarse_state).and.all(loaded%fine_state==before%fine_state).and. &
        all(loaded%coarse_temperature==before%coarse_temperature).and. &
        all(loaded%fine_temperature==before%fine_temperature),'rejected restart preserves hierarchy')
      call require(loaded_time==clock.and.loaded_steps==1.and.loaded_reflux==0.and. &
        all(loaded_integrals==0).and.all(ledger==0),'rejected restart preserves diagnostics')
      call require(all(loaded_distribution%coarse_counts==distribution%coarse_counts).and. &
        all(loaded_distribution%fine_counts==distribution%fine_counts).and. &
        loaded_patch%coarse_i_lower==patch%coarse_i_lower.and.loaded_patch%coarse_i_upper==patch%coarse_i_upper, &
        'rejected restart preserves geometry and ownership')
      if (rank==0) then
        open(newunit=unit,file=path,status='old',access='stream',form='unformatted',action='write')
        write(unit,pos=coarse_offset) first_byte
        close(unit)
      end if
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      if (rank==0) then
        open(newunit=unit,file=path,status='old',action='read')
        open(newunit=output_unit,file=trim(path)//'.spaced',status='replace',action='write')
        do
          read(unit,'(a)',iostat=io_status) line
          if (io_status/=0) exit
          write(output_unit,'(a)') trim(line)//' '
        end do
        close(unit);close(output_unit)
      end if
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call read_mpi_amr_checkpoint_3d(loaded_distribution,loaded,trim(path)//'.spaced',species,flow,amr,loaded_patch, &
        loaded_time,loaded_steps,loaded_integrals,loaded_reflux,ok,message,boundary_change=ledger)
      call require(ok,'streaming fallback accepts reformatted checkpoint: '//trim(message))
      call require(all(loaded%coarse_state==before%coarse_state).and.all(loaded%fine_state==before%fine_state).and. &
        all(loaded%coarse_temperature==before%coarse_temperature).and. &
        all(loaded%fine_temperature==before%fine_temperature),'reformatted restart preserves exact state')
    end if
    if (mode/=1) return
    if (distribution%fine_count>0) then
      saved_t=hierarchy%fine_temperature(1,1,1)
      hierarchy%fine_temperature(1,1,1)=saved_t*1.1_dp
    end if
    call write_mpi_amr_checkpoint_3d(distribution,hierarchy,path,species,flow,amr,patch, &
      clock,1,integrals,0.0_dp,ok,message,boundary_change=integrals)
    call require(.not.ok,'inconsistent distributed checkpoint temperature rejected')
    if (distribution%fine_count>0) hierarchy%fine_temperature(1,1,1)=saved_t
    if (distribution%fine_count>0) then
      parent=distribution%coarse_local_index(max(distribution%coarse_first,patch%coarse_i_lower))
      saved_cell=hierarchy%coarse_state(:,parent,patch%coarse_j_lower,patch%coarse_k_lower)
      hierarchy%coarse_state(:,parent,patch%coarse_j_lower,patch%coarse_k_lower)=1.01_dp*saved_cell
    end if
    call write_mpi_amr_checkpoint_3d(distribution,hierarchy,path,species,flow,amr,patch, &
      clock,1,integrals,0.0_dp,ok,message,boundary_change=integrals)
    call require(.not.ok.and.index(message,'not synchronized')>0,'unsynchronized distributed checkpoint rejected')
    if (distribution%fine_count>0) &
      hierarchy%coarse_state(:,parent,patch%coarse_j_lower,patch%coarse_k_lower)=saved_cell
    if (rank==0) clock=.2_dp
    call write_mpi_amr_checkpoint_3d(distribution,hierarchy,path,species,flow,amr,patch, &
      clock,1,integrals,0.0_dp,ok,message,boundary_change=integrals)
    call require(.not.ok,'inconsistent checkpoint metadata rejected')
    if (rank==0) call same_file(path)
  end subroutine check_checkpoint
  subroutine check_csv(check_rejection)
    logical, intent(in) :: check_rejection
    real(dp) :: xyz(8),fx(patch%fine_nx()),fy(patch%fine_ny()),fz(patch%fine_nz()),saved_density
    character(len=128) :: prefix,cp,fp,path
    character(len=256) :: message
    integer :: n
    xyz=[(real(n,dp),n=1,8)]
    fx=[(real(n,dp)/r,n=1,size(fx))];fy=[(real(n,dp)/r,n=1,size(fy))];fz=[(real(n,dp)/r,n=1,size(fz))]
    write(prefix,'(a,i0,a,i0)') 'distributed_csv_np',distribution%nranks,'_r',r
    cp=trim(prefix)//'_coarse.csv';fp=trim(prefix)//'_fine.csv'
    call write_mpi_amr_sparse_csv_3d(distribution,species,patch,hierarchy,cp,fp,xyz,xyz,xyz,fx,fy,fz,0.1_dp,ok)
    call require(ok,'distributed CSV including empty fine ranks')
    if (rank==0) then
      call write_reactive_3d_csv(trim(cp)//'.serial',species,xyz,xyz,xyz,c,t,8,8,8,0.1_dp,ok,message)
      call require(ok,'serial coarse CSV')
      call write_reactive_3d_csv(trim(fp)//'.serial',species,fx,fy,fz,f,ft, &
        size(fx),size(fy),size(fz),0.1_dp,ok,message)
      call require(ok,'serial fine CSV')
      call same_file(cp);call same_file(fp)
    end if
    if (.not.check_rejection) return
    ! Invalid fine data must be detected before either existing output is opened.
    if (distribution%fine_count>0) then
      saved_density=hierarchy%fine_state(1,1,1,1)
      hierarchy%fine_state(1,1,1,1)=-1.0_dp
    end if
    call write_mpi_amr_sparse_csv_3d(distribution,species,patch,hierarchy,cp,fp,xyz,xyz,xyz,fx,fy,fz,0.2_dp,ok)
    call require(.not.ok,'invalid fine CSV state rejected')
    if (distribution%fine_count>0) hierarchy%fine_state(1,1,1,1)=saved_density
    path=fp
    if (rank==0) path=trim(fp)//'.different'
    call write_mpi_amr_sparse_csv_3d(distribution,species,patch,hierarchy,cp,path,xyz,xyz,xyz,fx,fy,fz,0.2_dp,ok)
    call require(.not.ok,'rank-inconsistent fine path rejected')
    if (rank==0) then
      call same_file(cp);call same_file(fp)
    end if
  end subroutine check_csv
  subroutine same_file(path)
    character(len=*), intent(in) :: path
    character(len=:), allocatable :: actual,expected
    integer :: unit,bytes,reference_bytes
    inquire(file=trim(path),size=bytes)
    inquire(file=trim(path)//'.serial',size=reference_bytes)
    call require(bytes==reference_bytes,'CSV size and overwrite truncation match serial')
    allocate(character(len=bytes)::actual,expected)
    open(newunit=unit,file=trim(path),access='stream',form='unformatted',status='old')
    read(unit) actual
    close(unit)
    open(newunit=unit,file=trim(path)//'.serial',access='stream',form='unformatted',status='old')
    read(unit) expected
    close(unit)
    call require(actual==expected,'CSV bytes match independent serial writer')
  end subroutine same_file
  subroutine check_balance()
    real(dp) :: cost(8),best,peak
    integer :: peer,first
    cost=real(patch%coarse_ny*patch%coarse_nz,dp)
    cost(patch%coarse_i_lower:patch%coarse_i_upper)= &
      cost(patch%coarse_i_lower:patch%coarse_i_upper)+real(r*patch%fine_ny()*patch%fine_nz(),dp)
    best=huge(1.0_dp)
    call enumerate(1,distribution%nranks,0.0_dp,cost,best)
    call balance_mpi_amr_sparse_3d(distribution,patch,hierarchy,changed,ok)
    call require(ok,'automatic balance')
    peak=0
    do peer=1,distribution%nranks
      first=distribution%coarse_displacements(peer)+1
      peak=max(peak,sum(cost(first:first+distribution%coarse_counts(peer)-1)))
    end do
    call require(peak==best,'balance reaches exhaustive minimum maximum cell count')
    call gather_mpi_amr_sparse_hierarchy_3d(distribution,patch,0,hierarchy,gc,gt,gf,gft,ok)
    call require(ok,'gather balanced hierarchy')
    if (rank==0) call require(all(c==gc).and.all(t==gt).and.all(f==gf).and.all(ft==gft), &
      'automatic balance preserves serial state')
    call balance_mpi_amr_sparse_3d(distribution,patch,hierarchy,changed,ok)
    call require(ok.and..not.changed,'repeated balance does not transfer again')
  end subroutine check_balance
  recursive subroutine enumerate(first,remaining,previous,cost,best)
    integer, intent(in) :: first,remaining
    real(dp), intent(in) :: previous,cost(:)
    real(dp), intent(inout) :: best
    integer :: last
    if (remaining==1) then
      best=min(best,max(previous,sum(cost(first:))))
      return
    end if
    do last=first,size(cost)-remaining+1
      call enumerate(last+1,remaining-1,max(previous,sum(cost(first:last))),cost,best)
    end do
  end subroutine enumerate
  subroutine rollback()
    call require(.not.ok,'invalid remap rejected')
    call require(all(hierarchy%coarse_state==saved%coarse_state).and. &
      all(hierarchy%coarse_temperature==saved%coarse_temperature).and. &
      all(hierarchy%fine_state==saved%fine_state).and.all(hierarchy%fine_temperature==saved%fine_temperature), &
      'failure preserves hierarchy')
    call require(all(distribution%coarse_counts==saved_distribution%coarse_counts).and. &
      all(distribution%coarse_displacements==saved_distribution%coarse_displacements).and. &
      all(distribution%fine_counts==saved_distribution%fine_counts).and. &
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
