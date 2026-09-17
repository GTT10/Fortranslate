module spray_mpi_eb_particle_checkpoint_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64,iostat_end
  use, intrinsic :: iso_c_binding, only: c_int,c_char,c_null_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use spray_parcel_mod, only: spray_parcel,valid_spray_parcel
  use spray_mixture_mod, only: valid_component_inventory
  use spray_mpi_exchange_mod, only: exchange_spray_parcel_bins
  use spray_mpi_amr_tree_mod, only: migrate_eb_tree_spray_parcels
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d
  implicit none
  private
  public :: mpi_eb_particle_checkpoint
  interface
    function link_file(old,new) bind(C,name='link') result(status)
      import c_char,c_int
      character(c_char), intent(in) :: old(*),new(*)
      integer(c_int) :: status
    end function
  end interface
contains
  ! Particle payload only. Caller binds context to the gas generation and model.
  ! A closed .pending file is linked to an immutable final path. Reads are atomic.
  subroutine mpi_eb_particle_checkpoint(comm,writing,path,context,tree,owners,phi,lower,spacing,capacity, &
      parcels,components,time,steps,event,last_id,injected,ok)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: writing
    character(len=*), intent(in) :: path,context
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),capacity
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),spacing(3)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: components(:,:)
    real(dp), intent(inout) :: time,injected(:)
    integer, intent(inout) :: steps,event
    integer(int64), intent(inout) :: last_id
    logical, intent(out) :: ok
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:),ledger(:)
    integer, allocatable :: counts(:),offsets(:)
    integer :: rank,nranks,ierr,nc,nv,unit,ios,closed,np,p,info(2),dims(2),flag,lo,hi
    integer(int64) :: mark
    real(dp) :: at_time
    logical :: valid,agreed
    character(len=:), allocatable :: line
    ok=.false.
    allocate(character(len=131072)::line)
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=allocated(parcels).and.allocated(components)
    valid=valid.and.len_trim(path)>0.and.len(context)>0.and.len(context)<=len(line)
    valid=valid.and.index(context,achar(10))==0.and.index(context,achar(13))==0
    valid=valid.and.capacity>=0.and.capacity<=1000000.and.size(injected)>=6.and.size(injected)<=133
    call agree();if(.not.agreed) return
    nc=size(components,1);nv=size(injected)
    valid=nc>=1.and.nc<=16.and.size(components,2)==size(parcels)
    call agree();if(.not.agreed) return
    call same_string(path);if(.not.agreed) return
    call same_string(context);if(.not.agreed) return
    flag=merge(1,0,writing)
    call MPI_Allreduce(flag,lo,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(flag,hi,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.lo/=hi) return
    dims=[nc,nv]
    call MPI_Bcast(dims,2,MPI_INTEGER,0,comm,ierr)
    valid=ierr==MPI_SUCCESS.and.all(dims==[nc,nv])
    call agree();if(.not.agreed) return
    allocate(ledger(nv));at_time=time;info=[steps,event];mark=last_id;ledger=injected
    if(writing) then
      call broadcast_metadata()
      valid=valid.and.at_time==time.and.all(info==[steps,event]).and.mark==last_id.and.all(ledger==injected)
      valid=valid.and.ieee_is_finite(time).and.all(ieee_is_finite(ledger))
      valid=valid.and.time>=0.and.steps>=0.and.event>=1.and.last_id>=0
      call agree();if(.not.agreed) return
      drops=parcels;masses=components
      call migrate_eb_tree_spray_parcels(comm,tree,owners,phi,drops,lower,spacing,capacity,valid,masses)
      if(.not.valid) return
      allocate(counts(nranks),offsets(nranks))
      counts=0;counts(1)=size(drops);offsets=counts(1);offsets(1)=0
      call exchange_spray_parcel_bins(comm,drops,counts,offsets,capacity,valid,masses)
      if(.not.valid) return
      if(rank==0) then
        valid=all(drops%id<=mark)
        if(valid) then
          open(newunit=unit,file=path//'.pending',status='new',action='write',iostat=ios)
          valid=ios==0
          if(valid) then
            write(unit,'(a)',iostat=ios) 'PELEF_EB_PARTICLES_V1',context
            if(ios==0) write(unit,*,iostat=ios) nc,nv,size(drops),mark,info
            if(ios==0) write(unit,'(*(es27.18e3,1x))',iostat=ios) at_time,ledger
            do p=1,size(drops)
              if(ios/=0) exit
              write(unit,'(3(i0,1x),*(es27.18e3,1x))',iostat=ios) drops(p)%id,drops(p)%parent_id,drops(p)%material_id, &
                drops(p)%position,drops(p)%velocity,drops(p)%mass,drops(p)%temperature,drops(p)%multiplicity, &
                drops(p)%deformation,drops(p)%deformation_rate,drops(p)%rt_time,drops(p)%khrt_initial_count, &
                drops(p)%khrt_shed_volume,masses(:,p)
            end do
            if(ios==0) write(unit,'(a)',iostat=ios) 'END_EB_PARTICLES'
            close(unit,iostat=closed);valid=ios==0.and.closed==0
            if(valid) valid=link_file(path//'.pending'//c_null_char,path//c_null_char)==0
          end if
        end if
      end if
      call agree();ok=agreed;return
    end if
    allocate(drops(0),masses(nc,0))
    valid=.true.
    if(rank==0) then
      open(newunit=unit,file=path,status='old',action='read',iostat=ios)
      valid=ios==0
      if(valid) then
        call read_payload()
        close(unit,iostat=closed);valid=valid.and.closed==0
      end if
    end if
    call agree();if(.not.agreed) return
    call broadcast_metadata()
    call agree();if(.not.agreed) return
    call migrate_eb_tree_spray_parcels(comm,tree,owners,phi,drops,lower,spacing,capacity,valid,masses)
    if(.not.valid) return
    call move_alloc(drops,parcels);call move_alloc(masses,components)
    time=at_time;steps=info(1);event=info(2);last_id=mark;injected=ledger;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
    subroutine same_string(value)
      character(len=*), intent(in) :: value
      character(len=:), allocatable :: reference
      integer :: n
      n=len(value)
      call MPI_Bcast(n,1,MPI_INTEGER,0,comm,ierr)
      if(ierr/=MPI_SUCCESS) then
        agreed=.false.;return
      end if
      allocate(character(len=n)::reference);reference(:)=value
      call MPI_Bcast(reference,n,MPI_CHARACTER,0,comm,ierr)
      valid=ierr==MPI_SUCCESS.and.len(value)==n.and.reference==value
      call agree()
    end subroutine
    subroutine broadcast_metadata()
      valid=.false.
      call MPI_Bcast(at_time,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if(ierr/=MPI_SUCCESS) return
      call MPI_Bcast(info,2,MPI_INTEGER,0,comm,ierr)
      if(ierr/=MPI_SUCCESS) return
      call MPI_Bcast(mark,1,MPI_INTEGER8,0,comm,ierr)
      if(ierr/=MPI_SUCCESS) return
      call MPI_Bcast(ledger,nv,MPI_DOUBLE_PRECISION,0,comm,ierr)
      valid=ierr==MPI_SUCCESS
    end subroutine
    subroutine read_payload()
      valid=.false.
      read(unit,'(a)',iostat=ios) line
      if(ios/=0.or.line/='PELEF_EB_PARTICLES_V1') return
      read(unit,'(a)',iostat=ios) line
      if(ios/=0.or.line/=context) return
      read(unit,*,iostat=ios) dims,np,mark,info
      if(ios/=0) return
      if(any(dims/=[nc,nv]).or.np<0.or.np>capacity.or.mark<0.or.info(1)<0.or.info(2)<1) return
      read(unit,*,iostat=ios) at_time,ledger
      if(ios/=0) return
      if(.not.ieee_is_finite(at_time).or..not.all(ieee_is_finite(ledger))) return
      if(at_time<0) return
      deallocate(drops,masses);allocate(drops(np),masses(nc,np),stat=ios)
      if(ios/=0) return
      do p=1,np
        read(unit,*,iostat=ios) drops(p)%id,drops(p)%parent_id,drops(p)%material_id, &
          drops(p)%position,drops(p)%velocity,drops(p)%mass,drops(p)%temperature,drops(p)%multiplicity, &
          drops(p)%deformation,drops(p)%deformation_rate,drops(p)%rt_time,drops(p)%khrt_initial_count, &
          drops(p)%khrt_shed_volume,masses(:,p)
        if(ios/=0) return
        if(.not.valid_spray_parcel(drops(p)).or..not.valid_component_inventory(drops(p),masses(:,p))) return
        if(drops(p)%id>mark) return
        if(p>1) then
          if(drops(p)%id<=drops(p-1)%id) return
        end if
      end do
      read(unit,'(a)',iostat=ios) line
      if(ios/=0.or.line/='END_EB_PARTICLES') return
      read(unit,'(a)',iostat=ios) line
      valid=ios==iostat_end
    end subroutine
  end subroutine
end module
