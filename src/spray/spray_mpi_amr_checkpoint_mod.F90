module spray_mpi_amr_checkpoint_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: iso_c_binding, only: c_int,c_char,c_null_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use spray_parcel_mod, only: spray_parcel,valid_spray_parcel
  use spray_mixture_mod, only: valid_component_inventory
  use spray_mpi_exchange_mod, only: exchange_spray_parcel_bins
  use spray_wall_impact_mod, only: spray_wall_film,valid_wall_component_inventory
  use spray_mpi_amr_film_mod, only: migrate_mpi_amr_wall_films
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use amr_patch_tree_checkpoint_3d_mod, only: amr_tree_checkpoint_3d
  use mpi_amr_patch_tree_checkpoint_3d_mod, only: read_mpi_amr_tree_checkpoint_3d,write_mpi_amr_tree_checkpoint_3d
  implicit none
  private
  public :: mpi_amr_spray_checkpoint
  interface
    function link_file(old,new) bind(C,name='link') result(status)
      import c_char,c_int
      character(kind=c_char), intent(in) :: old(*),new(*)
      integer(c_int) :: status
    end function
  end interface
contains
  ! Immutable prefix: .pending reserves it, .gas/.particles are payloads, and
  ! the prefix itself is published last with POSIX link. Reads publish atomically.
  ! Shared path/context/model tables are required; at most one million total parcels.
  subroutine mpi_amr_spray_checkpoint(comm,writing,path,context,species,tree,owners,state,temperature, &
      parcels,component_mass,time,steps,event,regridded,last_id,injected,lower,spacing,maximum_parcels,ok,periodic,boundary, &
      wall_films,film_mass,maximum_films,source,reaction_source,heat_release,heat_available,frame,coordinate_energy, &
      internal_energy,energy_correction,internal_energy_pressure_work_change)
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    real(dp), optional, intent(inout) :: energy_correction,internal_energy_pressure_work_change
    real(dp) :: work_ledger,reference_work
    real(dp) :: held_ledger,reference_held
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(inout) :: coordinate_energy
    real(dp) :: coordinate_ledger,reference_coordinate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: writing
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(inout) :: tree
    integer, allocatable, intent(inout) :: owners(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(inout) :: time,injected(:)
    integer, intent(inout) :: steps,event
    logical, intent(inout) :: regridded
    integer(int64), intent(inout) :: last_id
    real(dp), intent(in) :: lower(3),spacing(3)
    integer, intent(in) :: maximum_parcels
    logical, intent(out) :: ok
    logical, optional, intent(in) :: periodic(3)
    real(dp), optional, intent(inout) :: boundary(:),source(:)
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source,heat_release
    logical, optional, intent(inout) :: heat_available
    type(amr_tree_state_3d) :: rates,heat
    logical :: with_heat
    real(dp), allocatable :: source_ledger(:),reference_source(:)
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    integer, optional, intent(in) :: maximum_films
    type(spray_wall_film), allocatable :: films(:)
    real(dp), allocatable :: film_components(:,:)
    character(len=32) :: payload_header
    integer :: nf
    logical :: wraps(3)
    real(dp), allocatable :: boundary_ledger(:),reference_boundary(:)
    character(len=3) :: boundary_key
    type(amr_tree_checkpoint_3d) :: snapshot
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:),ledger(:)
    integer, allocatable :: assigned(:),counts(:),offsets(:)
    integer :: rank,nranks,ierr,unit,ios,p,nc,np,nvar,info(3),date(8),low(13),high(13),header(13),node,n(3)
    integer(int64) :: mark,tick
    real(dp) :: at_time
    logical :: valid,agreed,exists
    character(len=256) :: token,line
    character(len=:), allocatable :: bound_context
    ok=.false.
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=allocated(owners).and.allocated(parcels).and.allocated(component_mass)
    valid=valid.and.tree%is_valid().and.len_trim(path)>0.and.len(context)>0
    valid=valid.and.(present(wall_films).eqv.present(film_mass))
    valid=valid.and.(present(wall_films).eqv.present(maximum_films))
    if (present(wall_films)) valid=valid.and.allocated(wall_films)
    if (present(film_mass)) valid=valid.and.allocated(film_mass)
    if (present(maximum_films)) valid=valid.and.maximum_films>=0.and.maximum_films<=1000000
    call agree();if (.not.agreed) return
    nc=size(component_mass,1);nvar=5+size(species)
    valid=nc>=1.and.nc<=16.and.all(shape(component_mass)==[nc,size(parcels)]).and.size(injected)==nvar
    if (present(boundary)) valid=valid.and.size(boundary)==nvar
    if (present(source)) valid=valid.and.size(source)==nvar
    call agree();if (.not.agreed) return
    wraps=.true.
    if (present(periodic)) wraps=periodic
    write(boundary_key,'(3l1)') wraps
    call same_string(boundary_key,valid);if (.not.valid) return
    with_heat=present(heat_release)
    if(present(heat_available).and.writing) with_heat=with_heat.and.heat_available
    header=[nc,nvar,merge(1,0,writing),merge(1,0,present(wall_films)),merge(1,0,present(source)), &
      merge(1,0,present(reaction_source)),merge(1,0,present(heat_release)),merge(1,0,with_heat), &
      merge(1,0,present(frame)),merge(1,0,present(coordinate_energy)), &
      merge(1,0,present(internal_energy)),merge(1,0,present(energy_correction)), &
      merge(1,0,present(internal_energy_pressure_work_change))]
    call MPI_Allreduce(header,low,13,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(header,high,13,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.any(low/=high)) return
    valid=present(frame).eqv.present(coordinate_energy)
    valid=valid.and.(present(internal_energy).eqv.present(energy_correction))
    valid=valid.and.(.not.present(internal_energy_pressure_work_change).or.present(internal_energy))
    call agree();if(.not.agreed) return
    held_ledger=0
    if(present(energy_correction)) held_ledger=energy_correction
    valid=ieee_is_finite(held_ledger);call agree();if(.not.agreed) return
    reference_held=held_ledger
    call MPI_Bcast(reference_held,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=held_ledger==reference_held;call agree();if(.not.agreed) return
    work_ledger=0
    if(present(internal_energy_pressure_work_change)) work_ledger=internal_energy_pressure_work_change
    valid=ieee_is_finite(work_ledger);call agree();if(.not.agreed) return
    reference_work=work_ledger
    call MPI_Bcast(reference_work,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=work_ledger==reference_work;call agree();if(.not.agreed) return
    coordinate_ledger=0
    if(present(coordinate_energy)) coordinate_ledger=coordinate_energy
    payload_header='PELEF_AMR_PARTICLES_V1'
    if (present(wall_films)) payload_header='PELEF_AMR_PARTICLES_V2'
    ! Verify shared strings before root performs any I/O.
    call same_string(path,valid);if (.not.valid) return
    call same_string(context,valid);if (.not.valid) return
    allocate(ledger(nvar),boundary_ledger(nvar),reference_boundary(nvar));token=''
    allocate(source_ledger(nvar),reference_source(nvar));source_ledger=0
    if (present(source)) source_ledger=source
    boundary_ledger=0
    if (present(boundary)) boundary_ledger=boundary
    if (writing) then
      valid=ieee_is_finite(time).and.all(ieee_is_finite(injected)).and.all(ieee_is_finite(boundary_ledger))
      valid=valid.and.all(ieee_is_finite(source_ledger)).and.ieee_is_finite(coordinate_ledger)
      valid=valid.and.time>=0.and.steps>=0.and.event>=0.and.last_id>=0
      call agree();if (.not.agreed) return
      info=[steps,event,merge(1,0,regridded)];at_time=time;mark=last_id;ledger=injected
      call MPI_Bcast(info,3,MPI_INTEGER,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Bcast(at_time,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Bcast(mark,1,MPI_INTEGER8,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Bcast(ledger,nvar,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      valid=all(info==[steps,event,merge(1,0,regridded)]).and.at_time==time.and.mark==last_id.and.all(ledger==injected)
      call agree();if (.not.agreed) return
      reference_boundary=boundary_ledger
      call MPI_Bcast(reference_boundary,nvar,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      valid=all(reference_boundary==boundary_ledger)
      call agree();if (.not.agreed) return
      reference_source=source_ledger
      call MPI_Bcast(reference_source,nvar,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      valid=all(reference_source==source_ledger)
      call agree();if (.not.agreed) return
      reference_coordinate=coordinate_ledger
      call MPI_Bcast(reference_coordinate,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if(ierr/=MPI_SUCCESS) return
      valid=reference_coordinate==coordinate_ledger
      call agree();if(.not.agreed) return
      drops=parcels;masses=component_mass
      allocate(counts(nranks),offsets(nranks));counts=0;counts(1)=size(drops);offsets=counts(1);offsets(1)=0
      call exchange_spray_parcel_bins(comm,drops,counts,offsets,1000000,valid,masses)
      if (.not.valid) return
      if (present(wall_films)) then
        films=wall_films;film_components=film_mass
        valid=size(film_components,1)==nc
        call agree();if (.not.agreed) return
        call migrate_mpi_amr_wall_films(comm,tree,[(0,p=1,size(tree%nodes))],films,film_components, &
          lower,spacing,1000000,valid)
        if (.not.valid) return
      end if
      valid=.true.
      if (rank==0) then
        do p=1,size(drops)
          if (drops(p)%id>mark) valid=.false.
        end do
        inquire(file=path,exist=exists);valid=valid.and..not.exists
        inquire(file=path//'.gas',exist=exists);valid=valid.and..not.exists
        inquire(file=path//'.particles',exist=exists);valid=valid.and..not.exists
        if (valid) then
          open(newunit=unit,file=path//'.pending',status='new',action='write',iostat=ios)
          valid=ios==0
          if (valid) then
            call date_and_time(values=date);call system_clock(tick)
            write(token,'(8(i0,":"),i0)') date,tick
            write(unit,'(a)',iostat=ios) 'PELEF_AMR_SPRAY_V1',trim(token),'END_AMR_SPRAY_MANIFEST'
            valid=ios==0
            close(unit,iostat=ios);valid=valid.and.ios==0
          end if
        end if
      end if
      call agree();if (.not.agreed) return
      call MPI_Bcast(token,len(token),MPI_CHARACTER,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      bound_context=context//'SPRAY_GENERATION '//trim(token)//achar(10)
      if (.not.all(wraps)) bound_context=bound_context//'SPRAY_PERIODIC '//boundary_key//achar(10)
      if (present(wall_films)) bound_context=bound_context//'SPRAY_FILMS_V1'//achar(10)
      if (present(source)) bound_context=bound_context//'SPRAY_SOURCE_V1'//achar(10)
      if(present(frame)) bound_context=bound_context//'SPRAY_COORDINATE_ENERGY_V1'//achar(10)
    if(present(internal_energy)) bound_context=bound_context//'SPRAY_ACTUAL_UEINT_V1'//achar(10)
    if(present(internal_energy_pressure_work_change)) &
      bound_context=bound_context//'SPRAY_UEINT_PRESSURE_WORK_V1'//achar(10)
      snapshot%tree=tree;snapshot%time=time;snapshot%steps=steps;snapshot%regrids=info(3)
      snapshot%initial=injected;snapshot%boundary=boundary_ledger
      if(present(reaction_source).or.present(heat_release)) snapshot%source=source_ledger
      if(present(reaction_source)) snapshot%reaction_source=reaction_source
      if(with_heat) snapshot%heat_release=heat_release
      if(present(internal_energy)) then
        snapshot%internal_energy=internal_energy;snapshot%energy_correction=energy_correction
      endif
      if(present(internal_energy_pressure_work_change)) snapshot%pressure_work_change=work_ledger
      call write_mpi_amr_tree_checkpoint_3d(path//'.gas',species,bound_context,snapshot,owners,state,temperature, &
        comm,valid,indexed=.true.,frame=frame,domain_lower=lower,root_spacing=spacing)
      if (.not.valid) return
      if (rank==0) then
        open(newunit=unit,file=path//'.particles',status='new',action='write',iostat=ios)
        valid=ios==0
        if (valid) then
          write(unit,'(a)',iostat=ios) trim(payload_header),trim(token)
          if (ios==0) write(unit,*,iostat=ios) nc,size(drops),mark,info
          if (ios==0) write(unit,'(*(es27.18e3,1x))',iostat=ios) time,ledger
          do p=1,size(drops)
            if (ios/=0) exit
            write(unit,'(3(i0,1x),*(es27.18e3,1x))',iostat=ios) drops(p)%id,drops(p)%parent_id,drops(p)%material_id, &
              drops(p)%position,drops(p)%velocity,drops(p)%mass,drops(p)%temperature,drops(p)%multiplicity, &
              drops(p)%deformation,drops(p)%deformation_rate,drops(p)%rt_time,drops(p)%khrt_initial_count, &
              drops(p)%khrt_shed_volume,masses(:,p)
          end do
          if (present(wall_films)) then
            if (ios==0) write(unit,*,iostat=ios) size(films)
            do p=1,size(films)
              if (ios/=0) exit
              write(unit,'(i0,1x,*(es27.18e3,1x))',iostat=ios) films(p)%material_id,films(p)%position, &
                films(p)%normal,films(p)%mass,films(p)%temperature,films(p)%diameter,films(p)%height,film_components(:,p)
            end do
          end if
          if (present(source)) then
            if (ios==0) write(unit,'(a)',iostat=ios) 'SPRAY_SOURCE_V1'
            if (ios==0) write(unit,'(*(es27.18e3,1x))',iostat=ios) source_ledger
          end if
          if(present(frame)) then
            if(ios==0) write(unit,'(a)',iostat=ios) 'SPRAY_COORDINATE_ENERGY_V1'
            if(ios==0) write(unit,'(es27.18e3)',iostat=ios) coordinate_ledger
          endif
          if (ios==0) write(unit,'(a)',iostat=ios) 'END_AMR_PARTICLES'
          valid=ios==0;close(unit,iostat=ios);valid=valid.and.ios==0
        end if
        if (valid) valid=link_file(path//'.pending'//c_null_char,path//c_null_char)==0
      end if
      call agree();ok=agreed;return
    end if
    valid=.true.
    if (rank==0) then
      open(newunit=unit,file=path,status='old',action='read',iostat=ios)
      valid=ios==0
      if (valid) then
        read(unit,'(a)',iostat=ios) line
        valid=ios==0.and.line=='PELEF_AMR_SPRAY_V1'
        if (valid) read(unit,'(a)',iostat=ios) token
        valid=valid.and.ios==0.and.len_trim(token)>0
        if (valid) read(unit,'(a)',iostat=ios) line
        valid=valid.and.ios==0.and.line=='END_AMR_SPRAY_MANIFEST'
        if (valid) read(unit,'(a)',iostat=ios) line
        valid=valid.and.ios<0
        close(unit)
      end if
    end if
    call agree();if (.not.agreed) return
    call MPI_Bcast(token,len(token),MPI_CHARACTER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    bound_context=context//'SPRAY_GENERATION '//trim(token)//achar(10)
    if (.not.all(wraps)) bound_context=bound_context//'SPRAY_PERIODIC '//boundary_key//achar(10)
    if (present(wall_films)) bound_context=bound_context//'SPRAY_FILMS_V1'//achar(10)
    if (present(source)) bound_context=bound_context//'SPRAY_SOURCE_V1'//achar(10)
    if(present(frame)) bound_context=bound_context//'SPRAY_COORDINATE_ENERGY_V1'//achar(10)
    if(present(internal_energy)) bound_context=bound_context//'SPRAY_ACTUAL_UEINT_V1'//achar(10)
    if(present(internal_energy_pressure_work_change)) &
      bound_context=bound_context//'SPRAY_UEINT_PRESSURE_WORK_V1'//achar(10)
    call read_mpi_amr_tree_checkpoint_3d(path//'.gas',species,bound_context,tree%nodes(1)%shape, &
      'round_robin',0,snapshot,assigned,comm,valid,frame=frame,domain_lower=lower,root_spacing=spacing, &
      allow_internal_energy=present(internal_energy))
    if (.not.valid) return
    valid=present(internal_energy).eqv.allocated(snapshot%internal_energy%fields)
    valid=valid.and.(present(energy_correction).eqv.allocated(snapshot%energy_correction))
    valid=valid.and.(present(internal_energy_pressure_work_change).eqv.allocated(snapshot%pressure_work_change))
    call agree();if(.not.agreed) return
    valid=present(boundary).or.all(snapshot%boundary==0)
    call agree();if (.not.agreed) return
    allocate(drops(0),masses(nc,0));info=0;mark=0
    if (present(wall_films)) allocate(films(0),film_components(nc,0))
    if (rank==0) then
      open(newunit=unit,file=path//'.particles',status='old',action='read',iostat=ios)
      valid=ios==0
      if (valid) then
        call read_particles()
        close(unit)
      end if
    end if
    call agree();if (.not.agreed) return
    call MPI_Bcast(info,3,MPI_INTEGER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Bcast(mark,1,MPI_INTEGER8,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Bcast(ledger,nvar,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Bcast(source_ledger,nvar,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Bcast(coordinate_ledger,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call migrate_amr_tree_spray_parcels(comm,snapshot%tree,assigned,drops,lower,spacing,wraps, &
      maximum_parcels,valid,masses)
    if (.not.valid) return
    if (present(wall_films)) then
      call migrate_mpi_amr_wall_films(comm,snapshot%tree,assigned,films,film_components,lower,spacing,maximum_films,valid)
      if (.not.valid) return
    end if
    ! Prepare legacy work fields on the restored lattice before publishing state.
    if(present(reaction_source)) then
      if(allocated(snapshot%reaction_source%fields)) then
        rates=snapshot%reaction_source
      else
        allocate(rates%fields(size(assigned)))
        do node=1,size(assigned)
          if(assigned(node)/=rank) cycle
          n=snapshot%tree%nodes(node)%shape
          allocate(rates%fields(node)%values(nvar-4,n(1),n(2),n(3)),source=0._dp)
        end do
      end if
    end if
    if(present(heat_release)) then
      if(allocated(snapshot%heat_release%fields)) then
        heat=snapshot%heat_release
      else
        heat=snapshot%temperature
        do node=1,size(assigned)
          if(assigned(node)/=rank) cycle
          heat%fields(node)%values=0
        end do
      end if
    end if
    if(present(reaction_source)) call move_alloc(rates%fields,reaction_source%fields)
    if(present(heat_release)) call move_alloc(heat%fields,heat_release%fields)
    if(present(heat_available)) heat_available=allocated(snapshot%heat_release%fields)
    if(present(internal_energy)) then
      call move_alloc(snapshot%internal_energy%fields,internal_energy%fields)
      energy_correction=snapshot%energy_correction
    endif
    if(present(internal_energy_pressure_work_change)) internal_energy_pressure_work_change=snapshot%pressure_work_change
    tree=snapshot%tree;call move_alloc(assigned,owners)
    call move_alloc(snapshot%state%fields,state%fields);call move_alloc(snapshot%temperature%fields,temperature%fields)
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    time=snapshot%time;steps=info(1);event=info(2);regridded=info(3)==1;last_id=mark;injected=ledger
    if (present(wall_films)) then
      call move_alloc(films,wall_films);call move_alloc(film_components,film_mass)
    end if
    if (present(boundary)) boundary=snapshot%boundary
    if (present(source)) source=source_ledger
    if(present(coordinate_energy)) coordinate_energy=coordinate_ledger
    ok=.true.
  contains
    subroutine read_particles()
      integer :: file_nc
      valid=.false.
      read(unit,'(a)',iostat=ios) line
      if (ios/=0.or.line/=payload_header) return
      read(unit,'(a)',iostat=ios) line
      if (ios/=0.or.line/=token) return
      read(unit,*,iostat=ios) file_nc,np,mark,info
      if (ios/=0) return
      if (file_nc/=nc.or.np<0.or.np>1000000.or.mark<0.or.any(info<0).or.info(3)>1) return
      read(unit,*,iostat=ios) at_time,ledger
      if (ios/=0) return
      if (.not.ieee_is_finite(at_time).or..not.all(ieee_is_finite(ledger))) return
      if (at_time/=snapshot%time.or.info(1)/=snapshot%steps.or.info(3)/=snapshot%regrids) return
      if (any(ledger/=snapshot%initial)) return
      deallocate(drops,masses);allocate(drops(np),masses(nc,np),stat=ios)
      if (ios/=0) return
      do p=1,np
        read(unit,*,iostat=ios) drops(p)%id,drops(p)%parent_id,drops(p)%material_id, &
          drops(p)%position,drops(p)%velocity,drops(p)%mass,drops(p)%temperature,drops(p)%multiplicity, &
          drops(p)%deformation,drops(p)%deformation_rate,drops(p)%rt_time,drops(p)%khrt_initial_count, &
          drops(p)%khrt_shed_volume,masses(:,p)
        if (ios/=0) return
        if (.not.valid_spray_parcel(drops(p)).or..not.valid_component_inventory(drops(p),masses(:,p))) return
        if (drops(p)%id>mark) return
        if (p>1) then
          if (drops(p)%id<=drops(p-1)%id) return
        end if
      end do
      if (present(wall_films)) then
        read(unit,*,iostat=ios) nf
        if (ios/=0) return
        if (nf<0.or.nf>1000000) return
        deallocate(films,film_components)
        allocate(films(nf),film_components(nc,nf),stat=ios)
        if (ios/=0) return
        do p=1,nf
          read(unit,*,iostat=ios) films(p)%material_id,films(p)%position,films(p)%normal, &
            films(p)%mass,films(p)%temperature,films(p)%diameter,films(p)%height,film_components(:,p)
          if (ios/=0) return
          if (.not.valid_wall_component_inventory(films(p),film_components(:,p))) return
        end do
      end if
      if (present(source)) then
        read(unit,'(a)',iostat=ios) line
        if (ios/=0.or.line/='SPRAY_SOURCE_V1') return
        read(unit,*,iostat=ios) source_ledger
        if (ios/=0) return
        if (.not.all(ieee_is_finite(source_ledger))) return
      end if
      if(present(frame)) then
        read(unit,'(a)',iostat=ios) line
        if(ios/=0.or.line/='SPRAY_COORDINATE_ENERGY_V1') return
        read(unit,*,iostat=ios) coordinate_ledger
        if(ios/=0.or..not.ieee_is_finite(coordinate_ledger)) return
      endif
      read(unit,'(a)',iostat=ios) line
      if (ios/=0.or.line/='END_AMR_PARTICLES') return
      read(unit,'(a)',iostat=ios) line
      valid=ios<0
    end subroutine
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
    subroutine same_string(value,success)
      character(len=*), intent(in) :: value
      logical, intent(out) :: success
      character(len=:), allocatable :: reference
      integer :: length
      length=len(value)
      call MPI_Bcast(length,1,MPI_INTEGER,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) then
        success=.false.;return
      end if
      allocate(character(len=length)::reference);reference(:)=value
      call MPI_Bcast(reference,length,MPI_CHARACTER,0,comm,ierr)
      valid=ierr==MPI_SUCCESS.and.value==reference.and.len(value)==length
      call agree();success=agreed
    end subroutine
  end subroutine
end module
