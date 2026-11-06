module spray_mpi_checkpoint_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use nasa7_thermo_mod, only: nasa7_species
  use spray_parcel_mod, only: spray_parcel
  use spray_partition_mod, only: spray_cartesian_block
  use spray_mpi_mod, only: migrate_spray_parcels,migrate_mpi_wall_films
  use spray_wall_impact_mod, only: spray_wall_film
  use spray_checkpoint_mod, only: write_spray_checkpoint,read_spray_checkpoint
  implicit none
  private
  public :: mpi_spray_checkpoint
contains
  ! Root owns file I/O and the temporary whole-grid buffer; no full grid is
  ! replicated on other ranks. The caller binds all model tables in context.
  subroutine mpi_spray_checkpoint(comm,writing,path,context,species,state,temperature,parcels,component_mass, &
      lower,spacing,cells,time,steps,event,last_id,reference,injected,ok,breakup_state,khrt_state,shedding_state, &
      wall_films,film_mass,wall_totals,retired_parent_state,external_totals,reaction_source,heat_release,frame,previous_timestep)
    type(MPI_Comm), intent(in) :: comm
    type(rotating_frame_3d), optional, intent(in) :: frame
    logical, intent(in) :: writing
    logical, intent(in), optional :: breakup_state,khrt_state,shedding_state,retired_parent_state
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),time,reference(:),injected(:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3)
    integer, intent(in) :: cells(3)
    integer, intent(inout) :: steps,event
    integer(int64), intent(inout) :: last_id
    logical, intent(out) :: ok
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    real(dp), optional, intent(inout) :: wall_totals(:),external_totals(:)
    real(dp), optional, intent(inout) :: reaction_source(:,:,:,:),heat_release(:,:,:)
    real(dp), optional, intent(inout) :: previous_timestep
    real(dp), allocatable :: history
    real(dp), allocatable :: potential(:,:,:)
    real(dp) :: frame_controls(10),root_controls(10)
    integer :: i,j,k,frame_axis
    real(dp), allocatable :: global_ir(:,:,:,:),global_heat(:,:,:)
    integer :: fields,ir_first,ir_last,heat_field
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:),wall(:)
    logical :: has_films
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:),global_state(:,:,:,:),global_t(:,:,:),local(:,:,:,:),packed(:),received(:)
    real(dp), allocatable :: ref(:),inj(:),external(:)
    integer, allocatable :: counts(:),offsets(:)
    integer :: rank,nranks,ierr,nvar,nc,r,first(3),extent(3),config(17),low(17),high(17),step_data(2)
    integer(int64) :: next_id
    real(dp) :: next_time
    logical :: valid,agreed,has_breakup,has_khrt,has_shedding,has_retired
    character(len=256) :: message
    ok=.false.;has_breakup=.false.;has_khrt=.false.;has_shedding=.false.
    if (present(breakup_state)) has_breakup=breakup_state
    if (present(khrt_state)) has_khrt=khrt_state
    if (present(shedding_state)) has_shedding=shedding_state
    has_retired=.false.
    if (present(retired_parent_state)) has_retired=retired_parent_state
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    has_films=present(wall_films)
    nvar=size(state,1);nc=0
    valid=allocated(parcels).and.allocated(component_mass)
    if (allocated(component_mass)) nc=size(component_mass,1)
    valid=valid.and.nc>=1.and.nc<=16.and.nvar==5+size(species)
    valid=valid.and.size(reference)==nvar.and.size(injected)==nvar
    valid=valid.and.(has_films.eqv.present(film_mass)).and.(has_films.eqv.present(wall_totals))
    if(present(external_totals)) valid=valid.and.size(external_totals)==nvar
    fields=nvar+1;ir_first=fields+1;ir_last=fields
    if(present(reaction_source)) then
      fields=fields+nvar-4;ir_last=fields
    endif
    heat_field=fields+1
    if(present(heat_release)) fields=fields+1
    call agree()
    if (.not.agreed) return
    frame_axis=0
    if(present(frame)) frame_axis=frame%axis
    config=[cells,nvar,nc,merge(1,0,writing),merge(1,0,has_breakup),merge(1,0,has_khrt), &
      merge(1,0,has_shedding),merge(1,0,has_films),merge(1,0,has_retired),merge(1,0,present(external_totals)), &
      merge(1,0,present(reaction_source)),merge(1,0,present(heat_release)),merge(1,0,present(frame)), &
      frame_axis,merge(1,0,present(previous_timestep))]
    call MPI_Allreduce(config,low,17,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,high,17,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(low/=high)) return
    ! Optional collectives require agreement on argument presence first.
    if(present(external_totals)) then
      valid=all(ieee_is_finite(external_totals))
      call agree()
      if(.not.agreed) return
      external=external_totals
    endif
    if(present(previous_timestep)) then
      allocate(history);history=previous_timestep
    endif
    if(present(frame)) then
      frame_controls(:6)=[lower,spacing]
      frame_controls(7:10)=[frame%angular_speed,frame%axis_point]
      root_controls=frame_controls
      call MPI_Bcast(root_controls,10,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if(ierr/=MPI_SUCCESS) return
      valid=frame%is_valid().and.all(ieee_is_finite(frame_controls)).and.all(frame_controls==root_controls)
      call agree();if(.not.agreed) return
    endif
    if (has_films) then
      valid=allocated(wall_films).and.allocated(film_mass).and.size(wall_totals)==nvar
      call agree()
      if (.not.agreed) return
      valid=all(shape(film_mass)==[nc,size(wall_films)]).and.all(ieee_is_finite(wall_totals))
      call agree()
      if (.not.agreed) return
      wall=wall_totals
    end if
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    valid=valid.and.all(cells>=2).and.all(shape(state)==[nvar,extent]).and.all(shape(temperature)==extent)
    if(present(reaction_source)) valid=valid.and.all(shape(reaction_source)==[nvar-4,extent])
    if(present(heat_release)) valid=valid.and.all(shape(heat_release)==extent)
    ! Bound the root buffer and MPI integer counts before products/allocations.
    if (valid) then
      valid=int(fields,int64)<=100000000_int64/cells(1)/cells(2)/cells(3)
    end if
    call agree()
    if (.not.agreed) return
    allocate(counts(nranks),offsets(nranks));offsets=0
    do r=0,nranks-1
      call spray_cartesian_block(cells,[nranks,1,1],r,first,extent,valid)
      counts(r+1)=fields*product(extent)
      if (r>0) offsets(r+1)=offsets(r)+counts(r)
    end do
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    allocate(local(fields,extent(1),extent(2),extent(3)),packed(counts(rank+1)))
    allocate(received(merge(sum(counts),1,rank==0)))
    if (rank==0) then
      allocate(global_state(nvar,cells(1),cells(2),cells(3)),global_t(cells(1),cells(2),cells(3)))
      if(present(reaction_source)) allocate(global_ir(nvar-4,cells(1),cells(2),cells(3)))
      if(present(heat_release)) allocate(global_heat(cells(1),cells(2),cells(3)))
    end if
    ! The root reconstructs the same global potential independently of the
    ! current rank layout. Only gas storage/EOS uses it; particle energy is ordinary.
    if(present(frame)) then
      valid=.true.
      if(rank==0) then
        allocate(potential(cells(1),cells(2),cells(3)))
        do k=1,cells(3)
          do j=1,cells(2)
            do i=1,cells(1)
              call rotating_potential_3d(frame,lower+([i,j,k]-.5_dp)*spacing,potential(i,j,k),valid)
              if(.not.valid) exit
            enddo
            if(.not.valid) exit
          enddo
          if(.not.valid) exit
        enddo
      endif
      call agree();if(.not.agreed) return
    endif
    ref=reference;inj=injected;next_time=time;step_data=[steps,event];next_id=last_id
    if (writing) then
      valid=all(ieee_is_finite([reference,injected,time])).and.last_id>=0
      if(present(previous_timestep)) then
        valid=valid.and.ieee_is_finite(previous_timestep)
        valid=valid.and.previous_timestep>=0.and.previous_timestep<=time
        valid=valid.and.((steps==0).eqv.(previous_timestep==0))
      endif
      call agree()
      if (.not.agreed) return
      call metadata()
      if (.not.agreed) return
      valid=all(ref==reference).and.all(inj==injected).and.next_time==time.and.next_id==last_id
      valid=valid.and.all(step_data==[steps,event])
      if(present(previous_timestep)) valid=valid.and.history==previous_timestep
      if (has_films) valid=valid.and.all(wall==wall_totals)
      if(present(external_totals)) valid=valid.and.all(external==external_totals)
      call agree()
      if (.not.agreed) return
      drops=parcels;masses=component_mass
      call migrate_spray_parcels(comm,drops,lower,spacing,cells,[nranks,1,1],[.true.,.true.,.true.], &
        1000000,valid,masses,gather_rank=0)
      if (.not.valid) return
      if (has_films) then
        patches=wall_films;patch_mass=film_mass
        call migrate_mpi_wall_films(comm,patches,patch_mass,lower,spacing,cells,[nranks,1,1],1000000,valid,gather_rank=0)
        if (.not.valid) return
      end if
      local(:nvar,:,:,:)=state;local(nvar+1,:,:,:)=temperature
      if(present(reaction_source)) local(ir_first:ir_last,:,:,:)=reaction_source
      if(present(heat_release)) local(heat_field,:,:,:)=heat_release
      valid=all(ieee_is_finite(local));call agree();if(.not.agreed) return
      packed=reshape(local,[size(packed)])
      call MPI_Gatherv(packed,size(packed),MPI_DOUBLE_PRECISION,received,counts,offsets, &
        MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      valid=.true.
      if (rank==0) then
        do r=0,nranks-1
          call spray_cartesian_block(cells,[nranks,1,1],r,first,extent,valid)
          local=reshape(received(offsets(r+1)+1:offsets(r+1)+counts(r+1)),[fields,extent])
          global_state(:,first(1)+1:first(1)+extent(1),:,:)=local(:nvar,:,:,:)
          global_t(first(1)+1:first(1)+extent(1),:,:)=local(nvar+1,:,:,:)
          if(present(reaction_source)) global_ir(:,first(1)+1:first(1)+extent(1),:,:)=local(ir_first:ir_last,:,:,:)
          if(present(heat_release)) global_heat(first(1)+1:first(1)+extent(1),:,:)=local(heat_field,:,:,:)
        end do
        if (has_films) then
          call write_spray_checkpoint(path,context,species,global_state,global_t,drops,next_time, &
              step_data(1),step_data(2),ref,inj,valid,message,component_mass=masses,last_id=next_id, &
            breakup_state=has_breakup,khrt_state=has_khrt,split_state=has_shedding,shedding_state=has_shedding, &
            retired_parent_state=has_retired, &
            film_state=.true.,wall_films=patches,film_component_mass=patch_mass,boundary_totals=wall,external_totals=external, &
            reaction_source=global_ir,heat_release=global_heat,potential=potential,previous_timestep=history)
        else
          call write_spray_checkpoint(path,context,species,global_state,global_t,drops,next_time, &
              step_data(1),step_data(2),ref,inj,valid,message,component_mass=masses,last_id=next_id, &
            breakup_state=has_breakup,khrt_state=has_khrt,split_state=has_shedding,shedding_state=has_shedding, &
            retired_parent_state=has_retired,external_totals=external, &
            reaction_source=global_ir,heat_release=global_heat,potential=potential,previous_timestep=history)
        end if
      end if
      call agree()
      ok=agreed;return
    end if
    allocate(drops(0),masses(nc,0));valid=.true.
    if (has_films) allocate(patches(0),patch_mass(nc,0))
    if (rank==0) then
      if (has_films) then
        call read_spray_checkpoint(path,context,species,global_state,global_t,drops,next_time, &
          step_data(1),step_data(2),ref,inj,valid,message,component_mass=masses,last_id=next_id, &
          breakup_state=has_breakup,khrt_state=has_khrt,split_state=has_shedding,shedding_state=has_shedding, &
            retired_parent_state=has_retired, &
          film_state=.true.,wall_films=patches,film_component_mass=patch_mass,boundary_totals=wall,external_totals=external, &
            reaction_source=global_ir,heat_release=global_heat,potential=potential,previous_timestep=history)
      else
        call read_spray_checkpoint(path,context,species,global_state,global_t,drops,next_time, &
          step_data(1),step_data(2),ref,inj,valid,message,component_mass=masses,last_id=next_id, &
          breakup_state=has_breakup,khrt_state=has_khrt,split_state=has_shedding,shedding_state=has_shedding, &
            retired_parent_state=has_retired,external_totals=external, &
            reaction_source=global_ir,heat_release=global_heat,potential=potential,previous_timestep=history)
      end if
      if (valid) valid=size(masses,1)==nc
      if (valid) then
        do r=0,nranks-1
          call spray_cartesian_block(cells,[nranks,1,1],r,first,extent,valid)
          if (allocated(local)) deallocate(local)
          allocate(local(fields,extent(1),extent(2),extent(3)))
          local(:nvar,:,:,:)=global_state(:,first(1)+1:first(1)+extent(1),:,:)
          local(nvar+1,:,:,:)=global_t(first(1)+1:first(1)+extent(1),:,:)
          if(present(reaction_source)) local(ir_first:ir_last,:,:,:)=global_ir(:,first(1)+1:first(1)+extent(1),:,:)
          if(present(heat_release)) local(heat_field,:,:,:)=global_heat(first(1)+1:first(1)+extent(1),:,:)
          received(offsets(r+1)+1:offsets(r+1)+counts(r+1))=reshape(local,[counts(r+1)])
        end do
      end if
    end if
    call agree()
    if (.not.agreed) return
    call MPI_Scatterv(received,counts,offsets,MPI_DOUBLE_PRECISION,packed,size(packed),MPI_DOUBLE_PRECISION,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call migrate_spray_parcels(comm,drops,lower,spacing,cells,[nranks,1,1],[.true.,.true.,.true.], &
      1000000,valid,masses)
    if (.not.valid) return
    if (has_films) then
      call migrate_mpi_wall_films(comm,patches,patch_mass,lower,spacing,cells,[nranks,1,1],1000000,valid)
      if (.not.valid) return
    end if
    call metadata()
    if (.not.agreed) return
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    local=reshape(packed,[fields,extent])
    valid=all(ieee_is_finite(local));call agree();if(.not.agreed) return
    state=local(:nvar,:,:,:);temperature=local(nvar+1,:,:,:)
    if(present(reaction_source)) reaction_source=local(ir_first:ir_last,:,:,:)
    if(present(heat_release)) heat_release=local(heat_field,:,:,:)
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    if (has_films) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass);wall_totals=wall
    end if
    if(present(external_totals)) external_totals=external
    if(present(previous_timestep)) previous_timestep=history
    reference=ref;injected=inj;time=next_time;steps=step_data(1);event=step_data(2);last_id=next_id;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
    subroutine metadata()
      agreed=.false.
      call MPI_Bcast(ref,nvar,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Bcast(inj,nvar,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      if(present(previous_timestep)) then
        call MPI_Bcast(history,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
        if(ierr/=MPI_SUCCESS) return
      endif
      call MPI_Bcast(next_time,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Bcast(step_data,2,MPI_INTEGER,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Bcast(next_id,1,MPI_INTEGER8,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      if (has_films) then
        call MPI_Bcast(wall,nvar,MPI_DOUBLE_PRECISION,0,comm,ierr)
        if (ierr/=MPI_SUCCESS) return
      end if
      if(present(external_totals)) then
        call MPI_Bcast(external,nvar,MPI_DOUBLE_PRECISION,0,comm,ierr)
        if(ierr/=MPI_SUCCESS) return
      endif
      agreed=.true.
    end subroutine
  end subroutine
end module
