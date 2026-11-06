module les_mpi_boundary_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use spray_partition_mod, only: spray_cartesian_block
  use les_mpi_filter_mod, only: sample_mpi_les_x_planes_3d
  use simulation_config_reactive_3d_mod, only: map_reactive_boundary_donor_3d, &
    valid_reactive_boundaries_3d,valid_wall_velocity_3d,fill_pressure_outlet_ghosts_3d,mach_inlet_ghost_3d
  implicit none
  private
  public :: fill_mpi_les_boundary_ghosts_3d
  public :: fill_mpi_sdc_source_ghosts_3d
contains
  ! Signed SDC predictors use periodic wrapping and physical extrapolation,
  ! matching reactive_3d_mod source halos. They are rates, not EOS states.
  subroutine fill_mpi_sdc_source_ghosts_3d(comm,source,cells,faces,width,ghost,ok)
    type(MPI_Comm), intent(in) :: comm
    real(dp), intent(in) :: source(:,:,:,:)
    integer, intent(in) :: cells(3),width
    character(len=*), intent(in) :: faces(6)
    real(dp), allocatable, intent(inout) :: ghost(:,:,:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: planes(:,:,:,:),candidate(:,:,:,:)
    integer, allocatable :: coordinates(:)
    integer :: rank,nranks,status,ios,nc,first(3),extent(3),dims(3),config(5),low(5),high(5)
    integer :: i,j,k,q(3),plane
    character(len=32) :: reference_faces(6)
    logical :: valid,agreed
    ok=.false.;nc=size(source,1)
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,status)
    if(status/=MPI_SUCCESS) return
    config=[cells,width,nc]
    call MPI_Allreduce(config,low,5,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(config,high,5,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    valid=all(low==high).and.all(cells>=2).and.all(cells<=huge(1)-256).and. &
      width>=1.and.width<=128.and.nc>=1.and.valid_reactive_boundaries_3d(faces)
    call agree();if(.not.agreed) return
    reference_faces=faces
    call MPI_Bcast(reference_faces,6*len(reference_faces),MPI_CHARACTER,0,comm,status)
    if(status/=MPI_SUCCESS) return
    valid=all(faces==reference_faces);call agree();if(.not.agreed) return
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    valid=valid.and.all(extent>=1).and.all(shape(source)==[nc,extent]).and.all(ieee_is_finite(source))
    call agree();if(.not.agreed) return
    dims=extent+2*width
    allocate(coordinates(2*width),candidate(nc,dims(1),dims(2),dims(3)),stat=ios)
    valid=ios==0;call agree();if(.not.agreed) return
    do i=1,2*width
      plane=i
      if(i>width) plane=plane+extent(1)
      coordinates(i)=mapped(first(1)+plane-width,1)
    enddo
    call sample_mpi_les_x_planes_3d(comm,source,cells,coordinates,planes,valid)
    if(.not.valid) return
    do k=1,dims(3)
    do j=1,dims(2)
      q(2)=mapped(j-width,2);q(3)=mapped(k-width,3)
      do i=1,dims(1)
        if(i<=width) then
          candidate(:,i,j,k)=planes(:,i,q(2),q(3))
        else if(i>width+extent(1)) then
          candidate(:,i,j,k)=planes(:,i-extent(1),q(2),q(3))
        else
          candidate(:,i,j,k)=source(:,i-width,q(2),q(3))
        endif
      enddo
    enddo
    enddo
    valid=all(ieee_is_finite(candidate));call agree();if(.not.agreed) return
    call move_alloc(candidate,ghost);ok=.true.
  contains
    integer function mapped(index,axis) result(donor)
      integer, intent(in) :: index,axis
      if(faces(2*axis-1)=='periodic') then
        donor=modulo(index-1,cells(axis))+1
      else
        donor=max(1,min(cells(axis),index))
      endif
    end function
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine fill_mpi_sdc_source_ghosts_3d

  subroutine fill_mpi_les_boundary_ghosts_3d(comm,species,state,temperature,cells,faces,width,ghost,gt,ok, &
      inflow_state,inflow_temperature,wall_velocity,outlet_pressure,inflow_mach,inflow_pressure,viscous_wall, &
      internal_energy_density,internal_energy_ghost)
    real(dp), optional, intent(in) :: internal_energy_density(:,:,:)
    real(dp), allocatable, optional, intent(inout) :: internal_energy_ghost(:,:,:)
    real(dp), allocatable :: candidate_h(:,:,:),held_cell
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:)
    integer, intent(in) :: cells(3),width
    character(len=*), intent(in) :: faces(6)
    real(dp), allocatable, intent(inout) :: ghost(:,:,:,:),gt(:,:,:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: viscous_wall
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_velocity(:,:),outlet_pressure(6)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_pressure(6)
    real(dp), allocatable :: packed(:,:,:,:),planes(:,:,:,:),donors(:,:,:,:),candidate(:,:,:,:),temps(:,:,:)
    real(dp), allocatable :: inlet(:,:),reference_inlet(:,:)
    integer, allocatable :: coordinates(:)
    integer :: rank,nranks,ierr,ios,nv,np,first(3),extent(3),dims(3),config(8),low(8),high(8)
    integer :: i,j,k,q(3),idx(3),inlet_face,a,bounds(2,3)
    real(dp) :: walls(3,6),pressure(6),controls(36),minimum(36),maximum(36),mach(6),pin(6),shift(3),momentum(3),rho
    character(len=32) :: reference_faces(6)
    logical :: valid,agreed,reflect(3),viscous
    ok=.false.;nv=5+size(species);walls=0;pressure=0;mach=-1;pin=0
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if(ierr/=MPI_SUCCESS) return
    viscous=.true.
    if(present(viscous_wall)) viscous=viscous_wall
    config=[cells,width,nv,merge(1,0,viscous),merge(1,0,present(internal_energy_density)), &
      merge(1,0,present(internal_energy_ghost))]
    call MPI_Allreduce(config,low,8,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,high,8,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=all(low==high).and.all(cells>=2).and.all(cells<=huge(1)-256).and.width>=1.and.width<=128
    valid=valid.and.valid_reactive_boundaries_3d(faces)
    valid=valid.and.(present(internal_energy_density).eqv.present(internal_energy_ghost))
    if(present(wall_velocity)) then
      valid=valid.and.valid_wall_velocity_3d(faces,wall_velocity)
      if(valid) walls=wall_velocity
    end if
    if(present(outlet_pressure)) pressure=outlet_pressure
    valid=valid.and.all(ieee_is_finite(pressure)).and.all(pressure>=0).and.all(pressure<=1e15_dp)
    valid=valid.and..not.any(pressure>0.and.faces/='outflow')
    if(present(inflow_mach)) mach=inflow_mach
    if(present(inflow_pressure)) then
      pin=inflow_pressure
      valid=valid.and.present(inflow_mach)
    end if
    valid=valid.and.all(ieee_is_finite(mach)).and.all(mach<=100)
    valid=valid.and..not.any(mach>=0.and.faces/='inflow')
    valid=valid.and.all(ieee_is_finite(pin)).and.all(pin>=0).and.all(pin<=1e15_dp)
    valid=valid.and..not.any(pin>0.and.mach<0)
    if(any(faces=='inflow')) then
      valid=valid.and.present(inflow_temperature)
      if(present(inflow_temperature)) valid=valid.and.size(inflow_temperature)==6
    end if
    if(any(faces=='inflow'.and.mach<0)) then
      valid=valid.and.present(inflow_state)
      if(present(inflow_state)) valid=valid.and.all(shape(inflow_state)==[nv,6])
      if(present(inflow_temperature)) valid=valid.and.size(inflow_temperature)==6
    end if
    call agree()
    if(.not.agreed) return
    reference_faces=faces
    call MPI_Bcast(reference_faces,6*len(reference_faces),MPI_CHARACTER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=all(faces==reference_faces);call agree()
    if(.not.agreed) return
    controls(:18)=reshape(walls,[18]);controls(19:24)=pressure;controls(25:30)=mach;controls(31:)=pin
    call MPI_Allreduce(controls,minimum,36,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,36,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum/=maximum)) return
    allocate(inlet(nv+1,6),reference_inlet(nv+1,6),stat=ios)
    valid=ios==0;call agree()
    if(.not.agreed) return
    inlet=0
    if(any(faces=='inflow')) then
      if(any(faces=='inflow'.and.mach<0)) inlet(:nv,:)=inflow_state
      inlet(nv+1,:)=inflow_temperature
    end if
    valid=all(ieee_is_finite(inlet));call agree()
    if(.not.agreed) return
    reference_inlet=inlet
    call MPI_Bcast(reference_inlet,size(inlet),MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=all(inlet==reference_inlet);call agree()
    if(.not.agreed) return
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    valid=valid.and.all(extent>=1).and.all(shape(state)==[nv,extent]).and.all(shape(temperature)==extent)
    valid=valid.and.all(ieee_is_finite(state)).and.all(ieee_is_finite(temperature)).and.all(temperature>0)
    call agree()
    if(.not.agreed) return
    if(present(internal_energy_density)) then
      valid=all(shape(internal_energy_density)==extent)
      if(valid) valid=all(ieee_is_finite(internal_energy_density))
      call agree();if(.not.agreed) return
    endif
    np=nv+1+merge(1,0,present(internal_energy_density))
    dims=extent+2*width
    if(present(internal_energy_density)) then
      allocate(candidate_h(dims(1),dims(2),dims(3)),held_cell,stat=ios)
      valid=ios==0;call agree();if(.not.agreed) return
    endif
    allocate(packed(np,extent(1),extent(2),extent(3)),coordinates(2*width), &
      donors(np,dims(1),extent(2),extent(3)),candidate(nv,dims(1),dims(2),dims(3)), &
      temps(dims(1),dims(2),dims(3)),stat=ios)
    valid=ios==0;call agree()
    if(.not.agreed) return
    packed(:nv,:,:,:)=state;packed(nv+1,:,:,:)=temperature
    if(present(internal_energy_density)) packed(np,:,:,:)=internal_energy_density
    do i=1,2*width
      j=i
      if(i>width) j=i+extent(1)
      idx=[first(1)+j-width,1,1]
      call map_reactive_boundary_donor_3d(idx,cells,faces,viscous,q,reflect,shift,inlet_face,valid,walls)
      if(.not.valid) exit
      coordinates(i)=q(1)
    end do
    call agree()
    if(.not.agreed) return
    call sample_mpi_les_x_planes_3d(comm,packed,cells,coordinates,planes,valid)
    if(.not.valid) return
    donors(:,:width,:,:)=planes(:,:width,:,:)
    donors(:,width+1:width+extent(1),:,:)=packed
    donors(:,width+extent(1)+1:,:,:)=planes(:,width+1:,:,:)
    do k=1,dims(3)
      do j=1,dims(2)
        do i=1,dims(1)
          idx=[first(1)+i-width,j-width,k-width]
          call map_reactive_boundary_donor_3d(idx,cells,faces,viscous,q,reflect,shift,inlet_face,valid,walls)
          if(.not.valid) exit
          candidate(:,i,j,k)=donors(:nv,i,q(2),q(3));temps(i,j,k)=donors(nv+1,i,q(2),q(3))
          if(allocated(held_cell)) then
            held_cell=donors(np,i,q(2),q(3));candidate_h(i,j,k)=held_cell
          endif
          if(inlet_face>0) then
            if(mach(inlet_face)>=0) then
              if(pin(inlet_face)>0) then
                call mach_inlet_ghost_3d(species,donors(:nv,i,q(2),q(3)),donors(nv+1,i,q(2),q(3)), &
                  inlet_face,inlet(nv+1,inlet_face),mach(inlet_face),candidate(:,i,j,k),temps(i,j,k),valid, &
                  inlet_pressure=pin(inlet_face),internal_energy_density=held_cell)
              else
                call mach_inlet_ghost_3d(species,donors(:nv,i,q(2),q(3)),donors(nv+1,i,q(2),q(3)), &
                  inlet_face,inlet(nv+1,inlet_face),mach(inlet_face),candidate(:,i,j,k),temps(i,j,k),valid, &
                  internal_energy_density=held_cell)
              end if
              if(.not.valid) exit
            else
              candidate(:,i,j,k)=inlet(:nv,inlet_face);temps(i,j,k)=inlet(nv+1,inlet_face)
            end if
          end if
          if(allocated(candidate_h).and.inlet_face>0) then
            rho=candidate(1,i,j,k)
            valid=rho>0
            if(.not.valid) exit
            candidate_h(i,j,k)=candidate(5,i,j,k)-.5_dp*sum(candidate(2:4,i,j,k)**2)/rho
            valid=ieee_is_finite(candidate_h(i,j,k))
            if(.not.valid) exit
          endif
          momentum=candidate(2:4,i,j,k)
          do a=1,3
            if(reflect(a)) candidate(a+1,i,j,k)=-candidate(a+1,i,j,k)
          end do
          if(any(abs(shift)>0)) then
            rho=candidate(1,i,j,k)
            valid=rho>0
            if(.not.valid) exit
            candidate(2:4,i,j,k)=candidate(2:4,i,j,k)+rho*shift
            candidate(5,i,j,k)=candidate(5,i,j,k)+(sum(candidate(2:4,i,j,k)**2)-sum(momentum**2))/(2*rho)
          end if
        end do
        if(.not.valid) exit
      end do
      if(.not.valid) exit
    end do
    if(valid) valid=all(ieee_is_finite(candidate)).and.all(ieee_is_finite(temps))
    call agree()
    if(.not.agreed) return
    bounds(1,:)=width+1;bounds(2,:)=width+extent
    if(first(1)<width) then
      bounds(1,1)=width+1-first(1)
    else
      pressure(1)=0
    end if
    if(cells(1)-first(1)-extent(1)<width) then
      bounds(2,1)=cells(1)-first(1)+width
    else
      pressure(2)=0
    end if
    call fill_pressure_outlet_ghosts_3d(species,extent,faces,pressure,candidate,temps,valid, &
      ghost_width=width,boundary_cells=bounds,internal_energy_ghost=candidate_h)
    call agree()
    if(.not.agreed) return
    if(present(internal_energy_ghost)) call move_alloc(candidate_h,internal_energy_ghost)
    call move_alloc(candidate,ghost);call move_alloc(temps,gt);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
end module les_mpi_boundary_mod
