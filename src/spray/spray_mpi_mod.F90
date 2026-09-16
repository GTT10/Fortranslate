module spray_mpi_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use spray_reactive_3d_mod, only: advance_spray_chemistry_3d
  use reactive_transport_3d_mod, only: reactive_transport_ghosted_fluxes_3d
  use les_smagorinsky_3d_mod, only: compute_les_ghosted_rhs_3d
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_conserved_to_primitive
  use reactive_3d_mod, only: recover_reactive_temperatures_3d,compute_reactive_plm_slab_face_fluxes_3d
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory,liquid_mixture_properties, &
    advance_mixture_frozen_gas,mixture_transfer_rates,mixture_trace_fraction,build_mixture_film
  use spray_breakup_mod, only: split_spray_parcels,advance_spray_breakup,advance_khrt_shedding, &
    construct_khrt_children,spray_khrt_release,spray_khrt_shedding_state
  use spray_coupling_3d_mod, only: particle_boundary_stencil,valid_spray_boundary_faces
  use spray_partition_mod
  implicit none
  private
  public :: migrate_spray_parcels,reserve_spray_ids,split_mpi_spray_parcels,inject_mpi_cone_parcels
  public :: advance_mpi_khrt_shedding
  public :: advance_mpi_spray_les
  public :: advance_mpi_spray_transport
  public :: advance_mpi_spray_hydro
  public :: advance_mpi_mixture_coupling,advance_mpi_mixture_reaction
  public :: sample_mpi_spray_field,deposit_mpi_spray_field,deposit_mpi_spray_gas
contains
  subroutine advance_mpi_spray_transport(comm,species,transport,state,temperature,spacing,cells,dt,enabled,minimum_theta,ok)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: enabled(4)
    real(dp), intent(out) :: minimum_theta
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),dt
    integer, intent(in) :: cells(3)
    logical, intent(out) :: ok
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),candidate(:,:,:,:),next_temp(:,:,:)
    real(dp), allocatable :: halo(:,:,:,:),halo_t(:,:,:),send(:,:,:),left(:,:,:),right(:,:,:)
    real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    real(dp), allocatable :: theta(:,:,:),theta_send(:,:),theta_left(:,:),theta_right(:,:)
    real(dp) :: controls(4),minimum(4),maximum(4),local_theta,global_theta,accepted_theta
    integer :: config(8),minconfig(8),maxconfig(8),rank,nranks,ierr,first(3),extent(3)
    integer :: nx,ny,nz,nvar,stage,i,j,k,jm,km,previous,next
    logical :: valid,agreed
    ok=.false.;minimum_theta=1;accepted_theta=1;nvar=size(state,1)
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(ieee_is_finite([spacing,dt])).and.size(transport)==size(species)
    valid=valid.and.all(spacing>=1e-12_dp).and.all(spacing<=1e9_dp).and.dt>=0.and.dt<=1e6_dp
    call agree()
    if (.not.agreed) return
    config=[cells,nvar,merge(1,0,enabled)];controls=[spacing,dt]
    call MPI_Allreduce(config,minconfig,8,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,maxconfig,8,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,minimum,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(minconfig/=maxconfig).or.any(minimum/=maximum)) return
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    valid=valid.and.all(extent>=1)
    valid=valid.and.all(shape(state)==[nvar,extent]).and.all(shape(temperature)==extent)
    ! Bound the packed halo MPI count before allocating or sending it.
    valid=valid.and.int(nvar+1,int64)*int(extent(2),int64)*int(extent(3),int64)<=huge(1)
    call agree()
    if (.not.agreed) return
    nx=extent(1);ny=extent(2);nz=extent(3)
    gas=state;candidate=state;allocate(temp(nx,ny,nz),next_temp(nx,ny,nz))
    call recover_reactive_temperatures_3d(species,gas,temperature,nx,ny,nz,temp,valid)
    call agree()
    if (.not.agreed) return
    if (dt<=tiny(1.0_dp).or..not.any(enabled(:3))) then
      ok=.true.;return
    end if
    allocate(halo(nvar,0:nx+1,0:ny+1,0:nz+1),halo_t(0:nx+1,0:ny+1,0:nz+1))
    allocate(send(nvar+1,ny,nz),left(nvar+1,ny,nz),right(nvar+1,ny,nz))
    allocate(fx(nvar,0:nx,ny,nz),fy(nvar,nx,0:ny,nz),fz(nvar,nx,ny,0:nz))
    allocate(theta(nx,ny,nz),theta_send(ny,nz),theta_left(ny,nz),theta_right(ny,nz))
    previous=modulo(rank-1,nranks);next=modulo(rank+1,nranks)
    do stage=1,2
      send(:nvar,:,:)=gas(:,nx,:,:);send(nvar+1,:,:)=temp(nx,:,:)
      call MPI_Sendrecv(send,size(send),MPI_DOUBLE_PRECISION,next,3311,left,size(left),MPI_DOUBLE_PRECISION, &
        previous,3311,comm,MPI_STATUS_IGNORE,ierr)
      if (ierr/=MPI_SUCCESS) return
      send(:nvar,:,:)=gas(:,1,:,:);send(nvar+1,:,:)=temp(1,:,:)
      call MPI_Sendrecv(send,size(send),MPI_DOUBLE_PRECISION,previous,3312,right,size(right),MPI_DOUBLE_PRECISION, &
        next,3312,comm,MPI_STATUS_IGNORE,ierr)
      if (ierr/=MPI_SUCCESS) return
      do k=0,nz+1
        km=modulo(k-1,nz)+1
        do j=0,ny+1
          jm=modulo(j-1,ny)+1
          halo(:,1:nx,j,k)=gas(:,:,jm,km);halo_t(1:nx,j,k)=temp(:,jm,km)
          halo(:,0,j,k)=left(:nvar,jm,km);halo_t(0,j,k)=left(nvar+1,jm,km)
          halo(:,nx+1,j,k)=right(:nvar,jm,km);halo_t(nx+1,j,k)=right(nvar+1,jm,km)
        end do
      end do
      call reactive_transport_ghosted_fluxes_3d(species,transport,halo,halo_t,nx,ny,nz, &
        spacing(1),spacing(2),spacing(3),dt,enabled(1),enabled(2),enabled(3),enabled(4), &
        fx,fy,fz,local_theta,valid,theta_cell_output=theta)
      call agree()
      if (.not.agreed) return
      theta_send=theta(nx,:,:)
      call MPI_Sendrecv(theta_send,size(theta_send),MPI_DOUBLE_PRECISION,next,3313, &
        theta_left,size(theta_left),MPI_DOUBLE_PRECISION,previous,3313,comm,MPI_STATUS_IGNORE,ierr)
      if (ierr/=MPI_SUCCESS) return
      theta_send=theta(1,:,:)
      call MPI_Sendrecv(theta_send,size(theta_send),MPI_DOUBLE_PRECISION,previous,3314, &
        theta_right,size(theta_right),MPI_DOUBLE_PRECISION,next,3314,comm,MPI_STATUS_IGNORE,ierr)
      if (ierr/=MPI_SUCCESS) return
      call reactive_transport_ghosted_fluxes_3d(species,transport,halo,halo_t,nx,ny,nz, &
        spacing(1),spacing(2),spacing(3),dt,enabled(1),enabled(2),enabled(3),enabled(4),fx,fy,fz,local_theta,valid, &
        exterior_theta_x_lower=theta_left,exterior_theta_x_upper=theta_right, &
        periodic_theta_y=.true.,periodic_theta_z=.true.)
      call agree()
      if (.not.agreed) return
      call MPI_Allreduce(local_theta,global_theta,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      accepted_theta=min(accepted_theta,global_theta)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            candidate(:,i,j,k)=gas(:,i,j,k)-dt/spacing(1)*(fx(:,i,j,k)-fx(:,i-1,j,k)) &
              -dt/spacing(2)*(fy(:,i,j,k)-fy(:,i,j-1,k))-dt/spacing(3)*(fz(:,i,j,k)-fz(:,i,j,k-1))
          end do
        end do
      end do
      call recover_reactive_temperatures_3d(species,candidate,temp,nx,ny,nz,next_temp,valid)
      call agree()
      if (.not.agreed) return
      gas=candidate;temp=next_temp
    end do
    candidate=.5_dp*state+.5_dp*gas;temp=.5_dp*temperature+.5_dp*next_temp
    call recover_reactive_temperatures_3d(species,candidate,temp,nx,ny,nz,next_temp,valid)
    call agree()
    if (.not.agreed) return
    state=candidate;temperature=next_temp;minimum_theta=accepted_theta;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine advance_mpi_spray_les(comm,species,state,temperature,spacing,cells,dt,cs,prt,sct,ok)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),dt
    integer, intent(in) :: cells(3)
    real(dp), intent(in) :: cs,prt,sct
    logical, intent(out) :: ok
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),candidate(:,:,:,:),next_temp(:,:,:)
    real(dp), allocatable :: halo(:,:,:,:),halo_t(:,:,:),send(:,:,:,:),left(:,:,:,:),right(:,:,:,:)
    real(dp), allocatable :: rhs(:,:,:,:)
    real(dp) :: controls(7),minimum(7),maximum(7),diffusivity
    integer :: config(4),minconfig(4),maxconfig(4),rank,nranks,ierr,first(3),extent(3)
    integer :: nx,ny,nz,nvar,stage,i,j,k,jm,km,previous,next
    logical :: valid,agreed
    ok=.false.;nvar=size(state,1)
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(ieee_is_finite([spacing,dt,cs,prt,sct]))
    valid=valid.and.all(spacing>=1e-12_dp).and.all(spacing<=1e9_dp).and.dt>=0.and.dt<=1e6_dp
    call agree()
    if (.not.agreed) return
    config=[cells,nvar];controls=[spacing,dt,cs,prt,sct]
    call MPI_Allreduce(config,minconfig,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,maxconfig,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,minimum,7,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,7,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(minconfig/=maxconfig).or.any(minimum/=maximum)) return
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    valid=valid.and.extent(1)>=2.and.all(extent(2:)>=3)
    valid=valid.and.all(shape(state)==[nvar,extent]).and.all(shape(temperature)==extent)
    ! Bound the packed halo MPI count before allocating or sending it.
    valid=valid.and.int(nvar+1,int64)*2*int(extent(2),int64)*int(extent(3),int64)<=huge(1)
    call agree()
    if (.not.agreed) return
    nx=extent(1);ny=extent(2);nz=extent(3)
    gas=state;allocate(temp(nx,ny,nz),next_temp(nx,ny,nz));candidate=gas
    call recover_reactive_temperatures_3d(species,gas,temperature,nx,ny,nz,temp,valid)
    call agree()
    if (.not.agreed) return
    allocate(halo(nvar,nx+4,ny+4,nz+4),halo_t(nx+4,ny+4,nz+4))
    allocate(send(nvar+1,2,ny,nz),left(nvar+1,2,ny,nz),right(nvar+1,2,ny,nz))
    allocate(rhs,mold=state)
    previous=modulo(rank-1,nranks);next=modulo(rank+1,nranks)
    do stage=1,2
      send(:nvar,:,:,:)=gas(:,nx-1:nx,:,:);send(nvar+1,:,:,:)=temp(nx-1:nx,:,:)
      call MPI_Sendrecv(send,size(send),MPI_DOUBLE_PRECISION,next,3321,left,size(left),MPI_DOUBLE_PRECISION, &
        previous,3321,comm,MPI_STATUS_IGNORE,ierr)
      if (ierr/=MPI_SUCCESS) return
      send(:nvar,:,:,:)=gas(:,:2,:,:);send(nvar+1,:,:,:)=temp(:2,:,:)
      call MPI_Sendrecv(send,size(send),MPI_DOUBLE_PRECISION,previous,3322,right,size(right),MPI_DOUBLE_PRECISION, &
        next,3322,comm,MPI_STATUS_IGNORE,ierr)
      if (ierr/=MPI_SUCCESS) return
      do k=1,nz+4
        km=modulo(k-3,nz)+1
        do j=1,ny+4
          jm=modulo(j-3,ny)+1
          halo(:,3:nx+2,j,k)=gas(:,:,jm,km);halo_t(3:nx+2,j,k)=temp(:,jm,km)
          halo(:,:2,j,k)=left(:nvar,:,jm,km);halo_t(:2,j,k)=left(nvar+1,:,jm,km)
          halo(:,nx+3:,j,k)=right(:nvar,:,jm,km);halo_t(nx+3:,j,k)=right(nvar+1,:,jm,km)
        end do
      end do
      call compute_les_ghosted_rhs_3d(species,halo,halo_t,[nx,ny,nz],spacing,cs,prt,sct,rhs,diffusivity,valid)
      if (valid) valid=dt*diffusivity*sum(1.0_dp/spacing**2)<=.2_dp
      call agree()
      if (.not.agreed) return
      candidate=gas+dt*rhs
      if (stage==2) candidate=.5_dp*(state+gas+dt*rhs)
      call recover_reactive_temperatures_3d(species,candidate,temp,nx,ny,nz,next_temp,valid)
      call agree()
      if (.not.agreed) return
      gas=candidate;temp=next_temp
    end do
    state=gas;temperature=temp;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine advance_mpi_spray_hydro(comm,species,state,temperature,spacing,cells,dt,solver,limiter,ok)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),dt
    integer, intent(in) :: cells(3)
    character(len=*), intent(in) :: solver,limiter
    logical, intent(out) :: ok
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),candidate(:,:,:,:),next_temp(:,:,:)
    real(dp), allocatable :: halo(:,:,:,:),halo_t(:,:,:),send(:,:,:,:),left(:,:,:,:),right(:,:,:,:)
    real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    real(dp) :: controls(4),minimum(4),maximum(4),rhs(size(state,1))
    integer :: config(6),minconfig(6),maxconfig(6),rank,nranks,ierr,first(3),extent(3)
    integer :: nx,ny,nz,nvar,stage,i,j,k,jm,km,previous,next,solver_code,limiter_code
    logical :: valid,agreed
    ok=.false.;nvar=size(state,1);solver_code=0;limiter_code=0
    select case(trim(solver))
    case('rusanov');solver_code=1
    case('hllc');solver_code=2
    case('pelec');solver_code=3
    end select
    select case(trim(limiter))
    case('minmod');limiter_code=1
    case('mc');limiter_code=2
    end select
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(ieee_is_finite([spacing,dt])).and.solver_code>0.and.limiter_code>0
    valid=valid.and.all(spacing>=1e-12_dp).and.all(spacing<=1e9_dp).and.dt>=0.and.dt<=1e6_dp
    call agree()
    if (.not.agreed) return
    config=[cells,nvar,solver_code,limiter_code];controls=[spacing,dt]
    call MPI_Allreduce(config,minconfig,6,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,maxconfig,6,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,minimum,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(minconfig/=maxconfig).or.any(minimum/=maximum)) return
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    valid=valid.and.extent(1)>=2.and.all(extent(2:)>=3)
    valid=valid.and.all(shape(state)==[nvar,extent]).and.all(shape(temperature)==extent)
    ! Bound the packed halo MPI count before allocating or sending it.
    valid=valid.and.int(nvar+1,int64)*2*int(extent(2),int64)*int(extent(3),int64)<=huge(1)
    call agree()
    if (.not.agreed) return
    nx=extent(1);ny=extent(2);nz=extent(3)
    gas=state;allocate(temp(nx,ny,nz),next_temp(nx,ny,nz));candidate=gas
    call recover_reactive_temperatures_3d(species,gas,temperature,nx,ny,nz,temp,valid)
    call agree()
    if (.not.agreed) return
    allocate(halo(nvar,nx+4,ny,nz),halo_t(nx+4,ny,nz))
    allocate(send(nvar+1,2,ny,nz),left(nvar+1,2,ny,nz),right(nvar+1,2,ny,nz))
    allocate(fx(nvar,nx+4,ny,nz),fy(nvar,nx+4,ny,nz),fz(nvar,nx+4,ny,nz))
    previous=modulo(rank-1,nranks);next=modulo(rank+1,nranks)
    do stage=1,2
      send(:nvar,:,:,:)=gas(:,nx-1:nx,:,:);send(nvar+1,:,:,:)=temp(nx-1:nx,:,:)
      call MPI_Sendrecv(send,size(send),MPI_DOUBLE_PRECISION,next,3301,left,size(left),MPI_DOUBLE_PRECISION, &
        previous,3301,comm,MPI_STATUS_IGNORE,ierr)
      if (ierr/=MPI_SUCCESS) return
      send(:nvar,:,:,:)=gas(:,:2,:,:);send(nvar+1,:,:,:)=temp(:2,:,:)
      call MPI_Sendrecv(send,size(send),MPI_DOUBLE_PRECISION,previous,3302,right,size(right),MPI_DOUBLE_PRECISION, &
        next,3302,comm,MPI_STATUS_IGNORE,ierr)
      if (ierr/=MPI_SUCCESS) return
      halo(:,3:nx+2,:,:)=gas;halo_t(3:nx+2,:,:)=temp
      halo(:,:2,:,:)=left(:nvar,:,:,:);halo_t(:2,:,:)=left(nvar+1,:,:,:)
      halo(:,nx+3:,:,:)=right(:nvar,:,:,:);halo_t(nx+3:,:,:)=right(nvar+1,:,:,:)
      call compute_reactive_plm_slab_face_fluxes_3d(species,halo,halo_t,nx+4,ny,nz,2,nx+2, &
        limiter,solver,fx,fy,fz,valid)
      call agree()
      if (.not.agreed) return
      do k=1,nz
        km=modulo(k-2,nz)+1
        do j=1,ny
          jm=modulo(j-2,ny)+1
          do i=1,nx
            rhs=-(fx(:,i+2,j,k)-fx(:,i+1,j,k))/spacing(1) &
              -(fy(:,i+2,j,k)-fy(:,i+2,jm,k))/spacing(2) &
              -(fz(:,i+2,j,k)-fz(:,i+2,j,km))/spacing(3)
            candidate(:,i,j,k)=gas(:,i,j,k)+dt*rhs
            if (stage==2) candidate(:,i,j,k)=.5_dp*state(:,i,j,k)+.5_dp*candidate(:,i,j,k)
          end do
        end do
      end do
      call recover_reactive_temperatures_3d(species,candidate,temp,nx,ny,nz,next_temp,valid)
      call agree()
      if (.not.agreed) return
      gas=candidate;temp=next_temp
    end do
    state=gas;temperature=temp;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine advance_mpi_mixture_reaction(comm,species,reactions,transport,components,options,state,temperature, &
      parcels,component_mass,lower,spacing,cells,parts,interval,rtol,atol,method,substeps,ok, &
      hydro_solver,hydro_limiter,transport_flags,last_id,les_settings)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),interval,rtol,atol
    integer, intent(in) :: cells(3),parts(3)
    character(len=*), intent(in) :: method
    character(len=*), intent(in), optional :: hydro_solver,hydro_limiter
    logical, intent(in), optional :: transport_flags(4)
    real(dp), intent(in), optional :: les_settings(3)
    real(dp) :: sgs(3)
    integer, intent(out) :: substeps
    logical, intent(out) :: ok
    integer(int64), optional, intent(inout) :: last_id
    integer(int64) :: candidate_id
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),masses(:,:)
    real(dp) :: controls(13),minimum(13),maximum(13),chemistry_interval,transport_theta
    integer :: ierr,first_steps,last_steps,method_code,hydro_code,limit_code,nranks
    logical :: valid,agreed,enabled(4)
    character(len=256) :: message
    ok=.false.;substeps=0;method_code=0
    select case(trim(method))
    case('implicit');method_code=1
    case('explicit');method_code=2
    case('cvode');method_code=3
    end select
    sgs=[0.0_dp,.7_dp,.7_dp]
    if (present(les_settings)) sgs=les_settings
    enabled=.false.
    if (present(transport_flags)) enabled=transport_flags
    hydro_code=0;limit_code=0
    if (present(hydro_solver)) then
      select case(trim(hydro_solver))
      case('rusanov');hydro_code=1
      case('hllc');hydro_code=2
      case('pelec');hydro_code=3
      end select
    end if
    if (present(hydro_limiter)) then
      select case(trim(hydro_limiter))
      case('minmod');limit_code=1
      case('mc');limit_code=2
      end select
    end if
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=allocated(parcels).and.allocated(component_mass).and.method_code>0
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    valid=valid.and.all(ieee_is_finite([interval,rtol,atol,sgs]))
    valid=valid.and.sgs(1)>=0.and.sgs(1)<=2.and.minval(sgs(2:))>=1e-6_dp.and.maxval(sgs(2:))<=1e6_dp
    valid=valid.and.(sgs(1)==0.or.present(hydro_solver))
    valid=valid.and.interval>=0.and.interval<=1e6_dp.and.min(rtol,atol)>0
    valid=valid.and.(.not.any(enabled).or.present(hydro_solver))
    valid=valid.and.(present(hydro_solver).eqv.present(hydro_limiter))
    if (present(hydro_solver)) then
      valid=valid.and.hydro_code>0.and.limit_code>0.and.all(parts==[nranks,1,1])
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (.not.agreed) return
    controls=[interval,rtol,atol,real(method_code,dp),real(hydro_code,dp),real(limit_code,dp), &
      merge(1.0_dp,0.0_dp,enabled),sgs]
    call MPI_Allreduce(controls,minimum,13,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,13,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(minimum/=maximum)) return
    candidate_id=0
    if (present(last_id)) candidate_id=last_id
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    call advance_mpi_mixture_coupling(comm,species,transport,components,options,gas,temp,drops,masses, &
      lower,spacing,cells,parts,.5_dp*interval,first_steps,valid,candidate_id)
    if (.not.valid) return
    chemistry_interval=interval
    if (present(hydro_solver)) chemistry_interval=.5_dp*interval
    call advance_spray_chemistry_3d(species,reactions,gas,temp,chemistry_interval,rtol,atol,method,valid,message)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (.not.agreed) return
    if (present(hydro_solver)) then
      if (any(enabled)) then
        call advance_mpi_spray_transport(comm,species,transport,gas,temp,spacing,cells,.5_dp*interval, &
          enabled,transport_theta,valid)
        if (.not.valid) return
      end if
      if (sgs(1)>0) then
        call advance_mpi_spray_les(comm,species,gas,temp,spacing,cells,.5_dp*interval,sgs(1),sgs(2),sgs(3),valid)
        if (.not.valid) return
      end if
      call advance_mpi_spray_hydro(comm,species,gas,temp,spacing,cells,interval,hydro_solver,hydro_limiter,valid)
      if (.not.valid) return
      if (sgs(1)>0) then
        call advance_mpi_spray_les(comm,species,gas,temp,spacing,cells,.5_dp*interval,sgs(1),sgs(2),sgs(3),valid)
        if (.not.valid) return
      end if
      if (any(enabled)) then
        call advance_mpi_spray_transport(comm,species,transport,gas,temp,spacing,cells,.5_dp*interval, &
          enabled,transport_theta,valid)
        if (.not.valid) return
      end if
      call advance_spray_chemistry_3d(species,reactions,gas,temp,chemistry_interval,rtol,atol,method,valid,message)
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      if (.not.agreed) return
    end if
    call advance_mpi_mixture_coupling(comm,species,transport,components,options,gas,temp,drops,masses, &
      lower,spacing,cells,parts,.5_dp*interval,last_steps,valid,candidate_id)
    if (.not.valid) return
    if (first_steps>huge(substeps)-last_steps) return
    state=gas;temperature=temp
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    if (present(last_id)) last_id=candidate_id
    substeps=first_steps+last_steps;ok=.true.
  end subroutine

  ! Advance KH shedding and construct local children with globally reserved IDs.
  ! gas_energy(p) is the gas energy increment at original parent p (not a rate).
  ! The caller must deposit it transactionally with these candidates.
  subroutine advance_mpi_khrt_shedding(comm,components,options,parcels,component_mass,gas_density, &
      gas_velocity,reynolds,dt,last_id,gas_energy,ok)
    type(MPI_Comm), intent(in) :: comm
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:),gas_energy(:)
    real(dp), intent(in) :: gas_density(:),gas_velocity(:,:),reynolds(:),dt
    integer(int64), intent(inout) :: last_id
    logical, intent(out) :: ok
    type(spray_parcel), allocatable :: candidate(:),children(:),result(:)
    type(spray_khrt_release), allocatable :: releases(:)
    type(spray_khrt_shedding_state) :: pending
    type(spray_liquid) :: effective
    real(dp), allocatable :: masses(:,:),fractions(:,:),energy(:),result_mass(:,:)
    real(dp) :: rho,cp,e,axis(3),bounds(1),low(1),high(1)
    integer :: np,nc,p,n,number,cursor,ierr
    integer(int64) :: next_id,first_id
    logical :: valid,agreed
    ok=.false.
    valid=allocated(parcels).and.allocated(component_mass).and.valid_spray_options(options)
    valid=valid.and.options%breakup_model=='khrt'.and.options%khrt_children
    call agree()
    if (.not.agreed) return
    np=size(parcels);nc=size(components)
    valid=np<=1000000.and.nc>=1.and.nc<=16.and.all(shape(component_mass)==[nc,np])
    valid=valid.and.size(gas_density)==np.and.size(reynolds)==np.and.all(shape(gas_velocity)==[3,np])
    valid=valid.and.ieee_is_finite(dt).and.dt>=0.and.dt<=1e6_dp
    call agree()
    if (.not.agreed) return
    bounds=[dt]
    call MPI_Allreduce(bounds,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.any(low/=high)) return
    candidate=parcels;masses=component_mass
    allocate(releases(np),fractions(nc,np),energy(np));energy=0;fractions=0;number=0
    do p=1,np
      valid=valid_component_inventory(candidate(p),masses(:,p)).and.candidate(p)%id<=last_id
      if (p>1) valid=valid.and.candidate(p)%id>candidate(p-1)%id
      if (.not.valid) exit
      if (candidate(p)%mass<=0) cycle
      fractions(:,p)=masses(:,p)/candidate(p)%mass
      call liquid_mixture_properties(components,fractions(:,p),candidate(p)%temperature,rho,cp,e,valid)
      if (.not.valid) exit
      effective=components(1);effective%density_model='constant';effective%density=rho
      effective%density_coefficients=0;effective%density_temperature_range=0
      effective%heat_capacity_model='constant';effective%heat_capacity=cp
      effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
      effective%reference_temperature=candidate(p)%temperature;effective%reference_energy=e
      effective%boiling_model='reject';effective%critical_temperature=0
      pending=spray_khrt_shedding_state(candidate(p)%khrt_initial_count,candidate(p)%khrt_shed_volume)
      call advance_khrt_shedding(candidate(p),effective,options,gas_density(p),gas_velocity(:,p),reynolds(p), &
        dt,options%khrt_reference_diameter_cubed,pending,releases(p),valid)
      if (.not.valid) exit
      candidate(p)%khrt_initial_count=pending%initial_count;candidate(p)%khrt_shed_volume=pending%shed_volume
      masses(:,p)=candidate(p)%mass*fractions(:,p)
      if (releases(p)%mass>0) number=number+4
    end do
    valid=valid.and.number<=options%maximum_parcels-np
    call agree()
    if (.not.agreed) return
    next_id=last_id
    call reserve_spray_ids(comm,number,next_id,first_id,valid)
    if (.not.valid) return
    allocate(result(np+number),result_mass(nc,np+number))
    result(:np)=candidate;result_mass(:,:np)=masses;cursor=np
    do p=1,np
      if (releases(p)%mass<=0) cycle
      axis=candidate(p)%velocity
      if (norm2(axis)<1e-12_dp) axis=gas_velocity(:,p)-candidate(p)%velocity
      call construct_khrt_children(candidate(p),releases(p),axis,first_id+int(cursor-np,int64),children,energy(p),valid)
      if (.not.valid) exit
      result(cursor+1:cursor+4)=children
      do n=1,4
        result_mass(:,cursor+n)=children(n)%mass*fractions(:,p)
      end do
      cursor=cursor+4
    end do
    call agree()
    if (.not.agreed) return
    call move_alloc(result,parcels);call move_alloc(result_mass,component_mass);call move_alloc(energy,gas_energy)
    last_id=next_id;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine advance_mpi_mixture_coupling(comm,species,transport,components,options,state,temperature, &
      parcels,component_mass,lower,spacing,cells,parts,interval,substeps,ok,last_id)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),interval
    integer, intent(in) :: cells(3),parts(3)
    integer, intent(out) :: substeps
    logical, intent(out) :: ok
    integer(int64), optional, intent(inout) :: last_id
    integer(int64) :: accepted_id,trial_id
    real(dp), allocatable :: reynolds(:),child_energy(:)
    integer :: child_flag,low_flag,high_flag
    type(spray_parcel), allocatable :: drops(:),trial(:),survivors(:)
    type(spray_liquid) :: effective
    type(spray_film) :: film
    real(dp) :: old_mass
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),next_gas(:,:,:,:),next_temp(:,:,:)
    real(dp), allocatable :: masses(:,:),next_mass(:,:),kept_mass(:,:),positions(:,:),field(:,:,:,:),samples(:,:)
    real(dp), allocatable :: primitive(:,:),tg(:),change(:,:)
    real(dp) :: elapsed,dt,local_dt,min_interval,max_interval,sound,invtau,heat,rho,cp,e,tdot,speed
    real(dp) :: y(size(components)),rates(size(components)),hv(size(components)),length(3)
    integer :: p,n,np,nvar,nc,ierr,attempt,count_steps,kept,min_nc,max_nc
    logical :: valid,agreed,accepted
    character(len=8), parameter :: faces(6)='periodic'
    ok=.false.;substeps=0
    valid=allocated(parcels).and.allocated(component_mass).and.valid_spray_options(options)
    valid=valid.and.ieee_is_finite(interval)
    valid=valid.and.interval>=0.and.interval<=1e6_dp
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    valid=valid.and.size(components)>=1.and.size(components)<=16
    call agree()
    if (.not.agreed) return
    child_flag=merge(1,0,options%khrt_children)
    call MPI_Allreduce(child_flag,low_flag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(child_flag,high_flag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_flag/=high_flag) return
    accepted_id=0
    if (present(last_id)) accepted_id=last_id
    np=size(parcels);nc=size(components);nvar=size(state,1)
    valid=size(component_mass,1)==nc.and.size(component_mass,2)==np.and.np<=125000
    valid=valid.and.size(transport)==size(species)
    call MPI_Allreduce(interval,min_interval,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(interval,max_interval,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nc,min_nc,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nc,max_nc,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=valid.and.min_interval==max_interval.and.min_nc==max_nc
    call agree()
    if (.not.agreed) return
    do p=1,np
      valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p))
      if (p>1) valid=valid.and.parcels(p)%id>parcels(p-1)%id
    end do
    do n=1,nc
      valid=valid.and.valid_spray_liquid(components(n)).and.components(n)%vapor_index<=size(species)
    end do
    call agree()
    if (.not.agreed) return
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    allocate(positions(3,np),change(nvar,np),primitive(nvar,np),tg(np))
    do p=1,np
      positions(:,p)=drops(p)%position
    end do
    change=0
    ! This collective also validates geometry and the original gas on every rank.
    call deposit_mpi_spray_gas(comm,species,gas,temp,positions,change,lower,spacing,cells,parts,faces,valid)
    if (.not.valid) return
    length=spacing*real(cells,dp);elapsed=0;count_steps=0
    do while(elapsed<interval)
      valid=count_steps<options%maximum_substeps
      call agree()
      if (.not.agreed) return
      np=size(drops)
      deallocate(positions,change,primitive,tg)
      if (allocated(reynolds)) deallocate(reynolds)
      allocate(positions(3,np),change(nvar,np),primitive(nvar,np),tg(np),reynolds(np))
      ! Sample the temperature guess alongside conserved variables in one exchange.
      if (allocated(field)) deallocate(field)
      allocate(field(nvar+1,size(gas,2),size(gas,3),size(gas,4)))
      field(:nvar,:,:,:)=gas;field(nvar+1,:,:,:)=temp
      do p=1,np
        positions(:,p)=drops(p)%position
      end do
      call sample_mpi_spray_field(comm,field,positions,lower,spacing,cells,parts,faces,samples,valid)
      if (.not.valid) return
      local_dt=interval-elapsed
      do p=1,np
        if (drops(p)%mass<=0) cycle
        call reactive_conserved_to_primitive(species,samples(:nvar,p),samples(nvar+1,p),primitive(:,p),tg(p),sound,valid)
        if (.not.valid) exit
        y=masses(:,p)/drops(p)%mass
        call mixture_transfer_rates(species,transport,primitive(:,p),tg(p),drops(p),components,y,options, &
          invtau,heat,rates,hv,valid)
        if (.not.valid) exit
        call liquid_mixture_properties(components,y,drops(p)%temperature,rho,cp,e,valid)
        if (.not.valid) exit
        if (invtau>0) local_dt=min(local_dt,options%max_fractional_change/invtau)
        tdot=heat*(tg(p)-drops(p)%temperature)
        do n=1,nc
          if (rates(n)>0.and.y(n)>mixture_trace_fraction) &
            local_dt=min(local_dt,options%max_fractional_change*masses(n,p)/rates(n))
          tdot=tdot-rates(n)*(hv(n)-liquid_energy(components(n),drops(p)%temperature))
        end do
        tdot=tdot/(drops(p)%mass*cp)
        if (abs(tdot)>0) local_dt=min(local_dt,min(2.0_dp,options%max_fractional_change*drops(p)%temperature)/abs(tdot))
        speed=max(maxval(abs(drops(p)%velocity)),maxval(abs(primitive(2:4,p))))
        if (speed>0) local_dt=min(local_dt,options%particle_cfl*minval(spacing)/speed)
      end do
      call agree()
      if (.not.agreed) return
      call MPI_Allreduce(local_dt,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      accepted=.false.
      do attempt=1,30
        if (dt<=0.or.elapsed+dt<=elapsed) exit
        trial=drops;next_mass=masses;change=0;valid=.true.;trial_id=accepted_id;reynolds=0
        do p=1,np
          if (drops(p)%mass<=0) cycle
          call advance_mixture_frozen_gas(species,transport,primitive(:,p),tg(p),trial(p),components, &
            next_mass(:,p),options,dt,change(:,p),valid)
          if (.not.valid) exit
          if (trial(p)%mass>0.and.options%breakup_model/='none') then
            old_mass=trial(p)%mass;y=next_mass(:,p)/old_mass
            call liquid_mixture_properties(components,y,trial(p)%temperature,rho,cp,e,valid)
            if (.not.valid) exit
            effective=components(1);effective%density=rho;effective%density_model='constant'
            effective%density_coefficients=0;effective%density_temperature_range=0
            effective%heat_capacity=cp;effective%heat_capacity_model='constant'
            effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
            effective%reference_temperature=trial(p)%temperature;effective%reference_energy=e
            effective%boiling_model='reject';effective%critical_temperature=0
            call build_mixture_film(species,transport,primitive(:,p),tg(p),trial(p),components,y,options, &
              film,rates,hv,valid)
            if (.not.valid) exit
            if (options%khrt_children) then
              reynolds(p)=film%density*norm2(primitive(2:4,p)-trial(p)%velocity)* &
                parcel_diameter(trial(p),effective)/film%viscosity
            else
              call advance_spray_breakup(trial(p),effective,options,primitive(1,p),primitive(2:4,p), &
                film%density,film%viscosity,dt,valid)
              if (.not.valid) exit
              next_mass(:,p)=next_mass(:,p)*(trial(p)%mass/old_mass)
            end if
          end if
          trial(p)%position=lower+modulo(trial(p)%position-lower,length)
        end do
        call agree()
        if (agreed.and.options%khrt_children) then
          call advance_mpi_khrt_shedding(comm,components,options,trial,next_mass,primitive(1,:), &
            primitive(2:4,:),reynolds,dt,trial_id,child_energy,valid)
          if (valid) valid=size(trial)<=125000
          call agree()
          if (agreed) change(5,:)=change(5,:)+child_energy
        end if
        if (agreed) then
          next_gas=gas;next_temp=temp
          call deposit_mpi_spray_gas(comm,species,next_gas,next_temp,positions,change,lower,spacing,cells,parts,faces,valid)
          if (valid) then
            accepted=.true.;exit
          end if
        end if
        dt=.5_dp*dt
      end do
      if (.not.accepted) return
      drops=trial;masses=next_mass;gas=next_gas;temp=next_temp;accepted_id=trial_id
      elapsed=elapsed+dt;count_steps=count_steps+1
    end do
    np=size(drops)
    kept=count(drops%mass>0);allocate(survivors(kept),kept_mass(nc,kept));kept=0
    do p=1,np
      if (drops(p)%mass<=0) cycle
      kept=kept+1;survivors(kept)=drops(p);kept_mass(:,kept)=masses(:,p)
    end do
    call migrate_spray_parcels(comm,survivors,lower,spacing,cells,parts,[.true.,.true.,.true.], &
      options%maximum_parcels,valid,kept_mass)
    if (.not.valid) return
    state=gas;temperature=temp
    call move_alloc(survivors,parcels);call move_alloc(kept_mass,component_mass)
    if (present(last_id)) last_id=accepted_id
    substeps=count_steps;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine deposit_mpi_spray_gas(comm,species,state,temperature,positions,contributions, &
      lower,spacing,cells,parts,faces,ok)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: positions(:,:),contributions(:,:),lower(3),spacing(3)
    integer, intent(in) :: cells(3),parts(3)
    character(len=*), intent(in) :: faces(6)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:,:,:,:),recovered(:,:,:)
    integer :: extent(3),ierr
    logical :: valid,agreed
    ! Keep both gas arrays private until every owner accepts its thermodynamics.
    ! The caller owns the identical ordered species table and particle transaction.
    ok=.false.;extent=shape(temperature)
    allocate(recovered(extent(1),extent(2),extent(3)))
    call recover_reactive_temperatures_3d(species,state,temperature,extent(1),extent(2),extent(3),recovered,valid)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (.not.agreed) return
    candidate=state
    call deposit_mpi_spray_field(comm,candidate,positions,contributions,lower,spacing,cells,parts,faces,valid)
    if (.not.valid) return
    call recover_reactive_temperatures_3d(species,candidate,temperature,extent(1),extent(2),extent(3),recovered,valid)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (.not.agreed) return
    state=candidate;temperature=recovered;ok=.true.
  end subroutine

  subroutine sample_mpi_spray_field(comm,field,positions,lower,spacing,cells,parts,faces,samples,ok)
    type(MPI_Comm), intent(in) :: comm
    real(dp), intent(in) :: field(:,:,:,:),positions(:,:),lower(3),spacing(3)
    integer, intent(in) :: cells(3),parts(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), allocatable, intent(inout) :: samples(:,:)
    logical, intent(out) :: ok
    integer :: rank,nranks,ierr,nvar,np,nrequests,nreceived,p,n,r,q,i,cell(3),first(3),local_cells(3)
    integer :: config(13),minimum(13),maximum(13),indices(3,8)
    integer, allocatable :: owners(:),requests(:,:),incoming(:,:),counts(:),offsets(:),received(:),displacements(:),cursor(:)
    integer, allocatable :: order(:)
    real(dp), allocatable :: weights(:),packed(:,:),answers(:,:),candidate(:,:)
    real(dp) :: w(8),bounds(6),minbounds(6),maxbounds(6),center(3)
    logical :: valid,agreed,periodic(3)
    ok=.false.;np=size(positions,2);nvar=size(field,1)
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=size(positions,1)==3.and.np<=125000.and.nvar>=1.and.nvar<=256
    valid=valid.and.valid_spray_boundary_faces(faces).and.all(ieee_is_finite([lower,spacing]))
    valid=valid.and.all(spacing>=1e-12_dp).and.all(spacing<=1e9_dp).and.maxval(abs(lower))<=1e9_dp
    valid=valid.and.all(ieee_is_finite(positions)).and.all(ieee_is_finite(field))
    call agree()
    if (.not.agreed) return
    config(1:6)=[cells,parts];config(13)=nvar
    do i=1,6
      select case(trim(faces(i)))
      case('periodic');config(6+i)=1
      case('inflow');config(6+i)=2
      case('outflow');config(6+i)=3
      case('slip_wall');config(6+i)=4
      case('no_slip_wall');config(6+i)=5
      case('isothermal_wall');config(6+i)=6
      end select
    end do
    bounds=[lower,spacing]
    call MPI_Allreduce(config,minimum,13,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,maximum,13,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,minbounds,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,maxbounds,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call spray_cartesian_block(cells,parts,rank,first,local_cells,valid)
    valid=valid.and.all(minimum==maximum).and.all(minbounds==maxbounds)
    if (valid) valid=product(parts)==nranks.and.all(shape(field(1,:,:,:))==local_cells).and.all(cells>=2)
    call agree()
    if (.not.agreed) return
    periodic=[faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']
    nrequests=8*np
    allocate(owners(nrequests),requests(3,nrequests),weights(nrequests),order(nrequests))
    allocate(counts(nranks),offsets(nranks),received(nranks),displacements(nranks),cursor(nranks));counts=0
    do p=1,np
      call particle_boundary_stencil(positions(:,p),lower,spacing,cells,faces,indices,w,valid)
      if (.not.valid) exit
      do n=1,8
        q=8*(p-1)+n;center=lower+(real(indices(:,n),dp)-.5_dp)*spacing
        r=spray_cartesian_owner(center,lower,spacing,cells,parts,periodic)
        if (r<0) then
          valid=.false.;exit
        end if
        owners(q)=r+1;weights(q)=w(n);counts(r+1)=counts(r+1)+1
      end do
      if (.not.valid) exit
    end do
    call agree()
    if (.not.agreed) return
    offsets(1)=0
    do r=2,nranks
      offsets(r)=offsets(r-1)+counts(r-1)
    end do
    cursor=offsets
    do p=1,np
      call particle_boundary_stencil(positions(:,p),lower,spacing,cells,faces,indices,w,valid)
      do n=1,8
        q=8*(p-1)+n;r=owners(q);cursor(r)=cursor(r)+1
        requests(:,cursor(r))=indices(:,n);order(q)=cursor(r)
      end do
    end do
    call MPI_Alltoall(counts,1,MPI_INTEGER,received,1,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=sum(int(received,int64))<=1000000_int64
    call agree()
    if (.not.agreed) return
    nreceived=sum(received);displacements(1)=0
    do r=2,nranks
      displacements(r)=displacements(r-1)+received(r-1)
    end do
    allocate(incoming(3,nreceived),packed(nvar,nreceived),answers(nvar,nrequests),candidate(nvar,np))
    call MPI_Alltoallv(requests,3*counts,3*offsets,MPI_INTEGER, &
      incoming,3*received,3*displacements,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    do q=1,nreceived
      cell=incoming(:,q)-first
      if (any(cell<1).or.any(cell>local_cells)) then
        valid=.false.;exit
      end if
      packed(:,q)=field(:,cell(1),cell(2),cell(3))
    end do
    call agree()
    if (.not.agreed) return
    call MPI_Alltoallv(packed,nvar*received,nvar*displacements,MPI_DOUBLE_PRECISION, &
      answers,nvar*counts,nvar*offsets,MPI_DOUBLE_PRECISION,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    candidate=0
    do p=1,np
      do n=1,8
        q=8*(p-1)+n;candidate(:,p)=candidate(:,p)+weights(q)*answers(:,order(q))
      end do
    end do
    valid=all(ieee_is_finite(candidate))
    call agree()
    if (.not.agreed) return
    call move_alloc(candidate,samples);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine deposit_mpi_spray_field(comm,field,positions,contributions,lower,spacing,cells,parts,faces,ok)
    type(MPI_Comm), intent(in) :: comm
    real(dp), intent(inout) :: field(:,:,:,:)
    real(dp), intent(in) :: positions(:,:),contributions(:,:),lower(3),spacing(3)
    integer, intent(in) :: cells(3),parts(3)
    character(len=*), intent(in) :: faces(6)
    logical, intent(out) :: ok
    integer :: rank,nranks,ierr,nvar,np,nrequests,nreceived,p,n,r,q,i,v,cell(3),first(3),local_cells(3)
    integer :: config(13),minimum(13),maximum(13),indices(3,8)
    integer, allocatable :: owners(:),requests(:,:),incoming(:,:),counts(:),offsets(:),received(:),displacements(:),cursor(:)
    integer, allocatable :: order(:)
    real(dp), allocatable :: weights(:),packed(:,:),answers(:,:),candidate(:,:,:,:)
    real(dp) :: w(8),bounds(6),minbounds(6),maxbounds(6),center(3),volume,old,delta
    logical :: valid,agreed,periodic(3)
    ok=.false.;np=size(positions,2);nvar=size(field,1)
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=size(positions,1)==3.and.np<=125000.and.nvar>=1.and.nvar<=256
    valid=valid.and.valid_spray_boundary_faces(faces).and.all(ieee_is_finite([lower,spacing]))
    valid=valid.and.all(spacing>=1e-12_dp).and.all(spacing<=1e9_dp).and.maxval(abs(lower))<=1e9_dp
    valid=valid.and.all(ieee_is_finite(positions)).and.all(ieee_is_finite(field))
    valid=valid.and.all(shape(contributions)==[nvar,np]).and.all(ieee_is_finite(contributions))
    call agree()
    if (.not.agreed) return
    config(1:6)=[cells,parts];config(13)=nvar
    do i=1,6
      select case(trim(faces(i)))
      case('periodic');config(6+i)=1
      case('inflow');config(6+i)=2
      case('outflow');config(6+i)=3
      case('slip_wall');config(6+i)=4
      case('no_slip_wall');config(6+i)=5
      case('isothermal_wall');config(6+i)=6
      end select
    end do
    bounds=[lower,spacing]
    call MPI_Allreduce(config,minimum,13,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,maximum,13,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,minbounds,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,maxbounds,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call spray_cartesian_block(cells,parts,rank,first,local_cells,valid)
    valid=valid.and.all(minimum==maximum).and.all(minbounds==maxbounds)
    if (valid) valid=product(parts)==nranks.and.all(shape(field(1,:,:,:))==local_cells).and.all(cells>=2)
    call agree()
    if (.not.agreed) return
    volume=product(spacing)
    valid=ieee_is_finite(volume).and.volume>0
    if (valid.and.volume<1) valid=all(abs(contributions)<=huge(1.0_dp)*volume)
    call agree()
    if (.not.agreed) return
    periodic=[faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']
    nrequests=8*np
    allocate(owners(nrequests),requests(3,nrequests),weights(nrequests),order(nrequests))
    allocate(counts(nranks),offsets(nranks),received(nranks),displacements(nranks),cursor(nranks));counts=0
    do p=1,np
      call particle_boundary_stencil(positions(:,p),lower,spacing,cells,faces,indices,w,valid)
      if (.not.valid) exit
      do n=1,8
        q=8*(p-1)+n;center=lower+(real(indices(:,n),dp)-.5_dp)*spacing
        r=spray_cartesian_owner(center,lower,spacing,cells,parts,periodic)
        if (r<0) then
          valid=.false.;exit
        end if
        owners(q)=r+1;weights(q)=w(n);counts(r+1)=counts(r+1)+1
      end do
      if (.not.valid) exit
    end do
    call agree()
    if (.not.agreed) return
    offsets(1)=0
    do r=2,nranks
      offsets(r)=offsets(r-1)+counts(r-1)
    end do
    cursor=offsets
    do p=1,np
      call particle_boundary_stencil(positions(:,p),lower,spacing,cells,faces,indices,w,valid)
      do n=1,8
        q=8*(p-1)+n;r=owners(q);cursor(r)=cursor(r)+1
        requests(:,cursor(r))=indices(:,n);order(q)=cursor(r)
      end do
    end do
    call MPI_Alltoall(counts,1,MPI_INTEGER,received,1,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=sum(int(received,int64))<=1000000_int64
    call agree()
    if (.not.agreed) return
    nreceived=sum(received);displacements(1)=0
    do r=2,nranks
      displacements(r)=displacements(r-1)+received(r-1)
    end do
    allocate(incoming(3,nreceived),packed(nvar,nrequests),answers(nvar,nreceived))
    allocate(candidate,source=field)
    call MPI_Alltoallv(requests,3*counts,3*offsets,MPI_INTEGER, &
      incoming,3*received,3*displacements,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    do p=1,np
      do n=1,8
        q=8*(p-1)+n;packed(:,order(q))=weights(q)*contributions(:,p)/volume
      end do
    end do
    call MPI_Alltoallv(packed,nvar*counts,nvar*offsets,MPI_DOUBLE_PRECISION, &
      answers,nvar*received,nvar*displacements,MPI_DOUBLE_PRECISION,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    do q=1,nreceived
      cell=incoming(:,q)-first
      if (any(cell<1).or.any(cell>local_cells)) then
        valid=.false.;exit
      end if
      do v=1,nvar
        old=candidate(v,cell(1),cell(2),cell(3));delta=answers(v,q)
        if (.not.ieee_is_finite(delta)) then
          valid=.false.;exit
        end if
        if ((old>0.and.delta>0).or.(old<0.and.delta<0)) then
          if (abs(old)>huge(1.0_dp)-abs(delta)) then
            valid=.false.;exit
          end if
        end if
        candidate(v,cell(1),cell(2),cell(3))=old+delta
      end do
      if (.not.valid) exit
    end do
    valid=valid.and.all(ieee_is_finite(candidate))
    call agree()
    if (.not.agreed) return
    field=candidate;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine inject_mpi_cone_parcels(comm,liquid,total_mass,number,center,axis,speed,half_angle,diameter,temperature, &
      maximum_parcels,last_id,parcels,injected_totals,ok,components,fractions,component_mass, &
      diameter_shape,diameter_min,diameter_max,diameter_weighting,diameter_distribution,diameter_std,nozzle_radius,nozzle_model)
    type(MPI_Comm), intent(in) :: comm
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: total_mass,center(3),axis(3),speed,half_angle,diameter,temperature
    integer, intent(in) :: number,maximum_parcels
    integer(int64), intent(inout) :: last_id
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(out) :: injected_totals(:)
    logical, intent(out) :: ok
    type(spray_liquid), optional, intent(in) :: components(:)
    real(dp), optional, intent(in) :: fractions(:)
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    real(dp), optional, intent(in) :: diameter_shape,diameter_min,diameter_max,diameter_std,nozzle_radius
    character(len=*), optional, intent(in) :: diameter_weighting,diameter_distribution
    type(spray_nozzle_model), optional, intent(in) :: nozzle_model
    type(spray_parcel), allocatable :: candidate(:),born(:)
    type(spray_liquid) :: material
    real(dp), allocatable :: masses(:,:)
    real(dp) :: rho,cp,e,ledger(size(injected_totals)),inventory(5)
    integer :: p,c,nc,lo,hi,old_count,ierr
    integer(int64) :: next_id,first_id
    logical :: valid,agreed,mixed
    ok=.false.;injected_totals=0;nc=0;mixed=present(components)
    valid=allocated(parcels).and.maximum_parcels>=1.and.maximum_parcels<=1000000
    valid=valid.and.number>=0.and.number<=1000000.and.last_id>=0.and.size(injected_totals)>=6
    valid=valid.and.(mixed.eqv.present(fractions)).and.(mixed.eqv.present(component_mass))
    if (present(component_mass)) valid=valid.and.allocated(component_mass)
    if (.not.ieee_is_finite(total_mass)) valid=.false.
    if (number==0.and.total_mass/=0) valid=.false.
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    old_count=size(parcels);valid=number<=maximum_parcels-old_count
    material=liquid
    if (mixed) then
      nc=size(components)
      valid=valid.and.nc>=1.and.nc<=16.and.size(fractions)==nc
      valid=valid.and.all(shape(component_mass)==[nc,old_count])
      if (valid) then
        call liquid_mixture_properties(components,fractions,temperature,rho,cp,e,valid)
        if (valid) then
          material=components(1);material%density=rho;material%heat_capacity=cp
          material%reference_temperature=temperature;material%reference_energy=e
          material%density_model='constant';material%density_coefficients=0;material%density_temperature_range=0
          material%heat_capacity_model='constant';material%heat_capacity_coefficients=0
          material%heat_capacity_temperature_range=0;material%boiling_model='reject';material%critical_temperature=0
        end if
      end if
    else
      valid=valid.and.valid_spray_liquid(material)
    end if
    if (valid) then
      do p=1,old_count
        valid=valid.and.valid_spray_parcel(parcels(p)).and.parcels(p)%id<=last_id
        if (p>1) valid=valid.and.parcels(p)%id>parcels(p-1)%id
        if (mixed) valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p))
      end do
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(nc,lo,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nc,hi,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.lo/=hi) return
    next_id=last_id;first_id=0
    call reserve_spray_ids(comm,number,next_id,first_id,valid)
    if (.not.valid) return
    allocate(born(0));ledger=0
    if (number>0) call inject_cone_parcels(material,total_mass,number,first_id,center,axis,speed,half_angle,diameter, &
      temperature,born,valid,diameter_shape,diameter_min,diameter_max,diameter_weighting,diameter_distribution, &
      diameter_std,nozzle_radius,nozzle_model)
    if (valid) then
      do p=1,size(born)
        inventory=parcel_totals(born(p),material);ledger(1:5)=ledger(1:5)+inventory
        if (mixed) then
          do c=1,nc
            if (components(c)%vapor_index>size(ledger)-5) then
              valid=.false.;exit
            end if
            ledger(5+components(c)%vapor_index)=ledger(5+components(c)%vapor_index)+inventory(1)*fractions(c)
          end do
        else
          if (material%vapor_index>size(ledger)-5) then
            valid=.false.;exit
          end if
          ledger(5+material%vapor_index)=ledger(5+material%vapor_index)+inventory(1)
        end if
      end do
    end if
    valid=valid.and.all(ieee_is_finite(ledger))
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    candidate=[parcels,born]
    if (mixed) then
      allocate(masses(nc,size(candidate)));masses(:,:old_count)=component_mass
      do p=1,size(born)
        masses(:,old_count+p)=born(p)%mass*fractions
      end do
      call move_alloc(masses,component_mass)
    end if
    call move_alloc(candidate,parcels);last_id=next_id;injected_totals=ledger;ok=.true.
  end subroutine

  subroutine split_mpi_spray_parcels(comm,parcels,maximum_multiplicity,maximum_parcels,last_id,ok,component_mass)
    type(MPI_Comm), intent(in) :: comm
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(in) :: maximum_multiplicity
    integer, intent(in) :: maximum_parcels
    integer(int64), intent(inout) :: last_id
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: masses(:,:)
    integer :: p,parent,old_count,nc,lo,hi,ierr
    integer(int64) :: next_id,first_id
    logical :: valid,agreed
    ok=.false.;valid=allocated(parcels).and.last_id>=0;nc=0
    if (present(component_mass)) then
      valid=valid.and.allocated(component_mass)
      if (allocated(component_mass)) nc=size(component_mass,1)
      valid=valid.and.nc>=1.and.nc<=16
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(nc,lo,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nc,hi,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.lo/=hi) return
    old_count=size(parcels);candidate=parcels
    if (present(component_mass)) valid=size(component_mass,2)==old_count
    if (valid) then
      do p=1,old_count
        valid=valid.and.valid_spray_parcel(parcels(p)).and.parcels(p)%id<=last_id
        if (p>1) valid=valid.and.parcels(p)%id>parcels(p-1)%id
        if (present(component_mass)) valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p))
      end do
    end if
    if (valid) call split_spray_parcels(candidate,maximum_multiplicity,maximum_parcels,valid)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    if (present(component_mass)) then
      allocate(masses(nc,size(candidate)));masses(:,:old_count)=component_mass;parent=1
      do p=old_count+1,size(candidate)
        do while(parent<old_count)
          if (parcels(parent)%id>=candidate(p)%parent_id) exit
          parent=parent+1
        end do
        masses(:,p)=component_mass(:,parent)
      end do
    end if
    next_id=last_id;first_id=0
    call reserve_spray_ids(comm,size(candidate)-old_count,next_id,first_id,valid)
    if (.not.valid) return
    do p=old_count+1,size(candidate)
      candidate(p)%id=first_id+int(p-old_count-1,int64)
    end do
    call move_alloc(candidate,parcels)
    if (present(component_mass)) call move_alloc(masses,component_mass)
    last_id=next_id;ok=.true.
  end subroutine

  ! All ranks retain the same high-water mark; rank order breaks allocation ties.
  subroutine reserve_spray_ids(comm,local_count,last_id,first_id,ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: local_count
    integer(int64), intent(inout) :: last_id,first_id
    logical, intent(out) :: ok
    integer(int64) :: minimum,maximum,requested,total,prefix,next_first
    integer :: rank,ierr
    logical :: valid,agreed
    ok=.false.;valid=local_count>=0.and.local_count<=1000000.and.last_id>=0
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(last_id,minimum,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(last_id,maximum,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.minimum/=maximum) return
    requested=int(local_count,int64)
    call MPI_Allreduce(requested,total,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (total>huge(last_id)-last_id) return
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Exscan(requested,prefix,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (rank==0) prefix=0
    next_first=0
    if (local_count>0) next_first=last_id+prefix+1_int64
    first_id=next_first;last_id=last_id+total;ok=.true.
  end subroutine

  subroutine migrate_spray_parcels(comm,parcels,lower,spacing,cells,parts,periodic,maximum_parcels,ok,component_mass,gather_rank)
    integer, intent(in), optional :: gather_rank
    type(MPI_Comm), intent(in) :: comm
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(in) :: lower(3),spacing(3)
    integer, intent(in) :: cells(3),parts(3),maximum_parcels
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    type(spray_parcel), allocatable :: outgoing(:),incoming(:)
    real(dp), allocatable :: masses(:,:),send_real(:,:),recv_real(:,:)
    integer(int64), allocatable :: send_int(:,:),recv_int(:,:)
    integer, allocatable :: counts(:),offsets(:),received(:),displacements(:)
    integer :: rank,nranks,ierr,nc,nreal,nrecv,p,r,config(11),lo(11),hi(11),destination
    real(dp) :: bounds(6),minbounds(6),maxbounds(6)
    logical :: valid,global_valid
    ok=.false.
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    destination=-1
    if (present(gather_rank)) destination=gather_rank
    valid=allocated(parcels).and.maximum_parcels>=0.and.maximum_parcels<=1000000;nc=0
    if (allocated(parcels)) valid=valid.and.size(parcels)<=1000000
    if (present(component_mass)) then
      valid=valid.and.allocated(component_mass)
      if (allocated(component_mass)) nc=size(component_mass,1)
      valid=valid.and.nc>=1.and.nc<=16
    end if
    valid=valid.and.destination>=-1.and.destination<nranks
    valid=valid.and.all(ieee_is_finite([lower,spacing]))
    call agree()
    if (.not.global_valid) return
    config=[cells,parts,merge(1,0,periodic),nc,destination];bounds=[lower,spacing]
    call MPI_Allreduce(config,lo,11,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,hi,11,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,minbounds,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,maxbounds,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(lo==hi).and.all(minbounds==maxbounds)
    valid=valid.and.spray_cartesian_owner(lower,lower,spacing,cells,parts,periodic)>=0
    call agree()
    if (.not.global_valid) return
    valid=product(parts)==nranks
    outgoing=parcels
    if (valid) then
      if (present(component_mass)) then
        masses=component_mass
        call partition_spray_parcels(outgoing,lower,spacing,cells,parts,periodic,counts,offsets,valid,masses)
      else
        call partition_spray_parcels(outgoing,lower,spacing,cells,parts,periodic,counts,offsets,valid)
      end if
    end if
    call agree()
    if (.not.global_valid) return
    if (destination>=0) then
      counts=0;counts(destination+1)=size(outgoing);offsets=0
      do r=2,nranks
        offsets(r)=offsets(r-1)+counts(r-1)
      end do
    end if
    allocate(received(nranks),displacements(nranks))
    call MPI_Alltoall(counts,1,MPI_INTEGER,received,1,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=sum(int(received,int64))<=int(maximum_parcels,int64)
    call agree()
    if (.not.global_valid) return
    nrecv=sum(received);displacements(1)=0
    do r=2,nranks
      displacements(r)=displacements(r-1)+received(r-1)
    end do
    nreal=14+nc
    allocate(send_int(3,size(outgoing)),send_real(nreal,size(outgoing)),recv_int(3,nrecv),recv_real(nreal,nrecv))
    do p=1,size(outgoing)
      send_int(:,p)=[outgoing(p)%id,outgoing(p)%parent_id,int(outgoing(p)%material_id,int64)]
      send_real(1:14,p)=[outgoing(p)%position,outgoing(p)%velocity,outgoing(p)%mass,outgoing(p)%temperature, &
        outgoing(p)%multiplicity,outgoing(p)%deformation,outgoing(p)%deformation_rate,outgoing(p)%rt_time, &
        outgoing(p)%khrt_initial_count,outgoing(p)%khrt_shed_volume]
      if (nc>0) send_real(15:,p)=masses(:,p)
    end do
    call MPI_Alltoallv(send_int,3*counts,3*offsets,MPI_INTEGER8,recv_int,3*received,3*displacements,MPI_INTEGER8,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Alltoallv(send_real,nreal*counts,nreal*offsets,MPI_DOUBLE_PRECISION, &
      recv_real,nreal*received,nreal*displacements,MPI_DOUBLE_PRECISION,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (nrecv>1) call sort_records(1,nrecv)
    allocate(incoming(nrecv))
    if (allocated(masses)) deallocate(masses)
    if (nc>0) allocate(masses(nc,nrecv))
    valid=.true.
    do p=1,nrecv
      incoming(p)=spray_parcel(recv_int(1,p),recv_real(1:3,p),recv_real(4:6,p),recv_real(7,p), &
        recv_real(8,p),recv_real(9,p),recv_real(10,p),recv_real(11,p),recv_int(2,p),recv_real(12,p), &
        int(recv_int(3,p)),recv_real(13,p),recv_real(14,p))
      valid=valid.and.valid_spray_parcel(incoming(p))
      if (p>1) valid=valid.and.incoming(p)%id>incoming(p-1)%id
      if (destination>=0) then
        valid=valid.and.rank==destination
      else
        valid=valid.and.spray_cartesian_owner(incoming(p)%position,lower,spacing,cells,parts,periodic)==rank
      end if
      if (nc>0) then
        masses(:,p)=recv_real(15:,p)
        valid=valid.and.valid_component_inventory(incoming(p),masses(:,p))
      end if
    end do
    call agree()
    if (.not.global_valid) return
    call move_alloc(incoming,parcels)
    if (present(component_mass)) call move_alloc(masses,component_mass)
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,global_valid,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) global_valid=.false.
    end subroutine
    recursive subroutine sort_records(first,last)
      integer, intent(in) :: first,last
      integer :: left,right
      integer(int64) :: pivot,tmp_int(3)
      real(dp) :: tmp_real(nreal)
      left=first;right=last;pivot=recv_int(1,first+(last-first)/2)
      do
        do while(recv_int(1,left)<pivot)
          left=left+1
        end do
        do while(recv_int(1,right)>pivot)
          right=right-1
        end do
        if (left<=right) then
          tmp_int=recv_int(:,left);recv_int(:,left)=recv_int(:,right);recv_int(:,right)=tmp_int
          tmp_real=recv_real(:,left);recv_real(:,left)=recv_real(:,right);recv_real(:,right)=tmp_real
          left=left+1;right=right-1
        end if
        if (left>right) exit
      end do
      if (first<right) call sort_records(first,right)
      if (left<last) call sort_records(left,last)
    end subroutine
  end subroutine
end module
