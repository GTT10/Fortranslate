module spray_mpi_mod
  use mpi_f08
  use spray_mpi_exchange_mod, only: exchange_spray_parcel_bins
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d, &
    valid_reactive_boundaries_3d,valid_wall_velocity_3d,fill_pressure_outlet_ghosts_3d
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use spray_reactive_3d_mod, only: advance_spray_chemistry_3d
  use reactive_transport_3d_mod, only: reactive_transport_ghosted_fluxes_3d,physical_transport_ghosted_fluxes_3d
  use les_smagorinsky_3d_mod, only: compute_les_ghosted_rhs_3d
  use les_mpi_boundary_mod, only: fill_mpi_les_boundary_ghosts_3d
  use les_mpi_dynamic_mod, only: advance_mpi_dynamic_les_ssprk2_3d
  use les_filter_3d_mod, only: les_filter_weights,les_filter_code
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_conserved_to_primitive
  use reactive_3d_mod, only: recover_reactive_temperatures_3d,compute_reactive_plm_slab_face_fluxes_3d, &
    compute_reactive_plm_face_fluxes_3d,compute_reactive_ppm_halo_ctu_fluxes_3d
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory,liquid_mixture_properties, &
    advance_mixture_frozen_gas,mixture_transfer_rates,mixture_trace_fraction,build_mixture_film
  use spray_breakup_mod, only: split_spray_parcels,advance_spray_breakup,advance_khrt_shedding, &
    construct_khrt_children,spray_khrt_release,spray_khrt_shedding_state
  use spray_coupling_3d_mod, only: particle_boundary_stencil,valid_spray_boundary_faces,apply_parcel_boundaries
  use spray_wall_impact_mod, only: spray_wall_film,valid_wall_component_inventory,spray_wall_options,valid_wall_film, &
    valid_wall_options,wall_liquid_conductivity_at_temperature,mixture_wall_film_rates,advance_mixture_wall_film_exchange, &
    classify_wall_impact,resolve_mixture_wall_impact,mixture_impact_boiling_temperature
  use spray_partition_mod
  implicit none
  private
  public :: migrate_spray_parcels,reserve_spray_ids,split_mpi_spray_parcels,inject_mpi_cone_parcels
  public :: advance_mpi_khrt_shedding
  public :: fill_mpi_spray_ghosts
  public :: impact_mpi_planar_walls
  public :: sample_mpi_wall_film_thickness
  public :: resolve_mpi_wall_impacts
  public :: advance_mpi_film_coupling
  public :: migrate_mpi_wall_films
  public :: advance_mpi_spray_les
  public :: advance_mpi_spray_transport
  public :: advance_mpi_spray_hydro
  public :: advance_mpi_mixture_coupling,advance_mpi_mixture_reaction
  public :: sample_mpi_spray_field,deposit_mpi_spray_field,deposit_mpi_spray_gas
contains
  ! Validate and agree controls before periodic/no-op dispatch in coupled entries.
  subroutine validate_mpi_mach_controls(comm,inflow_mach,inflow_pressure,faces,ok)
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: faces(6)
    logical, intent(out) :: ok
    real(dp) :: mach(6),pressure(6),controls(12),minimum(12),maximum(12)
    character(len=32) :: resolved(6)
    logical :: valid
    integer :: ierr
    ok=.false.;mach=-1;pressure=0;resolved='periodic'
    if(present(inflow_mach)) mach=inflow_mach
    if(present(inflow_pressure)) pressure=inflow_pressure
    if(present(faces)) resolved=faces
    valid=all(ieee_is_finite(mach)).and.all(mach<=100)
    valid=valid.and..not.any(mach>=0.and.resolved/='inflow')
    valid=valid.and.all(ieee_is_finite(pressure)).and.all(pressure>=0).and.all(pressure<=1e15_dp)
    valid=valid.and..not.any(pressure>0.and.mach<0)
    if(present(inflow_pressure)) valid=valid.and.present(inflow_mach)
    call MPI_Allreduce(valid,ok,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS) ok=.false.
    if(.not.ok) return
    ok=.false.;controls=[mach,pressure]
    call MPI_Allreduce(controls,minimum,12,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,12,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum/=maximum)) return
    ok=.true.
  end subroutine validate_mpi_mach_controls

  subroutine fill_mpi_spray_ghosts(comm,state,temperature,cells,faces,ghost,ghost_temperature,ok, &
      viscous_wall,inflow_state,inflow_temperature,wall_velocity,species,outlet_pressure,inflow_mach,inflow_pressure)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), optional, intent(in) :: species(:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp) :: pressure(6),local_pressure(6),minpressure(6),maxpressure(6)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:)
    integer, intent(in) :: cells(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), allocatable, intent(inout) :: ghost(:,:,:,:),ghost_temperature(:,:,:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: viscous_wall
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_velocity(:,:)
    real(dp), allocatable :: candidate(:,:,:,:),temp(:,:,:),send(:,:,:,:),left(:,:,:,:),right(:,:,:,:)
    real(dp) :: local_velocity(3,6),mach(6),pin(6),mach_controls(12),mach_min(12),mach_max(12)
    integer :: rank,nranks,ierr,first(3),extent(3),nx,ny,nz,nvar,previous,next,f,config(11),low(11),high(11)
    character(len=32) :: local_faces(6)
    logical :: valid,agreed,viscous
    ok=.false.;viscous=.false.
    if (present(viscous_wall)) viscous=viscous_wall
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    mach=-1;pin=0
    if(present(inflow_mach)) mach=inflow_mach
    if(present(inflow_pressure)) pin=inflow_pressure
    pressure=0
    if (present(outlet_pressure)) pressure=outlet_pressure
    nvar=size(state,1);config(:4)=[cells,nvar];config(11)=merge(1,0,viscous)
    valid=valid_reactive_boundaries_3d(faces).and.nvar>=5.and.nvar<=256
    valid=valid.and.all(ieee_is_finite(pressure))
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    valid=valid.and..not.any(pressure>0.and.faces/='outflow')
    if (any(pressure>0)) then
      valid=valid.and.present(species)
      if (present(species)) valid=valid.and.nvar==5+size(species)
    end if
    valid=valid.and.all(ieee_is_finite(mach)).and.all(mach<=100)
    valid=valid.and..not.any(mach>=0.and.faces/='inflow')
    valid=valid.and.all(ieee_is_finite(pin)).and.all(pin>=0).and.all(pin<=1e15_dp)
    valid=valid.and..not.any(pin>0.and.mach<0)
    if(present(inflow_pressure)) valid=valid.and.present(inflow_mach)
    if(any(mach>=0)) then
      valid=valid.and.present(species)
      if(present(species)) valid=valid.and.nvar==5+size(species)
    end if
    do f=1,6
      select case(trim(faces(f)))
      case('periodic');config(4+f)=1
      case('inflow');config(4+f)=2
      case('outflow');config(4+f)=3
      case('slip_wall');config(4+f)=4
      case('no_slip_wall');config(4+f)=5
      case('isothermal_wall');config(4+f)=6
      case default;valid=.false.;config(4+f)=0
      end select
    end do
    local_velocity=0
    if (present(wall_velocity)) then
      valid=valid.and.valid_wall_velocity_3d(faces,wall_velocity)
      if (valid) local_velocity=wall_velocity
    end if
    call agree()
    if (.not.agreed) return
    call MPI_Allreduce(pressure,minpressure,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(pressure,maxpressure,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(minpressure/=maxpressure)) return
    call MPI_Allreduce(config,low,11,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,high,11,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.any(low/=high)) return
    mach_controls=[mach,pin]
    call MPI_Allreduce(mach_controls,mach_min,12,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(mach_controls,mach_max,12,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(mach_min/=mach_max)) return
    if(any(mach>=0)) then
      call fill_mpi_les_boundary_ghosts_3d(comm,species,state,temperature,cells,faces,2,ghost,ghost_temperature,ok, &
        inflow_state,inflow_temperature,wall_velocity,outlet_pressure,mach,pin,viscous_wall=viscous)
      return
    end if
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    valid=valid.and.all(extent>=2).and.all(shape(state)==[nvar,extent]).and.all(shape(temperature)==extent)
    call agree()
    if (.not.agreed) return
    nx=extent(1);ny=extent(2);nz=extent(3)
    valid=2*(nvar+1)<=huge(1)/(ny+4)/(nz+4)
    call agree()
    if (.not.agreed) return
    local_faces=faces;local_faces(1:2)='outflow';local_velocity(:,1:2)=0
    previous=modulo(rank-1,nranks);next=modulo(rank+1,nranks)
    if (rank==0.and.faces(1)/='periodic') then
      previous=MPI_PROC_NULL;local_faces(1)=faces(1)
      if (present(wall_velocity)) local_velocity(:,1)=wall_velocity(:,1)
    end if
    if (rank==nranks-1.and.faces(2)/='periodic') then
      next=MPI_PROC_NULL;local_faces(2)=faces(2)
      if (present(wall_velocity)) local_velocity(:,2)=wall_velocity(:,2)
    end if
    allocate(candidate(nvar,nx+4,ny+4,nz+4),temp(nx+4,ny+4,nz+4))
    call fill_reactive_boundary_ghosts_3d(state,temperature,local_faces,candidate,temp,valid,viscous_wall=viscous, &
      inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_velocity=local_velocity)
    call agree()
    if (.not.agreed) return
    if (any(pressure>0)) then
      local_pressure=pressure
      if (rank/=0) local_pressure(1)=0
      if (rank/=nranks-1) local_pressure(2)=0
      call fill_pressure_outlet_ghosts_3d(species,extent,local_faces,local_pressure,candidate,temp,valid)
      call agree()
      if (.not.agreed) return
    end if
    allocate(send(nvar+1,2,ny+4,nz+4),left(nvar+1,2,ny+4,nz+4),right(nvar+1,2,ny+4,nz+4))
    send(:nvar,:,:,:)=candidate(:,nx+1:nx+2,:,:);send(nvar+1,:,:,:)=temp(nx+1:nx+2,:,:)
    call MPI_Sendrecv(send,size(send),MPI_DOUBLE_PRECISION,next,3331,left,size(left),MPI_DOUBLE_PRECISION, &
      previous,3331,comm,MPI_STATUS_IGNORE,ierr)
    if (ierr/=MPI_SUCCESS) return
    send(:nvar,:,:,:)=candidate(:,3:4,:,:);send(nvar+1,:,:,:)=temp(3:4,:,:)
    call MPI_Sendrecv(send,size(send),MPI_DOUBLE_PRECISION,previous,3332,right,size(right),MPI_DOUBLE_PRECISION, &
      next,3332,comm,MPI_STATUS_IGNORE,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (previous/=MPI_PROC_NULL) then
      candidate(:,:2,:,:)=left(:nvar,:,:,:);temp(:2,:,:)=left(nvar+1,:,:,:)
    end if
    if (next/=MPI_PROC_NULL) then
      candidate(:,nx+3:,:,:)=right(:nvar,:,:,:);temp(nx+3:,:,:)=right(nvar+1,:,:,:)
    end if
    call move_alloc(candidate,ghost);call move_alloc(temp,ghost_temperature);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine impact_mpi_planar_walls(comm,components,properties,settings,faces,molecular_weights,pressure, &
      previous,lower,upper,parcels,component_mass,films,film_mass,last_id,hit,wall_change,ok)
    type(MPI_Comm), intent(in) :: comm
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: properties
    type(spray_wall_options), intent(in) :: settings
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: molecular_weights(:),pressure(:),previous(:,:),lower(3),upper(3)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:),film_mass(:,:)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    integer(int64), intent(inout) :: last_id
    logical, allocatable, intent(inout) :: hit(:)
    real(dp), intent(out) :: wall_change(5)
    logical, intent(out) :: ok
    integer, allocatable :: indices(:),wall_faces(:)
    real(dp), allocatable :: location(:,:),normal(:,:),data(:,:),diameters(:),thickness(:)
    logical, allocatable :: crossed(:)
    real(dp) :: bounds(6),low(6),high(6),delta(3),earliest,fraction,boundary,rho,cp,e,boil,mu,un
    real(dp) :: y(size(components)),sigma
    integer :: np,nc,p,a,f,face,ne,n,ierr
    logical :: valid,agreed
    ok=.false.;wall_change=0
    valid=allocated(parcels).and.allocated(component_mass).and.allocated(films).and.allocated(film_mass)
    valid=valid.and.valid_spray_options(properties).and.valid_wall_options(settings,faces)
    valid=valid.and.all(ieee_is_finite([lower,upper])).and.all(upper>lower)
    valid=valid.and.maxval(abs([lower,upper]))<=1e9_dp
    call agree()
    if (.not.agreed) return
    np=size(parcels);nc=size(components)
    valid=np<=125000.and.nc>=1.and.nc<=16.and.size(molecular_weights)==nc.and.size(pressure)==np
    valid=valid.and.all(shape(previous)==[3,np]).and.all(shape(component_mass)==[nc,np])
    call agree()
    if (.not.agreed) return
    valid=all(ieee_is_finite(previous)).and.all(ieee_is_finite(pressure))
    valid=valid.and.all(pressure>=1).and.all(pressure<=1e10_dp)
    do p=1,np
      valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p))
      valid=valid.and.all(previous(:,p)>=lower).and.all(previous(:,p)<upper)
    end do
    call agree()
    if (.not.agreed) return
    bounds=[lower,upper]
    call MPI_Allreduce(bounds,low,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,high,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.any(low/=high)) return
    allocate(indices(np),wall_faces(np),location(3,np),normal(3,np),crossed(np));crossed=.false.;ne=0
    do p=1,np
      if (parcels(p)%mass<=0) cycle
      delta=parcels(p)%position-previous(:,p);earliest=2;face=0
      do a=1,3
        f=0;boundary=0
        if (parcels(p)%position(a)<lower(a)) then
          f=2*a-1;boundary=lower(a)
        else if (parcels(p)%position(a)>=upper(a)) then
          f=2*a;boundary=upper(a)
        end if
        if (f==0.or.delta(a)==0) cycle
        if (faces(f)=='periodic') cycle
        fraction=(boundary-previous(a,p))/delta(a)
        if (fraction<earliest) then
          earliest=fraction;face=f
        end if
      end do
      if (face==0) cycle
      if (.not.settings%enabled(face)) cycle
      ne=ne+1;indices(ne)=p;wall_faces(ne)=face;crossed(p)=.true.
      a=(face+1)/2;location(:,ne)=previous(:,p)+earliest*delta;normal(:,ne)=0
      location(a,ne)=merge(lower(a),upper(a),modulo(face,2)==1)
      normal(a,ne)=merge(1.0_dp,-1.0_dp,modulo(face,2)==1)
      do a=1,3
        if (faces(2*a)=='periodic') location(a,ne)=lower(a)+modulo(location(a,ne)-lower(a),upper(a)-lower(a))
      end do
    end do
    allocate(data(6,ne),diameters(ne));data=0;sigma=properties%surface_tension
    do n=1,ne
      p=indices(n);y=component_mass(:,p)/parcels(p)%mass
      call liquid_mixture_properties(components,y,parcels(p)%temperature,rho,cp,e,valid)
      if (.not.valid) exit
      call mixture_impact_boiling_temperature(components,y,molecular_weights,pressure(p),boil,valid)
      if (.not.valid) exit
      mu=liquid_viscosity_at_temperature(properties,parcels(p)%temperature)
      valid=mu>0.and.sigma>=1e-8_dp.and.sigma<=1e3_dp
      if (.not.valid) exit
      diameters(n)=(6*parcels(p)%mass/(acos(-1.0_dp)*rho))**(1.0_dp/3.0_dp)
      un=abs(dot_product(parcels(p)%velocity,normal(:,n)))
      data(5,n)=rho*diameters(n)*un**2/sigma;data(6,n)=un*diameters(n)*rho/mu
      data(1,n)=sqrt(data(5,n)*sqrt(data(6,n)));data(2,n)=settings%temperature(wall_faces(n))/boil
      data(3,n)=(26.788_dp*log(pressure(p)*1e-5_dp)+33.546_dp)* &
        exp(-.006_dp*(settings%temperature(wall_faces(n))-260))+36
    end do
    call agree()
    if (.not.agreed) return
    call sample_mpi_wall_film_thickness(comm,films,location(:,:ne),normal(:,:ne), &
      1e-10_dp*maxval(upper-lower),thickness,valid)
    if (.not.valid) return
    data(4,:)=thickness/diameters
    call resolve_mpi_wall_impacts(comm,components,parcels,component_mass,indices(:ne),location(:,:ne),normal(:,:ne),data, &
      settings%contact_angle,films,film_mass,properties%maximum_parcels,settings%maximum_films,last_id,wall_change,valid)
    if (.not.valid) return
    call move_alloc(crossed,hit);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  ! Query distributed circular stationary film footprints without replicating films.
  ! Bounded query batches visit every rank; rank-local film inventories stay local.
  subroutine sample_mpi_wall_film_thickness(comm,films,positions,normals,plane_tolerance,thickness,ok)
    type(MPI_Comm), intent(in) :: comm
    type(spray_wall_film), intent(in) :: films(:)
    real(dp), intent(in) :: positions(:,:),normals(:,:),plane_tolerance
    real(dp), allocatable, intent(inout) :: thickness(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: query(:,:),local(:),total(:),candidate(:)
    real(dp) :: delta(3),minimum,maximum,limit
    integer :: rank,nranks,ierr,source,nq,count_queries,q,p
    logical :: valid,agreed
    ok=.false.;nq=size(positions,2)
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=size(films)<=1000000.and.nq<=125000.and.size(positions,1)==3.and.all(shape(normals)==[3,nq])
    valid=valid.and.ieee_is_finite(plane_tolerance)
    valid=valid.and.plane_tolerance>=0.and.plane_tolerance<=1e9_dp
    call agree()
    if (.not.agreed) return
    valid=all(ieee_is_finite(positions)).and.all(ieee_is_finite(normals))
    valid=valid.and.maxval(abs(positions))<=1e9_dp.and.maxval(abs(normals))<=1
    do p=1,size(films)
      valid=valid.and.valid_wall_film(films(p))
    end do
    call agree()
    if (.not.agreed) return
    do q=1,nq
      valid=valid.and.abs(sum(normals(:,q)**2)-1)<1e-12_dp
    end do
    call agree()
    if (.not.agreed) return
    call MPI_Allreduce(plane_tolerance,minimum,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(plane_tolerance,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.minimum/=maximum) return
    allocate(candidate(nq));candidate=0;limit=huge(1.0_dp)/real(nranks+1,dp)
    do source=0,nranks-1
      count_queries=nq
      call MPI_Bcast(count_queries,1,MPI_INTEGER,source,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      if (count_queries==0) cycle
      allocate(query(6,count_queries),local(count_queries),total(count_queries));local=0
      if (rank==source) then
        query(:3,:)=positions;query(4:,:)=normals
      end if
      call MPI_Bcast(query,size(query),MPI_DOUBLE_PRECISION,source,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      valid=.true.
      do q=1,count_queries
        do p=1,size(films)
          if (maxval(abs(films(p)%normal-query(4:,q)))>1e-12_dp) cycle
          delta=query(:3,q)-films(p)%position
          if (abs(dot_product(delta,query(4:,q)))>plane_tolerance) cycle
          if (norm2(delta)>.5_dp*films(p)%diameter) cycle
          if (films(p)%height>limit-local(q)) then
            valid=.false.;exit
          end if
          local(q)=local(q)+films(p)%height
        end do
        if (.not.valid) exit
      end do
      call agree()
      if (.not.agreed) return
      call MPI_Reduce(local,total,count_queries,MPI_DOUBLE_PRECISION,MPI_SUM,source,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      if (rank==source) candidate=total
      deallocate(query,local,total)
    end do
    call move_alloc(candidate,thickness);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  ! Resolve supplied wall events. data rows: Kv,T/Tboil,Kcrit,film/d,We,Re.
  ! Geometry crossing and film-overlap queries remain the caller's responsibility.
  subroutine resolve_mpi_wall_impacts(comm,components,parcels,component_mass,indices,locations,normals,data, &
      contact_angle,films,film_mass,maximum_parcels,maximum_films,last_id,wall_change,ok)
    type(MPI_Comm), intent(in) :: comm
    type(spray_liquid), intent(in) :: components(:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:),film_mass(:,:)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    integer, intent(in) :: indices(:),maximum_parcels,maximum_films
    real(dp), intent(in) :: locations(:,:),normals(:,:),data(:,:),contact_angle
    integer(int64), intent(inout) :: last_id
    real(dp), intent(out) :: wall_change(5)
    logical, intent(out) :: ok
    type(spray_parcel), allocatable :: candidate(:),children(:)
    type(spray_wall_film), allocatable :: patches(:)
    type(spray_wall_film) :: deposited
    real(dp), allocatable :: masses(:,:),patch_mass(:,:),child_mass(:,:)
    real(dp) :: film_components(size(components)),budget(5),ledger(5),global_ledger(5),alpha,speed,un,splashed,deposit
    integer :: ne,np,nc,nf,n,p,q,number,cursor,film_cursor,kind,ierr
    integer(int64) :: next_id,first_id,child_id
    logical :: valid,agreed
    ok=.false.;wall_change=0
    valid=allocated(parcels).and.allocated(component_mass).and.allocated(films).and.allocated(film_mass)
    call agree()
    if (.not.agreed) return
    np=size(parcels);nf=size(films);nc=size(components);ne=size(indices)
    valid=np<=1000000.and.nf<=1000000.and.ne<=np.and.nc>=1.and.nc<=16
    valid=valid.and.all(shape(component_mass)==[nc,np]).and.all(shape(film_mass)==[nc,nf])
    valid=valid.and.all(shape(locations)==[3,ne]).and.all(shape(normals)==[3,ne]).and.all(shape(data)==[6,ne])
    valid=valid.and.maximum_parcels>=np.and.maximum_parcels<=1000000
    valid=valid.and.maximum_films>=nf.and.maximum_films<=1000000
    valid=valid.and.all(indices>=1).and.all(indices<=np)
    call agree()
    if (.not.agreed) return
    valid=all(ieee_is_finite(normals)).and.all(ieee_is_finite(locations)).and.all(ieee_is_finite(data))
    valid=valid.and.maxval(abs(normals))<=1
    do p=1,np
      valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p)).and.parcels(p)%id<=last_id
      if (p>1) valid=valid.and.parcels(p)%id>parcels(p-1)%id
    end do
    do p=1,nf
      valid=valid.and.valid_wall_component_inventory(films(p),film_mass(:,p))
    end do
    do n=2,ne
      valid=valid.and.indices(n)>indices(n-1)
    end do
    call agree()
    if (.not.agreed) return
    number=0
    do n=1,ne
      p=indices(n);speed=norm2(parcels(p)%velocity);un=dot_product(parcels(p)%velocity,normals(:,n))
      valid=speed>0.and.un<0.and.abs(sum(normals(:,n)**2)-1)<1e-12_dp
      if (.not.valid) exit
      alpha=max(acos(-1.0_dp)/6,asin(min(1.0_dp,abs(un)/speed)))
      call classify_wall_impact(data(1,n),data(2,n),alpha,data(3,n),data(4,n),kind,splashed,deposit,valid)
      if (.not.valid) exit
      if (splashed>0) number=number+4
    end do
    valid=valid.and.number<=maximum_parcels-np
    call agree()
    if (.not.agreed) return
    next_id=last_id
    call reserve_spray_ids(comm,number,next_id,first_id,valid)
    if (.not.valid) return
    allocate(candidate(np+number),masses(nc,np+number),patches(nf+ne),patch_mass(nc,nf+ne))
    candidate(:np)=parcels;masses(:,:np)=component_mass;patches(:nf)=films;patch_mass(:,:nf)=film_mass
    cursor=np;film_cursor=nf;ledger=0
    do n=1,ne
      p=indices(n);film_components=0;child_id=0
      if (cursor<np+number) child_id=first_id+int(cursor-np,int64)
      call resolve_mixture_wall_impact(candidate(p),components,masses(:,p),locations(:,n),normals(:,n), &
        data(1,n),data(2,n),data(3,n),data(4,n),data(5,n),data(6,n),contact_angle, &
        child_id,children,child_mass,deposited,film_components,budget,kind,valid)
      if (.not.valid) exit
      q=size(children)
      valid=q<=np+number-cursor
      if (.not.valid) exit
      candidate(cursor+1:cursor+q)=children;masses(:,cursor+1:cursor+q)=child_mass;cursor=cursor+q
      if (deposited%mass>0) then
        film_cursor=film_cursor+1;patches(film_cursor)=deposited;patch_mass(:,film_cursor)=film_components
      end if
      ledger=ledger+budget
    end do
    valid=valid.and.cursor==np+number.and.film_cursor<=maximum_films
    call agree()
    if (.not.agreed) return
    call assign_parent_child_ids(comm,candidate,np,last_id,valid)
    if (.not.valid) return
    call MPI_Allreduce(ledger,global_ledger,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call move_alloc(candidate,parcels);call move_alloc(masses,component_mass)
    films=patches(:film_cursor);film_mass=patch_mass(:,:film_cursor)
    last_id=next_id;wall_change=global_ledger;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine advance_mpi_film_coupling(comm,species,transport,components,options,wall_options,state,temperature, &
      films,film_mass,lower,spacing,cells,parts,interval,faces,substeps,boundary_change,ok)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_wall_options), intent(in) :: wall_options
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    real(dp), allocatable, intent(inout) :: film_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),interval
    integer, intent(in) :: cells(3),parts(3)
    character(len=*), intent(in) :: faces(6)
    integer, intent(out) :: substeps
    real(dp), intent(out) :: boundary_change(:)
    logical, intent(out) :: ok
    type(spray_wall_film), allocatable :: patches(:),trial(:),kept(:)
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),next_gas(:,:,:,:),next_temp(:,:,:),field(:,:,:,:)
    real(dp), allocatable :: masses(:,:),next_mass(:,:),kept_mass(:,:),positions(:,:),samples(:,:),change(:,:)
    real(dp), allocatable :: rates(:,:),hv(:,:),qg(:),qw(:)
    real(dp) :: primitive(size(state,1)),tg,sound,conductivity,tdot,elapsed,dt,local_dt,wall_energy
    real(dp) :: ledger,trial_ledger,global_ledger,controls(3),low(3),high(3)
    integer :: np,nc,nvar,p,c,axis,face,attempt,count_steps,ierr
    logical :: valid,agreed,accepted
    ok=.false.;substeps=0;boundary_change=0
    valid=allocated(films).and.allocated(film_mass).and.valid_spray_options(options)
    valid=valid.and.valid_wall_options(wall_options,faces).and.ieee_is_finite(interval)
    valid=valid.and.interval>=0.and.interval<=1e6_dp.and.size(transport)==size(species)
    call agree()
    if (.not.agreed) return
    np=size(films);nc=size(components);nvar=size(state,1)
    valid=nc>=1.and.nc<=min(16,size(species)).and.np<=min(125000,wall_options%maximum_films)
    valid=valid.and.all(shape(film_mass)==[nc,np]).and.size(boundary_change)==nvar
    call agree()
    if (.not.agreed) return
    controls=[interval,real(nc,dp),merge(1.0_dp,0.0_dp,wall_options%evolve)]
    call MPI_Allreduce(controls,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.any(low/=high)) return
    do p=1,np
      valid=valid.and.valid_wall_component_inventory(films(p),film_mass(:,p))
    end do
    do c=1,nc
      valid=valid.and.valid_spray_liquid(components(c)).and.components(c)%vapor_index<=size(species)
    end do
    call agree()
    if (.not.agreed) return
    gas=state;temp=temperature;patches=films;masses=film_mass
    allocate(positions(3,np),change(nvar,np),rates(nc,np),hv(nc,np),qg(np),qw(np))
    do p=1,np
      positions(:,p)=patches(p)%position
    end do
    change=0
    call deposit_mpi_spray_gas(comm,species,gas,temp,positions,change,lower,spacing,cells,parts,faces,valid)
    if (.not.valid) return
    elapsed=0;ledger=0;count_steps=0
    do while(wall_options%evolve.and.elapsed<interval)
      valid=count_steps<options%maximum_substeps
      call agree()
      if (.not.agreed) return
      if (allocated(field)) deallocate(field)
      allocate(field(nvar+1,size(gas,2),size(gas,3),size(gas,4)))
      field(:nvar,:,:,:)=gas;field(nvar+1,:,:,:)=temp
      call sample_mpi_spray_field(comm,field,positions,lower,spacing,cells,parts,faces,samples,valid)
      if (.not.valid) return
      local_dt=interval-elapsed;rates=0;hv=0;qg=0;qw=0
      do p=1,np
        if (patches(p)%mass<=0) cycle
        axis=maxloc(abs(patches(p)%normal),dim=1)
        face=2*axis-merge(1,0,patches(p)%normal(axis)>0)
        valid=wall_options%enabled(face)
        if (.not.valid) exit
        call reactive_conserved_to_primitive(species,samples(:nvar,p),samples(nvar+1,p),primitive,tg,sound,valid)
        if (.not.valid) exit
        conductivity=wall_liquid_conductivity_at_temperature(wall_options,patches(p)%material_id,patches(p)%temperature)
        call mixture_wall_film_rates(patches(p),components,masses(:,p),species,transport,primitive,tg,options, &
          conductivity,wall_options%temperature(face),rates(:,p),hv(:,p),qg(p),qw(p),tdot,valid)
        if (.not.valid) exit
        do c=1,nc
          if (rates(c,p)>0.and.masses(c,p)>mixture_trace_fraction*patches(p)%mass) &
            local_dt=min(local_dt,options%max_fractional_change*masses(c,p)/rates(c,p))
        end do
        if (abs(tdot)>0) local_dt=min(local_dt,min(2.0_dp,options%max_fractional_change*patches(p)%temperature)/abs(tdot))
      end do
      call agree()
      if (.not.agreed) return
      call MPI_Allreduce(local_dt,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      accepted=.false.
      do attempt=1,30
        if (dt<=0.or.elapsed+dt<=elapsed) exit
        trial=patches;next_mass=masses;change=0;trial_ledger=0;valid=.true.
        do p=1,np
          if (patches(p)%mass<=0) cycle
          call advance_mixture_wall_film_exchange(trial(p),components,next_mass(:,p),rates(:,p)*dt,hv(:,p), &
            qg(p)*dt,qw(p)*dt,wall_options%minimum_height,change(:,p),wall_energy,valid,deplete_traces=.true.)
          if (valid.and.trial(p)%mass>0) valid= &
            wall_liquid_conductivity_at_temperature(wall_options,trial(p)%material_id,trial(p)%temperature)>0
          if (.not.valid) exit
          trial_ledger=trial_ledger-wall_energy
        end do
        call agree()
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
      gas=next_gas;temp=next_temp;patches=trial;masses=next_mass;ledger=ledger+trial_ledger
      elapsed=elapsed+dt;count_steps=count_steps+1
    end do
    kept=pack(patches,patches%mass>0)
    allocate(kept_mass(nc,size(kept)))
    do c=1,nc
      kept_mass(c,:)=pack(masses(c,:),patches%mass>0)
    end do
    call migrate_mpi_wall_films(comm,kept,kept_mass,lower,spacing,cells,parts,wall_options%maximum_films,valid)
    if (.not.valid) return
    call MPI_Allreduce(ledger,global_ledger,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    state=gas;temperature=temp
    call move_alloc(kept,films);call move_alloc(kept_mass,film_mass)
    boundary_change(5)=global_ledger;substeps=count_steps;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  ! Stationary patches are owned by their containing Cartesian cell. Exact upper
  ! boundary coordinates select the final cell; physical film positions stay exact.
  ! Receive order is source-rank then original local order; films have no global ID.
  subroutine migrate_mpi_wall_films(comm,films,component_mass,lower,spacing,cells,parts,maximum_films,ok,gather_rank)
    type(MPI_Comm), intent(in) :: comm
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3)
    integer, intent(in) :: cells(3),parts(3),maximum_films
    logical, intent(out) :: ok
    integer, optional, intent(in) :: gather_rank
    integer :: destination
    type(spray_wall_film), allocatable :: candidate(:)
    real(dp), allocatable :: sent(:,:),received(:,:),masses(:,:)
    integer, allocatable :: counts(:),offsets(:),incoming(:),starts(:),cursor(:),owners(:),material(:),got_material(:)
    real(dp) :: bounds(6),lowbounds(6),highbounds(6),upper(3)
    integer :: config(9),lowconfig(9),highconfig(9),nranks,ierr,nc,np,p,r,q,nreceived,width
    logical :: valid,agreed
    ok=.false.;nc=0;destination=-1
    if (present(gather_rank)) destination=gather_rank
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=allocated(films).and.allocated(component_mass)
    if (allocated(component_mass)) nc=size(component_mass,1)
    call agree()
    if (.not.agreed) return
    np=size(films)
    valid=nc>=1.and.nc<=16.and.size(component_mass,2)==np.and.np<=1000000
    valid=valid.and.destination>=-1.and.destination<nranks
    valid=valid.and.maximum_films>=0.and.maximum_films<=1000000
    valid=valid.and.spray_cartesian_owner(lower,lower,spacing,cells,parts,[.false.,.false.,.false.])>=0
    call agree()
    if (.not.agreed) return
    config=[cells,parts,nc,maximum_films,destination];bounds=[lower,spacing]
    call MPI_Allreduce(config,lowconfig,9,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,highconfig,9,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,lowbounds,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,highbounds,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(lowconfig/=highconfig).or.any(lowbounds/=highbounds).or.product(parts)/=nranks) return
    upper=lower+spacing*real(cells,dp)
    allocate(counts(nranks),offsets(nranks),incoming(nranks),starts(nranks),cursor(nranks),owners(np));counts=0
    do p=1,np
      valid=valid_wall_component_inventory(films(p),component_mass(:,p))
      if (.not.valid) exit
      r=owner(films(p)%position);valid=r>=0
      if (.not.valid) exit
      owners(p)=r+1;counts(r+1)=counts(r+1)+1
    end do
    call agree()
    if (.not.agreed) return
    call MPI_Alltoall(counts,1,MPI_INTEGER,incoming,1,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=sum(int(incoming,int64))<=maximum_films
    call agree()
    if (.not.agreed) return
    offsets(1)=0;starts(1)=0
    do r=2,nranks
      offsets(r)=offsets(r-1)+counts(r-1);starts(r)=starts(r-1)+incoming(r-1)
    end do
    cursor=offsets;nreceived=sum(incoming);width=10+nc
    allocate(sent(width,np),received(width,nreceived),material(np),got_material(nreceived))
    do p=1,np
      r=owners(p);cursor(r)=cursor(r)+1;q=cursor(r)
      sent(:,q)=[films(p)%position,films(p)%normal,films(p)%mass,films(p)%temperature, &
        films(p)%diameter,films(p)%height,component_mass(:,p)]
      material(q)=films(p)%material_id
    end do
    call MPI_Alltoallv(sent,width*counts,width*offsets,MPI_DOUBLE_PRECISION,received,width*incoming, &
      width*starts,MPI_DOUBLE_PRECISION,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Alltoallv(material,counts,offsets,MPI_INTEGER,got_material,incoming,starts,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    allocate(candidate(nreceived),masses(nc,nreceived))
    call MPI_Comm_rank(comm,r,ierr)
    if (ierr/=MPI_SUCCESS) return
    do p=1,nreceived
      candidate(p)=spray_wall_film(received(1:3,p),received(4:6,p),received(7,p),received(8,p), &
        received(9,p),received(10,p),got_material(p))
      masses(:,p)=received(11:,p)
      valid=valid_wall_component_inventory(candidate(p),masses(:,p)).and.owner(candidate(p)%position)==r
      if (.not.valid) exit
    end do
    call agree()
    if (.not.agreed) return
    call move_alloc(candidate,films);call move_alloc(masses,component_mass);ok=.true.
  contains
    integer function owner(position) result(owner_rank)
      real(dp), intent(in) :: position(3)
      real(dp) :: probe(3)
      integer :: axis
      owner_rank=-1
      if (any(position<lower).or.any(position>upper)) return
      probe=position
      do axis=1,3
        if (position(axis)==upper(axis)) probe(axis)=lower(axis)+(real(cells(axis),dp)-.5_dp)*spacing(axis)
      end do
      owner_rank=spray_cartesian_owner(probe,lower,spacing,cells,parts,[.false.,.false.,.false.])
      if (owner_rank>=0.and.destination>=0) owner_rank=destination
    end function
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine advance_mpi_spray_transport(comm,species,transport,state,temperature,spacing,cells,dt,enabled,minimum_theta,ok, &
      faces,boundary_change,inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure, &
      inflow_mach,inflow_pressure)
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp) :: pressure(6)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: enabled(4)
    real(dp), intent(out) :: minimum_theta
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),dt
    integer, intent(in) :: cells(3)
    character(len=*), optional, intent(in) :: faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(out) :: boundary_change(:)
    character(len=32) :: local_faces(6)
    real(dp) :: local_change(size(state,1)),global_change(size(state,1)),rate(size(state,1)),velocity(3,6)
    logical, intent(out) :: ok
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),candidate(:,:,:,:),next_temp(:,:,:)
    real(dp), allocatable :: halo(:,:,:,:),halo_t(:,:,:),send(:,:,:),left(:,:,:),right(:,:,:)
    real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    real(dp), allocatable :: theta(:,:,:),theta_send(:,:),theta_left(:,:),theta_right(:,:)
    real(dp) :: controls(10),minimum(10),maximum(10),local_theta,global_theta,accepted_theta
    integer :: config(9),minconfig(9),maxconfig(9),rank,nranks,ierr,first(3),extent(3)
    integer :: nx,ny,nz,nvar,stage,i,j,k,jm,km,previous,next
    logical :: valid,agreed
    local_change=0
    pressure=0
    if (present(outlet_pressure)) pressure=outlet_pressure
    if (present(boundary_change)) boundary_change=0
    ok=.false.;minimum_theta=1;accepted_theta=1;nvar=size(state,1)
    call validate_mpi_mach_controls(comm,inflow_mach,inflow_pressure,faces,agreed)
    if(.not.agreed) return
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(ieee_is_finite([spacing,dt])).and.size(transport)==size(species)
    valid=valid.and.all(spacing>=1e-12_dp).and.all(spacing<=1e9_dp).and.dt>=0.and.dt<=1e6_dp
    if (present(boundary_change)) valid=valid.and.present(faces).and.size(boundary_change)==nvar
    if (present(faces)) then
      if (any(faces=='isothermal_wall')) then
        valid=valid.and.enabled(2).and.present(wall_temperature)
        if (present(wall_temperature)) valid=valid.and.size(wall_temperature)==6
      end if
      if (present(wall_velocity)) valid=valid.and.(enabled(1).or.all(wall_velocity==0))
    end if
    valid=valid.and.all(ieee_is_finite(pressure))
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    if (any(pressure>0)) valid=valid.and.present(faces)
    if (present(faces)) valid=valid.and..not.any(pressure>0.and.faces/='outflow')
    call agree()
    if (.not.agreed) return
    config=[cells,nvar,merge(1,0,enabled),merge(1,0,present(faces))];controls=[spacing,dt,pressure]
    call MPI_Allreduce(config,minconfig,9,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,maxconfig,9,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,minimum,10,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,10,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
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
    if (present(faces)) then
      call fill_mpi_spray_ghosts(comm,gas,temp,cells,faces,halo,halo_t,valid,viscous_wall=.true., &
        inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_velocity=wall_velocity, &
          species=species,outlet_pressure=pressure, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      call agree()
      if (.not.agreed) return
      local_faces=faces;velocity=0
      if (present(wall_velocity)) velocity=wall_velocity
      if (rank/=0) then
        local_faces(1)='outflow';velocity(:,1)=0
      end if
      if (rank/=nranks-1) then
        local_faces(2)='outflow';velocity(:,2)=0
      end if
    end if
    if (dt<=tiny(1.0_dp).or..not.any(enabled(:3))) then
      ok=.true.;return
    end if
    if (.not.present(faces)) allocate(halo(nvar,0:nx+1,0:ny+1,0:nz+1),halo_t(0:nx+1,0:ny+1,0:nz+1))
    allocate(send(nvar+1,ny,nz),left(nvar+1,ny,nz),right(nvar+1,ny,nz))
    allocate(fx(nvar,0:nx,ny,nz),fy(nvar,nx,0:ny,nz),fz(nvar,nx,ny,0:nz))
    allocate(theta(nx,ny,nz),theta_send(ny,nz),theta_left(ny,nz),theta_right(ny,nz))
    previous=modulo(rank-1,nranks);next=modulo(rank+1,nranks)
    do stage=1,2
      if (present(faces)) then
        call fill_mpi_spray_ghosts(comm,gas,temp,cells,faces,halo,halo_t,valid,viscous_wall=.true., &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_velocity=wall_velocity, &
          species=species,outlet_pressure=pressure, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
        call agree()
        if (.not.agreed) return
        call physical_transport_ghosted_fluxes_3d(species,transport,gas,temp,halo,halo_t,extent,spacing,dt, &
          local_faces,enabled(1),enabled(2),enabled(3),enabled(4),fx,fy,fz,local_theta,valid,wall_temperature,velocity)
        call agree()
        if (.not.agreed) return
        rate=0
        do k=1,nz
          do j=1,ny
            if (rank==0) rate=rate+fx(:,0,j,k)*spacing(2)*spacing(3)
            if (rank==nranks-1) rate=rate-fx(:,nx,j,k)*spacing(2)*spacing(3)
          end do
        end do
        do k=1,nz
          do i=1,nx
            rate=rate+(fy(:,i,0,k)-fy(:,i,ny,k))*spacing(1)*spacing(3)
          end do
        end do
        do j=1,ny
          do i=1,nx
            rate=rate+(fz(:,i,j,0)-fz(:,i,j,nz))*spacing(1)*spacing(2)
          end do
        end do
        local_change=local_change+.5_dp*dt*rate
        valid=all(ieee_is_finite(local_change))
        call agree()
        if (.not.agreed) return
      else
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
      end if
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
    if (present(faces)) then
      call MPI_Allreduce(local_change,global_change,nvar,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      if (.not.all(ieee_is_finite(global_change))) return
      if (present(boundary_change)) boundary_change=global_change
    end if
    state=candidate;temperature=next_temp;minimum_theta=accepted_theta;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine advance_mpi_spray_les(comm,species,state,temperature,spacing,cells,dt,cs,prt,sct,ok, &
      faces,boundary_change,inflow_state,inflow_temperature,wall_temperature,wall_velocity,model,ci,outlet_pressure, &
      inflow_mach,inflow_pressure)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),dt
    integer, intent(in) :: cells(3)
    real(dp), intent(in) :: cs,prt,sct
    character(len=*), optional, intent(in) :: model
    character(len=16) :: selected_model
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp) :: pressure(6)
    real(dp), optional, intent(in) :: ci
    real(dp) :: isotropic
    integer :: model_id
    character(len=*), optional, intent(in) :: faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(out) :: boundary_change(:)
    character(len=32) :: local_faces(6)
    logical :: include_face(6)
    real(dp) :: rate(size(state,1)),change(size(state,1)),global_change(size(state,1))
    logical, intent(out) :: ok
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),candidate(:,:,:,:),next_temp(:,:,:)
    real(dp), allocatable :: halo(:,:,:,:),halo_t(:,:,:),send(:,:,:,:),left(:,:,:,:),right(:,:,:,:)
    real(dp), allocatable :: rhs(:,:,:,:)
    real(dp) :: controls(14),minimum(14),maximum(14),diffusivity
    integer :: config(6),minconfig(6),maxconfig(6),rank,nranks,ierr,first(3),extent(3)
    integer :: nx,ny,nz,nvar,stage,i,j,k,jm,km,previous,next
    logical :: valid,agreed
    pressure=0
    if (present(outlet_pressure)) pressure=outlet_pressure
    isotropic=0
    if (present(ci)) isotropic=ci
    change=0
    if (present(boundary_change)) boundary_change=0
    ok=.false.;nvar=size(state,1)
    call validate_mpi_mach_controls(comm,inflow_mach,inflow_pressure,faces,agreed)
    if(.not.agreed) return
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    selected_model='smagorinsky';model_id=1
    if (present(model)) then
      select case(model)
      case('smagorinsky')
        model_id=1
      case('wale')
        model_id=2
      case('vreman')
        model_id=3
      case default
        model_id=0
      end select
      if (model_id>0) selected_model=model
    end if
    valid=model_id>0.and.all(ieee_is_finite([spacing,dt,cs,prt,sct,isotropic,pressure]))
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    if (any(pressure>0)) valid=valid.and.present(faces)
    if (present(faces)) valid=valid.and..not.any(pressure>0.and.faces/='outflow')
    valid=valid.and.isotropic>=0.and.isotropic<=1
    valid=valid.and.all(spacing>=1e-12_dp).and.all(spacing<=1e9_dp).and.dt>=0.and.dt<=1e6_dp
    if (present(boundary_change)) valid=valid.and.present(faces).and.size(boundary_change)==nvar
    call agree()
    if (.not.agreed) return
    config=[cells,nvar,merge(1,0,present(faces)),model_id];controls=[spacing,dt,cs,prt,sct,isotropic,pressure]
    call MPI_Allreduce(config,minconfig,6,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,maxconfig,6,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,minimum,14,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,14,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
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
    include_face=.true.
    include_face(1)=rank==0;include_face(2)=rank==nranks-1
    if (present(faces)) then
      local_faces=faces
      if (rank/=0) local_faces(1)='outflow'
      if (rank/=nranks-1) local_faces(2)='outflow'
    end if
    do stage=1,2
      if (present(faces)) then
        call fill_mpi_spray_ghosts(comm,gas,temp,cells,faces,halo,halo_t,valid,viscous_wall=.true., &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_velocity=wall_velocity, &
          species=species,outlet_pressure=pressure, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
        call agree()
        if (.not.agreed) return
        call compute_les_ghosted_rhs_3d(species,halo,halo_t,[nx,ny,nz],spacing,cs,prt,sct,rhs,diffusivity,valid, &
          local_faces,rate,wall_temperature,include_face,model=selected_model,ci=isotropic)
        change=change+.5_dp*dt*rate
        valid=valid.and.all(ieee_is_finite(change))
      else
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
        call compute_les_ghosted_rhs_3d(species,halo,halo_t,[nx,ny,nz],spacing,cs,prt,sct,rhs,diffusivity,valid, &
          model=selected_model,ci=isotropic)
      end if
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
    if (present(faces)) then
      call MPI_Allreduce(change,global_change,nvar,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      if (.not.all(ieee_is_finite(global_change))) return
      if (present(boundary_change)) boundary_change=global_change
    end if
    state=gas;temperature=temp;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine advance_mpi_spray_hydro(comm,species,state,temperature,spacing,cells,dt,solver,limiter,ok, &
      faces,boundary_change,inflow_state,inflow_temperature,outlet_pressure,inflow_mach,inflow_pressure,reconstruction,difmag)
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: difmag
    real(dp) :: viscosity
    real(dp) :: pressure(6)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),dt
    integer, intent(in) :: cells(3)
    character(len=*), intent(in) :: solver,limiter
    character(len=*), optional, intent(in) :: reconstruction
    character(len=32) :: hydro,ppm_faces(6)
    character(len=*), optional, intent(in) :: faces(6)
    real(dp), optional, intent(out) :: boundary_change(:)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),candidate(:,:,:,:),next_temp(:,:,:)
    real(dp), allocatable :: halo(:,:,:,:),halo_t(:,:,:),send(:,:,:,:),left(:,:,:,:),right(:,:,:,:)
    real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    real(dp) :: controls(11),minimum(11),maximum(11),rhs(size(state,1))
    real(dp) :: local_change(size(state,1)),global_change(size(state,1)),rate(size(state,1))
    integer :: config(8),minconfig(8),maxconfig(8),rank,nranks,ierr,first(3),extent(3)
    integer :: nx,ny,nz,nvar,stage,i,j,k,jm,km,previous,next,solver_code,limiter_code,reconstruction_code
    logical :: valid,agreed
    hydro='characteristic_plm'
    if(present(reconstruction)) hydro=reconstruction
    reconstruction_code=0
    if(hydro=='characteristic_plm') reconstruction_code=1
    if(hydro=='characteristic_ppm') reconstruction_code=2
    viscosity=0
    if(hydro=='characteristic_ppm'.and.solver=='pelec') viscosity=.1_dp
    if(present(difmag)) viscosity=difmag
    pressure=0
    if (present(outlet_pressure)) pressure=outlet_pressure
    if (present(boundary_change)) boundary_change=0
    local_change=0
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
    call validate_mpi_mach_controls(comm,inflow_mach,inflow_pressure,faces,agreed)
    if(.not.agreed) return
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(ieee_is_finite([spacing,dt])).and.solver_code>0.and.limiter_code>0.and.reconstruction_code>0
    valid=valid.and.all(spacing>=1e-12_dp).and.all(spacing<=1e9_dp).and.dt>=0.and.dt<=1e6_dp
    if (present(boundary_change)) valid=valid.and.size(boundary_change)==nvar.and.present(faces)
    valid=valid.and.all(ieee_is_finite(pressure))
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    if (any(pressure>0)) valid=valid.and.present(faces)
    if (present(faces)) valid=valid.and..not.any(pressure>0.and.faces/='outflow')
    valid=valid.and.ieee_is_finite(viscosity)
    valid=valid.and.viscosity>=0.and.(viscosity==0.or.reconstruction_code==2)
    call agree()
    if (.not.agreed) return
    config=[cells,nvar,solver_code,limiter_code,merge(1,0,present(faces)),reconstruction_code]
    controls=[spacing,dt,pressure,viscosity]
    call MPI_Allreduce(config,minconfig,8,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,maxconfig,8,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,minimum,11,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,11,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
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
    if (hydro=='characteristic_ppm') then
      call ppm_step()
      return
    end if
    if (.not.present(faces)) then
      allocate(halo(nvar,nx+4,ny,nz),halo_t(nx+4,ny,nz))
      allocate(send(nvar+1,2,ny,nz),left(nvar+1,2,ny,nz),right(nvar+1,2,ny,nz))
      allocate(fx(nvar,nx+4,ny,nz),fy(nvar,nx+4,ny,nz),fz(nvar,nx+4,ny,nz))
      previous=modulo(rank-1,nranks);next=modulo(rank+1,nranks)
    else
      allocate(fx(nvar,nx+4,ny+4,nz+4),fy(nvar,nx+4,ny+4,nz+4),fz(nvar,nx+4,ny+4,nz+4))
    end if
    do stage=1,2
      if (present(faces)) then
        call physical_stage()
        if (.not.agreed) return
      else
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
      end if
      call recover_reactive_temperatures_3d(species,candidate,temp,nx,ny,nz,next_temp,valid)
      call agree()
      if (.not.agreed) return
      gas=candidate;temp=next_temp
    end do
    if (present(faces)) then
      call MPI_Allreduce(local_change,global_change,nvar,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      if (.not.all(ieee_is_finite(global_change))) return
      if (present(boundary_change)) boundary_change=global_change
    end if
    state=gas;temperature=temp;ok=.true.
  contains
    subroutine ppm_step()
      real(dp) :: theta
      logical :: inflow_mask(6),wall_mask(6)
      ppm_faces='periodic'
      if(present(faces)) ppm_faces=faces
      call fill_mpi_les_boundary_ghosts_3d(comm,species,gas,temp,cells,ppm_faces,4,halo,halo_t,valid, &
        inflow_state,inflow_temperature,outlet_pressure=pressure,inflow_mach=inflow_mach, &
        inflow_pressure=inflow_pressure,viscous_wall=.false.)
      call agree()
      if (.not.agreed) return
      allocate(fx(nvar,0:nx,ny,nz),fy(nvar,nx,0:ny,nz),fz(nvar,nx,ny,0:nz))
      wall_mask=ppm_faces=='no_slip_wall'.or.ppm_faces=='isothermal_wall'
      wall_mask(1)=wall_mask(1).and.rank==0
      wall_mask(2)=wall_mask(2).and.rank==nranks-1
      inflow_mask=ppm_faces=='inflow'
      inflow_mask(1)=inflow_mask(1).and.rank==0
      inflow_mask(2)=inflow_mask(2).and.rank==nranks-1
      call compute_reactive_ppm_halo_ctu_fluxes_3d(species,halo,halo_t,extent,spacing,dt,.false.,.true., &
        solver,fx,fy,fz,theta,valid,inflow_faces=inflow_mask,no_slip_faces=wall_mask,difmag=viscosity)
      call agree()
      if (.not.agreed) return
      do k=1,nz
      do j=1,ny
      do i=1,nx
        candidate(:,i,j,k)=gas(:,i,j,k)-(dt/spacing(1))*(fx(:,i,j,k)-fx(:,i-1,j,k))- &
          (dt/spacing(2))*(fy(:,i,j,k)-fy(:,i,j-1,k))-(dt/spacing(3))*(fz(:,i,j,k)-fz(:,i,j,k-1))
      end do
      end do
      end do
      call recover_reactive_temperatures_3d(species,candidate,temp,nx,ny,nz,next_temp,valid)
      call agree()
      if (.not.agreed) return
      rate=0
      do k=1,nz
      do j=1,ny
        if(rank==0) rate=rate+fx(:,0,j,k)*spacing(2)*spacing(3)
        if(rank==nranks-1) rate=rate-fx(:,nx,j,k)*spacing(2)*spacing(3)
      end do
      end do
      do k=1,nz
      do i=1,nx
        rate=rate+(fy(:,i,0,k)-fy(:,i,ny,k))*spacing(1)*spacing(3)
      end do
      end do
      do j=1,ny
      do i=1,nx
        rate=rate+(fz(:,i,j,0)-fz(:,i,j,nz))*spacing(1)*spacing(2)
      end do
      end do
      local_change=dt*rate
      call MPI_Allreduce(local_change,global_change,nvar,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS) return
      if(.not.all(ieee_is_finite(global_change))) return
      state=candidate;temperature=next_temp
      if(present(boundary_change)) boundary_change=global_change
      ok=.true.
    end subroutine ppm_step
    subroutine physical_stage()
      call fill_mpi_spray_ghosts(comm,gas,temp,cells,faces,halo,halo_t,valid, &
        inflow_state=inflow_state,inflow_temperature=inflow_temperature,species=species,outlet_pressure=pressure, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      call agree()
      if (.not.agreed) return
      call compute_reactive_plm_face_fluxes_3d(species,halo,halo_t,nx+4,ny+4,nz+4, &
        limiter,solver,fx,fy,fz,valid)
      call agree()
      if (.not.agreed) return
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rhs=-(fx(:,i+2,j+2,k+2)-fx(:,i+1,j+2,k+2))/spacing(1) &
              -(fy(:,i+2,j+2,k+2)-fy(:,i+2,j+1,k+2))/spacing(2) &
              -(fz(:,i+2,j+2,k+2)-fz(:,i+2,j+2,k+1))/spacing(3)
            candidate(:,i,j,k)=gas(:,i,j,k)+dt*rhs
            if (stage==2) candidate(:,i,j,k)=.5_dp*state(:,i,j,k)+.5_dp*candidate(:,i,j,k)
          end do
        end do
      end do
      rate=0
      do k=3,nz+2
        do j=3,ny+2
          if (rank==0) rate=rate+fx(:,2,j,k)*spacing(2)*spacing(3)
          if (rank==nranks-1) rate=rate-fx(:,nx+2,j,k)*spacing(2)*spacing(3)
        end do
      end do
      do k=3,nz+2
        do i=3,nx+2
          rate=rate+(fy(:,i,2,k)-fy(:,i,ny+2,k))*spacing(1)*spacing(3)
        end do
      end do
      do j=3,ny+2
        do i=3,nx+2
          rate=rate+(fz(:,i,j,2)-fz(:,i,j,nz+2))*spacing(1)*spacing(2)
        end do
      end do
      local_change=local_change+.5_dp*dt*rate
      valid=all(ieee_is_finite(local_change))
      call agree()
    end subroutine
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine advance_mpi_mixture_reaction(comm,species,reactions,transport,components,options,state,temperature, &
      parcels,component_mass,lower,spacing,cells,parts,interval,rtol,atol,method,substeps,ok, &
      hydro_solver,hydro_limiter,transport_flags,last_id,les_settings, &
      boundary_faces,boundary_change,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
      wall_films,film_mass,film_settings,les_model,les_test_filter,les_test_ratio,les_grid_ratio, &
      les_harmonic_mean,les_ci,outlet_pressure,inflow_mach,inflow_pressure,hydro_reconstruction,difmag)
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
    character(len=*), optional, intent(in) :: les_model,les_test_filter
    integer, optional, intent(in) :: les_test_ratio,les_grid_ratio
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: difmag
    real(dp) :: viscosity
    real(dp) :: pressure(6)
    real(dp), optional, intent(in) :: les_ci
    real(dp) :: isotropic
    logical, optional, intent(in) :: les_harmonic_mean
    logical :: harmonic
    character(len=32) :: selected_filter
    integer :: filter_code,test_ratio,grid_ratio
    real(dp), allocatable :: filter_weights(:)
    logical :: filter_valid
    integer :: les_code
    character(len=*), intent(in) :: method
    character(len=*), intent(in), optional :: hydro_solver,hydro_limiter,hydro_reconstruction
    logical, intent(in), optional :: transport_flags(4)
    real(dp), intent(in), optional :: les_settings(3)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(out) :: boundary_change(:)
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    integer :: first_film_steps,last_film_steps
    real(dp), allocatable :: stage_change(:)
    real(dp) :: ledger(size(state,1))
    real(dp) :: sgs(3)
    integer, intent(out) :: substeps
    logical, intent(out) :: ok
    integer(int64), optional, intent(inout) :: last_id
    integer(int64) :: candidate_id
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),masses(:,:)
    real(dp) :: controls(29),minimum(29),maximum(29),chemistry_interval,transport_theta
    integer :: ierr,first_steps,last_steps,method_code,hydro_code,limit_code,nranks,reconstruction_code
    logical :: valid,agreed,enabled(4)
    character(len=256) :: message
    ok=.false.;substeps=0;method_code=0;ledger=0
    reconstruction_code=1
    if(present(hydro_reconstruction)) then
      reconstruction_code=0
      if(hydro_reconstruction=='characteristic_plm') reconstruction_code=1
      if(hydro_reconstruction=='characteristic_ppm') reconstruction_code=2
    end if
    first_film_steps=0;last_film_steps=0
    if (present(boundary_change)) boundary_change=0
    if (present(boundary_faces)) allocate(stage_change(size(state,1)))
    select case(trim(method))
    case('implicit');method_code=1
    case('explicit');method_code=2
    case('cvode');method_code=3
    end select
    les_code=1
    if (present(les_model)) then
      select case(les_model)
      case('smagorinsky');les_code=1
      case('wale');les_code=2
      case('vreman');les_code=3
      case('dynamic');les_code=4
      case default;les_code=0
      end select
    end if
    selected_filter='box';test_ratio=2;grid_ratio=1;filter_code=2;filter_valid=.true.
    pressure=0
    if (present(outlet_pressure)) pressure=outlet_pressure
    isotropic=0
    if (present(les_ci)) isotropic=les_ci
    harmonic=.false.
    if (present(les_harmonic_mean)) harmonic=les_harmonic_mean
    if (les_code==4) then
      if (present(les_test_filter)) selected_filter=les_test_filter
      if (present(les_test_ratio)) test_ratio=les_test_ratio
      if (present(les_grid_ratio)) grid_ratio=les_grid_ratio
      filter_code=les_filter_code(selected_filter)
      call les_filter_weights(selected_filter,test_ratio,filter_weights,filter_valid)
      filter_valid=filter_valid.and.grid_ratio>=1.and.grid_ratio<=64
    end if
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
    viscosity=0
    if(reconstruction_code==2.and.hydro_code==3) viscosity=.1_dp
    if(present(difmag)) viscosity=difmag
    if (present(hydro_limiter)) then
      select case(trim(hydro_limiter))
      case('minmod');limit_code=1
      case('mc');limit_code=2
      end select
    end if
    call validate_mpi_mach_controls(comm,inflow_mach,inflow_pressure,boundary_faces,agreed)
    if(.not.agreed) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=allocated(parcels).and.allocated(component_mass).and.method_code>0.and.les_code>0.and.filter_valid
    valid=valid.and.(.not.harmonic.or.les_code==4)
    if(present(inflow_mach)) then
      if(any(inflow_mach>=0)) valid=valid.and.present(hydro_solver)
    end if
    valid=valid.and.all(ieee_is_finite(pressure))
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    if (any(pressure>0)) valid=valid.and.present(boundary_faces).and.present(hydro_solver)
    if (present(boundary_faces)) valid=valid.and..not.any(pressure>0.and.boundary_faces/='outflow')
    valid=valid.and.ieee_is_finite(isotropic)
    valid=valid.and.isotropic>=0.and.isotropic<=1.and.(isotropic==0.or.les_code/=4)
    if (present(boundary_change)) valid=valid.and.size(boundary_change)==size(state,1)
    valid=valid.and.(present(wall_films).eqv.present(film_mass)).and.(present(wall_films).eqv.present(film_settings))
    if (present(wall_films)) valid=valid.and.present(boundary_faces)
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    valid=valid.and.all(ieee_is_finite([interval,rtol,atol,sgs]))
    valid=valid.and.sgs(1)>=0.and.sgs(1)<=2.and.minval(sgs(2:))>=1e-6_dp.and.maxval(sgs(2:))<=1e6_dp
    valid=valid.and.((sgs(1)==0.and.isotropic==0).or.present(hydro_solver))
    valid=valid.and.interval>=0.and.interval<=1e6_dp.and.min(rtol,atol)>0
    valid=valid.and.(.not.any(enabled).or.present(hydro_solver))
    valid=valid.and.(present(hydro_solver).eqv.present(hydro_limiter))
    if (present(hydro_solver)) then
      valid=valid.and.hydro_code>0.and.limit_code>0.and.all(parts==[nranks,1,1])
    end if
    valid=valid.and.reconstruction_code>0
    valid=valid.and.ieee_is_finite(viscosity)
    valid=valid.and.viscosity>=0.and.(viscosity==0.or.reconstruction_code==2)
    if(present(difmag)) valid=valid.and.present(hydro_solver)
    if(present(hydro_reconstruction)) valid=valid.and.present(hydro_solver)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (.not.agreed) return
    controls=[interval,rtol,atol,real(method_code,dp),real(hydro_code,dp),real(limit_code,dp), &
      merge(1.0_dp,0.0_dp,enabled),sgs,merge(1.0_dp,0.0_dp,present(boundary_faces)),merge(1.0_dp,0.0_dp,present(wall_films)), &
      real(les_code,dp),real(filter_code,dp),real(test_ratio,dp),real(grid_ratio,dp),&
      merge(1._dp,0._dp,harmonic),isotropic,pressure,real(reconstruction_code,dp),viscosity]
    call MPI_Allreduce(controls,minimum,29,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,29,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(minimum/=maximum)) return
    candidate_id=0
    if (present(last_id)) candidate_id=last_id
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    if (present(wall_films)) then
      valid=allocated(wall_films).and.allocated(film_mass)
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS.or..not.agreed) return
      patches=wall_films;patch_mass=film_mass
      call advance_mpi_film_coupling(comm,species,transport,components,options,film_settings,gas,temp, &
        patches,patch_mass,lower,spacing,cells,parts,.5_dp*interval,boundary_faces,first_film_steps,stage_change,valid)
      if (.not.valid) return
      ledger=ledger+stage_change
    end if
    if (present(wall_films)) then
      call advance_mpi_mixture_coupling(comm,species,transport,components,options,gas,temp,drops,masses, &
        lower,spacing,cells,parts,.5_dp*interval,first_steps,valid,candidate_id,boundary_faces,stage_change, &
        patches,patch_mass,film_settings)
    else
      call advance_mpi_mixture_coupling(comm,species,transport,components,options,gas,temp,drops,masses, &
        lower,spacing,cells,parts,.5_dp*interval,first_steps,valid,candidate_id,boundary_faces,stage_change)
    end if
    if (.not.valid) return
    if (allocated(stage_change)) ledger=ledger+stage_change
    chemistry_interval=interval
    if (present(hydro_solver)) chemistry_interval=.5_dp*interval
    call advance_spray_chemistry_3d(species,reactions,gas,temp,chemistry_interval,rtol,atol,method,valid,message)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (.not.agreed) return
    if (present(hydro_solver)) then
      if (any(enabled)) then
        call advance_mpi_spray_transport(comm,species,transport,gas,temp,spacing,cells,.5_dp*interval, &
          enabled,transport_theta,valid,boundary_faces,stage_change,inflow_state,inflow_temperature, &
          wall_temperature,wall_velocity,outlet_pressure=pressure, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
        if (.not.valid) return
        if (allocated(stage_change)) ledger=ledger+stage_change
      end if
      if (sgs(1)>0.or.isotropic>0) then
        if (les_code==4) then
          if (allocated(stage_change)) stage_change=0
          call advance_mpi_dynamic_les_ssprk2_3d(comm,species,gas,temp,spacing,cells,selected_filter,test_ratio, &
            .5_dp*interval,valid,grid_ratio,harmonic,boundary_faces,stage_change, &
            inflow_state,inflow_temperature,wall_temperature,wall_velocity,pressure,inflow_mach,inflow_pressure)
        else
        call advance_mpi_spray_les(comm,species,gas,temp,spacing,cells,.5_dp*interval,sgs(1),sgs(2),sgs(3),valid, &
          boundary_faces,stage_change,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
          model=les_model,ci=isotropic,outlet_pressure=pressure, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
        end if
        if (.not.valid) return
        if (allocated(stage_change)) ledger=ledger+stage_change
      end if
      call advance_mpi_spray_hydro(comm,species,gas,temp,spacing,cells,interval,hydro_solver,hydro_limiter,valid, &
        boundary_faces,stage_change,inflow_state,inflow_temperature,outlet_pressure=pressure, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure,reconstruction=hydro_reconstruction,difmag=viscosity)
      if (.not.valid) return
      if (allocated(stage_change)) ledger=ledger+stage_change
      if (sgs(1)>0.or.isotropic>0) then
        if (les_code==4) then
          if (allocated(stage_change)) stage_change=0
          call advance_mpi_dynamic_les_ssprk2_3d(comm,species,gas,temp,spacing,cells,selected_filter,test_ratio, &
            .5_dp*interval,valid,grid_ratio,harmonic,boundary_faces,stage_change, &
            inflow_state,inflow_temperature,wall_temperature,wall_velocity,pressure,inflow_mach,inflow_pressure)
        else
        call advance_mpi_spray_les(comm,species,gas,temp,spacing,cells,.5_dp*interval,sgs(1),sgs(2),sgs(3),valid, &
          boundary_faces,stage_change,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
          model=les_model,ci=isotropic,outlet_pressure=pressure, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
        end if
        if (.not.valid) return
        if (allocated(stage_change)) ledger=ledger+stage_change
      end if
      if (any(enabled)) then
        call advance_mpi_spray_transport(comm,species,transport,gas,temp,spacing,cells,.5_dp*interval, &
          enabled,transport_theta,valid,boundary_faces,stage_change,inflow_state,inflow_temperature, &
          wall_temperature,wall_velocity,outlet_pressure=pressure, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
        if (.not.valid) return
        if (allocated(stage_change)) ledger=ledger+stage_change
      end if
      call advance_spray_chemistry_3d(species,reactions,gas,temp,chemistry_interval,rtol,atol,method,valid,message)
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      if (.not.agreed) return
    end if
    if (present(wall_films)) then
      call advance_mpi_mixture_coupling(comm,species,transport,components,options,gas,temp,drops,masses, &
        lower,spacing,cells,parts,.5_dp*interval,last_steps,valid,candidate_id,boundary_faces,stage_change, &
        patches,patch_mass,film_settings)
    else
      call advance_mpi_mixture_coupling(comm,species,transport,components,options,gas,temp,drops,masses, &
        lower,spacing,cells,parts,.5_dp*interval,last_steps,valid,candidate_id,boundary_faces,stage_change)
    end if
    if (.not.valid) return
    if (allocated(stage_change)) ledger=ledger+stage_change
    if (present(wall_films)) then
      call advance_mpi_film_coupling(comm,species,transport,components,options,film_settings,gas,temp, &
        patches,patch_mass,lower,spacing,cells,parts,.5_dp*interval,boundary_faces,last_film_steps,stage_change,valid)
      if (.not.valid) return
      ledger=ledger+stage_change
    end if
    if (first_steps>huge(substeps)-last_steps) return
    substeps=first_steps+last_steps
    if (first_film_steps>huge(substeps)-substeps) then
      substeps=0;return
    end if
    substeps=substeps+first_film_steps
    if (last_film_steps>huge(substeps)-substeps) then
      substeps=0;return
    end if
    substeps=substeps+last_film_steps
    if (.not.all(ieee_is_finite(ledger))) then
      substeps=0;return
    end if
    if (present(boundary_change)) boundary_change=ledger
    state=gas;temperature=temp
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    if (present(last_id)) last_id=candidate_id
    if (present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    ok=.true.
  end subroutine

  ! Advance KH shedding and construct local children with globally reserved IDs.
  ! gas_energy(p) is the gas energy increment at original parent p (not a rate).
  ! The caller must deposit it transactionally with these candidates.
  subroutine advance_mpi_khrt_shedding(comm,components,options,parcels,component_mass,gas_density, &
      gas_velocity,reynolds,dt,last_id,gas_energy,ok,eligible)
    type(MPI_Comm), intent(in) :: comm
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:),gas_energy(:)
    real(dp), intent(in) :: gas_density(:),gas_velocity(:,:),reynolds(:),dt
    integer(int64), intent(inout) :: last_id
    logical, intent(out) :: ok
    logical, optional, intent(in) :: eligible(:)
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
    if(present(eligible)) valid=valid.and.size(eligible)==np
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
      if(present(eligible)) then
        if(.not.eligible(p)) cycle
      end if
      fractions(:,p)=masses(:,p)/candidate(p)%mass
      call liquid_mixture_properties(components,fractions(:,p),candidate(p)%temperature,rho,cp,e,valid)
      if (.not.valid) exit
      effective=components(1);effective%density_model='constant';effective%density=rho
      effective%density_coefficients=0;effective%density_temperature_range=0
      effective%heat_capacity_model='constant';effective%heat_capacity=cp
      effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
      effective%reference_temperature=candidate(p)%temperature;effective%reference_energy=e
      effective%boiling_model='reject';effective%critical_temperature=0
      effective%normal_boiling_temperature=0
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
    call assign_parent_child_ids(comm,result,np,last_id,valid)
    if (.not.valid) return
    call move_alloc(result,parcels);call move_alloc(result_mass,component_mass);call move_alloc(energy,gas_energy)
    last_id=next_id;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine advance_mpi_mixture_coupling(comm,species,transport,components,options,state,temperature, &
      parcels,component_mass,lower,spacing,cells,parts,interval,substeps,ok,last_id,boundary_faces,boundary_change, &
      wall_films,film_mass,film_settings)
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
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(out) :: boundary_change(:)
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    type(spray_wall_film), allocatable :: patches(:),trial_patches(:)
    real(dp), allocatable :: patch_mass(:,:),trial_patch_mass(:,:),impact_previous(:,:),impact_pressure(:)
    logical, allocatable :: hit(:)
    real(dp) :: impact_wall(5),wall_ledger(5),mw(size(components))
    real(dp) :: ledger(size(state,1)),trial_ledger(size(state,1)),global_ledger(size(state,1)),exchange(5)
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
    character(len=32) :: faces(6)
    faces='periodic';ledger=0;wall_ledger=0
    if (present(boundary_faces)) faces=boundary_faces
    if (present(boundary_change)) boundary_change=0
    ok=.false.;substeps=0
    valid=allocated(parcels).and.allocated(component_mass).and.valid_spray_options(options)
    valid=valid.and.ieee_is_finite(interval)
    valid=valid.and.interval>=0.and.interval<=1e6_dp
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    valid=valid.and.size(components)>=1.and.size(components)<=16
    valid=valid.and.valid_spray_boundary_faces(faces)
    if (present(boundary_change)) valid=valid.and.size(boundary_change)==size(state,1)
    valid=valid.and.(present(wall_films).eqv.present(film_mass)).and.(present(wall_films).eqv.present(film_settings))
    if (present(wall_films)) valid=valid.and.present(boundary_faces).and.present(last_id)
    call agree()
    if (.not.agreed) return
    child_flag=merge(1,0,options%khrt_children)+2*merge(1,0,present(boundary_faces))+4*merge(1,0,present(wall_films))
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
    if (present(wall_films)) then
      valid=allocated(wall_films).and.allocated(film_mass)
      call agree()
      if (.not.agreed) return
      patches=wall_films;patch_mass=film_mass
      do n=1,nc
        mw(n)=species(components(n)%vapor_index)%molecular_weight
      end do
      call migrate_mpi_wall_films(comm,patches,patch_mass,lower,spacing,cells,parts,film_settings%maximum_films,valid)
      if (.not.valid) return
    end if
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
        trial=drops;next_mass=masses;change=0;valid=.true.;trial_id=accepted_id;reynolds=0;trial_ledger=0
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
            effective%normal_boiling_temperature=0
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
          if (.not.present(boundary_faces)) trial(p)%position=lower+modulo(trial(p)%position-lower,length)
        end do
        call agree()
        if (agreed.and.options%khrt_children) then
          call advance_mpi_khrt_shedding(comm,components,options,trial,next_mass,primitive(1,:), &
            primitive(2:4,:),reynolds,dt,trial_id,child_energy,valid)
          if (valid) valid=size(trial)<=125000
          call agree()
          if (agreed) change(5,:)=change(5,:)+child_energy
        end if
        impact_wall=0
        if (agreed.and.present(wall_films)) then
          trial_patches=patches;trial_patch_mass=patch_mass
          call prepare_impacts()
          call agree()
          if (agreed) then
            call impact_mpi_planar_walls(comm,components,options,film_settings,faces,mw,impact_pressure,impact_previous, &
              lower,lower+length,trial,next_mass,trial_patches,trial_patch_mass,trial_id,hit,impact_wall,valid)
            if (valid) valid=size(trial)<=125000
            call agree()
          end if
        end if
        if (agreed.and.present(boundary_faces)) then
          call apply_boundaries()
          call agree()
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
      ledger=ledger+trial_ledger
      if (present(wall_films)) then
        patches=trial_patches;patch_mass=trial_patch_mass;wall_ledger=wall_ledger+impact_wall
        kept=count(drops%mass>0);survivors=pack(drops,drops%mass>0)
        allocate(kept_mass(nc,kept));kept=0
        do p=1,size(drops)
          if (drops(p)%mass<=0) cycle
          kept=kept+1;kept_mass(:,kept)=masses(:,p)
        end do
        call move_alloc(survivors,drops);call move_alloc(kept_mass,masses)
      end if
      elapsed=elapsed+dt;count_steps=count_steps+1
    end do
    np=size(drops)
    kept=count(drops%mass>0);allocate(survivors(kept),kept_mass(nc,kept));kept=0
    do p=1,np
      if (drops(p)%mass<=0) cycle
      kept=kept+1;survivors(kept)=drops(p);kept_mass(:,kept)=masses(:,p)
    end do
    call migrate_spray_parcels(comm,survivors,lower,spacing,cells,parts, &
      [faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic'], &
      options%maximum_parcels,valid,kept_mass)
    if (.not.valid) return
    if (present(wall_films)) then
      call migrate_mpi_wall_films(comm,patches,patch_mass,lower,spacing,cells,parts,film_settings%maximum_films,valid)
      if (.not.valid) return
    end if
    call MPI_Allreduce(ledger,global_ledger,nvar,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    global_ledger(:5)=global_ledger(:5)-wall_ledger
    if (.not.all(ieee_is_finite(global_ledger))) return
    if (present(boundary_change)) boundary_change=global_ledger
    state=gas;temperature=temp
    call move_alloc(survivors,parcels);call move_alloc(kept_mass,component_mass)
    if (present(last_id)) last_id=accepted_id
    if (present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    substeps=count_steps;ok=.true.
  contains
    subroutine prepare_impacts()
      integer :: q,parent,lo,hi,mid
      if (allocated(impact_previous)) deallocate(impact_previous,impact_pressure)
      allocate(impact_previous(3,size(trial)),impact_pressure(size(trial)))
      valid=.true.
      do q=1,size(trial)
        parent=q
        if (q>np) then
          lo=1;hi=np;parent=0
          do while(lo<=hi)
            mid=lo+(hi-lo)/2
            if (drops(mid)%id==trial(q)%parent_id) then
              parent=mid;exit
            else if (drops(mid)%id<trial(q)%parent_id) then
              lo=mid+1
            else
              hi=mid-1
            end if
          end do
          if (parent==0) then
            valid=.false.;return
          end if
        end if
        impact_previous(:,q)=positions(:,parent);impact_pressure(q)=1
        if (drops(parent)%mass>0) impact_pressure(q)=primitive(5,parent)
      end do
    end subroutine
    subroutine apply_boundaries()
      integer :: q,c
      valid=.true.
      do q=1,size(trial)
        if (trial(q)%mass<=0) cycle
        y=next_mass(:,q)/trial(q)%mass
        call liquid_mixture_properties(components,y,trial(q)%temperature,rho,cp,e,valid)
        if (.not.valid) return
        effective=components(1);effective%density=rho;effective%density_model='constant'
        effective%density_coefficients=0;effective%density_temperature_range=0
        effective%heat_capacity=cp;effective%heat_capacity_model='constant'
        effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
        effective%reference_temperature=trial(q)%temperature;effective%reference_energy=e
        effective%boiling_model='reject';effective%critical_temperature=0
        effective%normal_boiling_temperature=0
        call apply_parcel_boundaries(trial(q),effective,lower,lower+length,faces,exchange,valid)
        if (.not.valid) return
        trial_ledger(:5)=trial_ledger(:5)+exchange
        do c=1,nc
          trial_ledger(5+components(c)%vapor_index)=trial_ledger(5+components(c)%vapor_index)+exchange(1)*y(c)
        end do
        if (trial(q)%mass<=0) next_mass(:,q)=0
      end do
      valid=all(ieee_is_finite(trial_ledger))
    end subroutine
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
          material%normal_boiling_temperature=0
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
    call assign_parent_child_ids(comm,candidate,old_count,last_id,valid)
    if (.not.valid) return
    call move_alloc(candidate,parcels)
    if (present(component_mass)) call move_alloc(masses,component_mass)
    last_id=next_id;ok=.true.
  end subroutine

  ! Order children by global parent ID, then deterministic local sibling order.
  ! Broadcast one owner's compact parent/count metadata at a time, never parcel payloads.
  subroutine assign_parent_child_ids(comm,candidate,old_count,base,ok)
    type(MPI_Comm), intent(in) :: comm
    type(spray_parcel), intent(inout) :: candidate(:)
    integer, intent(in) :: old_count
    integer(int64), intent(in) :: base
    logical, intent(out) :: ok
    integer(int64), allocatable :: groups(:,:),remote(:,:),offset(:)
    integer(int64) :: prefix,sibling
    integer :: rank,nranks,ierr,n,p,j,k,sender,count,ios
    logical :: valid,agreed
    ok=.false.;n=0
    do p=old_count+1,size(candidate)
      if (p==old_count+1) then
        n=n+1
      else if (candidate(p)%parent_id/=candidate(p-1)%parent_id) then
        n=n+1
      end if
    end do
    allocate(groups(2,n),offset(n),stat=ios)
    valid=ios==0
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    groups=0;offset=0;j=0
    do p=old_count+1,size(candidate)
      if (p==old_count+1) then
        j=j+1
      else if (candidate(p)%parent_id/=candidate(p-1)%parent_id) then
        j=j+1
      end if
      groups(1,j)=candidate(p)%parent_id;groups(2,j)=groups(2,j)+1
    end do
    do j=2,n
      offset(j)=offset(j-1)+groups(2,j-1)
    end do
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    do sender=0,nranks-1
      count=n
      call MPI_Bcast(count,1,MPI_INTEGER,sender,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      allocate(remote(2,count),stat=ios)
      valid=ios==0
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS.or..not.agreed) return
      if (rank==sender) remote=groups
      call MPI_Bcast(remote,2*count,MPI_INTEGER8,sender,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      valid=.true.
      if (rank/=sender) then
        k=1;prefix=0
        do j=1,n
          do while(k<=count)
            if (remote(1,k)>=groups(1,j)) exit
            prefix=prefix+remote(2,k);k=k+1
          end do
          if (k<=count) then
            if (remote(1,k)==groups(1,j)) valid=.false.
          end if
          offset(j)=offset(j)+prefix
        end do
      end if
      deallocate(remote)
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS.or..not.agreed) return
    end do
    j=0;sibling=0
    do p=old_count+1,size(candidate)
      if (p==old_count+1) then
        j=j+1;sibling=0
      else if (candidate(p)%parent_id/=candidate(p-1)%parent_id) then
        j=j+1;sibling=0
      end if
      sibling=sibling+1;candidate(p)%id=base+offset(j)+sibling
    end do
    ok=.true.
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
    type(spray_parcel), allocatable :: outgoing(:)
    real(dp), allocatable :: masses(:,:)
    integer, allocatable :: counts(:),offsets(:)
    integer :: rank,nranks,ierr,nc,p,r,config(11),lo(11),hi(11),destination
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
    if (present(component_mass)) then
      call exchange_spray_parcel_bins(comm,outgoing,counts,offsets,maximum_parcels,valid,masses)
    else
      call exchange_spray_parcel_bins(comm,outgoing,counts,offsets,maximum_parcels,valid)
    end if
    if (.not.valid) return
    do p=1,size(outgoing)
      if (destination>=0) then
        valid=valid.and.rank==destination
      else
        valid=valid.and.spray_cartesian_owner(outgoing(p)%position,lower,spacing,cells,parts,periodic)==rank
      end if
    end do
    call agree()
    if (.not.global_valid) return
    call move_alloc(outgoing,parcels)
    if (present(component_mass)) call move_alloc(masses,component_mass)
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,global_valid,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) global_valid=.false.
    end subroutine
  end subroutine
end module
