program test_spray_mpi
  use mpi_f08
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use les_smagorinsky_3d_mod, only: advance_les_ssprk2_3d
  use reactive_transport_3d_mod, only: advance_reactive_transport_3d
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d,advance_reactive_euler_ssprk2_plm_3d, &
    advance_reactive_boundary_hydro_3d
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use spray_parcel_mod
  use spray_mpi_checkpoint_mod, only: mpi_spray_checkpoint
  use spray_mpi_mod
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,mixture_impact_at_planar_wall
  use spray_partition_mod, only: spray_cartesian_block
  use spray_coupling_3d_mod, only: particle_boundary_stencil,mixture_system_totals
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:),bad_reactions(:)
  type(spray_liquid) :: liquids(2)
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_options) :: options
  real(dp) :: local_totals(15),before(15),after(15)
  integer :: steps,checkpoint_event,checkpoint_unit
  real(dp) :: checkpoint_time
  character(len=80) :: checkpoint_path
  real(dp) :: fractions(2),ledger(15),direction(3)
  type(spray_parcel), allocatable :: parcels(:),saved(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:)
  integer :: rank,nranks,ierr,p,target,first(3),dims(3),cells(3),i,j,k,n,indices(3,8)
  real(dp), allocatable :: field(:,:,:,:),positions(:,:),samples(:,:),saved_samples(:,:)
  real(dp), allocatable :: contributions(:,:),reference(:,:,:,:),saved_field(:,:,:,:)
  real(dp), allocatable :: temperature(:,:,:),saved_temperature(:,:,:),reference_temperature(:,:,:)
  real(dp) :: primitive(15),gas_cell(15),gas_t,sound
  integer :: source
  real(dp) :: weights(8),expected_value,serial_theta,mpi_theta
  character(len=16) :: faces(6)
  integer(int64) :: last_id,first_id,saved_last
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call check_film_migration()
  call check_film_query()
  allocate(parcels(2),masses(2,2))
  do p=1,2
    parcels(p)=spray_parcel(int(2*rank+p,int64),[.25_dp,.25_dp,.25_dp], &
      [1.0_dp,2.0_dp,3.0_dp],1e-9_dp,350.0_dp,100.0_dp)
    parcels(p)%khrt_initial_count=100;parcels(p)%khrt_shed_volume=real(parcels(p)%id,dp)*1e-15_dp
    masses(:,p)=parcels(p)%mass*[.25_dp,.75_dp]
  end do
  saved=parcels;saved_mass=masses
  call migrate_spray_parcels(MPI_COMM_WORLD,parcels,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
    [nranks,1,1],[nranks,1,1],[.true.,.false.,.false.],2,ok,masses)
  call require(.not.ok.and.all(parcels%id==saved%id).and.all(masses==saved_mass),'capacity collective rollback')
  call migrate_spray_parcels(MPI_COMM_WORLD,parcels,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
    [nranks,1,1],[nranks,1,1],[.true.,.false.,.false.],100,ok,masses)
  call require(ok,'gather migration')
  call require(size(parcels)==merge(2*nranks,0,rank==0),'empty destinations')
  do p=1,size(parcels)
    call require(parcels(p)%id==p,'received ID ordering')
    call require(parcels(p)%khrt_shed_volume==real(p,dp)*1e-15_dp,'KH payload')
    target=modulo(p-1,nranks);parcels(p)%position(1)=target+.25_dp
  end do
  call migrate_spray_parcels(MPI_COMM_WORLD,parcels,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
    [nranks,1,1],[nranks,1,1],[.true.,.false.,.false.],100,ok,masses)
  call require(ok.and.size(parcels)==2,'redistribute from empty sources')
  do p=1,2
    call require(parcels(p)%id==rank+1+(p-1)*nranks,'destination ID ordering')
    call require(all(masses(:,p)==parcels(p)%mass*[.25_dp,.75_dp]),'component payload')
    call require(all(parcels(p)%velocity==[1.0_dp,2.0_dp,3.0_dp]),'velocity payload')
  end do
  saved=parcels;saved_mass=masses
  if (rank==0) masses(1,1)=-1
  call migrate_spray_parcels(MPI_COMM_WORLD,parcels,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
    [nranks,1,1],[nranks,1,1],[.true.,.false.,.false.],100,ok,masses)
  call require(.not.ok.and.all(parcels%id==saved%id),'one-rank invalid input rejects collectively')
  if (rank/=0) call require(all(masses==saved_mass),'other ranks preserve mass')
  masses=saved_mass;last_id=int(2*nranks,int64)
  call split_mpi_spray_parcels(MPI_COMM_WORLD,parcels,60.0_dp,merge(2,100,rank==0),last_id,ok,masses)
  call require(.not.ok.and.size(parcels)==2.and.last_id==2*nranks,'split capacity collective rollback')
  call require(all(parcels%id==saved%id).and.all(masses==saved_mass),'split failure preserves payload')
  call split_mpi_spray_parcels(MPI_COMM_WORLD,parcels,60.0_dp,100,last_id,ok,masses)
  call require(ok.and.size(parcels)==4.and.last_id==4*nranks,'collective split ID allocation')
  call require(all(parcels(3:4)%id==[2*nranks+2*rank+1,2*nranks+2*rank+2]),'split globally disjoint IDs')
  call require(all(parcels(3:4)%parent_id==saved%id),'split lineage')
  call require(sum(parcels%multiplicity)==sum(saved%multiplicity),'split multiplicity conserved')
  call require(sum(parcels%khrt_initial_count)==sum(saved%khrt_initial_count),'split KH count conserved')
  call require(abs(sum(parcels%khrt_shed_volume)-sum(saved%khrt_shed_volume))<1e-28_dp,'split KH volume conserved')
  do p=1,4
    call require(all(masses(:,p)==parcels(p)%mass*[.25_dp,.75_dp]),'split component inheritance')
    call require(parcels(p)%multiplicity<=60,'split count cap')
  end do
  last_id=100;first_id=-7
  call reserve_spray_ids(MPI_COMM_WORLD,rank,last_id,first_id,ok)
  call require(ok.and.last_id==100+nranks*(nranks-1)/2,'collective ID high-water mark')
  if (rank==0) then
    call require(first_id==0,'zero request has no ID range')
  else
    call require(first_id==101+rank*(rank-1)/2,'disjoint ID prefix range')
  end if
  saved_last=last_id
  call reserve_spray_ids(MPI_COMM_WORLD,1,last_id,first_id,ok)
  call require(ok.and.first_id==saved_last+rank+1.and.last_id==saved_last+nranks,'successive ID reservations')
  saved_last=last_id;first_id=-9
  call reserve_spray_ids(MPI_COMM_WORLD,merge(-1,1,rank==0),last_id,first_id,ok)
  call require(.not.ok.and.last_id==saved_last.and.first_id==-9,'invalid ID request collective rollback')
  last_id=huge(last_id)-1;first_id=-11
  call reserve_spray_ids(MPI_COMM_WORLD,1,last_id,first_id,ok)
  call require(.not.ok.and.last_id==huge(last_id)-1.and.first_id==-11,'ID overflow rollback')
  last_id=int(rank,int64);first_id=-13
  call reserve_spray_ids(MPI_COMM_WORLD,1,last_id,first_id,ok)
  call require(.not.ok.and.last_id==rank.and.first_id==-13,'inconsistent ID history rejected')
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'injection thermo')
  call initialize_spray_liquid(species,'H2O',1000.0_dp,4184.0_dp,373.15_dp,101325.0_dp,2.257e6_dp,liquids(1),ok)
  call require(ok,'injection water')
  call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp,liquids(2),ok)
  call require(ok,'injection peroxide')
  deallocate(parcels,masses);allocate(parcels(0),masses(2,0));last_id=0
  fractions=[.2_dp+.1_dp*rank,.8_dp-.1_dp*rank]
  call inject_mpi_cone_parcels(MPI_COMM_WORLD,liquids(1),1e-10_dp*rank,rank, &
    [rank+.25_dp,.25_dp,.25_dp],[0.0_dp,0.0_dp,1.0_dp],10.0_dp,0.0_dp,2e-5_dp,330.0_dp, &
    100,last_id,parcels,ledger,ok,liquids,fractions,masses)
  call require(ok.and.size(parcels)==rank.and.last_id==nranks*(nranks-1)/2,'collective mixture injection')
  call require(abs(ledger(1)-1e-10_dp*rank)<1e-24_dp,'injected total mass')
  call require(abs(ledger(5+liquids(1)%vapor_index)-ledger(1)*fractions(1))<1e-24_dp,'injected component ledger')
  do p=1,size(parcels)
    call require(all(masses(:,p)==parcels(p)%mass*fractions),'injected composition')
    call require(parcels(p)%id==rank*(rank-1)/2+p,'injected unique IDs')
  end do
  saved=parcels;saved_mass=masses;saved_last=last_id;direction=[0.0_dp,0.0_dp,1.0_dp]
  if (rank==nranks-1) direction=0
  call inject_mpi_cone_parcels(MPI_COMM_WORLD,liquids(1),1e-10_dp,1, &
    [rank+.25_dp,.25_dp,.25_dp],direction,10.0_dp,0.0_dp,2e-5_dp,330.0_dp, &
    100,last_id,parcels,ledger,ok,liquids,fractions,masses)
  call require(.not.ok.and.last_id==saved_last.and.all(ledger==0),'injection failure restores global IDs and ledger')
  call require(all(parcels%id==saved%id).and.all(masses==saved_mass),'injection failure preserves local inventories')
  cells=[2*nranks+1,3,2];faces='slip_wall';faces(1:2)='periodic'
  call spray_cartesian_block(cells,[nranks,1,1],rank,first,dims,ok)
  call require(ok,'distributed field block')
  allocate(field(2,dims(1),dims(2),dims(3)),positions(3,3))
  do k=1,dims(3)
    do j=1,dims(2)
      do i=1,dims(1)
        field(1,i,j,k)=real(i+first(1)+10*j+100*k,dp);field(2,i,j,k)=-field(1,i,j,k)
      end do
    end do
  end do
  positions(:,1)=[2.9_dp,1.3_dp,.7_dp]
  positions(:,2)=[real(cells(1),dp)-.25_dp,2.5_dp,1.0_dp]
  positions(:,3)=[-.25_dp,.1_dp,.1_dp]
  call sample_mpi_spray_field(MPI_COMM_WORLD,field,positions,[0.0_dp,0.0_dp,0.0_dp], &
    [1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,samples,ok)
  call require(ok,'distributed field sampling')
  do p=1,3
    call particle_boundary_stencil(positions(:,p),[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
      cells,faces,indices,weights,ok)
    call require(ok,'serial reference stencil');expected_value=0
    do n=1,8
      expected_value=expected_value+weights(n)*real(indices(1,n)+10*indices(2,n)+100*indices(3,n),dp)
    end do
    call require(samples(1,p)==expected_value.and.samples(2,p)==-expected_value,'distributed versus serial interpolation')
  end do
  if (rank==0) positions=positions(:,:0)
  call sample_mpi_spray_field(MPI_COMM_WORLD,field,positions,[0.0_dp,0.0_dp,0.0_dp], &
    [1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,samples,ok)
  call require(ok.and.size(samples,2)==merge(0,3,rank==0),'empty sampling source')
  saved_samples=samples
  if (rank==nranks-1) positions(2,1)=-1
  call sample_mpi_spray_field(MPI_COMM_WORLD,field,positions,[0.0_dp,0.0_dp,0.0_dp], &
    [1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,samples,ok)
  call require(.not.ok.and.all(samples==saved_samples),'one-rank invalid sample preserves every result')
  deallocate(positions);allocate(positions(3,3),contributions(2,3),reference(2,cells(1),cells(2),cells(3)))
  positions(:,1)=[2.9_dp,1.3_dp,.7_dp]
  positions(:,2)=[real(cells(1),dp)-.25_dp,2.5_dp,1.0_dp]
  positions(:,3)=[-.25_dp,.1_dp,.1_dp]
  contributions(1,:)=real(rank+1,dp);contributions(2,:)=-2*real(rank+1,dp)
  reference=0;field=0
  do source=0,nranks-1
    do p=1,3
      call particle_boundary_stencil(positions(:,p),[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
        cells,faces,indices,weights,ok)
      call require(ok,'serial deposition stencil')
      do n=1,8
        i=indices(1,n);j=indices(2,n);k=indices(3,n)
        reference(:,i,j,k)=reference(:,i,j,k)+weights(n)*[real(source+1,dp),-2*real(source+1,dp)]
      end do
    end do
  end do
  call deposit_mpi_spray_field(MPI_COMM_WORLD,field,positions,contributions,[0.0_dp,0.0_dp,0.0_dp], &
    [1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,ok)
  call require(ok,'distributed deposition')
  call require(all(field==reference(:,first(1)+1:first(1)+dims(1),:,:)),'distributed versus serial deposition')
  if (rank==0) then
    positions=positions(:,:0);contributions=contributions(:,:0)
  end if
  call deposit_mpi_spray_field(MPI_COMM_WORLD,field,positions,contributions,[0.0_dp,0.0_dp,0.0_dp], &
    [1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,ok)
  call require(ok,'empty depositing rank')
  saved_field=field
  if (rank==nranks-1) positions(2,1)=-1
  call deposit_mpi_spray_field(MPI_COMM_WORLD,field,positions,contributions,[0.0_dp,0.0_dp,0.0_dp], &
    [1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,ok)
  call require(.not.ok.and.all(field==saved_field),'invalid deposition collective rollback')
  ! A physical gas candidate commits only after all owner blocks recover.
  deallocate(field,positions,contributions)
  allocate(field(15,dims(1),dims(2),dims(3)),temperature(dims(1),dims(2),dims(3)))
  allocate(positions(3,1),contributions(15,1))
  primitive=0;primitive(1)=1;primitive(5)=101325;primitive(15)=1
  call reactive_primitive_to_conserved(species,primitive,gas_cell,gas_t,sound,ok)
  call require(ok,'gas deposition initial thermodynamics')
  do k=1,dims(3)
    do j=1,dims(2)
      do i=1,dims(1)
        field(:,i,j,k)=gas_cell
      end do
    end do
  end do
  temperature=gas_t;positions(:,1)=[real(first(1),dp)+.5_dp,.5_dp,.5_dp]
  contributions=0;contributions(5,1)=100
  reference=field
  call deposit_mpi_spray_field(MPI_COMM_WORLD,reference,positions,contributions,[0.0_dp,0.0_dp,0.0_dp], &
    [1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,ok)
  call require(ok,'gas deposition reference')
  allocate(reference_temperature(dims(1),dims(2),dims(3)))
  call recover_reactive_temperatures_3d(species,reference,temperature,dims(1),dims(2),dims(3),reference_temperature,ok)
  call require(ok,'serial thermodynamic reference')
  call deposit_mpi_spray_gas(MPI_COMM_WORLD,species,field,temperature,positions,contributions, &
    [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,ok)
  call require(ok.and.all(field==reference).and.all(temperature==reference_temperature),'gas and temperature commit')
  call require(temperature(1,1,1)>gas_t,'deposited heat raises temperature')
  saved_field=field;saved_temperature=temperature
  if (rank==nranks-1) contributions(5,1)=-1e12_dp
  call deposit_mpi_spray_gas(MPI_COMM_WORLD,species,field,temperature,positions,contributions, &
    [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,ok)
  call require(.not.ok.and.all(field==saved_field).and.all(temperature==saved_temperature), &
    'one owner thermal failure rolls back all owners')
  contributions=0
  if (rank==nranks-1) contributions(6,1)=-.1_dp
  call deposit_mpi_spray_gas(MPI_COMM_WORLD,species,field,temperature,positions,contributions, &
    [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,ok)
  call require(.not.ok.and.all(field==saved_field).and.all(temperature==saved_temperature), &
    'negative species rolls back all owners')
  call load_h2o2_full_transport(transport,ok);call require(ok,'coupling transport')
  call check_kh_release()
  call check_film_coupling()
  call check_wall_impacts()
  call check_wall_crossing()
  call check_coupled_impact()
  call check_particle_boundaries()
  ! Existing injected parcels include an empty rank and move to their spatial owners.
  call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],parcels,liquids,masses,local_totals,ok)
  call require(ok,'initial coupled totals')
  call MPI_Allreduce(local_totals,before,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  options%heat=.false.;options%evaporation=.false.
  saved=parcels;saved_mass=masses;saved_field=field;saved_temperature=temperature
  options%maximum_substeps=1
  call advance_mpi_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,field,temperature, &
    parcels,masses,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],.01_dp,steps,ok)
  call require(.not.ok.and.steps==0,'coupling substep exhaustion rejects interval')
  call require(all(field==saved_field).and.all(temperature==saved_temperature).and.all(masses==saved_mass), &
    'coupling exhaustion preserves gas and composition')
  do p=1,size(parcels)
    call require(all(parcels(p)%position==saved(p)%position).and.all(parcels(p)%velocity==saved(p)%velocity), &
      'coupling exhaustion preserves particles')
  end do
  options%maximum_substeps=100000
  call advance_mpi_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,field,temperature, &
    parcels,masses,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],1e-4_dp,steps,ok)
  call require(ok.and.steps>0,'coupled MPI mixture advancement')
  call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],parcels,liquids,masses,local_totals,ok)
  call require(ok,'final coupled totals')
  call MPI_Allreduce(local_totals,after,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(all(abs(after-before)<2e-12_dp*max(1e-12_dp,abs(before))),'MPI gas particle inventory conservation')
  do p=1,size(parcels)
    call require(parcels(p)%velocity(3)<10.0_dp,'coupled drag changes velocity')
  end do
  ! Representative breakup changes drop mass/multiplicity, retaining full component inventory.
  saved=parcels;saved_mass=masses;saved_field=field;saved_temperature=temperature
  options%drag=.false.;options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
  do n=1,2
    parcels=saved;masses=saved_mass;field=saved_field;temperature=saved_temperature
    options%breakup_model='etab'
    if (n==2) options%breakup_model='khrt'
    do p=1,size(parcels)
      parcels(p)%velocity=[1000.0_dp,0.0_dp,0.0_dp];parcels(p)%deformation=1.1_dp
    end do
    call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],parcels,liquids,masses,local_totals,ok)
    call require(ok,'breakup initial totals')
    call MPI_Allreduce(local_totals,before,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call advance_mpi_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,field,temperature, &
      parcels,masses,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],1e-8_dp,steps,ok)
    call require(ok.and.steps>0,'MPI representative breakup accepts')
    call require(size(parcels)==size(saved),'representative breakup preserves parcel count')
    do p=1,size(parcels)
      call require(parcels(p)%mass<saved(p)%mass,'breakup reduces physical drop mass')
      call require(parcels(p)%multiplicity>saved(p)%multiplicity,'breakup increases multiplicity')
      call require(all(abs(masses(:,p)*parcels(p)%multiplicity-saved_mass(:,p)*saved(p)%multiplicity) &
        <2e-12_dp*saved_mass(:,p)*saved(p)%multiplicity),'breakup retains each component inventory')
    end do
    call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],parcels,liquids,masses,local_totals,ok)
    call require(ok,'breakup final totals')
    call MPI_Allreduce(local_totals,after,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(abs(after-before)<2e-12_dp*max(1e-12_dp,abs(before))),'breakup preserves coupled totals')
  end do
  parcels=saved;masses=saved_mass;field=saved_field;temperature=saved_temperature
  options%breakup_model='none';options%drag=.true.;options%surface_tension=0;options%liquid_viscosity=0
  call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],parcels,liquids,masses,local_totals,ok)
  call require(ok,'restore pre-breakup totals')
  call MPI_Allreduce(local_totals,after,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  options%heat=.true.;options%evaporation=.true.
  before=after
  call advance_mpi_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,field,temperature, &
    parcels,masses,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],1e-7_dp,steps,ok)
  call require(ok.and.steps>0,'coupled mixture heating and evaporation')
  call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],parcels,liquids,masses,local_totals,ok)
  call require(ok,'evaporating coupled totals')
  call MPI_Allreduce(local_totals,after,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(all(abs(after-before)<2e-12_dp*max(1e-12_dp,abs(before))),'MPI evaporating inventory conservation')
  local_totals=0;local_totals(1)=sum(field(5+liquids(1)%vapor_index,:,:,:))
  call MPI_Allreduce(local_totals(1),after(1),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(after(1)>0,'liquid vapor reaches distributed gas')
  call load_h2o2_full_mechanism(reactions,ok);call require(ok,'MPI reaction mechanism')
  primitive=0;primitive(1)=1;primitive(6)=.03_dp;primitive(9)=.23_dp;primitive(15)=.74_dp
  primitive(5)=1400*mixture_specific_gas_constant(species,primitive(6:),ok)
  call require(ok,'reacting gas constant')
  call reactive_primitive_to_conserved(species,primitive,gas_cell,gas_t,sound,ok)
  call require(ok,'reacting initial state')
  do k=1,dims(3)
    do j=1,dims(2)
      do i=1,dims(1)
        field(:,i,j,k)=gas_cell
      end do
    end do
  end do
  temperature=gas_t;saved_field=field;saved_temperature=temperature;saved=parcels;saved_mass=masses
  bad_reactions=reactions
  if (rank==nranks-1) bad_reactions(1)%forward_rate%pre_exponential=-1
  call advance_mpi_mixture_reaction(MPI_COMM_WORLD,species,bad_reactions,transport,liquids,options, &
    field,temperature,parcels,masses,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1], &
    1e-7_dp,1e-7_dp,1e-12_dp,'implicit',steps,ok)
  call require(.not.ok.and.steps==0,'one-rank chemistry failure rejects split update')
  call require(all(field==saved_field).and.all(temperature==saved_temperature).and.all(masses==saved_mass), &
    'reaction failure rolls back gas and liquid composition')
  do p=1,size(parcels)
    call require(parcels(p)%id==saved(p)%id.and.all(parcels(p)%position==saved(p)%position).and. &
      all(parcels(p)%velocity==saved(p)%velocity).and.parcels(p)%temperature==saved(p)%temperature, &
      'reaction failure rolls back already advanced particles')
  end do
  call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],parcels,liquids,masses,local_totals,ok)
  call require(ok,'reacting initial budgets')
  call MPI_Allreduce(local_totals,before,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call advance_mpi_mixture_reaction(MPI_COMM_WORLD,species,reactions,transport,liquids,options, &
    field,temperature,parcels,masses,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1], &
    1e-7_dp,1e-7_dp,1e-12_dp,'implicit',steps,ok)
  call require(ok.and.steps>=2,'MPI split particle reaction update')
  call require(maxval(abs(field(6,:,:,:)-saved_field(6,:,:,:)))>1e-12_dp,'gas reaction advances hydrogen')
  call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],parcels,liquids,masses,local_totals,ok)
  call require(ok,'reacting final budgets')
  call MPI_Allreduce(local_totals,after,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(all(abs(after(:5)-before(:5))<2e-12_dp*max(1e-12_dp,abs(before(:5)))), &
    'reaction coupling conserves global hydrodynamic inventories')
  call require(all(abs(element_totals(after)-element_totals(before))< &
    2e-10_dp*max(1e-12_dp,abs(element_totals(before)))),'reaction coupling conserves element inventories')
  cells=[2*nranks+1,3,3]
  call spray_cartesian_block(cells,[nranks,1,1],rank,first,dims,ok)
  call require(ok,'hydro block')
  deallocate(reference,reference_temperature)
  allocate(reference(15,cells(1),3,3),reference_temperature(cells(1),3,3))
  do k=1,3
    do j=1,3
      do i=1,cells(1)
        primitive=0;primitive(15)=1
        primitive(1)=1+.1_dp*sin(2*acos(-1.0_dp)*real(i-1,dp)/cells(1))
        primitive(2:4)=[2.0_dp+.1_dp*j,-1.0_dp+.1_dp*k,.5_dp+.1_dp*i]
        primitive(5)=101325+100*j*k
        call reactive_primitive_to_conserved(species,primitive,reference(:,i,j,k),reference_temperature(i,j,k),sound,ok)
        call require(ok,'hydro initial state')
      end do
    end do
  end do
  field=reference(:,first(1)+1:first(1)+dims(1),:,:)
  temperature=reference_temperature(first(1)+1:first(1)+dims(1),:,:)
  saved_field=field;saved_temperature=temperature
  call advance_reactive_euler_ssprk2_plm_3d(species,reference,reference_temperature,cells(1),3,3, &
    1.0_dp,1.0_dp,1.0_dp,1e-5_dp,'mc','hllc',ok)
  call require(ok,'serial hydro reference')
  call check_boundary_ghosts()
  call check_boundary_split()
  call advance_mpi_spray_hydro(MPI_COMM_WORLD,species,field,temperature,[1.0_dp,1.0_dp,1.0_dp], &
    cells,1e-5_dp,'hllc','mc',ok)
  call require(ok,'distributed hydro')
  call require(maxval(abs(field-reference(:,first(1)+1:first(1)+dims(1),:,:)))<1e-10_dp, &
    'distributed hydro matches serial field')
  call require(maxval(abs(temperature-reference_temperature(first(1)+1:first(1)+dims(1),:,:)))<1e-10_dp, &
    'distributed hydro matches serial temperatures')
  call require(maxval(abs(field-saved_field))>1e-6_dp,'hydro actually advances')
  call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],parcels,liquids,masses,local_totals,ok)
  call require(ok,'full MPI split initial budgets')
  call MPI_Allreduce(local_totals,before,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call advance_mpi_mixture_reaction(MPI_COMM_WORLD,species,reactions,transport,liquids,options, &
    field,temperature,parcels,masses,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1], &
    1e-7_dp,1e-7_dp,1e-12_dp,'implicit',steps,ok,'hllc','mc')
  call require(ok.and.steps>=2,'MPI particle reaction hydro split')
  call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],parcels,liquids,masses,local_totals,ok)
  call require(ok,'full MPI split final budgets')
  call MPI_Allreduce(local_totals,after,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(all(abs(after(:5)-before(:5))<2e-12_dp*max(1e-12_dp,abs(before(:5)))), &
    'MPI particle reaction hydro conservation')
  do k=1,3
    do j=1,3
      do i=1,cells(1)
        primitive=0;primitive(1)=1;primitive(2:4)=[2.0_dp*j,.5_dp*k,.25_dp*i]
        primitive(5+liquids(1)%vapor_index)=.1_dp+.02_dp*sin(2*acos(-1.0_dp)*real(i-1,dp)/cells(1))
        primitive(15)=1-primitive(5+liquids(1)%vapor_index)
        primitive(5)=(500+20*i+10*j)*mixture_specific_gas_constant(species,primitive(6:),ok)
        call require(ok,'diffusing gas constant')
        call reactive_primitive_to_conserved(species,primitive,reference(:,i,j,k),reference_temperature(i,j,k),sound,ok)
        call require(ok,'transport initial state')
      end do
    end do
  end do
  field=reference(:,first(1)+1:first(1)+dims(1),:,:)
  temperature=reference_temperature(first(1)+1:first(1)+dims(1),:,:)
  saved_field=field
  saved_temperature=temperature
  call advance_mpi_spray_les(MPI_COMM_WORLD,species,field,temperature,[.001_dp,.001_dp,.001_dp], &
    cells,1e6_dp,.1_dp,.7_dp,.7_dp,ok)
  call require(.not.ok.and.all(field==saved_field).and.all(temperature==saved_temperature), &
    'MPI SGS excessive diffusion step rolls back')
  before=sum(sum(sum(reference,dim=4),dim=3),dim=2)
  call advance_les_ssprk2_3d(species,reference,reference_temperature,[.001_dp,.001_dp,.001_dp], &
    .1_dp,.7_dp,.7_dp,1e-6_dp,ok)
  call require(ok,'serial SGS reference')
  after=sum(sum(sum(reference,dim=4),dim=3),dim=2)
  call require(all(abs(after-before)<2e-12_dp*max(1e-12_dp,abs(before))),'SGS global conservation')
  call advance_mpi_spray_les(MPI_COMM_WORLD,species,field,temperature,[.001_dp,.001_dp,.001_dp], &
    cells,1e-6_dp,.1_dp,.7_dp,.7_dp,ok)
  call require(ok,'MPI SGS update')
  call require(all(abs(field-reference(:,first(1)+1:first(1)+dims(1),:,:))< &
    2e-12_dp*max(1.0_dp,abs(reference(:,first(1)+1:first(1)+dims(1),:,:)))),'MPI SGS matches serial fields')
  call require(maxval(abs(temperature-reference_temperature(first(1)+1:first(1)+dims(1),:,:)))<1e-8_dp, &
    'MPI SGS matches serial temperature')
  call require(maxval(abs(field-saved_field))>1e-6_dp,'SGS changes the gas')
  saved_field=field
  call advance_reactive_transport_3d(species,transport,reference,reference_temperature,cells(1),3,3, &
    .001_dp,.001_dp,.001_dp,1e-6_dp,.true.,.true.,.true.,.true.,serial_theta,ok)
  call require(ok,'serial transport reference')
  call advance_mpi_spray_transport(MPI_COMM_WORLD,species,transport,field,temperature,[.001_dp,.001_dp,.001_dp], &
    cells,1e-6_dp,[.true.,.true.,.true.,.true.],mpi_theta,ok)
  call require(ok.and.abs(mpi_theta-serial_theta)<1e-14_dp,'MPI transport limiter matches serial')
  call require(all(abs(field-reference(:,first(1)+1:first(1)+dims(1),:,:))< &
    2e-12_dp*max(1.0_dp,abs(reference(:,first(1)+1:first(1)+dims(1),:,:)))),'MPI transport matches serial fields')
  call require(maxval(abs(temperature-reference_temperature(first(1)+1:first(1)+dims(1),:,:)))<1e-8_dp, &
    'MPI transport matches serial temperatures')
  call require(maxval(abs(field-saved_field))>1e-6_dp,'transport changes the gas')
  call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],parcels,liquids,masses,local_totals,ok)
  call require(ok,'transport split initial budgets')
  call MPI_Allreduce(local_totals,before,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call advance_mpi_mixture_reaction(MPI_COMM_WORLD,species,reactions,transport,liquids,options, &
    field,temperature,parcels,masses,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1], &
    1e-7_dp,1e-7_dp,1e-12_dp,'implicit',steps,ok,'hllc','mc',[.true.,.true.,.true.,.true.])
  call require(ok,'MPI particles chemistry transport hydro split')
  call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],parcels,liquids,masses,local_totals,ok)
  call require(ok,'transport split final budgets')
  call MPI_Allreduce(local_totals,after,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(all(abs(after(:5)-before(:5))<2e-12_dp*max(1e-12_dp,abs(before(:5)))), &
    'MPI particles chemistry transport hydro conservation')
  if (rank==nranks-1) field(1,1,1,1)=-1
  saved_field=field;saved_temperature=temperature
  call advance_mpi_spray_hydro(MPI_COMM_WORLD,species,field,temperature,[1.0_dp,1.0_dp,1.0_dp], &
    cells,1e-5_dp,'hllc','mc',ok)
  call require(.not.ok.and.all(field==saved_field).and.all(temperature==saved_temperature), &
    'one rank invalid hydro state rolls back collectively')
  call advance_mpi_spray_transport(MPI_COMM_WORLD,species,transport,field,temperature,[.001_dp,.001_dp,.001_dp], &
    cells,1e-6_dp,[.true.,.true.,.true.,.true.],mpi_theta,ok)
  call require(.not.ok.and.all(field==saved_field).and.all(temperature==saved_temperature), &
    'one rank invalid transport state rolls back collectively')
  field=reference(:,first(1)+1:first(1)+dims(1),:,:)
  temperature=reference_temperature(first(1)+1:first(1)+dims(1),:,:)
  checkpoint_time=1e-8_dp;steps=3;checkpoint_event=2;last_id=huge(0_int64)-19
  write(checkpoint_path,'(a,i0,a)') 'mpi-spray-checkpoint-',nranks,'.chk'
  call mpi_spray_checkpoint(MPI_COMM_WORLD,.true.,trim(checkpoint_path),'mpi-checkpoint-context',species, &
    field,temperature,parcels,masses,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells, &
    checkpoint_time,steps,checkpoint_event,last_id,before,after,ok)
  call require(ok,'collective checkpoint write')
  saved_field=field;saved_temperature=temperature;saved=parcels;saved_mass=masses
  call mpi_spray_checkpoint(MPI_COMM_WORLD,.false.,trim(checkpoint_path),'wrong-context',species, &
    field,temperature,parcels,masses,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells, &
    checkpoint_time,steps,checkpoint_event,last_id,before,after,ok)
  call require(.not.ok.and.all(field==saved_field).and.all(masses==saved_mass).and.all(parcels%id==saved%id), &
    'checkpoint read failure preserves all ranks')
  ! Each rank now reads the same distributed checkpoint on a one-rank communicator.
  ! This tests a changed decomposition without requiring a second mpiexec process.
  field=reference*0;temperature=reference_temperature*0;checkpoint_time=0;steps=0;checkpoint_event=0;last_id=0
  deallocate(parcels,masses);allocate(parcels(0),masses(2,0))
  call mpi_spray_checkpoint(MPI_COMM_SELF,.false.,trim(checkpoint_path),'mpi-checkpoint-context',species, &
    field,temperature,parcels,masses,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells, &
    checkpoint_time,steps,checkpoint_event,last_id,before,after,ok)
  call require(ok.and.all(field==reference).and.all(temperature==reference_temperature), &
    'changed-rank checkpoint exact gas recovery')
  call require(size(parcels)==nranks*(nranks-1)/2.and.last_id==huge(0_int64)-19, &
    'changed-rank particle collection and ID reservation recovery')
  call require(checkpoint_time==1e-8_dp.and.steps==3.and.checkpoint_event==2,'checkpoint metadata recovery')
  do p=1,size(parcels)
    call require(parcels(p)%id==p,'checkpoint global particle ordering')
    call require(abs(sum(masses(:,p))-parcels(p)%mass)<1e-24_dp,'checkpoint component recovery')
  end do
  do p=1,size(saved)
    i=int(saved(p)%id)
    call require(all(parcels(i)%position==saved(p)%position).and.all(parcels(i)%velocity==saved(p)%velocity).and. &
      parcels(i)%mass==saved(p)%mass.and.parcels(i)%temperature==saved(p)%temperature.and. &
      parcels(i)%multiplicity==saved(p)%multiplicity.and.all(masses(:,i)==saved_mass(:,p)), &
      'changed-rank checkpoint exact local particle payload')
  end do
  field=saved_field*0;temperature=saved_temperature*0
  deallocate(parcels,masses);allocate(parcels(0),masses(2,0))
  call mpi_spray_checkpoint(MPI_COMM_WORLD,.false.,trim(checkpoint_path),'mpi-checkpoint-context',species, &
    field,temperature,parcels,masses,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells, &
    checkpoint_time,steps,checkpoint_event,last_id,before,after,ok)
  call require(ok.and.all(field==saved_field).and.all(temperature==saved_temperature), &
    'checkpoint scatters exact owner blocks')
  call require(all(parcels%id==saved%id).and.all(masses==saved_mass),'checkpoint redistributes owner particles')
  call check_film_checkpoint()
  call check_retired_parent_checkpoint()
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (rank==0) then
    open(newunit=checkpoint_unit,file=trim(checkpoint_path),status='old');close(checkpoint_unit,status='delete')
  end if
  call MPI_Finalize(ierr)
contains
  subroutine check_boundary_split()
    type(spray_parcel), allocatable :: drops(:)
    type(spray_options) :: settings
    real(dp), allocatable :: u(:,:,:,:),t(:,:,:),u0(:,:,:,:),t0(:,:,:),cm(:,:)
    real(dp) :: ledger(15),local_before(15),local_after(15),total_before(15),total_after(15),wall_t(6)
    character(len=32) :: boundary(6)
    integer :: np,nstep
    integer(int64) :: mark
    u=reference(:,first(1)+1:first(1)+dims(1),:,:)
    t=reference_temperature(first(1)+1:first(1)+dims(1),:,:);u0=u;t0=t
    np=merge(1,0,rank==0);allocate(drops(np),cm(2,np));mark=17
    if (rank==0) then
      drops(1)=spray_parcel(1_int64,[1e-6_dp,.5_dp,.5_dp],[-1.0_dp,0.0_dp,0.0_dp],1e-9_dp,330.0_dp,100.0_dp)
      cm(:,1)=drops(1)%mass*[.25_dp,.75_dp]
    end if
    settings%drag=.false.;settings%heat=.false.;settings%evaporation=.false.
    boundary='periodic';boundary(1:2)='outflow';boundary(3)='no_slip_wall';boundary(4)='isothermal_wall'
    wall_t=600
    if (rank==nranks-1) wall_t(4)=-1
    call advance_mpi_mixture_reaction(MPI_COMM_WORLD,species,reactions,transport,liquids,settings,u,t,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],1e-5_dp,1e-7_dp,1e-12_dp, &
      'implicit',nstep,ok,'hllc','mc',[.true.,.true.,.true.,.true.],mark,[.1_dp,.7_dp,.7_dp], &
      boundary,ledger,wall_temperature=wall_t)
    call require(.not.ok.and.nstep==0.and.all(u==u0).and.all(t==t0).and.mark==17.and.all(ledger==0), &
      'late boundary transport failure rolls back whole split gas ID and ledger')
    call require(size(drops)==np,'late boundary failure restores escaped particle')
    if (rank==0) call require(drops(1)%position(1)==1e-6_dp.and.drops(1)%mass==1e-9_dp.and. &
      all(cm(:,1)==1e-9_dp*[.25_dp,.75_dp]),'late boundary failure restores parcel inventory')
    call mixture_system_totals(u,[1.0_dp,1.0_dp,1.0_dp],drops,liquids,cm,local_before,ok)
    call require(ok,'boundary split initial inventory')
    call MPI_Allreduce(local_before,total_before,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    wall_t=600
    call advance_mpi_mixture_reaction(MPI_COMM_WORLD,species,reactions,transport,liquids,settings,u,t,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],1e-5_dp,1e-7_dp,1e-12_dp, &
      'implicit',nstep,ok,'hllc','mc',[.true.,.true.,.true.,.true.],mark,[.1_dp,.7_dp,.7_dp], &
      boundary,ledger,wall_temperature=wall_t)
    call require(ok.and.nstep>=2.and.size(drops)==0.and.mark==17,'complete nonperiodic split with escaped parcel')
    call mixture_system_totals(u,[1.0_dp,1.0_dp,1.0_dp],drops,liquids,cm,local_after,ok)
    call require(ok,'boundary split final inventory')
    call MPI_Allreduce(local_after,total_after,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(abs(total_after-total_before-ledger)<2e-12_dp*max(1e-12_dp,abs(total_before))), &
      'whole nonperiodic split inventory follows summed surface and particle ledger')
    call require(maxval(abs(u-u0))>1e-5_dp,'nonperiodic split advances gas')
  end subroutine

  subroutine check_boundary_ghosts()
    real(dp), allocatable :: global_ghost(:,:,:,:),global_t(:,:,:),local_ghost(:,:,:,:),local_t(:,:,:),u(:,:,:,:),t(:,:,:)
    real(dp), allocatable :: saved_ghost(:,:,:,:),saved_t(:,:,:)
    real(dp), allocatable :: serial_u(:,:,:,:),serial_t(:,:,:),local_u(:,:,:,:),local_temp(:,:,:)
    real(dp) :: velocity(3,6),inlet(15,6),inlet_t(6),serial_change(15),change(15),delta(15),total_delta(15),inventory(15)
    real(dp) :: theta_serial,theta_mpi
    real(dp) :: wall_t(6)
    real(dp), allocatable :: transport_initial(:,:,:,:)
    character(len=32) :: boundaries(6)
    real(dp), parameter :: spacing(3)=[1.0_dp,1.0_dp,1.0_dp]
    integer :: mode
    allocate(global_ghost(15,cells(1)+4,7,7),global_t(cells(1)+4,7,7))
    u=reference(:,first(1)+1:first(1)+dims(1),:,:)
    t=reference_temperature(first(1)+1:first(1)+dims(1),:,:)
    inlet=spread(reference(:,1,1,1),2,6);inlet_t=reference_temperature(1,1,1)
    do mode=1,3
      boundaries='periodic';velocity=0
      if (mode/=1) then
        boundaries(1)='slip_wall';boundaries(2)='outflow'
        boundaries(3)='no_slip_wall';boundaries(4)='isothermal_wall'
        velocity(1,3)=2;velocity(1,4)=-1
      end if
      if (mode==3) boundaries(1)='inflow'
      call fill_reactive_boundary_ghosts_3d(reference,reference_temperature,boundaries,global_ghost,global_t,ok, &
        viscous_wall=.true.,inflow_state=inlet,inflow_temperature=inlet_t,wall_velocity=velocity)
      call require(ok,'serial boundary ghost reference')
      call fill_mpi_spray_ghosts(MPI_COMM_WORLD,u,t,cells,boundaries,local_ghost,local_t,ok, &
        viscous_wall=.true.,inflow_state=inlet,inflow_temperature=inlet_t,wall_velocity=velocity)
      call require(ok,'MPI boundary ghost exchange')
      call require(all(local_ghost==global_ghost(:,first(1)+1:first(1)+dims(1)+4,:,:)), &
        'MPI ghost states exactly match serial including corners')
      call require(all(local_t==global_t(first(1)+1:first(1)+dims(1)+4,:,:)), &
        'MPI ghost temperatures exactly match serial')
      serial_u=reference;serial_t=reference_temperature;local_u=u;local_temp=t
      call advance_reactive_boundary_hydro_3d(species,serial_u,serial_t,spacing,1e-5_dp,boundaries, &
        'characteristic_plm','mc','hllc',serial_change,ok,inlet,inlet_t)
      call require(ok,'serial physical-boundary hydro reference')
      call advance_mpi_spray_hydro(MPI_COMM_WORLD,species,local_u,local_temp,spacing,cells,1e-5_dp, &
        'hllc','mc',ok,boundaries,change,inlet,inlet_t)
      call require(ok,'MPI physical-boundary hydro')
      call require(maxval(abs(local_u-serial_u(:,first(1)+1:first(1)+dims(1),:,:))/ &
        max(1.0_dp,abs(serial_u(:,first(1)+1:first(1)+dims(1),:,:))))<2e-12_dp, &
        'MPI physical hydro matches serial states')
      call require(maxval(abs(local_temp-serial_t(first(1)+1:first(1)+dims(1),:,:)))<1e-5_dp, &
        'MPI physical hydro matches serial temperature')
      call require(maxval(abs(change-serial_change)/max(1e-12_dp,abs(serial_change)))<1e-10_dp, &
        'MPI physical hydro global boundary ledger matches serial')
      delta=sum(sum(sum(local_u-u,dim=4),dim=3),dim=2)*product(spacing)
      call MPI_Allreduce(delta,total_delta,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      inventory=sum(sum(sum(abs(reference),dim=4),dim=3),dim=2)*product(spacing)
      call require(all(abs(total_delta-change)<2e-12_dp*max(1e-12_dp,inventory)), &
        'MPI boundary flux accounts for total inventory change')
      serial_u=reference;serial_t=reference_temperature;wall_t=600
      do k=1,3
        do j=1,3
          do i=1,cells(1)
            primitive=0;primitive(1)=1;primitive(2:4)=[2.0_dp*j,-1.0_dp*k,.5_dp*i]
            primitive(6)=.1_dp+.08_dp*sin(2*acos(-1.0_dp)*real(i-1,dp)/cells(1))
            primitive(9)=.2_dp;primitive(15)=1-sum(primitive(6:14))
            primitive(5)=800*mixture_specific_gas_constant(species,primitive(6:15),ok)
            call reactive_primitive_to_conserved(species,primitive,serial_u(:,i,j,k),serial_t(i,j,k),sound,ok)
            call require(ok,'multispecies physical transport initial state')
          end do
        end do
      end do
      local_u=serial_u(:,first(1)+1:first(1)+dims(1),:,:)
      local_temp=serial_t(first(1)+1:first(1)+dims(1),:,:);transport_initial=local_u
      inventory=sum(sum(sum(abs(serial_u),dim=4),dim=3),dim=2)
      call advance_reactive_transport_3d(species,transport,serial_u,serial_t,cells(1),3,3, &
        .001_dp,.001_dp,.001_dp,1e-5_dp,.true.,.true.,.true.,.true.,theta_serial,ok, &
        boundaries,serial_change,inlet,inlet_t,wall_t,velocity)
      call require(ok,'serial physical-boundary transport reference')
      call advance_mpi_spray_transport(MPI_COMM_WORLD,species,transport,local_u,local_temp, &
        [.001_dp,.001_dp,.001_dp],cells,1e-5_dp,[.true.,.true.,.true.,.true.],theta_mpi,ok, &
        boundaries,change,inlet,inlet_t,wall_t,velocity)
      call require(ok,'MPI physical-boundary transport')
      call require(maxval(abs(local_u-serial_u(:,first(1)+1:first(1)+dims(1),:,:))/ &
        max(1.0_dp,abs(serial_u(:,first(1)+1:first(1)+dims(1),:,:))))<2e-12_dp, &
        'MPI physical transport matches serial states')
      call require(maxval(abs(local_temp-serial_t(first(1)+1:first(1)+dims(1),:,:)))<1e-8_dp, &
        'MPI physical transport matches serial temperature')
      call require(all(abs(change-serial_change)<1e-10_dp*max(1e-12_dp,abs(serial_change))), &
        'MPI transport boundary ledger matches serial')
      delta=sum(sum(sum(local_u-transport_initial,dim=4),dim=3),dim=2)*1e-9_dp
      call MPI_Allreduce(delta,total_delta,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call require(all(abs(total_delta-change)<2e-12_dp*max(1e-12_dp,inventory*1e-9_dp)), &
        'MPI transport boundary flux accounts for inventory change')
      transport_initial=local_u
      inventory=sum(sum(sum(abs(serial_u),dim=4),dim=3),dim=2)
      call advance_les_ssprk2_3d(species,serial_u,serial_t,[.001_dp,.001_dp,.001_dp],.1_dp,.7_dp,.7_dp, &
        1e-5_dp,ok,boundaries,serial_change,inlet,inlet_t,wall_t,velocity)
      call require(ok,'serial physical-boundary SGS reference')
      call advance_mpi_spray_les(MPI_COMM_WORLD,species,local_u,local_temp,[.001_dp,.001_dp,.001_dp],cells, &
        1e-5_dp,.1_dp,.7_dp,.7_dp,ok,boundaries,change,inlet,inlet_t,wall_t,velocity)
      call require(ok,'MPI physical-boundary SGS')
      call require(maxval(abs(local_u-serial_u(:,first(1)+1:first(1)+dims(1),:,:))/ &
        max(1.0_dp,abs(serial_u(:,first(1)+1:first(1)+dims(1),:,:))))<2e-12_dp, &
        'MPI physical SGS matches serial states')
      call require(maxval(abs(local_temp-serial_t(first(1)+1:first(1)+dims(1),:,:)))<1e-8_dp, &
        'MPI physical SGS matches serial temperature')
      call require(all(abs(change-serial_change)<1e-10_dp*max(1e-12_dp,abs(serial_change))), &
        'MPI SGS boundary ledger matches serial')
      delta=sum(sum(sum(local_u-transport_initial,dim=4),dim=3),dim=2)*1e-9_dp
      call MPI_Allreduce(delta,total_delta,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call require(all(abs(total_delta-change)<2e-12_dp*max(1e-12_dp,inventory*1e-9_dp)), &
        'MPI SGS boundary flux accounts for inventory change')
    end do
    local_u=u;local_temp=t
    call advance_mpi_spray_hydro(MPI_COMM_WORLD,species,local_u,local_temp,spacing,cells,1e3_dp, &
      'hllc','mc',ok,boundaries,change,inlet,inlet_t)
    call require(.not.ok.and.all(local_u==u).and.all(local_temp==t).and.all(change==0), &
      'MPI physical hydro failed stage rolls back state and boundary ledger')
    local_u=u;local_temp=t;wall_t=600
    if (rank==0) wall_t(4)=-1
    call advance_mpi_spray_transport(MPI_COMM_WORLD,species,transport,local_u,local_temp, &
      [.001_dp,.001_dp,.001_dp],cells,1e-5_dp,[.true.,.true.,.true.,.true.],theta_mpi,ok, &
      boundaries,change,inlet,inlet_t,wall_t,velocity)
    call require(.not.ok.and.all(local_u==u).and.all(local_temp==t).and.all(change==0), &
      'MPI invalid wall thermal data rolls back transport')
    call advance_mpi_spray_les(MPI_COMM_WORLD,species,local_u,local_temp,[.001_dp,.001_dp,.001_dp],cells, &
      1e-5_dp,.1_dp,.7_dp,.7_dp,ok,boundaries,change,inlet,inlet_t,wall_t,velocity)
    call require(.not.ok.and.all(local_u==u).and.all(local_temp==t).and.all(change==0), &
      'MPI invalid wall thermal data rolls back SGS')
    saved_ghost=local_ghost;saved_t=local_t
    if (rank==0) boundaries(2)='slip_wall'
    call fill_mpi_spray_ghosts(MPI_COMM_WORLD,u,t,cells,boundaries,local_ghost,local_t,ok, &
      viscous_wall=.true.,inflow_state=inlet,inflow_temperature=inlet_t,wall_velocity=velocity)
    call require(.not.ok.and.all(local_ghost==saved_ghost).and.all(local_t==saved_t), &
      'MPI boundary disagreement rejects without changing ghost output')
    local_u=u;local_temp=t
    call advance_mpi_spray_hydro(MPI_COMM_WORLD,species,local_u,local_temp,spacing,cells,1e-5_dp, &
      'hllc','mc',ok,boundaries,change,inlet,inlet_t)
    call require(.not.ok.and.all(local_u==u).and.all(local_temp==t).and.all(change==0), &
      'MPI physical hydro face disagreement rolls back')
  end subroutine

  subroutine check_coupled_impact()
    type(spray_parcel), allocatable :: drops(:),original(:)
    type(spray_wall_film), allocatable :: films(:)
    type(spray_options) :: settings
    type(spray_wall_options) :: walls
    real(dp), allocatable :: cm(:,:),fm(:,:),gas(:,:,:,:),temp(:,:,:)
    real(dp) :: local(15),initial(15),final(15),ledger(15),scale(15)
    character(len=32) :: boundary(6)
    integer :: np,nstep
    integer(int64) :: mark
    np=merge(1,0,rank<nranks-1);allocate(drops(np),cm(2,np),films(0),fm(2,0));mark=int(nranks,int64)
    if (np>0) then
      drops(1)=spray_parcel(int(rank+1,int64),[real(cells(1),dp)-.00025_dp,.25_dp,.001_dp], &
        [1.0_dp,0.0_dp,-1.0_dp],1e-9_dp,330.0_dp,100.0_dp)
      cm(:,1)=drops(1)%mass*[.25_dp,.75_dp]
    end if
    original=drops;gas=field;temp=temperature
    settings%drag=.false.;settings%heat=.false.;settings%evaporation=.false.
    settings%surface_tension=.072_dp;settings%liquid_viscosity=.001_dp;settings%maximum_substeps=1
    boundary='periodic';boundary(5:6)='slip_wall';walls%enabled(5)=.true.
    call advance_mpi_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,settings,gas,temp,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],.5_dp,nstep,ok,mark, &
      boundary,ledger,films,fm,walls)
    call require(.not.ok.and.nstep==0.and.size(films)==0.and.mark==nranks.and.all(ledger==0), &
      'coupled impact exhaustion rolls back films IDs and ledger')
    call require(all(gas==field).and.all(temp==temperature).and.all(drops%mass==original%mass), &
      'coupled impact exhaustion rolls back gas and parent mass')
    if (np>0) call require(all(drops(1)%position==original(1)%position).and. &
      all(cm(:,1)==1e-9_dp*[.25_dp,.75_dp]),'coupled impact exhaustion restores parent position/composition')
    settings%maximum_substeps=100
    call mixture_system_totals(gas,[1.0_dp,1.0_dp,1.0_dp],drops,liquids,cm,local,ok,films,fm)
    call require(ok,'coupled impact initial inventory')
    call MPI_Allreduce(local,initial,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call advance_mpi_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,settings,gas,temp,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],.002_dp,nstep,ok,mark, &
      boundary,ledger,films,fm,walls)
    call require(ok.and.nstep==1.and.size(drops)==0.and.mark==nranks,'coupled deposition consumes parents without new IDs')
    call require(size(films)==merge(nranks-1,0,rank==0),'new films migrate to spatial owner')
    call require(all(abs(fm(1,:)-.25e-7_dp)<1e-20_dp).and.all(abs(fm(2,:)-.75e-7_dp)<1e-20_dp), &
      'deposited film inherits whole-parcel component inventories')
    call mixture_system_totals(gas,[1.0_dp,1.0_dp,1.0_dp],drops,liquids,cm,local,ok,films,fm)
    call require(ok,'coupled impact final inventory')
    call MPI_Allreduce(local,final,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(abs(final-initial-ledger)<2e-12_dp*max(1e-12_dp,abs(initial))), &
      'coupled gas particle film inventory includes wall exchange exactly once')
    drops=original;deallocate(cm,films,fm);allocate(cm(2,np),films(0),fm(2,0));gas=field;temp=temperature
    mark=int(nranks,int64)
    if (np>0) then
      drops(1)%position(1)=.25_dp;drops(1)%velocity(3)=-40;cm(:,1)=drops(1)%mass*[.25_dp,.75_dp]
    end if
    call mixture_system_totals(gas,[1.0_dp,1.0_dp,1.0_dp],drops,liquids,cm,local,ok,films,fm)
    call require(ok,'coupled splash initial inventory')
    call MPI_Allreduce(local,initial,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call advance_mpi_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,settings,gas,temp,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],1e-4_dp,nstep,ok,mark, &
      boundary,ledger,films,fm,walls)
    call require(ok.and.mark==nranks+4*(nranks-1),'coupled splash reserves global child IDs')
    call require(size(drops)==merge(4*(nranks-1),0,rank==0),'splash children migrate to owner')
    call require(all(drops%parent_id>=1).and.all(drops%parent_id<nranks),'splash lineage survives migration')
    call mixture_system_totals(gas,[1.0_dp,1.0_dp,1.0_dp],drops,liquids,cm,local,ok,films,fm)
    call require(ok,'coupled splash final inventory')
    call MPI_Allreduce(local,final,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    scale=max(1e-12_dp,abs(initial));scale(2:4)=max(1e-12_dp,norm2(initial(2:4)))
    call require(all(abs(final-initial-ledger)<2e-12_dp*scale), &
      'coupled splash gas particle film and wall budgets close')
  end subroutine

  subroutine check_particle_boundaries()
    type(spray_parcel), allocatable :: drops(:),initial(:)
    type(spray_options) :: settings
    real(dp), allocatable :: cm(:,:),initial_mass(:,:)
    real(dp) :: u(15,1,2,2),t(1,2,2),u0(15,1,2,2),t0(1,2,2),q(15),cell(15),tc,ac,ledger(15)
    character(len=32) :: boundary(6)
    integer :: np,nstep,c
    np=merge(3,0,rank==0);allocate(drops(np),cm(2,np))
    do c=1,np
      drops(c)=spray_parcel(int(c,int64),[.5_dp,.5_dp,.5_dp],[0.0_dp,0.0_dp,0.0_dp],1e-9_dp,330.0_dp,100.0_dp)
      cm(:,c)=drops(c)%mass*[.25_dp,.75_dp]
    end do
    if (rank==0) then
      drops(1)%position(1)=.01_dp;drops(1)%velocity(1)=-2
      drops(2)%position(3)=.01_dp;drops(2)%velocity(3)=-2
      drops(3)%position(1)=.99_dp;drops(3)%velocity(1)=2
    end if
    q=0;q(1)=1;q(5)=101325;q(15)=1
    call reactive_primitive_to_conserved(species,q,cell,tc,ac,ok)
    call require(ok,'particle boundary gas initialization')
    u=spread(spread(spread(cell,2,1),3,2),4,2);t=tc;u0=u;t0=t;initial=drops;initial_mass=cm
    boundary='periodic';boundary(1:2)='outflow';boundary(5:6)='slip_wall'
    settings%drag=.false.;settings%heat=.false.;settings%evaporation=.false.;settings%maximum_substeps=1
    call advance_mpi_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,settings,u,t,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],[nranks,2,2],[nranks,1,1],.8_dp,nstep,ok, &
      boundary_faces=boundary,boundary_change=ledger)
    call require(.not.ok.and.all(u==u0).and.all(t==t0).and.all(ledger==0).and.all(cm==initial_mass), &
      'nonperiodic particle substep exhaustion rolls back gas mass and ledger')
    do c=1,np
      call require(all(drops(c)%position==initial(c)%position).and.all(drops(c)%velocity==initial(c)%velocity), &
        'nonperiodic particle substep exhaustion rolls back positions and velocities')
    end do
    settings%maximum_substeps=100
    call advance_mpi_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,settings,u,t,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],[nranks,2,2],[nranks,1,1],.01_dp,nstep,ok, &
      boundary_faces=boundary,boundary_change=ledger)
    call require(ok.and.nstep==1,'nonperiodic particle coupling')
    call require(size(drops)==merge(1,0,rank<2),'outflow removed and survivors migrated')
    if (rank==0) then
      call require(drops(1)%id==2.and.abs(drops(1)%position(3)-.01_dp)<1e-14_dp.and.drops(1)%velocity(3)==2, &
        'stationary wall reflection matches ballistic path')
    else if (rank==1) then
      call require(drops(1)%id==3.and.abs(drops(1)%position(1)-1.01_dp)<1e-14_dp, &
        'nonperiodic survivor crosses owner interface')
    end if
    call require(abs(ledger(1)+1e-7_dp)<1e-20_dp.and.abs(ledger(2)-2e-7_dp)<1e-20_dp.and. &
      abs(ledger(4)-4e-7_dp)<1e-20_dp,'outflow mass and reflection impulse ledger')
    call require(abs(ledger(5+liquids(1)%vapor_index)+.25e-7_dp)<1e-20_dp.and. &
      abs(ledger(5+liquids(2)%vapor_index)+.75e-7_dp)<1e-20_dp,'outflow component ledger')
    tc=.25_dp*liquid_energy(liquids(1),330.0_dp)+.75_dp*liquid_energy(liquids(2),330.0_dp)+2.0_dp
    call require(abs(ledger(5)+1e-7_dp*tc)<1e-12_dp*max(1e-12_dp,abs(1e-7_dp*tc)), &
      'outflow energy includes both component internal energies and kinetic energy')
  end subroutine

  subroutine check_wall_crossing()
    type(spray_parcel), allocatable :: drops(:),children(:)
    type(spray_parcel) :: expected_parent
    type(spray_wall_film), allocatable :: patches(:),empty_films(:)
    type(spray_wall_film) :: expected_film
    type(spray_options) :: settings
    type(spray_wall_options) :: walls
    real(dp), allocatable :: cm(:,:),fm(:,:),previous(:,:),press(:),child_mass(:,:)
    logical, allocatable :: hits(:)
    real(dp) :: mw(2),parent_mass(2),film_components(2),wall(5),serial_wall(5),global_wall(5)
    logical :: serial_hit
    integer :: np,q
    integer(int64) :: mark
    character(len=16) :: face_names(6)
    np=merge(2,0,rank<nranks-1)
    allocate(drops(np),cm(2,np),fm(2,0),patches(0),empty_films(0),previous(3,np),press(np),hits(1))
    hits=.true.;press=101325;mark=int(2*nranks,int64)
    do q=1,np
      drops(q)=spray_parcel(int(2*rank+q,int64),[.25_dp,.25_dp,-.1_dp], &
        [.1_dp,0.0_dp,-1.0_dp],1e-9_dp,330.0_dp,100.0_dp)
      previous(:,q)=[.25_dp,.25_dp,.1_dp];cm(:,q)=drops(q)%mass*[.25_dp,.75_dp]
    end do
    if (np>0) drops(2)%position(3)=.05_dp
    mw=[species(liquids(1)%vapor_index)%molecular_weight,species(liquids(2)%vapor_index)%molecular_weight]
    settings%surface_tension=.072_dp;settings%liquid_viscosity=.001_dp
    face_names='periodic';face_names(5:6)='slip_wall';walls%enabled(5)=.true.
    if (rank==0) previous(3,1)=-.01_dp
    call impact_mpi_planar_walls(MPI_COMM_WORLD,liquids,settings,walls,face_names,mw,press,previous, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],drops,cm,patches,fm,mark,hits,wall,ok)
    call require(.not.ok.and.size(hits)==1.and.all(hits).and.mark==2*nranks,'invalid previous position preserves hit/ID')
    if (rank==0) previous(3,1)=.1_dp
    serial_wall=0
    if (np>0) then
      expected_parent=drops(1);parent_mass=cm(:,1);film_components=0
      call mixture_impact_at_planar_wall(expected_parent,previous(:,1),liquids,parent_mass,settings,mw,press(1), &
        [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],walls,empty_films,mark+1,children,child_mass, &
        expected_film,film_components,serial_wall,serial_hit,ok)
      call require(ok.and.serial_hit,'serial wall crossing reference')
    end if
    call MPI_Allreduce(serial_wall,global_wall,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call impact_mpi_planar_walls(MPI_COMM_WORLD,liquids,settings,walls,face_names,mw,press,previous, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],drops,cm,patches,fm,mark,hits,wall,ok)
    call require(ok.and.size(hits)==np.and.mark==2*nranks,'MPI wall crossing accepts without unused splash IDs')
    call require(all(abs(wall-global_wall)<1e-13_dp*max(1e-12_dp,abs(global_wall))),'MPI crossing matches serial wall budget')
    if (np>0) then
      call require(hits(1).and..not.hits(2).and.size(patches)==1,'only crossed particle deposits')
      call require(drops(1)%mass==expected_parent%mass.and.all(cm(:,1)==parent_mass),'MPI parent matches serial impact')
      call require(all(patches(1)%position==expected_film%position).and.patches(1)%height==expected_film%height, &
        'MPI film crossing geometry matches serial')
      call require(all(fm(:,1)==film_components),'MPI deposited components match serial')
      call require(drops(2)%mass==1e-9_dp.and.drops(2)%position(3)==.05_dp,'uncrossed particle unchanged')
    end if
  end subroutine

  subroutine check_film_query()
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: pos(:,:),normal(:,:),thickness(:),saved_thickness(:)
    real(dp) :: expected,tolerance
    integer :: np,nq
    np=merge(2,0,rank<nranks-1);nq=merge(4,0,rank>0)
    allocate(patches(np),pos(3,nq),normal(3,nq),thickness(1));thickness=123
    if (np>0) then
      patches(1)=spray_wall_film([.25_dp,.25_dp,0.0_dp],[0.0_dp,0.0_dp,1.0_dp], &
        1e-9_dp,330.0_dp,.2_dp,real(rank+1,dp)*1e-6_dp,1)
      patches(2)=patches(1);patches(2)%normal=-patches(1)%normal;patches(2)%height=.01_dp
    end if
    if (nq>0) then
      pos(:,1)=[.25_dp,.25_dp,0.0_dp];pos(:,2)=[.75_dp,.25_dp,0.0_dp]
      pos(:,3)=[.25_dp,.25_dp,1e-5_dp];pos(:,4)=pos(:,1)
      normal=0;normal(3,:)=1;normal(3,4)=-1
    end if
    call sample_mpi_wall_film_thickness(MPI_COMM_WORLD,patches,pos,normal,1e-8_dp,thickness,ok)
    call require(ok.and.size(thickness)==nq,'MPI film query supports empty sources and owners')
    if (nq>0) then
      expected=real(nranks*(nranks-1)/2,dp)*1e-6_dp
      call require(abs(thickness(1)-expected)<1e-14_dp*expected,'MPI overlap sums remote film heights')
      call require(thickness(2)==0.and.thickness(3)==0,'MPI film footprint and plane exclusion')
      call require(abs(thickness(4)-.01_dp*(nranks-1))<1e-15_dp,'MPI film normal filtering')
    end if
    saved_thickness=thickness;tolerance=1e-8_dp
    if (rank==0) tolerance=2e-8_dp
    call sample_mpi_wall_film_thickness(MPI_COMM_WORLD,patches,pos,normal,tolerance,thickness,ok)
    call require(.not.ok.and.all(thickness==saved_thickness),'MPI film tolerance disagreement preserves output')
    if (rank==1) normal(:,1)=0
    call sample_mpi_wall_film_thickness(MPI_COMM_WORLD,patches,pos,normal,1e-8_dp,thickness,ok)
    call require(.not.ok.and.all(thickness==saved_thickness),'MPI invalid query rejects collectively')
  end subroutine

  subroutine check_wall_impacts()
    type(spray_parcel), allocatable :: drops(:),original(:)
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: cm(:,:),original_cm(:,:),fm(:,:),loc(:,:),normal(:,:),data(:,:)
    integer, allocatable :: events(:)
    real(dp) :: wall(5),initial(15),final(15),global_initial(15),global_final(15)
    integer :: np,q
    integer(int64) :: mark
    np=merge(1,0,rank<nranks-1)
    allocate(drops(np),cm(2,np),patches(0),fm(2,0),loc(3,np),normal(3,np),data(6,np),events(np))
    if (np>0) then
      drops(1)=spray_parcel(int(rank+1,int64),[.25_dp,.25_dp,0.0_dp], &
        [10.0_dp,0.0_dp,-20.0_dp],1e-9_dp,330.0_dp,100.0_dp)
      cm(:,1)=drops(1)%mass*[.25_dp,.75_dp];events(1)=1
      loc(:,1)=drops(1)%position;normal(:,1)=[0.0_dp,0.0_dp,1.0_dp]
      data(:,1)=[150.0_dp,1.0_dp,80.0_dp,0.0_dp,100.0_dp,0.0_dp]
    end if
    original=drops;original_cm=cm;mark=int(nranks,int64)
    call resolve_mpi_wall_impacts(MPI_COMM_WORLD,liquids,drops,cm,events,loc,normal,data,acos(-1.0_dp)/3, &
      patches,fm,100,100,mark,wall,ok)
    call require(.not.ok.and.mark==nranks.and.all(wall==0),'late impact failure rolls back IDs and wall ledger')
    call require(all(cm==original_cm).and.all(drops%mass==original%mass).and.size(patches)==0, &
      'late impact failure rolls back particle and film candidates')
    data(6,:)=1000
    call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],drops,liquids,cm,initial,ok,patches,fm)
    call require(ok,'wall event initial budget')
    call MPI_Allreduce(initial,global_initial,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call resolve_mpi_wall_impacts(MPI_COMM_WORLD,liquids,drops,cm,events,loc,normal,data,acos(-1.0_dp)/3, &
      patches,fm,100,100,mark,wall,ok)
    call require(ok.and.mark==nranks+4*(nranks-1),'MPI impact reserves child IDs')
    call require(size(drops)==5*np.and.size(patches)==np,'MPI dry impact splashes and deposits')
    do q=2,size(drops)
      call require(drops(q)%id==nranks+4*rank+q-1.and.drops(q)%parent_id==rank+1,'splash IDs and lineage')
      call require(all(abs(cm(:,q)/drops(q)%mass-[.25_dp,.75_dp])<1e-14_dp),'splash composition inheritance')
    end do
    do q=1,2
      call require(abs(sum(cm(q,:)*drops%multiplicity)+sum(fm(q,:))-sum(original_cm(q,:)*original%multiplicity)) &
        <2e-12_dp*max(1e-30_dp,sum(original_cm(q,:)*original%multiplicity)),'impact component conservation')
    end do
    call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],drops,liquids,cm,final,ok,patches,fm)
    call require(ok,'wall event final budget')
    call MPI_Allreduce(final,global_final,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    global_final(:5)=global_final(:5)+wall
    call require(all(abs(global_final-global_initial)<2e-12_dp*max(1e-12_dp,abs(global_initial))), &
      'wall event mass momentum energy species conservation')
  end subroutine

  subroutine check_retired_parent_checkpoint()
    type(spray_parcel), allocatable :: drops(:)
    type(spray_wall_film), allocatable :: films(:)
    real(dp), allocatable :: cm(:,:),fm(:,:),gas(:,:,:,:),temp(:,:,:)
    real(dp) :: ref(15),inj(15),wall(15),time
    integer :: steps,event,unit,np
    integer(int64) :: mark
    character(len=100) :: path
    write(path,'(a,i0,a)') 'mpi_retired_parent_',nranks,'.chk'
    np=merge(1,0,rank==0);allocate(drops(np),cm(2,np),films(0),fm(2,0))
    if (np>0) then
      drops(1)=spray_parcel(2_int64,[.25_dp,.25_dp,.25_dp],[0.0_dp,0.0_dp,0.0_dp],1e-9_dp,330.0_dp,100.0_dp)
      drops(1)%parent_id=1;cm(:,1)=drops(1)%mass*[.25_dp,.75_dp]
    end if
    gas=field;temp=temperature;ref=0;inj=0;wall=0;time=0;steps=0;event=0;mark=2
    call mpi_spray_checkpoint(MPI_COMM_WORLD,.true.,trim(path),'retired-parent-test',species,gas,temp,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,time,steps,event,mark,ref,inj,ok, &
      wall_films=films,film_mass=fm,wall_totals=wall)
    call require(.not.ok,'strict checkpoint still rejects a missing parent')
    if (np>0) drops(1)%parent_id=2
    call mpi_spray_checkpoint(MPI_COMM_WORLD,.true.,trim(path),'retired-parent-test',species,gas,temp,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,time,steps,event,mark,ref,inj,ok, &
      wall_films=films,film_mass=fm,wall_totals=wall,retired_parent_state=.true.)
    call require(.not.ok,'retired lineage still rejects self parent')
    if (np>0) drops(1)%parent_id=1
    call mpi_spray_checkpoint(MPI_COMM_WORLD,.true.,trim(path),'retired-parent-test',species,gas,temp,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,time,steps,event,mark,ref,inj,ok, &
      wall_films=films,film_mass=fm,wall_totals=wall,retired_parent_state=.true.)
    call require(ok,'retired lineage checkpoint write')
    call mpi_spray_checkpoint(MPI_COMM_WORLD,.false.,trim(path),'retired-parent-test',species,gas,temp,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,time,steps,event,mark,ref,inj,ok, &
      wall_films=films,film_mass=fm,wall_totals=wall)
    call require(.not.ok,'retired lineage schema cannot be read as strict parent schema')
    call mpi_spray_checkpoint(MPI_COMM_WORLD,.false.,trim(path),'retired-parent-test',species,gas,temp,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,time,steps,event,mark,ref,inj,ok, &
      wall_films=films,film_mass=fm,wall_totals=wall,retired_parent_state=.true.)
    call require(ok.and.mark==2.and.all(drops%parent_id==1),'retired lineage survives read')
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if (rank==0) then
      open(newunit=unit,file=trim(path),status='old');close(unit,status='delete')
    end if
  end subroutine

  subroutine check_film_checkpoint()
    type(spray_wall_film), allocatable :: patches(:)
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: cm(:,:),pm(:,:),gas(:,:,:,:),temp(:,:,:)
    real(dp) :: wall(15),ref(15),inj(15),time
    integer :: step,event,q,unit
    integer(int64) :: mark
    character(len=80) :: path
    allocate(patches(1),cm(2,1),drops(0),pm(2,0))
    patches(1)=spray_wall_film([real(first(1),dp)+.5_dp,0.0_dp,.5_dp],[0.0_dp,1.0_dp,0.0_dp], &
      1e-9_dp,330.0_dp+rank,.001_dp,1e-6_dp,rank+1)
    cm(:,1)=patches(1)%mass*[.25_dp,.75_dp]
    gas=field;temp=temperature;wall=0;wall(5)=1.25_dp;ref=0;inj=0;time=1e-8_dp;step=3;event=2;mark=19
    write(path,'(a,i0,a)') 'mpi-film-checkpoint-',nranks,'.chk'
    call mpi_spray_checkpoint(MPI_COMM_WORLD,.true.,trim(path),'film-context',species,gas,temp,drops,pm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,time,step,event,mark,ref,inj,ok, &
      wall_films=patches,film_mass=cm,wall_totals=wall)
    call require(ok.and.size(patches)==1.and.patches(1)%temperature==330+rank,'film checkpoint write preserves ownership')
    if (rank==0) then
      call mpi_spray_checkpoint(MPI_COMM_WORLD,.false.,trim(path),'film-context',species,gas,temp,drops,pm, &
        [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,time,step,event,mark,ref,inj,ok, &
        wall_films=patches,film_mass=cm,wall_totals=wall)
    else
      call mpi_spray_checkpoint(MPI_COMM_WORLD,.false.,trim(path),'film-context',species,gas,temp,drops,pm, &
        [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,time,step,event,mark,ref,inj,ok)
    end if
    call require(.not.ok.and.wall(5)==1.25_dp.and.size(patches)==1,'film mode mismatch rejects collectively')
    deallocate(patches,cm);allocate(patches(0),cm(2,0));wall=0
    gas=reference*0;temp=reference_temperature*0
    call mpi_spray_checkpoint(MPI_COMM_SELF,.false.,trim(path),'film-context',species,gas,temp,drops,pm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,time,step,event,mark,ref,inj,ok, &
      wall_films=patches,film_mass=cm,wall_totals=wall)
    call require(ok.and.size(patches)==nranks.and.wall(5)==1.25_dp,'film changed-rank restart and wall budget')
    call require(all(gas==reference).and.all(temp==reference_temperature),'film checkpoint exact gas')
    do q=1,nranks
      call require(patches(q)%temperature==329+q.and.patches(q)%material_id==q,'film restored source order/material')
      call require(all(cm(:,q)==patches(q)%mass*[.25_dp,.75_dp]),'film restored component inventory')
    end do
    gas=field*0;temp=temperature*0
    call mpi_spray_checkpoint(MPI_COMM_WORLD,.false.,trim(path),'film-context',species,gas,temp,drops,pm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,time,step,event,mark,ref,inj,ok, &
      wall_films=patches,film_mass=cm,wall_totals=wall)
    call require(ok.and.size(patches)==1,'film restart redistributes from root')
    call require(patches(1)%temperature==330+rank,'film restart recovers original owner')
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if (rank==0) then
      open(newunit=unit,file=trim(path),status='old');close(unit,status='delete')
    end if
  end subroutine

  subroutine check_film_coupling()
    type(spray_wall_film), allocatable :: patches(:),original(:)
    type(spray_parcel), allocatable :: empty(:)
    type(spray_options) :: settings
    type(spray_wall_options) :: walls
    type(elementary_reaction), allocatable :: faulty(:),film_reactions(:)
    real(dp), allocatable :: cm(:,:),original_cm(:,:),empty_cm(:,:),gas(:,:,:,:),temp(:,:,:)
    real(dp) :: start(15),finish(15),global_start(15),global_finish(15),boundary(15),film_total
    integer :: count_steps,np
    character(len=16) :: film_faces(6)
    np=merge(1,0,rank==0)
    allocate(patches(np),cm(2,np),empty(0),empty_cm(2,0))
    if (np==1) then
      patches(1)=spray_wall_film([0.0_dp,.5_dp,.5_dp],[1.0_dp,0.0_dp,0.0_dp], &
        1e-9_dp,330.0_dp,.001_dp,1e-6_dp,1)
      cm(:,1)=patches(1)%mass*[.25_dp,.75_dp]
    end if
    original=patches;original_cm=cm;gas=field;temp=temperature
    film_faces='periodic';film_faces(1:2)='slip_wall'
    walls%enabled(1)=.true.;walls%evolve=.true.;walls%temperature=400
    walls%liquid_conductivity=.6_dp;walls%minimum_height=1e-10_dp
    settings%maximum_substeps=1
    call advance_mpi_film_coupling(MPI_COMM_WORLD,species,transport,liquids,settings,walls,gas,temp,patches,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],1e-3_dp,film_faces, &
      count_steps,boundary,ok)
    call require(.not.ok.and.count_steps==0.and.all(boundary==0),'film substep failure clears ledger')
    call require(all(gas==field).and.all(temp==temperature).and.all(cm==original_cm), &
      'film substep failure rolls back gas/composition')
    call require(all(patches%temperature==original%temperature).and.all(patches%mass==original%mass), &
      'film substep failure rolls back patch state')
    settings%maximum_substeps=100000
    call mixture_system_totals(gas,[1.0_dp,1.0_dp,1.0_dp],empty,liquids,empty_cm,start,ok,patches,cm)
    call require(ok,'film initial totals')
    call MPI_Allreduce(start,global_start,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call advance_mpi_film_coupling(MPI_COMM_WORLD,species,transport,liquids,settings,walls,gas,temp,patches,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],1e-8_dp,film_faces, &
      count_steps,boundary,ok)
    call require(ok.and.count_steps>0.and.boundary(5)>0,'MPI film heats with wall energy input')
    if (rank==0) call require(patches(1)%temperature>330,'MPI film temperature increases')
    film_total=sum(patches%mass)
    call MPI_Allreduce(MPI_IN_PLACE,film_total,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(film_total<1e-9_dp,'MPI film evaporation reduces liquid mass')
    call mixture_system_totals(gas,[1.0_dp,1.0_dp,1.0_dp],empty,liquids,empty_cm,finish,ok,patches,cm)
    call require(ok,'film final totals')
    call MPI_Allreduce(finish,global_finish,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(abs(global_finish-global_start-boundary)<2e-12_dp*max(1e-12_dp,abs(global_start))), &
      'MPI gas film wall budgets close')
    call load_h2o2_full_mechanism(film_reactions,ok);call require(ok,'film split mechanism')
    patches=original;cm=original_cm;gas=field;temp=temperature;faulty=film_reactions
    if (rank==nranks-1) faulty(1)%forward_rate%pre_exponential=-1
    call advance_mpi_mixture_reaction(MPI_COMM_WORLD,species,faulty,transport,liquids,settings,gas,temp,empty,empty_cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],1e-8_dp,1e-7_dp,1e-12_dp, &
      'implicit',count_steps,ok,boundary_faces=film_faces,boundary_change=boundary, &
      wall_films=patches,film_mass=cm,film_settings=walls)
    call require(.not.ok.and.count_steps==0.and.all(boundary==0),'late chemistry failure clears film split outputs')
    call require(all(gas==field).and.all(temp==temperature).and.all(cm==original_cm).and. &
      all(patches%mass==original%mass).and.all(patches%temperature==original%temperature), &
      'late chemistry failure rolls back earlier film heating and evaporation')
    call advance_mpi_mixture_reaction(MPI_COMM_WORLD,species,film_reactions,transport,liquids,settings,gas,temp,empty,empty_cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],1e-8_dp,1e-7_dp,1e-12_dp, &
      'implicit',count_steps,ok,boundary_faces=film_faces,boundary_change=boundary, &
      wall_films=patches,film_mass=cm,film_settings=walls)
    call require(ok.and.count_steps>=4.and.boundary(5)>0,'symmetric film particle chemistry split')
    if (rank==0) call require(patches(1)%temperature>330.and.patches(1)%mass<1e-9_dp,'split film heats and evaporates')
    call mixture_system_totals(gas,[1.0_dp,1.0_dp,1.0_dp],empty,liquids,empty_cm,finish,ok,patches,cm)
    call require(ok,'film split final totals')
    call MPI_Allreduce(finish,global_finish,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(abs(global_finish(:5)-global_start(:5)-boundary(:5))< &
      2e-12_dp*max(1e-12_dp,abs(global_start(:5)))),'film split mass momentum and wall energy budgets close')
    call require(all(abs(element_totals(global_finish-global_start-boundary))< &
      2e-10_dp*max(1e-12_dp,abs(element_totals(global_start)))),'film split element inventory closes')
  end subroutine

  subroutine check_film_migration()
    type(spray_wall_film), allocatable :: patches(:),original(:)
    real(dp), allocatable :: cm(:,:),original_cm(:,:)
    real(dp) :: total(2),initial(2)
    integer :: q,source_rank
    allocate(patches(2),cm(2,2))
    patches(1)=spray_wall_film([real(nranks,dp),.25_dp,.25_dp],[-1.0_dp,0.0_dp,0.0_dp], &
      1e-9_dp,330.0_dp+rank,.001_dp,1e-6_dp,rank+1)
    patches(2)=spray_wall_film([real(rank,dp)+.25_dp,2.0_dp,.25_dp],[0.0_dp,-1.0_dp,0.0_dp], &
      2e-9_dp,400.0_dp+rank,.002_dp,2e-6_dp,16)
    cm(:,1)=patches(1)%mass*[.25_dp,.75_dp];cm(:,2)=patches(2)%mass*[.5_dp,.5_dp]
    original=patches;original_cm=cm
    total=sum(cm,dim=2)
    call MPI_Allreduce(total,initial,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call migrate_mpi_wall_films(MPI_COMM_WORLD,patches,cm,[0.0_dp,0.0_dp,0.0_dp], &
      [1.0_dp,1.0_dp,1.0_dp],[nranks,2,2],[nranks,1,1],1,ok)
    call require(.not.ok.and.all(cm==original_cm),'film capacity rollback preserves components')
    call require(all(patches%temperature==original%temperature),'film capacity rollback preserves patches')
    call migrate_mpi_wall_films(MPI_COMM_WORLD,patches,cm,[0.0_dp,0.0_dp,0.0_dp], &
      [1.0_dp,1.0_dp,1.0_dp],[nranks,2,2],[nranks,1,1],100,ok)
    call require(ok.and.size(patches)==merge(nranks+1,1,rank==nranks-1),'film closed-upper-boundary owner')
    do q=1,size(patches)
      if (patches(q)%temperature<400) then
        source_rank=nint(patches(q)%temperature)-330
        call require(rank==nranks-1.and.patches(q)%material_id==source_rank+1,'film material routing')
        call require(all(patches(q)%position==[real(nranks,dp),.25_dp,.25_dp]),'film exact wall position')
        call require(all(patches(q)%normal==[-1.0_dp,0.0_dp,0.0_dp]),'film wall normal')
        call require(patches(q)%diameter==.001_dp.and.patches(q)%height==1e-6_dp,'film patch geometry')
        call require(all(cm(:,q)==original_cm(:,1)),'film first composition exact')
      else
        call require(patches(q)%temperature==400+rank.and.patches(q)%material_id==16,'film local wall patch')
        call require(all(cm(:,q)==original_cm(:,2)),'film second composition exact')
      end if
      patches(q)%position(1)=real(nranks,dp)
    end do
    call migrate_mpi_wall_films(MPI_COMM_WORLD,patches,cm,[0.0_dp,0.0_dp,0.0_dp], &
      [1.0_dp,1.0_dp,1.0_dp],[nranks,2,2],[nranks,1,1],100,ok)
    call require(ok.and.size(patches)==merge(2*nranks,0,rank==nranks-1),'film migration permits empty ranks')
    total=sum(cm,dim=2)
    call MPI_Allreduce(MPI_IN_PLACE,total,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(abs(total-initial)<2e-12_dp*initial),'film component inventory preserved')
  end subroutine

  subroutine check_kh_release()
    type(spray_parcel), allocatable :: drops(:),original(:)
    real(dp), allocatable :: cm(:,:),original_cm(:,:),energy(:),density(:),velocity(:,:),reynolds(:)
    type(spray_options) :: settings
    real(dp) :: initial(15),final(15),y(2),global_initial(15),global_final(15)
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:)
    integer :: coupled_steps,global_parcels
    integer :: np,q
    integer(int64) :: mark
    np=merge(1,0,rank<nranks-1)
    allocate(drops(np),cm(2,np),density(np),velocity(3,np),reynolds(np),energy(1));energy=123
    y=[.25_dp,.75_dp]
    if (np>0) then
      drops(1)=spray_parcel(int(rank+1,int64),[.25_dp,.25_dp,.25_dp], &
        [0.0_dp,0.0_dp,0.0_dp],1e-9_dp,330.0_dp,100.0_dp)
      cm(:,1)=drops(1)%mass*y
    end if
    density=10;velocity=0;velocity(1,:)=300;reynolds=2000
    settings%breakup_model='khrt';settings%khrt_children=.true.
    settings%surface_tension=.072_dp;settings%liquid_viscosity=.001_dp
    settings%khrt_reference_diameter_cubed=1e-25_dp
    settings%maximum_parcels=1;mark=int(nranks,int64);original=drops;original_cm=cm
    call advance_mpi_khrt_shedding(MPI_COMM_WORLD,liquids,settings,drops,cm,density,velocity,reynolds, &
      1e-6_dp,mark,energy,ok)
    call require(.not.ok.and.mark==nranks.and.all(energy==123),'KH capacity rollback includes IDs/energy')
    call require(all(drops%mass==original%mass).and.all(cm==original_cm),'KH capacity rollback includes composition')
    settings%maximum_parcels=100
    call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],drops,liquids,cm,initial,ok)
    call require(ok,'KH initial inventory')
    call advance_mpi_khrt_shedding(MPI_COMM_WORLD,liquids,settings,drops,cm,density,velocity,reynolds, &
      1e-6_dp,mark,energy,ok)
    call require(ok.and.size(drops)==5*np,'KH four-child release with empty rank')
    call require(mark==nranks+4*(nranks-1),'KH collective ID watermark')
    call require(size(energy)==np,'KH parent energy shape')
    if (np>0) then
      call require(energy(1)<0,'KH transverse kinetic energy debits gas')
      do q=2,5
        call require(drops(q)%id==nranks+4*rank+q-1,'KH rank-prefix child IDs')
        call require(drops(q)%parent_id==rank+1,'KH child lineage')
        call require(all(abs(cm(:,q)/drops(q)%mass-y)<1e-14_dp),'KH inherited composition')
      end do
    end if
    do q=1,2
      call require(abs(sum(cm(q,:)*drops%multiplicity)-sum(original_cm(q,:)*original%multiplicity)) &
        <2e-12_dp*max(1e-30_dp,sum(original_cm(q,:)*original%multiplicity)),'KH component inventory conservation')
    end do
    call mixture_system_totals(field,[1.0_dp,1.0_dp,1.0_dp],drops,liquids,cm,final,ok)
    call require(ok,'KH final inventory');final(5)=final(5)+sum(energy)
    call require(all(abs(final-initial)<2e-12_dp*max(1e-12_dp,abs(initial))),'KH mass momentum energy species conservation')
    drops=original;cm=original_cm;mark=int(nranks,int64)
    do q=1,size(drops)
      drops(q)%velocity=[1000.0_dp,0.0_dp,0.0_dp]
    end do
    original=drops;gas=field;temp=temperature
    settings%drag=.false.;settings%heat=.false.;settings%evaporation=.false.
    settings%maximum_substeps=1
    call advance_mpi_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,settings,gas,temp,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],.001_dp,coupled_steps,ok,mark)
    call require(.not.ok.and.coupled_steps==0.and.mark==nranks,'KH whole-interval rollback restores IDs')
    call require(all(gas==field).and.all(temp==temperature).and.all(cm==original_cm), &
      'KH whole-interval rollback restores gas and components')
    call require(all(drops%mass==original%mass).and.all(drops%khrt_shed_volume==original%khrt_shed_volume), &
      'KH whole-interval rollback restores parent state')
    settings%maximum_substeps=100000
    call mixture_system_totals(gas,[1.0_dp,1.0_dp,1.0_dp],drops,liquids,cm,initial,ok)
    call require(ok,'coupled KH initial inventory')
    call MPI_Allreduce(initial,global_initial,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call advance_mpi_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,settings,gas,temp,drops,cm, &
      [0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],1e-6_dp,coupled_steps,ok,mark)
    call require(ok.and.coupled_steps>0.and.mark>nranks,'coupled KH accepts children and IDs')
    q=size(drops)
    call MPI_Allreduce(q,global_parcels,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(global_parcels>=4*(nranks-1),'coupled KH children survive migration')
    call mixture_system_totals(gas,[1.0_dp,1.0_dp,1.0_dp],drops,liquids,cm,final,ok)
    call require(ok,'coupled KH final inventory')
    call MPI_Allreduce(final,global_final,15,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(abs(global_final-global_initial)<2e-12_dp*max(1e-12_dp,abs(global_initial))), &
      'coupled KH conserves gas and particle inventories')

  end subroutine

  function element_totals(inventory) result(elements)
    real(dp), intent(in) :: inventory(15)
    real(dp) :: elements(4),atoms(4)
    integer :: q
    elements=0
    do q=1,size(species)
      atoms=0
      select case(trim(species(q)%name))
      case('H2');atoms(1)=2
      case('H');atoms(1)=1
      case('O');atoms(2)=1
      case('O2');atoms(2)=2
      case('OH');atoms(:2)=[1,1]
      case('H2O');atoms(:2)=[2,1]
      case('HO2');atoms(:2)=[1,2]
      case('H2O2');atoms(:2)=[2,2]
      case('AR');atoms(3)=1
      case('N2');atoms(4)=2
      case default;call require(.false.,'unknown element inventory species')
      end select
      elements=elements+atoms*inventory(5+q)/species(q)%molecular_weight
    end do
  end function
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (condition) return
    print *,rank,': ',message
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine
end program
