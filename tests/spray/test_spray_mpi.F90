program test_spray_mpi
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use spray_parcel_mod
  use spray_mpi_mod
  use spray_partition_mod, only: spray_cartesian_block
  use spray_coupling_3d_mod, only: particle_boundary_stencil,mixture_system_totals
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:),bad_reactions(:)
  type(spray_liquid) :: liquids(2)
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_options) :: options
  real(dp) :: local_totals(15),before(15),after(15)
  integer :: steps
  real(dp) :: fractions(2),ledger(15),direction(3)
  type(spray_parcel), allocatable :: parcels(:),saved(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:)
  integer :: rank,nranks,ierr,p,target,first(3),dims(3),cells(3),i,j,k,n,indices(3,8)
  real(dp), allocatable :: field(:,:,:,:),positions(:,:),samples(:,:),saved_samples(:,:)
  real(dp), allocatable :: contributions(:,:),reference(:,:,:,:),saved_field(:,:,:,:)
  real(dp), allocatable :: temperature(:,:,:),saved_temperature(:,:,:),reference_temperature(:,:,:)
  real(dp) :: primitive(15),gas_cell(15),gas_t,sound
  integer :: source
  real(dp) :: weights(8),expected_value
  character(len=16) :: faces(6)
  integer(int64) :: last_id,first_id,saved_last
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
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
  call MPI_Finalize(ierr)
contains
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
