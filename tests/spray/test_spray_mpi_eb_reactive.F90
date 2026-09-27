program test_spray_mpi_eb_reactive
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: build_eb_tree_extended_geometry_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_subcycle_3d
  use spray_mpi_eb_reactive_mod, only: advance_mpi_eb_spray_reactive,advance_mpi_eb_spray_sdc, &
    advance_mpi_eb_spray_sdc_subcycle
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_eb_tree_halo_3d_mod, only: fill_sdc_mpi_eb_tree_source_halo_3d
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options
  implicit none
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp), allocatable :: boundary(:),chemical(:)
  real(dp), allocatable :: forcing(:),exchange(:)
  integer(int64) :: last_id=1_int64,created=0_int64
  real(dp) :: accel(3)=0
  real(dp) :: theta
  character(len=32) :: solver,test_mode
  character(len=16), allocatable :: labels(:)
  real(dp), allocatable :: inlet(:,:),inlet_t(:)
  logical, allocatable :: open_boundary
  integer :: post_calls=0
  logical :: fail_post=.false.
  integer :: predictor_calls=0
  logical :: fail_predictor=.false.
  logical :: recursive_steps
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  type(spray_parcel), allocatable :: drops(:),saved_drops(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:),before(:),after(:),ledger(:)
  type(spray_wall_film), allocatable :: films(:),saved_films(:)
  type(spray_wall_options), allocatable :: film_settings
  real(dp), allocatable :: film_mass(:,:),saved_film_mass(:,:)
  integer :: film_steps
  integer :: steps,p
  real(dp) :: rho,cp,e,amount
  type(amr_patch_3d) :: patch(1)
  type(amr_patch_tree_3d) :: tree
  type(eb_geometry_3d) :: finest,g(2),extended_finest
  type(eb_geometry_3d), allocatable :: extended(:)
  real(dp) :: extended_phi(0:32,0:32,0:32)
  type(amr_tree_state_3d) :: field,temp,saved,saved_t
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: positions(:,:),feedback(:,:),cell(:),primitive(:),y(:),total(:),parent(:),fine(:)
  real(dp), allocatable :: wall_temperature(:)
  real(dp) :: phi(0:8,0:8,0:8),zero(3)=0,spacing(3)=.001_dp,t,sound,r,volume
  integer :: rank,nranks,ierr,owners(2),dims(3),node,i,j,k,nv,nitrogen
  logical :: ok
  call get_command_argument(1,test_mode)
  recursive_steps=index(test_mode,'subcycle')>0
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'rank count')
  owners=[0,nranks-1]
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patch(1),ok)
  call require(ok,'patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1],patch,tree,ok)
  call require(ok,'topology')
  do j=0,8
    phi(:,j,:)=real(j,dp)*.0005_dp-.00125_dp
  end do
  call build_level_set_eb_geometry_3d(phi,zero,[.004_dp,.004_dp,.004_dp],finest,ok)
  call require(ok,'fine geometry')
  call coarsen_eb_geometry_3d(finest,2,g(1),ok)
  call require(ok,'root geometry')
  call extract_eb_geometry_3d(finest,[1,1,1],[4,8,8],g(2),ok)
  call require(ok,'child geometry')
  if(index(test_mode,'ctu')==1.or.index(test_mode,'sdc')==1) then
    do j=0,32
      extended_phi(:,j,:)=real(max(0,min(8,j-12)),dp)*.0005_dp-.00125_dp
    end do
    call build_level_set_eb_geometry_3d(extended_phi,[-.006_dp,-.006_dp,-.006_dp], &
      [.010_dp,.010_dp,.010_dp],extended_finest,ok)
    call require(ok,'extended nodal geometry')
    call build_eb_tree_extended_geometry_3d(tree,g,extended_finest,extended,ok)
    call require(ok,'CTU patch geometry')
  end if
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermo');nv=reactive_nvar(size(species))
  allocate(cell(nv),primitive(nv),y(size(species)),total(nv),parent(nv),fine(nv))
  nitrogen=0
  do i=1,size(species)
    if(species(i)%name=='N2') nitrogen=i
  end do
  call require(nitrogen>0,'species')
  y=0;y(nitrogen)=1;r=mixture_specific_gas_constant(species,y,ok)
  call require(ok,'gas constant')
  primitive=0;primitive(1)=1;primitive(2:4)=[1._dp,2._dp,3._dp];primitive(5)=r*1000;primitive(6:)=y
  call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
  call require(ok,'gas state')
  allocate(field%fields(2),temp%fields(2),positions(3,merge(1,0,rank==0)),feedback(nv,merge(1,0,rank==0)))
  do node=1,2
    if(owners(node)/=rank) cycle
    dims=tree%nodes(node)%shape
    allocate(field%fields(node)%values(nv,dims(1),dims(2),dims(3)),temp%fields(node)%values(1,dims(1),dims(2),dims(3)))
    do k=1,dims(3)
    do j=1,dims(2)
    do i=1,dims(1)
      field%fields(node)%values(:,i,j,k)=cell;temp%fields(node)%values(1,i,j,k)=t
      if(g(node)%volume_fraction(i,j,k)==0) then
        field%fields(node)%values(:,i,j,k)=0;temp%fields(node)%values(1,i,j,k)=42
      end if
    end do
    end do
    end do
  end do
  call load_h2o2_full_transport(transport,ok)
  call require(ok,'transport')
  call load_h2o2_full_mechanism(reactions,ok)
  call require(ok,'mechanism')
  allocate(boundary(nv),chemical(nv),forcing(nv),exchange(nv))
  forcing=46;exchange=47
  if(index(test_mode,'sdc')==1) then
    if(test_mode=='sdc'.or.test_mode=='sdc-force'.or.test_mode=='sdc-recursive'.or. &
      test_mode=='sdc-recursive-crossing'.or.test_mode=='sdc-recursive-entering'.or. &
      test_mode=='sdc-recursive-breakup'.or.test_mode=='sdc-recursive-khrt'.or. &
      test_mode=='sdc-recursive-children'.or.test_mode=='sdc-recursive-film') then
      allocate(labels(6));labels='outflow'
    end if
  end if
  call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquids(1),ok)
  call require(ok,'water')
  call initialize_spray_liquid(species,'H2O2',800._dp,2000._dp,350._dp,101325._dp,1e6_dp,liquids(2),ok)
  call require(ok,'peroxide')
  allocate(drops(merge(1,0,rank==nranks-1)),masses(2,merge(1,0,rank==nranks-1)),before(nv),after(nv),ledger(nv))
  if(rank==nranks-1) then
    drops(1)=spray_parcel(1_int64,[.00199_dp,.001250001_dp,.0015_dp], &
      [20._dp,-20._dp,5._dp],1e-12_dp,350._dp,100._dp)
    masses(:,1)=drops(1)%mass*[.25_dp,.75_dp]
  end if
  if(index(test_mode,'inflow')>0) then
    allocate(labels(6),inlet(nv,6),inlet_t(6),open_boundary)
    labels='outflow';labels(1)='inflow';open_boundary=.true.
    do i=1,6
      inlet(:,i)=cell
    end do
    inlet_t=t
  end if
  if(index(test_mode,'periodic')>0) then
    allocate(labels(6),open_boundary)
    labels='periodic';labels(3:4)='outflow';open_boundary=.true.
    if(rank==nranks-1) drops(1)%position(1)=.00399999_dp
  end if
  if(index(test_mode,'wall')>0.or.index(test_mode,'isothermal')>0) then
    allocate(labels(6),open_boundary)
    labels='outflow';labels(1:2)='no_slip_wall';open_boundary=.true.
    if(rank==nranks-1) then
      drops(1)%position(1)=1.e-8_dp;drops(1)%velocity(1)=-20._dp
    end if
  end if
  if(index(test_mode,'isothermal')>0) then
    labels(1:2)='isothermal_wall';allocate(wall_temperature(6));wall_temperature=2000
  end if
  if(test_mode=='sdc-force') accel=[1000._dp,-2000._dp,3000._dp]
  saved=field;saved_t=temp;saved_drops=drops;saved_mass=masses
  if(index(test_mode,'sdc-recursive')==1) then
    call check_recursive_spray()
    if(test_mode=='sdc-recursive') call check_post_level_callback()
    call MPI_Finalize(ierr)
    stop
  end if
  ledger=42;boundary=43;chemical=44;theta=45;solver='invalid'
  call advance()
  call require(.not.ok.and.steps==0.and.all(ledger==42),'whole interval failure and ledger rollback')
  call require(all(masses==saved_mass).and.size(drops)==size(saved_drops),'particle population rollback')
  if(rank==nranks-1) then
    call require(all(drops(1)%position==saved_drops(1)%position),'particle position rollback')
    call require(all(drops(1)%velocity==saved_drops(1)%velocity),'particle velocity rollback')
    call require(drops(1)%temperature==saved_drops(1)%temperature,'particle temperature rollback')
    call require(drops(1)%mass==saved_drops(1)%mass,'particle mass rollback')
  end if
  do node=1,2
    if(owners(node)/=rank) cycle
    call require(all(field%fields(node)%values==saved%fields(node)%values),'gas rollback')
    call require(all(temp%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
  end do
  call require(all(boundary==43).and.all(chemical==44).and.theta==45,'gas diagnostics rollback')
  if(index(test_mode,'sdc')==1) call require(last_id==1.and.all(forcing==46).and.all(exchange==47), &
    'SDC ID and exchange rollback')
  solver='hllc'
  if(index(test_mode,'ctu')==1) solver='pelec'
  if(allocated(inlet)) then
    if(rank==0) inlet(1,1)=-1
    call advance()
    call require(.not.ok.and.steps==0.and.all(boundary==43).and.all(ledger==42),'inlet failure rollback')
    call require(all(masses==saved_mass).and.size(drops)==size(saved_drops),'inlet population rollback')
    do node=1,2
      if(owners(node)/=rank) cycle
      call require(all(field%fields(node)%values==saved%fields(node)%values),'inlet gas rollback')
    end do
    inlet(1,1)=cell(1)
  end if
  if(index(test_mode,'sdc')==1) then
    options%breakup_model='etab';options%surface_tension=1e3_dp;options%liquid_viscosity=.001_dp
    options%maximum_multiplicity=50;options%maximum_parcels=1
    call advance()
    call require(.not.ok.and.steps==0.and.created==0.and.last_id==1,'late split capacity rollback')
    call require(all(masses==saved_mass).and.size(drops)==size(saved_drops),'split population rollback')
    call require(all(boundary==43).and.all(ledger==42).and.all(chemical==44).and.theta==45, &
      'split budget rollback')
    call require(all(forcing==46).and.all(exchange==47),'split source rollback')
    do node=1,2
      if(owners(node)/=rank) cycle
      call require(all(field%fields(node)%values==saved%fields(node)%values),'split gas rollback')
      call require(all(temp%fields(node)%values==saved_t%fields(node)%values),'split temperature rollback')
    end do
    options%maximum_parcels=3
    if(rank==0) options%maximum_multiplicity=0
    call advance()
    call require(.not.ok.and.steps==0.and.created==0.and.last_id==1,'split controls agree across ranks')
    call require(all(masses==saved_mass).and.all(forcing==46),'split disagreement rollback')
    options%maximum_multiplicity=50
  end if
  call inventory(before)
  call advance()
  call require(ok.and.steps>1,'distributed EB interval')
  call inventory(after)
  if(index(test_mode,'sdc')==1) after=after-forcing
  call require(all(abs(after-before-ledger-boundary-chemical)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'composite gas parcel wall conservation')
  call require(any(chemical(6:)/=0),'chemical species exchange')
  if(index(test_mode,'sdc')==1) then
    call require(size(drops)==merge(2,0,rank==nranks-1),'split final fine owner')
    call require(created==1.and.last_id==2,'accepted split ID')
  else
    call require(size(drops)==merge(1,0,rank==nranks-1),'final fine owner')
  end if
  if(rank==nranks-1) then
    call require(drops(1)%position(1)<.002_dp.and.drops(1)%position(2)>.00125_dp,'fine cell and fluid position')
    if(index(test_mode,'wall')>0.or.index(test_mode,'isothermal')>0) call require(drops(1)%velocity(1)>0,'outer wall reflection')
    call require(drops(1)%temperature>350.and.drops(1)%mass<1e-12_dp,'heating and evaporation')
    call require(drops(1)%velocity(2)>0.and.drops(1)%velocity(1)<20,'reflection and drag')
  end if
  if(index(test_mode,'sdc')==1) then
    if(test_mode=='sdc-force') then
      call require(any(forcing(2:5)/=0).and.forcing(1)==0.and.all(forcing(6:)==0),'external impulse and work')
    else
      call require(all(forcing==0),'zero external forcing')
    end if
    call require(any(exchange/=0),'internal particle exchange')
    if(rank==0) print *, 'SDC combined residual: ', &
      maxval(abs(after-before-ledger-boundary-chemical)/max(abs(before),abs(after),1e-30_dp))
  end if
  if(test_mode=='sdc-force') then
    call check_ballistic_force()
    call check_gas_acceleration()
  end if
  if(rank==0) print *, 'PASS: MPI EB coupled reactive step and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine check_post_level_callback()
    type(amr_tree_state_3d) :: old,old_t
    integer :: attempt
    old=field;old_t=temp
    do attempt=1,2
      post_calls=0;fail_post=attempt==1
      boundary=17;ledger=19;chemical=23;forcing=29;theta=31
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,g,field,temp, &
        1e-8_dp,2,.false.,[.false.,.false.,.false.,.false.],zero,2e-7_dp,1e-12_dp,extended, &
        labels,boundary,ledger,chemical,forcing,theta,MPI_COMM_WORLD,ok,post_level=check_post_context)
      call require(post_calls==3,'post-level callback follows both child steps and root synchronization')
      if(attempt==1) then
        call require(.not.ok,'one-rank final post-level failure rejects entire gas step')
        call require(all(boundary==17).and.all(ledger==19).and.all(chemical==23).and.all(forcing==29).and. &
          theta==31,'post-level failure retains gas ledgers')
        do node=1,2
          if(owners(node)/=rank) cycle
          call require(all(field%fields(node)%values==old%fields(node)%values),'post-level gas rollback')
          call require(all(temp%fields(node)%values==old_t%fields(node)%values),'post-level temperature rollback')
        end do
      else
        call require(ok,'post-level accepted retry')
      end if
    end do
  end subroutine check_post_level_callback
  subroutine check_post_context(depth,time,interval,iteration,ncycle,comm,accepted)
    integer, intent(in) :: depth,iteration,ncycle
    real(dp), intent(in) :: time,interval
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: accepted
    integer :: callback_rank,status
    call MPI_Comm_rank(comm,callback_rank,status)
    call require(status==MPI_SUCCESS,'post-level communicator')
    post_calls=post_calls+1
    if(post_calls<=2) then
      call require(depth==1.and.iteration==post_calls.and.ncycle==2,'child post-level iteration context')
      call require(abs(time-real(post_calls-1,dp)*5e-9_dp)<1e-22_dp.and. &
        abs(interval-5e-9_dp)<1e-22_dp,'child post-level clock')
    else
      call require(post_calls==3.and.depth==0.and.iteration==1.and.ncycle==1,'root post-level context')
      call require(time==0.and.interval==1e-8_dp,'root post-level clock')
    end if
    accepted=.not.(fail_post.and.depth==0.and.callback_rank==0)
  end subroutine check_post_context
  subroutine check_recursive_spray()
    integer(int64) :: counts(2)
    integer :: attempt,final_owner,attempts,final_count
    logical :: breakup_test,children_test
    real(dp) :: coefficient,gas_before(nv),gas_after(nv),particle_change(nv)
    if(rank==nranks-1) drops(1)%position=[.0012_dp,.0025_dp,.0015_dp]
    final_owner=nranks-1
    if(test_mode=='sdc-recursive-crossing') then
      if(rank==nranks-1) drops(1)%position(1)=.00199999_dp
      final_owner=0
    end if
    if(test_mode=='sdc-recursive-entering') then
      if(rank==nranks-1) then
        drops(1)%position(1)=.00200001_dp;drops(1)%velocity(1)=-20._dp
      end if
    end if
    children_test=test_mode=='sdc-recursive-children'
    breakup_test=children_test.or.test_mode=='sdc-recursive-breakup'.or.test_mode=='sdc-recursive-khrt'
    attempts=2
    if(breakup_test) then
      attempts=3;options%breakup_model='etab';options%surface_tension=.072_dp
      options%liquid_viscosity=.001_dp;options%maximum_multiplicity=50
      if(test_mode=='sdc-recursive-khrt') options%breakup_model='khrt'
      if(rank==final_owner) drops(1)%deformation=1.1_dp
    end if
    if(children_test) then
      options%breakup_model='khrt';options%khrt_children=.true.;options%maximum_multiplicity=0
      options%khrt_reference_diameter_cubed=1e-25_dp
      options%drag=.false.;options%heat=.false.;options%evaporation=.false.
      if(rank==final_owner) drops(1)%velocity=[1000._dp,10._dp,5._dp]
    end if
    if(test_mode=='sdc-recursive-film') then
      attempts=3
      allocate(film_settings,films(merge(1,0,rank==0.or.rank==nranks-1)))
      allocate(film_mass(2,size(films)))
      film_settings%enabled=.true.;film_settings%evolve=.true.;film_settings%temperature=400
      film_settings%liquid_conductivity=.6_dp;film_settings%minimum_height=1e-15_dp
      if(size(films)>0) then
        films(1)=spray_wall_film([merge(.0025_dp,.0015_dp,rank==0),.00125_dp,.0015_dp], &
          [0._dp,1._dp,0._dp],1e-10_dp,350._dp,1e-4_dp,1e-5_dp,merge(1,2,rank==0))
        film_mass(:,1)=[.6_dp,.4_dp]*films(1)%mass
      end if
      saved_films=films;saved_film_mass=film_mass
    end if
    saved_drops=drops;saved_mass=masses
    call inventory(before)
    call inventory(gas_before,.false.)
    do attempt=1,attempts
      if(breakup_test) options%maximum_parcels=merge(1,10000,attempt==2)
      if(allocated(film_settings)) then
        film_settings%port_conductivity_models(2)=''
        film_settings%port_conductivity_coefficients(:,2)=0
        film_settings%port_conductivity_temperature_ranges(:,2)=0
        if(attempt==2) then
          film_settings%port_conductivity_models(2)='polynomial'
          film_settings%port_conductivity_coefficients(1,2)=.6_dp
          film_settings%port_conductivity_temperature_ranges(:,2)=[300._dp,340._dp]
        end if
      end if
      coefficient=merge(-1._dp,.1_dp,attempt==1)
      boundary=17;ledger=19;chemical=23;forcing=29;theta=31;exchange=37;counts=41
      if(allocated(film_settings)) then
      call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,tree,owners,g,extended,phi, &
        species,reactions,transport,liquids,options,field,temp,drops,masses,zero,spacing, &
        1e-8_dp,2,.true.,[.true.,.true.,.true.,.true.],zero,2e-7_dp,1e-12_dp,labels, &
        boundary,ledger,chemical,forcing,theta,steps,last_id,ok,difmag=coefficient, &
        gas_particle_change=exchange,level_steps=counts, &
        inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=wall_temperature,split_created=created, &
        films=films,film_mass=film_mass,film_settings=film_settings,film_steps=film_steps)
      else
      call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,tree,owners,g,extended,phi, &
        species,reactions,transport,liquids,options,field,temp,drops,masses,zero,spacing, &
        1e-8_dp,2,.true.,[.true.,.true.,.true.,.true.],zero,2e-7_dp,1e-12_dp,labels, &
        boundary,ledger,chemical,forcing,theta,steps,last_id,ok,difmag=coefficient, &
        gas_particle_change=exchange,level_steps=counts, &
        inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=wall_temperature,split_created=created)
      end if
      if(attempt<attempts) then
        call require(.not.ok.and.steps==0.and.last_id==1,'recursive early or late rejection retains ID')
        call require(all(boundary==17).and.all(ledger==19).and.all(chemical==23).and.all(forcing==29).and. &
          all(exchange==37).and.created==0.and.theta==31.and.all(counts==41),'recursive whole-step ledger rollback')
        call require(size(drops)==size(saved_drops).and.all(masses==saved_mass),'recursive population rollback')
        if(allocated(films)) then
          call require(film_steps==0.and.size(films)==size(saved_films).and.all(film_mass==saved_film_mass), &
            'recursive film inventory rollback')
          call require(all(films%mass==saved_films%mass).and.all(films%temperature==saved_films%temperature), &
            'recursive film thermal rollback')
        end if
        do node=1,2
          if(owners(node)/=rank) cycle
          call require(all(field%fields(node)%values==saved%fields(node)%values),'recursive gas rollback')
          call require(all(temp%fields(node)%values==saved_t%fields(node)%values),'recursive temperature rollback')
        end do
      else
        call require(ok.and.steps==6.and.all(counts==[1_int64,2_int64]),'coupled recursive gas particle intervals')
      end if
    end do
    if(allocated(films)) then
      call require(film_steps==6.and.size(films)==size(saved_films),'recursive physical film collection and stages')
      if(size(films)>0) call require(films(1)%mass<saved_films(1)%mass.and. &
        films(1)%temperature>saved_films(1)%temperature,'recursive films evaporate and heat on both levels')
    end if
    final_count=1
    if(children_test) then
      final_count=size(drops)
      call MPI_Allreduce(MPI_IN_PLACE,final_count,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      call require(last_id>1.and.final_count>1.and.created==0,'recursive KHRT shedding IDs, separate from splitting')
      call require(exchange(5)<0.and.all(exchange(:4)==0).and.all(exchange(6:)==0), &
        'recursive child transverse energy debited from gas only')
      if(rank==0) print *, 'Recursive KHRT live count / watermark / gas energy: ',final_count,last_id,exchange(5)
      if(rank==final_owner) then
        call require(maxval(drops%id)==last_id,'recursive child watermark')
        do p=2,size(drops)
          call require(drops(p)%id>drops(p-1)%id,'recursive child IDs unique and sorted')
        end do
      end if
    else if(breakup_test) then
      final_count=int(last_id)
      call require(created==last_id-1.and.created>=1,'recursive computational split IDs')
      if(rank==final_owner) then
        call require(all(drops%multiplicity<=50._dp),'recursive split multiplicity bound')
        if(test_mode=='sdc-recursive-breakup') call require(all(drops%mass<.9_dp*saved_drops(1)%mass), &
          'ETAB changes droplet mass beyond evaporation')
      end if
    end if
    call inventory(after)
    call inventory(gas_after,.false.)
    particle_change=(after-gas_after)-(before-gas_before)
    call require((children_test.or.last_id==int(final_count,int64)).and.size(drops)==merge(final_count,0,rank==final_owner), &
      'recursive physical collection')
    if(rank==final_owner.and.test_mode=='sdc-recursive-crossing') &
      call require(drops(1)%position(1)>.002_dp,'crossed fine boundary')
    if(rank==final_owner.and.test_mode=='sdc-recursive-entering') &
      call require(drops(1)%position(1)<.002_dp,'entered fine patch')
    if(rank==0.and.(test_mode=='sdc-recursive-crossing'.or.test_mode=='sdc-recursive-entering')) then
      print *, 'Crossing composite residual: ', &
        maxval(abs(after-before-ledger-boundary-chemical-forcing)/max(abs(before),abs(after),1e-30_dp))
      print *, 'Gas source closure: ',maxval(abs(gas_after-gas_before-ledger-boundary-chemical-forcing-exchange))
      print *, 'Particle plus deposited exchange: ',particle_change+exchange
    end if
    call require(all(abs(after-before-ledger-boundary-chemical-forcing)<= &
      3e-12_dp*max(abs(before),abs(after),1e-30_dp)),'recursive composite gas particle inventory')
    call require((children_test.or.last_id==int(final_count,int64)).and.size(drops)==merge(final_count,0,rank==final_owner), &
      'recursive physical collection')
    call require(any(exchange/=0),'recursive particle feedback reaches gas')
    if(rank==final_owner) then
      if(test_mode=='sdc-recursive-crossing') call require(drops(1)%position(1)>.002_dp,'crossed fine boundary')
      if(.not.children_test) call require(all(drops%temperature>350._dp).and.all(drops%mass<1e-12_dp), &
        'recursive heating and evaporation')
    end if
    if(rank==0) print *, 'Recursive coupled residual: ', &
      maxval(abs(after-before-ledger-boundary-chemical-forcing)/max(abs(before),abs(after),1e-30_dp))
  end subroutine check_recursive_spray
  subroutine check_halo_context(state,temperature,interval,stage,depth,time,rate,comm,accepted, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(amr_tree_state_3d), intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(in) :: interval,time,ancestor_fraction(:)
    integer, intent(in) :: stage,depth
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: accepted
    real(dp) :: expected_fraction
    integer :: a
    call require(size(ancestor_fraction)==1,'two-level ancestor fraction shape')
    expected_fraction=0
    if(depth==1) expected_fraction=(time+merge(interval,0._dp,stage>0))/1e-8_dp
    call require(abs(ancestor_fraction(1)-expected_fraction)<1e-14_dp,'stage-start/end ancestor fraction')
    do a=1,2
      call require(allocated(ancestor_start%fields(a)%values).eqv.(owners(a)==rank),'ancestor snapshot ownership')
      if(owners(a)/=rank) cycle
      call require(all(shape(ancestor_start%fields(a)%values)==shape(state%fields(a)%values)), &
        'ancestor conserved start shape')
      call require(all(shape(ancestor_end%fields(a)%values)==shape(state%fields(a)%values)), &
        'ancestor conserved end shape')
      call require(all(shape(temperature_start%fields(a)%values)==shape(temperature%fields(a)%values)).and. &
        all(shape(temperature_end%fields(a)%values)==shape(temperature%fields(a)%values)), &
        'ancestor temperature shapes')
      if(a==1.and.depth==1) call require(any(ancestor_end%fields(a)%values(2,:,:,:)/= &
        ancestor_start%fields(a)%values(2,:,:,:)),'ancestor endpoint carries completed coarse acceleration')
    end do
    rate=state
    do a=1,2
      if(owners(a)==rank) rate%fields(a)%values=0
    end do
    accepted=.true.
  end subroutine check_halo_context
  subroutine recursive_predictor(rate,interval,iteration,depth,time,halo,comm,accepted, &
      source_start,source_end,source_fraction)
    type(amr_tree_state_3d), intent(in) :: rate,source_start,source_end
    real(dp), intent(in) :: interval,time,source_fraction(:)
    integer, intent(in) :: iteration,depth
    type(amr_tree_state_3d), intent(inout) :: halo
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: accepted
    type(amr_tree_state_3d) :: metrics
    integer :: a,n(3)
    predictor_calls=predictor_calls+1
    call require(iteration==modulo(predictor_calls-1,2)+1,'recursive predictor iteration')
    call require(size(source_fraction)==1,'predictor ancestor fraction shape')
    if(depth==0) then
      call require(time==0.and.interval==1e-8_dp.and.all(source_fraction==0),'root predictor clock')
    else
      call require(depth==1.and.interval==5e-9_dp,'fine predictor interval')
      call require(abs(source_fraction(1)-(time+.5_dp*interval)/1e-8_dp)<1e-14_dp, &
        'predictor parent midpoint fraction')
    end if
    allocate(metrics%fields(2))
    do a=1,2
      call require(allocated(source_start%fields(a)%values).eqv.(owners(a)==rank),'predictor donor ownership')
      if(owners(a)/=rank) cycle
      n=tree%nodes(a)%shape
      allocate(metrics%fields(a)%values(1,n(1)+6,n(2)+6,n(3)+6))
      metrics%fields(a)%values(1,:,:,:)=extended(a)%volume_fraction(4:n(1)+9,4:n(2)+9,4:n(3)+9)
    end do
    call fill_sdc_mpi_eb_tree_source_halo_3d(tree,owners,g,rate,nv,3,halo,comm,accepted, &
      boundary_faces=labels,target_volume_fraction=metrics,ancestor_start=source_start, &
      ancestor_end=source_end,ancestor_fraction=source_fraction)
    if(fail_predictor.and.predictor_calls==6.and.rank==nranks-1) accepted=.false.
  end subroutine recursive_predictor
  subroutine check_gas_acceleration()
    real(dp) :: acceleration(3),force_density(3),dv(3),velocity(3),expected(nv),interval,density
    integer :: density_case,advance_mode
    type(amr_tree_state_3d) :: old,old_t
    acceleration=[1000._dp,0._dp,3000._dp];force_density=[20._dp,0._dp,40._dp]
    velocity=[1._dp,0._dp,3._dp];interval=1e-8_dp
    labels='outflow'
    deallocate(drops,masses);allocate(drops(0),masses(2,0))
    do advance_mode=1,3
    do density_case=1,2
      density=2._dp*density_case;field=saved;temp=saved_t
      primitive=0;primitive(1)=density;primitive(2:4)=velocity;primitive(5)=density*r*1000;primitive(6:)=y
      call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
      call require(ok,'non-unit density gas fixture')
      do node=1,2
        if(owners(node)/=rank) cycle
        dims=tree%nodes(node)%shape
        do k=1,dims(3)
        do j=1,dims(2)
        do i=1,dims(1)
          if(g(node)%volume_fraction(i,j,k)<=0) cycle
          field%fields(node)%values(:,i,j,k)=cell
          temp%fields(node)%values(1,i,j,k)=t
        end do
        end do
        end do
      end do
      call inventory(before)
      dv=interval*(acceleration+force_density/density)
      expected=0;expected(2:4)=before(1)*dv
      expected(5)=before(1)*dot_product(velocity+.5_dp*dv,dv)
      predictor_calls=0
      if(advance_mode==1) then
      call advance_mpi_eb_spray_sdc(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
        liquids,options,field,temp,drops,masses,zero,spacing,interval,2,.false., &
        [.false.,.false.,.false.,.false.],force_density,2e-7_dp,1e-12_dp,labels,boundary,ledger,chemical, &
        forcing,theta,steps,last_id,ok,gas_particle_change=exchange,gas_acceleration=acceleration)
      else if(advance_mode==2) then
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,g,field,temp, &
          interval,2,.false.,[.false.,.false.,.false.,.false.],force_density,2e-7_dp,1e-12_dp,extended, &
          labels,boundary,ledger,chemical,forcing,theta,MPI_COMM_WORLD,ok,gas_acceleration=acceleration, &
          predictor_builder=recursive_predictor)
      else
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,g,field,temp, &
          interval,2,.false.,[.false.,.false.,.false.,.false.],force_density,2e-7_dp,1e-12_dp,extended, &
          labels,boundary,ledger,chemical,forcing,theta,MPI_COMM_WORLD,ok,gas_acceleration=acceleration, &
          halo_source=check_halo_context,auxiliary_change=exchange,predictor_builder=recursive_predictor)
      end if
      call require(ok,'non-unit density accelerated gas')
      if(advance_mode>1) call require(predictor_calls==6,'recursive predictor covers all intervals')
      call require(all(abs(forcing-expected)<=3e-13_dp*max(abs(expected),1e-30_dp)),'density-weighted impulse/work')
      call require(all(exchange==0),'acceleration stays outside particle exchange')
      do node=1,2
        if(owners(node)/=rank) cycle
        dims=tree%nodes(node)%shape
        do k=1,dims(3)
        do j=1,dims(2)
        do i=1,dims(1)
          if(g(node)%volume_fraction(i,j,k)<=0) cycle
          call require(maxval(abs(field%fields(node)%values(2:4,i,j,k)/density-velocity-dv))<1e-12_dp, &
            'force density and acceleration yield analytic velocity')
          call require(abs(temp%fields(node)%values(1,i,j,k)-1000)<1e-8_dp,'external work retains temperature')
        end do
        end do
        end do
      end do
      if(advance_mode==3.and.density_case==2) then
        old=field;old_t=temp;predictor_calls=0;fail_predictor=.true.
        boundary=17;ledger=19;chemical=23;forcing=29;theta=31;exchange=37
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,g,field,temp, &
          interval,2,.false.,[.false.,.false.,.false.,.false.],force_density,2e-7_dp,1e-12_dp,extended, &
          labels,boundary,ledger,chemical,forcing,theta,MPI_COMM_WORLD,ok,gas_acceleration=acceleration, &
          halo_source=check_halo_context,auxiliary_change=exchange,predictor_builder=recursive_predictor)
        call require(.not.ok.and.predictor_calls==6,'last descendant predictor failure')
        do node=1,2
          if(owners(node)/=rank) cycle
          call require(all(field%fields(node)%values==old%fields(node)%values),'recursive predictor gas rollback')
          call require(all(temp%fields(node)%values==old_t%fields(node)%values),'recursive predictor T rollback')
        end do
        call require(all(boundary==17).and.all(ledger==19).and.all(chemical==23).and.all(forcing==29).and. &
          theta==31.and.all(exchange==37),'recursive predictor inventory rollback')
        fail_predictor=.false.
      end if
    end do
    end do
  end subroutine check_gas_acceleration
  subroutine check_ballistic_force()
    real(dp) :: dv(3),expected(nv),initial_velocity(3),initial_position(3),interval
    integer :: iteration_count
    interval=1e-8_dp;initial_velocity=[20._dp,-20._dp,5._dp]
    initial_position=[.0015_dp,.002_dp,.0015_dp];dv=interval*accel
    options=spray_options();options%drag=.false.;options%heat=.false.;options%evaporation=.false.
    do iteration_count=1,3,2
      field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass;last_id=1
      if(rank==nranks-1) drops(1)%position=initial_position
      expected=0;expected(2:4)=1e-10_dp*dv
      expected(5)=1e-10_dp*dot_product(initial_velocity+.5_dp*dv,dv)
      call advance_mpi_eb_spray_sdc(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
        liquids,options,field,temp,drops,masses,zero,spacing,interval,iteration_count,.false., &
        [.false.,.false.,.false.,.false.],zero,2e-7_dp,1e-12_dp,labels,boundary,ledger,chemical, &
        forcing,theta,steps,last_id,ok,gas_particle_change=exchange,particle_acceleration=accel)
      call require(ok,'ballistic particle forcing')
      call require(all(abs(forcing-expected)<=1e-13_dp*max(abs(expected),1e-30_dp)), &
        'analytic particle external impulse and work')
      call require(all(exchange==0),'particle external force is not gas exchange')
      if(rank==nranks-1) then
        call require(maxval(abs(drops(1)%velocity-initial_velocity-dv))<1e-14_dp,'analytic velocity')
        call require(maxval(abs(drops(1)%position-initial_position- &
          interval*(initial_velocity+.5_dp*dv)))<1e-18_dp,'analytic position without repeated drift')
      end if
    end do
    if(rank==0) accel(1)=0
    forcing=46;exchange=47
    call advance()
    call require(.not.ok.and.steps==0.and.all(forcing==46).and.all(exchange==47),'particle force agreement rollback')
    accel=[1000._dp,-2000._dp,3000._dp]
  end subroutine check_ballistic_force
  subroutine advance()
    real(dp) :: coefficient
    if(index(test_mode,'sdc')==1) then
      coefficient=.1_dp
      if(solver=='invalid') coefficient=-1
      call advance_mpi_eb_spray_sdc(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
        liquids,options,field,temp,drops,masses,zero,spacing,1e-8_dp,2,.true., &
        [.true.,.true.,.true.,.true.],accel,2e-7_dp,1e-12_dp,labels,boundary,ledger,chemical, &
        forcing,theta,steps,last_id,ok,difmag=coefficient,gas_particle_change=exchange, &
        inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=wall_temperature,split_created=created, &
        particle_acceleration=accel)
      return
    end if
    call advance_mpi_eb_spray_reactive(MPI_COMM_WORLD,tree,owners,g,phi,species,reactions,transport, &
      liquids,options,field,temp,drops,masses,zero,spacing,1e-8_dp,solver, &
      .true.,.true.,.true.,.true.,.true.,.true.,2e-7_dp,1e-12_dp,steps,theta,boundary,ledger,chemical,ok,'implicit', &
      subcycling=recursive_steps,extended_geometry=extended,particle_outflow=open_boundary, &
      boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=wall_temperature)
  end subroutine
  subroutine inventory(result,with_particles)
    real(dp), intent(out) :: result(:)
    logical, optional, intent(in) :: with_particles
    logical :: include_particles
    include_particles=.true.
    if(present(with_particles)) include_particles=with_particles
    result=0
    do node=1,2
      if(owners(node)/=rank) cycle
      dims=tree%nodes(node)%shape
      do k=1,dims(3)
      do j=1,dims(2)
      do i=1,dims(1)
        if(node==1.and.i<=2) cycle
        volume=g(node)%volume_fraction(i,j,k)*g(node)%dx*g(node)%dy*g(node)%dz
        result=result+volume*field%fields(node)%values(:,i,j,k)
      end do
      end do
      end do
    end do
    if(include_particles) then
    do p=1,size(drops)
      call liquid_mixture_properties(liquids,masses(:,p)/drops(p)%mass,drops(p)%temperature,rho,cp,e,ok)
      call require(ok,'liquid properties')
      amount=drops(p)%mass*drops(p)%multiplicity
      result(1)=result(1)+amount;result(2:4)=result(2:4)+amount*drops(p)%velocity
      result(5)=result(5)+amount*(e+.5_dp*sum(drops(p)%velocity**2))
      do i=1,2
        result(5+liquids(i)%vapor_index)=result(5+liquids(i)%vapor_index)+masses(i,p)*drops(p)%multiplicity
      end do
    end do
    if(allocated(films)) then
      do p=1,size(films)
        call liquid_mixture_properties(liquids,film_mass(:,p)/films(p)%mass,films(p)%temperature,rho,cp,e,ok)
        call require(ok,'film liquid properties')
        result(1)=result(1)+films(p)%mass;result(5)=result(5)+films(p)%mass*e
        do i=1,2
          result(5+liquids(i)%vapor_index)=result(5+liquids(i)%vapor_index)+film_mass(i,p)
        end do
      end do
    end if
    end if
    call MPI_Allreduce(MPI_IN_PLACE,result,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',rank,': ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine
end program
