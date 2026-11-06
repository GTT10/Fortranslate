program test_spray_mpi_eb_reactive
  use dual_energy_mod, only: dual_energy_policy,initialize_nasa7_internal_energy_patch
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
  use spray_mpi_amr_step_mod, only: convert_mpi_amr_tree_rotating_energy
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use elementary_kinetics_mod, only: elementary_reaction
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi_eb_tree_halo_3d_mod, only: fill_sdc_mpi_eb_tree_source_halo_3d,fill_mpi_eb_tree_halo_3d
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options
  use les_mpi_eb_sdc_mod, only: compute_mpi_regular_eb_les_sdc_3d,compute_mpi_regular_eb_les_spatial_sdc_3d
  use mpi_eb_tree_full_3d_mod, only: mpi_eb_sdc_diffusion_source
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
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
  logical :: recursive_steps,les_failure=.false.
  integer :: les_calls=0
  integer :: rejected_reactor_calls=0
  real(dp) :: les_largest=0
  character(len=16) :: les_model='smagorinsky'
  type(eb_geometry_3d) :: les_cut(2)
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
  real(dp), allocatable :: extended_phi(:,:,:)
  integer :: ctu_padding
  character(len=16) :: recursive_redistribution
  type(amr_tree_state_3d) :: field,temp,saved,saved_t
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: positions(:,:),feedback(:,:),cell(:),primitive(:),y(:),total(:),parent(:),fine(:)
  real(dp), allocatable :: wall_temperature(:)
  real(dp) :: phi(0:8,0:8,0:8),zero(3)=0,spacing(3)=.001_dp,t,sound,r,volume
  integer :: rank,nranks,ierr,owners(2),dims(3),node,i,j,k,nv,nitrogen
  logical :: ok
  call get_command_argument(1,test_mode)
  ctu_padding=6;recursive_redistribution='none'
  if(test_mode=='sdc-recursive-state-redist') then
    ctu_padding=10;recursive_redistribution='state_redist'
  endif
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
  if(index(test_mode,"sdc-les")==1.or.test_mode=="sdc-mol-sources") phi=1
  call build_level_set_eb_geometry_3d(phi,zero,[.004_dp,.004_dp,.004_dp],finest,ok)
  call require(ok,'fine geometry')
  call coarsen_eb_geometry_3d(finest,2,g(1),ok)
  call require(ok,'root geometry')
  call extract_eb_geometry_3d(finest,[1,1,1],[4,8,8],g(2),ok)
  call require(ok,'child geometry')
  if(index(test_mode,'sdc-les')==1) then
    block
      type(eb_geometry_3d) :: cut
      real(dp) :: cut_phi(0:8,0:8,0:8)
      do j=0,8
        cut_phi(:,j,:)=real(j,dp)*.0005_dp-.00125_dp
      end do
      call build_level_set_eb_geometry_3d(cut_phi,zero,[.004_dp,.004_dp,.004_dp],cut,ok)
      call require(ok,'LES rejection geometry')
      call coarsen_eb_geometry_3d(cut,2,les_cut(1),ok)
      call require(ok,'LES rejection root')
      call extract_eb_geometry_3d(cut,[1,1,1],[4,8,8],les_cut(2),ok)
      call require(ok,'LES rejection child')
    end block
  end if
  if(index(test_mode,'ctu')==1.or.index(test_mode,'sdc')==1) then
    allocate(extended_phi(0:8+4*ctu_padding,0:8+4*ctu_padding,0:8+4*ctu_padding))
    do j=0,8+4*ctu_padding
      extended_phi(:,j,:)=real(max(0,min(8,j-2*ctu_padding)),dp)*.0005_dp-.00125_dp
    end do
    if(index(test_mode,"sdc-les")==1.or.test_mode=="sdc-mol-sources") extended_phi=1
    call build_level_set_eb_geometry_3d(extended_phi,spread(-.001_dp*ctu_padding,1,3), &
      spread(.004_dp+.001_dp*ctu_padding,1,3),extended_finest,ok)
    call require(ok,'extended nodal geometry')
    call build_eb_tree_extended_geometry_3d(tree,g,extended_finest,extended,ok,halo_width=ctu_padding)
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
      if(index(test_mode,'sdc-les')==1) then
        primitive(2)=1+sin(2*acos(-1._dp)*(g(node)%y_lower+(j-.5_dp)*g(node)%dy)/.004_dp)
        primitive(3)=2+cos(2*acos(-1._dp)*(g(node)%x_lower+(i-.5_dp)*g(node)%dx)/.004_dp)
        call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
        call require(ok,'LES nonuniform state')
      end if
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
    if(test_mode=='sdc'.or.test_mode=='sdc-fixed'.or.test_mode=='sdc-reaction-rollback'.or.test_mode=='sdc-force'.or. &
      test_mode=='sdc-force-euler'.or.test_mode=='sdc-recursive'.or. &
      test_mode=='sdc-recursive-state-redist'.or.test_mode=='sdc-recursive-force'.or.test_mode=='sdc-recursive-fixed'.or. &
      test_mode=='sdc-recursive-crossing'.or.test_mode=='sdc-recursive-entering'.or. &
      test_mode=='sdc-recursive-breakup'.or.test_mode=='sdc-recursive-khrt'.or. &
      test_mode=='sdc-recursive-etab'.or.test_mode=='sdc-recursive-etab-film'.or. &
      test_mode=='sdc-recursive-rotating'.or.test_mode=='sdc-recursive-rotating-film'.or. &
      test_mode=='sdc-recursive-children'.or.test_mode=='sdc-recursive-film'.or.test_mode=='sdc-recursive-impact') then
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
  if((test_mode=='sdc-force'.or.test_mode=='sdc-force-euler')) accel=[1000._dp,-2000._dp,3000._dp]
  if(test_mode=='sdc-force-euler') options%mixture_time_scheme='euler'
  saved=field;saved_t=temp;saved_drops=drops;saved_mass=masses
  if(test_mode=='sdc-held') then
    call check_held_coordinator()
    call MPI_Finalize(ierr)
    stop
  endif
  if(test_mode=='sdc-mol-sources') then
    allocate(labels(6));labels='outflow'
    call check_mol_sources()
    call MPI_Finalize(ierr)
    stop
  endif
  if(test_mode=='sdc-reaction-rollback') then
    call check_reaction_rollback()
    call MPI_Finalize(ierr)
    stop
  end if
  if(test_mode=='sdc-fixed') then
    call check_fixed_common()
    call MPI_Finalize(ierr)
    stop
  end if
  if(test_mode=='sdc-recursive-rotating') then
    call check_rotating_recursive_spray()
    call MPI_Finalize(ierr)
    stop
  end if
  if(test_mode=='sdc-recursive-rotating-film') then
    call check_rotating_recursive_film()
    call MPI_Finalize(ierr)
    stop
  end if
  if(index(test_mode,'sdc-les')==1) then
    allocate(labels(6));labels='outflow'
    if(test_mode=='sdc-les-spatial') then
      call check_les_spatial()
    else
      call check_les_spray()
    endif
    call MPI_Finalize(ierr)
    stop
  end if
  if(index(test_mode,'sdc-recursive')==1) then
    if(index(test_mode,'sdc-recursive-impact')==1) then
      call check_recursive_impact()
    else
      call check_recursive_spray()
    end if
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
    if((test_mode=='sdc-force'.or.test_mode=='sdc-force-euler')) then
      call require(any(forcing(2:5)/=0).and.forcing(1)==0.and.all(forcing(6:)==0),'external impulse and work')
    else
      call require(all(forcing==0),'zero external forcing')
    end if
    call require(any(exchange/=0),'internal particle exchange')
    if(rank==0) print *, 'SDC combined residual: ', &
      maxval(abs(after-before-ledger-boundary-chemical)/max(abs(before),abs(after),1e-30_dp))
  end if
  if((test_mode=='sdc-force'.or.test_mode=='sdc-force-euler')) then
    call check_ballistic_force()
    call check_gas_acceleration()
  end if
  if(rank==0) print *, 'PASS: MPI EB coupled reactive step and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine check_held_coordinator()
    type(amr_tree_state_3d) :: held,initial_held,held_temperature
    type(dual_energy_policy) :: policy
    integer(int64) :: counts(2)
    integer :: mode,attempt,p,d(3)
    real(dp) :: prefix,previous_prefix,pressure_work,residual,corrected(nv)
    logical :: mol,common,equal
    character(len=16) :: faces(6)='outflow'
    initial_held=saved_t;held_temperature=saved_t
    ! Independent700K scalar and1000K Eden force native dual-energy policy
    ! to retain UEINT and account for the resulting numerical energy change.
    policy%allow_negative_energy=.false.;policy%eta2=2
    do p=1,2
      if(owners(p)/=rank) cycle
      held_temperature%fields(p)%values=700._dp
      call initialize_nasa7_internal_energy_patch(species,saved%fields(p)%values,held_temperature%fields(p)%values(1,:,:,:), &
        g(p)%volume_fraction,initial_held%fields(p)%values(1,:,:,:),ok)
      call require(ok,'actual UEINT initialized independently from total energy')
    enddo
    allocate(film_settings,films(merge(1,0,rank==owners(1).or.rank==owners(2))))
    allocate(film_mass(2,size(films)))
    film_settings%enabled=.true.;film_settings%evolve=.true.;film_settings%temperature=400
    film_settings%liquid_conductivity=.6_dp;film_settings%minimum_height=1e-15_dp
    if(size(films)>0) then
      films(1)=spray_wall_film([merge(.0025_dp,.0015_dp,rank==owners(1)),.00125_dp,.0015_dp], &
        [0._dp,1._dp,0._dp],1e-10_dp,350._dp,1e-4_dp,1e-5_dp,merge(1,2,rank==owners(1)))
      film_mass(:,1)=[.6_dp,.4_dp]*films(1)%mass
    endif
    saved_films=films;saved_film_mass=film_mass
    do mode=0,3
      mol=mode>=2;common=modulo(mode,2)==0
      do attempt=1,3
        field=saved;temp=saved_t;held=initial_held;drops=saved_drops;masses=saved_mass
        films=saved_films;film_mass=saved_film_mass;last_id=1
        options%breakup_model='etab';options%surface_tension=1e3_dp;options%liquid_viscosity=.001_dp
        options%maximum_multiplicity=50;options%maximum_parcels=100
        options%etab_parcel_factor=-1
        if(attempt==2) options%maximum_parcels=1
        prefix=3.5_dp;pressure_work=4.5_dp;counts=41;boundary=43;ledger=42;chemical=44
        forcing=46;exchange=47;theta=45
        if(attempt==1.and.rank==0) prefix=prefix+1
        previous_prefix=prefix
        call inventory(before)
        if(common) then
          call advance_mpi_eb_spray_sdc(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
            liquids,options,field,temp,drops,masses,zero,spacing,1e-8_dp,2,.false., &
            [.false.,.false.,.false.,.false.],zero,1e-8_dp,1e-12_dp,faces,boundary,ledger,chemical,forcing, &
            theta,steps,last_id,ok,gas_particle_change=exchange,mol_mode=mol,films=films,film_mass=film_mass, &
            film_settings=film_settings,internal_energy=held,energy_policy=policy,energy_correction=prefix, &
            internal_energy_pressure_work_change=pressure_work)
        else
          call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
            liquids,options,field,temp,drops,masses,zero,spacing,1e-8_dp,2,.false., &
            [.false.,.false.,.false.,.false.],zero,1e-8_dp,1e-12_dp,faces,boundary,ledger,chemical,forcing, &
            theta,steps,last_id,ok,gas_particle_change=exchange,mol_mode=mol,films=films,film_mass=film_mass, &
            film_settings=film_settings,level_steps=counts,internal_energy=held,energy_policy=policy, &
            energy_correction=prefix,internal_energy_pressure_work_change=pressure_work)
        endif
        if(attempt<3) then
          call require(.not.ok.and.steps==0.and.last_id==1,'held early/late coupled failure retains physical ID and step count')
          call require(prefix==previous_prefix.and.pressure_work==4.5_dp.and.theta==45.and.all(counts==41), &
            'held source refusal retains pressure/reset/count ledgers')
          call require(all(boundary==43).and.all(ledger==42).and.all(chemical==44).and. &
            all(forcing==46).and.all(exchange==47),'held source refusal retains gas/particle inventory ledgers')
          equal=size(drops)==size(saved_drops).and.size(films)==size(saved_films)
          equal=equal.and.all(masses==saved_mass).and.all(film_mass==saved_film_mass)
          do p=1,size(drops)
            equal=equal.and.drops(p)%temperature==saved_drops(p)%temperature.and. &
              all(drops(p)%position==saved_drops(p)%position).and.all(drops(p)%velocity==saved_drops(p)%velocity)
          enddo
          do p=1,size(films)
            equal=equal.and.films(p)%temperature==saved_films(p)%temperature
          enddo
          call require(equal,'held source refusal retains physical particles and films')
          do p=1,2
            if(owners(p)/=rank) cycle
            call require(all(field%fields(p)%values==saved%fields(p)%values).and. &
              all(temp%fields(p)%values==saved_t%fields(p)%values).and. &
              all(held%fields(p)%values==initial_held%fields(p)%values),'held coupled failure full gas/T/UEINT rollback')
          enddo
          cycle
        endif
        call require(ok,'common/recursive native H particle+film source assembly')
        call require(steps==merge(4,6,common),'actual physical old/final half-kick clock count')
        if(.not.common) call require(all(counts==[1_int64,2_int64]),'held recursive level clock counts')
        call require(ieee_is_finite(prefix).and.ieee_is_finite(pressure_work).and.pressure_work/=4.5_dp, &
          'held pressure/reset ledgers publish only on outer coupled success')
        do p=1,2
          if(owners(p)/=rank) cycle
          call require(all(ieee_is_finite(held%fields(p)%values)),'accepted actual UEINT remains finite')
        enddo
        call inventory(after)
        corrected=after-before-ledger-boundary-chemical-forcing
        corrected(5)=corrected(5)-(prefix-previous_prefix)
        residual=maxval(abs(corrected)/max(abs(before),abs(after),1e-30_dp))
        if(rank==0.and.residual>3e-12_dp) then
          print *, 'Held coupled residual:',common,mol,residual,corrected
          print *, 'Held reset prefix:',previous_prefix,prefix
          print *, 'Held energy inventory/changes:',before(5),after(5),ledger(5),boundary(5),chemical(5),forcing(5)
        endif
        call require(residual<=3e-12_dp,'held coupled gas/particle/film inventory includes separate numerical reset')
        call require(any(exchange/=0),'held fixture has nontrivial physical gas exchange')
        call require(abs(prefix-previous_prefix)>1e-10_dp,'independent held scalar produces separately recorded reset')
        if(rank==0) print *, 'Held particle/film common, MOL, composite residual:',common,mol,residual
      enddo
    enddo
    if(rank==0) print *, 'PASS: actual H coupled gas/particle/film source assembly and early/late rollback'
  end subroutine check_held_coordinator
  subroutine check_mol_sources()
    real(dp) :: gas_before(nv),gas_after(nv),residual
    integer :: mode,iterations
    ! A source-only forced reactor isolates the reacting MOL control path
    ! and completed particle exchange from reaction-model error. This checks
    ! source budgets, not reacting fuel chemistry.
    do mode=0,1
    do iterations=1,3
      field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass;last_id=1
      call inventory(before);call inventory(gas_before,.false.)
      if(mode==0) then
        call advance_mpi_eb_spray_sdc(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
          liquids,options,field,temp,drops,masses,zero,spacing,1e-8_dp,iterations,.true., &
          [.false.,.false.,.false.,.false.],zero,1e-8_dp,1e-12_dp,labels,boundary,ledger,chemical,forcing, &
          theta,steps,last_id,ok,gas_particle_change=exchange,mol_mode=.true.,evaluate_reaction=source_only_reactor)
      else
        call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
          liquids,options,field,temp,drops,masses,zero,spacing,1e-8_dp,iterations,.true., &
          [.false.,.false.,.false.,.false.],zero,1e-8_dp,1e-12_dp,labels,boundary,ledger,chemical,forcing, &
          theta,steps,last_id,ok,gas_particle_change=exchange,mol_mode=.true.,evaluate_reaction=source_only_reactor)
      endif
      call require(ok,'reacting-control MOL particle source step')
      call inventory(after);call inventory(gas_after,.false.)
      residual=maxval(abs(after-before-ledger-boundary-chemical-forcing)/max(abs(before),abs(after),1e-30_dp))
      call require(residual<=3e-12_dp,'MOL completed particle exchange preserves composite inventory')
      call require(all(abs(gas_after-gas_before-ledger-boundary-chemical-forcing-exchange)<= &
        3e-12_dp*max(abs(gas_before),abs(gas_after),1e-30_dp)),'MOL gas particle exchange ledger')
      call require(any(exchange(6:)>0),'MOL source test has nontrivial evaporation')
      if(rank==0) print *, 'MOL common/recursive, iterations, composite residual:',mode,iterations,residual
    enddo
    enddo
  end subroutine check_mol_sources
  subroutine source_only_reactor(sp,rx,e0,external_y,external_e,duration,rtol,atol,densities,temp,count,good,maximum_steps)
    type(nasa7_species), intent(in) :: sp(:)
    type(elementary_reaction), intent(in) :: rx(:)
    real(dp), intent(in) :: e0,external_y(:),external_e,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),temp
    integer, intent(out) :: count
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    type(elementary_reaction) :: no_reactions(0)
    real(dp) :: initial_density(size(densities)),final_density(size(densities))
    integer :: limiting
    initial_density=densities
    call advance_forced_reactor_interval(sp,no_reactions,e0,external_y,external_e,duration,rtol,atol, &
      densities,temp,count,good,maximum_steps=maximum_steps)
    if(.not.good) then
      final_density=initial_density+duration*external_y
      limiting=minloc(final_density,dim=1)
      print *, 'Source-only reactor rejected rank/species/start/change/final:',rank,sp(limiting)%name, &
        initial_density(limiting),duration*external_y(limiting),final_density(limiting)
    endif
  end subroutine
  subroutine check_reaction_rollback()
    type(amr_tree_state_3d) :: heat,ir
    integer :: mode,q,n(3),global_calls
    do mode=1,4
      field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass;last_id=1
      boundary=43;ledger=42;chemical=44;forcing=46;exchange=47;theta=45
      allocate(heat%fields(size(owners)),ir%fields(size(owners)))
      do q=1,size(owners)
        if(owners(q)/=rank) cycle
        n=tree%nodes(q)%shape
        allocate(heat%fields(q)%values(1,n(1),n(2),n(3)),source=53._dp)
        allocate(ir%fields(q)%values(size(species)+1,n(1),n(2),n(3)),source=0._dp)
      end do
      rejected_reactor_calls=0
      if(mode<=2) then
        call advance_mpi_eb_spray_sdc(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
          liquids,options,field,temp,drops,masses,zero,spacing,1e-8_dp,2,.true., &
          [.false.,.false.,.false.,.false.],zero,2e-7_dp,1e-12_dp,labels,boundary,ledger,chemical, &
          forcing,theta,steps,last_id,ok,gas_particle_change=exchange,reaction_source=ir,heat_release=heat, &
          evaluate_reaction=reject_reaction,mol_mode=mod(mode,2)==0)
      else
        call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
          liquids,options,field,temp,drops,masses,zero,spacing,1e-8_dp,2,.true., &
          [.false.,.false.,.false.,.false.],zero,2e-7_dp,1e-12_dp,labels,boundary,ledger,chemical, &
          forcing,theta,steps,last_id,ok,gas_particle_change=exchange,reaction_source=ir,heat_release=heat, &
          evaluate_reaction=reject_reaction,mol_mode=mod(mode,2)==0)
      endif
      call MPI_Allreduce(rejected_reactor_calls,global_calls,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      call require(global_calls>0,'reactor rejection reached after particle callback')
      call require(.not.ok.and.steps==0.and.last_id==1,'gas rejection propagates through particle wrapper')
      call require(all(boundary==43).and.all(ledger==42).and.all(chemical==44).and. &
        all(forcing==46).and.all(exchange==47).and.theta==45,'reactor rejection preserves budgets')
      call require(size(drops)==size(saved_drops).and.all(masses==saved_mass),'reactor rejection preserves population')
      if(rank==nranks-1) then
        call require(all(drops(1)%position==saved_drops(1)%position).and. &
          all(drops(1)%velocity==saved_drops(1)%velocity).and.drops(1)%mass==saved_drops(1)%mass.and. &
          drops(1)%temperature==saved_drops(1)%temperature,'reactor rejection preserves parcel state')
      endif
      do q=1,size(owners)
        if(owners(q)/=rank) cycle
        call require(all(field%fields(q)%values==saved%fields(q)%values).and. &
          all(temp%fields(q)%values==saved_t%fields(q)%values),'reactor rejection preserves gas')
        call require(all(heat%fields(q)%values==53).and.all(ir%fields(q)%values==0), &
          'reactor rejection preserves rates and heat')
      enddo
      deallocate(heat%fields,ir%fields)
    enddo
    if(rank==0) print *, 'PASS: common and recursive CTU/MOL reactor rejection rollback'
  end subroutine

  subroutine reject_reaction(sp,rx,e0,external_y,external_e,duration,rtol,atol,densities,temp,count,good,maximum_steps)
    type(nasa7_species),intent(in) :: sp(:)
    type(elementary_reaction),intent(in) :: rx(:)
    real(dp),intent(in) :: e0,external_y(:),external_e,duration,rtol,atol
    real(dp),intent(inout) :: densities(:),temp
    integer,intent(out) :: count
    logical,intent(out) :: good
    integer,optional,intent(in) :: maximum_steps
    rejected_reactor_calls=rejected_reactor_calls+1
    count=0;good=.false.
  end subroutine

  subroutine check_fixed_common()
    type(amr_tree_state_3d) :: heat
    integer :: scheme,iterations,q,n(3)
    real(dp) :: kinetic
    options%fixed_parts=.true.;accel=[1000._dp,-2000._dp,3000._dp]
    do scheme=1,2
      options%mixture_time_scheme='conservative'
      if(scheme==2) options%mixture_time_scheme='euler'
      iterations=merge(1,3,scheme==1)
      field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass;last_id=1
      allocate(heat%fields(size(owners)))
      do q=1,size(owners)
        if(owners(q)/=rank) cycle
        n=tree%nodes(q)%shape
        allocate(heat%fields(q)%values(1,n(1),n(2),n(3)),source=53._dp)
      end do
      call inventory(before)
      call advance_mpi_eb_spray_sdc(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
        liquids,options,field,temp,drops,masses,zero,spacing,1e-8_dp,iterations,.true., &
        [.true.,.true.,.true.,.true.],zero,2e-7_dp,1e-12_dp,labels,boundary,ledger,chemical, &
        forcing,theta,steps,last_id,ok,gas_particle_change=exchange,particle_acceleration=accel,heat_release=heat)
      call require(ok.and.steps>0,'fixed common SDC step')
      do q=1,size(owners)
        if(owners(q)/=rank) cycle
        call require(all(ieee_is_finite(heat%fields(q)%values)),'common spray finite stage heat')
        call require(any(heat%fields(q)%values/=53._dp),'common spray publishes stage heat')
        call require(all(pack(heat%fields(q)%values(1,:,:,:),g(q)%volume_fraction==0)==0), &
          'common spray covered heat is zero')
      end do
      deallocate(heat%fields)
      call inventory(after)
      call require(all(abs(after-before-ledger-boundary-chemical-forcing)<= &
        3e-12_dp*max(abs(before),abs(after),1e-30_dp)),'fixed common SDC composite inventory')
      kinetic=-.5_dp*1e-10_dp*sum([20._dp,-20._dp,5._dp]**2)
      call require(abs(forcing(5)/kinetic-1)<2e-13_dp,'fixed common SDC constraint counted once')
      if(rank==nranks-1) then
        call require(all(drops(1)%velocity==0).and.all(drops(1)%position==saved_drops(1)%position), &
          'fixed common SDC physical stationary state')
        call require(drops(1)%temperature>350.and.drops(1)%mass<saved_drops(1)%mass,'fixed common SDC evaporation')
      end if
    end do
  end subroutine check_fixed_common
  subroutine check_rotating_recursive_spray()
    type(rotating_frame_3d) :: frame
    type(amr_tree_state_3d) :: before_gas,before_t
    real(dp) :: initial(nv),final(nv),frame_change
    integer :: a
    frame_change=31
    frame%axis=3;frame%angular_speed=1000._dp;frame%axis_point=[.001_dp,.001_dp,.001_dp]
    field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass;last_id=1
    call convert_mpi_amr_tree_rotating_energy(MPI_COMM_WORLD,tree,owners,field,zero,spacing,frame,.true.,ok,g)
    call require(ok,'recursive rotating initial energy')
    before_gas=field;before_t=temp
    boundary=11;ledger=12;chemical=13;forcing=14;exchange=16;theta=15
    if(rank==0) then
      call advance_mpi_eb_spray_sdc(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
        liquids,options,field,temp,drops,masses,zero,spacing,1.e-8_dp,2,.false., &
        [.false.,.false.,.false.,.false.],zero,2.e-7_dp,1.e-12_dp,labels,boundary,ledger,chemical, &
        forcing,theta,steps,last_id,ok,mol_mode=.true.,frame=frame,gas_particle_change=exchange,frame_energy_change=frame_change)
    else
      call advance_mpi_eb_spray_sdc(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
        liquids,options,field,temp,drops,masses,zero,spacing,1.e-8_dp,2,.false., &
        [.false.,.false.,.false.,.false.],zero,2.e-7_dp,1.e-12_dp,labels,boundary,ledger,chemical, &
        forcing,theta,steps,last_id,ok,mol_mode=.true.,gas_particle_change=exchange,frame_energy_change=frame_change)
    end if
    call require(.not.ok.and.steps==0.and.last_id==1.and.all(boundary==11).and.all(ledger==12).and. &
      all(chemical==13).and.all(forcing==14).and.all(exchange==16).and.theta==15.and.frame_change==31, &
      'common rotating frame presence rejects collectively')
    do a=1,size(owners)
      if(owners(a)/=rank) cycle
      call require(all(field%fields(a)%values==before_gas%fields(a)%values).and. &
        all(temp%fields(a)%values==before_t%fields(a)%values),'common rotating rejection preserves gas')
    end do
    call rotating_inventory(frame,initial)
    call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
      liquids,options,field,temp,drops,masses,zero,spacing,1.e-8_dp,2,.false., &
      [.false.,.false.,.false.,.false.],zero,2.e-7_dp,1.e-12_dp,labels,boundary,ledger,chemical, &
      forcing,theta,steps,last_id,ok,mol_mode=.true.,frame=frame,gas_particle_change=exchange,frame_energy_change=frame_change)
    call require(ok.and.steps>0,'recursive rotating particle step')
    call rotating_inventory(frame,final)
    final(5)=final(5)-frame_change
    call require(ieee_is_finite(frame_change).and.frame_change/=31,"rotating coordinate ledger published")
    call require(all(abs(final-initial-boundary-ledger-chemical-forcing)<= &
      3.e-12_dp*max(abs(initial),abs(final),1.e-30_dp)),'recursive rotating composite inventory')
    call require(any(exchange/=0),'recursive rotating particle exchange')
    do a=1,size(owners)
      if(owners(a)/=rank) cycle
      call require(all(ieee_is_finite(field%fields(a)%values)),'recursive rotating finite gas')
      call require(any(field%fields(a)%values/=before_gas%fields(a)%values),'recursive rotating gas evolved')
      call require(all(ieee_is_finite(temp%fields(a)%values)),'recursive rotating finite temperature')
    end do
    if(rank==0) print *, 'PASS: MPI recursive rotating particle coupling'
  end subroutine check_rotating_recursive_spray

  subroutine rotating_inventory(frame,result)
    type(rotating_frame_3d), intent(in) :: frame
    real(dp), intent(out) :: result(:)
    real(dp) :: local_potential,global_potential,potential
    integer :: parcel,status
    logical :: good
    call inventory(result)
    local_potential=0
    do parcel=1,size(drops)
      call rotating_potential_3d(frame,drops(parcel)%position,potential,good)
      call require(good,'particle rotating potential')
      local_potential=local_potential+drops(parcel)%mass*drops(parcel)%multiplicity*potential
    end do
    if(allocated(films)) then
      do parcel=1,size(films)
        call rotating_potential_3d(frame,films(parcel)%position,potential,good)
        call require(good,'film rotating potential')
        local_potential=local_potential+films(parcel)%mass*potential
      end do
    end if
    call MPI_Allreduce(local_potential,global_potential,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS,'particle rotating potential reduction')
    result(5)=result(5)-global_potential
  end subroutine rotating_inventory

  subroutine check_rotating_recursive_film()
    type(rotating_frame_3d) :: frame
    real(dp) :: initial(nv),final(nv),frame_change
    frame_change=31
    frame%axis=3;frame%angular_speed=1000._dp;frame%axis_point=[.001_dp,.001_dp,.001_dp]
    field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass;last_id=1
    if(rank==nranks-1) drops(1)%position=[.0012_dp,.0025_dp,.0015_dp]
    allocate(film_settings,films(merge(1,0,rank==0.or.rank==nranks-1)))
    allocate(film_mass(2,size(films)))
    film_settings%enabled=.true.;film_settings%evolve=.true.;film_settings%temperature=400
    film_settings%liquid_conductivity=.6_dp;film_settings%minimum_height=1.e-15_dp
    if(size(films)>0) then
      films(1)=spray_wall_film([merge(.0025_dp,.0015_dp,rank==0),.00125_dp,.0015_dp], &
        [0._dp,1._dp,0._dp],1.e-10_dp,350._dp,1.e-4_dp,1.e-5_dp,merge(1,2,rank==0))
      film_mass(:,1)=[.6_dp,.4_dp]*films(1)%mass
    end if
    call convert_mpi_amr_tree_rotating_energy(MPI_COMM_WORLD,tree,owners,field,zero,spacing,frame,.true.,ok,g)
    call require(ok,'recursive film initial rotation')
    call rotating_inventory(frame,initial)
    call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,tree,owners,g,extended,phi,species,reactions,transport, &
      liquids,options,field,temp,drops,masses,zero,spacing,1.e-8_dp,2,.false., &
      [.false.,.false.,.false.,.false.],zero,2.e-7_dp,1.e-12_dp,labels,boundary,ledger,chemical, &
      forcing,theta,steps,last_id,ok,mol_mode=.true.,frame=frame,gas_particle_change=exchange, &
      frame_energy_change=frame_change, &
      films=films,film_mass=film_mass,film_settings=film_settings,film_steps=film_steps)
    call require(ok.and.steps>0.and.film_steps>0,'recursive rotating particle and film step')
    call rotating_inventory(frame,final)
    final(5)=final(5)-frame_change
    call require(ieee_is_finite(frame_change).and.frame_change/=31,"rotating coordinate ledger published")
    call require(all(abs(final-initial-boundary-ledger-chemical-forcing)<= &
      3.e-12_dp*max(abs(initial),abs(final),1.e-30_dp)),'recursive rotating film composite inventory')
    call require(any(exchange/=0),'recursive rotating film source')
    if(rank==0) print *, 'PASS: MPI recursive rotating particle and film coupling'
  end subroutine check_rotating_recursive_film

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
  subroutine check_recursive_impact()
    integer :: attempt,live,film_count
    integer(int64) :: counts(2)
    logical :: rejected
    options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    if(rank==nranks-1) then
      drops(1)=spray_parcel(1_int64,[.0012_dp,.001250001_dp,.0015_dp], &
        [0._dp,-100._dp,0._dp],1e-9_dp,350._dp,1._dp)
      masses(:,1)=[.25_dp,.75_dp]*drops(1)%mass
    end if
    if(index(test_mode,'periodic')>0.and.rank==nranks-1) then
      drops(1)%position(3)=.003999999_dp;drops(1)%velocity(3)=100._dp
    end if
    allocate(film_settings,films(merge(1,0,rank==0.or.rank==nranks-1)))
    allocate(film_mass(2,size(films)))
    film_settings%enabled=.true.;film_settings%evolve=.true.;film_settings%temperature=350
    film_settings%liquid_conductivity=.6_dp;film_settings%minimum_height=1e-15_dp
    if(size(films)>0) then
      films(1)=spray_wall_film([merge(.0025_dp,.0015_dp,rank==0),.00125_dp,.0015_dp], &
        [0._dp,1._dp,0._dp],1e-10_dp,350._dp,1e-4_dp,1e-5_dp,1)
      film_mass(:,1)=[.6_dp,.4_dp]*films(1)%mass
    end if
    saved_drops=drops;saved_mass=masses;saved_films=films;saved_film_mass=film_mass
    call inventory(before)
    do attempt=1,2
      options%maximum_parcels=merge(1,10000,attempt==1)
      boundary=17;ledger=19;chemical=23;forcing=29;theta=31;exchange=37;counts=41
      call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,tree,owners,g,extended,phi, &
        species,reactions,transport,liquids,options,field,temp,drops,masses,zero,spacing, &
        1e-8_dp,2,.true.,[.true.,.true.,.true.,.true.],zero,2e-7_dp,1e-12_dp,labels, &
        boundary,ledger,chemical,forcing,theta,steps,last_id,ok,gas_particle_change=exchange,level_steps=counts, &
        films=films,film_mass=film_mass,film_settings=film_settings,film_steps=film_steps,film_impacts=.true., &
        model_rejected=rejected)
      if(attempt==1) then
        call require(.not.ok.and..not.rejected.and.steps==0.and.film_steps==0.and.last_id==1, &
          'recursive splash capacity failure rolls back IDs and clocks')
        call require(size(drops)==size(saved_drops).and.all(masses==saved_mass).and. &
          size(films)==size(saved_films).and.all(film_mass==saved_film_mass),'recursive impact inventory rollback')
        call require(all(films%temperature==saved_films%temperature).and.all(films%mass==saved_films%mass), &
          'recursive impact film thermal rollback')
        call require(all(boundary==17).and.all(ledger==19).and.all(chemical==23).and.all(forcing==29).and. &
          theta==31.and.all(exchange==37).and.all(counts==41),'recursive impact ledger rollback')
        do node=1,2
          if(owners(node)/=rank) cycle
          call require(all(field%fields(node)%values==saved%fields(node)%values).and. &
            all(temp%fields(node)%values==saved_t%fields(node)%values),'recursive impact gas rollback')
        end do
      else
        call require(ok.and..not.rejected.and.steps==6.and.film_steps==6.and.all(counts==[1_int64,2_int64]), &
          'recursive impact thermal coupled interval')
      end if
    end do
    live=size(drops);film_count=size(films)
    call MPI_Allreduce(MPI_IN_PLACE,live,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,film_count,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(live>1.and.last_id>1.and.film_count>2,'recursive splash children and newborn films')
    call inventory(after)
    if(rank==0) print *, 'Recursive impact residual / particles / films: ', &
      maxval(abs(after-before-ledger-boundary-chemical-forcing)/max(abs(before),abs(after),1e-30_dp)),live,film_count
    call require(all(abs(after-before-ledger-boundary-chemical-forcing)<= &
      3e-12_dp*max(abs(before),abs(after),1e-30_dp)),'recursive gas particle film wall conservation')
  end subroutine check_recursive_impact
  subroutine check_les_spatial()
    type(amr_tree_state_3d) :: halo,halo_t,rate,spatial,old_rate,old_spatial,volume_halo,metric,short_halo
    type(amr_tree_state_3d) :: rotated_halo,rotation_before,rotated_rate,rotated_spatial
    type(amr_tree_state_3d) :: rotated_state,temporal_rate,rotated_temporal_rate
    type(amr_tree_state_3d) :: ancestors,rotated_ancestors,ancestor_t
    type(rotating_frame_3d) :: rotation
    type(eb_integrated_flux_3d), allocatable :: flux(:),rotated_flux(:),temporal_flux(:),rotated_temporal_flux(:)
    type(eb_integrated_flux_3d), allocatable :: previous_temporal_flux(:)
    character(len=16) :: models(4)=['smagorinsky     ','wale            ','vreman          ','dynamic         ']
    real(dp) :: diffusion,old_diffusion,rotated_diffusion,divergence(nv),error,scale,h(3)
    integer :: mode,q,n(3),x,z,w
    allocate(metric%fields(2))
    do q=1,2
      if(owners(q)/=rank) cycle
      n=tree%nodes(q)%shape
      allocate(metric%fields(q)%values(1,-9:n(1)+10,-9:n(2)+10,-9:n(3)+10),source=1._dp)
    enddo
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,g,field,temp,halo,halo_t,MPI_COMM_WORLD,ok, &
      halo_width=10,boundary_faces=labels,ghost_volume_fraction=volume_halo,target_volume_fraction=metric)
    call require(ok,'LES spatial supplied grow10')
    rotated_halo=halo;rotation%axis=2;rotation%angular_speed=600
    rotation%axis_point=[.001_dp,.0007_dp,-.0003_dp]
    call convert_mpi_amr_tree_rotating_energy(MPI_COMM_WORLD,tree,owners,rotated_halo,zero,spacing, &
      rotation,.true.,ok,g,halo_width=10)
    call require(ok,'LES published rotating halo setup')
    rotation_before=rotated_halo
    do mode=1,4
      diffusion=-1
      call compute_mpi_regular_eb_les_spatial_sdc_3d(tree,owners,extended,species,halo,1e-7_dp,-1, &
        .12_dp,.7_dp,.7_dp,rate,spatial,flux,diffusion,MPI_COMM_WORLD,ok,model=models(mode))
      call require(ok.and.diffusion>0,'LES spatial models')
      call compute_mpi_regular_eb_les_spatial_sdc_3d(tree,owners,extended,species,rotated_halo,1e-7_dp,-1, &
        .12_dp,.7_dp,.7_dp,rotated_rate,rotated_spatial,rotated_flux,rotated_diffusion,MPI_COMM_WORLD,ok, &
        model=models(mode),frame=rotation)
      call require(ok,'LES rotating spatial models')
      call require(abs(rotated_diffusion-diffusion)<2e-11_dp*max(diffusion,tiny(1._dp)), &
        'LES rotating diffusion bound')
      do q=1,2
        if(owners(q)/=rank) cycle
        n=tree%nodes(q)%shape;h=[g(q)%dx,g(q)%dy,g(q)%dz]
        call require(all(rate%fields(q)%values==spatial%fields(q)%values(:,7:n(1)+6,7:n(2)+6,7:n(3)+6)), &
          'LES spatial exact valid crop')
        call require(all(rotated_halo%fields(q)%values==rotation_before%fields(q)%values), &
          'LES published halo immutable')
        scale=max(1._dp,maxval(abs(spatial%fields(q)%values)))
        call require(maxval(abs(rotated_spatial%fields(q)%values-spatial%fields(q)%values))<2e-11_dp*scale, &
          'LES rotating spatial rates match physical control')
        scale=max(tiny(1._dp),maxval(abs(flux(q)%x)),maxval(abs(flux(q)%y)),maxval(abs(flux(q)%z)))
        call require(maxval(abs(rotated_flux(q)%x-flux(q)%x))<2e-11_dp*scale.and. &
          maxval(abs(rotated_flux(q)%y-flux(q)%y))<2e-11_dp*scale.and. &
          maxval(abs(rotated_flux(q)%z-flux(q)%z))<2e-11_dp*scale,'LES rotating integrated flux control')
        error=0;scale=max(1._dp,maxval(abs(rate%fields(q)%values)))
        do w=1,n(3)
        do z=1,n(2)
        do x=1,n(1)
          divergence=-(flux(q)%x(:,x,z,w)-flux(q)%x(:,x-1,z,w))/h(1)/1e-7_dp &
            -(flux(q)%y(:,x,z,w)-flux(q)%y(:,x,z-1,w))/h(2)/1e-7_dp &
            -(flux(q)%z(:,x,z,w)-flux(q)%z(:,x,z,w-1))/h(3)/1e-7_dp
          error=max(error,maxval(abs(divergence-rate%fields(q)%values(:,x,z,w))))
        enddo
        enddo
        enddo
        call require(error<2e-12_dp*scale,'LES spatial flux divergence')
        if(mode==1) call require(any(spatial%fields(q)%values(:,4:6,4:n(2)+9,4:n(3)+9)/=0), &
          'LES active grow3 source')
        ! The clamped x-outflow slab has only du_x/dy; Vreman's
        ! rank-one velocity-gradient invariant and WALE's squared-gradient
        ! deviator are both zero in this region.
        if(mode==2.or.mode==3) call require(all(spatial%fields(q)%values(:,4:6,4:n(2)+9,4:n(3)+9)==0), &
          'LES WALE/Vreman rank-one halo control')
      enddo
      old_rate=rate;old_spatial=spatial;old_diffusion=diffusion
      call compute_mpi_regular_eb_les_spatial_sdc_3d(tree,owners,extended,species,halo,-1._dp,-1, &
        .12_dp,.7_dp,.7_dp,rate,spatial,flux,diffusion,MPI_COMM_WORLD,ok,model=models(mode))
      call require(.not.ok.and.diffusion==old_diffusion,'LES spatial failure')
      do q=1,2
        if(owners(q)/=rank) cycle
        call require(all(rate%fields(q)%values==old_rate%fields(q)%values).and. &
          all(spatial%fields(q)%values==old_spatial%fields(q)%values),'LES spatial failure rollback')
      enddo
    enddo
    rotated_state=field
    call convert_mpi_amr_tree_rotating_energy(MPI_COMM_WORLD,tree,owners,rotated_state,zero,spacing, &
      rotation,.true.,ok,g)
    call require(ok,'LES rotating temporal state setup')
    call compute_mpi_regular_eb_les_sdc_3d(tree,owners,g,species,field,temp,spacing,labels, &
      1e-7_dp,-1,.12_dp,.7_dp,.7_dp,temporal_rate,temporal_flux,diffusion,MPI_COMM_WORLD,ok,model='wale')
    call require(ok,'LES temporal physical control')
    call compute_mpi_regular_eb_les_sdc_3d(tree,owners,g,species,rotated_state,temp,spacing,labels, &
      1e-7_dp,-1,.12_dp,.7_dp,.7_dp,rotated_temporal_rate,rotated_temporal_flux,rotated_diffusion, &
      MPI_COMM_WORLD,ok,model='wale',frame=rotation,domain_lower=zero)
    call require(ok,'LES temporal rotating state')
    call require(abs(rotated_diffusion-diffusion)<2e-11_dp*max(diffusion,tiny(1._dp)), &
      'LES temporal rotating diffusivity')
    do q=1,2
      if(owners(q)/=rank) cycle
      scale=max(1._dp,maxval(abs(temporal_rate%fields(q)%values)))
      call require(maxval(abs(rotated_temporal_rate%fields(q)%values-temporal_rate%fields(q)%values))<2e-11_dp*scale, &
        'LES temporal rotating rates match physical control')
      scale=max(tiny(1._dp),maxval(abs(temporal_flux(q)%x)),maxval(abs(temporal_flux(q)%y)), &
        maxval(abs(temporal_flux(q)%z)))
      call require(maxval(abs(rotated_temporal_flux(q)%x-temporal_flux(q)%x))<2e-11_dp*scale.and. &
        maxval(abs(rotated_temporal_flux(q)%y-temporal_flux(q)%y))<2e-11_dp*scale.and. &
        maxval(abs(rotated_temporal_flux(q)%z-temporal_flux(q)%z))<2e-11_dp*scale, &
        'LES temporal rotating flux matches physical control')
    enddo
    allocate(ancestors%fields(2),ancestor_t%fields(2))
    if(owners(1)==rank) then
      ancestors%fields(1)%values=field%fields(1)%values
      ancestor_t%fields(1)%values=temp%fields(1)%values
    endif
    rotated_ancestors=ancestors
    call convert_mpi_amr_tree_rotating_energy(MPI_COMM_WORLD,tree,owners,rotated_ancestors,zero,spacing, &
      rotation,.true.,ok,g,maximum_depth=0)
    call require(ok,'LES rotating temporal ancestor setup')
    call compute_mpi_regular_eb_les_sdc_3d(tree,owners,g,species,field,temp,spacing,labels, &
      1e-7_dp,1,.12_dp,.7_dp,.7_dp,temporal_rate,temporal_flux,diffusion,MPI_COMM_WORLD,ok, &
      ancestors,ancestors,ancestor_t,ancestor_t,[.5_dp],model='wale')
    call require(ok,'LES temporal level-one physical control')
    call compute_mpi_regular_eb_les_sdc_3d(tree,owners,g,species,rotated_state,temp,spacing,labels, &
      1e-7_dp,1,.12_dp,.7_dp,.7_dp,rotated_temporal_rate,rotated_temporal_flux,rotated_diffusion, &
      MPI_COMM_WORLD,ok,rotated_ancestors,rotated_ancestors,ancestor_t,ancestor_t,[.5_dp], &
      model='wale',frame=rotation,domain_lower=zero)
    call require(ok,'LES temporal level-one rotating state')
    if(owners(2)==rank) then
      scale=max(1._dp,maxval(abs(temporal_rate%fields(2)%values)))
      call require(maxval(abs(rotated_temporal_rate%fields(2)%values-temporal_rate%fields(2)%values))<2e-11_dp*scale, &
        'LES temporal level-one rotating rates match physical control')
      scale=max(tiny(1._dp),maxval(abs(temporal_flux(2)%x)),maxval(abs(temporal_flux(2)%y)), &
        maxval(abs(temporal_flux(2)%z)))
      call require(maxval(abs(rotated_temporal_flux(2)%x-temporal_flux(2)%x))<2e-11_dp*scale.and. &
        maxval(abs(rotated_temporal_flux(2)%y-temporal_flux(2)%y))<2e-11_dp*scale.and. &
        maxval(abs(rotated_temporal_flux(2)%z-temporal_flux(2)%z))<2e-11_dp*scale, &
        'LES temporal level-one rotating flux matches physical control')
    endif
    if(rank==1)rotation%angular_speed=601
    old_rate=rotated_temporal_rate;old_diffusion=rotated_diffusion
    previous_temporal_flux=rotated_temporal_flux
    call compute_mpi_regular_eb_les_sdc_3d(tree,owners,g,species,rotated_state,temp,spacing,labels, &
      1e-7_dp,-1,.12_dp,.7_dp,.7_dp,rotated_temporal_rate,rotated_temporal_flux,rotated_diffusion, &
      MPI_COMM_WORLD,ok,model='wale',frame=rotation,domain_lower=zero)
    call require(.not.ok.and.rotated_diffusion==old_diffusion, &
      'LES temporal rank-disagreeing rotating frame rejected')
    do q=1,2
      if(owners(q)/=rank) cycle
      call require(all(rotated_temporal_rate%fields(q)%values==old_rate%fields(q)%values), &
        'LES temporal frame rejection preserves rates')
      call require(all(rotated_temporal_flux(q)%x==previous_temporal_flux(q)%x).and. &
        all(rotated_temporal_flux(q)%y==previous_temporal_flux(q)%y).and. &
        all(rotated_temporal_flux(q)%z==previous_temporal_flux(q)%z), &
        'LES temporal frame rejection preserves flux')
    enddo
    call compute_mpi_regular_eb_les_spatial_sdc_3d(tree,owners,extended,species,rotated_halo,1e-7_dp,-1, &
      .12_dp,.7_dp,.7_dp,rotated_rate,rotated_spatial,rotated_flux,rotated_diffusion,MPI_COMM_WORLD,ok, &
      model='dynamic',frame=rotation)
    call require(.not.ok,'LES rank-disagreeing rotating frame rejected')
    rotation%angular_speed=600
    allocate(short_halo%fields(2))
    do q=1,2
      if(owners(q)/=rank) cycle
      n=tree%nodes(q)%shape
      short_halo%fields(q)%values=halo%fields(q)%values(:,5:n(1)+16,5:n(2)+16,5:n(3)+16)
    enddo
    old_rate=rate;old_spatial=spatial;old_diffusion=diffusion
    call compute_mpi_regular_eb_les_spatial_sdc_3d(tree,owners,extended,species,short_halo,1e-7_dp,-1, &
      .12_dp,.7_dp,.7_dp,rate,spatial,flux,diffusion,MPI_COMM_WORLD,ok,model='dynamic')
    call require(.not.ok.and.diffusion==old_diffusion,'LES dynamic missing filter halo rejected')
    call compute_mpi_regular_eb_les_spatial_sdc_3d(tree,owners,extended,species,halo,1e-7_dp,1, &
      .12_dp,.7_dp,.7_dp,rate,spatial,flux,diffusion,MPI_COMM_WORLD,ok,model='wale')
    call require(ok,'LES selected level')
    if(rank==owners(1))then
      call require(all(rate%fields(1)%values==0).and..not.allocated(spatial%fields(1)%values), &
        'LES inactive level rate contract')
    endif
  end subroutine check_les_spatial
  subroutine check_les_spray()
    character(len=16) :: models(4)
    integer :: m,peer
    real(dp) :: normalized,scale(nv),delta,global_delta
    type(amr_tree_state_3d) :: no_les
    models=[character(len=16)::'smagorinsky','wale','vreman','dynamic']
    call inventory(before)
    call les_step()
    call require(ok,'gas particle reference without LES')
    no_les=field
    do m=1,4
      les_model=models(m);field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass;last_id=1
      les_calls=0;les_largest=0
      call les_step(les_source)
      call require(ok.and.les_largest>0,'active gas particle LES')
      delta=0
      do peer=1,size(owners)
        if(owners(peer)/=rank) cycle
        delta=max(delta,maxval(abs(field%fields(peer)%values-no_les%fields(peer)%values)))
      end do
      call MPI_Allreduce(delta,global_delta,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call require(global_delta>0,'LES changes coupled gas state')
      call inventory(after)
      scale=max(abs(before),abs(after),1e-20_dp)
      scale(2:4)=max(scale(2:4),sqrt(before(1)*before(5)))
      normalized=maxval(abs(after-before-boundary-ledger-chemical-forcing)/scale)
      call require(normalized<5e-12_dp,'gas particle LES composite inventory')
      if(rank==0) print *, 'LES particle residual ',trim(les_model),normalized
    end do
    field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass;last_id=1
    boundary=17;ledger=19;chemical=23;forcing=29;theta=31;exchange=37
    les_failure=.true.;les_calls=0
    call les_step(les_source)
    call require(.not.ok,'cut-cell LES rejected during particle step')
    call require(all(masses==saved_mass).and.size(drops)==size(saved_drops),'LES particle mass rollback')
    do peer=1,size(drops)
      call require(all(drops(peer)%position==saved_drops(peer)%position),'LES particle position rollback')
      call require(all(drops(peer)%velocity==saved_drops(peer)%velocity),'LES particle velocity rollback')
      call require(drops(peer)%temperature==saved_drops(peer)%temperature,'LES particle temperature rollback')
    end do
    do peer=1,size(owners)
      if(owners(peer)/=rank) cycle
      call require(all(field%fields(peer)%values==saved%fields(peer)%values),'LES gas rollback')
      call require(all(temp%fields(peer)%values==saved_t%fields(peer)%values),'LES temperature rollback')
    end do
    call require(last_id==1.and.all(boundary==17).and.all(ledger==19).and.all(chemical==23).and. &
      all(forcing==29).and.theta==31.and.all(exchange==37),'LES diagnostics rollback')
    if(rank==0) print *, 'PASS: four LES models, particle conservation and cut-cell rejection rollback'
  end subroutine
  subroutine les_step(evaluator)
    procedure(mpi_eb_sdc_diffusion_source), optional :: evaluator
    if(index(test_mode,'subcycle')>0) then
      call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,tree,owners,g,extended,phi, &
        species,reactions,transport,liquids,options,field,temp,drops,masses,zero,spacing, &
        1e-8_dp,2,.false.,[.true.,.true.,.true.,.true.],zero,2e-7_dp,1e-12_dp,labels, &
        boundary,ledger,chemical,forcing,theta,steps,last_id,ok,gas_particle_change=exchange,evaluate_diffusion=evaluator)
    else
      call advance_mpi_eb_spray_sdc(MPI_COMM_WORLD,tree,owners,g,extended,phi, &
        species,reactions,transport,liquids,options,field,temp,drops,masses,zero,spacing, &
        1e-8_dp,2,.false.,[.true.,.true.,.true.,.true.],zero,2e-7_dp,1e-12_dp,labels, &
        boundary,ledger,chemical,forcing,theta,steps,last_id,ok,gas_particle_change=exchange,evaluate_diffusion=evaluator)
    end if
  end subroutine
  subroutine les_source(state,temperature,interval,depth,rate,flux,limiter,comm,good, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,fractions)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: interval
    integer, intent(in) :: depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: limiter
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: fractions(:)
    type(eb_geometry_3d) :: geometry(2)
    real(dp) :: diffusivity
    les_calls=les_calls+1;geometry=g
    if(les_failure.and.les_calls==2) geometry=les_cut
    call compute_mpi_regular_eb_les_sdc_3d(tree,owners,geometry,species,state,temperature,spacing,labels, &
      interval,depth,.17_dp,.7_dp,.7_dp,rate,flux,diffusivity,comm,good, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,fractions,model=les_model)
    if(good) then
      limiter=1;les_largest=max(les_largest,diffusivity)
    end if
  end subroutine
  subroutine check_recursive_spray()
    type(amr_tree_state_3d) :: heat
    integer :: n(3)
    integer(int64) :: counts(2)
    integer :: attempt,final_owner,attempts,final_count
    logical :: breakup_test,children_test,dispersed,etab_test
    real(dp) :: coefficient,gas_before(nv),gas_after(nv),particle_change(nv),recursive_accel(3),accepted_mass
    if(rank==nranks-1) drops(1)%position=[.0012_dp,.0025_dp,.0015_dp]
    recursive_accel=0
    if(test_mode=='sdc-recursive-force'.or.test_mode=='sdc-recursive-fixed') then
      options%mixture_time_scheme='euler';recursive_accel=[10000._dp,-20000._dp,5000._dp]
    end if
    if(test_mode=='sdc-recursive-fixed') options%fixed_parts=.true.
    final_owner=nranks-1
    if(test_mode=='sdc-recursive-crossing') then
      if(rank==nranks-1) drops(1)%position(1)=.00199999_dp
      final_owner=0
    end if
    if(test_mode=='sdc-recursive-periodic') then
      if(rank==nranks-1) drops(1)%position(1)=.00399999_dp
    end if
    if(test_mode=='sdc-recursive-entering') then
      if(rank==nranks-1) then
        drops(1)%position(1)=.00200001_dp;drops(1)%velocity(1)=-20._dp
      end if
    end if
    children_test=test_mode=='sdc-recursive-children'
    etab_test=test_mode=='sdc-recursive-etab'.or.test_mode=='sdc-recursive-etab-film'
    breakup_test=children_test.or.test_mode=='sdc-recursive-breakup'.or.test_mode=='sdc-recursive-khrt'.or. &
      etab_test
    attempts=2
    if(test_mode=='sdc-recursive-force'.or.test_mode=='sdc-recursive-fixed') attempts=3
    if(breakup_test) then
      attempts=3;options%breakup_model='etab';options%surface_tension=.072_dp
      options%liquid_viscosity=.001_dp;options%maximum_multiplicity=50
      if(test_mode=='sdc-recursive-khrt') options%breakup_model='khrt'
      if(etab_test) options%etab_parcel_factor=.5_dp
      if(rank==final_owner) drops(1)%deformation=1.1_dp
    end if
    if(children_test) then
      options%breakup_model='khrt';options%khrt_children=.true.;options%maximum_multiplicity=0
      options%khrt_reference_diameter_cubed=1e-25_dp
      options%drag=.false.;options%heat=.false.;options%evaporation=.false.
      if(rank==final_owner) drops(1)%velocity=[1000._dp,10._dp,5._dp]
    end if
    if(test_mode=='sdc-recursive-film'.or.test_mode=='sdc-recursive-etab-film') then
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
    allocate(heat%fields(size(owners)))
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(heat%fields(node)%values(1,n(1),n(2),n(3)),source=53._dp)
    end do
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
      accel=recursive_accel
      if((test_mode=='sdc-recursive-force'.or.test_mode=='sdc-recursive-fixed').and.attempt==2.and.rank==0) accel(1)=accel(1)+1
      coefficient=merge(-1._dp,.1_dp,attempt==1)
      boundary=17;ledger=19;chemical=23;forcing=29;theta=31;exchange=37;counts=41
      if(allocated(film_settings)) then
      call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,tree,owners,g,extended,phi, &
        species,reactions,transport,liquids,options,field,temp,drops,masses,zero,spacing, &
        1e-8_dp,2,.true.,[.true.,.true.,.true.,.true.],zero,2e-7_dp,1e-12_dp,labels, &
        boundary,ledger,chemical,forcing,theta,steps,last_id,ok,difmag=coefficient, &
        gas_particle_change=exchange,level_steps=counts,particle_acceleration=accel, &
        redistribution=recursive_redistribution, &
        inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=wall_temperature,split_created=created,heat_release=heat, &
        films=films,film_mass=film_mass,film_settings=film_settings,film_steps=film_steps)
      else
      call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,tree,owners,g,extended,phi, &
        species,reactions,transport,liquids,options,field,temp,drops,masses,zero,spacing, &
        1e-8_dp,2,.true.,[.true.,.true.,.true.,.true.],zero,2e-7_dp,1e-12_dp,labels, &
        boundary,ledger,chemical,forcing,theta,steps,last_id,ok,difmag=coefficient, &
        gas_particle_change=exchange,level_steps=counts,particle_acceleration=accel, &
        redistribution=recursive_redistribution, &
        inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=wall_temperature,split_created=created,heat_release=heat)
      end if
      if(attempt<attempts) then
        call require(.not.ok.and.steps==0.and.last_id==1,'recursive early or late rejection retains ID')
        call require(all(boundary==17).and.all(ledger==19).and.all(chemical==23).and.all(forcing==29).and. &
          all(exchange==37).and.created==0.and.theta==31.and.all(counts==41),'recursive whole-step ledger rollback')
        call require(size(drops)==size(saved_drops).and.all(masses==saved_mass),'recursive population rollback')
        do p=1,size(drops)
          call require(all(drops(p)%position==saved_drops(p)%position).and. &
            all(drops(p)%velocity==saved_drops(p)%velocity),'recursive particle motion rollback')
        end do
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
          call require(all(heat%fields(node)%values==53._dp),'recursive late failure retains stage heat')
        end do
      else
        if(test_mode=='sdc-recursive-state-redist') print *, 'StateRedist recursive result: ',rank,ok,steps,counts
        call require(ok.and.steps==6.and.all(counts==[1_int64,2_int64]),'coupled recursive gas particle intervals')
        do node=1,size(owners)
          if(owners(node)/=rank) cycle
          call require(all(ieee_is_finite(heat%fields(node)%values)),'recursive spray finite stage heat')
          call require(any(heat%fields(node)%values/=53._dp),'recursive spray publishes stage heat')
        end do
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
      if(etab_test) then
        final_count=size(drops)
        call MPI_Allreduce(MPI_IN_PLACE,final_count,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call require(final_count>1.and.last_id>int(final_count,int64), &
          'recursive ETAB replaces retired parents')
        if(rank==final_owner) then
          dispersed=.false.
          do p=2,size(drops)
            dispersed=dispersed.or.any(abs(drops(p)%velocity-drops(1)%velocity)>1e-6_dp)
          enddo
          call require(dispersed,'recursive ETAB transverse dispersion')
        endif
      endif
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
    call require((children_test.or.etab_test.or.last_id==int(final_count,int64)).and. &
      size(drops)==merge(final_count,0,rank==final_owner), &
      'recursive physical collection')
    if(rank==final_owner.and.test_mode=='sdc-recursive-crossing') &
      call require(drops(1)%position(1)>.002_dp,'crossed fine boundary')
    if(rank==final_owner.and.test_mode=='sdc-recursive-periodic') &
      call require(drops(1)%position(1)>=0.and.drops(1)%position(1)<2e-7_dp,'recursive wrap enters fine patch')
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
    call require((children_test.or.etab_test.or.last_id==int(final_count,int64)).and. &
      size(drops)==merge(final_count,0,rank==final_owner), &
      'recursive physical collection')
    call require(any(exchange/=0),'recursive particle feedback reaches gas')
    if(test_mode=='sdc-recursive-force') then
      call require(all(forcing(6:)==0).and.forcing(1)==0,'particle force changes momentum and energy only')
      call require(all(forcing(2:4)*recursive_accel>0).and.forcing(5)>0,'recursive physical force ledger')
      accepted_mass=sum(drops%mass*drops%multiplicity)
      call MPI_Allreduce(MPI_IN_PLACE,accepted_mass,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call require(ierr==MPI_SUCCESS,'recursive accepted particle mass')
      call require(all(forcing(2:4)/(1e-8_dp*recursive_accel)>=accepted_mass-1e-22_dp).and. &
        all(forcing(2:4)/(1e-8_dp*recursive_accel)<=1e-10_dp+1e-22_dp), &
        'recursive force impulse bounded by physical initial and accepted mass')
    end if
    if(test_mode=='sdc-recursive-fixed') then
      coefficient=-.5_dp*1e-10_dp*sum([20._dp,-20._dp,5._dp]**2)
      call require(abs(forcing(5)/coefficient-1)<2e-13_dp,'fixed recursive constraint removes initial kinetic energy')
      if(rank==final_owner) then
        call require(all(drops(1)%velocity==0).and.all(drops(1)%position==saved_drops(1)%position), &
          'fixed recursive physical position and velocity')
      end if
    end if
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
    if(test_mode=='sdc-force-euler') options%mixture_time_scheme='euler'
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
