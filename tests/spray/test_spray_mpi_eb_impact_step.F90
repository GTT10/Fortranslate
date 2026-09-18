program test_spray_mpi_eb_impact_step
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use spray_mpi_amr_step_mod, only: advance_mpi_amr_mixture_coupling
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  type(spray_wall_options) :: walls
  type(spray_wall_film), allocatable :: films(:)
  real(dp), allocatable :: film_mass(:,:)
  integer(int64) :: last_id
  integer :: scenario,local_count,global_count
  real(dp) :: interval
  logical :: rejected
  type(spray_parcel), allocatable :: drops(:),saved_drops(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:),before(:),after(:),ledger(:),exterior(:)
  integer :: steps,p
  real(dp) :: rho,cp,e,amount
  type(amr_patch_3d) :: patch(1)
  type(amr_patch_tree_3d) :: tree
  type(eb_geometry_3d) :: finest,g(2)
  type(amr_tree_state_3d) :: field,temp,saved,saved_t
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: positions(:,:),feedback(:,:),cell(:),primitive(:),y(:),total(:),parent(:),fine(:)
  real(dp) :: phi(0:8,0:8,0:8),zero(3)=0,spacing(3)=.001_dp,t,sound,r,volume
  integer :: rank,nranks,ierr,owners(2),dims(3),node,i,j,k,nv,nitrogen
  logical :: ok
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
  saved=field;saved_t=temp;saved_drops=drops;saved_mass=masses
  allocate(films(0),film_mass(2,0))
  options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
  walls%enabled=.true.;walls%temperature=300
  do scenario=1,2
  field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass
  deallocate(films,film_mass);allocate(films(0),film_mass(2,0))
  last_id=1;interval=2e-5_dp
  if(scenario==2) then
    options%breakup_model='khrt';options%khrt_children=.true.
    options%khrt_reference_diameter_cubed=1e-25_dp
    interval=1e-8_dp
    if(rank==nranks-1) drops(1)%velocity=[20._dp,-200._dp,5._dp]
  end if
  if(scenario==2) then
    ! A preexisting thick wet film produces inadmissible empirical child fractions.
    ! The collective diagnostic must escape the timestep-halving loop unchanged.
    deallocate(films,film_mass)
    allocate(films(merge(1,0,rank==nranks-1)),film_mass(2,merge(1,0,rank==nranks-1)))
    if(rank==nranks-1) then
      films(1)=spray_wall_film([.00199_dp,.00125_dp,.0015_dp],[0._dp,1._dp,0._dp], &
        1e-10_dp,350._dp,.001_dp,2e-6_dp,1)
      film_mass(:,1)=films(1)%mass*[.25_dp,.75_dp]
    end if
    options%maximum_substeps=100000;ledger=42
    call advance()
    call require(.not.ok.and.rejected.and.steps==0.and.all(ledger==42),'collective model rejection')
    call require(last_id==1.and.all(masses==saved_mass),'model rejection particle and ID rollback')
    if(rank==nranks-1) then
      call require(all(drops(1)%position==saved_drops(1)%position),'model rejection position rollback')
      call require(films(1)%mass==1e-10_dp.and.films(1)%height==2e-6_dp,'model rejection film rollback')
    end if
    do node=1,2
      if(owners(node)/=rank) cycle
      call require(all(field%fields(node)%values==saved%fields(node)%values),'model rejection gas rollback')
    end do
    deallocate(films,film_mass);allocate(films(0),film_mass(2,0))
  end if
  ledger=42;options%maximum_substeps=1
  if(scenario==2) options%maximum_parcels=1
  call advance()
  call require(.not.ok.and.steps==0.and.all(ledger==42),'whole interval failure and ledger rollback')
  call require(all(masses==saved_mass).and.size(drops)==size(saved_drops),'particle population rollback')
  call require(size(films)==0.and.size(film_mass,2)==0.and.last_id==1,'film and watermark rollback')
  if(rank==nranks-1) call require(all(drops(1)%position==saved_drops(1)%position),'particle position rollback')
  do node=1,2
    if(owners(node)/=rank) cycle
    call require(all(field%fields(node)%values==saved%fields(node)%values),'gas rollback')
    call require(all(temp%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
  end do
  options%maximum_substeps=100000;options%maximum_parcels=100
  call inventory(before)
  call advance()
  call require(ok.and.steps>=1,'distributed EB interval')
  call inventory(after)
  call require(all(abs(after-before-ledger)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'composite gas parcel wall conservation')
  call require(ledger(3)>0,'wall impulse returned with system sign')
  local_count=size(films)
  call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(global_count==1,'one deposited film')
  call require(size(films)==merge(1,0,rank==nranks-1),'film fine ownership')
  if(size(films)>0) call require(films(1)%mass>0.and.films(1)%temperature>350,'heated liquid deposition')
  local_count=size(drops)
  call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(scenario==1) then
    call require(global_count==0.and.last_id==1,'deposition consumes parent without IDs')
  else
    call require(global_count==4.and.last_id==5,'four splash children and shared watermark')
    do p=1,size(drops)
      call require(drops(p)%id>=2.and.drops(p)%id<=5.and.drops(p)%parent_id==1,'child lineage')
      call require(drops(p)%position(2)>.00125_dp.and.drops(p)%velocity(2)>0,'remaining child flight')
      call require((drops(p)%position(1)<.002_dp).eqv.(rank==nranks-1),'final finest particle owner')
    end do
  end if
  end do
  ! No wall crossing: KH release conserves composite gas/particle inventory,
  ! advances children and migrates their full state to the final finest owners.
  field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass
  deallocate(films,film_mass);allocate(films(0),film_mass(2,0))
  if(rank==nranks-1) then
    drops(1)%position=[.001995_dp,.0015_dp,.0015_dp]
    drops(1)%velocity=[1000._dp,10._dp,5._dp]
  end if
  saved_drops=drops;interval=1e-8_dp;last_id=1;ledger=42
  options%drag=.false.;options%heat=.false.;options%evaporation=.false.
  options%maximum_parcels=1;options%maximum_substeps=1
  call advance()
  call require(.not.ok.and.steps==0.and.last_id==1.and.all(ledger==42),'KH capacity interval rollback')
  call require(size(drops)==size(saved_drops).and.all(masses==saved_mass),'KH payload rollback')
  if(rank==nranks-1) call require(all(drops(1)%position==saved_drops(1)%position),'KH position rollback')
  options%maximum_parcels=100;options%maximum_substeps=100000
  call inventory(before)
  call advance()
  call require(ok.and.last_id>1,'EB KH release')
  call inventory(after)
  call require(all(abs(after-before-ledger)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'KH transverse energy and composite conservation')
  local_count=size(drops)
  call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(global_count==5.and.last_id==5.and.size(films)==0,'one KH four-child release')
  do p=1,size(drops)
    call require(drops(p)%position(1)>=.002_dp.and.rank==0,'KH final coarse owner and flight')
    if(drops(p)%id>1) call require(drops(p)%parent_id==1,'KH child lineage')
    call require(all(abs(masses(:,p)/drops(p)%mass-[.25_dp,.75_dp])<1e-13_dp),'KH component preservation')
  end do
  field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass;last_id=1;ledger=0
  call inventory(before)
  call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,field,temp, &
    tree,owners,drops,masses,zero,spacing,interval,steps,ok,last_id=last_id,boundary_change=ledger, &
    eb_geometry=g,nodal_phi=phi)
  call require(ok.and.last_id==5,'KH release with specular EB flight')
  call inventory(after)
  call require(all(abs(after-before-ledger)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'specular KH composite conservation')
  field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass
  deallocate(films,film_mass);allocate(films(0),film_mass(2,0),exterior(nv))
  if(rank==nranks-1) then
    drops(1)%position=[.00199_dp,.001250001_dp,.0015_dp]
    drops(1)%velocity=[20._dp,-200._dp,5._dp]
  end if
  saved_drops=drops;last_id=1;ledger=42;exterior=43;interval=1e-3_dp
  options%breakup_model='none';options%khrt_children=.false.;options%khrt_reference_diameter_cubed=0
  options%maximum_substeps=1
  call outflow_advance()
  call require(.not.ok.and.steps==0.and.last_id==1,'outflow interval rollback')
  call require(size(drops)==size(saved_drops).and.size(films)==0.and.all(masses==saved_mass),'outflow arrays rollback')
  call require(all(ledger==42).and.all(exterior==43),'outflow separate ledger rollback')
  if(rank==nranks-1) call require(all(drops(1)%position==saved_drops(1)%position),'outflow position rollback')
  options%maximum_substeps=100000
  call inventory(before)
  call outflow_advance()
  call require(ok.and.size(drops)==0.and.last_id==5,'source splash children escape')
  call require(exterior(1)<0,'source global exterior inventory')
  call require(abs(exterior(5+liquids(1)%vapor_index)-.25_dp*exterior(1))<1e-23_dp,'first exterior species')
  call require(abs(exterior(5+liquids(2)%vapor_index)-.75_dp*exterior(1))<1e-23_dp,'second exterior species')
  call inventory(after)
  call require(all(abs(after-before-ledger-exterior)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'MPI source gas film wall exterior conservation')
  field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass
  deallocate(films,film_mass);allocate(films(0),film_mass(2,0))
  last_id=1;ledger=42;exterior=43;options%maximum_substeps=1
  call specular_outflow_advance()
  call require(.not.ok.and.steps==0.and.last_id==1,'specular exterior interval rollback')
  call require(size(drops)==size(saved_drops).and.all(masses==saved_mass),'specular exterior payload rollback')
  call require(all(ledger==42).and.all(exterior==43),'specular exterior budgets rollback')
  options%maximum_substeps=100000
  call inventory(before)
  call specular_outflow_advance()
  call require(ok.and.size(drops)==0.and.last_id==1,'specular reflection then escape')
  call require(ledger(3)>0.and.exterior(1)<0,'specular wall and exterior signs')
  call require(abs(exterior(5+liquids(1)%vapor_index)-.25_dp*exterior(1))<1e-23_dp,'specular first species')
  call require(abs(exterior(5+liquids(2)%vapor_index)-.75_dp*exterior(1))<1e-23_dp,'specular second species')
  call inventory(after)
  call require(all(abs(after-before-ledger-exterior)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'MPI specular source exterior conservation')
  if(rank==0) print *, 'PASS: MPI EB impact source coupling and global gas parcel film wall inventory'
  call MPI_Finalize(ierr)
contains
  subroutine specular_outflow_advance()
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,field,temp, &
      tree,owners,drops,masses,zero,spacing,interval,steps,ok,last_id=last_id,boundary_change=ledger, &
      eb_geometry=g,nodal_phi=phi,eb_outflow_change=exterior)
  end subroutine
  subroutine outflow_advance()
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,field,temp, &
      tree,owners,drops,masses,zero,spacing,interval,steps,ok,last_id=last_id,boundary_change=ledger, &
      wall_films=films,film_mass=film_mass,film_settings=walls,eb_geometry=g,nodal_phi=phi, &
      model_rejected=rejected,eb_outflow_change=exterior)
  end subroutine
  subroutine advance()
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,field,temp, &
      tree,owners,drops,masses,zero,spacing,interval,steps,ok,last_id=last_id,boundary_change=ledger, &
      wall_films=films,film_mass=film_mass,film_settings=walls,eb_geometry=g,nodal_phi=phi,model_rejected=rejected)
  end subroutine
  subroutine inventory(result)
    real(dp), intent(out) :: result(:)
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
    do p=1,size(films)
      call liquid_mixture_properties(liquids,film_mass(:,p)/films(p)%mass,films(p)%temperature,rho,cp,e,ok)
      call require(ok,'film properties')
      result(1)=result(1)+films(p)%mass;result(5)=result(5)+films(p)%mass*e
      do i=1,2
        result(5+liquids(i)%vapor_index)=result(5+liquids(i)%vapor_index)+film_mass(i,p)
      end do
    end do
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
