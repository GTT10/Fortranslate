program test_spray_mpi_eb_step
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
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  type(spray_parcel), allocatable :: drops(:),saved_drops(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:),before(:),after(:),ledger(:)
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
  character(len=16) :: mode
  character(len=32), allocatable :: faces(:)
  real(dp), allocatable :: exterior(:)
  integer :: initial_rank,final_rank
  call get_command_argument(1,mode)
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'rank count')
  owners=[0,nranks-1]
  initial_rank=nranks-1;final_rank=0
  if(mode=='periodic') then
    initial_rank=0;final_rank=nranks-1
    allocate(faces(6));faces='periodic';faces(3:4)='outflow'
  end if
  if(mode=='wall') then
    initial_rank=0;final_rank=0
    allocate(faces(6));faces='outflow';faces(1:2)='no_slip_wall'
  end if
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
  allocate(drops(merge(1,0,rank==initial_rank)),masses(2,merge(1,0,rank==initial_rank)),before(nv),after(nv),ledger(nv))
  if(rank==initial_rank) then
    drops(1)=spray_parcel(1_int64,[.00199_dp,.001250001_dp,.0015_dp], &
      [20._dp,-20._dp,5._dp],1e-12_dp,350._dp,100._dp)
    if(allocated(faces)) drops(1)%position(1)=.00399_dp
    masses(:,1)=drops(1)%mass*[.25_dp,.75_dp]
  end if
  if(allocated(faces)) then
    allocate(exterior(nv));exterior=23
  end if
  saved=field;saved_t=temp;saved_drops=drops;saved_mass=masses
  ledger=42;options%maximum_substeps=1
  call advance()
  call require(.not.ok.and.steps==0.and.all(ledger==42),'whole interval failure and ledger rollback')
  call require(all(masses==saved_mass).and.size(drops)==size(saved_drops),'particle population rollback')
  if(rank==initial_rank) call require(all(drops(1)%position==saved_drops(1)%position),'particle position rollback')
  do node=1,2
    if(owners(node)/=rank) cycle
    call require(all(field%fields(node)%values==saved%fields(node)%values),'gas rollback')
    call require(all(temp%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
  end do
  options%maximum_substeps=100000
  call inventory(before)
  call advance()
  call require(ok.and.steps>1,'distributed EB interval')
  call inventory(after)
  if(mode=='wall') ledger=ledger+exterior
  call require(all(abs(after-before-ledger)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'composite gas parcel wall conservation')
  call require(ledger(3)>0,'reflected wall impulse')
  if(mode=='periodic') call require(all(exterior==0),'periodic crossing has no exterior loss')
  if(mode=='wall') then
    call require(exterior(2)<0,'outer wall impulse')
    call require(exterior(1)==0.and.all(exterior(5:)==0),'outer wall has no mass or energy loss')
  end if
  call require(size(drops)==merge(1,0,rank==final_rank),'final gas patch owner')
  if(rank==final_rank) then
    if(mode=='wall') then
      call require(drops(1)%position(1)<.004_dp.and.drops(1)%velocity(1)<0,'outer wall specular reflection')
    else if(allocated(faces)) then
      call require(drops(1)%position(1)<.002_dp.and.drops(1)%position(1)>=0,'periodic crossing to refined owner')
    else
      call require(drops(1)%position(1)>.002_dp,'coarse crossing')
    end if
    call require(drops(1)%position(2)>.00125_dp,'fluid position')
    call require(drops(1)%temperature>350.and.drops(1)%mass<1e-12_dp,'heating and evaporation')
    call require(drops(1)%velocity(2)>0.and.drops(1)%velocity(1)<20,'reflection and drag')
  end if
  call check_kick_stages()
  if(rank==0) print *, 'PASS: MPI EB mixed-particle interval and global wall inventory'
  call MPI_Finalize(ierr)
contains
  subroutine check_kick_stages()
    type(amr_tree_state_3d) :: reference_gas,reference_t,reference_source,stage_gas,stage_t,stage_source
    type(spray_parcel), allocatable :: initial(:),reference(:),candidate(:)
    real(dp), allocatable :: initial_mass(:,:),reference_mass(:,:),candidate_mass(:,:)
    character(len=16) :: stage_mode
    real(dp) :: interval,wall(nv),wall_reference(nv),expected_position(3)
    integer :: mode_index,q,substeps
    logical :: good
    initial=saved_drops;initial_mass=saved_mass;interval=2.e-10_dp
    do q=1,size(initial)
      initial(q)%position=[.0015_dp,.002_dp,.0015_dp]
    end do
    reference_gas=saved;reference_t=saved_t;reference=initial;reference_mass=initial_mass
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options, &
      reference_gas,reference_t,tree,owners,reference,reference_mass,zero,spacing,.5_dp*interval,substeps,good, &
      boundary_change=wall_reference,eb_geometry=g,nodal_phi=phi,gas_change=reference_source)
    call require(good.and.substeps==1,'reference half kick')
    do mode_index=1,2
      stage_mode='kick'
      if(mode_index==2) stage_mode='kick_drift'
      stage_gas=saved;stage_t=saved_t;candidate=initial;candidate_mass=initial_mass
      call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options, &
        stage_gas,stage_t,tree,owners,candidate,candidate_mass,zero,spacing,interval,substeps,good, &
        boundary_change=wall,eb_geometry=g,nodal_phi=phi,gas_change=stage_source,motion_mode=stage_mode)
      call require(good.and.substeps==1,'kick stage interval')
      call require(all(candidate_mass==reference_mass),'half-interval component transfer')
      call require(size(candidate)==size(reference),'stage parcel ownership')
      call require(all(wall==wall_reference),'stage wall ledger')
      do q=1,size(candidate)
        call require(all(candidate(q)%velocity==reference(q)%velocity),'half-interval velocity kick')
        call require(candidate(q)%mass==reference(q)%mass.and.candidate(q)%temperature==reference(q)%temperature, &
          'half-interval heat and mass update')
        expected_position=[.0015_dp,.002_dp,.0015_dp]
        if(mode_index==2) expected_position=expected_position+interval*candidate(q)%velocity
        call require(all(candidate(q)%position==expected_position),'full drift or stationary kick')
      end do
      do q=1,size(owners)
        if(owners(q)/=rank) cycle
        call require(all(stage_gas%fields(q)%values==reference_gas%fields(q)%values),'half-interval gas feedback')
        call require(all(stage_source%fields(q)%values==reference_source%fields(q)%values),'half-interval source export')
      end do
    end do
    ! A kick after a completed drift must not move the particles again.
    initial=candidate;initial_mass=candidate_mass
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options, &
      stage_gas,stage_t,tree,owners,candidate,candidate_mass,zero,spacing,interval,substeps,good, &
      boundary_change=wall,eb_geometry=g,nodal_phi=phi,gas_change=stage_source,motion_mode='kick')
    call require(good,'second kick')
    call require(size(candidate)==size(initial),'second kick ownership')
    do q=1,size(candidate)
      call require(all(candidate(q)%position==initial(q)%position),'second kick leaves drift endpoint unchanged')
    end do
    initial=candidate;initial_mass=candidate_mass;reference_gas=stage_gas;reference_source=stage_source
    stage_mode='kick'
    if(rank==0) stage_mode='kick_drift'
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options, &
      stage_gas,stage_t,tree,owners,candidate,candidate_mass,zero,spacing,interval,substeps,good, &
      boundary_change=wall,eb_geometry=g,nodal_phi=phi,gas_change=stage_source,motion_mode=stage_mode)
    call require(.not.good,'rank-dependent kick mode rejected')
    call require(all(candidate_mass==initial_mass),'mode disagreement component rollback')
    do q=1,size(candidate)
      call require(all(candidate(q)%position==initial(q)%position),'mode disagreement position rollback')
    end do
    do q=1,size(owners)
      if(owners(q)/=rank) cycle
      call require(all(stage_gas%fields(q)%values==reference_gas%fields(q)%values),'mode disagreement gas rollback')
      call require(all(stage_source%fields(q)%values==reference_source%fields(q)%values),'mode disagreement source rollback')
    end do
    field=saved;temp=saved_t;drops=saved_drops;masses=saved_mass
    call inventory(before)
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,field,temp, &
      tree,owners,drops,masses,zero,spacing,2.e-5_dp,substeps,good, &
      boundary_change=wall,eb_geometry=g,nodal_phi=phi,gas_change=stage_source,motion_mode='kick_drift')
    call require(good.and.substeps>1,'subcycled EB kick-drift')
    call inventory(after)
    call require(all(abs(after-before-wall)<=3.e-12_dp*max(abs(before),abs(after),1.e-30_dp)), &
      'subcycled kick-drift gas liquid wall conservation')
    if(rank==0) print *, 'kick-drift relative inventory residual: ', &
      maxval(abs(after-before-wall)/max(abs(before),abs(after),1.e-30_dp))
  end subroutine check_kick_stages
  subroutine advance()
    type(amr_tree_state_3d) :: previous,source
    integer :: a
    previous=field;source=field
    do a=1,size(owners)
      if(owners(a)==rank) source%fields(a)%values=17
    end do
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,field,temp, &
      tree,owners,drops,masses,zero,spacing,2e-5_dp,steps,ok,boundary_change=ledger,eb_geometry=g,nodal_phi=phi, &
      boundary_faces=faces,eb_outflow_change=exterior,gas_change=source)
    do a=1,size(owners)
      call require(allocated(source%fields(a)%values).eqv.(owners(a)==rank),'owner-local particle source')
      if(owners(a)/=rank) cycle
      if(ok) then
        call require(all(abs(field%fields(a)%values-previous%fields(a)%values-source%fields(a)%values)<= &
          3.e-13_dp*max(1._dp,abs(field%fields(a)%values))),'source reconstructs accepted particle feedback')
      else
        call require(all(source%fields(a)%values==17),'particle source rollback')
      end if
    end do
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
