program test_spray_mpi_eb_checkpoint
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
  use spray_mpi_eb_particle_checkpoint_mod, only: mpi_eb_particle_checkpoint
  implicit none
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  type(spray_parcel), allocatable :: drops(:),saved_drops(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:),before(:),after(:),ledger(:)
  character(len=1024) :: path,mode
  real(dp) :: time
  integer :: events
  integer(int64) :: mark
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
  if(rank==nranks-1) then
    drops(1)%id=3;drops(1)%parent_id=1;drops(1)%material_id=2
    drops(1)%deformation=.2_dp;drops(1)%deformation_rate=.3_dp;drops(1)%rt_time=.4_dp
    drops(1)%khrt_initial_count=100;drops(1)%khrt_shed_volume=1e-20_dp
  end if
  call get_command_argument(1,mode)
  call get_command_argument(2,path)
  time=1e-8_dp;steps=2;events=3;mark=7;ledger=5
  if(mode=='write') then
    call checkpoint(.true.,'test generation')
    call require(ok,'write')
    call checkpoint(.true.,'test generation')
    call require(.not.ok,'immutable archive')
  else
    deallocate(drops,masses);allocate(drops(0),masses(2,0))
    time=42;steps=43;events=44;mark=45;ledger=46
    call checkpoint(.false.,'wrong generation')
    call require(.not.ok.and.size(drops)==0.and.time==42,'context rollback')
    call checkpoint(.false.,'test generation')
    if(mode=='reject') then
      call require(.not.ok.and.size(drops)==0.and.time==42.and.steps==43.and.events==44, 'truncated rollback')
      call require(mark==45.and.all(ledger==46),'metadata rollback')
    else
      call require(ok,'changed-rank restore')
      call require(time==1e-8_dp.and.steps==2.and.events==3.and.mark==7.and.all(ledger==5),'metadata')
      call require(size(drops)==merge(1,0,rank==nranks-1),'new owner')
      if(rank==nranks-1) then
        call require(drops(1)%id==3.and.drops(1)%mass==1e-12_dp.and.drops(1)%temperature==350, 'parcel')
        call require(all(drops(1)%position==[.00199_dp,.001250001_dp,.0015_dp]),'position')
        call require(all(drops(1)%velocity==[20._dp,-20._dp,5._dp]).and.drops(1)%multiplicity==100,'velocity/count')
        call require(all(masses(:,1)==1e-12_dp*[.25_dp,.75_dp]),'components')
        call require(drops(1)%parent_id==1.and.drops(1)%material_id==2,'lineage/material')
        call require(drops(1)%deformation==.2_dp.and.drops(1)%deformation_rate==.3_dp,'deformation')
        call require(drops(1)%rt_time==.4_dp.and.drops(1)%khrt_initial_count==100.and. &
          drops(1)%khrt_shed_volume==1e-20_dp,'KHRT state')
      end if
    end if
  end if
  if(rank==0) print *, 'PASS: EB particle checkpoint ',trim(mode)
  call MPI_Finalize(ierr)
contains
  subroutine checkpoint(writing,context)
    logical, intent(in) :: writing
    character(len=*), intent(in) :: context
    call mpi_eb_particle_checkpoint(MPI_COMM_WORLD,writing,trim(path),context,tree,owners,phi,zero,spacing, &
      100,drops,masses,time,steps,events,mark,ledger,ok)
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
