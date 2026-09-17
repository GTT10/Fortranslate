program test_spray_mpi_amr_gas
  use mpi_f08
  use spray_mpi_amr_field_mod, only: deposit_mpi_amr_tree_spray_gas
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use spray_amr_tree_coupling_mod, only: deposit_amr_tree_spray_gas
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: patches(2)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: gas,temp,saved,saved_t,expected,expected_t,local,local_t
  real(dp), allocatable :: root(:,:,:,:),primitive(:),cell(:),y(:),feedback(:,:),before(:),after(:),scale(:)
  real(dp) :: positions(3,3),spacing(3)=.001_dp,zero(3)=0._dp,t,sound,r
  integer :: nv,node,n(3),i,j,k,p,water,nitrogen,mode
  real(dp), allocatable :: local_positions(:,:),local_feedback(:,:),ordered_positions(:,:),ordered_feedback(:,:)
  integer :: rank,nranks,ierr,owners(3),count,q,offset,source
  logical :: ok,periodic(3)=[.true.,.true.,.true.]
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'ranks')
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermodynamics')
  nv=reactive_nvar(size(species))
  allocate(root(nv,4,4,4),primitive(nv),cell(nv),y(size(species)),feedback(nv,3),before(nv),after(nv),scale(nv))
  water=0;nitrogen=0
  do p=1,size(species)
    if (species(p)%name=='H2O') water=p
    if (species(p)%name=='N2') nitrogen=p
  end do
  call require(water>0.and.nitrogen>0,'species identity')
  y=0._dp;y(nitrogen)=1._dp
  r=mixture_specific_gas_constant(species,y,ok)
  call require(ok,'gas constant')
  primitive=0._dp;primitive(1)=1._dp;primitive(5)=r*1000._dp;primitive(6:)=y
  call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
  call require(ok,'initial gas')
  do k=1,4
  do j=1,4
  do i=1,4
    root(:,i,j,k)=cell
  end do
  end do
  end do
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),ok)
  call require(ok,'fine patch')
  call initialize_amr_patch_3d(4,8,8,1,2,1,2,1,2,2,patches(2),ok)
  call require(ok,'nested patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
  call prolong_amr_patch_tree_3d(tree,root,gas,ok)
  call require(ok,'gas fields')
  allocate(temp%fields(3))
  do node=1,3
    n=tree%nodes(node)%shape
    allocate(temp%fields(node)%values(1,n(1),n(2),n(3)))
    temp%fields(node)%values=1000._dp
  end do
  y=0._dp;y(water)=1._dp;r=mixture_specific_gas_constant(species,y,ok)
  call require(ok,'vapor gas constant')
  primitive(1)=1._dp;primitive(2:4)=[10._dp,-20._dp,30._dp];primitive(5)=r*1500._dp;primitive(6:)=y
  call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
  call require(ok,'vapor feedback state')
  positions(:,1)=spacing*[1.9_dp,1.1_dp,1.1_dp]
  positions(:,2)=spacing*[.9_dp,.9_dp,.9_dp]
  positions(:,3)=spacing*[3.9_dp,1.1_dp,1.1_dp]
  do p=1,3
    feedback(:,p)=real(p,dp)*1.e-11_dp*cell
  end do
  owners=[0,1,nranks-1]
  count=0
  do p=1,3
    if (mod(p,nranks)==rank) count=count+1
  end do
  allocate(local_positions(3,count),local_feedback(nv,count),ordered_positions(3,3),ordered_feedback(nv,3))
  q=0
  do p=1,3
    if (mod(p,nranks)/=rank) cycle
    q=q+1;local_positions(:,q)=positions(:,p);local_feedback(:,q)=feedback(:,p)
  end do
  offset=0
  do source=0,nranks-1
    do p=1,3
      if (mod(p,nranks)/=source) cycle
      offset=offset+1;ordered_positions(:,offset)=positions(:,p);ordered_feedback(:,offset)=feedback(:,p)
    end do
  end do
  expected=gas;expected_t=temp
  call deposit_amr_tree_spray_gas(species,tree,expected,expected_t,ordered_positions,ordered_feedback, &
    zero,spacing,periodic,ok)
  call require(ok,'serial feedback reference')
  call composite_integrals_amr_patch_tree_3d(tree,gas,spacing,before,ok)
  call require(ok,'initial budget')
  call composite_integrals_amr_patch_tree_3d(tree,expected,spacing,after,ok)
  scale=max(abs(before),sum(abs(feedback),dim=2),tiny(1._dp))
  call require(ok.and.all(abs(after-before-sum(feedback,dim=2))<=3.e-13_dp*scale),'reference feedback budget')
  call distribute()
  call deposit()
  call require(ok,'distributed thermodynamic feedback');call compare_expected()
  saved=local;saved_t=local_t
  do mode=1,3
    local_feedback=0._dp
    if (rank==nranks-1) then
      select case(mode)
      case(1)
        local_feedback(1,1)=-1._dp
      case(2)
        local_feedback(5,1)=-1._dp
      case(3)
        local_feedback(5+nitrogen,1)=-1._dp
      end select
    end if
    call deposit()
    call require(.not.ok,'remote inadmissible feedback rejected');call unchanged()
  end do
  local_feedback=0._dp
  call deposit_mpi_amr_tree_spray_gas(MPI_COMM_WORLD,species,tree,owners,local,local_t, &
    local_positions,local_feedback,zero,spacing,periodic,ok,maximum_received=0)
  call require(.not.ok,'capacity rejection');call unchanged()
  if (rank==owners(3)) local%fields(3)%values(1,1,1,1)=-1._dp
  saved=local
  call deposit()
  call require(.not.ok,'invalid initial owner state rejected');call unchanged()
  call distribute();saved=local;saved_t=local_t
  if (rank==owners(3)) deallocate(local_t%fields(3)%values)
  call deposit()
  call require(.not.ok,'missing owner temperature rejected')
  do node=1,3
    if (owners(node)/=rank) cycle
    call require(all(local%fields(node)%values==saved%fields(node)%values),'missing-temperature conserved rollback')
  end do
  owners=modulo(owners+1,nranks);call distribute()
  q=0
  do p=1,3
    if (mod(p,nranks)/=rank) cycle
    q=q+1;local_feedback(:,q)=feedback(:,p)
  end do
  call deposit()
  call require(ok,'changed gas ownership');call compare_expected()
  saved=local;saved_t=local_t
  call deposit_mpi_amr_tree_spray_gas(MPI_COMM_WORLD,species,tree,owners,local,local_t, &
    local_positions(:,:0),local_feedback(:,:0),zero,spacing,periodic,ok)
  call require(ok,'globally empty feedback');call unchanged()
  if (rank==0) print *, 'PASS: transactional MPI AMR spray gas feedback'
  call MPI_Finalize(ierr)
contains
  subroutine deposit()
    call deposit_mpi_amr_tree_spray_gas(MPI_COMM_WORLD,species,tree,owners,local,local_t, &
      local_positions,local_feedback,zero,spacing,periodic,ok)
  end subroutine
  subroutine distribute()
    integer :: q
    if (allocated(local%fields)) deallocate(local%fields)
    if (allocated(local_t%fields)) deallocate(local_t%fields)
    allocate(local%fields(3),local_t%fields(3))
    do q=1,3
      if (owners(q)/=rank) cycle
      local%fields(q)%values=gas%fields(q)%values
      local_t%fields(q)%values=temp%fields(q)%values
    end do
  end subroutine
  subroutine compare_expected()
    integer :: q
    do q=1,3
      if (owners(q)==rank) then
        call require(all(local%fields(q)%values==expected%fields(q)%values),'exact serial conserved fields')
        call require(all(local_t%fields(q)%values==expected_t%fields(q)%values),'exact serial EOS temperatures')
      else
        call require(.not.allocated(local%fields(q)%values).and..not.allocated(local_t%fields(q)%values), &
          'nonowner U/T storage absent')
      end if
    end do
  end subroutine
  subroutine unchanged()
    integer :: q
    do q=1,3
      if (owners(q)/=rank) cycle
      call require(all(local%fields(q)%values==saved%fields(q)%values),'conserved rollback')
      call require(all(local_t%fields(q)%values==saved_t%fields(q)%values),'temperature rollback')
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL rank ',rank,': ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine
end program
