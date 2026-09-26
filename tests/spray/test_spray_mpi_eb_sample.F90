program test_spray_mpi_eb_sample
  use mpi_f08
  use spray_mpi_amr_field_mod, only: sample_mpi_amr_tree_spray_field,deposit_mpi_amr_tree_spray_field
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use spray_eb_tree_stencil_mod, only: eb_tree_spray_stencil
  implicit none
  type(amr_tree_state_3d) :: field,change,initial_field
  real(dp), allocatable :: positions(:,:),samples(:,:),kept(:,:),contribution(:,:)
  integer :: rank,nranks,ierr,owners(2),dims(3)
  type(amr_patch_3d) :: patch(1)
  type(amr_patch_tree_3d) :: tree
  type(eb_geometry_3d) :: finest,g(2)
  real(dp) :: phi(0:8,0:8,0:8),position(3),zero(3)=0,spacing(3)=1,expected,total
  real(dp), allocatable :: w(:),v(:),saved_w(:),saved_v(:)
  integer, allocatable :: cells(:,:),saved(:,:)
  integer :: j,q,node
  logical :: ok,periodic(3)=.false.
  character(len=16) :: mode
  call get_command_argument(1,mode)
  if(mode=='periodic') periodic=[.true.,.false.,.true.]
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
    phi(:,j,:)=real(j,dp)*.5_dp-1.25_dp
  end do
  call build_level_set_eb_geometry_3d(phi,zero,[4._dp,4._dp,4._dp],finest,ok)
  call require(ok,'finest geometry')
  call coarsen_eb_geometry_3d(finest,2,g(1),ok)
  call require(ok,'root geometry')
  call extract_eb_geometry_3d(finest,[1,1,1],[4,8,8],g(2),ok)
  call require(ok,'child geometry')
  position=[1.9_dp,1.3_dp,1.5_dp]
  if(any(periodic)) position=[3.9_dp,1.3_dp,3.9_dp]
  call run()
  call require(ok.and.all(w>0).and.abs(sum(w)-1)<1e-14_dp,'normalized cut-interface weights')
  call require(any(cells(1,:)==1).and.any(cells(1,:)==2),'coarse and fine support')
  allocate(field%fields(2),positions(3,merge(1,0,rank==0)),samples(1,merge(1,0,rank==0)))
  do node=1,2
    if(owners(node)/=rank) cycle
    dims=tree%nodes(node)%shape
    allocate(field%fields(node)%values(1,dims(1),dims(2),dims(3)))
    field%fields(node)%values=real(10*node,dp)
  end do
  expected=sum(w*real(10*cells(1,:),dp))
  call require(expected>10.and.expected<20,'two-owner reference')
  if(rank==0) positions(:,1)=position
  samples=99;kept=samples
  call sample(0)
  call require(.not.ok.and.all(samples==kept),'request capacity rollback')
  if(rank==nranks-1) phi(0,0,0)=phi(0,0,0)+.01_dp
  call sample(100)
  call require(.not.ok.and.all(samples==kept),'rank geometry mismatch rollback')
  phi(0,0,0)=-1.25_dp
  if(rank==0) positions(2,1)=1._dp
  call sample(100)
  call require(.not.ok.and.all(samples==kept),'solid query rollback')
  if(rank==0) positions(:,1)=position
  call sample(100)
  call require(ok,'owner-routed EB sampling')
  if(rank==0) call require(abs(samples(1,1)-expected)<1e-13_dp,'weighted remote values')
  if(rank/=0) call require(size(samples,2)==0,'empty query ranks')
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,field,1,positions,zero,spacing, &
    periodic,samples,ok,100,g,phi,maximum_level=0)
  call require(ok.and.all(abs(samples-10._dp)<1e-13_dp),'covered coarse EB sampling')
  initial_field=field;allocate(contribution(1,size(positions,2)));contribution=2
  call deposit_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,field,1,positions,contribution, &
    zero,spacing,periodic,ok,eb_geometry=g,nodal_phi=phi,gas_change=change,maximum_level=0)
  call require(ok,'coarse EB source deposition')
  total=0
  if(rank==owners(1)) total=sum(change%fields(1)%values(1,:,:,:)*g(1)%volume_fraction)*g(1)%dx*g(1)%dy*g(1)%dz
  call MPI_Allreduce(MPI_IN_PLACE,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(abs(total-2._dp)<1e-13_dp,'coarse EB fluid-volume source inventory')
  if(rank==owners(2)) call require(all(field%fields(2)%values==initial_field%fields(2)%values).and. &
    all(change%fields(2)%values==0),'fine EB field unchanged by coarse source')
  deallocate(positions);allocate(positions(3,0))
  call sample(100)
  call require(ok.and.size(samples,2)==0,'globally empty queries')
  if(any(periodic)) then
    phi(0,8,8)=phi(0,8,8)+.1_dp
    call sample(100)
    call require(.not.ok.and.size(samples,2)==0,'empty periodic query still validates nodal endpoints')
  end if
  if(rank==0) print *, 'PASS: owner-local MPI EB composite particle sampling'
  call MPI_Finalize(ierr)
contains
  subroutine sample(limit)
    integer, intent(in) :: limit
    call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,field,1,positions,zero,spacing, &
      periodic,samples,ok,limit,g,phi)
  end subroutine
  subroutine run()
    call eb_tree_spray_stencil(tree,g,phi,position,zero,spacing,cells,w,v,ok,periodic=periodic)
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine
end program
