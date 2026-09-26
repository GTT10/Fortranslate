program test_spray_mpi_eb_deposit
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use spray_mpi_amr_field_mod, only: deposit_mpi_amr_tree_spray_gas
  implicit none
  type(amr_patch_3d) :: patch(1)
  type(amr_patch_tree_3d) :: tree
  type(eb_geometry_3d) :: finest,g(2)
  type(amr_tree_state_3d) :: field,temp,saved,saved_t,source,saved_source
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: positions(:,:),feedback(:,:),cell(:),primitive(:),y(:),total(:),parent(:),fine(:)
  real(dp) :: phi(0:8,0:8,0:8),zero(3)=0,spacing(3)=1,t,sound,r,volume
  integer :: rank,nranks,ierr,owners(2),dims(3),node,i,j,k,nv,nitrogen
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
  if(rank==0) then
    positions(:,1)=[1.9_dp,1.3_dp,1.5_dp];feedback(:,1)=.01_dp*cell
    if(any(periodic)) positions(:,1)=[3.9_dp,1.3_dp,3.9_dp]
  end if
  call deposit()
  call require(ok,'EB gas deposit')
  total=0;parent=0;fine=0
  do node=1,2
    if(owners(node)/=rank) cycle
    dims=tree%nodes(node)%shape
    do k=1,dims(3)
    do j=1,dims(2)
    do i=1,dims(1)
      volume=g(node)%volume_fraction(i,j,k)*g(node)%dx*g(node)%dy*g(node)%dz
      if(volume==0) then
        call require(all(field%fields(node)%values(:,i,j,k)==0),'solid state untouched')
        call require(temp%fields(node)%values(1,i,j,k)==42,'solid EOS skipped')
      end if
      if(node==1.and.i<=2) then
        parent=parent+volume*field%fields(node)%values(:,i,j,k)
      else
        total=total+volume*field%fields(node)%values(:,i,j,k)
      end if
      if(node==2) fine=fine+volume*field%fields(node)%values(:,i,j,k)
    end do
    end do
    end do
  end do
  call MPI_Allreduce(MPI_IN_PLACE,total,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(MPI_IN_PLACE,parent,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(MPI_IN_PLACE,fine,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(all(abs(total-44.01_dp*cell)<1e-13_dp*max(1._dp,abs(total))),'composite inventory')
  call require(all(abs(parent-fine)<1e-13_dp*max(1._dp,abs(fine))),'fluid weighted parent restriction')
  saved=field;saved_t=temp;saved_source=source
  if(rank==0) feedback(:,1)=-1e8_dp*cell
  call deposit()
  call require(.not.ok,'nonphysical final gas rejection')
  do node=1,2
    if(owners(node)/=rank) cycle
    call require(all(source%fields(node)%values==saved_source%fields(node)%values),'source EOS rollback')
    call require(all(field%fields(node)%values==saved%fields(node)%values),'gas rollback')
    call require(all(temp%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
  end do
  ! Below gas-state ulp: the source must still retain the transferred inventory.
  if(rank==0) feedback(:,1)=1.e-25_dp*cell
  call deposit()
  call require(ok,'tiny conservative source deposition')
  total=0
  do node=1,2
    call require(allocated(source%fields(node)%values).eqv.(owners(node)==rank),'owner-local source field')
    if(owners(node)/=rank) cycle
    dims=tree%nodes(node)%shape
    do k=1,dims(3)
    do j=1,dims(2)
    do i=1,dims(1)
      if(node==1.and.i<=2) cycle
      volume=g(node)%volume_fraction(i,j,k)*g(node)%dx*g(node)%dy*g(node)%dz
      total=total+volume*source%fields(node)%values(:,i,j,k)
    end do
    end do
    end do
    call require(all(field%fields(node)%values==saved%fields(node)%values),'tiny feedback below gas-state ulp')
  end do
  call MPI_Allreduce(MPI_IN_PLACE,total,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(all(abs(total-1.e-25_dp*cell)<3.e-14_dp*max(abs(1.e-25_dp*cell),tiny(1._dp))), &
    'source retains sub-ulp particle inventory')
  if(rank==0) print *, 'tiny source relative error: ',maxval(abs(total-1.e-25_dp*cell)/max(abs(1.e-25_dp*cell),tiny(1._dp)))
  if(rank==0) print *, 'PASS: conservative MPI EB feedback, fluid restriction and EOS rollback'
  call MPI_Finalize(ierr)
contains
  subroutine deposit()
    call deposit_mpi_amr_tree_spray_gas(MPI_COMM_WORLD,species,tree,owners,field,temp,positions,feedback, &
      zero,spacing,periodic,ok,eb_geometry=g,nodal_phi=phi,gas_change=source)
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
