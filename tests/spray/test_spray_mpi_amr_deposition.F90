program test_spray_mpi_amr_deposition
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    composite_integrals_amr_patch_tree_3d
  use spray_amr_tree_coupling_mod, only: deposit_amr_tree_spray_field
  use spray_mpi_amr_field_mod, only: deposit_mpi_amr_tree_spray_field,sum_mpi_amr_particle_source_halo
  implicit none
  type(amr_patch_3d) :: patches(2)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: full,local,expected,saved,change
  real(dp), allocatable :: positions(:,:),contributions(:,:),all_positions(:,:),all_contributions(:,:)
  real(dp) :: before(2),local_integral(2),total_integral(2),volume
  integer :: r,p,count,offset,child,idx(3)
  logical :: covered
  real(dp) :: points(3,4),zero(3)=0._dp,spacing(3)=1._dp
  integer :: rank,nranks,ierr,node,n(3),i,j,k,np,owners(3),level
  logical :: ok,periodic(3)=[.true.,.true.,.true.]
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'rank count')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),ok)
  call require(ok,'fine patch')
  call initialize_amr_patch_3d(4,8,8,1,2,1,2,1,2,2,patches(2),ok)
  call require(ok,'nested patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
  allocate(full%fields(3))
  do node=1,3
    n=tree%nodes(node)%shape
    allocate(full%fields(node)%values(2,n(1),n(2),n(3)))
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      full%fields(node)%values(:,i,j,k)=[real(10*node+i+3*j+7*k,dp),real(node*i*j+k,dp)]
    end do
    end do
    end do
  end do
  owners=[0,1,nranks-1];call distribute()
  points(:,1)=[1.9_dp,1.1_dp,1.1_dp]
  points(:,2)=[2.1_dp,1.1_dp,1.1_dp]
  points(:,3)=[.9_dp,.9_dp,.9_dp]
  points(:,4)=[3.9_dp,1.1_dp,1.1_dp]
  np=2
  if (rank==0) np=0
  if (rank==1) np=4
  positions=points(:,:np)
  count=4+2*(nranks-2)
  allocate(contributions(2,np),all_positions(3,count),all_contributions(2,count))
  offset=0
  do r=1,nranks-1
    count=2
    if (r==1) count=4
    all_positions(:,offset+1:offset+count)=points(:,:count)
    do p=1,count
      all_contributions(:,offset+p)=[real((r+1)*p,dp),-.25_dp*(r+1)*p]
      if (rank==r) contributions(:,p)=all_contributions(:,offset+p)
    end do
    offset=offset+count
  end do
  do level=0,1
    call distribute()
    call deposit_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,contributions, &
      zero,spacing,periodic,ok,gas_change=change,maximum_level=level)
    call require(ok,'bounded source deposition')
    local_integral=0
    do node=1,3
      if(owners(node)/=rank) cycle
      call require(all(abs(local%fields(node)%values-full%fields(node)%values-change%fields(node)%values)<1e-12_dp), &
        'raw source increment reconstructs field')
      if(tree%nodes(node)%level>level) then
        call require(all(change%fields(node)%values==0).and. &
          all(local%fields(node)%values==full%fields(node)%values),'deeper field unchanged')
        cycle
      end if
      n=tree%nodes(node)%shape;volume=1._dp/8._dp**tree%nodes(node)%level
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        idx=[i,j,k];covered=.false.
        do child=node+1,3
          if(tree%nodes(child)%parent/=node.or.tree%nodes(child)%level>level) cycle
          associate(box=>tree%nodes(child)%patch)
          if(all(idx>=[box%coarse_i_lower,box%coarse_j_lower,box%coarse_k_lower]).and. &
            all(idx<=[box%coarse_i_upper,box%coarse_j_upper,box%coarse_k_upper])) covered=.true.
          end associate
        end do
        if(.not.covered) local_integral=local_integral+volume*change%fields(node)%values(:,i,j,k)
      end do
      end do
      end do
    end do
    call MPI_Allreduce(local_integral,total_integral,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(abs(total_integral-sum(all_contributions,dim=2))<1e-12_dp),'capped source inventory')
  end do
  saved=local
  call deposit_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,contributions, &
    zero,spacing,periodic,ok,maximum_received=0,maximum_level=0)
  call require(.not.ok,'bounded receive capacity rejection')
  do node=1,3
    if(owners(node)/=rank) cycle
    call require(all(local%fields(node)%values==saved%fields(node)%values),'bounded source rollback')
  end do
  call distribute()
  call composite_integrals_amr_patch_tree_3d(tree,full,spacing,before,ok)
  call require(ok,'initial composite integral')
  expected=full
  call deposit_amr_tree_spray_field(tree,expected,all_positions,all_contributions,zero,spacing,periodic,ok)
  call require(ok,'serial reference')
  call deposit()
  call require(ok,'owner deposition')
  call compare_expected()
  local_integral=0._dp
  do node=1,3
    if (owners(node)/=rank) cycle
    n=tree%nodes(node)%shape;volume=1._dp/8._dp**tree%nodes(node)%level
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      idx=[i,j,k];covered=.false.
      do child=node+1,3
        if (tree%nodes(child)%parent/=node) cycle
        associate(box=>tree%nodes(child)%patch)
        if (all(idx>=[box%coarse_i_lower,box%coarse_j_lower,box%coarse_k_lower]).and. &
          all(idx<=[box%coarse_i_upper,box%coarse_j_upper,box%coarse_k_upper])) covered=.true.
        end associate
      end do
      if (.not.covered) local_integral=local_integral+volume*local%fields(node)%values(:,i,j,k)
    end do
    end do
    end do
  end do
  call MPI_Allreduce(local_integral,total_integral,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(ierr==MPI_SUCCESS,'inventory reduction')
  call require(all(abs(total_integral-before-sum(all_contributions,dim=2))<=3.e-13_dp*max(1._dp,abs(before))), &
    'independent composite conservation')
  saved=local
  call deposit_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,contributions, &
    zero,spacing,periodic,ok,maximum_received=merge(0,1000000,rank==owners(2)))
  call require(.not.ok,'receive capacity rejection');call unchanged()
  if (rank==nranks-1) contributions(1,1)=ieee_value(0._dp,ieee_quiet_nan)
  call deposit()
  call require(.not.ok,'nonfinite local feedback rejection');call unchanged()
  contributions=1._dp
  if (rank==nranks-1) contributions=huge(1._dp)
  call deposit_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions*1.e-12_dp,contributions, &
    zero,spacing*1.e-12_dp,periodic,ok)
  call require(.not.ok,'density overflow collective rejection');call unchanged()
  if (rank==owners(1)) local%fields(1)%values(:,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
  call deposit()
  call require(.not.ok,'nonfinite initial payload rejection')
  call distribute()
  owners=modulo(owners+1,nranks);call distribute()
  offset=0
  do r=1,nranks-1
    count=2
    if (r==1) count=4
    if (rank==r) contributions=all_contributions(:,offset+1:offset+count)
    offset=offset+count
  end do
  call deposit()
  call require(ok,'changed ownership deposition');call compare_expected()
  call distribute();saved=local
  deallocate(positions,contributions);allocate(positions(3,0),contributions(2,0))
  call deposit_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,contributions, &
    zero,spacing,periodic,ok,maximum_requests=0,maximum_received=0)
  call require(ok,'globally empty feedback');call unchanged()
  if (rank==0) print *, 'PASS: owner-routed AMR spray field deposition'
  call check_source_halo()
  call MPI_Finalize(ierr)
contains
  subroutine check_source_halo()
    type(amr_patch_3d) :: boxes(2)
    type(amr_patch_tree_3d) :: peers
    type(amr_tree_state_3d) :: halo,source,remainder,accepted_source,accepted_remainder
    integer :: map(3),p,n(3)
    real(dp) :: integral(2),remaining(2),pair(2)=[1._dp,-2._dp]
    call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,boxes(1),ok)
    call require(ok,'halo source left patch')
    call initialize_amr_patch_3d(4,4,4,3,4,1,4,1,4,2,boxes(2),ok)
    call require(ok,'halo source right patch')
    call initialize_amr_patch_tree_3d([4,4,4],[1,1],boxes,peers,ok)
    call require(ok,'halo source tree')
    map=[0,0,nranks-1];allocate(halo%fields(3),source%fields(3))
    do p=1,3
      if(map(p)/=rank) cycle
      n=peers%nodes(p)%shape
      allocate(source%fields(p)%values(2,n(1),n(2),n(3)));source%fields(p)%values=0
      if(p==1) then
        source%fields(p)%values=17;cycle
      end if
      allocate(halo%fields(p)%values(2,0:n(1)+1,0:n(2)+1,0:n(3)+1));halo%fields(p)%values=0
      if(p==2) then
        halo%fields(p)%values(:,1,2,2)=pair
        halo%fields(p)%values(:,5,2,2)=2*pair
        halo%fields(p)%values(:,1,0,2)=5*pair
      else
        halo%fields(p)%values(:,0,2,2)=3*pair
        halo%fields(p)%values(:,5,2,2)=4*pair
      end if
    end do
    call sum_mpi_amr_particle_source_halo(MPI_COMM_WORLD,peers,map,halo,2,1,1,[.true.,.false.,.false.], &
      source,remainder,ok)
    call require(ok,'same-level halo source sum')
    integral=0;remaining=0
    do p=1,3
      if(map(p)/=rank) cycle
      if(p==1) then
        call require(all(source%fields(p)%values==17),'other-level source unchanged');cycle
      end if
      integral=integral+sum(sum(sum(source%fields(p)%values,dim=4),dim=3),dim=2)
      remaining=remaining+sum(sum(sum(remainder%fields(p)%values,dim=4),dim=3),dim=2)
      if(p==2) then
        call require(all(source%fields(p)%values(:,1,2,2)==5*pair),'periodic source routes to peer')
        call require(all(source%fields(p)%values(:,4,2,2)==3*pair),'remote left source')
        call require(all(remainder%fields(p)%values(:,1,0,2)==5*pair),'physical-boundary source retained')
      else
        call require(all(source%fields(p)%values(:,1,2,2)==2*pair),'remote right source')
      end if
    end do
    call MPI_Allreduce(MPI_IN_PLACE,integral,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,remaining,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(integral==10*pair).and.all(remaining==5*pair),'source plus remainder preserves every contribution')
    accepted_source=source;accepted_remainder=remainder
    call sum_mpi_amr_particle_source_halo(MPI_COMM_WORLD,peers,map,halo,2,1,1,[.true.,.false.,.false.], &
      source,remainder,ok,maximum_received=0)
    call require(.not.ok,'halo source capacity rejection')
    do p=1,3
      if(map(p)/=rank) cycle
      call require(all(source%fields(p)%values==accepted_source%fields(p)%values),'source capacity rollback')
      if(p>1) call require(all(remainder%fields(p)%values==accepted_remainder%fields(p)%values), &
        'remainder capacity rollback')
    end do
    call sum_mpi_amr_particle_source_halo(MPI_COMM_WORLD,peers,map,halo,2,merge(0,1,rank==0),1, &
      [.true.,.false.,.false.],source,remainder,ok)
    call require(.not.ok,'source width disagreement rejection')
  end subroutine check_source_halo
  subroutine deposit()
    call deposit_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,contributions, &
      zero,spacing,periodic,ok)
  end subroutine
  subroutine compare_expected()
    integer :: q
    do q=1,3
      if (owners(q)==rank) then
        call require(all(local%fields(q)%values==expected%fields(q)%values),'exact serial owner field')
      else
        call require(.not.allocated(local%fields(q)%values),'nonowner storage remains absent')
      end if
    end do
  end subroutine
  subroutine unchanged()
    integer :: q
    do q=1,3
      if (owners(q)/=rank) cycle
      call require(all(local%fields(q)%values==saved%fields(q)%values),'distributed field rollback')
    end do
  end subroutine
  subroutine distribute()
    integer :: q
    if (allocated(local%fields)) deallocate(local%fields)
    allocate(local%fields(3))
    do q=1,3
      if (owners(q)==rank) local%fields(q)%values=full%fields(q)%values
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
