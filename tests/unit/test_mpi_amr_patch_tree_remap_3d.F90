program test_mpi_amr_patch_tree_remap_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use amr_patch_tree_remap_3d_mod, only: remap_amr_patch_tree_3d,remap_amr_patch_tree_reactive_3d
  use mpi_amr_patch_tree_remap_3d_mod, only: remap_mpi_amr_tree_3d,remap_mpi_amr_tree_reactive_3d
  use mpi_amr_patch_tree_distribution_3d_mod, only: gather_mpi_amr_tree_3d
  implicit none
  type(amr_patch_tree_3d) :: old_tree,new_tree,root_tree,target,changed
  type(amr_patch_3d) :: old_patches(2),new_patches(3)
  type(amr_tree_state_3d) :: old,old_t,local,local_t,out,out_t,oracle,oracle_t,saved,saved_t,gathered,source
  type(nasa7_species), allocatable :: species(:)
  integer, allocatable :: old_owners(:),new_owners(:),bad_owners(:)
  real(dp) :: root(2,2,2,2),before(2),after(2),gas(12,2,2,2),temp(1,2,2,2),q(12),sound
  integer :: rank,nranks,ierr,mode,geometry,node,i,j,k,n(3),c
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks>=2,'rank count')
  call initialize_amr_patch_3d(2,2,2,1,1,1,2,1,2,2,old_patches(1),ok)
  call require(ok,'old child')
  call initialize_amr_patch_3d(2,4,4,1,1,1,1,1,1,2,old_patches(2),ok)
  call require(ok,'old grandchild')
  call initialize_amr_patch_tree_3d([2,2,2],[1,2],old_patches,old_tree,ok)
  call require(ok,'old tree')
  call initialize_amr_patch_3d(2,2,2,2,2,1,2,1,2,3,new_patches(1),ok)
  call require(ok,'new child')
  call initialize_amr_patch_3d(2,2,2,1,1,1,1,1,1,3,new_patches(2),ok)
  call require(ok,'new sibling')
  call initialize_amr_patch_3d(3,6,6,1,1,1,1,1,1,2,new_patches(3),ok)
  call require(ok,'new grandchild')
  call initialize_amr_patch_tree_3d([2,2,2],[1,1,2],new_patches,new_tree,ok)
  call require(ok,'new tree')
  call initialize_amr_patch_tree_3d([2,2,2],[integer::],old_patches(:0),root_tree,ok)
  call require(ok,'root tree')
  root=1
  call prolong_amr_patch_tree_3d(old_tree,root,old,ok)
  call require(ok,'old fields')
  do node=1,size(old_tree%nodes)
    n=old_tree%nodes(node)%shape
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          old%fields(node)%values(:,i,j,k)=[real(node*10+i+j+2*k,dp),real(node-i+j-k,dp)]/7.0_dp
        end do
      end do
    end do
  end do
  call composite_integrals_amr_patch_tree_3d(old_tree,old,[.5_dp,.5_dp,.5_dp],before,ok)
  call require(ok,'initial inventory')
  do mode=1,2
    old_owners=mod([0,1,2],nranks)
    if (mode==2) old_owners=0
    call localize(old,old_owners,local)
    source=local
    do geometry=1,3
      target=new_tree
      if (geometry==2) target=root_tree
      if (geometry==3) target=old_tree
      new_owners=[(mod(node,nranks),node=1,size(target%nodes))]
      if (mode==2) new_owners=nranks-1
      call remap_amr_patch_tree_3d(old_tree,old,target,oracle,ok)
      call require(ok,'serial remap')
      call remap_mpi_amr_tree_3d(old_tree,old_owners,local,target,new_owners,2,out,MPI_COMM_WORLD,ok)
      call require(ok,'distributed remap')
      call compare_oracle(out,oracle,new_owners)
      call same(local,source,'source unchanged')
      call gather_mpi_amr_tree_3d(target,new_owners,out,gathered,MPI_COMM_WORLD,0,ok)
      call require(ok,'gather for inventory check')
      if (rank==0) then
        call composite_integrals_amr_patch_tree_3d(target,gathered,[.5_dp,.5_dp,.5_dp],after,ok)
        call require(ok.and.maxval(abs(after-before))<2e-12_dp,'conservation')
      end if
      saved=out
      bad_owners=new_owners
      if (rank==0) bad_owners(1)=mod(new_owners(1)+1,nranks)
      call remap_mpi_amr_tree_3d(old_tree,old_owners,local,target,bad_owners,2,out,MPI_COMM_WORLD,ok)
      call require(.not.ok,'inconsistent target ownership rejected')
      call same(out,saved,'metadata rollback')
      if (rank==old_owners(3)) local%fields(3)%values(1,1,1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
      call remap_mpi_amr_tree_3d(old_tree,old_owners,local,target,new_owners,2,out,MPI_COMM_WORLD,ok)
      call require(.not.ok,'nonfinite input rejected collectively')
      call same(out,saved,'nonfinite rollback')
      local=source
    end do
  end do
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  q=0;q(1)=.3_dp;q(5)=101325;q(6)=.03_dp;q(9)=.22_dp;q(12)=.75_dp
  call reactive_primitive_to_conserved(species,q,gas(:,1,1,1),temp(1,1,1,1),sound,ok)
  call require(ok,'mixture')
  do c=1,12
    gas(c,:,:,:)=gas(c,1,1,1)
  end do
  temp=temp(1,1,1,1)
  call prolong_amr_patch_tree_3d(old_tree,gas,old,ok)
  call require(ok,'gas')
  call prolong_amr_patch_tree_3d(old_tree,temp,old_t,ok)
  call require(ok,'temperature')
  q(5)=202650;q(6)=.04_dp;q(12)=.74_dp
  call reactive_primitive_to_conserved(species,q,old%fields(3)%values(:,1,1,1), &
    old_t%fields(3)%values(1,1,1,1),sound,ok)
  call require(ok,'leaf perturbation')
  call remap_amr_patch_tree_reactive_3d(old_tree,old,old_t,new_tree,species,oracle,oracle_t,ok)
  call require(ok,'serial reactive remap')
  old_owners=mod([0,1,2],nranks);new_owners=mod([1,2,0,1],nranks)
  call localize(old,old_owners,local)
  call localize(old_t,old_owners,local_t)
  call remap_mpi_amr_tree_reactive_3d(old_tree,old_owners,local,local_t,new_tree,new_owners, &
    species,out,out_t,MPI_COMM_WORLD,ok)
  call require(ok,'distributed reactive remap')
  call compare_oracle(out,oracle,new_owners)
  call compare_oracle(out_t,oracle_t,new_owners)
  saved=out;saved_t=out_t
  if (rank==old_owners(3)) local%fields(3)%values(5,:,:,:)=-1e12_dp
  call remap_mpi_amr_tree_reactive_3d(old_tree,old_owners,local,local_t,new_tree,new_owners, &
    species,out,out_t,MPI_COMM_WORLD,ok)
  call require(.not.ok,'late EOS failure rejected collectively')
  call same(out,saved,'EOS state rollback')
  call same(out_t,saved_t,'EOS temperature rollback')
  if (rank==0) print *, 'PASS: sparse remap, changing ratios/owners, empty ranks, conservation, exact serial parity and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine localize(full,owners,part)
    type(amr_tree_state_3d), intent(in) :: full
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(out) :: part
    integer :: peer
    allocate(part%fields(size(full%fields)))
    do peer=1,size(owners)
      if (owners(peer)==rank) part%fields(peer)=full%fields(peer)
    end do
  end subroutine
  subroutine compare_oracle(part,full,owners)
    type(amr_tree_state_3d), intent(in) :: part,full
    integer, intent(in) :: owners(:)
    integer :: peer
    call require(size(part%fields)==size(owners),'field count')
    do peer=1,size(owners)
      call require(allocated(part%fields(peer)%values).eqv.(owners(peer)==rank),'owner-only allocation')
      if (owners(peer)/=rank) cycle
      call require(all(part%fields(peer)%values==full%fields(peer)%values),'exact serial parity')
    end do
  end subroutine
  subroutine same(actual,expected,label)
    type(amr_tree_state_3d), intent(in) :: actual,expected
    character(len=*), intent(in) :: label
    integer :: peer
    call require(size(actual%fields)==size(expected%fields),label//' count')
    do peer=1,size(actual%fields)
      call require(allocated(actual%fields(peer)%values).eqv.allocated(expected%fields(peer)%values),label//' allocation')
      if (.not.allocated(expected%fields(peer)%values)) cycle
      call require(all(actual%fields(peer)%values==expected%fields(peer)%values),label//' values')
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: status
    if (condition) return
    print *, 'FAIL rank ',rank,': '//label
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
end program
