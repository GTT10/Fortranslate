program test_amr_patch_tree_remap_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,average_down_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d, &
    amr_tree_integer_layout_3d
  use amr_patch_tree_remap_3d_mod, only: remap_amr_patch_tree_3d,remap_amr_patch_tree_reactive_3d
  implicit none
  type(amr_patch_tree_3d) :: old_tree,new_tree,root_tree,bad_tree
  type(amr_patch_3d) :: old_patches(2),new_patches(3)
  type(amr_tree_state_3d) :: old,new,saved,sync,old_t,new_t,saved_t,guess
  type(nasa7_species), allocatable :: species(:)
  integer(int64), allocatable :: origin(:,:),domain(:,:)
  integer, allocatable :: ratios(:)
  real(dp) :: root(2,2,2,2),dense(2,24,24,24),expected(2),before(2),after(2)
  real(dp) :: gas(12,2,2,2),temp(1,2,2,2),q(12),sound,gas_before(12),gas_after(12)
  integer :: node,i,j,k,c,n(3),lo(3),hi(3),factor(3)
  logical :: ok
  call initialize_amr_patch_3d(2,2,2,1,1,1,2,1,2,2,old_patches(1),ok)
  call require(ok,'old child')
  call initialize_amr_patch_3d(2,4,4,1,1,1,1,1,1,2,old_patches(2),ok)
  call require(ok,'old grandchild')
  call initialize_amr_patch_tree_3d([2,2,2],[1,2],old_patches,old_tree,ok)
  call require(ok,'old tree')
  call initialize_amr_patch_3d(2,2,2,2,2,1,2,1,2,3,new_patches(1),ok)
  call initialize_amr_patch_3d(2,2,2,1,1,1,1,1,1,3,new_patches(2),ok)
  call initialize_amr_patch_3d(3,6,6,1,1,1,1,1,1,2,new_patches(3),ok)
  call initialize_amr_patch_tree_3d([2,2,2],[1,1,2],new_patches,new_tree,ok)
  call require(ok,'new tree with changed ratios and branches')
  call initialize_amr_patch_tree_3d([2,2,2],[integer::],old_patches(:0),root_tree,ok)
  call require(ok,'root-only tree')
  root=1
  call prolong_amr_patch_tree_3d(old_tree,root,old,ok)
  call require(ok,'old fields')
  call amr_tree_integer_layout_3d(old_tree,origin,domain,ratios,ok)
  call require(ok,'old layout')
  ! Independent dense raster: parent values are overwritten by descendants.
  do node=1,size(old_tree%nodes)
    n=old_tree%nodes(node)%shape;factor=int(24/domain(:,old_tree%nodes(node)%level+1))
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          old%fields(node)%values(:,i,j,k)=[real(node*10+i+j+2*k,dp),real(node-i+j-k,dp)]
          lo=(int(origin(:,node))+[i-1,j-1,k-1])*factor+1;hi=lo+factor-1
          do c=1,2
            dense(c,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=old%fields(node)%values(c,i,j,k)
          end do
        end do
      end do
    end do
  end do
  call composite_integrals_amr_patch_tree_3d(old_tree,old,[.5_dp,.5_dp,.5_dp],before,ok)
  call require(ok,'old inventory')
  call remap_amr_patch_tree_3d(old_tree,old,new_tree,new,ok)
  call require(ok,'moving/refining/coarsening remap')
  call amr_tree_integer_layout_3d(new_tree,origin,domain,ratios,ok)
  call require(ok,'new layout')
  do node=1,size(new_tree%nodes)
    n=new_tree%nodes(node)%shape;factor=int(24/domain(:,new_tree%nodes(node)%level+1))
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          lo=(int(origin(:,node))+[i-1,j-1,k-1])*factor+1;hi=lo+factor-1
          do c=1,2
            expected(c)=sum(dense(c,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))/real(product(factor),dp)
          end do
          call require(maxval(abs(new%fields(node)%values(:,i,j,k)-expected))<2e-12_dp,'dense overlap oracle')
        end do
      end do
    end do
  end do
  call composite_integrals_amr_patch_tree_3d(new_tree,new,[.5_dp,.5_dp,.5_dp],after,ok)
  call require(ok.and.maxval(abs(after-before))<2e-12_dp,'remap conservation')
  saved=new
  call remap_amr_patch_tree_3d(old_tree,old,root_tree,new,ok)
  call require(ok,'remove all refined levels')
  call composite_integrals_amr_patch_tree_3d(root_tree,new,[.5_dp,.5_dp,.5_dp],after,ok)
  call require(ok.and.maxval(abs(after-before))<2e-12_dp,'coarsening conservation')
  sync=old
  call average_down_amr_patch_tree_3d(old_tree,sync,ok)
  call require(ok,'identity reference')
  call remap_amr_patch_tree_3d(old_tree,old,old_tree,new,ok)
  call require(ok,'identity remap')
  do node=1,size(old_tree%nodes)
    call require(all(new%fields(node)%values==sync%fields(node)%values),'identity exact')
  end do
  new=saved;bad_tree=new_tree;bad_tree%nodes(2)%parent=2
  call remap_amr_patch_tree_3d(old_tree,old,bad_tree,new,ok)
  call require(.not.ok,'invalid topology rejection')
  call same(new,saved,'generic rollback')
  call overflow_rejection()
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
  call require(ok,'gas fields')
  call prolong_amr_patch_tree_3d(old_tree,temp,old_t,ok)
  call require(ok,'temperature fields')
  q(5)=202650;q(6)=.04_dp;q(12)=.74_dp
  call reactive_primitive_to_conserved(species,q,old%fields(3)%values(:,1,1,1), &
    old_t%fields(3)%values(1,1,1,1),sound,ok)
  call require(ok,'leaf perturbation')
  call composite_integrals_amr_patch_tree_3d(old_tree,old,[.5_dp,.5_dp,.5_dp],gas_before,ok)
  call require(ok,'old gas inventory')
  call remap_amr_patch_tree_reactive_3d(old_tree,old,old_t,new_tree,species,new,new_t,ok)
  call require(ok,'reactive remap')
  call composite_integrals_amr_patch_tree_3d(new_tree,new,[.5_dp,.5_dp,.5_dp],gas_after,ok)
  call require(ok.and.maxval(abs(gas_after-gas_before)/max(1.0_dp,abs(gas_before)))<2e-12_dp,'gas conservation')
  call remap_amr_patch_tree_3d(old_tree,old_t,new_tree,guess,ok)
  call require(ok,'temperature guess')
  call require(maxval(abs(new_t%fields(1)%values-guess%fields(1)%values))>1e-5_dp,'EOS recovery changes guess')
  do node=1,size(new_tree%nodes)
    call require(all(new_t%fields(node)%values>0),'positive recovered temperature')
    call require(maxval(abs(sum(new%fields(node)%values(6:,:,:,:),dim=1)- &
      new%fields(node)%values(1,:,:,:)))<1e-12_dp,'species closure')
  end do
  saved=new;saved_t=new_t
  old%fields(3)%values(5,:,:,:)=-1e12_dp
  call remap_amr_patch_tree_reactive_3d(old_tree,old,old_t,new_tree,species,new,new_t,ok)
  call require(.not.ok,'nonphysical candidate rejection')
  call same(new,saved,'reactive rollback')
  call same(new_t,saved_t,'temperature rollback')
contains
  subroutine overflow_rejection()
    type(amr_patch_tree_3d) :: binary,ternary
    type(amr_patch_3d) :: p2(39),p3(23)
    type(amr_tree_state_3d) :: source
    integer :: level
    do level=1,39
      call initialize_amr_patch_3d(2,2,2,1,1,1,1,1,1,2,p2(level),ok)
      call require(ok,'binary deep patch')
    end do
    do level=1,23
      c=3
      if (level==1) c=2
      call initialize_amr_patch_3d(c,c,c,1,1,1,1,1,1,3,p3(level),ok)
      call require(ok,'ternary deep patch')
    end do
    call initialize_amr_patch_tree_3d([2,2,2],[(level,level=1,39)],p2,binary,ok)
    call require(ok,'binary deep tree')
    call initialize_amr_patch_tree_3d([2,2,2],[(level,level=1,23)],p3,ternary,ok)
    call require(ok,'ternary deep tree')
    call prolong_amr_patch_tree_3d(binary,root,source,ok)
    call require(ok,'small local deep state')
    call remap_amr_patch_tree_3d(binary,source,ternary,new,ok)
    call require(.not.ok,'common-coordinate overflow rejection')
    call same(new,saved,'overflow rollback')
  end subroutine
  subroutine same(a,b,label)
    type(amr_tree_state_3d), intent(in) :: a,b
    character(len=*), intent(in) :: label
    integer :: peer
    call require(size(a%fields)==size(b%fields),label)
    do peer=1,size(a%fields)
      call require(all(shape(a%fields(peer)%values)==shape(b%fields(peer)%values)),label)
      call require(all(a%fields(peer)%values==b%fields(peer)%values),label)
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
