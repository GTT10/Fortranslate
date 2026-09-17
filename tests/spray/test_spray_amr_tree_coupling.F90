program test_spray_amr_tree_coupling
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use spray_amr_tree_coupling_mod
  use spray_coupling_3d_mod, only: periodic_particle_stencil
  implicit none
  type(amr_patch_3d) :: patches(2)
  type(amr_patch_tree_3d) :: tree,uniform
  type(amr_tree_state_3d) :: field,saved
  integer, allocatable :: cells(:,:),saved_cells(:,:)
  real(dp), allocatable :: weights(:),volumes(:),samples(:,:),saved_samples(:,:),saved_weights(:)
  real(dp) :: root(2,4,4,4),positions(3,4),contributions(2,4),before(2),after(2)
  real(dp) :: zero(3)=0._dp,spacing(3)=1._dp,legacy_weights(8),reference(2),deposited(2,4,4,4)
  integer :: i,j,k,p,q,n,legacy_cells(3,8)
  logical :: ok,periodic(3)=[.true.,.true.,.true.]
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches(:0),uniform,ok)
  call require(ok,'uniform topology')
  do k=1,4
  do j=1,4
  do i=1,4
    root(:,i,j,k)=[real(i+3*j+7*k,dp),real(i*j+k,dp)]
  end do
  end do
  end do
  call prolong_amr_patch_tree_3d(uniform,root,field,ok)
  call require(ok,'uniform field')
  positions(:,1)=[1.9_dp,1.1_dp,1.1_dp]
  positions(:,2)=[2.1_dp,1.1_dp,1.1_dp]
  positions(:,3)=[.9_dp,.9_dp,.9_dp]
  positions(:,4)=[3.9_dp,1.1_dp,1.1_dp]
  call sample_amr_tree_spray_field(uniform,field,positions,zero,spacing,periodic,samples,ok)
  call require(ok,'uniform sample')
  do p=1,4
    call periodic_particle_stencil(positions(:,p),zero,spacing,[4,4,4],legacy_cells,legacy_weights,ok)
    call require(ok,'legacy stencil')
    reference=0._dp
    do q=1,8
      reference=reference+legacy_weights(q)*root(:,legacy_cells(1,q),legacy_cells(2,q),legacy_cells(3,q))
    end do
    call require(all(abs(samples(:,p)-reference)<3.e-14_dp*max(1._dp,abs(reference))),'uniform CIC parity')
    contributions(:,p)=[real(p,dp),-.25_dp*p]
  end do
  deposited=root
  do p=1,4
    call periodic_particle_stencil(positions(:,p),zero,spacing,[4,4,4],legacy_cells,legacy_weights,ok)
    do q=1,8
      i=legacy_cells(1,q);j=legacy_cells(2,q);k=legacy_cells(3,q)
      deposited(:,i,j,k)=deposited(:,i,j,k)+legacy_weights(q)*contributions(:,p)
    end do
  end do
  call deposit_amr_tree_spray_field(uniform,field,positions,contributions,zero,spacing,periodic,ok)
  call require(ok.and.all(abs(field%fields(1)%values-deposited)<3.e-14_dp*max(1._dp,abs(deposited))), &
    'uniform deposit parity')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),ok)
  call require(ok,'fine patch')
  call initialize_amr_patch_3d(4,8,8,1,2,1,2,1,2,2,patches(2),ok)
  call require(ok,'nested patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'nested topology')
  root=1._dp
  call prolong_amr_patch_tree_3d(tree,root,field,ok)
  call require(ok,'constant field')
  call sample_amr_tree_spray_field(tree,field,positions,zero,spacing,periodic,samples,ok)
  call require(ok.and.all(abs(samples-1._dp)<3.e-14_dp),'constant preservation across levels and periodic seam')
  do n=1,3
    field%fields(n)%values=10._dp*n
  end do
  call sample_amr_tree_spray_field(tree,field,positions,zero,spacing,periodic,samples,ok)
  call require(ok,'piecewise field sample')
  call require(all(abs(samples(1,:)-[17._dp,14._dp,27.29_dp,14.64_dp])<3.e-13_dp),'analytic overlap weights')
  call amr_tree_spray_stencil(tree,positions(:,3),zero,spacing,periodic,cells,weights,volumes,ok)
  call require(ok.and.abs(sum(weights)-1._dp)<3.e-14_dp.and.all(weights>0._dp),'positive normalized stencil')
  do q=1,size(weights)
    call require(volumes(q)==1._dp/8._dp**(cells(1,q)-1),'leaf cell volume')
    if (cells(1,q)==1) call require(cells(2,q)>2,'covered root excluded')
    if (cells(1,q)==2) call require(any(cells(2:4,q)>2),'covered fine excluded')
  end do
  saved_cells=cells;saved_weights=weights
  call amr_tree_spray_stencil(tree,positions(:,3),zero,spacing,periodic,cells,weights,volumes,ok,maximum_cells=1)
  call require(.not.ok.and.all(cells==saved_cells).and.all(weights==saved_weights),'stencil capacity rollback')
  call composite_integrals_amr_patch_tree_3d(tree,field,spacing,before,ok)
  call require(ok,'initial composite inventory')
  call deposit_amr_tree_spray_field(tree,field,positions,contributions,zero,spacing,periodic,ok)
  call require(ok,'composite deposit')
  call composite_integrals_amr_patch_tree_3d(tree,field,spacing,after,ok)
  call require(ok.and.all(abs(after-before-sum(contributions,dim=2))<3.e-14_dp*max(1._dp,abs(before))), &
    'volume conservation without covered-cell double counting')
  saved=field;saved_samples=samples;positions(1,4)=5._dp
  call deposit_amr_tree_spray_field(tree,field,positions,contributions,zero,spacing, &
    [.false.,.false.,.false.],ok)
  call require(.not.ok,'late bad position rejected')
  do n=1,3
    call require(all(field%fields(n)%values==saved%fields(n)%values),'complete deposit rollback')
  end do
  call sample_amr_tree_spray_field(tree,field,positions,zero,spacing,[.false.,.false.,.false.],samples,ok)
  call require(.not.ok.and.all(samples==saved_samples),'sample rollback')
  positions(:,1)=0._dp
  call amr_tree_spray_stencil(tree,positions(:,1),zero,spacing,[.false.,.false.,.false.],cells,weights,volumes,ok)
  call require(ok.and.size(weights)==1.and.weights(1)==1._dp,'clipped physical corner normalization')
  call deposit_amr_tree_spray_field(tree,field,positions(:,:0),contributions(:,:0),zero,spacing,periodic,ok)
  call require(ok,'empty deposit')
  do n=1,3
    call require(all(field%fields(n)%values==saved%fields(n)%values),'empty deposit is unchanged')
  end do
  call sample_amr_tree_spray_field(tree,field,positions(:,:0),zero,zero,periodic,samples,ok)
  call require(.not.ok.and.all(samples==saved_samples),'empty sample still validates geometry')
  call sample_amr_tree_spray_field(tree,field,positions(:,:0),zero,spacing,periodic,samples,ok)
  call require(ok.and.size(samples,2)==0,'empty sample')
  contributions=huge(1._dp)
  call deposit_amr_tree_spray_field(tree,field,positions*1.e-12_dp,contributions,zero,spacing*1.e-12_dp,periodic,ok)
  call require(.not.ok,'density overflow rejected before division')
  do n=1,3
    call require(all(field%fields(n)%values==saved%fields(n)%values),'overflow rollback')
  end do
  print *, 'PASS: composite AMR spray sampling and deposition'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
