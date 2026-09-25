program test_spray_eb_tree_stencil
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use spray_eb_tree_stencil_mod, only: eb_tree_spray_stencil
  implicit none
  type(amr_patch_3d) :: patch(1)
  type(amr_patch_tree_3d) :: tree
  type(eb_geometry_3d) :: finest,g(2)
  real(dp) :: phi(0:8,0:8,0:8),position(3),zero(3)=0,spacing(3)=1,expected,total
  real(dp), allocatable :: w(:),v(:),saved_w(:),saved_v(:)
  integer, allocatable :: cells(:,:),saved(:,:)
  integer :: i,j,q,node
  logical :: ok
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
  call run()
  call require(ok.and.all(w>0).and.abs(sum(w)-1)<1e-14_dp,'normalized cut-interface weights')
  call require(any(cells(1,:)==1).and.any(cells(1,:)==2),'coarse and fine support')
  total=0
  do q=1,size(w)
    node=cells(1,q)
    if(node==1) then
      call require(cells(2,q)>=3,'hidden parent excluded')
      expected=.75_dp
    else
      expected=merge(.0625_dp,.125_dp,cells(3,q)==3)
    end if
    call require(abs(v(q)-expected)<1e-14_dp,'physical cut-cell volume')
    total=total+w(q)*7._dp
  end do
  call require(abs(total-7)<1e-14_dp,'constant preservation')
  saved=cells;saved_w=w;saved_v=v
  position(2)=1._dp
  call run()
  call require(.not.ok.and.all(cells==saved).and.all(w==saved_w).and.all(v==saved_v),'solid-point rollback')
  position=[3.9_dp,3.9_dp,3.9_dp]
  call run()
  call require(ok.and.abs(sum(w)-1)<1e-14_dp.and.all(cells(1,:)==1),'truncated root-boundary support')
  position=[3.9_dp,3.9_dp,3.9_dp]
  call eb_tree_spray_stencil(tree,g,phi,position,zero,spacing,cells,w,v,ok,periodic=[.true.,.false.,.true.])
  call require(ok.and.abs(sum(w)-1)<1e-14_dp,'periodic corner normalized support')
  call require(any(cells(1,:)==1).and.any(cells(1,:)==2),'periodic seam reaches refined opposite side')
  saved=cells;saved_w=w;saved_v=v
  position=position+[4._dp,0._dp,-4._dp]
  call eb_tree_spray_stencil(tree,g,phi,position,zero,spacing,cells,w,v,ok,periodic=[.true.,.false.,.true.])
  call require(ok.and.all(cells==saved).and.maxval(abs(w-saved_w))<1e-14_dp.and.all(v==saved_v), &
    'periodic image invariance')
  call eb_tree_spray_stencil(tree,g,phi,position,zero,spacing,cells,w,v,ok,periodic=[.true.,.true.,.true.])
  call require(.not.ok.and.all(cells==saved).and.maxval(abs(w-saved_w))<1e-14_dp.and.all(v==saved_v), &
    'incompatible periodic normal rejects without output changes')
  do j=0,8
  do i=0,8
    phi(i,j,:)=abs(real(i,dp)*.5_dp-2._dp)-.6_dp
  end do
  end do
  call build_level_set_eb_geometry_3d(phi,zero,[4._dp,4._dp,4._dp],finest,ok)
  call require(ok,'periodic geometry with central solid slab')
  call coarsen_eb_geometry_3d(finest,2,g(1),ok)
  call require(ok,'periodic slab root')
  call extract_eb_geometry_3d(finest,[1,1,1],[4,8,8],g(2),ok)
  call require(ok,'periodic slab child')
  position=[3.9_dp,3.9_dp,3.9_dp]
  call eb_tree_spray_stencil(tree,g,phi,position,zero,spacing,cells,w,v,ok,periodic=[.true.,.false.,.true.])
  call require(ok.and.any(cells(1,:)==2),'visibility follows seam instead of crossing central solid slab')
  saved=cells;saved_w=w;saved_v=v
  phi(0,8,8)=phi(0,8,8)+.1_dp
  call eb_tree_spray_stencil(tree,g,phi,position,zero,spacing,cells,w,v,ok,periodic=[.true.,.false.,.true.])
  call require(.not.ok.and.all(cells==saved).and.all(w==saved_w).and.all(v==saved_v), &
    'nonperiodic nodal field rolls back outputs')
  print *, 'PASS: EB composite particle stencil at a cut coarse/fine interface'
contains
  subroutine run()
    call eb_tree_spray_stencil(tree,g,phi,position,zero,spacing,cells,w,v,ok)
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
