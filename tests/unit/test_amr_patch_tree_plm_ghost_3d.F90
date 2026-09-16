program test_amr_patch_tree_plm_ghost_3d
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_ghost_3d_mod, only: fill_amr_patch_tree_ghosts_3d
  implicit none
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patch(1)
  type(amr_tree_state_3d) :: state,future
  real(dp), allocatable :: ghost(:,:,:,:),saved(:,:,:,:)
  real(dp) :: expected(8),average(8),scale
  integer :: node,n,i,j,k,mode
  character(len=32) :: faces(6)
  logical :: ok
  call initialize_amr_patch_3d(6,6,6,3,4,3,4,3,4,2,patch(1),ok)
  call require(ok,'patch')
  call initialize_amr_patch_tree_3d([6,6,6],[1],patch,tree,ok)
  call require(ok,'tree')
  allocate(state%fields(2),future%fields(2))
  do node=1,2
    n=tree%nodes(node)%shape(1)
    allocate(state%fields(node)%values(8,n,n,n))
    do k=1,n
      do j=1,n
        do i=1,n
          if (node==1) then
            expected=affine(real([i,j,k],dp)-.5_dp)
          else
            expected=affine(2._dp+(real([i,j,k],dp)-.5_dp)/2)
          end if
          state%fields(node)%values(:,i,j,k)=expected
        end do
      end do
    end do
  end do
  future%fields(1)%values=1.2_dp*state%fields(1)%values
  faces='outflow'
  do mode=1,2
    if (mode==1) then
      call fill_amr_patch_tree_ghosts_3d(tree,state,2,faces,ghost,ok,coarse_limiter='mc',species_count=2)
    else
      call fill_amr_patch_tree_ghosts_3d(tree,state,2,faces,ghost,ok,coarse_future=future, &
        coarse_fractions=[.5_dp],coarse_limiter='mc',species_count=2)
    end if
    call require(ok,'PLM ghost fill')
    do k=1,8
      do j=1,8
        do i=1,8
          scale=1
          if (mode==2.and.(any([i,j,k]<3).or.any([i,j,k]>6))) scale=1.1_dp
          expected=scale*affine(2._dp+(real([i,j,k]-2,dp)-.5_dp)/2)
          call require(maxval(abs(ghost(:,i,j,k)-expected))<2e-12_dp,'3D affine exactness and temporal interpolation')
          call require(abs(sum(ghost(6:7,i,j,k))-ghost(1,i,j,k))<2e-14_dp,'species closure')
        end do
      end do
    end do
    average=0
    do k=1,2
      do j=1,2
        do i=1,2
          average=average+ghost(:,i,j,k)/8
        end do
      end do
    end do
    scale=1
    if (mode==2) scale=1.1_dp
    call require(maxval(abs(average-scale*affine([1.5_dp,1.5_dp,1.5_dp])))<2e-12_dp,'donor subdivision mean')
  end do
  saved=ghost
  call fill_amr_patch_tree_ghosts_3d(tree,state,2,faces,ghost,ok,coarse_limiter='unknown',species_count=2)
  call require(.not.ok,'invalid limiter')
  call require(all(ghost==saved),'invalid limiter rollback')
  call fill_amr_patch_tree_ghosts_3d(tree,state,2,faces,ghost,ok,coarse_limiter='mc',species_count=4)
  call require(.not.ok,'invalid species count')
  call require(all(ghost==saved),'invalid species count rollback')
  print *, 'PASS: affine coarse donors, temporal slopes, conservative means, species closure and rollback'
contains
  function affine(x) result(value)
    real(dp), intent(in) :: x(3)
    real(dp) :: value(8),rho
    rho=2+dot_product([.1_dp,.2_dp,.3_dp],x)
    value=[rho,.1_dp*rho,.2_dp*rho,-.3_dp*rho,100*rho,rho/3,2*rho/3,300+dot_product([2._dp,3._dp,4._dp],x)]
  end function
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: '//label
      error stop 1
    end if
  end subroutine
end program
