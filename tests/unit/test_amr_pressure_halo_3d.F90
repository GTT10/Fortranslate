program test_amr_pressure_halo_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d, &
    initialize_amr_patch_tree_3d,prolong_amr_patch_tree_3d
  use amr_patch_tree_ghost_3d_mod, only: fill_amr_patch_tree_ghosts_3d,apply_amr_tree_pressure_ghosts_3d
  use simulation_config_reactive_3d_mod, only: characteristic_pressure_ghost_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: patch(1)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state
  real(dp) :: q(12),u(12),t,sound,root(13,4,4,4),pressure(6),expected(12),et
  real(dp), allocatable :: ghost(:,:,:,:),saved(:,:,:,:)
  character(len=32) :: faces(6)
  integer :: i,j,k,face,layer
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  q=0;q(1)=.3_dp;q(2)=20;q(5)=101325;q(12)=1
  call reactive_primitive_to_conserved(species,q,u,t,sound,ok)
  call require(ok,'gas')
  do k=1,4
    do j=1,4
      do i=1,4
        root(:12,i,j,k)=u;root(13,i,j,k)=t
      end do
    end do
  end do
  ! Child occupies fine coordinates [2,5] in an eight-cell domain.
  ! Its seven-cell halo crosses both boundaries without touching either face.
  call initialize_amr_patch_3d(4,4,4,2,3,2,3,2,3,2,patch(1),ok)
  call require(ok,'interior child')
  call initialize_amr_patch_tree_3d([4,4,4],[1],patch,tree,ok)
  call require(ok,'tree')
  call prolong_amr_patch_tree_3d(tree,root,state,ok)
  call require(ok,'fields')
  faces='outflow';pressure=0;pressure(1:2)=[102000._dp,100000._dp]
  call fill_amr_patch_tree_ghosts_3d(tree,state,2,faces,ghost,ok,ghost_width=7)
  call require(ok,'wide halo')
  call apply_amr_tree_pressure_ghosts_3d(tree,species,2,faces,pressure,ghost,ok,ghost_width=7)
  call require(ok,'wide pressure correction')
  do i=1,18
    expected=u;et=t;face=0;layer=0
    if(i<6) then
      face=1;layer=6-i
    else if(i>13) then
      face=2;layer=i-13
    end if
    if(face>0) then
      call characteristic_pressure_ghost_3d(species,u,t,face,pressure(face),real(2*layer,dp),expected,et,ok)
      call require(ok,'independent acoustic state')
    end if
    do k=1,18
      do j=1,18
        call require(all(ghost(:12,i,j,k)==expected).and.ghost(13,i,j,k)==et, &
          'physical distance and in-domain donor preservation')
      end do
    end do
  end do
  saved=ghost
  call apply_amr_tree_pressure_ghosts_3d(tree,species,2,faces,pressure,ghost,ok,ghost_width=0)
  call require(.not.ok.and.all(ghost==saved),'invalid width preserves private buffer')
  print '(a)','test_amr_pressure_halo_3d: PASS'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(condition) return
    print '(a)','FAIL: '//label
    error stop 1
  end subroutine
end program
