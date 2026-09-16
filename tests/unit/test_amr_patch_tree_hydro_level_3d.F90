program test_amr_patch_tree_hydro_level_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use amr_patch_tree_hydro_3d_mod, only: advance_amr_patch_tree_hydro_3d,advance_amr_patch_tree_hydro_level_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: state,temp,reference,reference_t,start,finish,initial,initial_t,saved,saved_t
  type(amr_tree_flux_3d), allocatable :: flux(:),saved_flux(:)
  real(dp) :: root(12,4,4,4),rt(1,4,4,4),q(12),sound,budget(12),correction,expected(12),h(3)
  real(dp), parameter :: spacing(3)=[.001_dp,.001_dp,.001_dp],dt=1e-9_dp
  integer :: i,j,k,node,mode
  character(len=32) :: faces(6),reconstruction
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  do k=1,4
    do j=1,4
      do i=1,4
        q=0;q(1)=.3_dp+.002_dp*i;q(2:4)=[2._dp,-1._dp,3._dp];q(5)=101325+100*j
        q(6)=.03_dp;q(9)=.22_dp;q(12)=.75_dp
        call reactive_primitive_to_conserved(species,q,root(:,i,j,k),rt(1,i,j,k),sound,ok)
        call require(ok,'initial primitive')
      end do
    end do
  end do
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches(:0),tree,ok)
  call require(ok,'root tree')
  do mode=1,2
    faces='periodic';reconstruction='pcm'
    if (mode==2) then
      faces='outflow';faces(3:4)='slip_wall';reconstruction='characteristic_plm'
    end if
    call initialize()
    reference=state;reference_t=temp
    call advance_amr_patch_tree_hydro_3d(tree,species,reference,reference_t,spacing,dt,faces, &
      reconstruction,'mc','hllc',budget,correction,ok)
    call require(ok,'synchronous root reference')
    call advance_amr_patch_tree_hydro_level_3d(tree,species,state,temp,0,spacing,dt,faces, &
      reconstruction,'mc','hllc',flux,ok)
    call require(ok,'root level step')
    call require(maxval(abs(state%fields(1)%values-reference%fields(1)%values))<1e-11_dp,'root state oracle')
    call require(maxval(abs(temp%fields(1)%values-reference_t%fields(1)%values))<1e-10_dp,'root temperature oracle')
    call check_flux(1,spacing)
  end do
  call initialize_amr_patch_3d(4,4,4,1,2,2,3,2,3,2,patches(1),ok)
  call require(ok,'first patch')
  call initialize_amr_patch_3d(4,4,4,3,4,2,3,2,3,2,patches(2),ok)
  call require(ok,'second patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1],patches,tree,ok)
  call require(ok,'peer tree')
  call initialize()
  allocate(start%fields(3),finish%fields(3))
  allocate(start%fields(1)%values(13,4,4,4))
  start%fields(1)%values(:12,:,:,:)=root;start%fields(1)%values(13,:,:,:)=rt(1,:,:,:)
  finish=start;finish%fields(1)%values(:12,:,:,:)=1.01_dp*root
  faces='periodic'
  call advance_amr_patch_tree_hydro_level_3d(tree,species,state,temp,1,spacing,dt,faces, &
    'characteristic_plm','mc','hllc',flux,ok,start,finish,[.25_dp],[.5_dp])
  call require(ok,'fine interval with coarse endpoints')
  call require(all(state%fields(1)%values==initial%fields(1)%values),'inactive root unchanged')
  call require(all(temp%fields(1)%values==initial_t%fields(1)%values),'inactive temperature unchanged')
  call require(.not.allocated(flux(1)%x),'inactive flux absent')
  do node=2,3
    call check_flux(node,spacing/2)
  end do
  call require(maxval(abs(flux(2)%x(:,5,:,:)-flux(3)%x(:,1,:,:)))<1e-12_dp,'synchronous peer flux')
  reference=initial;reference_t=initial_t
  call advance_amr_patch_tree_hydro_level_3d(tree,species,reference,reference_t,1,spacing,dt,faces, &
    'characteristic_plm','mc','hllc',saved_flux,ok,start,finish,[0._dp],[0._dp])
  call require(ok,'frozen donor comparison')
  call require(maxval(abs(reference%fields(2)%values-state%fields(2)%values))>1e-8_dp,'coarse time affects fine solution')
  saved=state;saved_t=temp;saved_flux=flux
  call advance_amr_patch_tree_hydro_level_3d(tree,species,state,temp,1,spacing,dt,faces, &
    'pcm','mc','hllc',flux,ok,start,finish,[.5_dp],[.25_dp])
  call require(.not.ok,'backward time interval rejected')
  call unchanged()
  call advance_amr_patch_tree_hydro_level_3d(tree,species,state,temp,1,spacing,1._dp,faces, &
    'pcm','mc','hllc',flux,ok,start,finish,[0._dp],[1._dp])
  call require(.not.ok,'failed stage rejected')
  call unchanged()
  print *, 'PASS: level RK2 root oracle, peer updates, time donors, integrated flux divergence and rollback'
contains
  subroutine initialize()
    call prolong_amr_patch_tree_3d(tree,root,state,ok)
    call require(ok,'initialize state')
    call prolong_amr_patch_tree_3d(tree,rt,temp,ok)
    call require(ok,'initialize temperature')
    initial=state;initial_t=temp
  end subroutine
  subroutine check_flux(peer,cell_h)
    integer, intent(in) :: peer
    real(dp), intent(in) :: cell_h(3)
    integer :: a,b,c,n(3)
    n=tree%nodes(peer)%shape
    do c=1,n(3)
      do b=1,n(2)
        do a=1,n(1)
          expected=initial%fields(peer)%values(:,a,b,c) &
            -(flux(peer)%x(:,a+1,b,c)-flux(peer)%x(:,a,b,c))/cell_h(1) &
            -(flux(peer)%y(:,a,b+1,c)-flux(peer)%y(:,a,b,c))/cell_h(2) &
            -(flux(peer)%z(:,a,b,c+1)-flux(peer)%z(:,a,b,c))/cell_h(3)
          call require(maxval(abs(expected-state%fields(peer)%values(:,a,b,c))/max(1._dp,abs(expected)))<3e-15_dp, &
            'integrated flux reproduces complete RK update')
        end do
      end do
    end do
  end subroutine
  subroutine unchanged()
    integer :: peer
    do peer=1,size(tree%nodes)
      call require(all(state%fields(peer)%values==saved%fields(peer)%values),'state rollback')
      call require(all(temp%fields(peer)%values==saved_t%fields(peer)%values),'temperature rollback')
      if (.not.allocated(saved_flux(peer)%x)) cycle
      call require(all(flux(peer)%x==saved_flux(peer)%x),'x flux rollback')
      call require(all(flux(peer)%y==saved_flux(peer)%y),'y flux rollback')
      call require(all(flux(peer)%z==saved_flux(peer)%z),'z flux rollback')
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: '//label
      error stop 1
    end if
  end subroutine
end program
