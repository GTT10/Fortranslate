program test_les_dynamic_boundaries
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d, &
    initialize_amr_patch_tree_3d,prolong_amr_patch_tree_3d
  use les_amr_patch_tree_mod, only: advance_amr_patch_tree_les_3d
  use les_dynamic_smagorinsky_mod, only: advance_dynamic_les_ssprk2_3d, &
    compute_dynamic_les_boundary_rhs_3d,compute_dynamic_les_rhs_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(0)
  type(amr_tree_state_3d) :: state,temps
  real(dp) :: root(12,4,4,4),t(4,4,4),q(12),sound,h(3),dt,u(12,4,4,4),ut(4,4,4)
  real(dp) :: rhs(12,4,4,4),old_rhs(12,4,4,4),diffusion,old_diffusion,rate(12)
  real(dp) :: budget(12),reference(12),correction,wall_t(6),wall_v(3,6),pressure(6),inlet(12,6),inlet_t(6)
  real(dp) :: saved(12,4,4,4),saved_t(4,4,4),saved_budget(12),delta(12),scale
  integer :: i,j,k,mode
  character(len=32) :: faces(6)
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  do k=1,4
    do j=1,4
      do i=1,4
        q=0;q(1)=.3_dp;q(2)=2*j;q(3)=3*i;q(5)=101325;q(12)=1
        call reactive_primitive_to_conserved(species,q,root(:,i,j,k),t(i,j,k),sound,ok)
        call require(ok,'initial gas')
      end do
    end do
  end do
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches,tree,ok)
  call require(ok,'root geometry')
  h=.001_dp;dt=1e-9_dp;inlet=spread(root(:,1,1,1),2,6);inlet_t=t(1,1,1)
  call compute_dynamic_les_rhs_3d(species,root,t,h,'box',2,old_rhs,ok,maximum_diffusivity=old_diffusion)
  call require(ok,'periodic oracle')
  call compute_dynamic_les_boundary_rhs_3d(species,root,t,h,'box',2,rhs,diffusion,ok,boundary_rate=rate)
  call require(ok,'periodic wrapper')
  call require(all(rhs==old_rhs).and.diffusion==old_diffusion.and.all(rate==0),'unchanged periodic path')
  do mode=1,4
    faces='outflow';wall_t=1100;wall_v=0;pressure=0
    select case(mode)
    case(1)
      faces(1)='inflow'
    case(2)
      faces(1:2)='slip_wall';faces(3:4)='no_slip_wall';wall_v(1,3:4)=2
    case(3)
      faces(3:4)='isothermal_wall';wall_v(1,3)=1;wall_v(1,4)=3
    case(4)
      pressure(2)=100000
    end select
    call prolong_amr_patch_tree_3d(tree,root,state,ok)
    call require(ok,'oracle state')
    call prolong_amr_patch_tree_3d(tree,reshape(t,[1,4,4,4]),temps,ok)
    call require(ok,'oracle temperature')
    call advance_amr_patch_tree_les_3d(tree,species,state,temps,h,dt,faces, &
      .17_dp,.7_dp,.7_dp,reference,correction,ok,inflow_state=inlet,inflow_temperature=inlet_t, &
      wall_temperature=wall_t,wall_velocity=wall_v,outlet_pressure=pressure,model='dynamic')
    call require(ok,'AMR root oracle')
    u=root;ut=t;budget=123
    call advance_dynamic_les_ssprk2_3d(species,u,ut,h,'box',2,dt,ok, &
      boundary_faces=faces,boundary_change=budget,inflow_state=inlet,inflow_temperature=inlet_t, &
      wall_temperature=wall_t,wall_velocity=wall_v,outlet_pressure=pressure)
    call require(ok,'regular physical Dynamic step')
    ! AMR averages two Euler states; the regular integrator adds the weighted RHS.
    call require(all(abs(u-state%fields(1)%values)<=8*epsilon(1._dp)*max(abs(u),abs(root))), &
      'AMR root state within floating-point accumulation error')
    call require(all(ut==temps%fields(1)%values(1,:,:,:)),'exact AMR root temperature')
    scale=max(maxval(abs(reference)),1e-30_dp)
    call require(maxval(abs(budget-reference))<=1e-12_dp*scale,'AMR boundary budget')
    delta=sum(sum(sum(u-root,dim=4),dim=3),dim=2)*product(h)
    scale=maxval(abs(root))*product(h)*64
    call require(maxval(abs(delta-budget))<=2e-14_dp*scale,'integrated boundary conservation')
    saved=u;saved_t=ut;saved_budget=budget
    call advance_dynamic_les_ssprk2_3d(species,u,ut,h,'box',2,1._dp,ok, &
      boundary_faces=faces,boundary_change=budget,inflow_state=inlet,inflow_temperature=inlet_t, &
      wall_temperature=wall_t,wall_velocity=wall_v,outlet_pressure=pressure)
    call require(.not.ok,'excessive timestep rejects')
    call require(all(u==saved).and.all(ut==saved_t).and.all(budget==saved_budget),'step rollback')
  end do
  rhs=123;diffusion=456;rate=789;wall_t=-1;faces='isothermal_wall'
  call compute_dynamic_les_boundary_rhs_3d(species,root,t,h,'box',2,rhs,diffusion,ok, &
    boundary_faces=faces,boundary_rate=rate,wall_temperature=wall_t)
  call require(.not.ok.and.all(rhs==123).and.diffusion==456.and.all(rate==789),'invalid wall RHS rollback')
  print *, 'PASS: regular Dynamic physical boundaries, AMR parity, conservation and rollback'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(condition) return
    print *, 'FAIL: '//label
    error stop 1
  end subroutine
end program
