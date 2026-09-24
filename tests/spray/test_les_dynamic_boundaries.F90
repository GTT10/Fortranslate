program test_les_dynamic_boundaries
  use simulation_config_reactive_3d_mod, only: mach_inlet_ghost_3d
  use les_smagorinsky_3d_mod, only: compute_les_rhs_3d,advance_les_ssprk2_3d
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
  call check_mach_les()
  print *, 'PASS: regular Dynamic physical boundaries, AMR parity, conservation and rollback'
contains
  subroutine check_mach_les()
    character(len=16), parameter :: models(4)=[character(len=16)::'smagorinsky','wale','vreman','dynamic']
    real(dp) :: mach(6),pin(6),target_t(6),uniform(12,4,4,4),uniform_t(4,4,4),ug(12),tg,ref_rate(12)
    integer :: m,face
    q=0;q(1)=.3_dp;q(2:4)=[3._dp,4._dp,5._dp];q(5)=101325;q(12)=1
    call reactive_primitive_to_conserved(species,q,ug,tg,sound,ok)
    call require(ok,'Mach LES initial state')
    do k=1,4
      do j=1,4
        do i=1,4
          uniform(:,i,j,k)=ug;uniform_t(i,j,k)=tg
        end do
      end do
    end do
    target_t=1100
    do m=1,4
      face=1
      if(m==4) face=6
      faces='outflow';faces(face)='inflow';mach=-1;mach(face)=.2_dp;pin=0
      if(m==4) pin(face)=120000
      if(pin(face)>0) then
        call mach_inlet_ghost_3d(species,uniform(:,1,1,1),tg,face,target_t(face),mach(face), &
          ug,inlet_t(face),ok,inlet_pressure=pin(face))
      else
        call mach_inlet_ghost_3d(species,uniform(:,1,1,1),tg,face,target_t(face),mach(face),ug,inlet_t(face),ok)
      end if
      call require(ok,'Mach LES reservoir oracle')
      inlet=spread(ug,2,6);inlet_t=inlet_t(face)
      if(m<4) then
        call compute_les_rhs_3d(species,uniform,uniform_t,h,.17_dp,.7_dp,.7_dp,old_rhs,old_diffusion,ok, &
          boundary_faces=faces,boundary_rate=ref_rate,inflow_state=inlet,inflow_temperature=inlet_t,model=models(m))
      else
        call compute_dynamic_les_boundary_rhs_3d(species,uniform,uniform_t,h,'box',2,old_rhs,old_diffusion,ok, &
          boundary_faces=faces,boundary_rate=ref_rate,inflow_state=inlet,inflow_temperature=inlet_t)
      end if
      call require(ok,'LES reservoir RHS')
      if(m<4) then
        call compute_les_rhs_3d(species,uniform,uniform_t,h,.17_dp,.7_dp,.7_dp,rhs,diffusion,ok, &
          boundary_faces=faces,boundary_rate=rate,inflow_temperature=target_t,model=models(m), &
          inflow_mach=mach,inflow_pressure=pin)
      else
        call compute_dynamic_les_boundary_rhs_3d(species,uniform,uniform_t,h,'box',2,rhs,diffusion,ok, &
          boundary_faces=faces,boundary_rate=rate,inflow_temperature=target_t,inflow_mach=mach,inflow_pressure=pin)
      end if
      call require(ok,'Mach LES RHS')
      call require(all(rhs==old_rhs).and.diffusion==old_diffusion.and.all(rate==ref_rate),'Mach LES initial oracle')
      u=root;ut=t;budget=123
      if(m<4) then
        call advance_les_ssprk2_3d(species,u,ut,h,.17_dp,.7_dp,.7_dp,dt,ok, &
          boundary_faces=faces,boundary_change=budget,inflow_temperature=target_t,model=models(m), &
          inflow_mach=mach,inflow_pressure=pin)
      else
        call advance_dynamic_les_ssprk2_3d(species,u,ut,h,'box',2,dt,ok, &
          boundary_faces=faces,boundary_change=budget,inflow_temperature=target_t,inflow_mach=mach,inflow_pressure=pin)
      end if
      call require(ok,'Mach LES SSPRK2')
      call require(any(u/=root),'Mach LES evolves nonuniform gas')
      delta=sum(sum(sum(u-root,dim=4),dim=3),dim=2)*product(h)
      scale=maxval(abs(root))*product(h)*64
      call require(maxval(abs(delta-budget))<=2e-14_dp*scale,'Mach LES boundary inventory')
      saved=u;saved_t=ut;target_t(face)=-1
      if(m<4) then
        call advance_les_ssprk2_3d(species,u,ut,h,.17_dp,.7_dp,.7_dp,dt,ok, &
          boundary_faces=faces,boundary_change=budget,inflow_temperature=target_t,model=models(m),inflow_mach=mach)
      else
        saved_budget=budget
        call advance_dynamic_les_ssprk2_3d(species,u,ut,h,'box',2,dt,ok, &
          boundary_faces=faces,boundary_change=budget,inflow_temperature=target_t,inflow_mach=mach)
        call require(all(budget==saved_budget),'Mach Dynamic diagnostic rollback')
      end if
      call require(.not.ok.and.all(u==saved).and.all(ut==saved_t),'Mach LES invalid temperature rollback')
      target_t=1100
    end do
    faces='periodic'
    call compute_les_rhs_3d(species,uniform,uniform_t,h,.17_dp,.7_dp,.7_dp,rhs,diffusion,ok, &
      boundary_faces=faces,inflow_mach=mach)
    call require(.not.ok,'static LES rejects Mach on periodic face')
    rhs=123;diffusion=456;rate=789
    call compute_dynamic_les_boundary_rhs_3d(species,uniform,uniform_t,h,'box',2,rhs,diffusion,ok, &
      boundary_faces=faces,boundary_rate=rate,inflow_mach=mach)
    call require(.not.ok.and.all(rhs==123).and.diffusion==456.and.all(rate==789),'Dynamic periodic Mach rollback')
  end subroutine check_mach_les

  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(condition) return
    print *, 'FAIL: '//label
    error stop 1
  end subroutine
end program
