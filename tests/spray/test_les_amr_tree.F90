program test_les_amr_tree
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use les_amr_patch_tree_mod, only: advance_amr_patch_tree_les_3d,compute_amr_patch_tree_les_timestep_3d
  use les_smagorinsky_3d_mod, only: advance_les_ssprk2_3d,compute_les_rhs_3d
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_rhs_3d,advance_dynamic_les_ssprk2_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,temperature,saved,saved_t
  type(amr_patch_3d) :: patches(5)
  real(dp) :: root(12,4,4,4),rt(1,4,4,4),serial(12,4,4,4),st(4,4,4),q(12),sound
  real(dp) :: spacing(3),dt,budget(12),serial_budget(12),largest,initial(12),final(12),wall_t(6),wall_v(3,6),outlet(6)
  real(dp) :: chosen,diffusivity,reference_diff,reference_rhs(12,4,4,4),half_dt,half_diff
  real(dp) :: inlet(12,6),inlet_t(6)
  integer :: i,j,k,mode,geometry,node,m
  character(len=32) :: faces(6),model
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  spacing=[.001_dp,.0015_dp,.002_dp];dt=1e-9_dp;wall_t=1100
  do k=1,4
    do j=1,4
      do i=1,4
        q=0;q(1)=.3_dp*(1+.02_dp*(i+j+k));q(5)=101325*(1+.02_dp*k)
        q(2:4)=[10.0_dp+i,20.0_dp+j,-4.0_dp+.5_dp*k]
        q(6)=.025_dp+.003_dp*i;q(9)=.22_dp;q(12)=1-q(6)-q(9)
        call reactive_primitive_to_conserved(species,q,root(:,i,j,k),rt(1,i,j,k),sound,ok)
        call require(ok,'initial state')
      end do
    end do
  end do
  call initialize_amr_patch_3d(4,4,4,1,2,2,3,1,4,2,patches(1),ok)
  call initialize_amr_patch_3d(4,4,4,3,4,2,3,1,4,2,patches(2),ok)
  call initialize_amr_patch_3d(4,4,8,3,4,1,2,1,2,2,patches(3),ok)
  call initialize_amr_patch_3d(4,4,8,1,2,1,2,1,2,2,patches(4),ok)
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,1,2,2,patches(5),ok)
  inlet=spread(root(:,1,1,1),2,6);inlet_t=rt(1,1,1,1)
  do geometry=1,2
    if (geometry==1) then
      call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches(:0),tree,ok)
    else
      call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3,4],patches,tree,ok)
    end if
    call require(ok,'topology')
    do m=1,4
      select case(m)
      case(1)
        model='smagorinsky'
      case(2)
        model='wale'
      case(3)
        model='vreman'
      case(4)
        model='dynamic'
      end select
      do mode=1,4
        if(mode==4.and.m/=4) cycle
        faces='periodic';wall_v=0;outlet=0
        if (mode>=2) faces='outflow'
        if (mode==3) then
          faces(3:4)='isothermal_wall';wall_v(1,3:4)=2
          outlet(2)=102325
        end if
        if(m==4.and.mode==3) then
          ! Keep the seven-layer acoustic extrapolation thermodynamically admissible.
          faces='outflow';faces(1)='inflow';outlet=0;outlet(2)=110000;wall_v=0
        end if
        if(mode==4) then
          faces(1:2)='slip_wall';faces(3:4)='no_slip_wall';faces(5:6)='isothermal_wall'
          wall_v(1,3:6)=2;outlet=0
        end if
        call prolong_amr_patch_tree_3d(tree,root,state,ok)
        call require(ok,'state initialization')
        call prolong_amr_patch_tree_3d(tree,rt,temperature,ok)
        call require(ok,'temperature initialization')
        call composite_integrals_amr_patch_tree_3d(tree,state,spacing,initial,ok)
        call require(ok,'initial inventory')
        chosen=-1;diffusivity=-2
        call compute_amr_patch_tree_les_timestep_3d(tree,species,state,temperature,spacing,.05_dp,faces, &
          .17_dp,.7_dp,.7_dp,chosen,diffusivity,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
          coarse_limiter='mc',outlet_pressure=outlet,model=model,ci=.1_dp,inflow_state=inlet,inflow_temperature=inlet_t)
        call require(ok.and.chosen>0.and.diffusivity>0,'hierarchy LES interval')
        call compute_amr_patch_tree_les_timestep_3d(tree,species,state,temperature,spacing,.025_dp,faces, &
          .17_dp,.7_dp,.7_dp,half_dt,half_diff,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
          coarse_limiter='mc',outlet_pressure=outlet,model=model,ci=.1_dp,inflow_state=inlet,inflow_temperature=inlet_t)
        call require(ok.and.half_dt==.5_dp*chosen.and.half_diff==diffusivity,'CFL scales chosen interval')
        if (geometry==1.and.(m/=4.or.mode==1)) then
          if(m==4) then
            call compute_dynamic_les_rhs_3d(species,root,rt(1,:,:,:),spacing,'box',2,reference_rhs,ok, &
              maximum_diffusivity=reference_diff)
            call require(ok.and.diffusivity>=reference_diff,'conservative dynamic diffusion bound')
            call require(chosen*reference_diff*sum(1._dp/spacing**2)<=.05_dp*(1+1e-12_dp), &
              'dynamic root interval bound')
          else
          call compute_les_rhs_3d(species,root,rt(1,:,:,:),spacing,.17_dp,.7_dp,.7_dp,reference_rhs,reference_diff,ok, &
            boundary_faces=faces,wall_temperature=wall_t,wall_velocity=wall_v,model=model,ci=.1_dp,outlet_pressure=outlet)
          call require(ok,'root interval oracle')
          call require(abs(chosen*reference_diff*sum(1._dp/spacing**2)-.05_dp)<1e-12_dp,'root interval parity')
          end if
        end if
        dt=min(chosen,1e-9_dp)
        saved=state
        call advance_amr_patch_tree_les_3d(tree,species,state,temperature,spacing,dt,faces, &
          .17_dp,.7_dp,.7_dp,budget,largest,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
          coarse_limiter='mc',outlet_pressure=outlet,model=model,ci=.1_dp,inflow_state=inlet,inflow_temperature=inlet_t)
        call require(ok.and.largest>=0,'AMR LES step')
        if(m==4.and.mode>1) call require(maxval(abs(budget))>1e-25_dp,'dynamic open boundary transfers inventory')
        call require(any(state%fields(1)%values/=saved%fields(1)%values),'AMR LES evolves')
        if (geometry==1.and.(m/=4.or.mode==1)) then
          serial=root;st=rt(1,:,:,:)
          if(m==4) then
            call advance_dynamic_les_ssprk2_3d(species,serial,st,spacing,'box',2,dt,ok)
            serial_budget=0
          else
          call advance_les_ssprk2_3d(species,serial,st,spacing,.17_dp,.7_dp,.7_dp,dt,ok, &
            boundary_faces=faces,boundary_change=serial_budget,wall_temperature=wall_t,wall_velocity=wall_v, &
            model=model,ci=.1_dp,outlet_pressure=outlet)
          end if
          call require(ok,'regular LES oracle')
          call require(maxval(abs(state%fields(1)%values-serial)/max(1._dp,abs(serial)))<2e-13_dp,'root state parity')
          call require(maxval(abs(temperature%fields(1)%values(1,:,:,:)-st))<1e-9_dp,'root temperature parity')
          call require(maxval(abs(budget-serial_budget))<1e-14_dp,'root ledger parity')
        end if
        call composite_integrals_amr_patch_tree_3d(tree,state,spacing,final,ok)
        call require(ok,'final inventory')
        call require(maxval(abs(final-initial-budget)/max(1e-12_dp,abs(initial)))<2e-12_dp,'composite conservation')
        chosen=-7;diffusivity=-9
        call compute_amr_patch_tree_les_timestep_3d(tree,species,state,temperature,spacing,.21_dp,faces, &
          .17_dp,.7_dp,.7_dp,chosen,diffusivity,ok,model=model)
        call require(.not.ok.and.chosen==-7.and.diffusivity==-9,'invalid CFL preserves outputs')
        saved=state;saved_t=temperature;budget=-17;largest=-19
        call advance_amr_patch_tree_les_3d(tree,species,state,temperature,spacing,1._dp,faces, &
          .17_dp,.7_dp,.7_dp,budget,largest,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
          coarse_limiter='mc',outlet_pressure=outlet,model=model,ci=.1_dp,inflow_state=inlet,inflow_temperature=inlet_t)
        call require(.not.ok,'excessive interval rejected')
        do node=1,size(tree%nodes)
          call require(all(state%fields(node)%values==saved%fields(node)%values),'state rollback')
          call require(all(temperature%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
        end do
        call require(all(budget==-17).and.largest==-19,'diagnostics rollback')
      end do
    end do
  end do
  print *, 'PASS: AMR static/dynamic LES, root parity, composite conservation and rollback'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (condition) return
    print *, 'FAIL: ',label
    error stop 1
  end subroutine
end program
