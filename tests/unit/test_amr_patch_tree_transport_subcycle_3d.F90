program test_amr_patch_tree_transport_subcycle_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,average_down_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use amr_patch_tree_transport_3d_mod, only: advance_amr_patch_tree_transport_3d, &
    compute_amr_patch_tree_transport_timestep_3d
  use amr_patch_tree_transport_subcycle_3d_mod, only: advance_amr_patch_tree_subcycle_transport_3d, &
    compute_amr_patch_tree_subcycle_transport_timestep_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(3)
  type(amr_tree_state_3d) :: state,temp,reference,reference_t,saved,saved_t,synchronized
  real(dp) :: root(12,4,4,4),rt(1,4,4,4),q(12),sound,initial(12),final(12),budget(12),maximum
  real(dp) :: reference_budget(12),reference_max,dt,common_dt,diffusivity,theta,reference_theta,wall_t(6),wall_v(3,6)
  real(dp), parameter :: spacing(3)=[.001_dp,.0015_dp,.002_dp]
  real(dp) :: outlet_pressure(6),mach(6),tin(6),pin(6)
  integer(int64), allocatable :: counts(:)
  integer :: geometry,mode,node,i,j,k
  character(len=32) :: faces(6),reconstruction
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  wall_t=1100;wall_v=0
  do k=1,4
    do j=1,4
      do i=1,4
        q=0;q(1)=.3_dp+.001_dp*(i+j+k);q(2:4)=[2._dp,-1._dp,3._dp];q(5)=101325+100*j
        q(6)=.03_dp;q(9)=.22_dp;q(12)=.75_dp
        call reactive_primitive_to_conserved(species,q,root(:,i,j,k),rt(1,i,j,k),sound,ok)
        call require(ok,'primitive')
      end do
    end do
  end do
  call initialize_amr_patch_3d(4,4,4,1,2,2,3,1,2,2,patches(1),ok)
  call require(ok,'left patch')
  call initialize_amr_patch_3d(4,4,4,3,4,2,3,1,2,2,patches(2),ok)
  call require(ok,'right patch')
  call initialize_amr_patch_3d(4,4,4,4,4,1,1,1,1,3,patches(3),ok)
  call require(ok,'ratio-three grandchild')
  do geometry=1,2
    if (geometry==1) then
      call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches(:0),tree,ok)
    else
      call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],patches,tree,ok)
    end if
    call require(ok,'tree')
    if (allocated(counts)) deallocate(counts)
    allocate(counts(tree%level_count()))
    do mode=1,5
      faces='periodic';reconstruction='pcm';wall_v=0
      if (mode>=2) then
        faces='outflow';faces(3:4)='isothermal_wall';reconstruction='characteristic_plm'
        wall_v(1,3:4)=2
      end if
      outlet_pressure=0;mach=-1;tin=700;pin=0
      if (mode==3) outlet_pressure([2,6])=102325
      if(mode>3) then
        outlet_pressure=0;faces(1)='inflow';mach(1)=.2_dp
        if(mode==5) pin(1)=120000
      end if
      call prolong_amr_patch_tree_3d(tree,root,state,ok)
      call require(ok,'initial state')
      call prolong_amr_patch_tree_3d(tree,rt,temp,ok)
      call require(ok,'initial temperature')
      call composite_integrals_amr_patch_tree_3d(tree,state,spacing,initial,ok)
      call require(ok,'initial inventory')
      call compute_amr_patch_tree_subcycle_transport_timestep_3d(tree,species,transport,state,temp,spacing, &
        .2_dp,.true.,.true.,.true.,dt,diffusivity,ok,faces,wall_t)
      call require(ok,'subcycle CFL')
      call compute_amr_patch_tree_transport_timestep_3d(tree,species,transport,state,temp,spacing, &
        .2_dp,.true.,.true.,.true.,common_dt,diffusivity,ok,faces,wall_t)
      call require(ok,'common CFL')
      if (geometry==1) then
        call require(dt==common_dt,'root timestep equality')
      else
        call require(dt>=common_dt,'scaled coarse diffusion bound')
        if (mode==1) call require(dt>35*common_dt.and.dt<37*common_dt,'ratio-squared diffusion scaling')
      end if
      dt=min(dt,1e-8_dp)
      reference=state;reference_t=temp
      if (geometry==1) then
        call advance_amr_patch_tree_transport_3d(tree,species,transport,reference,reference_t,spacing,dt,faces, &
          .true.,.true.,.true.,.true.,reference_budget,reference_max,reference_theta,ok, &
          wall_temperature=wall_t,wall_velocity=wall_v,coarse_limiter='mc',outlet_pressure=outlet_pressure, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
        call require(ok,'root oracle')
      end if
      call advance_amr_patch_tree_subcycle_transport_3d(tree,species,transport,state,temp,spacing,dt,faces, &
        .true.,.true.,.true.,.true.,budget,maximum,theta,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
        coarse_limiter='mc',level_steps=counts,outlet_pressure=outlet_pressure, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
      call require(ok,'recursive transport step')
      call require(theta>=0.and.theta<=1,'transport limiter bounds')
      if (geometry==1) then
        call require(all(counts==1),'root step count')
        call require(maxval(abs(state%fields(1)%values-reference%fields(1)%values))<2e-11_dp,'root state equality')
        call require(maxval(abs(temp%fields(1)%values-reference_t%fields(1)%values))<1e-9_dp,'root temperature')
        call require(maxval(abs(budget-reference_budget))<1e-18_dp,'root boundary equality')
      else
        call require(all(counts==[1_int64,4_int64,36_int64]),'nested substep counts')
      end if
      call composite_integrals_amr_patch_tree_3d(tree,state,spacing,final,ok)
      call require(ok,'final inventory')
      call require(maxval(abs(final-initial-budget)/max(1e-20_dp,abs(initial)))<2e-10_dp,'composite conservation')
      synchronized=state
      call average_down_amr_patch_tree_3d(tree,synchronized,ok)
      call require(ok,'restriction oracle')
      do node=1,size(tree%nodes)
        call require(all(synchronized%fields(node)%values==state%fields(node)%values),'synchronized parent cells')
        call require(all(temp%fields(node)%values>0),'positive temperatures')
        call require(maxval(abs(sum(state%fields(node)%values(6:,:,:,:),dim=1)- &
          state%fields(node)%values(1,:,:,:)))<1e-12_dp,'species closure')
      end do
      saved=state;saved_t=temp;budget=-17;maximum=-19;counts=-23;theta=-29
      call advance_amr_patch_tree_subcycle_transport_3d(tree,species,transport,state,temp,spacing,1._dp,faces, &
        .true.,.true.,.true.,.true.,budget,maximum,theta,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
        coarse_limiter='mc',level_steps=counts,outlet_pressure=outlet_pressure, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
      call require(.not.ok,'nonphysical step rejected')
      do node=1,size(tree%nodes)
        call require(all(state%fields(node)%values==saved%fields(node)%values),'state rollback')
        call require(all(temp%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
      end do
      call require(all(budget==-17).and.maximum==-19.and.all(counts==-23).and.theta==-29,'diagnostic rollback')
    end do
  end do
  print *, 'PASS: recursive transport, 1/4/36 steps, root oracle, physical/periodic budgets and rollback'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: '//label
      error stop 1
    end if
  end subroutine
end program
