program test_les_amr_subcycle
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,average_down_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use les_amr_patch_tree_mod, only: advance_amr_patch_tree_les_3d,advance_amr_patch_tree_les_level_3d, &
    convert_amr_les_tree_energy_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use les_amr_subcycle_mod, only: advance_amr_patch_tree_subcycle_les_3d, &
    compute_amr_patch_tree_subcycle_les_timestep_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(3)
  type(amr_tree_state_3d) :: state,temp,reference,reference_t,saved,saved_t,synchronized
  type(amr_tree_state_3d) :: coarse_plain,coarse_rotated
  type(rotating_frame_3d) :: frame
  type(amr_tree_flux_3d), allocatable :: plain_flux(:),rotated_flux(:)
  real(dp) :: root(12,4,4,4),rt(1,4,4,4),q(12),sound,initial(12),final(12),budget(12),maximum
  real(dp) :: reference_budget(12),reference_max,dt,diffusivity,wall_t(6),wall_v(3,6)
  real(dp), parameter :: spacing(3)=[.001_dp,.0015_dp,.002_dp]
  real(dp) :: outlet_pressure(6),mach(6),tin(6),pin(6),selected_dt,selected_d
  integer(int64), allocatable :: counts(:)
  integer :: geometry,mode,node,i,j,k,m
  character(len=32) :: faces(6),model
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  wall_t=1100;wall_v=0
  do k=1,4
    do j=1,4
      do i=1,4
        q=0;q(1)=.3_dp+.001_dp*(i+j+k);q(2:4)=[2._dp*j,-1._dp*i,3._dp*k];q(5)=101325+100*j
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
  do m=1,4
    select case(m)
    case(1)
      model="smagorinsky"
    case(2)
      model="wale"
    case(3)
      model="vreman"
    case(4)
      model="dynamic"
    end select
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
      faces='periodic';wall_v=0
      if (mode>=2) then
        faces='outflow';faces(3:4)='isothermal_wall'
        wall_v(1,3:4)=2
      end if
      outlet_pressure=0;mach=-1;tin=700;pin=0
      if (mode==3) outlet_pressure([2,6])=102325
      ! Avoid repeated seven-layer corner extrapolation producing negative pressure.
      if (mode==3.and.m==4) outlet_pressure(6)=0
      if(mode>=4) then
        faces(1)='inflow';mach(1)=.2_dp
        if(mode==5) pin(1)=120000
      end if
      call prolong_amr_patch_tree_3d(tree,root,state,ok)
      call require(ok,'initial state')
      call prolong_amr_patch_tree_3d(tree,rt,temp,ok)
      call require(ok,'initial temperature')
      call composite_integrals_amr_patch_tree_3d(tree,state,spacing,initial,ok)
      call require(ok,'initial inventory')
      if(mode>=4) then
        call compute_amr_patch_tree_subcycle_les_timestep_3d(tree,species,state,temp,spacing,.05_dp,faces, &
          .17_dp,.7_dp,.7_dp,selected_dt,selected_d,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
          coarse_limiter='mc',model=model,inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
        call require(ok.and.selected_dt>0.and.selected_d>=0,'Mach subcycle timestep')
      end if
      dt=1e-9_dp
      reference=state;reference_t=temp
      if (geometry==1) then
        call advance_amr_patch_tree_les_3d(tree,species,reference,reference_t,spacing,dt,faces, &
          .17_dp,.7_dp,.7_dp,reference_budget,reference_max,ok, &
          wall_temperature=wall_t,wall_velocity=wall_v,coarse_limiter='mc',outlet_pressure=outlet_pressure,model=model, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
        if(.not.ok) print *, 'root oracle failure: ',trim(model),mode
        call require(ok,'root oracle')
      end if
      call advance_amr_patch_tree_subcycle_les_3d(tree,species,state,temp,spacing,dt,faces, &
        .17_dp,.7_dp,.7_dp,budget,maximum,diffusivity,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
        coarse_limiter='mc',level_steps=counts,outlet_pressure=outlet_pressure,model=model, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
      call require(ok,'recursive LES step')
      call require(diffusivity>=0,'LES diffusivity')
      if (geometry==1) then
        call require(all(counts==1),'root step count')
        call require(maxval(abs(state%fields(1)%values-reference%fields(1)%values))<2e-11_dp,'root state equality')
        if(maxval(abs(temp%fields(1)%values-reference_t%fields(1)%values))>=1e-9_dp) &
          print *, 'root temperature difference: ',trim(model),mode,maxval(abs(temp%fields(1)%values-reference_t%fields(1)%values))
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
      saved=state;saved_t=temp;budget=-17;maximum=-19;counts=-23;diffusivity=-29
      call advance_amr_patch_tree_subcycle_les_3d(tree,species,state,temp,spacing,1._dp,faces, &
        .17_dp,.7_dp,.7_dp,budget,maximum,diffusivity,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
        coarse_limiter='mc',level_steps=counts,outlet_pressure=outlet_pressure,model=model, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
      call require(.not.ok,'nonphysical step rejected')
      do node=1,size(tree%nodes)
        call require(all(state%fields(node)%values==saved%fields(node)%values),'state rollback')
        call require(all(temp%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
      end do
      call require(all(budget==-17).and.maximum==-19.and.all(counts==-23).and.diffusivity==-29,'diagnostic rollback')
  end do
  end do
  end do
  frame%axis=3;frame%angular_speed=400._dp;frame%axis_point=[.001_dp,-.002_dp,0._dp]
  faces='outflow'
  call prolong_amr_patch_tree_3d(tree,root,reference,ok)
  call require(ok,'rotating recursive LES initial state')
  call prolong_amr_patch_tree_3d(tree,rt,reference_t,ok)
  call require(ok,'rotating recursive LES initial temperature')
  state=reference;temp=reference_t
  call convert_amr_les_tree_energy_3d(tree,state,spacing,frame,[0._dp,0._dp,0._dp],.true.,ok)
  call require(ok,'rotating recursive LES energy storage')
  selected_dt=-1;selected_d=-1
  call compute_amr_patch_tree_subcycle_les_timestep_3d(tree,species,reference,reference_t,spacing,.05_dp,faces, &
    .17_dp,.7_dp,.7_dp,selected_dt,selected_d,ok)
  call require(ok.and.selected_dt>0,'recursive LES timestep oracle')
  dt=-1;diffusivity=-1
  call compute_amr_patch_tree_subcycle_les_timestep_3d(tree,species,state,temp,spacing,.05_dp,faces, &
    .17_dp,.7_dp,.7_dp,dt,diffusivity,ok,frame=frame,domain_lower=[0._dp,0._dp,0._dp])
  call require(ok.and.dt==selected_dt.and.diffusivity==selected_d,'rotating recursive LES timestep parity')
  dt=min(selected_dt,1e-9_dp)
  call advance_amr_patch_tree_subcycle_les_3d(tree,species,reference,reference_t,spacing,dt,faces, &
    .17_dp,.7_dp,.7_dp,reference_budget,reference_max,selected_d,ok)
  call require(ok,'recursive LES step oracle')
  call advance_amr_patch_tree_subcycle_les_3d(tree,species,state,temp,spacing,dt,faces, &
    .17_dp,.7_dp,.7_dp,budget,maximum,diffusivity,ok,frame=frame,domain_lower=[0._dp,0._dp,0._dp])
  call require(ok,'rotating recursive LES step')
  call convert_amr_les_tree_energy_3d(tree,state,spacing,frame,[0._dp,0._dp,0._dp],.false.,ok)
  call require(ok,'rotating recursive LES energy recovery')
  do node=1,size(tree%nodes)
    call require(maxval(abs(state%fields(node)%values-reference%fields(node)%values)/ &
      max(1._dp,abs(reference%fields(node)%values)))<2e-12_dp,'rotating recursive LES state parity')
  enddo
  call require(maxval(abs(budget-reference_budget))<1e-12_dp,'rotating recursive LES budget parity')
  call prolong_amr_patch_tree_3d(tree,root,reference,ok)
  call require(ok,'rotating LES level state')
  call prolong_amr_patch_tree_3d(tree,rt,reference_t,ok)
  call require(ok,'rotating LES level temperature')
  state=reference;temp=reference_t
  call convert_amr_les_tree_energy_3d(tree,state,spacing,frame,[0._dp,0._dp,0._dp],.true.,ok)
  call require(ok,'rotating LES level storage')
  call advance_amr_patch_tree_les_level_3d(tree,species,reference,reference_t,0,spacing,1e-9_dp,faces, &
    .17_dp,.7_dp,.7_dp,plain_flux,selected_d,ok)
  call require(ok,'ordinary LES level oracle')
  call advance_amr_patch_tree_les_level_3d(tree,species,state,temp,0,spacing,1e-9_dp,faces, &
    .17_dp,.7_dp,.7_dp,rotated_flux,diffusivity,ok,frame=frame,domain_lower=[0._dp,0._dp,0._dp])
  call require(ok,'rotating LES level step')
  call convert_amr_les_tree_energy_3d(tree,state,spacing,frame,[0._dp,0._dp,0._dp],.false.,ok)
  call require(ok,'rotating LES level recovery')
  call require(maxval(abs(state%fields(1)%values-reference%fields(1)%values)/ &
    max(1._dp,abs(reference%fields(1)%values)))<2e-12_dp,'rotating LES level state parity')
  call require(maxval(abs(rotated_flux(1)%x-plain_flux(1)%x))<1e-12_dp,'rotating LES level flux parity')
  call prolong_amr_patch_tree_3d(tree,root,reference,ok)
  call require(ok,'rotating LES fine level state')
  call prolong_amr_patch_tree_3d(tree,rt,reference_t,ok)
  call require(ok,'rotating LES fine level temperature')
  state=reference;temp=reference_t
  call convert_amr_les_tree_energy_3d(tree,state,spacing,frame,[0._dp,0._dp,0._dp],.true.,ok)
  call require(ok,'rotating LES fine level storage')
  allocate(coarse_plain%fields(size(tree%nodes)))
  allocate(coarse_plain%fields(1)%values(13,4,4,4))
  coarse_plain%fields(1)%values(:12,:,:,:)=reference%fields(1)%values
  coarse_plain%fields(1)%values(13,:,:,:)=reference_t%fields(1)%values(1,:,:,:)
  coarse_rotated=coarse_plain
  call convert_amr_les_tree_energy_3d(tree,coarse_rotated,spacing,frame,[0._dp,0._dp,0._dp],.true.,ok, &
    maximum_depth=0)
  call require(ok,'rotating LES coarse endpoint storage')
  call advance_amr_patch_tree_les_level_3d(tree,species,reference,reference_t,1,spacing,1e-9_dp,faces, &
    .17_dp,.7_dp,.7_dp,plain_flux,selected_d,ok,coarse_start=coarse_plain,coarse_end=coarse_plain, &
    fractions_start=[0._dp],fractions_end=[1._dp])
  call require(ok,'ordinary LES fine level oracle')
  call advance_amr_patch_tree_les_level_3d(tree,species,state,temp,1,spacing,1e-9_dp,faces, &
    .17_dp,.7_dp,.7_dp,rotated_flux,diffusivity,ok,coarse_start=coarse_rotated,coarse_end=coarse_rotated, &
    fractions_start=[0._dp],fractions_end=[1._dp],frame=frame,domain_lower=[0._dp,0._dp,0._dp])
  call require(ok,'rotating LES fine level step')
  call convert_amr_les_tree_energy_3d(tree,state,spacing,frame,[0._dp,0._dp,0._dp],.false.,ok)
  call require(ok,'rotating LES fine level recovery')
  do node=1,size(tree%nodes)
    if(tree%nodes(node)%level/=1) cycle
    call require(maxval(abs(state%fields(node)%values-reference%fields(node)%values)/ &
      max(1._dp,abs(reference%fields(node)%values)))<2e-12_dp,'rotating LES fine level state parity')
    call require(maxval(abs(rotated_flux(node)%x-plain_flux(node)%x))<1e-12_dp, &
      'rotating LES fine level flux parity')
  enddo
  print *, 'PASS: recursive LES, 1/4/36 steps, root oracle, physical/periodic budgets and rollback'
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
