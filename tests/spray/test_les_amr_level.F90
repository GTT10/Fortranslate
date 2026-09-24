program test_les_amr_level
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d, &
    initialize_amr_patch_tree_3d,prolong_amr_patch_tree_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use les_amr_patch_tree_mod, only: advance_amr_patch_tree_les_level_3d,advance_amr_patch_tree_les_3d, &
    compute_amr_patch_tree_les_timestep_3d
  use les_amr_subcycle_mod, only: compute_amr_patch_tree_subcycle_les_timestep_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patch(1)
  type(amr_tree_state_3d) :: state,temps,before,before_t,oracle,oracle_t,start,finish
  type(amr_tree_flux_3d), allocatable :: flux(:),saved_flux(:),control_flux(:)
  real(dp) :: root(12,4,4,4),temperature(1,4,4,4),q(12),sound,h(3),dt,diffusion,budget(12),correction
  real(dp) :: common_dt,sub_dt,half_dt,common_diff,sub_diff,half_diff
  real(dp) :: change(12,4,4,4),fractions0(1),fractions1(1)
  integer :: geometry,m,i,j,k,node,active
  character(len=32) :: faces(6),model
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  do k=1,4
    do j=1,4
      do i=1,4
        q=0;q(1)=.3_dp;q(2)=2*j;q(3)=3*i;q(5)=101325;q(12)=1
        call reactive_primitive_to_conserved(species,q,root(:,i,j,k),temperature(1,i,j,k),sound,ok)
        call require(ok,'initial gas')
      end do
    end do
  end do
  call initialize_amr_patch_3d(4,4,4,2,3,2,3,2,3,2,patch(1),ok)
  call require(ok,'child geometry')
  h=.001_dp;dt=1e-9_dp;faces='periodic';fractions0=.25_dp;fractions1=.75_dp
  do geometry=0,1
    if(geometry==0) then
      call initialize_amr_patch_tree_3d([4,4,4],[integer::],patch(:0),tree,ok)
    else
      call initialize_amr_patch_tree_3d([4,4,4],[1],patch,tree,ok)
    end if
    call require(ok,'tree')
    call prolong_amr_patch_tree_3d(tree,root,before,ok)
    call require(ok,'state')
    call prolong_amr_patch_tree_3d(tree,temperature,before_t,ok)
    call require(ok,'temperatures')
    if(allocated(start%fields)) deallocate(start%fields,finish%fields)
    allocate(start%fields(size(tree%nodes)),finish%fields(size(tree%nodes)))
    allocate(start%fields(1)%values(13,4,4,4))
    start%fields(1)%values(:12,:,:,:)=root;start%fields(1)%values(13,:,:,:)=temperature(1,:,:,:)
    finish%fields(1)=start%fields(1)
    finish%fields(1)%values(:12,:,:,:)=1.01_dp*root
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
      state=before;temps=before_t
      call compute_amr_patch_tree_les_timestep_3d(tree,species,state,temps,h,.05_dp,faces, &
        .17_dp,.7_dp,.7_dp,common_dt,common_diff,ok,model=model)
      call require(ok,'common interval')
      call compute_amr_patch_tree_subcycle_les_timestep_3d(tree,species,state,temps,h,.05_dp,faces, &
        .17_dp,.7_dp,.7_dp,sub_dt,sub_diff,ok,model=model)
      call require(ok.and.sub_dt>=common_dt.and.sub_dt<=4*common_dt,'subcycle interval bounds')
      call require(sub_diff==common_diff,'subcycling preserves diffusivity diagnostic')
      if(geometry==0) call require(sub_dt==common_dt,'root interval equality')
      if(geometry==1.and.m==4) call require(sub_dt>common_dt,'dynamic fine-grid bound is relaxed by subcycling')
      call compute_amr_patch_tree_subcycle_les_timestep_3d(tree,species,state,temps,h,.025_dp,faces, &
        .17_dp,.7_dp,.7_dp,half_dt,half_diff,ok,model=model)
      call require(ok.and.half_dt==.5_dp*sub_dt.and.half_diff==sub_diff,'subcycle CFL scaling')
      half_dt=-7;half_diff=-9
      call compute_amr_patch_tree_subcycle_les_timestep_3d(tree,species,state,temps,h,.3_dp,faces, &
        .17_dp,.7_dp,.7_dp,half_dt,half_diff,ok,model=model)
      call require(.not.ok.and.half_dt==-7.and.half_diff==-9,'invalid subcycle CFL preserves outputs')

      if(geometry==0) then
        call advance_amr_patch_tree_les_level_3d(tree,species,state,temps,0,h,dt,faces, &
          .17_dp,.7_dp,.7_dp,flux,diffusion,ok,model=model)
      else
        call advance_amr_patch_tree_les_level_3d(tree,species,state,temps,1,h,dt,faces, &
          .17_dp,.7_dp,.7_dp,flux,diffusion,ok,start,finish,fractions0,fractions1,model=model,coarse_limiter='mc')
      end if
      call require(ok.and.diffusion>=0,'level SSPRK2')
      active=geometry+1
      call require(any(state%fields(active)%values/=before%fields(active)%values),'active level evolves')
      do node=1,size(tree%nodes)
        call require(allocated(flux(node)%x).eqv.(node==active),'only active integrated fluxes')
        if(node==active) cycle
        call require(all(state%fields(node)%values==before%fields(node)%values), 'inactive state preserved')
        call require(all(temps%fields(node)%values==before_t%fields(node)%values), 'inactive temperature preserved')
      end do
      change=(flux(active)%x(:,:4,:,:)-flux(active)%x(:,2:,:,:))/(h(1)/2**geometry)+ &
        (flux(active)%y(:,:,:4,:)-flux(active)%y(:,:,2:,:))/(h(2)/2**geometry)+ &
        (flux(active)%z(:,:,:,:4)-flux(active)%z(:,:,:,2:))/(h(3)/2**geometry)
      call require(maxval(abs(state%fields(active)%values-before%fields(active)%values-change))<= &
        2e-12_dp*maxval(abs(before%fields(active)%values)),'integrated flux accounts for accepted state')
      if(geometry==0) then
        oracle=before;oracle_t=before_t
        call advance_amr_patch_tree_les_3d(tree,species,oracle,oracle_t,h,dt,faces, &
          .17_dp,.7_dp,.7_dp,budget,correction,ok,model=model)
        call require(ok.and.all(oracle%fields(1)%values==state%fields(1)%values), 'root common-step state parity')
        call require(maxval(abs(oracle_t%fields(1)%values-temps%fields(1)%values))<1e-8_dp,'root temperature parity')
      else
        oracle=before;oracle_t=before_t
        call advance_amr_patch_tree_les_level_3d(tree,species,oracle,oracle_t,1,h,dt,faces, &
          .17_dp,.7_dp,.7_dp,control_flux,diffusion,ok,start,start,fractions0,fractions1, &
          model=model,coarse_limiter='mc')
        call require(ok.and.any(oracle%fields(2)%values/=state%fields(2)%values), &
          'time-varying coarse donors affect the active level')
        oracle=state;oracle_t=temps;saved_flux=flux;diffusion=-7
        call advance_amr_patch_tree_les_level_3d(tree,species,state,temps,1,h,dt,faces, &
          .17_dp,.7_dp,.7_dp,flux,diffusion,ok,model=model)
        call require(.not.ok.and.diffusion==-7,'missing coarse endpoints rejected')
        call require(all(state%fields(2)%values==oracle%fields(2)%values).and. &
          all(temps%fields(2)%values==oracle_t%fields(2)%values),'failed level step rollback')
        call require(all(flux(2)%x==saved_flux(2)%x).and.all(flux(2)%y==saved_flux(2)%y).and. &
          all(flux(2)%z==saved_flux(2)%z),'failed level step preserves flux output')
        call advance_amr_patch_tree_les_level_3d(tree,species,state,temps,1,h,1._dp,faces, &
          .17_dp,.7_dp,.7_dp,flux,diffusion,ok,start,finish,fractions0,fractions1,model=model,coarse_limiter='mc')
        call require(.not.ok.and.diffusion==-7,'excessive level interval rejected')
        call require(all(state%fields(2)%values==oracle%fields(2)%values).and. &
          all(temps%fields(2)%values==oracle_t%fields(2)%values),'stage failure preserves gas state')
        call require(all(flux(2)%x==saved_flux(2)%x).and.all(flux(2)%y==saved_flux(2)%y).and. &
          all(flux(2)%z==saved_flux(2)%z),'stage failure preserves integrated flux')
      end if
    end do
  end do
  print '(a)','test_les_amr_level: PASS'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(condition) return
    print '(a)','FAIL: '//label
    error stop 1
  end subroutine
end program
