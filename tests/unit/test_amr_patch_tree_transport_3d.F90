program test_amr_patch_tree_transport_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_transport_3d_mod, only: advance_reactive_transport_3d,reactive_transport_timestep_3d, &
    reactive_isothermal_timestep_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_transport_3d_mod, only: advance_amr_patch_tree_transport_3d,compute_amr_patch_tree_transport_timestep_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,temperature,saved,saved_t
  type(amr_patch_3d) :: patches(5)
  real(dp) :: root(12,4,4,4),rt(1,4,4,4),serial(12,4,4,4),st(4,4,4),q(12),sound
  real(dp) :: spacing(3),dt,root_dt,expected_dt,local_dt,h(3),budget(12),serial_budget(12),largest, &
      initial(12),final(12),theta,serial_theta,diffusivity,local_diffusivity
  real(dp) :: wall_t(6),wall_v(3,6),wall_dt
  integer :: i,j,k,mode,geometry,node,n(3),axis,level
  integer(int64), allocatable :: origin(:,:),domain(:,:)
  integer, allocatable :: ratios(:)
  logical :: ok
  character(len=32) :: faces(6),local_faces(6)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  wall_t=1100;wall_v=0;wall_v(1,3:4)=2
  spacing=[.001_dp,.0015_dp,.002_dp]
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
  do geometry=1,2
    if (geometry==1) then
      call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches(:0),tree,ok)
    else
      call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3,4],patches,tree,ok)
    end if
    call require(ok,'topology')
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    call require(ok,'global layout')
    do mode=1,4
      faces='periodic';wall_v=0
      if (mode==2) faces='outflow'
      if (mode>2) then
        wall_v(1,3:4)=2
        faces='outflow';faces(3:4)='no_slip_wall'
        if (mode==4) faces(3:4)='isothermal_wall'
      end if

      call prolong_amr_patch_tree_3d(tree,root,state,ok)
      call require(ok,'state initialization')
      call prolong_amr_patch_tree_3d(tree,rt,temperature,ok)
      call require(ok,'temperature initialization')
      call compute_amr_patch_tree_transport_timestep_3d(tree,species,transport,state,temperature,spacing, &
      .15_dp,.true.,.true.,.true.,dt,diffusivity,ok)
      call require(ok.and.dt>0,'CFL')
      if (geometry==1) root_dt=dt
      expected_dt=huge(1.0_dp)
      do node=1,size(tree%nodes)
        n=tree%nodes(node)%shape;h=spacing/2**tree%nodes(node)%level
        call reactive_transport_timestep_3d(species,transport,state%fields(node)%values, &
          temperature%fields(node)%values(1,:,:,:),n(1),n(2),n(3),h(1),h(2),h(3),.15_dp,.true.,.true.,.true., &
      local_dt,local_diffusivity,ok)
        call require(ok,'local CFL oracle')
        expected_dt=min(expected_dt,local_dt)
      end do
      call require(dt==expected_dt,'hierarchy CFL minimum')
      if (geometry==2) call require(dt<root_dt/16,'finest grid CFL controls step')
      if (mode==4) then
        do node=1,size(tree%nodes)
          n=tree%nodes(node)%shape;h=spacing/2**tree%nodes(node)%level
          local_faces='outflow';level=tree%nodes(node)%level+1
          do axis=1,3
            if (origin(axis,node)==0.and.faces(2*axis-1)=='isothermal_wall') &
              local_faces(2*axis-1)='isothermal_wall'
            if (origin(axis,node)+int(n(axis),int64)==domain(axis,level).and.faces(2*axis)=='isothermal_wall') &
              local_faces(2*axis)='isothermal_wall'
          end do
          call reactive_isothermal_timestep_3d(species,transport,state%fields(node)%values, &
            temperature%fields(node)%values(1,:,:,:),h,local_faces,wall_t,.15_dp,wall_dt,ok)
          call require(ok,'local wall timestep')
          expected_dt=min(expected_dt,wall_dt)
        end do
        call compute_amr_patch_tree_transport_timestep_3d(tree,species,transport,state,temperature, &
          spacing,.15_dp,.true.,.true.,.true.,dt,diffusivity,ok,faces,wall_t)
        call require(ok.and.dt==expected_dt,'wall-aware timestep')
      end if
      call composite_integrals_amr_patch_tree_3d(tree,state,spacing,initial,ok)
      call require(ok,'initial inventory')
      dt=min(dt,1e-8_dp)
      call advance_amr_patch_tree_transport_3d(tree,species,transport,state,temperature,spacing,dt,faces, &
        .true.,.true.,.true.,.true.,budget,largest,theta,ok,wall_temperature=wall_t,wall_velocity=wall_v)
      call require(ok,'hydro step')
      if (geometry==1) then
        serial=root;st=rt(1,:,:,:)
        call advance_reactive_transport_3d(species,transport,serial,st,4,4,4,spacing(1),spacing(2),spacing(3),dt, &
          .true.,.true.,.true.,.true.,serial_theta,ok,faces,serial_budget,wall_temperature=wall_t,wall_velocity=wall_v)
        call require(ok,'serial oracle')
        call require(maxval(abs(state%fields(1)%values-serial)/max(1.0_dp,abs(serial)))<2e-13_dp,'serial state')
        call require(maxval(abs(temperature%fields(1)%values(1,:,:,:)-st))<1e-9_dp,'serial temperature')
        call require(maxval(abs(budget-serial_budget))<1e-14_dp,'serial boundary ledger')
      end if
      call composite_integrals_amr_patch_tree_3d(tree,state,spacing,final,ok)
      call require(ok,'final inventory')
      call require(maxval(abs(final-initial-budget)/max(1e-12_dp,abs(initial)))<2e-12_dp,'composite conservation')
      do node=1,size(tree%nodes)
        call require(all(temperature%fields(node)%values>0),'positive temperature')
        call require(maxval(abs(sum(state%fields(node)%values(6:,:,:,:),dim=1)- &
          state%fields(node)%values(1,:,:,:)))<1e-12_dp,'species closure')
      end do
      saved=state;saved_t=temperature;budget=-17;largest=-19;theta=-23
      call advance_amr_patch_tree_transport_3d(tree,species,transport,state,temperature,spacing,1.0_dp,faces, &
        .true.,.true.,.true.,.true.,budget,largest,theta,ok,wall_temperature=wall_t,wall_velocity=wall_v)
      call require(.not.ok,'nonphysical candidate rejection')
      do node=1,size(tree%nodes)
        call require(all(state%fields(node)%values==saved%fields(node)%values),'state rollback')
        call require(all(temperature%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
      end do
      call require(all(budget==-17).and.largest==-19.and.theta==-23,'diagnostics rollback')
    end do
  end do
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
