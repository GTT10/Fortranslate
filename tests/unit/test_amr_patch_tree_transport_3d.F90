program test_amr_patch_tree_transport_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_transport_3d_mod, only: advance_reactive_transport_3d,reactive_transport_timestep_3d, &
    reactive_isothermal_timestep_3d,reactive_transport_euler_update_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_transport_3d_mod, only: advance_amr_patch_tree_transport_3d,compute_amr_patch_tree_transport_timestep_3d
  use amr_patch_tree_transport_3d_mod, only: compute_sdc_amr_patch_tree_transport_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,reflux_amr_patch_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,temperature,saved,saved_t
  type(amr_patch_3d) :: patches(5)
  real(dp) :: root(12,4,4,4),rt(1,4,4,4),serial(12,4,4,4),st(4,4,4),q(12),sound
  real(dp) :: spacing(3),dt,root_dt,expected_dt,local_dt,h(3),budget(12),serial_budget(12),largest, &
      initial(12),final(12),theta,serial_theta,diffusivity,local_diffusivity
  real(dp) :: wall_t(6),wall_v(3,6),wall_dt,outlet(6),mach(6),tin(6),pin(6)
  integer :: i,j,k,mode,geometry,node,n(3),axis,level,boundary_mode
  character(len=:), allocatable :: coarse_limiter
  type(amr_tree_state_3d) :: pcm_results(4)
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
    do mode=1,12
      outlet=0;mach=-1;tin=700;pin=0
      if (mode>8) then
        outlet(2)=102325;outlet(6)=102325
      end if
      boundary_mode=mod(mode-1,4)+1
      if (mode==9) boundary_mode=2
      if (mode==10) boundary_mode=4
      if (allocated(coarse_limiter)) deallocate(coarse_limiter)
      if (mode>4) coarse_limiter="mc"
      faces='periodic';wall_v=0
      if (boundary_mode==2) faces='outflow'
      if (boundary_mode>2) then
        wall_v(1,3:4)=2
        faces='outflow';faces(3:4)='no_slip_wall'
        if (boundary_mode==4) faces(3:4)='isothermal_wall'
      end if

      if(mode>10) then
        outlet=0;faces(1)='inflow';mach(1)=.2_dp
        if(mode==12) pin(1)=120000
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
      if (boundary_mode==4) then
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
      if(mode==1) call check_sdc_transport_stage()
      call advance_amr_patch_tree_transport_3d(tree,species,transport,state,temperature,spacing,dt,faces, &
        .true.,.true.,.true.,.true.,budget,largest,theta,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
        coarse_limiter=coarse_limiter,outlet_pressure=outlet, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
      call require(ok,'transport step')
      if (geometry==2.and.mode<=8) then
        if (mode<=4) then
          pcm_results(mode)=state
        else
          call require(maxval(abs(state%fields(2)%values-pcm_results(boundary_mode)%fields(2)%values))>1e-12_dp, &
            'linear donors change coarse/fine transport')
        end if
      end if
      if (geometry==1) then
        serial=root;st=rt(1,:,:,:)
        call advance_reactive_transport_3d(species,transport,serial,st,4,4,4,spacing(1),spacing(2),spacing(3),dt, &
          .true.,.true.,.true.,.true.,serial_theta,ok,faces,serial_budget,wall_temperature=wall_t,wall_velocity=wall_v, &
          outlet_pressure=outlet, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
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
        .true.,.true.,.true.,.true.,budget,largest,theta,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
        coarse_limiter=coarse_limiter,outlet_pressure=outlet, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
      call require(.not.ok,'nonphysical candidate rejection')
      do node=1,size(tree%nodes)
        call require(all(state%fields(node)%values==saved%fields(node)%values),'state rollback')
        call require(all(temperature%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
      end do
      call require(all(budget==-17).and.largest==-19.and.theta==-23,'diagnostics rollback')
    end do
    call check_thermal_slip_stage()
  end do
contains
  subroutine check_thermal_slip_stage()
    type(amr_tree_state_3d) :: reference,actual,retained,sample,sample_t
    type(amr_tree_flux_3d), allocatable :: reference_flux(:),actual_flux(:)
    character(len=32) :: bc(6)
    logical :: mask(6)
    real(dp) :: minimum,reference_dt,actual_dt,maximum,wall_values(6)
    real(dp) :: heat_q(12),heat_u(12),heat_t,heat_sound
    integer :: face,peer,a,b,cidx,dims(3)
    ! With zero velocity the slip/no-slip ghost states coincide. This oracle
    ! isolates thermal face placement; nonuniform moving ghosts need not agree.
    sample=state;sample_t=temperature
    heat_q=0;heat_q(1)=.3_dp;heat_q(5)=101325
    heat_q(6)=.03_dp;heat_q(9)=.22_dp;heat_q(12)=.75_dp
    call reactive_primitive_to_conserved(species,heat_q,heat_u,heat_t,heat_sound,ok)
    call require(ok,'stationary thermal face oracle')
    do peer=1,size(tree%nodes)
      dims=tree%nodes(peer)%shape
      do cidx=1,dims(3)
      do b=1,dims(2)
      do a=1,dims(1)
        sample%fields(peer)%values(:,a,b,cidx)=heat_u
        sample_t%fields(peer)%values(1,a,b,cidx)=heat_t
      enddo
      enddo
      enddo
    enddo
    wall_values=1100
    do face=1,6
      bc='slip_wall';bc(face)='isothermal_wall'
      minimum=1
      call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,sample,sample_t,spacing,1e-8_dp,bc, &
        [.false.,.true.,.false.,.false.],reference,reference_flux,minimum,ok,wall_temperature=wall_values)
      call require(ok,'reference AMR wall heat')
      call compute_amr_patch_tree_transport_timestep_3d(tree,species,transport,sample,sample_t,spacing, &
        .15_dp,.false.,.true.,.false.,reference_dt,maximum,ok,bc,wall_values)
      call require(ok,'reference AMR wall timestep')
      bc='slip_wall';mask=.false.;mask(face)=.true.
      call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,sample,sample_t,spacing,1e-8_dp,bc, &
        [.false.,.true.,.false.,.false.],actual,actual_flux,minimum,ok,wall_temperature=wall_values, &
        isothermal_slip_faces=mask)
      call require(ok,'explicit AMR thermal slip')
      do peer=1,size(tree%nodes)
        call require(all(actual%fields(peer)%values==reference%fields(peer)%values),'thermal slip rates exact')
        call require(all(actual_flux(peer)%x==reference_flux(peer)%x).and. &
          all(actual_flux(peer)%y==reference_flux(peer)%y).and. &
          all(actual_flux(peer)%z==reference_flux(peer)%z),'physical-face thermal slip flux exact')
      enddo
      call compute_amr_patch_tree_transport_timestep_3d(tree,species,transport,sample,sample_t,spacing, &
        .15_dp,.false.,.true.,.false.,actual_dt,maximum,ok,bc,wall_values,isothermal_slip_faces=mask)
      call require(ok.and.actual_dt==reference_dt,'thermal slip AMR timestep exact')
      retained=actual;minimum=-23;bc(face)='outflow'
      call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,sample,sample_t,spacing,1e-8_dp,bc, &
        [.false.,.true.,.false.,.false.],actual,actual_flux,minimum,ok,wall_temperature=wall_values, &
        isothermal_slip_faces=mask)
      call require(.not.ok.and.minimum==-23,'invalid thermal slip mask rejects transaction')
      do peer=1,size(tree%nodes)
        call require(all(actual%fields(peer)%values==retained%fields(peer)%values),'thermal slip failed rate rollback')
        call require(all(actual_flux(peer)%x==reference_flux(peer)%x).and. &
          all(actual_flux(peer)%y==reference_flux(peer)%y).and. &
          all(actual_flux(peer)%z==reference_flux(peer)%z),'thermal slip failed flux rollback')
      enddo
    enddo
    print *, 'PASS: six-face AMR thermal slip, physical-face clipping, timestep and rollback'
  end subroutine check_thermal_slip_stage
  subroutine check_sdc_transport_stage()
    type(amr_tree_state_3d) :: rates,candidate,prior,selected,packed,future,retained,held_probe,total_rates,held_rates
    type(amr_tree_flux_3d), allocatable :: flux(:),prior_flux(:),selected_flux(:)
    real(dp) :: before(12),after(12),correction(12),largest,minimum,expected(12,4,4,4),et(4,4,4),unused
    integer :: peer,depth,dims(3)
    real(dp), allocatable :: fractions(:)
    minimum=-1
    call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,state,temperature,spacing,dt,faces, &
      [.true.,.true.,.true.,.true.],rates,flux,minimum,ok)
    call require(ok.and.minimum>=0.and.minimum<=1,'SDC hierarchy transport rate and integrated flux')
    if(geometry==1) then
      held_probe=temperature
      held_probe%fields(1)%values(1,:,:,:)=1.001_dp*state%fields(1)%values(5,:,:,:)
      call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,state,temperature,spacing,dt,faces, &
        [.true.,.true.,.true.,.true.],total_rates,selected_flux,minimum,ok, &
        internal_energy=held_probe,total_energy_primitives=.true.)
      call require(ok,'MOL transport with distinct retained internal energy')
      call require(all(total_rates%fields(1)%values==rates%fields(1)%values), &
        'MOL transport uses total energy despite distinct retained internal energy')
      call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,state,temperature,spacing,dt,faces, &
        [.true.,.true.,.true.,.true.],held_rates,selected_flux,minimum,ok,internal_energy=held_probe)
      call require(ok,'default held transport retained')
      call require(maxval(abs(held_rates%fields(1)%values-rates%fields(1)%values))>1e-8_dp, &
        'held transport control distinguishes thermodynamic input')
    endif
    candidate=state
    do peer=1,size(tree%nodes)
      candidate%fields(peer)%values=state%fields(peer)%values+dt*rates%fields(peer)%values
    end do
    if(geometry==1) then
      call reactive_transport_euler_update_3d(species,transport,root,rt(1,:,:,:),4,4,4, &
        spacing(1),spacing(2),spacing(3),dt,.true.,.true.,.true.,.true.,expected,et,unused,ok)
      call require(ok,'SDC transport independent single-grid update')
      call require(maxval(abs(candidate%fields(1)%values-expected)/max(1._dp,abs(expected)))<1e-12_dp, &
        'SDC root transport matches existing Euler update')
    end if
    call composite_integrals_amr_patch_tree_3d(tree,state,spacing,before,ok)
    call require(ok,'SDC transport initial inventory')
    call reflux_amr_patch_tree_3d(tree,candidate,flux,spacing,faces,correction,largest,ok)
    call require(ok,'SDC transport flux-register compatibility')
    call composite_integrals_amr_patch_tree_3d(tree,candidate,spacing,after,ok)
    call require(ok.and.maxval(abs(after-before)/max(1e-12_dp,abs(before)))<1e-12_dp, &
      'SDC periodic composite transport conservation')
    prior=rates;prior_flux=flux;unused=minimum
    call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,state,temperature,spacing,-dt,faces, &
      [.true.,.true.,.true.,.true.],rates,flux,minimum,ok)
    call require(.not.ok.and.minimum==unused,'SDC transport invalid step preserves theta')
    do peer=1,size(tree%nodes)
      call require(all(rates%fields(peer)%values==prior%fields(peer)%values).and. &
        all(flux(peer)%x==prior_flux(peer)%x).and.all(flux(peer)%y==prior_flux(peer)%y).and. &
        all(flux(peer)%z==prior_flux(peer)%z),'SDC transport stage output rollback')
    end do
    do depth=0,tree%level_count()-1
      call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,state,temperature,spacing,dt,faces, &
        [.true.,.true.,.true.,.true.],selected,selected_flux,minimum,ok,target_depth=depth)
      call require(ok,'level-selective transport')
      do peer=1,size(tree%nodes)
        if(tree%nodes(peer)%level==depth) then
          call require(all(selected%fields(peer)%values==rates%fields(peer)%values),'selected transport rate exact')
          call require(all(selected_flux(peer)%x==flux(peer)%x),'selected transport flux exact')
        else
          call require(all(selected%fields(peer)%values==0),'inactive transport rate zero')
          call require(.not.allocated(selected_flux(peer)%x),'inactive transport flux absent')
        end if
      end do
    end do
    depth=tree%level_count()-1
    if(depth==0) return
    allocate(packed%fields(size(tree%nodes)),fractions(depth))
    do peer=1,size(tree%nodes)
      dims=tree%nodes(peer)%shape
      allocate(packed%fields(peer)%values(13,dims(1),dims(2),dims(3)))
      packed%fields(peer)%values(:12,:,:,:)=state%fields(peer)%values
      packed%fields(peer)%values(13,:,:,:)=temperature%fields(peer)%values(1,:,:,:)
    end do
    future=packed
    do peer=1,size(tree%nodes)
      if(tree%nodes(peer)%level>=depth) cycle
      future%fields(peer)%values(13,:,:,:)=future%fields(peer)%values(13,:,:,:)+1
    end do
    fractions=0
    call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,state,temperature,spacing,dt,faces, &
      [.true.,.true.,.true.,.true.],selected,selected_flux,minimum,ok,target_depth=depth, &
      coarse_start=packed,coarse_end=future,coarse_fractions=fractions)
    call require(ok,'transport old coarse endpoint')
    do peer=1,size(tree%nodes)
      if(tree%nodes(peer)%level/=depth) cycle
      call require(all(selected%fields(peer)%values==rates%fields(peer)%values),'transport old endpoint exact')
    end do
    retained=selected;fractions=.5_dp
    call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,state,temperature,spacing,dt,faces, &
      [.true.,.true.,.true.,.true.],selected,selected_flux,minimum,ok,target_depth=depth, &
      coarse_start=packed,coarse_end=future,coarse_fractions=fractions)
    call require(ok,'transport intermediate coarse endpoint')
    unused=0
    do peer=1,size(tree%nodes)
      unused=max(unused,maxval(abs(selected%fields(peer)%values-retained%fields(peer)%values)))
    end do
    call require(unused>1e-10_dp,'parent temperature time changes transport')
  end subroutine check_sdc_transport_stage
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
