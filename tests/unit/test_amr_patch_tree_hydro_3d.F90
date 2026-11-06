program test_amr_patch_tree_hydro_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: advance_reactive_boundary_hydro_3d,compute_reactive_cfl_timestep_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_hydro_3d_mod, only: advance_amr_patch_tree_hydro_3d,compute_amr_patch_tree_cfl_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_sdc_amr_patch_tree_hydro_3d,compute_mol_amr_patch_tree_hydro_3d, &
    convert_amr_patch_tree_rotating_energy_3d
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_sdc_3d,advance_amr_patch_tree_sdc_subcycle_3d
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use amr_patch_tree_ghost_3d_mod, only: fill_amr_patch_tree_ghosts_3d,apply_amr_tree_pressure_ghosts_3d
  use reactive_mol_3d_mod, only: compute_reactive_mol_halo_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,reflux_amr_patch_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,temperature,saved,saved_t
  type(amr_patch_3d) :: patches(5)
  real(dp) :: root(12,4,4,4),rt(1,4,4,4),serial(12,4,4,4),st(4,4,4),q(12),sound
  real(dp) :: spacing(3),dt,root_dt,expected_dt,local_dt,h(3),budget(12),serial_budget(12),largest,initial(12),final(12)
  real(dp) :: outlet(6),mach(6),tin(6),pin(6)
  integer :: i,j,k,mode,geometry,node,n(3)
  logical :: ok,ppm
  character(len=32) :: faces(6),reconstruction,arg
  call get_command_argument(1,arg)
  ppm=trim(arg)=='ppm'
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
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
    do mode=1,8
      if (ppm.and.(mode==2.or.mode==4.or.mode==6)) cycle
      outlet=0;mach=-1;tin=700;pin=0
      if (mode>4) then
        outlet(2)=102325;outlet(6)=102325
      end if
      faces='periodic'
      if (mode>2) then
        faces='outflow';faces(3:4)='slip_wall'
      end if
      reconstruction='pcm'
      if (mod(mode,2)==0) reconstruction='characteristic_plm'
      if(mode>6) then
        outlet=0;faces(1)='inflow';mach(1)=.2_dp
        if(mode==8) pin(1)=120000
        reconstruction='pcm'
        if(mode==8) reconstruction='characteristic_plm'
      end if
      if (ppm) reconstruction='characteristic_ppm'
      call prolong_amr_patch_tree_3d(tree,root,state,ok)
      call require(ok,'state initialization')
      call prolong_amr_patch_tree_3d(tree,rt,temperature,ok)
      call require(ok,'temperature initialization')
      if(trim(arg)=='rotation_energy') then
        if(mode==1) call check_rotation_energy()
        exit
      endif
      call compute_amr_patch_tree_cfl_3d(tree,species,state,temperature,spacing,.15_dp,dt,ok)
      call require(ok.and.dt>0,'CFL')
      if(trim(arg)=='rotating_mol') then
        if(mode==1) call check_rotation_mol_stage()
        exit
      endif
      if(trim(arg)=='rotating_step') then
        if(mode==1) call check_rotation_step()
        exit
      endif
      if(trim(arg)=='wall_mol') then
        if(geometry==1) call check_no_slip_mol()
        exit
      endif
      if(mode==1) call check_auxiliary_mol_stage()
      if(trim(arg)=='auxiliary') exit
      if(ppm.and.mode==1) call check_sdc_hydro_stage()
      if (geometry==1) root_dt=dt
      expected_dt=huge(1.0_dp)
      do node=1,size(tree%nodes)
        n=tree%nodes(node)%shape;h=spacing/2**tree%nodes(node)%level
        call compute_reactive_cfl_timestep_3d(species,state%fields(node)%values, &
          temperature%fields(node)%values(1,:,:,:),n(1),n(2),n(3),h(1),h(2),h(3),.15_dp,local_dt,ok)
        call require(ok,'local CFL oracle')
        expected_dt=min(expected_dt,local_dt)
      end do
      call require(dt==expected_dt,'hierarchy CFL minimum')
      if (geometry==2) call require(dt<root_dt/4,'finest grid CFL controls step')
      call composite_integrals_amr_patch_tree_3d(tree,state,spacing,initial,ok)
      call require(ok,'initial inventory')
      call advance_amr_patch_tree_hydro_3d(tree,species,state,temperature,spacing,dt,faces, &
        reconstruction,'mc','hllc',budget,largest,ok,outlet_pressure=outlet, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
      call require(ok,'hydro step')
      if (geometry==1) then
        serial=root;st=rt(1,:,:,:)
        call advance_reactive_boundary_hydro_3d(species,serial,st,spacing,dt,faces, &
          reconstruction,'mc','hllc',serial_budget,ok,outlet_pressure=outlet, &
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
      saved=state;saved_t=temperature;budget=-17;largest=-19
      call advance_amr_patch_tree_hydro_3d(tree,species,state,temperature,spacing,1.0_dp,faces, &
        reconstruction,'mc','hllc',budget,largest,ok,outlet_pressure=outlet, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
      call require(.not.ok,'nonphysical candidate rejection')
      do node=1,size(tree%nodes)
        call require(all(state%fields(node)%values==saved%fields(node)%values),'state rollback')
        call require(all(temperature%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
      end do
      call require(all(budget==-17).and.largest==-19,'diagnostics rollback')
    end do
  end do
contains
  subroutine check_rotation_step()
    type(rotating_frame_3d) :: frame
    type(amr_tree_state_3d) :: initial_state,initial_temperature
    type(elementary_reaction), allocatable :: reactions(:)
    type(gas_transport_species), allocatable :: transport(:)
    real(dp) :: lower(3),start(12),finish(12),boundary(12),source(12),largest,theta,residual
    integer :: variant
    frame%axis=3;frame%angular_speed=400._dp
    frame%axis_point=[.0004_dp,-.0007_dp,0._dp]
    lower=[-.002_dp,.003_dp,.001_dp]
    faces=['outflow ','outflow ','outflow ','outflow ','periodic','periodic']
    allocate(reactions(0),transport(0))
    call convert_amr_patch_tree_rotating_energy_3d(tree,state,lower,spacing,frame,.true.,ok)
    call require(ok,'rotating MOL step initial state')
    initial_state=state;initial_temperature=temperature
    call composite_integrals_amr_patch_tree_3d(tree,state,spacing,start,ok)
    call require(ok,'rotating MOL initial composite inventory')
    do variant=1,2
      state=initial_state;temperature=initial_temperature
      if(variant==1) then
        call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,state,temperature,spacing,dt*.001_dp,faces, &
          2,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          boundary,source,largest,theta,ok,mol_mode=.true.,mol_order=2,frame=frame,domain_lower=lower)
      else
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,state,temperature,spacing,dt*.001_dp, &
          faces,2,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          boundary,source,largest,theta,ok,mol_mode=.true.,mol_order=2,frame=frame,domain_lower=lower)
      endif
      if(.not.ok) print *, 'rotating MOL interval variant:',variant
      call require(ok,'rotating MOL interval')
      call composite_integrals_amr_patch_tree_3d(tree,state,spacing,finish,ok)
      call require(ok,'rotating MOL final composite inventory')
      residual=maxval(abs(finish-start-boundary-source)/max(1e-18_dp,abs(start)))
      call require(residual<3e-12_dp,'rotating MOL composite conservation')
    enddo
  end subroutine check_rotation_step

  subroutine check_rotation_mol_stage()
    type(rotating_frame_3d) :: frame,bad_frame
    type(amr_tree_state_3d) :: rate,saved_rate
    type(amr_tree_flux_3d), allocatable :: flux(:)
    real(dp) :: lower(3),h(3),balance,maximum
    integer :: peer,depth,x,y,z,dims(3)
    frame%axis=3;frame%angular_speed=400._dp
    frame%axis_point=[.0004_dp,-.0007_dp,0._dp]
    lower=[-.002_dp,.003_dp,.001_dp]
    faces='periodic'
    call convert_amr_patch_tree_rotating_energy_3d(tree,state,lower,spacing,frame,.true.,ok)
    call require(ok,'rotating MOL initial state')
    call compute_mol_amr_patch_tree_hydro_3d(tree,species,state,temperature,spacing,dt,faces,rate,flux,ok, &
      frame=frame,domain_lower=lower)
    call require(ok,'rotating MOL tree spatial rates')
    do peer=1,size(tree%nodes)
      dims=tree%nodes(peer)%shape;h=spacing
      do depth=1,tree%nodes(peer)%level
        h=h/real(tree%nodes(peer)%patch%refinement_ratio,dp)
      enddo
      call require(all(ieee_is_finite(rate%fields(peer)%values)), 'rotating MOL finite rate')
      maximum=0
      do z=1,dims(3)
      do y=1,dims(2)
      do x=1,dims(1)
        balance=-(flux(peer)%x(1,x+1,y,z)-flux(peer)%x(1,x,y,z))/(dt*h(1))- &
          (flux(peer)%y(1,x,y+1,z)-flux(peer)%y(1,x,y,z))/(dt*h(2))- &
          (flux(peer)%z(1,x,y,z+1)-flux(peer)%z(1,x,y,z))/(dt*h(3))
        maximum=max(maximum,abs(rate%fields(peer)%values(1,x,y,z)-balance))
      enddo
      enddo
      enddo
      call require(maximum<1e-10_dp,'rotating MOL mass flux divergence')
    enddo
    saved_rate=rate;bad_frame=frame;bad_frame%axis=0
    call compute_mol_amr_patch_tree_hydro_3d(tree,species,state,temperature,spacing,dt,faces,rate,flux,ok, &
      frame=bad_frame,domain_lower=lower)
    call require(.not.ok,'rotating MOL invalid frame rejection')
    do peer=1,size(tree%nodes)
      call require(all(rate%fields(peer)%values==saved_rate%fields(peer)%values), &
        'rotating MOL rejected rate rollback')
    enddo
  end subroutine check_rotation_mol_stage

  subroutine check_rotation_energy()
    type(rotating_frame_3d) :: frame
    type(amr_tree_state_3d) :: ordinary,rotated,bad
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: lower(3),h(3),position(3),potential,expected
    integer :: peer,depth,x,y,z
    frame%axis=3;frame%angular_speed=4000._dp
    frame%axis_point=[.0004_dp,-.0007_dp,0._dp]
    lower=[-.002_dp,.003_dp,.001_dp]
    ordinary=state
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    call require(ok,'rotation tree layout')
    call convert_amr_patch_tree_rotating_energy_3d(tree,state,lower,spacing,frame,.true.,ok)
    call require(ok,'rotation to stored energy')
    rotated=state
    do peer=1,size(tree%nodes)
      h=spacing
      do depth=2,tree%nodes(peer)%level+1
        h=h/real(ratios(depth),dp)
      enddo
      do z=1,tree%nodes(peer)%shape(3)
      do y=1,tree%nodes(peer)%shape(2)
      do x=1,tree%nodes(peer)%shape(1)
        position=lower+(real(origin(:,peer),dp)+real([x,y,z],dp)-.5_dp)*h
        call rotating_potential_3d(frame,position,potential,ok)
        call require(ok,'rotation potential')
        expected=ordinary%fields(peer)%values(5,x,y,z)-ordinary%fields(peer)%values(1,x,y,z)*potential
        call require(abs(state%fields(peer)%values(5,x,y,z)-expected)<2e-15_dp*abs(expected), &
          'rotation at patch cell position')
      enddo
      enddo
      enddo
    enddo
    call convert_amr_patch_tree_rotating_energy_3d(tree,state,lower,spacing,frame,.false.,ok)
    call require(ok,'rotation to ordinary energy')
    do peer=1,size(tree%nodes)
      call require(maxval(abs(state%fields(peer)%values(5,:,:,:)-ordinary%fields(peer)%values(5,:,:,:)))< &
        2e-15_dp*maxval(abs(ordinary%fields(peer)%values(5,:,:,:))), 'rotation round trip')
    enddo
    bad=rotated
    peer=size(tree%nodes)
    bad%fields(peer)%values(1,1,1,1)=0
    rotated=bad
    call convert_amr_patch_tree_rotating_energy_3d(tree,rotated,lower,spacing,frame,.false.,ok)
    call require(.not.ok,'rotation rejects zero density')
    do peer=1,size(tree%nodes)
      call require(all(rotated%fields(peer)%values==bad%fields(peer)%values),'rotation failure rollback')
    enddo
  end subroutine check_rotation_energy

  subroutine check_no_slip_mol()
    type(amr_tree_state_3d) :: u,packed,rates
    type(amr_tree_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: ghost(:,:,:,:),expected(:,:,:,:),slip_rate(:,:,:,:), &
      fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    real(dp) :: kinetic,rho,b(12)
    character(len=32) :: walls(6)
    integer :: x,y,z,n0(3),kind
    u=state;n0=tree%nodes(1)%shape
    do z=1,n0(3)
    do y=1,n0(2)
    do x=1,n0(1)
      rho=u%fields(1)%values(1,x,y,z)
      kinetic=.5_dp*sum(u%fields(1)%values(2:4,x,y,z)**2)/rho
      u%fields(1)%values(2:4,x,y,z)=rho*[1.5_dp,1._dp+.2_dp*x*x,0._dp]
      u%fields(1)%values(5,x,y,z)=u%fields(1)%values(5,x,y,z)-kinetic+ &
        .5_dp*sum(u%fields(1)%values(2:4,x,y,z)**2)/rho
    enddo
    enddo
    enddo
    packed=u;deallocate(packed%fields(1)%values)
    allocate(packed%fields(1)%values(13,n0(1),n0(2),n0(3)))
    packed%fields(1)%values(:12,:,:,:)=u%fields(1)%values
    packed%fields(1)%values(13,:,:,:)=temperature%fields(1)%values(1,:,:,:)
    allocate(expected(12,n0(1),n0(2),n0(3)),slip_rate(12,n0(1),n0(2),n0(3)), &
      fx(12,0:n0(1),n0(2),n0(3)),fy(12,n0(1),0:n0(2),n0(3)),fz(12,n0(1),n0(2),0:n0(3)))
    walls='periodic';walls(1:2)='no_slip_wall'
    do kind=1,2
      call fill_amr_patch_tree_ghosts_3d(tree,packed,1,walls,ghost,ok,viscous_wall=kind==2, &
        coarse_limiter='linear',species_count=size(species),ghost_width=3,species=species,limit_internal_energy=.true.)
      call require(ok,'no-slip MOL oracle halo')
      call compute_reactive_mol_halo_3d(species,ghost(:12,:,:,:),ghost(13,:,:,:),n0,spacing,dt, &
        expected,fx,fy,fz,b,ok)
      call require(ok,'no-slip MOL oracle flux')
      if(kind==1) slip_rate=expected
    enddo
    call require(maxval(abs(expected(3,:,:,:)-slip_rate(3,:,:,:)))>1e-6_dp, &
      'tangential gradient distinguishes no-slip from normal-only reflection')
    call compute_mol_amr_patch_tree_hydro_3d(tree,species,u,temperature,spacing,dt,walls,rates,flux,ok)
    call require(ok,'no-slip MOL tree rate')
    call require(all(rates%fields(1)%values==expected),'no-slip MOL tangential ghost reflection matches explicit halo')
  end subroutine check_no_slip_mol
  subroutine check_auxiliary_mol_stage()
    type(amr_tree_state_3d) :: aux,rates,arates,packed,saved_rates,saved_aux,level_rates,level_aux
    type(amr_tree_flux_3d), allocatable :: flux(:),aflux(:),level_flux(:),level_aflux(:)
    type(eb_integrated_flux_3d) :: reference_flux
    real(dp), allocatable :: ghost(:,:,:,:),reference(:,:,:,:),aux_reference(:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    real(dp), allocatable :: before(:,:,:,:),fractions(:)
    real(dp) :: b(12),pressure(6),divergence
    character(len=32) :: physical(6)
    integer :: peer,dims(3),x,y,z,depth
    aux=temperature;packed=state
    do peer=1,size(tree%nodes)
      ! Retained rhoe deliberately differs from E-KE. It must survive packing
      ! and participate in coarse interpolation, rather than being rebuilt.
      aux%fields(peer)%values=10000
      dims=tree%nodes(peer)%shape
      deallocate(packed%fields(peer)%values)
      allocate(packed%fields(peer)%values(14,dims(1),dims(2),dims(3)))
      packed%fields(peer)%values(:12,:,:,:)=state%fields(peer)%values
      packed%fields(peer)%values(13,:,:,:)=temperature%fields(peer)%values(1,:,:,:)
      packed%fields(peer)%values(14,:,:,:)=aux%fields(peer)%values(1,:,:,:)
    enddo
    call compute_mol_amr_patch_tree_hydro_3d(tree,species,state,temperature,spacing,dt,faces,rates,flux,ok, &
      internal_energy=aux,internal_energy_rate=arates,internal_energy_integrated_flux=aflux)
    call require(ok,'stored auxiliary MOL tree stage')
    do peer=1,size(tree%nodes)
      dims=tree%nodes(peer)%shape;h=spacing/2**tree%nodes(peer)%level
      call fill_amr_patch_tree_ghosts_3d(tree,packed,peer,faces,ghost,ok,coarse_limiter='linear', &
        species_count=size(species),ghost_width=3,species=species,limit_internal_energy=.true.,stored_internal_energy=.true.)
      call require(ok,'stored auxiliary MOL oracle halo')
      allocate(reference(12,dims(1),dims(2),dims(3)),aux_reference(dims(1),dims(2),dims(3)), &
        fx(12,0:dims(1),dims(2),dims(3)),fy(12,dims(1),0:dims(2),dims(3)),fz(12,dims(1),dims(2),0:dims(3)))
      call compute_reactive_mol_halo_3d(species,ghost(:12,:,:,:),ghost(13,:,:,:),dims,h,dt, &
        reference,fx,fy,fz,b,ok,internal_energy_rate=aux_reference,internal_energy_integrated_flux=reference_flux)
      call require(ok,'stored auxiliary MOL regular operator oracle')
      call require(all(rates%fields(peer)%values==reference),'stored MOL conserved rate matches kernel')
      call require(all(arates%fields(peer)%values(1,:,:,:)==aux_reference),'stored MOL auxiliary rate matches kernel')
      call require(all(aflux(peer)%x==reference_flux%x).and.all(aflux(peer)%y==reference_flux%y).and. &
        all(aflux(peer)%z==reference_flux%z),'stored MOL scalar face index conversion')
      do z=1,dims(3)
      do y=1,dims(2)
      do x=1,dims(1)
        divergence=-(aflux(peer)%x(1,x+1,y,z)-aflux(peer)%x(1,x,y,z))/(dt*h(1))- &
          (aflux(peer)%y(1,x,y+1,z)-aflux(peer)%y(1,x,y,z))/(dt*h(2))- &
          (aflux(peer)%z(1,x,y,z+1)-aflux(peer)%z(1,x,y,z))/(dt*h(3))
        call require(abs(divergence-aux_reference(x,y,z))<1e-10_dp*max(1._dp,abs(divergence)), &
          'auxiliary integrated flux divergence')
      enddo
      enddo
      enddo
      deallocate(reference,aux_reference,fx,fy,fz)
    enddo
    depth=tree%level_count()-1
    if(depth>0) then
      allocate(fractions(depth));fractions=.5_dp
      call compute_mol_amr_patch_tree_hydro_3d(tree,species,state,temperature,spacing,dt,faces, &
        level_rates,level_flux,ok,target_depth=depth,coarse_start=packed,coarse_end=packed,coarse_fractions=fractions, &
        internal_energy=aux,internal_energy_rate=level_aux,internal_energy_integrated_flux=level_aflux)
      call require(ok,'stored auxiliary coarse endpoint stage')
      do peer=1,size(tree%nodes)
        if(tree%nodes(peer)%level==depth) then
          call require(all(level_rates%fields(peer)%values==rates%fields(peer)%values).and. &
            all(level_aux%fields(peer)%values==arates%fields(peer)%values),'stored endpoint sampling exact')
          call require(all(level_aflux(peer)%x==aflux(peer)%x),'stored selected level scalar flux exact')
        else
          call require(all(level_aux%fields(peer)%values==0),'stored inactive level scalar rate zero')
          call require(.not.allocated(level_aflux(peer)%x),'stored inactive scalar flux absent')
        endif
      enddo
    endif
    saved_rates=rates;saved_aux=arates
    call compute_mol_amr_patch_tree_hydro_3d(tree,species,state,temperature,spacing,dt,faces,rates,flux,ok, &
      internal_energy=aux,internal_energy_rate=arates)
    call require(.not.ok,'missing auxiliary output rejection')
    do peer=1,size(tree%nodes)
      call require(all(rates%fields(peer)%values==saved_rates%fields(peer)%values).and. &
        all(arates%fields(peer)%values==saved_aux%fields(peer)%values),'auxiliary stage rollback')
    enddo
    if(size(tree%nodes)/=1) return
    physical='outflow';pressure=0;pressure(2)=102325
    call fill_amr_patch_tree_ghosts_3d(tree,packed,1,physical,ghost,ok,coarse_limiter='linear', &
      species_count=size(species),ghost_width=3,species=species,limit_internal_energy=.true.,stored_internal_energy=.true.)
    call require(ok,'stored pressure halo')
    before=ghost
    call apply_amr_tree_pressure_ghosts_3d(tree,species,1,physical,pressure,ghost,ok,ghost_width=3, &
      stored_internal_energy=.true.)
    call require(ok,'stored pressure ghost update')
    call require(all(ghost(14,:7,:,:)==before(14,:7,:,:)),'pressure leaves other retained scalars unchanged')
    do z=1,10
    do y=1,10
    do x=8,10
      divergence=ghost(5,x,y,z)-.5_dp*sum(ghost(2:4,x,y,z)**2)/ghost(1,x,y,z)
      call require(ghost(14,x,y,z)==divergence,'pressure ghost auxiliary follows prescribed state')
    enddo
    enddo
    enddo
  end subroutine check_auxiliary_mol_stage

  subroutine check_sdc_hydro_stage()
    type(amr_tree_state_3d) :: source,rates,candidate,expected,et,forced,level_rates,packed,future,source_future,retained
    type(amr_tree_flux_3d), allocatable :: flux(:),forced_flux(:),level_flux(:)
    real(dp) :: boundary(12),correction(12),maximum,change
    integer :: peer,depth,nv,dims(3)
    real(dp), allocatable :: fractions(:)
    source=state
    do peer=1,size(tree%nodes)
      source%fields(peer)%values=0
    end do
    call compute_sdc_amr_patch_tree_hydro_3d(tree,species,state,temperature,source,spacing,dt,faces, &
      rates,flux,ok)
    call require(ok,'SDC hierarchy CTU stage')
    candidate=state
    do peer=1,size(tree%nodes)
      candidate%fields(peer)%values=state%fields(peer)%values+dt*rates%fields(peer)%values
    end do
    call reflux_amr_patch_tree_3d(tree,candidate,flux,spacing,faces,correction,maximum,ok)
    call require(ok,'SDC stage integrated flux reflux')
    expected=state;et=temperature
    call advance_amr_patch_tree_hydro_3d(tree,species,expected,et,spacing,dt,faces, &
      'characteristic_ppm','mc','pelec',boundary,maximum,ok)
    call require(ok,'SDC hierarchy CTU reference')
    do peer=1,size(tree%nodes)
      call require(maxval(abs(candidate%fields(peer)%values-expected%fields(peer)%values)/ &
        max(1._dp,abs(expected%fields(peer)%values)))<1e-12_dp,'SDC hydro stage plus reflux matches existing tree step')
      source%fields(peer)%values(2,:,:,:)=20000
    end do
    call compute_sdc_amr_patch_tree_hydro_3d(tree,species,state,temperature,source,spacing,dt,faces, &
      forced,forced_flux,ok)
    call require(ok,'SDC hierarchy signed source interpolation')
    change=0
    do peer=1,size(tree%nodes)
      change=max(change,maxval(abs(forced%fields(peer)%values-rates%fields(peer)%values)))
    end do
    call require(change>1e-10_dp,'SDC predictor source changes CTU rate')
    do depth=0,tree%level_count()-1
      call compute_sdc_amr_patch_tree_hydro_3d(tree,species,state,temperature,source,spacing,dt,faces, &
        level_rates,level_flux,ok,target_depth=depth)
      call require(ok,'level-selective SDC hydro')
      do peer=1,size(tree%nodes)
        if(tree%nodes(peer)%level==depth) then
          call require(all(level_rates%fields(peer)%values==forced%fields(peer)%values),'selected level rates exact')
          call require(all(level_flux(peer)%x==forced_flux(peer)%x),'selected level flux exact')
        else
          call require(all(level_rates%fields(peer)%values==0),'inactive rate zero')
          call require(.not.allocated(level_flux(peer)%x),'inactive flux absent')
        end if
      end do
    end do
    depth=tree%level_count()-1
    if(depth==0) return
    allocate(packed%fields(size(tree%nodes)),fractions(depth))
    nv=size(state%fields(1)%values,1)
    do peer=1,size(tree%nodes)
      dims=tree%nodes(peer)%shape
      allocate(packed%fields(peer)%values(nv+1,dims(1),dims(2),dims(3)))
      packed%fields(peer)%values(:nv,:,:,:)=state%fields(peer)%values
      packed%fields(peer)%values(nv+1,:,:,:)=temperature%fields(peer)%values(1,:,:,:)
    end do
    future=packed;source_future=source
    do peer=1,size(tree%nodes)
      if(tree%nodes(peer)%level>=depth) cycle
      future%fields(peer)%values(5,:,:,:)=1.0001_dp*packed%fields(peer)%values(5,:,:,:)
      source_future%fields(peer)%values=2*source%fields(peer)%values
    end do
    fractions=0
    call compute_sdc_amr_patch_tree_hydro_3d(tree,species,state,temperature,source,spacing,dt,faces, &
      level_rates,level_flux,ok,target_depth=depth,coarse_start=packed,coarse_end=future,coarse_fractions=fractions, &
      source_start=source,source_end=source_future,source_fractions=fractions)
    call require(ok,'coarse old endpoint sampling')
    do peer=1,size(tree%nodes)
      if(tree%nodes(peer)%level/=depth) cycle
      call require(all(level_rates%fields(peer)%values==forced%fields(peer)%values),'old endpoint exact')
    end do
    fractions=.5_dp
    retained=level_rates
    call compute_sdc_amr_patch_tree_hydro_3d(tree,species,state,temperature,source,spacing,dt,faces, &
      level_rates,level_flux,ok,target_depth=depth,coarse_start=packed,coarse_end=future,coarse_fractions=fractions, &
      source_start=source,source_end=source_future,source_fractions=fractions)
    call require(ok,'coarse intermediate state and signed source sampling')
    change=0
    do peer=1,size(tree%nodes)
      change=max(change,maxval(abs(level_rates%fields(peer)%values-retained%fields(peer)%values)))
    end do
    call require(change>1e-10_dp,'coarse temporal donors change fine rate')
    retained=level_rates;fractions=1.1_dp
    call compute_sdc_amr_patch_tree_hydro_3d(tree,species,state,temperature,source,spacing,dt,faces, &
      level_rates,level_flux,ok,target_depth=depth,coarse_start=packed,coarse_end=future,coarse_fractions=fractions)
    call require(.not.ok,'invalid coarse time rejected')
    do peer=1,size(tree%nodes)
      call require(all(level_rates%fields(peer)%values==retained%fields(peer)%values),'invalid coarse time rollback')
    end do
  end subroutine check_sdc_hydro_stage
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
