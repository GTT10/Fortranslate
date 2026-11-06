program test_reactive_mol_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_finite
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use dual_energy_mod, only: dual_energy_policy
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved, &
    reactive_conserved_to_primitive
  use reactive_mol_3d_mod, only: compute_reactive_mol_3d,advance_reactive_mol_3d
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use mixture_thermo_mod, only: mixture_mass_properties
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:),q(:),u(:),rate(:,:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
  real(dp), allocatable :: boundary(:),expected(:,:)
  real(dp) :: inlet_t(6),mach(6),pressure(6)
  real(dp) :: t,c,spacing(3)=[.01_dp,.02_dp,.03_dp],dt=1e-7_dp,rho,residual,volume
  integer :: nv,n(3),i,j,k,mode,order,halo_calls,agreement_calls
  logical :: ok,bad_halo,inert_probe,fail_probe
  integer :: probe_calls,source_calls,water_index,birth_mode
  real(dp) :: birth_energy,birth_rate,birth_velocity(3)
  integer :: reject_window_mode=0
  character(len=32) :: faces(6)
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  nv=reactive_nvar(size(species));allocate(q(reactive_nprim(size(species))),u(nv),boundary(nv),expected(nv,3))
  q=0;q(1)=1;q(2:4)=[2._dp,-3._dp,4._dp];q(5)=1e5_dp;q(size(q))=1
  call reactive_primitive_to_conserved(species,q,u,t,c,ok);call require(ok,'uniform gas')
  call check_admitted_species_ratios()
  do mode=1,2
    n=[1,1,1]
    if(mode==2) n=[5,3,2]
    allocate(state(nv,n(1),n(2),n(3)),temperature(n(1),n(2),n(3)),rate(nv,n(1),n(2),n(3)))
    allocate(fx(nv,0:n(1),n(2),n(3)),fy(nv,n(1),0:n(2),n(3)),fz(nv,n(1),n(2),0:n(3)))
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      state(:,i,j,k)=u
    enddo
    enddo
    enddo
    temperature=t
    do i=1,3
      expected(:,i)=q(i+1)*u;expected(i+1,i)=expected(i+1,i)+q(5)
      expected(5,i)=q(i+1)*(u(5)+q(5))
    enddo
    do order=1,2
    do j=0,1
      call compute_reactive_mol_3d(species,state,temperature,spacing,dt,0._dp,rate,fx,fy,fz,boundary,ok, &
        order=order,lax_friedrichs=j==1)
      call require(ok,'regular MOL uniform pass')
      call require(maxval(abs(rate))<1e-8_dp.and.maxval(abs(boundary))<1e-20_dp,'uniform spatial balance')
      do i=1,nv
        call require(maxval(abs(fx(i,:,:,:)-dt*expected(i,1)))<1e-12_dp*max(1._dp,abs(dt*expected(i,1))), &
          'analytic x Euler flux')
        call require(maxval(abs(fy(i,:,:,:)-dt*expected(i,2)))<1e-12_dp*max(1._dp,abs(dt*expected(i,2))), &
          'analytic y Euler flux')
        call require(maxval(abs(fz(i,:,:,:)-dt*expected(i,3)))<1e-12_dp*max(1._dp,abs(dt*expected(i,3))), &
          'analytic z Euler flux')
      enddo
    enddo
    enddo
    ! A genuinely nonuniform periodic field exercises reconstruction, divergence
    ! and seam fluxes; composite volume changes must equal boundary inventories.
    if(mode==2) then
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        rho=1+.1_dp*cos(2*acos(-1._dp)*real(i-1,dp)/real(n(1),dp))
        q(1)=rho
        call reactive_primitive_to_conserved(species,q,state(:,i,j,k),temperature(i,j,k),c,ok)
        call require(ok,'nonuniform gas')
      enddo
      enddo
      enddo
      volume=product(spacing)
      do j=0,1
        call compute_reactive_mol_3d(species,state,temperature,spacing,dt,0._dp,rate,fx,fy,fz,boundary,ok, &
          lax_friedrichs=j==1)
        call require(ok.and.maxval(abs(rate))>1e-3_dp,'nontrivial periodic MOL rate')
        residual=0
        do i=1,nv
          residual=max(residual,abs(dt*volume*sum(rate(i,:,:,:))-boundary(i))/ &
            max(1e-25_dp,dt*volume*sum(abs(rate(i,:,:,:)))))
        enddo
        call require(residual<1e-12_dp,'periodic MOL finite-volume conservation')
        print *, 'Regular MOL LF flag, relative divergence residual:',j,residual
      enddo
    endif
    if(mode==1) then
      ! One-cell physical reflection and characteristic inlet need the same
      ! multi-layer donor mapping as ordinary distributed MOL tiles.
      q(2:4)=0
      call reactive_primitive_to_conserved(species,q,state(:,1,1,1),temperature(1,1,1),c,ok)
      call require(ok,'wall gas')
      faces='slip_wall'
      call compute_reactive_mol_3d(species,state,temperature,spacing,dt,0._dp,rate,fx,fy,fz,boundary,ok, &
        boundary_faces=faces)
      call require(ok.and.maxval(abs(rate))<1e-8_dp.and.maxval(abs(boundary))<1e-20_dp,'one-cell wall balance')
      q(2)=2
      call reactive_primitive_to_conserved(species,q,state(:,1,1,1),temperature(1,1,1),c,ok)
      call require(ok,'Mach inlet gas')
      faces='periodic';faces(1)='inflow';faces(2)='outflow'
      mach=-1;mach(1)=2/c;inlet_t=temperature(1,1,1);pressure=0;pressure(1)=q(5)
      call compute_reactive_mol_3d(species,state,temperature,spacing,dt,0._dp,rate,fx,fy,fz,boundary,ok, &
        boundary_faces=faces,inflow_mach=mach,inflow_temperature=inlet_t,inflow_pressure=pressure)
      call require(ok,'one-cell Mach inlet spatial evaluation')
      residual=maxval(abs(rate))*minval(spacing)/maxval(abs(state))
      call require(residual<1e-12_dp,'balanced one-cell Mach inlet')
      q(2:4)=[2._dp,-3._dp,4._dp]
    endif
    rate=17;fx=19;fy=23;fz=29;boundary=31
    halo_calls=0;agreement_calls=0;bad_halo=.true.
    call compute_reactive_mol_3d(species,state,temperature,spacing,dt,0._dp,rate,fx,fy,fz,boundary,ok, &
      fill_halos=halos,agree_stage=agreement)
    call require(.not.ok.and.halo_calls==1.and.agreement_calls>=2,'late halo EOS rejection')
    call require(all(rate==17).and.all(fx==19).and.all(fy==23).and.all(fz==29).and.all(boundary==31), &
      'failed spatial evaluation preserves all outputs')
    deallocate(state,temperature,rate,fx,fy,fz)
    q(1)=1
  enddo
  call check_auxiliary_interval()
  call check_mol_interval()
  call check_window_recovery()
  call check_rotating_interval()
  print *, 'PASS: regular MOL spatial flux, interval reaction/forcing, conservation and rollback'
contains
  subroutine check_admitted_species_ratios()
    real(dp) :: trial(size(u)),original(size(u)),primitive(size(q)),temp,sound,expected_mass
    real(dp) :: uniform(size(u),1,1,1),seed(1,1,1),rhs(size(u),1,1,1),ledger(size(u))
    real(dp) :: xflux(size(u),0:1,1,1),yflux(size(u),1,0:1,1),zflux(size(u),1,1,0:1)
    logical :: valid
    trial=u;trial(size(u))=trial(size(u))*(1._dp+2.e-12_dp);original=trial
    call reactive_conserved_to_primitive(species,trial,t,primitive,temp,sound,valid, &
      preserve_species_ratios=.true.)
    call require(valid.and.all(primitive(6:)==trial(6:)/trial(1)), &
      'admitted primitive retains rhoY/rho without a second normalization')
    call require(all(trial==original),'primitive recovery never projects conserved species')
    call reactive_conserved_to_primitive(species,trial,t,primitive,temp,sound,valid)
    call require(valid.and.abs(sum(primitive(6:))-1._dp)<epsilon(1._dp), &
      'default primitive composition remains normalized')
    uniform(:,1,1,1)=trial;seed=t
    call compute_reactive_mol_3d(species,uniform,seed,spacing,dt,0._dp,rhs, &
      xflux,yflux,zflux,ledger,valid)
    call require(valid,'MOL admits the original strict species closure bound')
    expected_mass=dt*(trial(2)/trial(1))*sum(trial(6:))
    call require(abs(xflux(1,0,1,1)-expected_mass)<8*epsilon(1._dp)*abs(expected_mass), &
      'uniform MOL face reconstructs the admitted partial density sum')
    trial(size(u))=trial(size(u))+1.e-7_dp*trial(1)
    call reactive_conserved_to_primitive(species,trial,t,primitive,temp,sound,valid, &
      preserve_species_ratios=.true.)
    call require(.not.valid,'preserving ratios still rejects species closure failure')
  end subroutine check_admitted_species_ratios
  subroutine check_auxiliary_interval()
    real(dp) :: gas(nv,2,1,1),temp(2,1,1),aux(2,1,1),kept(nv,2,1,1),kept_t(2,1,1), &
      kept_aux(2,1,1),ledger,change(nv),chemical(nv),b(nv),prim(nv),cs,tt,expected_aux
    type(dual_energy_policy) :: policy
    type(eb_integrated_flux_3d) :: aux_flux,registered
    integer :: cell,lax,interpolation
    prim=q;prim(1)=1;prim(2:4)=[2._dp,-3._dp,4._dp]
    prim(5)=3e5_dp
    do cell=1,2
      call reactive_primitive_to_conserved(species,prim,gas(:,cell,1,1),tt,cs,ok)
      call require(ok,'auxiliary interval initial EOS')
      temp(cell,1,1)=tt
    enddo
    aux=gas(5,:,:,:)-.5_dp*sum(gas(2:4,:,:,:)**2,dim=1)/gas(1,:,:,:)
    kept=gas;kept_t=temp;kept_aux=aux
    do lax=0,1
    do interpolation=1,2
      gas=kept;temp=kept_t;aux=kept_aux;ledger=7
      call advance_reactive_mol_3d(species,gas,temp,[1._dp,1._dp,1._dp],dt,0._dp,1,.false., &
        1e-7_dp,1e-10_dp,ok,lax_friedrichs=lax==1,order=interpolation, &
        internal_energy_density=aux,energy_correction=ledger,internal_energy_integrated_flux=aux_flux)
      call require(ok.and.ledger==7,'default auxiliary MOL interval')
      call require(maxval(abs(gas-kept))<1e-9_dp,'uniform auxiliary interval preserves conserved state')
      call require(maxval(abs(aux-(gas(5,:,:,:)-.5_dp*sum(gas(2:4,:,:,:)**2,dim=1)/gas(1,:,:,:))))<1e-9_dp, &
        'default reset retains auxiliary from accepted total')
      expected_aux=dt*prim(2)*kept_aux(1,1,1)
      call require(maxval(abs(aux_flux%x-expected_aux))<1e-12_dp,'auxiliary trapezoidal integrated face flux')
    enddo
    enddo
    policy%allow_negative_energy=.false.;policy%eta2=1
    gas=kept;temp=kept_t;aux=kept_aux;ledger=7;fail_probe=.false.
    call advance_reactive_mol_3d(species,gas,temp,[1._dp,1._dp,1._dp],dt,0._dp,1,.false., &
      1e-7_dp,1e-10_dp,ok,source_count=1,evaluate_sources=energy_sink,source_change=change,chemical_change=chemical, &
      boundary_change=b,internal_energy_density=aux,energy_correction=ledger,energy_policy=policy, &
      internal_energy_integrated_flux=aux_flux,internal_energy_registered_flux=registered)
    call require(ok,'retained auxiliary policy interval')
    ! Two unit-volume cells each receive -10 J and a +10 J reset correction.
    call require(abs(ledger-27)<1e-8_dp.and.abs(change(5)+20)<1e-8_dp,'weighted predictor/final reset ledger')
    call require(maxval(abs(gas-kept))<1e-8_dp.and.maxval(abs(aux-kept_aux))<1e-8_dp,'retained auxiliary energy repair')
    call require(abs(sum(gas(5,:,:,:)-kept(5,:,:,:))-change(5)-(ledger-7)-b(5))<1e-8_dp, &
      'total-energy conservation includes policy correction')
    call require(abs(chemical(5))<1e-8_dp,'reset correction is not chemical energy')
    call require(all(aux_flux%x==registered%x),'one-iteration auxiliary register equals physical flux')
    gas=kept;temp=kept_t;aux=kept_aux;ledger=31;fail_probe=.true.
    call advance_reactive_mol_3d(species,gas,temp,[1._dp,1._dp,1._dp],dt,0._dp,1,.false., &
      1e-7_dp,1e-10_dp,ok,source_count=1,evaluate_sources=energy_sink, &
      internal_energy_density=aux,energy_correction=ledger,energy_policy=policy)
    call require(.not.ok.and.all(gas==kept).and.all(temp==kept_t).and.all(aux==kept_aux).and.ledger==31, &
      'late source failure preserves auxiliary field and correction ledger')
    fail_probe=.false.
  end subroutine check_auxiliary_interval
  subroutine energy_sink(gas,temp,stamp,families,good)
    real(dp), intent(in) :: gas(:,:,:,:),temp(:,:,:),stamp
    real(dp), intent(out) :: families(:,:,:,:,:)
    logical, intent(out) :: good
    families=0;families(5,:,:,:,1)=-1e8_dp
    good=all(ieee_is_finite(gas)).and.all(ieee_is_finite(temp))
    if(fail_probe.and.stamp>0) good=.false.
  end subroutine energy_sink
  subroutine check_mol_interval()
    type(elementary_reaction), allocatable :: reactions(:)
    real(dp) :: mass_y(size(species)),prim(reactive_nprim(size(species))),v(nv),mw,r,cp,cv,gamma,h,e,entropy
    real(dp) :: gas(nv,2,2,2),temp(2,2,2),ir(nv-4,2,2,2),heat(2,2,2),reference_y(size(species)),reference_t
    real(dp) :: boundary(nv),source(nv),chemical(nv),initial_t,weight,error,volume
    real(dp) :: held(2,2,2),held_ledger
    type(eb_integrated_flux_3d) :: physical,registered
    integer :: s,it,steps,x,y,z
    call load_h2o2_elementary_mechanism(reactions,ok);call require(ok,'interval mechanism')
    mass_y=0
    do s=1,size(species)
      select case(species(s)%name)
      case('H2');mass_y(s)=.1_dp
      case('O2');mass_y(s)=.2_dp
      case('N2');mass_y(s)=.699_dp
      case('H');mass_y(s)=.001_dp
      end select
    enddo
    call mixture_mass_properties(species,mass_y,1200._dp,mw,r,cp,cv,gamma,h,e,entropy,ok)
    call require(ok,'interval properties')
    prim=0;prim(1)=1;prim(5)=r*1200;prim(6:)=mass_y
    call reactive_primitive_to_conserved(species,prim,v,initial_t,c,ok);call require(ok,'interval gas')
    reference_y=v(6:);reference_t=initial_t
    call advance_forced_reactor_interval(species,reactions,v(5),0*mass_y,0._dp,1e-10_dp,1e-8_dp,1e-12_dp, &
      reference_y,reference_t,steps,ok)
    call require(ok.and.maxval(abs(reference_y-v(6:)))>1e-10_dp,'independent reacting oracle')
    volume=8*product(spacing);fail_probe=.false.
    do it=1,3
      call initialize_gas(gas,temp,v,initial_t)
      ir=0;heat=37;probe_calls=0;inert_probe=.false.
      call advance_reactive_mol_3d(species,gas,temp,spacing,1e-10_dp,0._dp,it,.true.,1e-8_dp,1e-12_dp,ok, &
        reactions=reactions,reaction_source=ir,heat_release=heat,evaluate_reaction=probe_reactor, &
        boundary_change=boundary,source_change=source,chemical_change=chemical, &
        integrated_flux=physical,registered_flux=registered)
      call require(ok.and.probe_calls==8*it,'MOL reaction iteration count')
      error=0
      do z=1,2
      do y=1,2
      do x=1,2
        error=max(error,maxval(abs(gas(6:,x,y,z)-reference_y)),abs(temp(x,y,z)-reference_t)/reference_t)
      enddo
      enddo
      enddo
      call require(error<1e-11_dp,'MOL homogeneous independent reactor agreement')
      call require(maxval(abs(volume*(gas(:,1,1,1)-v)-boundary-source-chemical))<1e-12_dp, &
        'MOL homogeneous chemical inventory')
      print *, 'Regular MOL iterations and reacting oracle error:',it,error
      call initialize_gas(gas,temp,v,initial_t)
      ir=0;heat=37;probe_calls=0;inert_probe=.false.;held=v(5);held_ledger=7
      call advance_reactive_mol_3d(species,gas,temp,spacing,1e-10_dp,0._dp,it,.true.,1e-8_dp,1e-12_dp,ok, &
        reactions=reactions,reaction_source=ir,heat_release=heat,evaluate_reaction=probe_reactor, &
        boundary_change=boundary,source_change=source,chemical_change=chemical, &
        internal_energy_density=held,energy_correction=held_ledger)
      call require(ok.and.probe_calls==8*it.and.held_ledger==7,'held MOL reaction and numerical-energy ledger')
      do z=1,2
      do y=1,2
      do x=1,2
        call require(maxval(abs(gas(6:,x,y,z)-reference_y))<1e-11_dp, &
          'held MOL chemical trajectory matches independent native reactor')
        call require(held(x,y,z)==v(5),'held MOL retains forced internal-energy endpoint')
      enddo
      enddo
      enddo
      call initialize_gas(gas,temp,v,initial_t)
      ir=0;heat=37;probe_calls=0;source_calls=0;inert_probe=.true.
      call advance_reactive_mol_3d(species,gas,temp,spacing,1e-10_dp,0._dp,it,.true.,1e-8_dp,1e-12_dp,ok, &
        reactions=reactions,reaction_source=ir,heat_release=heat,evaluate_reaction=probe_reactor, &
        source_count=1,evaluate_sources=interval_sources,boundary_change=boundary,source_change=source)
      call require(ok.and.source_calls==2,'MOL old/new nonspatial sources evaluated once')
      weight=merge(1._dp,.5_dp,it==1)
      call require(maxval(abs(gas(2,:,:,:)-weight*1e-10_dp*1e7_dp))<1e-14_dp, &
        'later MOL passes retain only half old nonspatial source')
      call require(abs(source(2)-volume*weight*1e-10_dp*1e7_dp)<1e-20_dp,'MOL source inventory matches actual impulse')
    enddo
    ! A full first reacted pass accepts privately, then the next reactor fails.
    call initialize_gas(gas,temp,v,initial_t)
    ir=0;heat=37;probe_calls=0;inert_probe=.false.;fail_probe=.true.
    boundary=41;source=42;chemical=43
    physical%x=17;physical%y=17;physical%z=17;registered%x=19;registered%y=19;registered%z=19
    call advance_reactive_mol_3d(species,gas,temp,spacing,1e-10_dp,0._dp,2,.true.,1e-8_dp,1e-12_dp,ok, &
      reactions=reactions,reaction_source=ir,heat_release=heat,evaluate_reaction=probe_reactor, &
      boundary_change=boundary,source_change=source,chemical_change=chemical, &
      integrated_flux=physical,registered_flux=registered)
    call require(.not.ok.and.probe_calls==9,'MOL later reactor rejection')
    call require(all(gas(:,1,1,1)==v).and.all(temp==initial_t).and.all(ir==0).and.all(heat==37), &
      'MOL gas/IR/heat rollback')
    call require(all(boundary==41).and.all(source==42).and.all(chemical==43),'MOL inventory rollback')
    call require(all(physical%x==17).and.all(physical%y==17).and.all(physical%z==17).and. &
      all(registered%x==19).and.all(registered%y==19).and.all(registered%z==19),'MOL flux ledger rollback')
    fail_probe=.false.
  end subroutine check_mol_interval
  subroutine check_window_recovery()
    type(elementary_reaction), allocatable :: rx(:)
    type(eb_integrated_flux_3d) :: flow_flux,flow_register,held_flux,held_register,kept_flux,kept_held_flux
    real(dp) :: gas(nv,3,1,1),temp(3,1,1),held(3,1,1),ir(nv-4,3,1,1),heat(3,1,1)
    real(dp) :: base(nv,3,1,1),base_t(3,1,1),base_h(3,1,1),held_prefix,blend
    real(dp) :: b(nv),src(nv),chem(nv),prim(nv),v(nv),tt,cs,composition(size(species))
    real(dp) :: mw,rg,cp,cv,gamma,enthalpy,e,entropy,delta,scale
    integer :: cell,it,variant,component
    call load_h2o2_elementary_mechanism(rx,ok);call require(ok,'recovery real mechanism')
    water_index=0
    do component=1,size(species)
      if(species(component)%name=='H2O') water_index=component
    enddo
    call require(water_index>0,'water tracer identity')
    prim=0;prim(1)=.3_dp;prim(2)=20;prim(5)=101325;prim(nv)=1
    call reactive_primitive_to_conserved(species,prim,v,tt,cs,ok);call require(ok,'dry recovery gas')
    do cell=1,3
      base(:,cell,1,1)=v;base_t(cell,1,1)=tt
    enddo
    base_h=base(5,:,:,:)-.5_dp*sum(base(2:4,:,:,:)**2,dim=1)/base(1,:,:,:)
    composition=0;composition(water_index)=1
    call mixture_mass_properties(species,composition,tt,mw,rg,cp,cv,gamma,enthalpy,e,entropy,ok)
    call require(ok,'birth internal energy')
    birth_velocity=prim(2:4);birth_energy=e+.5_dp*sum(birth_velocity**2);birth_rate=1000
    ! The first new source creates a wet cell. Later pinned MOL passes omit
    ! that source, leaving a negative water flux from a dry old donor. A
    ! conservative whole-window fallback must freeze the two source calls.
    reject_window_mode=0
    do variant=0,1
    do it=1,2
      birth_mode=variant;gas=base;temp=base_t;held=base_h;held_prefix=11;ir=0;heat=37;source_calls=0;blend=31
      call advance_reactive_mol_3d(species,gas,temp,spacing,dt,0._dp,it,.true.,1e-6_dp,1e-12_dp,ok, &
        reactions=rx,reaction_source=ir,heat_release=heat,source_count=1,evaluate_sources=window_birth, &
        boundary_change=b,source_change=src,chemical_change=chem,integrated_flux=flow_flux,registered_flux=flow_register, &
        internal_energy_density=held,energy_correction=held_prefix,internal_energy_integrated_flux=held_flux, &
        internal_energy_registered_flux=held_register,minimum_flux_blend=blend)
      call require(ok.and.source_calls==2,'natural MOL recovery without source replay')
      if(variant==1.and.it==2) then
        call require(blend==0,'dry tracer selects conservative old spatial family')
        call require(all(gas(6:,:,:,:)>=0),'recovered gas species remain physical')
        call require(all(flow_flux%x==flow_register%x).and.all(held_flux%x==held_register%x), &
          'recovered physical and registered gas/UEINT families agree')
        call require(maxval(abs(held_flux%x-dt*prim(2)*base_h(1,1,1)))<1e-12_dp, &
          'UEINT face carries the same zero-blend donor family')
      else
        call require(blend==1,'admissible MOL interval retains original spatial family')
      endif
      do component=1,nv
        delta=product(spacing)*sum(gas(component,:,:,:)-base(component,:,:,:))
        scale=max(1._dp,product(spacing)*sum(abs(base(component,:,:,:))))
        call require(abs(delta-b(component)-src(component)-chem(component))<=128*epsilon(1._dp)*scale, &
          'recovered full gas inventory equals face/source/chemical ledger')
      enddo
      call require(held_prefix==11,'recovery does not fabricate numerical reset history')
      if(variant==1.and.it==2) then
        call require(maxval(abs(gas-base))<=128*epsilon(1._dp)*max(1._dp,maxval(abs(base))), &
          'dry recovery endpoint equals old uniform physical family')
        call require(all(held==base_h),'dry recovery retains actual old UEINT')
      endif
    enddo
    enddo
    kept_flux=flow_flux;kept_held_flux=held_flux
    do variant=1,2
    gas=base;temp=base_t;held=base_h;held_prefix=11;ir=0;heat=37;source_calls=0;blend=31
    b=41;src=42;chem=43;birth_mode=1;reject_window_mode=variant
    call advance_reactive_mol_3d(species,gas,temp,spacing,dt,0._dp,2,.true.,1e-6_dp,1e-12_dp,ok, &
      reactions=rx,reaction_source=ir,heat_release=heat,source_count=1,evaluate_sources=window_birth, &
      boundary_change=b,source_change=src,chemical_change=chem,integrated_flux=flow_flux, &
      internal_energy_density=held,energy_correction=held_prefix,internal_energy_integrated_flux=held_flux, &
      minimum_flux_blend=blend,agree_flux_blend=reject_recovery)
    call require(.not.ok.and.source_calls==2,'late window coefficient rejection without source replay')
    call require(all(gas==base).and.all(temp==base_t).and.all(held==base_h).and.held_prefix==11.and. &
      all(ir==0).and.all(heat==37).and.blend==31,'late recovery gas/UEINT/IR/heat/prefix rollback')
    call require(all(b==41).and.all(src==42).and.all(chem==43),'late recovery inventory rollback')
    call require(all(flow_flux%x==kept_flux%x).and.all(flow_flux%y==kept_flux%y).and. &
      all(flow_flux%z==kept_flux%z).and.all(held_flux%x==kept_held_flux%x).and. &
      all(held_flux%y==kept_held_flux%y).and.all(held_flux%z==kept_held_flux%z),'late recovery paired face rollback')
    enddo
    reject_window_mode=0
    ! An old birth and a new sink cancel as sources, but the new wet face
    ! operator would remove more water than the old donor holds. Exercise a
    ! strict fractional face coefficient, not only the zero-budget case.
    prim(6:)=0;prim(nv)=1-1e-10_dp/prim(1);prim(5+water_index)=1e-10_dp/prim(1)
    call reactive_primitive_to_conserved(species,prim,base(:,2,1,1),base_t(2,1,1),cs,ok)
    call require(ok,'fractional recovery old wet donor')
    gas=base;temp=base_t;birth_mode=2;source_calls=0;blend=31
    call advance_reactive_mol_3d(species,gas,temp,spacing,dt,0._dp,1,.false.,1e-6_dp,1e-12_dp,ok, &
      source_count=1,evaluate_sources=window_birth,boundary_change=b,source_change=src,chemical_change=chem, &
      integrated_flux=flow_flux,registered_flux=flow_register,minimum_flux_blend=blend)
    call require(ok.and.source_calls==2.and.blend>0.and.blend<1,'fractional whole-window recovery without source replay')
    call require(all(gas(6:,:,:,:)>=0),'fractional recovery retains physical species')
    delta=sum(gas(5+water_index,:,:,:)-base(5+water_index,:,:,:))
    scale=sum(base(5+water_index,:,:,:))
    call require(abs(delta)<=128*epsilon(1._dp)*scale,'fractional repair conserves represented trace water')
    call require(all(src==0),'opposite frozen kicks cancel without source clipping')
    print *, 'Fractional MOL face blend:',blend
    print *, 'PASS: natural MOL dry-tracer whole-window repair, actual UEINT faces, frozen sources and late rollback'
  end subroutine check_window_recovery
  subroutine window_birth(gas,temp,stamp,families,good)
    real(dp), intent(in) :: gas(:,:,:,:),temp(:,:,:),stamp
    real(dp), intent(out) :: families(:,:,:,:,:)
    logical, intent(out) :: good
    source_calls=source_calls+1;families=0
    if((birth_mode==1.and.stamp>0).or.birth_mode==2) then
      families(1,2,1,1,1)=birth_rate
      families(2:4,2,1,1,1)=birth_rate*birth_velocity
      families(5,2,1,1,1)=birth_rate*birth_energy
      families(5+water_index,2,1,1,1)=birth_rate
      if(birth_mode==2.and.stamp>0) families=-families
    endif
    good=all(ieee_is_finite(gas)).and.all(ieee_is_finite(temp))
  end subroutine window_birth
  subroutine reject_recovery(local_blend,global_blend,good)
    real(dp), intent(in) :: local_blend
    real(dp), intent(out) :: global_blend
    logical, intent(out) :: good
    good=reject_window_mode/=1
    if(.not.good) return
    global_blend=local_blend
    if(reject_window_mode==2) global_blend=ieee_value(0._dp,ieee_quiet_nan)
  end subroutine reject_recovery
  subroutine check_rotating_interval()
    type(rotating_frame_3d) :: rotation,zero_rotation
    type(elementary_reaction), allocatable :: rx(:)
    type(gas_transport_species), allocatable :: db(:)
    real(dp) :: gas(nv,1,1,1),temp(1,1,1),reference(nv,1,1,1),reference_t(1,1,1)
    real(dp) :: seed(nv),initial_t,sound,phi,old_momentum(3),old_source(3),new_source(3),predicted(3)
    real(dp) :: b(nv),ext(nv),chem(nv),before(nv),error,interval,mw,r,cp,cv,gamma,h,e,entropy
    real(dp) :: reference_y(size(species)),oracle_t
    integer :: reactor_steps
    real(dp) :: small_h(3),center(3),offset(3),probe_ir(nv-4,1,1,1),probe_heat(1,1,1)
    type(eb_integrated_flux_3d) :: physical
    integer :: pass
    bad_halo=.false.
    q=0;q(1)=1;q(5)=1e5_dp;q(size(q))=1
    call reactive_primitive_to_conserved(species,q,seed,initial_t,sound,ok)
    call require(ok,'rotating regular initial gas')
    small_h=[.001_dp,.002_dp,.003_dp];offset=[.02_dp,-.03_dp,.01_dp];center=offset+.5_dp*small_h
    interval=1e-8_dp;rotation%angular_speed=300;rotation%axis=3
    phi=.5_dp*rotation%angular_speed**2*sum(center(:2)**2)
    gas(:,1,1,1)=seed;gas(5,1,1,1)=seed(5)-seed(1)*phi;temp=initial_t
    before=gas(:,1,1,1);old_momentum=before(2:4)
    old_source=before(1)*rotation%angular_speed**2*[center(1),center(2),0._dp]
    predicted=old_momentum+interval*old_source
    new_source=old_source+2*rotation%angular_speed*[predicted(2),-predicted(1),0._dp]
    call advance_reactive_mol_3d(species,gas,temp,small_h,interval,0._dp,2,.false.,1e-8_dp,1e-12_dp,ok, &
      frame=rotation,domain_lower=offset,boundary_change=b,source_change=ext,chemical_change=chem, &
      integrated_flux=physical)
    call require(ok,'regular rotating MOL interval')
    error=maxval(abs(gas(2:4,1,1,1)-(old_momentum+.5_dp*interval*(old_source+new_source))))
    call require(error<1e-13_dp,'regular rotating trapezoidal Coriolis/centrifugal impulse')
    print *, 'Regular rotating MOL momentum oracle error:',error
    call require(ext(5)==0.and.all(chem==0),'rotation adds no energy source')
    call require(maxval(abs(product(small_h)*(gas(:,1,1,1)-before)-b-ext-chem))<1e-16_dp, &
      'regular rotating published energy/flux inventory')
    call load_h2o2_elementary_transport(db,ok);call require(ok,'regular rotating transport database')
    reference=gas;reference_t=temp
    do pass=1,2
      gas(:,1,1,1)=before;temp=initial_t
      if(pass==1) then
        call advance_reactive_mol_3d(species,gas,temp,small_h,interval,0._dp,2,.false.,1e-8_dp,1e-12_dp,ok, &
          frame=rotation,domain_lower=offset,transport=db,transport_flags=[.true.,.true.,.true.,.false.])
      else
        call advance_reactive_mol_3d(species,gas,temp,small_h,interval,0._dp,2,.false.,1e-8_dp,1e-12_dp,ok, &
          frame=rotation,domain_lower=offset,transport=db,transport_flags=[.true.,.true.,.true.,.false.], &
          fill_halos=rotating_halos)
      endif
      call require(ok,'regular rotating molecular transport with local/callback halos')
      call require(maxval(abs(gas-reference))<1e-9_dp.and.maxval(abs(temp-reference_t))<1e-9_dp, &
        'homogeneous ordinary-energy transport is zero in rotating storage')
    enddo
    ! A real homogeneous reactor at nonzero potential: cancel rotation's
    ! momentum source, leaving chemistry as an independent cell oracle.
    interval=1e-10_dp
    q=0;q(1)=1;q(6)=.1_dp;q(7)=.001_dp;q(9)=.2_dp;q(size(q))=.699_dp
    call mixture_mass_properties(species,q(6:),1200._dp,mw,r,cp,cv,gamma,h,e,entropy,ok)
    call require(ok,'rotating reacting properties');q(5)=r*1200
    call reactive_primitive_to_conserved(species,q,seed,initial_t,sound,ok)
    call require(ok,'rotating regular reacting initial gas')
    call load_h2o2_elementary_mechanism(rx,ok);call require(ok,'rotating regular chemistry')
    do pass=1,3
      gas(:,1,1,1)=seed;gas(5,1,1,1)=seed(5)-seed(1)*phi;temp=initial_t
      reference_y=seed(6:);oracle_t=initial_t
      call advance_forced_reactor_interval(species,rx,seed(5),0*reference_y,0._dp,interval,1e-8_dp,1e-12_dp, &
        reference_y,oracle_t,reactor_steps,ok)
      call require(ok.and.maxval(abs(reference_y-seed(6:)))>1e-10_dp,'independent rotating chemistry oracle')
      reference(:,1,1,1)=seed;reference(6:,1,1,1)=reference_y;reference_t=oracle_t
      probe_ir=0;probe_heat=37
      call advance_reactive_mol_3d(species,gas,temp,small_h,interval,0._dp,pass,.true.,1e-8_dp,1e-12_dp,ok, &
        reactions=rx,frame=rotation,domain_lower=offset,reaction_source=probe_ir,heat_release=probe_heat, &
        source_count=1,evaluate_sources=cancel_rotation,boundary_change=b,source_change=ext,chemical_change=chem)
      call require(ok,'regular rotating real chemistry')
      error=maxval(abs(gas(6:,1,1,1)-reference(6:,1,1,1)))
      call require(error<1e-12_dp.and.abs(temp(1,1,1)-reference_t(1,1,1))<1e-8_dp, &
        'nonzero-potential regular MOL independent homogeneous reactor')
      print *, 'Regular rotating MOL iterations, species/temperature oracle error:',pass,error, &
        abs(temp(1,1,1)-reference_t(1,1,1))
    enddo
    ! Zero omega retains the ordinary interval; invalid frame rejects privately.
    reference(:,1,1,1)=seed;reference_t=initial_t
    gas=reference;temp=reference_t
    call advance_reactive_mol_3d(species,reference,reference_t,small_h,interval,0._dp,1,.false., &
      1e-8_dp,1e-12_dp,ok)
    call require(ok,'ordinary zero-frame control')
    call advance_reactive_mol_3d(species,gas,temp,small_h,interval,0._dp,1,.false.,1e-8_dp,1e-12_dp,ok, &
      frame=zero_rotation,domain_lower=offset)
    call require(ok.and.all(gas==reference).and.all(temp==reference_t),'zero rotation interval exact identity')
    before=gas(:,1,1,1);b=41;ext=42;chem=43;probe_ir=31;probe_heat=37;bad_halo=.true.
    call advance_reactive_mol_3d(species,gas,temp,small_h,interval,0._dp,1,.false.,1e-8_dp,1e-12_dp,ok, &
      frame=rotation,domain_lower=offset,boundary_change=b,source_change=ext,chemical_change=chem, &
      fill_halos=rotating_halos,reaction_source=probe_ir,heat_release=probe_heat)
    call require(.not.ok.and.all(gas(:,1,1,1)==before).and.all(temp==reference_t), &
      'late rotating halo failure preserves state/temperature')
    call require(all(b==41).and.all(ext==42).and.all(chem==43).and.all(probe_ir==31).and.all(probe_heat==37), &
      'late rotating halo failure preserves ledgers/IR/heat')
    bad_halo=.false.
    rotation%axis=4;before=gas(:,1,1,1);b=41;ext=42;chem=43
    call advance_reactive_mol_3d(species,gas,temp,small_h,interval,0._dp,1,.false.,1e-8_dp,1e-12_dp,ok, &
      frame=rotation,domain_lower=offset,boundary_change=b,source_change=ext,chemical_change=chem)
    call require(.not.ok.and.all(gas(:,1,1,1)==before).and.all(temp==reference_t), &
      'invalid regular rotating frame rollback')
    call require(all(b==41).and.all(ext==42).and.all(chem==43),'invalid rotation ledger rollback')
    print *, 'PASS: regular rotating MOL impulse, nonzero-potential chemistry and rollback'
  end subroutine check_rotating_interval
  subroutine rotating_halos(gas,temp,time,width,viscous,ghost,gt,good, &
      internal_energy_density,internal_energy_ghost)
    real(dp), optional, intent(in) :: internal_energy_density(:,:,:)
    real(dp), optional, intent(out) :: internal_energy_ghost(:,:,:)
    real(dp), intent(in) :: gas(:,:,:,:),temp(:,:,:),time
    integer, intent(in) :: width
    logical, intent(in) :: viscous
    real(dp), intent(out) :: ghost(:,:,:,:),gt(:,:,:)
    logical, intent(out) :: good
    real(dp) :: ordinary(size(gas,1),1,1,1),center(3),position(3),h(3),lower(3),phi
    integer :: x,y,z
    lower=[.02_dp,-.03_dp,.01_dp];h=[.001_dp,.002_dp,.003_dp];center=lower+.5_dp*h
    ordinary=gas;ordinary(5,1,1,1)=gas(5,1,1,1)+gas(1,1,1,1)*.5_dp*300._dp**2*sum(center(:2)**2)
    call fill_reactive_boundary_ghosts_3d(ordinary,temp,[character(len=8):: &
      'periodic','periodic','periodic','periodic','periodic','periodic'],ghost,gt,good,ghost_width=width, &
      internal_energy_density=internal_energy_density,internal_energy_ghost=internal_energy_ghost)
    if(.not.good) return
    do z=1,size(gt,3)
    do y=1,size(gt,2)
    do x=1,size(gt,1)
      position=lower+(real([x,y,z]-width,dp)-.5_dp)*h;phi=.5_dp*300._dp**2*sum(position(:2)**2)
      ghost(5,x,y,z)=ghost(5,x,y,z)-ghost(1,x,y,z)*phi
    enddo
    enddo
    enddo
    if(bad_halo) ghost(5,size(ghost,2),size(ghost,3),size(ghost,4))=ieee_value(0._dp,ieee_quiet_nan)
  end subroutine rotating_halos
  subroutine cancel_rotation(gas,temp,time,rates,good)
    real(dp), intent(in) :: gas(:,:,:,:),temp(:,:,:),time
    real(dp), intent(out) :: rates(:,:,:,:,:)
    logical, intent(out) :: good
    real(dp) :: center(3),momentum(3),omega
    center=[.0205_dp,-.029_dp,.0115_dp];omega=300
    momentum=gas(2:4,1,1,1)
    rates=0
    rates(2:4,1,1,1,1)=-gas(1,1,1,1)*omega**2*[center(1),center(2),0._dp]- &
      2*omega*[momentum(2),-momentum(1),0._dp]
    good=.true.
  end subroutine cancel_rotation
  subroutine initialize_gas(gas,temp,v,initial_t)
    real(dp), intent(out) :: gas(:,:,:,:),temp(:,:,:)
    real(dp), intent(in) :: v(:),initial_t
    integer :: x,y,z
    do z=1,size(temp,3)
    do y=1,size(temp,2)
    do x=1,size(temp,1)
      gas(:,x,y,z)=v
    enddo
    enddo
    enddo
    temp=initial_t
  end subroutine initialize_gas
  subroutine interval_sources(gas,temp,time,rates,good)
    real(dp), intent(in) :: gas(:,:,:,:),temp(:,:,:),time
    real(dp), intent(out) :: rates(:,:,:,:,:)
    logical, intent(out) :: good
    source_calls=source_calls+1;rates=0;rates(2,:,:,:,1)=1e7_dp;good=.true.
  end subroutine interval_sources
  subroutine probe_reactor(sp,rx,rhoe,external_y,external_e,duration,rtol,atol,densities,temp,count,good,maximum_steps)
    type(nasa7_species), intent(in) :: sp(:)
    type(elementary_reaction), intent(in) :: rx(:)
    real(dp), intent(in) :: rhoe,external_y(:),external_e,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),temp
    integer, intent(out) :: count
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    type(elementary_reaction) :: empty(0)
    probe_calls=probe_calls+1;good=.false.;count=0
    if(fail_probe.and.probe_calls>8) return
    if(inert_probe) then
      call advance_forced_reactor_interval(sp,empty,rhoe,external_y,external_e,duration,rtol,atol, &
        densities,temp,count,good,maximum_steps=maximum_steps)
    else
      call advance_forced_reactor_interval(sp,rx,rhoe,external_y,external_e,duration,rtol,atol, &
        densities,temp,count,good,maximum_steps=maximum_steps)
    endif
  end subroutine probe_reactor
  subroutine require(good,label)
    logical, intent(in) :: good
    character(len=*), intent(in) :: label
    if(.not.good) then
      print *, 'FAIL ',label
      error stop 1
    endif
  end subroutine require
  subroutine halos(gas,temp,time,width,viscous,ghost,gt,good, &
      internal_energy_density,internal_energy_ghost)
    real(dp), optional, intent(in) :: internal_energy_density(:,:,:)
    real(dp), optional, intent(out) :: internal_energy_ghost(:,:,:)
    real(dp), intent(in) :: gas(:,:,:,:),temp(:,:,:),time
    integer, intent(in) :: width
    logical, intent(in) :: viscous
    real(dp), intent(out) :: ghost(:,:,:,:),gt(:,:,:)
    logical, intent(out) :: good
    halo_calls=halo_calls+1;faces='periodic'
    call fill_reactive_boundary_ghosts_3d(gas,temp,faces,ghost,gt,good,ghost_width=width, &
      internal_energy_density=internal_energy_density,internal_energy_ghost=internal_energy_ghost)
    if(bad_halo) ghost(5,size(ghost,2)-1,size(ghost,3)-1,size(ghost,4)-1)=ieee_value(0._dp,ieee_quiet_nan)
  end subroutine halos
  subroutine agreement(local_ok,global_ok)
    logical, intent(in) :: local_ok
    logical, intent(out) :: global_ok
    agreement_calls=agreement_calls+1;global_ok=local_ok
  end subroutine agreement
end program test_reactive_mol_3d
