program test_regular_population_sdc
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  use mpi_f08
  use precision_mod, only: dp
  use dual_energy_mod, only: dual_energy_policy
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,wall_film_totals
  use spray_mpi_amr_film_mod, only: spray_level_film_population,initialize_mpi_amr_film_levels
  use spray_mpi_amr_film_step_mod, only: spray_sdc_film_stage,advance_mpi_amr_film_sdc_stage, &
    advance_mpi_amr_recursive_film_stage
  use mpi_amr_patch_tree_thermo_3d_mod, only: compute_mpi_amr_tree_cfl_3d
  use mixture_thermo_mod, only: mixture_mass_properties
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use elementary_kinetics_mod, only: elementary_reaction
  use spray_mpi_amr_reactive_mod, only: advance_mpi_amr_spray_sdc
  use spray_mpi_amr_remap_mod, only: remap_mpi_amr_spray
  use spray_mpi_amr_balance_mod, only: sum_mpi_spray_rotating_potential,rebalance_mpi_amr_spray
  use spray_eb_drift_mod, only: drift_spray_level_set
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use transport_database_mod, only: load_h2o2_full_transport
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options,initialize_spray_liquid,parcel_totals,parcel_diameter
  use spray_mpi_amr_tree_mod, only: spray_level_population,redistribute_mpi_amr_level_particles
  use spray_mpi_amr_step_mod, only: spray_sdc_population_stage,advance_mpi_amr_population_sdc_stage
  use spray_mpi_amr_step_mod, only: advance_mpi_eb_population_sdc_stage,convert_mpi_amr_tree_rotating_energy
  use spray_mpi_amr_step_mod, only: kick_mpi_amr_mixture_halo,advance_mpi_amr_level_halo, &
    advance_mpi_amr_recursive_population_stage,advance_mpi_amr_ghost_halo
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d
  use spray_parcel_mod, only: spray_film
  use spray_mixture_mod, only: build_mixture_film,advance_mixture_frozen_gas
  use spray_breakup_mod, only: advance_spray_breakup,khrt_reference_uniforms,construct_reference_etab_children, &
    advance_tab_no_birth_khrt,advance_khrt_shedding,spray_khrt_shedding_state,spray_khrt_release
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquid(1)
  type(spray_options) :: options
  type(spray_parcel) :: original
  type(amr_patch_3d) :: patches(0)
  type(amr_patch_tree_3d) :: tree,callback_mesh
  type(amr_tree_state_3d) :: gas,temp,source,remainder,physical,old_source
  type(spray_level_population) :: groups(3)
  type(spray_sdc_population_stage) :: work
  real(dp), allocatable :: primitive(:),cell(:),boundary(:),old_integral(:),new_integral(:),expected(:)
  real(dp) :: lower(3),spacing(3),wall(5),temperature,sound,interval
  integer :: rank,status,nv,p,i,j,k,count,stage
  logical :: ok,callback_failure=.false.,callback_negative=.false.
  integer :: diffusion_calls=0,reaction_calls=0
  character(len=16) :: faces(6),mode
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call load_h2o2_full_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_full_transport(transport,ok);call require(ok,'transport')
  call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2257000._dp,liquid(1),ok)
  call require(ok,'liquid')
  call initialize_amr_patch_tree_3d([4,4,4],[integer ::],patches,tree,ok);call require(ok,'tree')
  nv=5+size(species);allocate(primitive(nv),cell(nv),boundary(nv),old_integral(nv),new_integral(nv),expected(nv))
  primitive=0;primitive(1)=.4_dp;primitive(5)=101325
  do p=1,size(species)
    if(species(p)%name=='N2') primitive(5+p)=1
  enddo
  call reactive_primitive_to_conserved(species,primitive,cell,temperature,sound,ok);call require(ok,'cell')
  allocate(gas%fields(1),temp%fields(1))
  if(rank==0) then
    allocate(gas%fields(1)%values(nv,0:5,0:5,0:5),temp%fields(1)%values(1,0:5,0:5,0:5))
    do k=0,5
    do j=0,5
    do i=0,5
      gas%fields(1)%values(:,i,j,k)=cell
    enddo
    enddo
    enddo
    temp%fields(1)%values=temperature
  endif
  do p=1,3
    count=0
    if(p==1.and.rank==0) count=1
    allocate(groups(p)%parcels(count),groups(p)%component_mass(1,count),groups(p)%nodes(count))
  enddo
  original%id=1;original%mass=1e-12_dp;original%multiplicity=10;original%temperature=330
  original%position=[.0015_dp,.0015_dp,.0015_dp];original%velocity=[1._dp,0._dp,0._dp]
  if(rank==0) then
    groups(1)%parcels(1)=original;groups(1)%component_mass(1,1)=original%mass;groups(1)%nodes(1)=1
  endif
  options%evaporation=.false.;lower=0;spacing=.001_dp;faces='periodic';interval=2e-8_dp
  call get_command_argument(1,mode)
  if(mode=='held-upper') then
    call check_held_upper()
    call MPI_Finalize(status)
    stop
  endif
  if(mode=='held') then
    call check_held_wrappers()
    call MPI_Finalize(status)
    stop
  endif
  if(mode=='held-topology') then
    call check_rotating_remap(.true.)
    call MPI_Finalize(status)
    stop
  endif
  if(mode=='rotating-remap') then
    call check_rotating_remap()
    call MPI_Finalize(status)
    stop
  endif
  if(mode=='rotating-mol') then
    call check_frame_ledger()
    call check_coupled(.false.,.false.,mol_enabled=.true.,rotating=.true.)
    call check_coupled(.true.,.true.,mol_enabled=.true.,rotating=.true.)
    call check_coupled(.true.,.true.,films_enabled=.true.,impact_enabled=.true.,mol_enabled=.true.,rotating=.true.)
    call check_coupled(.true.,.false.,callbacks=.true.,mol_enabled=.true.,rotating=.true.)
    call MPI_Finalize(status)
    stop
  endif
  if(mode=='mol') then
    call check_coupled(.false.,.false.,mol_enabled=.true.)
    call check_coupled(.true.,.true.,mol_enabled=.true.)
    call check_coupled(.true.,.true.,films_enabled=.true.,impact_enabled=.true.,mol_enabled=.true.)
    call check_coupled(.true.,.false.,callbacks=.true.,mol_enabled=.true.)
    call MPI_Finalize(status)
    stop
  endif
  do stage=0,2
    call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,[0], &
      gas,temp,1,0,1,100,groups,work,0._dp,interval,stage,2,faces,source,remainder,wall,boundary,ok, &
      lower,spacing,physical_source=physical)
    call require(ok,'regular population stage')
    if(stage==0) then
      old_source=source;call integrate(physical,old_integral)
    endif
    if(stage==1.and.rank==0) then
      call require(all(source%fields(1)%values==old_source%fields(1)%values),'intermediate source retained')
    endif
    if(stage<2.and.rank==0) then
      call require(groups(1)%parcels(1)%temperature==original%temperature,'predictor keeps physical particle private')
    endif
  enddo
  call integrate(physical,new_integral);expected=0
  if(rank==0) then
    expected(:5)=parcel_totals(groups(1)%parcels(1),liquid(1))-parcel_totals(original,liquid(1))
    call require(groups(1)%parcels(1)%temperature>original%temperature,'particle heating')
    call require(groups(1)%parcels(1)%position(1)>original%position(1),'particle drift')
    call require(groups(1)%parcels(1)%mass==original%mass,'inactive evaporation preserves mass')
  endif
  call MPI_Allreduce(MPI_IN_PLACE,expected,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
  call require(maxval(abs(.5_dp*interval*(old_integral+new_integral)+expected))<1e-19_dp,'gas particle exchange balance')
  call require(work%next_stage==0,'completed stage lifecycle')
  if(rank==0) print *, 'PASS: regular particle SDC predictor, skipped stage, drift, heating and conservative exchange'
  call check_film_stage()
  call check_breakup_clock(.false.)
  call check_breakup_clock(.true.)
  call check_primitive_halo_sampling()
  call check_tab_history()
  call check_population_flow_cycles()
  call check_flow_birth_rank_crossing()
  call check_coupled(.false.,.false.)
  call check_coupled(.true.,.false.)
  call check_coupled(.true.,.true.)
  call check_coupled(.true.,.false.,.true.)
  call check_coupled(.true.,.true.,.true.)
  call check_coupled(.false.,.false.,films_enabled=.true.)
  call check_coupled(.true.,.false.,films_enabled=.true.)
  call check_coupled(.true.,.true.,films_enabled=.true.)
  call check_coupled(.false.,.false.,films_enabled=.true.,impact_enabled=.true.)
  call check_coupled(.true.,.true.,films_enabled=.true.,impact_enabled=.true.)
  call MPI_Finalize(status)
contains
  subroutine check_held_upper()
    type :: trial
      type(amr_tree_state_3d) :: gas,temp,held
      type(spray_parcel), allocatable :: drops(:)
      real(dp), allocatable :: mass(:,:),film_mass(:,:)
      type(spray_wall_film), allocatable :: films(:)
      real(dp), allocatable :: boundary(:),external(:),exchange(:),outside(:)
      real(dp) :: wall(5),largest,theta,prefix,work_change
      integer(int64) :: watermark
    end type trial
    type(trial) :: run(2),saved
    type(amr_patch_tree_3d) :: mesh
    type(amr_patch_3d) :: box(1)
    type(dual_energy_policy) :: policy
    type(spray_wall_options) :: settings
    type(spray_options) :: local_options
    type(elementary_reaction) :: no_reactions(0)
    integer, allocatable :: owners(:)
    integer :: kind,a,node,d(3),q,nparcel,nfilm,rank_count,c,i,j,k,model
    real(dp) :: y(size(species)),mw,rg,cp,cv,gamma,h,e700,e1000,entropy,vector(nv),energy_gap,work_roundoff
    character(len=16) :: case_faces(6)
    logical :: good,recursive,held_mol
    call MPI_Comm_size(MPI_COMM_WORLD,rank_count,status)
    call require(rank_count==3,'held upper two owners and idle rank')
    y=0
    do c=1,size(species)
      if(species(c)%name=='N2') y(c)=1
    enddo
    call mixture_mass_properties(species,y,700._dp,mw,rg,cp,cv,gamma,h,e700,entropy,good)
    call require(good,'held upper actual700K')
    call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,h,e1000,entropy,good)
    call require(good,'held upper bulk1000K')
    ! Keep UEDEN different during sampling. A total-energy reset before the
    ! callback would otherwise hide a missing actualH forwarding route.
    policy%eta2=2;policy%update_total_from_internal=.false.;policy%allow_negative_energy=.false.
    local_options=options;local_options%evaporation=.false.;settings%evolve=.true.
    settings%enabled(1:2)=.true.;settings%temperature=400
    settings%liquid_conductivity=.6_dp;settings%minimum_height=1e-15_dp
    do model=1,2
    held_mol=model==1
    if(.not.held_mol) then
      policy%allow_small_energy=.false.;policy%allow_negative_energy=.true.
    endif
    do kind=1,5
      recursive=kind==2.or.kind==3
      case_faces=faces
      if(kind==3.or.kind==5) case_faces(1:2)='slip_wall'
      if(allocated(owners)) deallocate(owners)
      if(kind==1) then
        call initialize_amr_patch_tree_3d([4,4,4],[integer ::],patches,mesh,good)
        allocate(owners(1));owners=0
      else
        call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,box(1),good)
        call require(good,'held upper fine child')
        call initialize_amr_patch_tree_3d([4,4,4],[1],box,mesh,good)
        allocate(owners(2));owners=[0,1]
      endif
      call require(good,'held upper topology')
      do a=1,2
        if(.not.allocated(run(a)%boundary)) &
          allocate(run(a)%boundary(nv),run(a)%external(nv),run(a)%exchange(nv),run(a)%outside(nv))
        if(allocated(run(a)%gas%fields)) deallocate(run(a)%gas%fields,run(a)%temp%fields,run(a)%held%fields)
        allocate(run(a)%gas%fields(size(owners)),run(a)%temp%fields(size(owners)),run(a)%held%fields(size(owners)))
        vector=0;vector(1)=.4_dp;vector(6:)=.4_dp*y
        vector(5)=.4_dp*merge(e1000,e700,a==1)
        do node=1,size(owners)
          if(owners(node)/=rank) cycle
          d=mesh%nodes(node)%shape
          allocate(run(a)%gas%fields(node)%values(nv,d(1),d(2),d(3)))
          allocate(run(a)%temp%fields(node)%values(1,d(1),d(2),d(3)),source=700._dp)
          allocate(run(a)%held%fields(node)%values(1,d(1),d(2),d(3)),source=.4_dp*e700)
          do k=1,d(3)
          do j=1,d(2)
          do i=1,d(1)
            run(a)%gas%fields(node)%values(:,i,j,k)=vector
          enddo
          enddo
          enddo
        enddo
        nparcel=0;nfilm=0
        if(rank==merge(0,1,kind==1)) nparcel=1
        if((kind==3.or.kind==5).and.rank==1) nfilm=1
        if(allocated(run(a)%drops)) deallocate(run(a)%drops,run(a)%mass,run(a)%films,run(a)%film_mass)
        allocate(run(a)%drops(nparcel),run(a)%mass(1,nparcel),run(a)%films(nfilm),run(a)%film_mass(1,nfilm))
        if(nparcel>0) then
          run(a)%drops(1)=original
          run(a)%drops(1)%position=[.0019_dp,.0015_dp,.0015_dp]
          run(a)%mass(1,1)=original%mass
        endif
        if(nfilm>0) then
          run(a)%films(1)%mass=1e-12_dp;run(a)%films(1)%temperature=330
          run(a)%films(1)%diameter=1e-4_dp
          run(a)%films(1)%height=run(a)%films(1)%mass/(1000._dp*.25_dp*acos(-1._dp)*1e-8_dp)
          run(a)%films(1)%position=[0._dp,.0015_dp,.0015_dp]
          run(a)%films(1)%normal=[1._dp,0._dp,0._dp];run(a)%film_mass(1,1)=1e-12_dp
        endif
        run(a)%boundary=0;run(a)%external=0;run(a)%exchange=0;run(a)%outside=0;run(a)%wall=0
        run(a)%largest=0;run(a)%theta=1;run(a)%prefix=3.5_dp;run(a)%watermark=2;run(a)%work_change=123
        if(kind==3.or.kind==5) then
          call advance_mpi_amr_spray_sdc(MPI_COMM_WORLD,mesh,owners,species,no_reactions,transport,liquid,local_options, &
            run(a)%gas,run(a)%temp,run(a)%drops,run(a)%mass,lower,spacing,interval,2,.false., &
            [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp,case_faces, &
            run(a)%boundary,run(a)%external,run(a)%exchange,run(a)%wall,run(a)%outside, &
            run(a)%largest,run(a)%theta,run(a)%watermark,good,subcycling=recursive,mol_mode=held_mol, &
            wall_films=run(a)%films,film_mass=run(a)%film_mass,film_settings=settings, &
            internal_energy=run(a)%held,energy_correction=run(a)%prefix,energy_policy=policy, &
            internal_energy_pressure_work_change=run(a)%work_change)
        else
          call advance_mpi_amr_spray_sdc(MPI_COMM_WORLD,mesh,owners,species,no_reactions,transport,liquid,local_options, &
            run(a)%gas,run(a)%temp,run(a)%drops,run(a)%mass,lower,spacing,interval,2,.false., &
            [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp,case_faces, &
            run(a)%boundary,run(a)%external,run(a)%exchange,run(a)%wall,run(a)%outside, &
            run(a)%largest,run(a)%theta,run(a)%watermark,good,subcycling=recursive,mol_mode=held_mol, &
            internal_energy=run(a)%held,energy_correction=run(a)%prefix,energy_policy=policy, &
            internal_energy_pressure_work_change=run(a)%work_change)
        endif
        call require(good,'held regular upper native interval')
      enddo
      if(held_mol) call require(run(1)%work_change==0.and.run(2)%work_change==0,'MOL zero pressure work')
      if(rank==0.and..not.held_mol) print *, 'Godunov pressure work actual-H comparison:',kind, &
        run(1)%work_change,run(2)%work_change
      ! Different bulk energies need not produce bitwise identical restricted
      ! predictors. Scale this additional comparison by p*V*dt/dx*u, not by
      ! the small signed work left after leaf cancellations. Existing field
      ! and inventory gates retain their tolerances.
      work_roundoff=64*epsilon(1._dp)*(.4_dp*rg*700._dp)*product(4*spacing)*interval/minval(spacing)* &
        real(2**(mesh%level_count()-1),dp)*max(1._dp,maxval(abs(original%velocity)))
      call require(ieee_is_finite(run(1)%work_change).and.ieee_is_finite(run(2)%work_change).and. &
        abs(run(1)%work_change-run(2)%work_change)<=work_roundoff.and.run(1)%work_change/=123, &
        'accepted pressure work uses actual H and publishes once')
      call require(size(run(1)%drops)==size(run(2)%drops),'held upper physical parcel count')
      do q=1,size(run(1)%drops)
        call require(abs(run(1)%drops(q)%temperature-run(2)%drops(q)%temperature)<= &
          64*epsilon(1._dp)*max(1._dp,abs(run(2)%drops(q)%temperature)), &
          'actual700K UEINT particle agrees with700K reference despite1000K Eden')
        call require(run(1)%drops(q)%temperature>original%temperature,'held upper parcel heat source reached')
        call require(all(run(1)%drops(q)%position==run(2)%drops(q)%position),'held upper parcel drift agrees')
      enddo
      do q=1,size(run(1)%films)
        call require(abs(run(1)%films(q)%temperature-run(2)%films(q)%temperature)<= &
          64*epsilon(1._dp)*max(1._dp,abs(run(2)%films(q)%temperature)), &
          'actual700K UEINT film agrees with700K reference despite1000K Eden')
        call require(run(1)%films(q)%temperature>330,'held upper film heat source reached')
      enddo
      call require(all(abs(run(1)%exchange-run(2)%exchange)<= &
        64*epsilon(1._dp)*max(1._dp,abs(run(2)%exchange))),'held upper physical gas source agrees')
      call require(run(1)%prefix==3.5_dp.and.run(2)%prefix==3.5_dp,'sampling-only policy preserves numerical prefix')
      energy_gap=0
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        energy_gap=max(energy_gap,maxval(abs(run(1)%gas%fields(node)%values(5,:,:,:)- &
          run(2)%gas%fields(node)%values(5,:,:,:))))
        call require(all(abs(run(1)%held%fields(node)%values-run(2)%held%fields(node)%values)<= &
          64*epsilon(1._dp)*max(1._dp,abs(run(2)%held%fields(node)%values))),'held upper actual scalar endpoint')
        call require(all(abs(run(1)%temp%fields(node)%values-run(2)%temp%fields(node)%values)<= &
          64*epsilon(1._dp)*max(1._dp,abs(run(2)%temp%fields(node)%values))),'held upper actual EOS endpoint')
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,energy_gap,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
      call require(status==MPI_SUCCESS.and.energy_gap>10000,'upper sampling retains distinct bulk energy')
      ! A collective prefix disagreement on the idle rank must reject before
      ! native sources, populations, scalar or physical histories can publish.
      run(1)%work_change=-99._dp
      saved=run(1)
      if(rank==2) run(1)%prefix=run(1)%prefix+1
      call advance_mpi_amr_spray_sdc(MPI_COMM_WORLD,mesh,owners,species,no_reactions,transport,liquid,local_options, &
        run(1)%gas,run(1)%temp,run(1)%drops,run(1)%mass,lower,spacing,interval,2,.false., &
        [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp,case_faces, &
        run(1)%boundary,run(1)%external,run(1)%exchange,run(1)%wall,run(1)%outside, &
        run(1)%largest,run(1)%theta,run(1)%watermark,good,subcycling=recursive,mol_mode=held_mol, &
        internal_energy=run(1)%held,energy_correction=run(1)%prefix,energy_policy=policy, &
        internal_energy_pressure_work_change=run(1)%work_change)
      call require(.not.good,'held upper idle prefix disagreement refuses')
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        call require(all(run(1)%gas%fields(node)%values==saved%gas%fields(node)%values).and. &
          all(run(1)%held%fields(node)%values==saved%held%fields(node)%values).and. &
          all(run(1)%temp%fields(node)%values==saved%temp%fields(node)%values),'held upper prefix field rollback')
      enddo
      call require(all(run(1)%exchange==saved%exchange).and.all(run(1)%boundary==saved%boundary).and. &
        run(1)%watermark==saved%watermark.and.run(1)%work_change==saved%work_change,'held upper prefix budget rollback')
      run(1)=saved
      if(kind==1) then
        ! The initial population fits, but post-gas multiplicity splitting does
        ! not. The final gas interval and work must remain private on refusal.
        local_options%maximum_parcels=1;local_options%maximum_multiplicity=3;local_options%etab_parcel_factor=-1
        local_options%breakup_model="etab";local_options%surface_tension=.072_dp;local_options%liquid_viscosity=.001_dp
        callback_mesh=mesh;diffusion_calls=0
        call advance_mpi_amr_spray_sdc(MPI_COMM_WORLD,mesh,owners,species,no_reactions,transport,liquid,local_options, &
          run(1)%gas,run(1)%temp,run(1)%drops,run(1)%mass,lower,spacing,interval,2,.false., &
          [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp,case_faces, &
          run(1)%boundary,run(1)%external,run(1)%exchange,run(1)%wall,run(1)%outside, &
          run(1)%largest,run(1)%theta,run(1)%watermark,good,subcycling=recursive,mol_mode=held_mol, &
          evaluate_diffusion=custom_diffusion, &
          internal_energy=run(1)%held,energy_correction=run(1)%prefix,energy_policy=policy, &
          internal_energy_pressure_work_change=run(1)%work_change)
        call require(.not.good.and.diffusion_calls>=2,'post-gas split capacity rejects')
        call require(run(1)%work_change==saved%work_change.and.run(1)%prefix==saved%prefix, &
          'post-gas refusal retains pressure work and reset')
        do node=1,size(owners)
          if(owners(node)/=rank) cycle
          call require(all(run(1)%gas%fields(node)%values==saved%gas%fields(node)%values).and. &
            all(run(1)%held%fields(node)%values==saved%held%fields(node)%values).and. &
            all(run(1)%temp%fields(node)%values==saved%temp%fields(node)%values),'post-gas field rollback')
        enddo
        call require(all(run(1)%mass==saved%mass).and.run(1)%watermark==saved%watermark, &
          'post-gas liquid rollback')
        call require(all(run(1)%boundary==saved%boundary).and.all(run(1)%external==saved%external).and. &
          all(run(1)%exchange==saved%exchange).and.all(run(1)%outside==saved%outside).and. &
          all(run(1)%wall==saved%wall),'post-gas physical histories rollback')
        do q=1,size(run(1)%drops)
          call require(all(run(1)%drops(q)%position==saved%drops(q)%position).and. &
            run(1)%drops(q)%temperature==saved%drops(q)%temperature.and. &
            run(1)%drops(q)%multiplicity==saved%drops(q)%multiplicity,'post-gas parcel rollback')
        enddo
        local_options=options;local_options%evaporation=.false.
      endif
      if(rank==0) print *, 'PASS: regular held upper, MOL/model refinement/film case',model,kind
    enddo
    enddo
  end subroutine check_held_upper

  subroutine check_held_wrappers()
    type(amr_tree_state_3d) :: held,bad_held,actual_gas,control_gas,ah,ch,ap,cp,ar,cr,saved_source
    type(amr_tree_state_3d) :: fine_actual,fine_control,fine_temp,fine_held
    type(amr_patch_3d) :: fine_box(1)
    type(amr_patch_tree_3d) :: fine_tree
    type(spray_level_population) :: actual_groups(3),control_groups(3),drops_a,drops_c
    type(spray_level_population) :: levels(1),virtuals(1),ghosts(1)
    type(spray_sdc_population_stage) :: actual_work,control_work,hierarchy_work(1)
    type(spray_level_film_population), allocatable :: films_base(:)
    type(spray_level_film_population) :: films_a(3),films_c(3),film_levels(1),film_virtuals(1),film_ghosts(1)
    type(spray_sdc_film_stage) :: film_work_a,film_work_c,film_hierarchy(1)
    type(spray_wall_film), allocatable :: films(:)
    type(spray_wall_options) :: settings
    real(dp), allocatable :: film_mass(:,:)
    real(dp) :: mw,gas_r,heat_cp,heat_cv,gamma,enthalpy,energy,entropy,held_energy,bulk_energy
    real(dp) :: walls_a(5),walls_c(5),outside_a(nv),outside_c(nv)
    integer :: step,kind,c
    logical :: good
    call mixture_mass_properties(species,cell(6:)/cell(1),700._dp,mw,gas_r,heat_cp,heat_cv,gamma, &
      enthalpy,energy,entropy,good)
    call require(good,'held regular 700K reference energy');held_energy=cell(1)*energy
    call mixture_mass_properties(species,cell(6:)/cell(1),1000._dp,mw,gas_r,heat_cp,heat_cv,gamma, &
      enthalpy,energy,entropy,good)
    call require(good,'held regular 1000K bulk energy');bulk_energy=cell(1)*energy
    actual_gas=gas;control_gas=gas;allocate(held%fields(1))
    if(rank==0) then
      actual_gas%fields(1)%values(5,:,:,:)=bulk_energy
      control_gas%fields(1)%values(5,:,:,:)=held_energy
      allocate(held%fields(1)%values(1,0:5,0:5,0:5),source=held_energy)
    endif
    actual_groups=groups;control_groups=groups
    do step=0,2
      call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,[0], &
        actual_gas,temp,1,0,1,100,actual_groups,actual_work,0._dp,interval,step,2,faces, &
        ah,ar,walls_a,outside_a,good,lower,spacing,physical_source=ap,internal_energy_halo=held)
      call require(good,'regular held particle SDC stage')
      call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,[0], &
        control_gas,temp,1,0,1,100,control_groups,control_work,0._dp,interval,step,2,faces, &
        ch,cr,walls_c,outside_c,good,lower,spacing,physical_source=cp)
      call require(good,'regular particle energy reference')
      call compare_held_fields(ah,ch,'held regular particle source')
      call compare_held_fields(ar,cr,'held regular particle remainder')
      call compare_held_fields(ap,cp,'held regular particle physical export')
      call require(all(walls_a==walls_c).and.all(outside_a==outside_c),'held regular particle ledgers')
      if(step/=0) cycle
      saved_source=ah;bad_held=held
      if(rank==0) bad_held%fields(1)%values(1,0,0,0)=ieee_value(1._dp,ieee_quiet_nan)
      call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,[0], &
        actual_gas,temp,1,0,1,100,actual_groups,actual_work,0._dp,interval,1,2,faces, &
        ah,ar,walls_a,outside_a,good,lower,spacing,physical_source=ap,internal_energy_halo=bad_held)
      call require(.not.good.and.actual_work%next_stage==1,'held particle NaN rollback')
      call compare_held_fields(ah,saved_source,'held particle NaN retains source')
      call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,[0], &
        actual_gas,temp,1,0,1,100,actual_groups,actual_work,0._dp,interval,1,2,faces, &
        ah,ar,walls_a,outside_a,good,lower,spacing,physical_source=ap)
      call require(.not.good.and.actual_work%next_stage==1,'held particle cached mode rollback')
      call compare_held_fields(ah,saved_source,'held particle mode retains source')
      if(rank==2) then
        call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,[0], &
          actual_gas,temp,1,0,1,100,actual_groups,actual_work,0._dp,interval,1,2,faces, &
          ah,ar,walls_a,outside_a,good,lower,spacing,physical_source=ap)
      else
        call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,[0], &
          actual_gas,temp,1,0,1,100,actual_groups,actual_work,0._dp,interval,1,2,faces, &
          ah,ar,walls_a,outside_a,good,lower,spacing,physical_source=ap,internal_energy_halo=held)
      endif
      call require(.not.good.and.actual_work%next_stage==1,'held particle idle presence disagreement')
    enddo
    if(rank==0) then
      call require(actual_groups(1)%parcels(1)%temperature==control_groups(1)%parcels(1)%temperature, &
        'held regular particle endpoint matches 700K gas')
      call require(actual_gas%fields(1)%values(5,1,1,1)==bulk_energy,'particle sampling retains actual Eden')
    endif
    levels(1)=groups(1);virtuals(1)=groups(2);ghosts(1)=groups(3)
    do step=0,2
      call advance_mpi_amr_recursive_population_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,[0], &
        actual_gas,temp,1,0,1,100,levels,virtuals,ghosts,hierarchy_work,0._dp,interval,step,2,faces, &
        ah,ar,walls_a,outside_a,good,lower,spacing,physical_source=ap,internal_energy_halo=held)
      call require(good,'held regular recursive particle stage')
    enddo
    call compare_held_fields(ap,cp,'held recursive particle physical export')
    if(rank==0) call require(levels(1)%parcels(1)%temperature==control_groups(1)%parcels(1)%temperature, &
      'held recursive particle endpoint')
    call initialize_amr_patch_3d(4,4,4,1,2,1,2,1,2,2,fine_box(1),good)
    call require(good,'held ghost fine patch')
    call initialize_amr_patch_tree_3d([4,4,4],[1],fine_box,fine_tree,good)
    call require(good,'held ghost fine tree')
    allocate(fine_actual%fields(2),fine_control%fields(2),fine_temp%fields(2),fine_held%fields(2))
    do c=1,2
      fine_actual%fields(c)=actual_gas%fields(1);fine_control%fields(c)=control_gas%fields(1)
      fine_temp%fields(c)=temp%fields(1);fine_held%fields(c)=held%fields(1)
    enddo
    do kind=1,2
      drops_a=groups(1);drops_c=groups(1)
      if(kind==1) then
        call advance_mpi_amr_level_halo(MPI_COMM_WORLD,species,transport,liquid,options,tree,[0],actual_gas,temp, &
          1,0,1,100,drops_a%nodes,drops_a%parcels,drops_a%component_mass,lower,spacing, &
          [.true.,.true.,.true.],interval,0._dp,faces,ah,ar,walls_a,outside_a,good,internal_energy_halo=held)
      else
        if(rank==0) then
          drops_a%nodes=2;drops_c%nodes=2
        endif
        call advance_mpi_amr_ghost_halo(MPI_COMM_WORLD,species,transport,liquid,options,fine_tree,[0,0], &
          fine_actual,fine_temp,1,1,1,100,drops_a%nodes,drops_a%parcels,drops_a%component_mass,lower,spacing, &
          [.true.,.true.,.true.],interval,0._dp,faces,ah,ar,walls_a,outside_a,good,internal_energy_halo=fine_held)
      endif
      call require(good,'held regular level/ghost convenience')
      if(kind==1) then
        call advance_mpi_amr_level_halo(MPI_COMM_WORLD,species,transport,liquid,options,tree,[0],control_gas,temp, &
          1,0,1,100,drops_c%nodes,drops_c%parcels,drops_c%component_mass,lower,spacing, &
          [.true.,.true.,.true.],interval,0._dp,faces,ch,cr,walls_c,outside_c,good)
      else
        call advance_mpi_amr_ghost_halo(MPI_COMM_WORLD,species,transport,liquid,options,fine_tree,[0,0], &
          fine_control,fine_temp,1,1,1,100,drops_c%nodes,drops_c%parcels,drops_c%component_mass,lower,spacing, &
          [.true.,.true.,.true.],interval,0._dp,faces,ch,cr,walls_c,outside_c,good)
      endif
      call require(good,'regular level/ghost physical reference')
      if(rank==0) call require(drops_a%parcels(1)%temperature==drops_c%parcels(1)%temperature, &
        'held regular level/ghost temperature')
    enddo
    settings%enabled=.true.;settings%evolve=.true.;settings%temperature=400
    settings%liquid_conductivity=.6_dp;settings%minimum_height=1e-15_dp
    c=merge(1,0,rank==0);allocate(films(c),film_mass(1,c))
    if(c==1) then
      films(1)=spray_wall_film([0._dp,.0015_dp,.0015_dp],[1._dp,0._dp,0._dp], &
        1e-10_dp,330._dp,1e-4_dp,1e-5_dp,1)
      film_mass(1,1)=films(1)%mass
    endif
    call initialize_mpi_amr_film_levels(MPI_COMM_WORLD,tree,[0],films,film_mass,lower,spacing,100,films_base,good)
    call require(good,'held regular film initialization')
    films_a(1)=films_base(1);films_a(2)=films_base(1)
    allocate(films_a(3)%films(0),films_a(3)%component_mass(1,0),films_a(3)%nodes(0),films_a(3)%probes(3,0))
    films_c=films_a
    do step=0,2
      call advance_mpi_amr_film_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,settings,tree,[0], &
        actual_gas,temp,1,0,films_a,film_work_a,lower,spacing,[.false.,.false.,.false.],0._dp,interval,step,2, &
        ah,ar,walls_a,good,ap,internal_energy_halo=held)
      call require(good,'held regular film SDC')
      call advance_mpi_amr_film_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,settings,tree,[0], &
        control_gas,temp,1,0,films_c,film_work_c,lower,spacing,[.false.,.false.,.false.],0._dp,interval,step,2, &
        ch,cr,walls_c,good,cp)
      call require(good,'regular film energy reference')
      call compare_held_fields(ah,ch,'held regular film source')
      call compare_held_fields(ar,cr,'held regular film remainder')
      call require(all(walls_a==walls_c),'held regular film wall ledger')
    enddo
    film_levels=films_base;film_virtuals(1)=films_a(3);film_ghosts(1)=films_a(3)
    do step=0,2
      call advance_mpi_amr_recursive_film_stage(MPI_COMM_WORLD,species,transport,liquid,options,settings,tree,[0], &
        actual_gas,temp,1,0,1,100,film_levels,film_virtuals,film_ghosts,film_hierarchy,lower,spacing, &
        0._dp,interval,step,2,faces,ah,ar,walls_a,good,ap,internal_energy_halo=held)
      call require(good,'held regular recursive film stage')
    enddo
    call compare_held_fields(ap,cp,'held recursive film physical export')
    if(rank==0) call require(film_levels(1)%films(1)%temperature==films_c(1)%films(1)%temperature, &
      'held regular recursive film endpoint')
    if(rank==0) print *, 'PASS: regular UEINT wrappers, 700K/1000K particle-film sampling and mode/NaN rollback'
  end subroutine check_held_wrappers
  subroutine compare_held_fields(actual,expected,label)
    type(amr_tree_state_3d), intent(in) :: actual,expected
    character(len=*), intent(in) :: label
    integer :: node
    call require(allocated(actual%fields).and.allocated(expected%fields),label)
    call require(size(actual%fields)==size(expected%fields),label)
    do node=1,size(actual%fields)
      call require(allocated(actual%fields(node)%values).eqv.allocated(expected%fields(node)%values),label)
      if(allocated(actual%fields(node)%values)) &
        call require(all(actual%fields(node)%values==expected%fields(node)%values),label)
    enddo
  end subroutine compare_held_fields
  subroutine check_frame_ledger()
    type(rotating_frame_3d) :: frame,bad_frame
    type(amr_patch_tree_3d) :: mesh
    type(amr_patch_3d) :: child(1)
    type(amr_tree_state_3d) :: field,saved
    type(spray_parcel), allocatable :: drops(:)
    type(spray_wall_film), allocatable :: films(:)
    real(dp) :: analytic,local,ledger,inverse,sentinel,h(3),position(3),origin(3),rho,phi_value
    real(dp) :: before,after,expected_before,expected_after,outside(5),walls(5),outflow_conversion
    real(dp) :: nodal(0:1,0:1,0:1),duration,ordinary_energy,exit_position(3)
    real(dp), allocatable :: energies(:)
    character(len=16) :: boundaries(6)
    integer :: owners(2),node,x,y,z,c,n(3),hits
    logical :: good
    frame%axis=3;frame%angular_speed=200;frame%axis_point=[.002_dp,.002_dp,0._dp]
    call initialize_amr_patch_3d(4,4,4,1,1,1,2,2,2,2,child(1),good)
    call require(good,'coordinate ledger child')
    call initialize_amr_patch_tree_3d([4,4,4],[1],child,mesh,good)
    call require(good,'coordinate ledger tree')
    owners=[0,1];allocate(field%fields(2));local=0
    do node=1,2
      if(owners(node)/=rank) cycle
      n=mesh%nodes(node)%shape;h=spacing;origin=lower
      if(node==2) h=.5_dp*spacing
      allocate(field%fields(node)%values(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1),source=0._dp)
      ! Large halo density is deliberate: ghosts must not enter the ledger.
      field%fields(node)%values(1,:,:,:)=1e3_dp
      field%fields(node)%values(5,:,:,:)=1e6_dp
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        rho=(-1._dp)**x*real(x+y+z,dp)
        field%fields(node)%values(1,x,y,z)=rho
        if(covered(mesh,node,[x,y,z])) cycle
        position=origin+h*(real([x,y,z],dp)-.5_dp)
        phi_value=.5_dp*frame%angular_speed**2*sum((position(1:2)-frame%axis_point(1:2))**2)
        local=local-rho*phi_value*product(h)
      enddo
      enddo
      enddo
    enddo
    call MPI_Allreduce(local,analytic,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    saved=field;ledger=19
    call convert_mpi_amr_tree_rotating_energy(MPI_COMM_WORLD,mesh,owners,field,lower,spacing,frame,.true.,good, &
      halo_width=1,energy_change=ledger)
    call require(good.and.abs(ledger-analytic)<2e-23_dp,'signed source conversion excludes covered cells and halos')
    inverse=19
    call convert_mpi_amr_tree_rotating_energy(MPI_COMM_WORLD,mesh,owners,field,lower,spacing,frame,.false.,good, &
      halo_width=1,energy_change=inverse)
    call require(good.and.abs(inverse+ledger)<2e-23_dp,'inverse coordinate ledger')
    bad_frame=frame
    if(rank==2) bad_frame%angular_speed=201
    sentinel=31;field=saved
    call convert_mpi_amr_tree_rotating_energy(MPI_COMM_WORLD,mesh,owners,field,lower,spacing,bad_frame,.true.,good, &
      halo_width=1,energy_change=sentinel)
    call require(.not.good.and.sentinel==31,'coordinate conversion ledger rollback')
    do node=1,2
      if(owners(node)==rank) call require(all(field%fields(node)%values==saved%fields(node)%values), &
        'coordinate conversion field rollback')
    enddo
    field=saved
    if(rank==0) then
      field%fields(1)%values(1,1,1,1)=huge(1._dp)
      field%fields(1)%values(5,1,1,1)=huge(1._dp)
    endif
    saved=field;sentinel=31
    call convert_mpi_amr_tree_rotating_energy(MPI_COMM_WORLD,mesh,owners,field,lower,spacing,frame,.false.,good, &
      halo_width=1,energy_change=sentinel)
    call require(.not.good.and.sentinel==31,'overflow rejects before floating exception and ledger publication')
    do node=1,2
      if(owners(node)==rank) call require(all(field%fields(node)%values==saved%fields(node)%values), &
        'overflow coordinate field rollback')
    enddo
    c=merge(1,0,rank==0);allocate(drops(c),films(c),energies(c));energies=0
    if(c==1) then
      drops(1)=original
      films(1)=spray_wall_film([0._dp,.0015_dp,.0015_dp],[1._dp,0._dp,0._dp], &
        1e-10_dp,330._dp,1e-4_dp,1e-5_dp,1)
    endif
    expected_before=original%mass*original%multiplicity*.5_dp*200._dp**2* &
      sum((original%position(1:2)-frame%axis_point(1:2))**2)+ &
      1e-10_dp*.5_dp*200._dp**2*((0._dp-.002_dp)**2+(.0015_dp-.002_dp)**2)
    before=31
    call sum_mpi_spray_rotating_potential(MPI_COMM_WORLD,frame,drops,before,good,films)
    call require(good.and.abs(before-expected_before)<1e-26_dp,'owned particle and film coordinate inventory')
    nodal=1;boundaries='slip_wall';duration=2e-5_dp;outflow_conversion=31
    call drift_spray_level_set(nodal,lower,lower+4*spacing,duration,32,drops,walls,hits,good, &
      boundaries,energies,outside,frame=frame,energy_conversion=outflow_conversion)
    call require(good.and.all(walls==0).and.all(outside==0).and.outflow_conversion==0,'pure drift has no exterior work')
    if(c==1) then
      call require(all(drops(1)%velocity==original%velocity), 'rotation does not add particle acceleration')
      call require(all(drops(1)%position==original%position+duration*original%velocity),'native particle drift preserved')
    endif
    expected_after=expected_before+original%mass*original%multiplicity*.5_dp*200._dp**2* &
      ((original%position(1)+duration*original%velocity(1)-.002_dp)**2-(original%position(1)-.002_dp)**2)
    after=31
    call sum_mpi_spray_rotating_potential(MPI_COMM_WORLD,frame,drops,after,good,films)
    call require(good.and.abs(after-expected_after)<1e-26_dp,'pure drift coordinate change from analytic position')
    if(c==1) then
      drops(1)=original;drops(1)%position=[.0035_dp,.0015_dp,.0015_dp];drops(1)%velocity=[100._dp,0._dp,0._dp]
    endif
    boundaries(2)='outflow';duration=1e-4_dp;outflow_conversion=31
    call drift_spray_level_set(nodal,lower,lower+4*spacing,duration,32,drops,walls,hits,good, &
      boundaries,energies,outside,frame=frame,energy_conversion=outflow_conversion)
    call require(good,'rotating exterior coordinate conversion')
    exit_position=[.004_dp,.0015_dp,.0015_dp]
    phi_value=.5_dp*200._dp**2*sum((exit_position(1:2)-frame%axis_point(1:2))**2)
    if(c==1) then
      ordinary_energy=.5_dp*original%mass*original%multiplicity*100._dp**2
      call require(abs(outflow_conversion-original%mass*original%multiplicity*phi_value)<1e-26_dp, &
        'outflow coordinate energy uses exit point')
      call require(abs(outside(5)+ordinary_energy-outflow_conversion)<1e-23_dp, &
        'outflow ordinary and rotating storage energies')
      call require(drops(1)%mass==0,'outflow retires physical particle')
    else
      call require(outflow_conversion==0.and.all(outside==0),'idle outflow rank does not duplicate ledger')
    endif
    if(rank==0) print *, 'PASS: analytic rotating source, physical liquid, drift, exit and atomic rollback ledgers'
  end subroutine check_frame_ledger

  subroutine check_rotating_remap(held_mode)
    logical, optional, intent(in) :: held_mode
    logical :: retained,moved
    type(amr_tree_state_3d), allocatable :: held,saved_h
    type(amr_tree_state_3d) :: ordinary,reference
    real(dp) :: fractions(size(species)),mw,rg,cp,cv,gamma,enthalpy,e700,e1000,entropy,dt_a,dt_b
    real(dp), allocatable :: costs(:)
    integer :: parcel_count
    type(amr_patch_tree_3d) :: mesh,refined
    type(amr_patch_3d) :: child(1)
    type(amr_tree_state_3d) :: u,t,heat,rates,saved_u,saved_t
    type(rotating_frame_3d) :: frame,bad_frame
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    integer, allocatable :: owners(:)
    real(dp) :: before(nv),after(nv),volume
    integer :: node,d(3),x,y,z,ancestor
    retained=.false.
    if(present(held_mode)) retained=held_mode
    if(retained) then
      fractions=cell(6:)/cell(1)
      call mixture_mass_properties(species,fractions,700._dp,mw,rg,cp,cv,gamma,enthalpy,e700,entropy,ok)
      call require(ok,'topology actual700K')
      call mixture_mass_properties(species,fractions,1000._dp,mw,rg,cp,cv,gamma,enthalpy,e1000,entropy,ok)
      call require(ok,'topology bulk1000K')
      allocate(held,saved_h);allocate(held%fields(2))
    endif
    call initialize_amr_patch_3d(4,4,4,1,1,2,2,2,2,2,child(1),ok)
    call require(ok,'rotating remap patch')
    call initialize_amr_patch_tree_3d([4,4,4],[1],child,mesh,ok)
    call require(ok,'rotating remap tree')
    refined=mesh
    owners=[0,1]
    allocate(u%fields(2),t%fields(2),heat%fields(2),rates%fields(2))
    do node=1,2
      if(owners(node)/=rank) cycle
      d=mesh%nodes(node)%shape
      allocate(u%fields(node)%values(nv,d(1),d(2),d(3)),t%fields(node)%values(1,d(1),d(2),d(3)), &
        heat%fields(node)%values(1,d(1),d(2),d(3)),rates%fields(node)%values(nv-4,d(1),d(2),d(3)))
      do z=1,d(3)
      do y=1,d(2)
      do x=1,d(1)
        u%fields(node)%values(:,x,y,z)=cell
      enddo
      enddo
      enddo
      t%fields(node)%values=temperature;heat%fields(node)%values=23;rates%fields(node)%values=2
      if(retained) then
        allocate(held%fields(node)%values(1,d(1),d(2),d(3)),source=cell(1)*e700)
        t%fields(node)%values=700
        u%fields(node)%values(5,:,:,:)=cell(1)*e1000
      endif
    enddo
    allocate(drops(merge(1,0,rank==1)),masses(1,merge(1,0,rank==1)))
    if(rank==1) then
      drops(1)=original;drops(1)%position=[.0005_dp,.0015_dp,.0015_dp];masses(1,1)=original%mass
    endif
    frame%axis=3;frame%angular_speed=1000;frame%axis_point=.001_dp
    call convert_mpi_amr_tree_rotating_energy(MPI_COMM_WORLD,mesh,owners,u,lower,spacing,frame,.true.,ok)
    call require(ok,'rotating remap initial conversion')
    saved_u=u;saved_t=t;bad_frame=frame
    if(retained) saved_h=held
    if(rank==2) bad_frame%axis_point(1)=bad_frame%axis_point(1)+.001_dp
    call remap_mpi_amr_spray(MPI_COMM_WORLD,species,mesh,owners,u,t,drops,masses,tree,[2],lower,spacing, &
      [.false.,.false.,.false.],100,ok,reaction_source=rates,heat_release=heat,frame=bad_frame,internal_energy=held)
    call require(.not.ok.and.size(mesh%nodes)==2.and.all(owners==[0,1]),'rotating remap frame rollback')
    do node=1,2
      if(owners(node)/=rank) cycle
      call require(all(u%fields(node)%values==saved_u%fields(node)%values),'rotating remap gas rollback')
      call require(all(t%fields(node)%values==saved_t%fields(node)%values),'rotating remap temperature rollback')
      if(retained) call require(all(held%fields(node)%values==saved_h%fields(node)%values), &
        'actualH remap collective frame rollback')
    enddo
    before=0
    do node=1,2
      if(owners(node)/=rank) cycle
      d=mesh%nodes(node)%shape;volume=product(spacing);ancestor=node
      do while(ancestor>1)
        volume=volume/real(mesh%nodes(ancestor)%patch%refinement_ratio,dp)**3
        ancestor=mesh%nodes(ancestor)%parent
      enddo
      do z=1,d(3)
      do y=1,d(2)
      do x=1,d(1)
        if(.not.covered(mesh,node,[x,y,z])) before=before+volume*u%fields(node)%values(:,x,y,z)
      enddo
      enddo
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,before,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call remap_mpi_amr_spray(MPI_COMM_WORLD,species,mesh,owners,u,t,drops,masses,tree,[2],lower,spacing, &
      [.false.,.false.,.false.],100,ok,reaction_source=rates,heat_release=heat,frame=frame,internal_energy=held)
    call require(ok.and.size(mesh%nodes)==1.and.all(owners==[2]),'rotating coarsen and owner migration')
    after=0
    if(rank==2) then
      do x=1,nv
        after(x)=product(spacing)*sum(u%fields(1)%values(x,:,:,:))
      enddo
      call require(size(drops)==1.and.all(masses==original%mass),'rotating remap particle inventory')
      if(retained) then
        call require(all(abs(held%fields(1)%values-cell(1)*e700)<=64*epsilon(1._dp)*abs(cell(1)*e700)), &
          'rotating remap retains independent actualH')
        call require(all(abs(t%fields(1)%values-700)<=64*epsilon(1._dp)*700), &
          'rotating remap EOS uses actual700K instead of bulk1000K')
      endif
      call require(all(heat%fields(1)%values==23).and.all(rates%fields(1)%values==2),'rotating remap stage history')
      call require(all(ieee_is_finite(t%fields(1)%values)).and.all(t%fields(1)%values>0),'rotating remap EOS recovery')
    else
      call require(size(drops)==0,'rotating remap empty owners')
    endif
    call MPI_Allreduce(MPI_IN_PLACE,after,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call require(maxval(abs(after-before))<1e-16_dp,'rotating conserved energy and species remap')
    if(rank==0) print *, 'PASS: rotating conserved remap, owner migration, stage history and frame rollback'
    if(retained) then
      ordinary=u
      call convert_mpi_amr_tree_rotating_energy(MPI_COMM_WORLD,mesh,owners,ordinary,lower,spacing,frame,.false.,ok)
      call require(ok,'actualH topology CFL ordinary energy')
      reference=ordinary
      if(rank==2) reference%fields(1)%values(5,:,:,:)=held%fields(1)%values(1,:,:,:)+ &
        .5_dp*sum(reference%fields(1)%values(2:4,:,:,:)**2,dim=1)/reference%fields(1)%values(1,:,:,:)
      call compute_mpi_amr_tree_cfl_3d(mesh,owners,species,ordinary,t,spacing,.2_dp,dt_a,MPI_COMM_WORLD,ok,held)
      call require(ok,'actualH topology CFL')
      call compute_mpi_amr_tree_cfl_3d(mesh,owners,species,reference,t,spacing,.2_dp,dt_b,MPI_COMM_WORLD,ok,held)
      call require(ok.and.dt_a==dt_b,'CFL actualH sound speed independent of distinct bulkenergy')
      ! Two owned patches create a strict positive work gain; one patch
      ! cannot improve the maximum load by changing ranks.
      call remap_mpi_amr_spray(MPI_COMM_WORLD,species,mesh,owners,u,t,drops,masses,refined,[2,2],lower,spacing, &
        [.false.,.false.,.false.],100,ok,reaction_source=rates,heat_release=heat,frame=frame,internal_energy=held)
      call require(ok,'actualH rebalance two-patch preparation')
      allocate(costs(2));costs=0
      if(rank==2) costs=1
      call rebalance_mpi_amr_spray(MPI_COMM_WORLD,mesh,owners,nv,u,t,drops,masses,lower,spacing, &
        [.false.,.false.,.false.],costs,0._dp,100,moved,ok,reaction_source=rates,heat_release=heat,internal_energy=held)
      call require(ok.and.moved.and.all(owners==[0,1]),'actualH topology measured-cost owner migration')
      parcel_count=size(drops)
      call MPI_Allreduce(MPI_IN_PLACE,parcel_count,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
      call require(status==MPI_SUCCESS.and.parcel_count==1.and.all(masses==original%mass), &
        'actualH rebalance physicalparcel')
      do node=1,2
        if(owners(node)==rank) then
          call require(all(abs(held%fields(node)%values-cell(1)*e700)<=64*epsilon(1._dp)*abs(cell(1)*e700)), &
            'actualH rebalance independent scalar')
          call require(all(abs(t%fields(node)%values-700)<=64*epsilon(1._dp)*700).and. &
            all(heat%fields(node)%values==23).and.all(rates%fields(node)%values==2),'actualH rebalance histories')
        else
          call require(.not.allocated(held%fields(node)%values),'actualH rebalance idlefields')
        endif
      enddo
      if(rank==0) print *, 'PASS: actualH rotating remap, all-history rebalance and EOS/CFL independent of bulkenergy'
    endif
  end subroutine check_rotating_remap
  subroutine check_flow_birth_rank_crossing()
    type(amr_patch_3d) :: child(2)
    type(amr_patch_tree_3d) :: mesh
    type(amr_tree_state_3d) :: halo,halo_t,feedback,leftover,raw,saved_feedback
    type(spray_level_population) :: population(3),initial_population(3)
    type(spray_sdc_population_stage) :: stage_work
    type(spray_parcel), allocatable :: manual(:)
    type(spray_options) :: settings
    real(dp), allocatable :: mass(:,:),history(:)
    integer, allocatable :: nodes(:)
    integer :: owners(3),node,n(3),np,g,slice,x,y,z
    integer(int64) :: marker,manual_marker
    real(dp) :: duration,walls(5),outside(nv)
    call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,child(1),ok)
    call require(ok,'flow crossing left patch')
    call initialize_amr_patch_3d(4,4,4,3,4,1,4,1,4,2,child(2),ok)
    call require(ok,'flow crossing right patch')
    call initialize_amr_patch_tree_3d([4,4,4],[1,1],child,mesh,ok)
    call require(ok,'flow crossing tree')
    owners=[0,0,1];allocate(halo%fields(3),halo_t%fields(3))
    do node=1,3
      if(owners(node)/=rank) cycle
      n=mesh%nodes(node)%shape
      allocate(halo%fields(node)%values(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1), &
        halo_t%fields(node)%values(1,0:n(1)+1,0:n(2)+1,0:n(3)+1))
      do z=0,n(3)+1
      do y=0,n(2)+1
      do x=0,n(1)+1
        halo%fields(node)%values(:,x,y,z)=cell
      enddo
      enddo
      enddo
      halo_t%fields(node)%values=temperature
    enddo
    do g=1,3
      np=merge(1,0,g==1.and.rank==0)
      allocate(population(g)%parcels(np),population(g)%component_mass(1,np),population(g)%nodes(np))
      if(np>0) then
        population(g)%parcels(1)=original;population(g)%parcels(1)%id=37
        population(g)%parcels(1)%mass=1e-9_dp;population(g)%parcels(1)%multiplicity=.01_dp
        population(g)%parcels(1)%position=[.00175_dp,.00175_dp,.00175_dp]
        population(g)%parcels(1)%velocity=[100._dp,0._dp,0._dp];population(g)%parcels(1)%deformation=1.1_dp
        population(g)%component_mass(1,1)=1e-9_dp;population(g)%nodes=2
      endif
    enddo
    initial_population=population
    settings=options;settings%drag=.false.;settings%heat=.false.;settings%evaporation=.false.
    settings%breakup_model='etab';settings%surface_tension=.072_dp;settings%liquid_viscosity=.001_dp
    settings%etab_parcel_factor=1;settings%maximum_multiplicity=.001_dp
    manual=population(1)%parcels;mass=population(1)%component_mass;nodes=population(1)%nodes
    allocate(history(size(manual)),source=0._dp);manual_marker=37;duration=6e-6_dp
    do slice=1,3
      call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,settings,mesh,owners,halo,halo_t, &
        1,nodes,manual,mass,lower,spacing,[.true.,.true.,.true.],duration/6,raw,ok, &
        drift_interval=duration/3,boundary_faces=faces,wall_change=walls,boundary_change=outside, &
        last_id=manual_marker,tab_interval=duration/3,tab_tangent_history=history, &
        last_flow_cycle=slice==3,khrt_interval=.5_dp*duration)
      call require(ok,'explicit crossing slice with TAB history')
      if(slice<3) call require(manual_marker==37,'ETAB birth deferred across crossing slices')
      call redistribute_mpi_amr_level_particles(MPI_COMM_WORLD,mesh,owners,manual,mass,nodes,lower,spacing, &
        [.true.,.true.,.true.],1,1,100,ok,tab_tangent_history=history)
      call require(ok,'intermediate rank crossing with TAB history')
      if(slice==2) then
        call require(size(manual)==merge(1,0,rank==1),'parent migrated before final birth')
        if(rank==1) call require(history(1)>0,'new rank receives accumulated TAB tangent')
      endif
    enddo
    call require(manual_marker==38,'final ETAB birth after rank crossing')
    call require(all(history==0),'replacement starts fresh TAB history')
    marker=37
    call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,settings,mesh,owners, &
      halo,halo_t,1,1,1,100,population,stage_work,0._dp,duration,0,2,faces,feedback,leftover,walls,outside,ok, &
      lower,spacing,last_id=marker,flow_cfl=1.2_dp)
    call require(ok.and.marker==38,'SDC crossing publishes final birth watermark')
    call require(size(stage_work%midpoint(1)%parcels)==size(manual),'SDC crossing ownership matches slices')
    if(size(manual)>0) then
      call require(all(stage_work%midpoint(1)%parcels%id==manual%id),'crossing replacement identity')
      call require(maxval(abs(stage_work%midpoint(1)%parcels(1)%position-manual(1)%position))<1e-15_dp.and. &
        maxval(abs(stage_work%midpoint(1)%parcels(1)%velocity-manual(1)%velocity))<1e-12_dp, &
        'crossing TAB cumulative release matches explicit slices')
    endif
    saved_feedback=feedback;population=initial_population;stage_work=spray_sdc_population_stage();marker=huge(marker)
    call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,settings,mesh,owners, &
      halo,halo_t,1,1,1,100,population,stage_work,0._dp,duration,0,2,faces,feedback,leftover,walls,outside,ok, &
      lower,spacing,last_id=marker,flow_cfl=1.2_dp)
    call require(.not.ok.and.marker==huge(marker).and.stage_work%next_stage==0,'late birth failure rejects entire interval')
    call require(size(population(1)%parcels)==size(initial_population(1)%parcels),'late failure preserves original owner')
    if(rank==0) call require(all(population(1)%parcels(1)%position==initial_population(1)%parcels(1)%position), &
      'late failure preserves pre-crossing particle')
    do node=1,3
      if(.not.allocated(feedback%fields(node)%values)) cycle
      call require(all(feedback%fields(node)%values==saved_feedback%fields(node)%values),'late failure preserves gas source')
    enddo
  end subroutine check_flow_birth_rank_crossing
  subroutine check_population_flow_cycles()
    type(spray_parcel), allocatable :: automatic(:),manual(:)
    real(dp), allocatable :: mass(:,:),reference_mass(:,:)
    integer, allocatable :: nodes(:),reference_nodes(:)
    type(spray_options) :: settings
    type(spray_level_population) :: population(3)
    type(spray_sdc_population_stage) :: stage_work
    type(amr_tree_state_3d) :: feedback,leftover,raw,delta,delta_left,total
    real(dp) :: wall_auto(5),wall_ref(5),exterior_auto(nv),exterior_ref(nv),full_dt
    integer :: slice,n,g
    settings=options;settings%mixture_time_scheme='euler';settings%heat=.false.;settings%evaporation=.false.
    n=merge(1,0,rank==0)
    allocate(automatic(n),mass(1,n),nodes(n))
    if(n>0) then
      automatic(1)=original;automatic(1)%velocity=[1000._dp,0._dp,0._dp]
      mass(1,1)=automatic(1)%mass;nodes=1
    endif
    manual=automatic;reference_mass=mass;reference_nodes=nodes
    full_dt=1.2e-6_dp
    call advance_mpi_amr_level_halo(MPI_COMM_WORLD,species,transport,liquid,settings,tree,[0],gas,temp, &
      1,0,1,100,nodes,automatic,mass,lower,spacing,[.true.,.true.,.true.],.5_dp*full_dt,full_dt,faces, &
      feedback,leftover,wall_auto,exterior_auto,ok,flow_cfl=1.2_dp)
    call require(ok,'automatic CFL population flow cycles')
    do slice=1,3
      call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,settings,tree,[0],gas,temp, &
        1,reference_nodes,manual,reference_mass,lower,spacing,[.true.,.true.,.true.],full_dt/6,raw,ok, &
        drift_interval=full_dt/3,boundary_faces=faces,wall_change=wall_ref,boundary_change=exterior_ref, &
        target_level=0,level_source=delta,source_remainder=delta_left,route_to_coarse=.false., &
        tab_interval=full_dt/6,last_flow_cycle=slice==3,khrt_interval=.5_dp*full_dt)
      call require(ok,'manual CFL slice')
      if(slice==1) then
        total=delta
      else if(rank==0) then
        total%fields(1)%values=total%fields(1)%values+delta%fields(1)%values
      endif
    enddo
    if(rank==0) then
      call require(maxval(abs(automatic(1)%position-manual(1)%position))<1e-15_dp,'flow cycle trajectory')
      call require(maxval(abs(automatic(1)%velocity-manual(1)%velocity))<1e-12_dp,'flow cycle transfer')
      call require(all(mass==reference_mass),'flow cycle liquid inventory')
      call require(maxval(abs(feedback%fields(1)%values-total%fields(1)%values))<1e-15_dp, &
        'all flow cycle sources accumulated')
    endif
    if(rank==0) then
      automatic(1)=original;automatic(1)%velocity=[1000._dp,0._dp,0._dp]
      mass(1,1)=automatic(1)%mass
    endif
    manual=automatic;reference_mass=mass
    call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,settings,tree,[0],gas,temp, &
      1,reference_nodes,manual,reference_mass,lower,spacing,[.true.,.true.,.true.],.5_dp*full_dt,raw,ok, &
      drift_interval=full_dt,boundary_faces=faces,wall_change=wall_ref,boundary_change=exterior_ref, &
      target_level=0,level_source=delta,source_remainder=delta_left,route_to_coarse=.false.)
    call require(ok,'single native default reference')
    call advance_mpi_amr_level_halo(MPI_COMM_WORLD,species,transport,liquid,settings,tree,[0],gas,temp, &
      1,0,1,100,nodes,automatic,mass,lower,spacing,[.true.,.true.,.true.],.5_dp*full_dt,full_dt,faces, &
      feedback,leftover,wall_auto,exterior_auto,ok)
    call require(ok,'pinned default level population kick')
    if(rank==0) then
      call require(all(automatic(1)%position==manual(1)%position).and. &
        all(automatic(1)%velocity==manual(1)%velocity),'default CFL retains one native update')
      call require(all(feedback%fields(1)%values==delta%fields(1)%values),'default single source')
    endif
    call advance_mpi_amr_level_halo(MPI_COMM_WORLD,species,transport,liquid,settings,tree,[0],gas,temp, &
      1,0,1,100,nodes,automatic,mass,lower,spacing,[.true.,.true.,.true.],.5_dp*full_dt,full_dt,faces, &
      feedback,leftover,wall_auto,exterior_auto,ok,flow_cfl=merge(1.2_dp,.6_dp,rank==0))
    call require(.not.ok,'explicit CFL must agree on empty and populated ranks')
    if(rank==0) call require(all(automatic(1)%position==manual(1)%position),'CFL disagreement rollback')
    settings%drag=.false.
    do g=1,3
      n=merge(1,0,g==1.and.rank==0)
      allocate(population(g)%parcels(n),population(g)%component_mass(1,n),population(g)%nodes(n))
      if(n>0) then
        population(g)%parcels(1)=original;population(g)%parcels(1)%velocity=[1000._dp,0._dp,0._dp]
        population(g)%component_mass(1,1)=original%mass;population(g)%nodes=1
      endif
    enddo
    call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,settings,tree,[0], &
      gas,temp,1,0,1,100,population,stage_work,0._dp,full_dt,0,2,faces,feedback,leftover, &
      wall_auto,exterior_auto,ok,lower,spacing,flow_cfl=1.2_dp)
    call require(ok.and.stage_work%native_flow_cfl==1.2_dp,'SDC records explicit native CFL')
    call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,settings,tree,[0], &
      gas,temp,1,0,1,100,population,stage_work,0._dp,full_dt,1,2,faces,feedback,leftover, &
      wall_auto,exterior_auto,ok,lower,spacing,flow_cfl=.6_dp)
    call require(.not.ok.and.stage_work%next_stage==1,'cached stage rejects changed native CFL')
    call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,settings,tree,[0], &
      gas,temp,1,0,1,100,population,stage_work,0._dp,full_dt,1,2,faces,feedback,leftover, &
      wall_auto,exterior_auto,ok,lower,spacing,flow_cfl=1.2_dp)
    call require(ok.and.stage_work%next_stage==2,'cached stage keeps agreed native CFL')
    call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,settings,tree,[0], &
      gas,temp,1,0,1,100,population,stage_work,0._dp,full_dt,2,2,faces,feedback,leftover, &
      wall_auto,exterior_auto,ok,lower,spacing)
    call require(.not.ok.and.stage_work%next_stage==2,'final stage rejects dropped explicit CFL')
    if(rank==0) call require(all(population(1)%parcels(1)%position==original%position),'failed stage keeps physical state')
    call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,settings,tree,[0], &
      gas,temp,1,0,1,100,population,stage_work,0._dp,full_dt,2,2,faces,feedback,leftover, &
      wall_auto,exterior_auto,ok,lower,spacing,flow_cfl=1.2_dp)
    call require(ok.and.stage_work%next_stage==0,'final native CFL stage completes')
    if(rank==0) call require(abs(population(1)%parcels(1)%position(1)-original%position(1)-.0012_dp)<1e-15_dp, &
      'explicit CFL old drift and single final kick')
  end subroutine check_population_flow_cycles
  subroutine check_tab_history()
    type(spray_parcel), allocatable :: particles(:),children(:)
    real(dp), allocatable :: inventory(:,:),history(:)
    integer, allocatable :: node_ids(:)
    type(spray_parcel) :: reference,saved
    type(spray_options) :: controls
    type(spray_film) :: skin
    type(amr_tree_state_3d) :: feedback
    real(dp) :: cm(1),delta(nv),re,transfer(1),hv(1),increment,total_tangent,uniforms(1),compensation(4),old_mass
    real(dp) :: saved_history,step_dt
    integer(int64) :: marker
    integer :: np,step
    np=merge(1,0,rank==0);allocate(particles(np),inventory(1,np),history(np),node_ids(np))
    controls%drag=.false.;controls%heat=.false.;controls%evaporation=.false.
    controls%breakup_model='etab';controls%surface_tension=.072_dp;controls%liquid_viscosity=.001_dp
    controls%etab_parcel_factor=1;controls%maximum_multiplicity=3
    reference=original;reference%velocity=[1000._dp,0._dp,0._dp];reference%deformation=1.1_dp
    cm=reference%mass;history=0;node_ids=1;marker=37;total_tangent=0;step_dt=5e-9_dp
    if(rank==0) particles(1)=reference
    inventory=reference%mass
    do step=1,2
      if(rank==0) then
        call advance_mixture_frozen_gas(species,transport,primitive,temperature,reference,liquid,cm,controls, &
          step_dt,delta,ok,reynolds_number=re)
        call require(ok,'TAB history reference transfer')
        reference%position=original%position;old_mass=reference%mass
        call build_mixture_film(species,transport,primitive,temperature,reference,liquid,[1._dp],controls, &
          skin,transfer,hv,ok)
        call require(ok,'TAB history reference film')
        call advance_spray_breakup(reference,liquid(1),controls,primitive(1),primitive(2:4),skin%density, &
          skin%viscosity,2*step_dt,ok,increment,reynolds_number=re)
        call require(ok,'TAB history reference update')
        total_tangent=total_tangent+increment;cm=cm*(reference%mass/old_mass)
        if(step==2) then
          saved=particles(1);saved_history=history(1)
          call khrt_reference_uniforms(reference%id,37_int64,uniforms)
          call construct_reference_etab_children(reference,reference%velocity,1._dp,total_tangent,uniforms, &
            38_int64,children,compensation,ok,drift_dt=0._dp)
          call require(ok.and.size(children)==1,'TAB history reference replacement')
        endif
      endif
      if(step==2) then
        marker=huge(marker)
        call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,controls,tree,[0],gas,temp,1,node_ids, &
          particles,inventory,lower,spacing,[.true.,.true.,.true.],step_dt,feedback,ok,last_id=marker, &
          tab_interval=2*step_dt,tab_tangent_history=history,last_flow_cycle=.true.,khrt_interval=2*step_dt)
        call require(.not.ok.and.marker==huge(marker),'late TAB history birth ID failure')
        if(rank==0) call require(history(1)==saved_history.and.particles(1)%mass==saved%mass.and. &
          particles(1)%deformation==saved%deformation,'TAB history rollback with parcel state')
        marker=37
      endif
      call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,controls,tree,[0],gas,temp,1,node_ids, &
        particles,inventory,lower,spacing,[.true.,.true.,.true.],step_dt,feedback,ok,last_id=marker, &
        tab_interval=2*step_dt,tab_tangent_history=history,last_flow_cycle=step==2,khrt_interval=2*step_dt)
      call require(ok,'TAB history subinterval kick')
      if(step==1) then
        call require(marker==37,'TAB birth deferred until final flow cycle')
        if(rank==0) call require(history(1)==total_tangent.and.history(1)>0.and.particles(1)%id==1, &
          'TAB tangent retained across flow cycles')
      else
        call require(marker==38,'final flow cycle reserves replacement ID')
        if(rank==0) call require(particles(1)%id==38.and.history(1)==0.and. &
          maxval(abs(particles(1)%velocity-children(1)%velocity))<1e-12_dp, &
          'replacement uses cumulative tangent and starts its own history')
      endif
    enddo
  end subroutine check_tab_history
  subroutine check_primitive_halo_sampling()
    type(amr_tree_state_3d) :: u,t,feedback
    type(spray_parcel), allocatable :: particles(:)
    real(dp), allocatable :: inventory(:,:)
    integer, allocatable :: node_ids(:)
    type(spray_options) :: controls
    type(spray_parcel) :: reference
    real(dp) :: q(nv),state_cell(nv),tc,cs,cm(1),delta(nv)
    integer :: x,y,z,np
    allocate(u%fields(1),t%fields(1));np=0
    if(rank==0) np=1
    allocate(particles(np),inventory(1,np),node_ids(np))
    controls%evaporation=.false.;controls%drag=.false.
    reference=original;reference%position=.002_dp;reference%velocity=0;cm=reference%mass
    if(rank==0) then
      particles(1)=reference;inventory=reference%mass;node_ids=1
      allocate(u%fields(1)%values(nv,0:5,0:5,0:5),t%fields(1)%values(1,0:5,0:5,0:5))
      do z=0,5
      do y=0,5
      do x=0,5
        q=primitive;q(2)=merge(1000._dp,-1000._dp,mod(x,2)==0)
        call reactive_primitive_to_conserved(species,q,state_cell,tc,cs,ok)
        call require(ok,'opposing velocity halo cell')
        u%fields(1)%values(:,x,y,z)=state_cell;t%fields(1)%values(1,x,y,z)=tc
      enddo
      enddo
      enddo
      call advance_mixture_frozen_gas(species,transport,primitive,temperature,reference,liquid,cm,controls, &
        1e-8_dp,delta,ok)
      call require(ok,'primitive halo reference kick')
    endif
    call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,controls,tree,[0],u,t,1,node_ids, &
      particles,inventory,lower,spacing,[.true.,.true.,.true.],1e-8_dp,feedback,ok)
    call require(ok,'primitive halo kick with empty peer ranks')
    if(rank==0) then
      call require(abs(particles(1)%temperature-reference%temperature)<1e-10_dp, &
        'opposing gas velocities do not heat sampled gas through kinetic variance')
      call require(all(particles(1)%velocity==reference%velocity).and.all(inventory==cm(1)), &
        'primitive halo transfer velocity and liquid inventory')
    endif
  end subroutine check_primitive_halo_sampling
  subroutine check_breakup_clock(embedded)
    logical, intent(in) :: embedded
    type(spray_level_population) :: population(3)
    type(spray_sdc_population_stage) :: clock_work
    type(spray_options) :: breakup_options
    type(spray_parcel) :: reference
    type(spray_khrt_shedding_state) :: shedding
    type(spray_khrt_release) :: release
    type(spray_parcel), allocatable :: children(:)
    type(spray_film) :: film
    type(eb_geometry_3d) :: geometry(1)
    type(amr_tree_state_3d) :: u,t,centroids,metrics,rate,leftover
    real(dp) :: phi(0:4,0:4,0:4),transfer(1),enthalpy(1),breakup_dt,walls(5),outside(nv),local_tab
    real(dp) :: cm(1),delta(nv),transfer_re,position(3),old_mass,tangent,uniforms(1),compensation(4)
    integer(int64) :: marker
    integer :: model,group,x,y,z,s
    allocate(u%fields(1),t%fields(1),centroids%fields(1),metrics%fields(1))
    phi=1
    call build_level_set_eb_geometry_3d(phi,lower,lower+4*spacing,geometry(1),ok)
    call require(ok,'breakup clock geometry')
    if(rank==0) then
      allocate(u%fields(1)%values(nv,4,4,4),t%fields(1)%values(1,4,4,4), &
        centroids%fields(1)%values(3,0:5,0:5,0:5),metrics%fields(1)%values(1,0:5,0:5,0:5))
      do z=1,4
      do y=1,4
      do x=1,4
        u%fields(1)%values(:,x,y,z)=cell
      enddo
      enddo
      enddo
      t%fields(1)%values=temperature;metrics%fields(1)%values=1
      do z=0,5
      do y=0,5
      do x=0,5
        centroids%fields(1)%values(:,x,y,z)=lower+(real([x,y,z],dp)-.5_dp)*spacing
      enddo
      enddo
      enddo
    endif
    do group=1,3
      count=merge(1,0,rank==0.and.group==1)
      allocate(population(group)%parcels(count),population(group)%component_mass(1,count),population(group)%nodes(count))
      population(group)%nodes=1
    enddo
    breakup_options=options;breakup_options%drag=.false.;breakup_options%heat=.false.
    breakup_options%surface_tension=.072_dp;breakup_options%liquid_viscosity=.001_dp
    do model=1,6
      breakup_options%drag=model>2.and.model/=5
      breakup_options%etab_parcel_factor=-1;breakup_options%maximum_multiplicity=0;marker=37
      if(model>=5) then
        breakup_options%etab_parcel_factor=1;breakup_options%maximum_multiplicity=3
      endif
      breakup_options%breakup_model='etab';breakup_dt=interval
      if(model<=4.and.mod(model,2)==0) then
        breakup_options%breakup_model='khrt';breakup_dt=.5_dp*interval
      endif
      clock_work=spray_sdc_population_stage()
      reference=original;reference%velocity=[1000._dp,0._dp,0._dp];reference%deformation=1.1_dp
      cm=reference%mass
      if(rank==0) then
        population(1)%parcels(1)=reference;population(1)%component_mass(1,1)=reference%mass
      endif
      do s=0,2
        if(s/=1) then
          position=reference%position
          call advance_mixture_frozen_gas(species,transport,primitive,temperature,reference,liquid,cm, &
            breakup_options,.5_dp*interval,delta,ok,reynolds_number=transfer_re)
          call require(ok,'reference transfer for breakup Reynolds')
          reference%position=position;old_mass=reference%mass
          call build_mixture_film(species,transport,primitive,temperature,reference,liquid,[1._dp],breakup_options, &
            film,transfer,enthalpy,ok)
          call require(ok,'reference breakup clock film')
          if(model>=5) then
            call advance_spray_breakup(reference,liquid(1),breakup_options,primitive(1),primitive(2:4), &
              film%density,film%viscosity,breakup_dt,ok,tangent,reynolds_number=transfer_re)
          else
            call advance_spray_breakup(reference,liquid(1),breakup_options,primitive(1),primitive(2:4), &
              film%density,film%viscosity,breakup_dt,ok,reynolds_number=transfer_re)
          endif
          call require(ok,'reference breakup clock update')
          if(model>=5.and.s==2) then
            call advance_tab_no_birth_khrt(reference,liquid(1),breakup_options,primitive(1),primitive(2:4), &
              transfer_re,.5_dp*interval,ok)
            call require(ok,'reference no-birth TAB to KHRT fallback')
          endif
          cm=cm*(reference%mass/old_mass)
          if(model>=5.and.s==0) then
            call khrt_reference_uniforms(reference%id,37_int64,uniforms)
            call construct_reference_etab_children(reference,reference%velocity,1._dp,tangent,uniforms, &
              38_int64,children,compensation,ok,drift_dt=0._dp)
            call require(ok.and.size(children)==1,'reference old-stage ETAB replacement')
            reference=children(1)
          endif
          if(s==0) reference%position=modulo(reference%position+interval*reference%velocity,4*spacing)
        endif
        if(embedded) then
          call advance_mpi_eb_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,breakup_options,tree,[0], &
            geometry,u,t,1,0,1,100,population,clock_work,0._dp,interval,s,2,faces,rate,leftover,walls,outside,ok, &
            phi,centroids,metrics,last_id=marker)
        else
          call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,breakup_options,tree,[0], &
            gas,temp,1,0,1,100,population,clock_work,0._dp,interval,s,2,faces,rate,leftover,walls,outside,ok,lower,spacing, &
            last_id=marker)
        endif
        if(.not.ok.and.rank==0) print *, 'Breakup clock case: EB=',embedded,' model=',model,' stage=',s
        call require(ok,'regular or EB breakup clock stage')
        if(model>=5) call require(marker==38,'old replacement reserves one ID; final creates no new IDs')
        if(rank==0) then
          if(s<2) then
            call require(clock_work%midpoint(1)%parcels(1)%mass==reference%mass.and. &
              clock_work%midpoint(1)%parcels(1)%deformation==reference%deformation, &
              'old or cached breakup clock matches TAB full and KHRT half interval')
          else
            call require(population(1)%parcels(1)%mass==reference%mass.and. &
              population(1)%parcels(1)%deformation==reference%deformation, &
              'final breakup clock matches TAB full and KHRT half interval')
            if(model>=5) call require(size(population(1)%parcels)==1.and.population(1)%parcels(1)%id==38, &
              'final reference ETAB retains replacement parcel identity')
          endif
        endif
      enddo
    enddo
    if(.not.embedded) then
      if(rank==0) reference=population(1)%parcels(1)
      local_tab=interval
      call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,breakup_options,tree,[0],gas,temp,1, &
        population(1)%nodes,population(1)%parcels,population(1)%component_mass,lower,spacing, &
        [.true.,.true.,.true.],.5_dp*interval,rate,ok,last_id=marker,etab_births=rank==2)
      call require(.not.ok,'idle rank ETAB birth-policy disagreement rejects')
      if(rank==2) then
        call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,breakup_options,tree,[0],gas,temp,1, &
          population(1)%nodes,population(1)%parcels,population(1)%component_mass,lower,spacing, &
          [.true.,.true.,.true.],.5_dp*interval,rate,ok,last_id=marker)
      else
        call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,breakup_options,tree,[0],gas,temp,1, &
          population(1)%nodes,population(1)%parcels,population(1)%component_mass,lower,spacing, &
          [.true.,.true.,.true.],.5_dp*interval,rate,ok,last_id=marker,etab_births=.true.)
      endif
      call require(.not.ok,'idle rank ETAB birth-policy presence disagreement rejects')
      if(rank==2) local_tab=.5_dp*interval
      call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,breakup_options,tree,[0],gas,temp,1, &
        population(1)%nodes,population(1)%parcels,population(1)%component_mass,lower,spacing, &
        [.true.,.true.,.true.],.5_dp*interval,rate,ok,tab_interval=local_tab)
      call require(.not.ok,'idle rank TAB interval disagreement rejects')
      if(rank==2) then
        call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,breakup_options,tree,[0],gas,temp,1, &
          population(1)%nodes,population(1)%parcels,population(1)%component_mass,lower,spacing, &
          [.true.,.true.,.true.],.5_dp*interval,rate,ok)
      else
        call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,breakup_options,tree,[0],gas,temp,1, &
          population(1)%nodes,population(1)%parcels,population(1)%component_mass,lower,spacing, &
          [.true.,.true.,.true.],.5_dp*interval,rate,ok,tab_interval=interval)
      endif
      call require(.not.ok,'idle rank TAB interval presence disagreement rejects')
      if(rank==0) call require(population(1)%parcels(1)%mass==reference%mass.and. &
        population(1)%parcels(1)%deformation==reference%deformation,'TAB interval disagreement population rollback')
    endif
    breakup_options%breakup_model='khrt';breakup_options%etab_parcel_factor=-1
    breakup_options%maximum_multiplicity=0;breakup_options%drag=.false.
    breakup_options%khrt_children=.true.;breakup_options%khrt_parcel_factor=0
    breakup_options%khrt_reference_diameter_cubed=1._dp
    clock_work=spray_sdc_population_stage();marker=37
    reference=original;reference%mass=1e-9_dp;reference%multiplicity=.01_dp
    reference%velocity=[300._dp,0._dp,0._dp];reference%khrt_initial_count=reference%multiplicity
    cm=reference%mass
    if(rank==0) then
      population(1)%parcels=[reference];population(1)%component_mass(1,1)=reference%mass
    endif
    do s=0,2
      if(s==2) then
        ! Force the release threshold on MK: size and shedding evolve, births remain disabled.
        reference%khrt_shed_volume=.95_dp*reference%multiplicity*(.5_dp*parcel_diameter(reference,liquid(1)))**3
        if(rank==0) clock_work%midpoint(1)%parcels(1)%khrt_shed_volume=reference%khrt_shed_volume
      endif
      if(s/=1) then
        position=reference%position
        call advance_mixture_frozen_gas(species,transport,primitive,temperature,reference,liquid,cm, &
          breakup_options,.5_dp*interval,delta,ok,reynolds_number=transfer_re)
        call require(ok,'KHRT stage reference transfer')
        reference%position=position
        shedding=spray_khrt_shedding_state(reference%khrt_initial_count,reference%khrt_shed_volume)
        call advance_khrt_shedding(reference,liquid(1),breakup_options,primitive(1),primitive(2:4), &
          transfer_re,.5_dp*interval, &
          1._dp,shedding,release,ok,make_new_drops=s==0)
        call require(ok,'reference KHRT no-birth stage')
        reference%khrt_initial_count=shedding%initial_count;reference%khrt_shed_volume=shedding%shed_volume
        cm=reference%mass
        if(s==0) reference%position=modulo(reference%position+interval*reference%velocity,4*spacing)
      endif
      if(embedded) then
        call advance_mpi_eb_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,breakup_options,tree,[0], &
          geometry,u,t,1,0,1,100,population,clock_work,0._dp,interval,s,2,faces,rate,leftover,walls,outside,ok, &
          phi,centroids,metrics,last_id=marker)
      else
        call advance_mpi_amr_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,breakup_options,tree,[0], &
          gas,temp,1,0,1,100,population,clock_work,0._dp,interval,s,2,faces,rate,leftover,walls,outside,ok,lower,spacing, &
          last_id=marker)
      endif
      call require(ok.and.marker==37,'KHRT stage keeps watermark without release')
      if(rank==0.and.s==2) then
        call require(size(population(1)%parcels)==1.and.population(1)%parcels(1)%mass>0, &
          'KHRT MK forbids release above shed threshold')
        call require(population(1)%parcels(1)%mass==reference%mass.and. &
          population(1)%parcels(1)%khrt_shed_volume==reference%khrt_shed_volume, &
          'KHRT MK preserves reference coefficients and accumulated shedding')
      endif
    enddo
  end subroutine check_breakup_clock
  subroutine check_film_stage()
    type(spray_wall_film), allocatable :: films(:)
    real(dp), allocatable :: mass(:,:)
    type(spray_level_film_population), allocatable :: levels(:),virtual_films(:),ghost_films(:)
    type(spray_sdc_film_stage), allocatable :: hierarchy_work(:)
    type(spray_level_film_population) :: populations(3)
    type(spray_sdc_film_stage) :: film_work
    type(spray_wall_options) :: settings
    type(spray_wall_film) :: initial
    type(amr_tree_state_3d) :: rate,leftover,actual,old_rate
    real(dp) :: old(nv),new(nv),delta(nv),walls(5),saved_wall(5)
    integer :: p,c
    settings%enabled=.true.;settings%evolve=.true.;settings%temperature=400
    settings%liquid_conductivity=.6_dp;settings%minimum_height=1e-15_dp
    initial=spray_wall_film([0._dp,.0015_dp,.0015_dp],[1._dp,0._dp,0._dp],1e-10_dp,330._dp,1e-4_dp,1e-5_dp,1)
    c=merge(1,0,rank==0);allocate(films(c),mass(1,c))
    if(c==1) then
      films(1)=initial;mass(1,1)=initial%mass
    endif
    call initialize_mpi_amr_film_levels(MPI_COMM_WORLD,tree,[0],films,mass,lower,spacing,100,levels,ok)
    call require(ok,'regular film level initialization')
    populations(1)=levels(1)
    ! Deliberate copy verifies physical export and wall accounting exclude copies.
    populations(2)=levels(1)
    allocate(populations(3)%films(0),populations(3)%component_mass(1,0), &
      populations(3)%nodes(0),populations(3)%probes(3,0))
    do p=0,2
      call advance_mpi_amr_film_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,settings,tree,[0], &
        gas,temp,1,0,populations,film_work,lower,spacing,[.false.,.false.,.false.],0._dp,interval,p,2, &
        rate,leftover,walls,ok,actual)
      call require(ok,'regular film SDC stage')
      if(p==0) then
        old_rate=rate;call integrate(actual,old);saved_wall=walls
      endif
      if(p<2.and.rank==0) call require(populations(1)%films(1)%temperature==initial%temperature, &
        'film predictor private')
      if(p==1) then
        call require(all(walls==saved_wall),'film intermediate wall retained')
        if(rank==0) call require(all(rate%fields(1)%values==old_rate%fields(1)%values),'film intermediate source retained')
      endif
    enddo
    call integrate(actual,new);delta=0
    if(rank==0) then
      delta(:5)=wall_film_totals(populations(1)%films(1),liquid(1))-wall_film_totals(initial,liquid(1))
      call require(populations(1)%films(1)%temperature>initial%temperature,'film heating')
      call require(populations(1)%films(1)%mass==initial%mass,'film inactive evaporation mass')
    endif
    call MPI_Allreduce(MPI_IN_PLACE,delta,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    delta(:5)=delta(:5)-walls
    call require(maxval(abs(.5_dp*interval*(old+new)+delta))<1e-19_dp,'physical film gas and wall balance')
    call require(film_work%next_stage==0,'regular film lifecycle complete')
    allocate(virtual_films(1),ghost_films(1),hierarchy_work(1))
    virtual_films(1)=populations(3);ghost_films(1)=populations(3)
    do p=0,2
      call advance_mpi_amr_recursive_film_stage(MPI_COMM_WORLD,species,transport,liquid,options,settings,tree,[0], &
        gas,temp,1,0,1,100,levels,virtual_films,ghost_films,hierarchy_work,lower,spacing,0._dp,interval,p,2, &
        faces,rate,leftover,walls,ok,actual)
      call require(ok,'regular recursive film hierarchy stage')
    enddo
    if(rank==0) call require(levels(1)%films(1)%temperature==populations(1)%films(1)%temperature, &
      'hierarchy film matches physical stage')
    call require(hierarchy_work(1)%next_stage==0,'recursive film lifecycle complete')
    if(rank==0) print *, 'PASS: regular film SDC privacy, cached stage, physical export, wall balance and hierarchy'
  end subroutine check_film_stage
  subroutine check_coupled(refined,recursive,callbacks,films_enabled,impact_enabled,mol_enabled,rotating)
    logical, intent(in) :: refined,recursive
    logical, optional, intent(in) :: callbacks,films_enabled,impact_enabled,mol_enabled,rotating
    logical :: custom,with_films,impact,mol
    type(rotating_frame_3d), allocatable :: frame,bad_frame
    real(dp) :: potential,frame_change,saved_frame_change
    type(spray_parcel) :: starting
    type(spray_options) :: local_options
    type(spray_wall_film), allocatable :: films(:),saved_films(:)
    real(dp), allocatable :: film_mass(:,:),saved_film_mass(:,:)
    type(spray_wall_film) :: initial_film
    type(spray_wall_options) :: film_settings
    character(len=16) :: coupled_faces(6)
    integer :: film_steps
    type(amr_patch_tree_3d) :: mesh
    type(amr_patch_3d) :: child(1)
    integer, allocatable :: owners(:)
    real(dp) :: particle_force(5),accel(3),volume,saved_force(5),bad_accel(3)
    integer :: node,ancestor,d(3)
    type(amr_tree_state_3d) :: u,t,saved_u,saved_t,heat,saved_heat,initial_u
    type(spray_parcel), allocatable :: drops(:),saved_drops(:)
    type(elementary_reaction) :: reactions(0)
    real(dp), allocatable :: masses(:,:),saved_masses(:,:)
    real(dp) :: budget(nv),external(nv),exchange(nv),outside(nv),particle_wall(5),largest,theta,inventory(nv)
    real(dp) :: saved_budget(nv),saved_external(nv),saved_exchange(nv),saved_outside(nv),saved_wall(5)
    real(dp) :: saved_largest,saved_theta
    integer(kind=8) :: watermark,saved_watermark
    integer :: c,x,y,z
    custom=.false.;with_films=.false.;impact=.false.;coupled_faces=faces
    starting=original;local_options=options
    mol=.false.
    if(present(mol_enabled)) mol=mol_enabled
    if(present(callbacks)) custom=callbacks
    if(present(films_enabled)) with_films=films_enabled
    if(present(impact_enabled)) impact=impact_enabled
    if(impact) then
      starting%position(1)=.25_dp*interval;starting%velocity=[-1._dp,0._dp,0._dp]
      local_options%surface_tension=.072_dp;local_options%liquid_viscosity=.001_dp
    endif
    if(with_films) then
      coupled_faces(1:2)='slip_wall'
      film_settings%enabled(1:2)=.true.;film_settings%evolve=.true.;film_settings%temperature=400
      film_settings%liquid_conductivity=.6_dp;film_settings%minimum_height=1e-15_dp
      initial_film=spray_wall_film([0._dp,.0015_dp,.0015_dp],[1._dp,0._dp,0._dp],1e-10_dp,330._dp,1e-4_dp,1e-5_dp,1)
      if(impact) film_settings%temperature=300
      c=merge(1,0,rank==0.and..not.impact);allocate(films(c),film_mass(1,c))
      if(c==1) then
        films(1)=initial_film;film_mass(1,1)=initial_film%mass
      endif
    endif
    mesh=tree;owners=[0]
    if(refined) then
      call initialize_amr_patch_3d(4,4,4,1,1,1,2,2,2,2,child(1),ok)
      call require(ok,'refined patch')
      call initialize_amr_patch_tree_3d([4,4,4],[1],child,mesh,ok)
      call require(ok,'refined tree')
      owners=[0,1]
    endif
    allocate(u%fields(size(owners)),t%fields(size(owners)),heat%fields(size(owners)));count=0
    if(rank==0) count=1
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      d=mesh%nodes(node)%shape
      allocate(u%fields(node)%values(nv,d(1),d(2),d(3)),t%fields(node)%values(1,d(1),d(2),d(3)))
      do z=1,d(3)
      do y=1,d(2)
      do x=1,d(1)
        u%fields(node)%values(:,x,y,z)=cell
      enddo
      enddo
      enddo
      t%fields(node)%values=temperature
      allocate(heat%fields(node)%values(1,d(1),d(2),d(3)),source=53._dp)
    enddo
    if(present(rotating)) then
      if(rotating) then
        allocate(frame)
        ! Rotation about an off-center axis is not periodic across radial faces.
        ! Use impermeable radial walls to test a compatible rotating domain.
        coupled_faces='slip_wall'
        frame%axis=3;frame%angular_speed=1000._dp;frame%axis_point=[.001_dp,.001_dp,.001_dp]
        call convert_mpi_amr_tree_rotating_energy(MPI_COMM_WORLD,mesh,owners,u,lower,spacing,frame,.true.,ok)
        call require(ok,'coupled rotating initial gas energy')
      endif
    endif
    initial_u=u
    allocate(drops(count),masses(1,count))
    if(rank==0) then
      drops(1)=starting;masses(1,1)=starting%mass
    endif
    frame_change=0
    particle_force=0;accel=[1000._dp,-2000._dp,3000._dp]
    if(impact) accel=0
    watermark=1;budget=0;external=0;exchange=0;particle_wall=0;outside=0;largest=0;theta=1
    callback_mesh=mesh;diffusion_calls=0;reaction_calls=0;callback_failure=.false.
    if(custom) then
    call advance_mpi_amr_spray_sdc(MPI_COMM_WORLD,mesh,owners,species,reactions,transport,liquid,local_options,u,t, &
      drops,masses,lower,spacing,interval,2,.true.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp], &
      1e-6_dp,1e-14_dp,coupled_faces,budget,external,exchange,particle_wall,outside,largest,theta,watermark,ok, &
      particle_acceleration=accel,particle_force_change=particle_force,subcycling=recursive,mol_mode=mol, &
      frame=frame,frame_energy_change=frame_change, &
      evaluate_diffusion=custom_diffusion,evaluate_reaction=custom_reaction,heat_release=heat)
    call require(diffusion_calls>0,'custom diffusion reached from coupled SDC')
    c=reaction_calls
    call MPI_Allreduce(MPI_IN_PLACE,c,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
    call require(c>0,'custom forced chemistry reached from coupled SDC')
    else if(with_films) then
    call advance_mpi_amr_spray_sdc(MPI_COMM_WORLD,mesh,owners,species,reactions,transport,liquid,local_options,u,t, &
      drops,masses,lower,spacing,interval,2,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp], &
      1e-6_dp,1e-14_dp,coupled_faces,budget,external,exchange,particle_wall,outside,largest,theta,watermark,ok, &
      particle_acceleration=accel,particle_force_change=particle_force,subcycling=recursive,mol_mode=mol, &
      frame=frame,frame_energy_change=frame_change, &
      wall_films=films,film_mass=film_mass,film_settings=film_settings,film_steps=film_steps,film_impacts=impact,heat_release=heat)
    else
    call advance_mpi_amr_spray_sdc(MPI_COMM_WORLD,mesh,owners,species,reactions,transport,liquid,local_options,u,t, &
      drops,masses,lower,spacing,interval,2,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp], &
      1e-6_dp,1e-14_dp,coupled_faces,budget,external,exchange,particle_wall,outside,largest,theta,watermark,ok, &
      particle_acceleration=accel,particle_force_change=particle_force,subcycling=recursive,mol_mode=mol, &
      frame=frame,frame_energy_change=frame_change, &
      heat_release=heat)
    endif
    call require(ok,'coupled gas particle SDC')
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      call require(all(ieee_is_finite(heat%fields(node)%values)),'finite coupled stage heat')
      if(.not.custom) call require(all(heat%fields(node)%values==0),'nonreacting coupled heat zero')
    end do
    saved_heat=heat
    expected=0;inventory=0
    if(rank==0) then
      expected(:5)=-parcel_totals(starting,liquid(1))
      if(allocated(frame)) then
        call rotating_potential_3d(frame,starting%position,potential,ok)
        call require(ok,'initial particle potential')
        expected(5)=expected(5)+starting%mass*starting%multiplicity*potential
      endif
    endif
    do c=1,size(drops)
      if(.not.impact) call require(size(drops)==1.and.drops(c)%temperature>starting%temperature, &
        'coupled physical publication')
      expected(:5)=expected(:5)+parcel_totals(drops(c),liquid(1))
      if(allocated(frame)) then
        call rotating_potential_3d(frame,drops(c)%position,potential,ok)
        call require(ok,'final particle potential')
        expected(5)=expected(5)-drops(c)%mass*drops(c)%multiplicity*potential
      endif
    enddo
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      d=mesh%nodes(node)%shape;volume=product(spacing);ancestor=node
      do while(ancestor>1)
        volume=volume/real(mesh%nodes(ancestor)%patch%refinement_ratio,dp)**3
        ancestor=mesh%nodes(ancestor)%parent
      enddo
      do z=1,d(3)
      do y=1,d(2)
      do x=1,d(1)
        if(covered(mesh,node,[x,y,z])) cycle
        inventory=inventory+volume*(u%fields(node)%values(:,x,y,z)-initial_u%fields(node)%values(:,x,y,z))
      enddo
      enddo
      enddo
    enddo
    if(with_films) then
      do c=1,size(films)
        expected(:5)=expected(:5)+wall_film_totals(films(c),liquid(1))
        if(allocated(frame)) then
          call rotating_potential_3d(frame,films(c)%position,potential,ok)
          call require(ok,'final film potential')
          expected(5)=expected(5)-films(c)%mass*potential
        endif
        if(.not.impact) call require(films(c)%temperature>initial_film%temperature,'coupled film heating')
      enddo
      if(rank==0.and..not.impact) then
        expected(:5)=expected(:5)-wall_film_totals(initial_film,liquid(1))
        if(allocated(frame)) then
          call rotating_potential_3d(frame,initial_film%position,potential,ok)
          call require(ok,'initial film potential')
          expected(5)=expected(5)+initial_film%mass*potential
        endif
      endif
      if(impact) then
        c=size(films)
        call MPI_Allreduce(MPI_IN_PLACE,c,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
        call require(c>0,'physical impact creates accepted film')
        c=size(drops)
        call MPI_Allreduce(MPI_IN_PLACE,c,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
        call require(c==0,'deposited particle retired')
      endif
      call require(film_steps==merge(6,merge(4,2,refined),recursive),'coupled film accepted half kicks')
    endif
    call MPI_Allreduce(MPI_IN_PLACE,expected,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call MPI_Allreduce(MPI_IN_PLACE,inventory,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    expected(:5)=expected(:5)-particle_force
    expected(5)=expected(5)-frame_change
    if(.not.allocated(frame)) call require(frame_change==0, &
      "ordinary particle energy has zero coordinate ledger")
    if(with_films) expected(:5)=expected(:5)-particle_wall
    if(rank==0.and.maxval(abs(exchange+expected))>=1e-19_dp) then
      print *, 'exchange:',exchange
      print *, 'liquid minus ledgers:',expected
      print *, 'coordinate change:',frame_change
    endif
    call MPI_Barrier(MPI_COMM_WORLD,status)
    call require(maxval(abs(exchange+expected))<1e-19_dp,'coupled particle exchange and force ledgers')
    call require(maxval(abs(particle_force(2:4)-starting%mass*starting%multiplicity*interval*accel))<1e-24_dp, &
      'physical particle force excludes copy counts')
    call require(maxval(abs(inventory-budget-external-exchange))<1e-16_dp,'coupled gas exchange balance')
    call require(all(outside==0),'coupled particle outflow inventory')
    if(.not.with_films) call require(all(particle_wall==0),'periodic coupled wall inventory')
    if(custom) then
      saved_u=u;saved_t=t;saved_force=particle_force;saved_frame_change=frame_change;callback_failure=.true.;diffusion_calls=0
      saved_drops=drops;saved_masses=masses;saved_watermark=watermark
      saved_budget=budget;saved_external=external;saved_exchange=exchange;saved_outside=outside
      saved_wall=particle_wall;saved_largest=largest;saved_theta=theta
      call advance_mpi_amr_spray_sdc(MPI_COMM_WORLD,mesh,owners,species,reactions,transport,liquid,local_options,u,t, &
        drops,masses,lower,spacing,interval,2,.true.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp], &
        1e-6_dp,1e-14_dp,coupled_faces,budget,external,exchange,particle_wall,outside,largest,theta,watermark,ok, &
        particle_acceleration=accel,particle_force_change=particle_force,subcycling=recursive,mol_mode=mol, &
      frame=frame,frame_energy_change=frame_change, &
        evaluate_diffusion=custom_diffusion,evaluate_reaction=custom_reaction,heat_release=heat)
      call require(.not.ok.and.all(particle_force==saved_force).and.frame_change==saved_frame_change, &
        'late custom diffusion failure rejects')
      call require(all(masses==saved_masses).and.size(drops)==size(saved_drops),'diffusion failure particle inventory rollback')
      do c=1,size(drops)
        call require(all(drops(c)%position==saved_drops(c)%position).and. &
          all(drops(c)%velocity==saved_drops(c)%velocity).and.drops(c)%temperature==saved_drops(c)%temperature, &
          'diffusion failure particle state rollback')
      enddo
      call require(all(budget==saved_budget).and.all(external==saved_external).and.all(exchange==saved_exchange).and. &
        all(outside==saved_outside).and.all(particle_wall==saved_wall),'diffusion failure all ledgers rollback')
      call require(watermark==saved_watermark.and.largest==saved_largest.and.theta==saved_theta, &
        'diffusion failure controls rollback')
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        call require(all(u%fields(node)%values==saved_u%fields(node)%values),'diffusion failure gas rollback')
        call require(all(t%fields(node)%values==saved_t%fields(node)%values),'diffusion failure temperature rollback')
        call require(all(heat%fields(node)%values==saved_heat%fields(node)%values),'late failure heat rollback')
      enddo
      callback_failure=.false.
      callback_negative=.true.;diffusion_calls=0
      call advance_mpi_amr_spray_sdc(MPI_COMM_WORLD,mesh,owners,species,reactions,transport,liquid,local_options,u,t, &
        drops,masses,lower,spacing,interval,2,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
        1e-6_dp,1e-14_dp,coupled_faces,budget,external,exchange,particle_wall,outside,largest,theta,watermark,ok, &
        particle_acceleration=accel,particle_force_change=particle_force,subcycling=recursive,mol_mode=mol, &
      frame=frame,frame_energy_change=frame_change, &
        evaluate_diffusion=custom_diffusion,heat_release=heat)
      call require(.not.ok.and.diffusion_calls>=merge(2,3,mol),'signed roundoff species cannot be published')
      call require(all(masses==saved_masses).and.size(drops)==size(saved_drops),'negative endpoint particle inventory rollback')
      do c=1,size(drops)
        call require(all(drops(c)%position==saved_drops(c)%position).and. &
          all(drops(c)%velocity==saved_drops(c)%velocity).and.drops(c)%temperature==saved_drops(c)%temperature, &
          'negative endpoint particle state rollback')
      enddo
      call require(all(budget==saved_budget).and.all(external==saved_external).and.all(exchange==saved_exchange).and. &
        all(outside==saved_outside).and.all(particle_wall==saved_wall).and.all(particle_force==saved_force).and. &
        frame_change==saved_frame_change, &
        'negative endpoint all ledgers rollback')
      call require(watermark==saved_watermark.and.largest==saved_largest.and.theta==saved_theta, &
        'negative endpoint controls rollback')
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        call require(all(u%fields(node)%values==saved_u%fields(node)%values),'negative endpoint gas rollback')
        call require(all(t%fields(node)%values==saved_t%fields(node)%values),'negative endpoint temperature rollback')
      enddo
      callback_negative=.false.
    endif
    saved_u=u;saved_t=t;saved_force=particle_force;saved_frame_change=frame_change;bad_accel=accel
    if(rank==2) bad_accel(1)=bad_accel(1)+1
    if(with_films) then
      saved_films=films;saved_film_mass=film_mass
    call advance_mpi_amr_spray_sdc(MPI_COMM_WORLD,mesh,owners,species,reactions,transport,liquid,local_options,u,t, &
      drops,masses,lower,spacing,interval,2,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp], &
      1e-6_dp,1e-14_dp,coupled_faces,budget,external,exchange,particle_wall,outside,largest,theta,watermark,ok, &
      particle_acceleration=bad_accel,particle_force_change=particle_force,subcycling=recursive,mol_mode=mol, &
      frame=frame,frame_energy_change=frame_change, &
      wall_films=films,film_mass=film_mass,film_settings=film_settings,film_impacts=impact)
      call require(all(film_mass==saved_film_mass),'failed coupled film inventory rollback')
      do c=1,size(films)
        call require(films(c)%temperature==saved_films(c)%temperature,'failed coupled film temperature rollback')
      enddo
    else
    call advance_mpi_amr_spray_sdc(MPI_COMM_WORLD,mesh,owners,species,reactions,transport,liquid,local_options,u,t, &
      drops,masses,lower,spacing,interval,2,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp], &
      1e-6_dp,1e-14_dp,coupled_faces,budget,external,exchange,particle_wall,outside,largest,theta,watermark,ok, &
      particle_acceleration=bad_accel,particle_force_change=particle_force,subcycling=recursive,mol_mode=mol, &
      frame=frame,frame_energy_change=frame_change, &
      heat_release=heat)
    endif
    call require(.not.ok.and.all(particle_force==saved_force).and.frame_change==saved_frame_change, &
      'mismatched particle force rejects without ledger publication')
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      call require(all(u%fields(node)%values==saved_u%fields(node)%values),'force agreement gas rollback')
      call require(all(t%fields(node)%values==saved_t%fields(node)%values),'force agreement temperature rollback')
    enddo
    if(allocated(frame)) then
      allocate(bad_frame);bad_frame=frame
      if(rank==2) bad_frame%angular_speed=bad_frame%angular_speed+1
      call advance_mpi_amr_spray_sdc(MPI_COMM_WORLD,mesh,owners,species,reactions,transport,liquid,local_options,u,t, &
        drops,masses,lower,spacing,interval,2,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
        1e-6_dp,1e-14_dp,coupled_faces,budget,external,exchange,particle_wall,outside,largest,theta,watermark,ok, &
        subcycling=recursive,mol_mode=mol,frame=bad_frame,frame_energy_change=frame_change)
      call require(.not.ok.and.frame_change==saved_frame_change,'mismatched coupled rotating frame rejects collectively')
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        call require(all(u%fields(node)%values==saved_u%fields(node)%values),'frame agreement gas rollback')
        call require(all(t%fields(node)%values==saved_t%fields(node)%values),'frame agreement temperature rollback')
      enddo
    endif
    if(rank==0) print *, 'PASS: coupled regular SDC/MOL exchange/force, refined=',refined,' recursive=',recursive,' MOL=',mol
  end subroutine
  subroutine custom_diffusion(state,temps,depth,fractions,coarse_start,coarse_end,rates,flux,good, &
      internal_energy)
      type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(amr_tree_state_3d), intent(in) :: state,temps
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
    integer, intent(in) :: depth
    real(dp), optional, intent(in) :: fractions(:)
    type(amr_tree_state_3d), intent(inout) :: rates
    type(amr_tree_flux_3d), allocatable, intent(inout) :: flux(:)
    logical, intent(out) :: good
    integer :: peer,cells(3)
    rates=state
    if(allocated(flux)) deallocate(flux)
    allocate(flux(size(callback_mesh%nodes)))
    do peer=1,size(callback_mesh%nodes)
      if(.not.allocated(state%fields(peer)%values)) cycle
      rates%fields(peer)%values=0
      if(depth>=0.and.callback_mesh%nodes(peer)%level/=depth) cycle
      if(callback_negative) then
        rates%fields(peer)%values(6,:,:,:)=-1e-20_dp
        rates%fields(peer)%values(nv,:,:,:)=1e-20_dp
      endif
      cells=callback_mesh%nodes(peer)%shape
      allocate(flux(peer)%x(nv,cells(1)+1,cells(2),cells(3)), &
        flux(peer)%y(nv,cells(1),cells(2)+1,cells(3)),flux(peer)%z(nv,cells(1),cells(2),cells(3)+1))
      flux(peer)%x=0;flux(peer)%y=0;flux(peer)%z=0
    enddo
    diffusion_calls=diffusion_calls+1
    good=.not.(callback_failure.and.rank==2.and.diffusion_calls==2)
  end subroutine
  subroutine custom_reaction(records,reactions,energy,forcing,energy_forcing,duration,rtol,atol, &
      densities,temp,steps,good,maximum_steps)
    type(nasa7_species), intent(in) :: records(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: energy,forcing(:),energy_forcing,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),temp
    integer, intent(out) :: steps
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    reaction_calls=reaction_calls+1
    call advance_forced_reactor_interval(records,reactions,energy,forcing,energy_forcing,duration,rtol,atol, &
      densities,temp,steps,good,maximum_steps=maximum_steps)
  end subroutine
  subroutine integrate(field,integral)
    type(amr_tree_state_3d), intent(in) :: field
    real(dp), intent(out) :: integral(:)
    integer :: c
    integral=0
    if(rank==0) then
      do c=1,nv
        integral(c)=sum(field%fields(1)%values(c,:,:,:))*product(spacing)
      enddo
    endif
    call MPI_Allreduce(MPI_IN_PLACE,integral,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(condition) return
    print *, 'FAIL rank',rank,': ',label
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
end program
