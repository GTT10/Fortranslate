program test_mpi_eb_native_mol_replay_3d
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: load_h2o2_full_transport
  use gas_transport_mod, only: gas_transport_species
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,eb_mol_held_interval_3d, &
    initialize_reactive_eb_tree_3d,build_eb_tree_extended_geometry_3d,advance_reactive_eb_tree_sdc_subcycle_3d, &
    integrate_reactive_eb_tree_3d
  use mpi_eb_tree_full_3d_mod, only: mpi_eb_mol_interval_3d,advance_mpi_eb_tree_sdc_subcycle_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_held_funding_queue_3d_mod, only: eb_held_funding_context_3d,valid_eb_held_funding_context_3d
  use dual_energy_mod, only: dual_energy_policy
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: patches(2)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,padded
  type(eb_geometry_3d), allocatable :: geometry(:),extended(:)
  type(eb_reactive_tree_3d) :: initial,serial,trial
  type(amr_tree_state_3d) :: initial_h,h,ir,heat,th,tir,theat,gas,temp,mh,mir,mheat,ug,ut,uh,uir,uheat
  type(eb_mol_held_interval_3d), allocatable :: records(:),replayed(:),bad(:)
  type(mpi_eb_mol_interval_3d), allocatable :: mrecords(:),mreplayed(:),mbad(:)
  type(eb_integrated_flux_3d), allocatable :: gf(:),hf(:),tgf(:),thf(:),mgf(:),mhf(:),ugf(:),uhf(:)
  type(eb_held_funding_context_3d) :: context,tcontext,mcontext,ucontext
  type(dual_energy_policy) :: policy
  real(dp), allocatable :: root(:,:,:,:),root_t(:,:,:),phi(:,:,:),padded_phi(:,:,:),q(:),cell(:)
  real(dp), allocatable :: b(:),w(:),ch(:),ext(:),aux(:),tb(:),tw(:),tch(:),te(:),ta(:)
  real(dp) :: t,sound,rg,correction,tcorrection,mc,uc,pressure,tp,mp,up,theta,ttheta,mt,umt
  real(dp), parameter :: dt=2e-8_dp,force(3)=[1e6_dp,-2e6_dp,3e6_dp],prefix=3.5_dp
  integer(int64) :: steps(3),tsteps(3),msteps(3),usteps(3)
  integer :: rank,nranks,status,nv,p,i,j,k,n(3),iterations,layout,owners(3),source_calls,saved_calls,case_id
  logical :: good,equal,forbid_sources=.false.,refused,refuse_chemistry=.false.
  logical :: selected_momentum=.false.,selected_pressure=.false.
  logical :: weighted,subparent,refuse_after_first_root=.false.
  real(dp) :: audit_interval=dt
  character(len=16) :: mode
  integer :: root_reactor_calls=0,fine_refusals=0,global_root_calls,global_refusals
  character(len=16) :: faces(6)='outflow'
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call get_command_argument(1,mode);weighted=mode=='weighted';subparent=mode=='subparent'
  selected_pressure=mode=='pressure'
  selected_momentum=mode=='momentum'.or.selected_pressure
  call require(nranks==3,'two owners and an idle rank')
  call load_h2o2_full_thermo(species,good);call require(good,'thermo')
  call load_h2o2_full_mechanism(reactions,good);call require(good,'chemistry')
  call load_h2o2_full_transport(transport,good);call require(good,'transport')
  nv=reactive_nvar(size(species));allocate(q(nv),cell(nv));q=0;q(1)=1;q(2:4)=[1._dp,2._dp,3._dp]
  do p=1,size(species)
    if(species(p)%name=='N2') q(5+p)=1
  enddo
  rg=mixture_specific_gas_constant(species,q(6:),good);call require(good,'gas constant')
  q(5)=rg*1200
  call reactive_primitive_to_conserved(species,q,cell,t,sound,good);call require(good,'initial cell')
  allocate(root(nv,4,4,4),root_t(4,4,4),phi(0:16,0:16,0:16),padded_phi(0:64,0:64,0:64))
  do p=1,nv
    root(p,:,:,:)=cell(p)
  enddo
  root_t=t
  if(weighted.or.subparent) then
    do k=1,4
    do j=1,4
    do i=1,4
      q(1)=1+.02_dp*i;q(2:4)=[1+.1_dp*i,2+.05_dp*j,3+.1_dp*k]
      q(5)=q(1)*rg*(1200+5*i+3*j)
      call reactive_primitive_to_conserved(species,q,cell,t,sound,good);call require(good,'nonuniform initial cell')
      root(:,i,j,k)=cell;root_t(i,j,k)=t
    enddo
    enddo
    enddo
  endif
  do j=0,16
    phi(:,j,:)=real(j,dp)*.00025_dp-.00131_dp
  enddo
  do j=0,64
    padded_phi(:,j,:)=real(max(0,min(16,j-24)),dp)*.00025_dp-.00131_dp
  enddo
  call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[.004_dp,.004_dp,.004_dp],finest,good)
  call require(good,'cut geometry')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),good);call require(good,'child')
  call initialize_amr_patch_3d(4,8,8,1,2,1,4,1,4,2,patches(2),good);call require(good,'grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,topology,good);call require(good,'three-level topology')
  call initialize_reactive_eb_tree_3d(species,topology,finest,root,root_t,initial,good,internal_energy=initial_h)
  call require(good,'initial paired family')
  allocate(geometry(3))
  do p=1,3
    geometry(p)=initial%nodes(p)%geometry
  enddo
  call build_level_set_eb_geometry_3d(padded_phi,[-.006_dp,-.006_dp,-.006_dp],[.010_dp,.010_dp,.010_dp],padded,good)
  call require(good,'extended geometry')
  call build_eb_tree_extended_geometry_3d(topology,geometry,padded,extended,good,halo_width=6)
  call require(good,'native donor geometry')
  allocate(b(nv),w(nv),ch(nv),ext(nv),aux(nv+1),tb(nv),tw(nv),tch(nv),te(nv),ta(nv+1))
  policy%eta2=2
  do iterations=1,2
    serial=initial;h=initial_h;call zero_field(ir,nv-4);call zero_field(heat,1)
    context=eb_held_funding_context_3d();correction=prefix;pressure=0;source_calls=0;forbid_sources=.false.
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,iterations,.true., &
      [.false.,.false.,.false.,.false.],force,1e-8_dp,1e-12_dp,serial,extended,faces,b,w,ch,ext,theta,good,ir, &
      mol_mode=.true.,redistribution='none',internal_energy=h,internal_energy_flux=hf,energy_policy=policy, &
      energy_correction=correction,heat_release=heat,held_interval_history=records,integrated_flux=gf, &
      auxiliary_source=serial_source,auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
      internal_energy_pressure_work_change=pressure,level_steps=steps,source_funding_context=context, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
    call require(good,'native serial root accepts')
    call require(all(records%separate_momentum.eqv.selected_momentum),'serial window binds policy')
    call require(size(records)==7.and.all(steps==[1_int64,2_int64,4_int64]),'seven physical clocks')
    call check_ancestors()
    call check_inventory()
    if(subparent) then
      call check_subparent_windows()
      cycle
    endif
    if(weighted) then
      call check_weighted_window()
      cycle
    endif
    saved_calls=source_calls;forbid_sources=.true.
    ! Fabricated completed tails and prefix must be replaced by the full anchor.
    trial=serial;th=h;tir=ir;theat=heat;tcorrection=999;tp=999;tcontext=context
    trial%nodes(1)%state(2,:,:,:)=trial%nodes(1)%state(2,:,:,:)+.01_dp
    call serial_replay(records,good)
    call require(good.and.source_calls==saved_calls,'serial full native replay freezes all external sources')
    call require(serial_equal(),'serial root replay replaces tails with exact gas/T/H/IR/heat/faces/budgets')
    do layout=1,2
      owners=[0,2,0]
      if(layout==2) owners=0
      call localize(initial,initial_h,gas,temp,mh);call localize_field(records(1)%clock_reaction_source,mir)
      call localize_field(records(1)%clock_heat,mheat)
      mcontext=eb_held_funding_context_3d();mc=prefix;mp=0;source_calls=0;forbid_sources=.false.
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
        dt,iterations,.true.,[.false.,.false.,.false.,.false.],force,1e-8_dp,1e-12_dp,extended,faces, &
        tb,tw,tch,te,mt,MPI_COMM_WORLD,good,mir,mol_mode=.true.,redistribution='none', &
        internal_energy=mh,internal_energy_flux=mhf,energy_policy=policy,energy_correction=mc,heat_release=mheat, &
        interval_history=mrecords,integrated_flux=mgf,auxiliary_source=native_mpi_source,auxiliary_change=ta, &
        auxiliary_has_internal_energy=.true.,internal_energy_pressure_work_change=mp,level_steps=msteps, &
        source_funding_context=mcontext, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
      call require(good,'native MPI root accepts')
      call require(all(mrecords%separate_momentum.eqv.selected_momentum),'MPI window binds policy')
      call check_mpi_ancestors(records,mrecords)
      call require(owner_equal(),'native MPI and serial accepted fields/faces/budgets exact')
      saved_calls=source_calls;forbid_sources=.true.
      ug=gas;ut=temp;uh=mh;uir=mir;uheat=mheat;uc=999;up=999;ucontext=mcontext
      call mpi_replay(mrecords,good)
      call require(good.and.source_calls==saved_calls,'MPI root replay freezes callbacks including idle ranks')
      call require(mpi_equal(),'MPI native replay restores all fields/faces/budgets exactly')
      ! Whole-window validation must reject a late corrupt packet before any
      ! chemistry, including mutations on a remote owner or an idle rank.
      do case_id=1,5
        mbad=mrecords
        select case(case_id)
        case(1)
          if(rank==owners(3)) mbad(7)%ancestor_internal_energy_end%fields(3)%values(1,1,1,1)= &
            ieee_value(0._dp,ieee_quiet_nan)
        case(2)
          if(rank==1) mbad(7)%fractions_end(1)=.75_dp
        case(3)
          if(rank==owners(1)) deallocate(mbad(7)%ancestor_start%fields(1)%values)
        case(4)
          if(rank==1) mbad(7)%time=mbad(7)%time+dt
        case(5)
          if(rank==1) mbad(7)%separate_momentum=.not.selected_momentum
          if(selected_pressure.and.rank==1) then
            mbad(7)%separate_momentum=selected_momentum;mbad(7)%paired_pressure=.false.
          endif
        end select
        call mpi_replay(mbad,good)
        call require(.not.good.and.source_calls==saved_calls,'late owner/idle archive corruption rejects before sources')
        call require(mpi_equal(),'MPI archive refusal preserves every output field, face and ledger')
      enddo
      refuse_chemistry=.true.;root_reactor_calls=0;fine_refusals=0
      call mpi_replay(mrecords,good)
      call require(.not.good.and.refused,'MPI native-admitted child provider refusal aborts root replay')
      call MPI_Allreduce(root_reactor_calls,global_root_calls,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
      call MPI_Allreduce(fine_refusals,global_refusals,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
      call require(global_root_calls>0.and.global_refusals==1,'completed root chemistry then one refusing child, no retry')
      call require(mpi_equal(),'MPI late chemistry refusal preserves complete accepted caller state')
      refuse_chemistry=.false.
    enddo
    refuse_chemistry=.true.;root_reactor_calls=0;fine_refusals=0
    call serial_replay(records,good)
    call require(.not.good.and.refused.and.root_reactor_calls>0.and.fine_refusals==1, &
      'serial completed root then native-admitted child refusal without retry')
    call require(serial_equal(),'serial late chemistry refusal preserves complete prior result')
    refuse_chemistry=.false.
    do case_id=1,4
      bad=records
      select case(case_id)
      case(1)
        bad(7)%ancestor_internal_energy_end%fields(3)%values(1,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
      case(2)
        bad(7)%fractions_end(1)=.75_dp
      case(3)
        bad(7)%time=bad(7)%time+dt
      case(4)
        bad(7)%separate_momentum=.not.selected_momentum
      end select
      call serial_replay(bad,good)
      call require(.not.good,'serial late archive refusal')
      call require(serial_equal(),'serial refusal preserves all prior accepted outputs')
    enddo
    if(rank==0) print *, 'PASS_NATIVE_MOL_ROOT_REPLAY',iterations
  enddo
  call MPI_Finalize(status)
contains
  subroutine require(local,message)
    logical, intent(in) :: local
    character(len=*), intent(in) :: message
    logical :: accepted
    call MPI_Allreduce(local,accepted,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,status)
    if(.not.accepted) then
      if(.not.local) print *, 'FAIL:',rank,message
      call MPI_Abort(MPI_COMM_WORLD,1,status)
    endif
  end subroutine require
  subroutine zero_field(field,ncomp)
    type(amr_tree_state_3d), intent(out) :: field
    integer, intent(in) :: ncomp
    integer :: node,d(3)
    allocate(field%fields(3))
    do node=1,3
      d=topology%nodes(node)%shape
      allocate(field%fields(node)%values(ncomp,d(1),d(2),d(3)),source=0._dp)
    enddo
  end subroutine zero_field
  subroutine localize(state,stored,gas,temp,held)
    type(eb_reactive_tree_3d), intent(in) :: state
    type(amr_tree_state_3d), intent(in) :: stored
    type(amr_tree_state_3d), intent(out) :: gas,temp,held
    integer :: node,d(3)
    allocate(gas%fields(3),temp%fields(3),held%fields(3))
    do node=1,3
      if(owners(node)/=rank) cycle
      d=topology%nodes(node)%shape
      gas%fields(node)%values=state%nodes(node)%state
      allocate(temp%fields(node)%values(1,d(1),d(2),d(3)))
      temp%fields(node)%values(1,:,:,:)=state%nodes(node)%temperature
      held%fields(node)%values=stored%fields(node)%values
    enddo
  end subroutine localize
  subroutine localize_field(input,output)
    type(amr_tree_state_3d), intent(in) :: input
    type(amr_tree_state_3d), intent(out) :: output
    integer :: node
    allocate(output%fields(3))
    do node=1,3
      if(owners(node)==rank) output%fields(node)=input%fields(node)
    enddo
  end subroutine localize_field
  subroutine serial_source(state,interval,stage,depth,time,rate,ok)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: ok
    integer :: node
    source_calls=source_calls+1;ok=.not.forbid_sources
    do node=1,size(rate%fields)
      rate%fields(node)%values=0
      if(topology%nodes(node)%level==depth) rate%fields(node)%values(5,:,:,:)=1e6_dp*real(1+stage+depth,dp)
    enddo
  end subroutine serial_source
  subroutine native_mpi_source(state,temp,interval,stage,depth,time,rate,comm,ok)
    type(amr_tree_state_3d), intent(in) :: state,temp
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: node
    source_calls=source_calls+1;ok=.not.forbid_sources
    do node=1,size(rate%fields)
      if(.not.allocated(rate%fields(node)%values)) cycle
      rate%fields(node)%values=0
      if(topology%nodes(node)%level==depth) rate%fields(node)%values(5,:,:,:)=1e6_dp*real(1+stage+depth,dp)
    enddo
  end subroutine native_mpi_source
  subroutine serial_replay(saved,ok,weight)
    type(eb_mol_held_interval_3d), intent(in) :: saved(:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: weight
    audit_interval=saved(1)%interval
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,audit_interval,iterations,.true., &
      [.false.,.false.,.false.,.false.],force,1e-8_dp,1e-12_dp,trial,extended,faces,tb,tw,tch,te,ttheta,ok,tir, &
      mol_mode=.true.,redistribution='none',internal_energy=th,internal_energy_flux=thf,energy_policy=policy, &
      energy_correction=tcorrection,heat_release=theat,held_interval_history=replayed,integrated_flux=tgf, &
      auxiliary_source=serial_source,auxiliary_change=ta,auxiliary_has_internal_energy=.true., &
      internal_energy_pressure_work_change=tp,level_steps=tsteps,source_funding_context=tcontext,native_replay_history=saved, &
      evaluate_reaction=audited_reactor,provider_refused=refused,native_spatial_weight=weight, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
  end subroutine serial_replay
  subroutine mpi_replay(saved,ok,weight)
    type(mpi_eb_mol_interval_3d), intent(in) :: saved(:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: weight
    audit_interval=saved(1)%interval
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,ug,ut, &
      audit_interval,iterations,.true.,[.false.,.false.,.false.,.false.],force,1e-8_dp,1e-12_dp,extended,faces, &
      tb,tw,tch,te,umt,MPI_COMM_WORLD,ok,uir,mol_mode=.true.,redistribution='none', &
      internal_energy=uh,internal_energy_flux=uhf,energy_policy=policy,energy_correction=uc,heat_release=uheat, &
      interval_history=mreplayed,integrated_flux=ugf,auxiliary_source=native_mpi_source,auxiliary_change=ta, &
      auxiliary_has_internal_energy=.true.,internal_energy_pressure_work_change=up,level_steps=usteps, &
      source_funding_context=ucontext,native_replay_history=saved, &
      evaluate_reaction=audited_reactor,provider_refused=refused,native_spatial_weight=weight, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
  end subroutine mpi_replay
  subroutine audited_reactor(table,rxns,energy,source_species,source_energy,duration,rtol,atol, &
      density,temperature,completed_steps,ok,maximum_steps)
    type(nasa7_species), intent(in) :: table(:)
    type(elementary_reaction), intent(in) :: rxns(:)
    real(dp), intent(in) :: energy,source_species(:),source_energy,duration,rtol,atol
    real(dp), intent(inout) :: density(:),temperature
    integer, intent(out) :: completed_steps
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    if((refuse_chemistry.and.duration<audit_interval).or. &
      (refuse_after_first_root.and.duration==audit_interval.and.root_reactor_calls>0)) then
      fine_refusals=fine_refusals+1;completed_steps=0;ok=.false.;return
    endif
    if(duration==audit_interval) root_reactor_calls=root_reactor_calls+1
    call advance_forced_reactor_interval(table,rxns,energy,source_species,source_energy,duration,rtol,atol, &
      density,temperature,completed_steps,ok,maximum_steps=maximum_steps)
  end subroutine audited_reactor
  logical function face_equal(a,b) result(same)
    type(eb_integrated_flux_3d), intent(in) :: a,b
    same=all(a%x==b%x).and.all(a%y==b%y).and.all(a%z==b%z)
  end function face_equal
  logical function serial_equal() result(same)
    integer :: node
    same=tcorrection==correction.and.tp==pressure.and.all(tsteps==steps).and. &
      all(tb==b).and.all(tw==w).and.all(tch==ch).and.all(te==ext).and.all(ta==aux)
    do node=1,3
      same=same.and.all(trial%nodes(node)%state==serial%nodes(node)%state).and. &
        all(trial%nodes(node)%temperature==serial%nodes(node)%temperature).and. &
        all(th%fields(node)%values==h%fields(node)%values).and.all(tir%fields(node)%values==ir%fields(node)%values).and. &
        all(theat%fields(node)%values==heat%fields(node)%values).and. &
        face_equal(tgf(node),gf(node)).and.face_equal(thf(node),hf(node))
    enddo
  end function serial_equal
  logical function owner_equal() result(same)
    integer :: node
    same=mc==correction.and.mp==pressure.and.all(msteps==steps).and.all(tb==b).and. &
      all(tw==w).and.all(tch==ch).and.all(te==ext).and.all(ta==aux)
    do node=1,3
      if(owners(node)/=rank) cycle
      same=same.and.all(gas%fields(node)%values==serial%nodes(node)%state).and. &
        all(temp%fields(node)%values(1,:,:,:)==serial%nodes(node)%temperature).and. &
        all(mh%fields(node)%values==h%fields(node)%values).and.all(mir%fields(node)%values==ir%fields(node)%values).and. &
        all(mheat%fields(node)%values==heat%fields(node)%values).and. &
        face_equal(mgf(node),gf(node)).and.face_equal(mhf(node),hf(node))
    enddo
  end function owner_equal
  logical function mpi_equal() result(same)
    integer :: node
    same=uc==mc.and.up==mp.and.all(usteps==msteps).and.all(tb==b).and.all(tw==w).and.all(tch==ch).and. &
      all(te==ext).and.all(ta==aux)
    do node=1,3
      if(owners(node)/=rank) cycle
      same=same.and.all(ug%fields(node)%values==gas%fields(node)%values).and. &
        all(ut%fields(node)%values==temp%fields(node)%values).and.all(uh%fields(node)%values==mh%fields(node)%values).and. &
        all(uir%fields(node)%values==mir%fields(node)%values).and.all(uheat%fields(node)%values==mheat%fields(node)%values).and. &
        face_equal(ugf(node),mgf(node)).and.face_equal(uhf(node),mhf(node))
    enddo
  end function mpi_equal
  subroutine check_inventory()
    real(dp) :: before(nv),after(nv),residual(nv),error
    logical :: ok
    call integrate_reactive_eb_tree_3d(initial,before,ok);call require(ok,'initial independent inventory')
    call integrate_reactive_eb_tree_3d(serial,after,ok);call require(ok,'final independent inventory')
    residual=after-before-b-w-ch-ext-aux(:nv)
    residual(5)=residual(5)-(correction-prefix)
    error=maxval(abs(residual)/max(1._dp,abs(before)))
    call require(error<3e-14_dp,'independent gas/species/source/numerical inventory closes')
    if(rank==0) print *, 'NATIVE_MOL_ROOT_INVENTORY',iterations,error
  end subroutine check_inventory

  subroutine check_subparent_windows()
    type(eb_reactive_tree_3d) :: base
    type(amr_tree_state_3d) :: base_h,base_ir,base_heat
    type(eb_held_funding_context_3d) :: base_context
    type(eb_mol_held_interval_3d), allocatable :: window(:)
    type(mpi_eb_mol_interval_3d), allocatable :: mwindow(:)
    type(eb_integrated_flux_3d), allocatable :: unit_gf(:),unit_hf(:)
    real(dp), parameter :: weights(3)=[1._dp,.5_dp,0._dp]
    real(dp) :: base_c,base_p,weight,gas_difference,held_difference
    integer :: window_index,index,first,last,next,rec,stage,node,depth,bucket
    integer(int64) :: expected_steps(3)
    logical :: local(3),same
    base=serial;base_h=h;base_ir=ir;base_heat=heat;base_context=context;base_c=correction;base_p=pressure
    do layout=1,2
      owners=[0,2,0];if(layout==2) owners=0
      call localize(initial,initial_h,gas,temp,mh);call localize_field(records(1)%clock_reaction_source,mir)
      call localize_field(records(1)%clock_heat,mheat)
      mcontext=eb_held_funding_context_3d();mc=prefix;mp=0;source_calls=0;forbid_sources=.false.
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
        dt,iterations,.true.,[.false.,.false.,.false.,.false.],force,1e-8_dp,1e-12_dp,extended,faces, &
        tb,tw,tch,te,mt,MPI_COMM_WORLD,good,mir,mol_mode=.true.,redistribution='none', &
        internal_energy=mh,internal_energy_flux=mhf,energy_policy=policy,energy_correction=mc,heat_release=mheat, &
        interval_history=mrecords,integrated_flux=mgf,auxiliary_source=native_mpi_source,auxiliary_change=ta, &
        auxiliary_has_internal_energy=.true.,internal_energy_pressure_work_change=mp,level_steps=msteps, &
        source_funding_context=mcontext, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
      call require(good,'native MPI source windows captured')
      call check_mpi_ancestors(records,mrecords)
      saved_calls=source_calls;forbid_sources=.true.
      do window_index=1,2
        first=2;last=4;next=5;expected_steps=[0_int64,1_int64,2_int64]
        if(window_index==2) then
          first=6;last=6;next=7;expected_steps=[0_int64,0_int64,1_int64]
        endif
        window=records(first:last);mwindow=mrecords(first:last);depth=window(1)%depth
        do index=1,3
          weight=weights(index)
          ! A completed caller tail is deliberately unrelated to this open
          ! window. Only the saved entry/ancestor/context is authoritative.
          trial=base;th=base_h;tir=base_ir;theat=base_heat;tcontext=base_context;tcorrection=999;tp=999
          trial%nodes(1)%state(2,:,:,:)=trial%nodes(1)%state(2,:,:,:)+.01_dp
          call serial_replay(window,good,weight)
          call require(good.and.source_calls==saved_calls.and.all(tsteps==expected_steps), &
            'serial sub-parent replays precisely its physical clocks with frozen sources')
          if(weight==1) then
            same=tcorrection==records(next)%clock_energy_correction
            do node=1,3
              same=same.and.all(trial%nodes(node)%state==records(next)%clock_initial%nodes(node)%state).and. &
                all(trial%nodes(node)%temperature==records(next)%clock_initial%nodes(node)%temperature).and. &
                all(th%fields(node)%values==records(next)%clock_internal_energy%fields(node)%values).and. &
                all(tir%fields(node)%values==records(next)%clock_reaction_source%fields(node)%values).and. &
                all(theat%fields(node)%values==records(next)%clock_heat%fields(node)%values)
            enddo
            call require(same,'unit sub-parent result matches the next native clock anchor bitwise')
            call require(context_equal(tcontext,records(next)%clock_source_funding_context), &
              'unit sub-parent pending stock and every still-open parent face match native continuation')
            unit_gf=tgf;unit_hf=thf
          endif
          gas_difference=0;held_difference=0
          do node=1,3
            gas_difference=max(gas_difference,face_difference(tgf(node),unit_gf(node)))
            held_difference=max(held_difference,face_difference(thf(node),unit_hf(node)))
          enddo
          if(weight<1) call require(gas_difference>1e-12_dp.and.held_difference>1e-12_dp, &
            'sub-parent weights change both conservative face families')
          local=.true.
          call require(valid_eb_held_funding_context_3d(topology,geometry,nv,local,tcontext,depth), &
            'closed sub-window keeps exactly its external ancestors open')
          do rec=1,size(window)
            do stage=0,iterations
              do node=1,3
                call require(all(replayed(rec)%external_stages(stage)%fields(node)%values== &
                  window(rec)%external_stages(stage)%fields(node)%values).and. &
                  all(replayed(rec)%auxiliary_stages(stage)%fields(node)%values== &
                  window(rec)%auxiliary_stages(stage)%fields(node)%values),'sub-parent volume stages frozen')
              enddo
            enddo
          enddo
          call check_subparent_inventory(window(1))
          serial=trial;h=th;ir=tir;heat=theat;correction=tcorrection;pressure=tp;gf=tgf;hf=thf;context=tcontext
          b=tb;w=tw;ch=tch;ext=te;aux=ta;steps=tsteps
          ug=gas;ut=temp;uh=mh;uir=mir;uheat=mheat;ucontext=mcontext;uc=999;up=999
          call mpi_replay(mwindow,good,weight)
          call require(good.and.source_calls==saved_calls.and.all(usteps==expected_steps), &
            'MPI sub-parent physical clocks and frozen source providers agree')
          call check_mpi_ancestors(replayed,mreplayed)
          gas=ug;temp=ut;mh=uh;mir=uir;mheat=uheat;mc=uc;mp=up;mgf=ugf;mhf=uhf;msteps=usteps;mcontext=ucontext
          call require(owner_equal(),'sub-parent gas/T/H/IR/heat/face/budget outputs bitwise serial/MPI')
          call require(context_equal(ucontext,tcontext,owners==rank), &
            'sub-parent owner-local pending/open families bitwise serial/MPI')
          local=owners==rank
          call require(valid_eb_held_funding_context_3d(topology,geometry,nv,local,ucontext,depth), &
            'MPI open ancestors retained with the complete owner/idle layout')
          do case_id=1,5
            mbad=mwindow
            select case(case_id)
            case(1)
              if(rank==1) deallocate(mbad(size(mbad))%clock_interval)
            case(2)
              if(rank==1) mbad(size(mbad))%clock_start(1)=dt
            case(3)
              if(rank==owners(1)) mbad(size(mbad))%ancestor_internal_energy_end%fields(1)%values(1,1,3,1)= &
                ieee_value(0._dp,ieee_quiet_nan)
            case(4)
              mbad(1)%clock_interval(1)=1.1_dp*mbad(1)%clock_interval(1)
            case(5)
              if(rank==owners(1)) deallocate(mbad(1)%clock_source_funding_context%open(0)%gas(1)%x)
            end select
            root_reactor_calls=0;fine_refusals=0
            call mpi_replay(mbad,good,weight)
            call require(.not.good.and.root_reactor_calls==0.and.fine_refusals==0.and.source_calls==saved_calls, &
              'invalid late clock/ancestor/open-face metadata rejects collectively before physics')
            call require(mpi_equal().and.context_equal(ucontext,mcontext), &
              'sub-parent metadata refusal preserves fields, budgets and every open parent family')
          enddo
          bad=window;bad(1)%clock_interval(1)=1.1_dp*bad(1)%clock_interval(1)
          root_reactor_calls=0;fine_refusals=0
          call serial_replay(bad,good,weight)
          call require(.not.good.and.root_reactor_calls==0.and.serial_equal().and.context_equal(tcontext,context), &
            'serial invalid ancestor period rejects atomically before physics')
          refuse_chemistry=depth==1;refuse_after_first_root=depth==2
          root_reactor_calls=0;fine_refusals=0
          call mpi_replay(mwindow,good,weight)
          call require(.not.good.and.refused.and.mpi_equal().and.context_equal(ucontext,mcontext), &
            'sub-parent admitted provider refusal rolls back fields and pending/open parent families')
          call MPI_Allreduce(root_reactor_calls,global_root_calls,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
          call MPI_Allreduce(fine_refusals,global_refusals,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
          call require(global_root_calls>0.and.global_refusals==1,'one admitted sub-parent refusal without retry')
          root_reactor_calls=0;fine_refusals=0
          call serial_replay(window,good,weight)
          call require(.not.good.and.refused.and.root_reactor_calls>0.and.fine_refusals==1.and. &
            serial_equal().and.context_equal(tcontext,context),'serial admitted sub-parent refusal rolls back entire window')
          refuse_chemistry=.false.;refuse_after_first_root=.false.
          if(rank==0) print *, 'PASS_NATIVE_MOL_SUBPARENT_REPLAY',iterations,layout,depth,window(1)%time,weight
        enddo
      enddo
    enddo
  end subroutine check_subparent_windows

  real(dp) function face_difference(a,b) result(value)
    type(eb_integrated_flux_3d), intent(in) :: a,b
    value=max(maxval(abs(a%x-b%x)),maxval(abs(a%y-b%y)),maxval(abs(a%z-b%z)))
  end function face_difference

  logical function context_equal(a,b,local) result(same)
    type(eb_held_funding_context_3d), intent(in) :: a,b
    logical, optional, intent(in) :: local(:)
    integer :: bucket,node
    logical :: held_a,held_b
    same=all(a%open%clock==b%open%clock)
    do bucket=0,2
      do node=1,3
        if(present(local)) then
          if(.not.local(node)) cycle
        endif
        held_a=allocated(a%pending(bucket)%fields(node)%values);held_b=allocated(b%pending(bucket)%fields(node)%values)
        same=same.and.(held_a.eqv.held_b)
        if(held_a.and.held_b) same=same.and.all(a%pending(bucket)%fields(node)%values==b%pending(bucket)%fields(node)%values)
        if(a%open(bucket)%clock<0.or.b%open(bucket)%clock<0) cycle
        held_a=allocated(a%open(bucket)%gas(node)%x);held_b=allocated(b%open(bucket)%gas(node)%x)
        ! Comparing a serial context with an MPI context ignores absent idle
        ! payloads; owner layout validity is checked separately above.
        if(held_a.and.held_b) same=same.and.face_equal(a%open(bucket)%gas(node),b%open(bucket)%gas(node)).and. &
          face_equal(a%open(bucket)%scalar(node),b%open(bucket)%scalar(node))
      enddo
    enddo
  end function context_equal

  subroutine check_subparent_inventory(anchor)
    type(eb_mol_held_interval_3d), intent(in) :: anchor
    type(eb_reactive_tree_3d) :: ledger
    real(dp) :: before(nv),after(nv),faces(nv),residual(nv),hbefore(1),hafter(1),hfaces(1),volume,reset,error,herror
    integer :: node,x,y,z,rec
    logical :: ok
    call integrate_reactive_eb_tree_3d(anchor%clock_initial,before,ok);call require(ok,'sub-window entry inventory')
    call integrate_reactive_eb_tree_3d(trial,after,ok);call require(ok,'sub-window closing inventory')
    ledger=initial
    do node=1,3
      ledger%nodes(node)%state=0
      associate(g=>geometry(node),f=>tgf(node))
      do z=1,g%nz
      do y=1,g%ny
      do x=1,g%nx
        if(g%volume_fraction(x,y,z)<=0) cycle
        volume=g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz
        ledger%nodes(node)%state(:,x,y,z)=( &
          g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(:,x-1,y,z)-g%x_face_fraction(x,y,z)*f%x(:,x,y,z))+ &
          g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(:,x,y-1,z)-g%y_face_fraction(x,y,z)*f%y(:,x,y,z))+ &
          g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(:,x,y,z-1)-g%z_face_fraction(x,y,z)*f%z(:,x,y,z)))/volume
      enddo
      enddo
      enddo
      end associate
    enddo
    call integrate_reactive_eb_tree_3d(ledger,faces,ok);call require(ok,'open-window composite gas face impulse')
    residual=after-before-faces-tw-tch-te-ta(:nv);residual(5)=residual(5)-(tcorrection-anchor%clock_energy_correction)
    error=maxval(abs(residual)/max(1._dp,abs(before)))
    do node=1,3
      ledger%nodes(node)%state=anchor%clock_internal_energy%fields(node)%values
    enddo
    call integrate_reactive_eb_tree_3d(ledger,hbefore,ok);call require(ok,'sub-window entry actual UEINT inventory')
    do node=1,3
      ledger%nodes(node)%state=th%fields(node)%values
    enddo
    call integrate_reactive_eb_tree_3d(ledger,hafter,ok);call require(ok,'sub-window closing actual UEINT inventory')
    do node=1,3
      ledger%nodes(node)%state=0
      associate(g=>geometry(node),f=>thf(node))
      do z=1,g%nz
      do y=1,g%ny
      do x=1,g%nx
        if(g%volume_fraction(x,y,z)<=0) cycle
        volume=g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz
        ledger%nodes(node)%state(1,x,y,z)=( &
          g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(1,x-1,y,z)-g%x_face_fraction(x,y,z)*f%x(1,x,y,z))+ &
          g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(1,x,y-1,z)-g%y_face_fraction(x,y,z)*f%y(1,x,y,z))+ &
          g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(1,x,y,z-1)-g%z_face_fraction(x,y,z)*f%z(1,x,y,z)))/volume
      enddo
      enddo
      enddo
      end associate
    enddo
    call integrate_reactive_eb_tree_3d(ledger,hfaces,ok);call require(ok,'open-window composite actual UEINT face impulse')
    reset=0
    do rec=1,size(replayed)
      reset=reset+replayed(rec)%internal_energy_projection_change+replayed(rec)%internal_energy_closure_change+ &
        replayed(rec)%held_clock_start_change+replayed(rec)%held_clock_end_change
    enddo
    herror=abs(hafter(1)-hbefore(1)-hfaces(1)-ta(nv+1)-reset)/max(1._dp,abs(hbefore(1)))
    if(rank==0) print *, 'NATIVE_MOL_SUBPARENT_INVENTORY',iterations,anchor%depth,error,herror
    call require(error<3e-14_dp.and.herror<3e-14_dp,'open-window gas/species/actual UEINT face/source/reset inventories close')
  end subroutine check_subparent_inventory

  subroutine check_weighted_window()
    type(eb_reactive_tree_3d) :: base
    type(amr_tree_state_3d) :: base_h,base_ir,base_heat
    type(eb_integrated_flux_3d), allocatable :: base_gf(:),base_hf(:)
    type(eb_held_funding_context_3d) :: base_context
    real(dp) :: base_b(nv),base_w(nv),base_ch(nv),base_ext(nv),base_aux(nv+1)
    real(dp) :: base_c,base_p,base_theta,weight,bad_weight,gas_difference,held_difference
    real(dp), parameter :: weights(3)=[1._dp,.5_dp,0._dp]
    integer :: index,rec,stage,node
    base=serial;base_h=h;base_ir=ir;base_heat=heat;base_context=context;base_gf=gf;base_hf=hf
    base_c=correction;base_p=pressure;base_theta=theta
    base_b=b;base_w=w;base_ch=ch;base_ext=ext;base_aux=aux
    do layout=1,2
      owners=[0,2,0];if(layout==2) owners=0
      call localize(initial,initial_h,gas,temp,mh);call localize_field(records(1)%clock_reaction_source,mir)
      call localize_field(records(1)%clock_heat,mheat)
      mcontext=eb_held_funding_context_3d();mc=prefix;mp=0;source_calls=0;forbid_sources=.false.
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
        dt,iterations,.true.,[.false.,.false.,.false.,.false.],force,1e-8_dp,1e-12_dp,extended,faces, &
        tb,tw,tch,te,mt,MPI_COMM_WORLD,good,mir,mol_mode=.true.,redistribution='none', &
        internal_energy=mh,internal_energy_flux=mhf,energy_policy=policy,energy_correction=mc,heat_release=mheat, &
        interval_history=mrecords,integrated_flux=mgf,auxiliary_source=native_mpi_source,auxiliary_change=ta, &
        auxiliary_has_internal_energy=.true.,internal_energy_pressure_work_change=mp,level_steps=msteps, &
        source_funding_context=mcontext, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
      call require(good,'nonuniform native MPI archive accepts')
      saved_calls=source_calls;forbid_sources=.true.
      do index=1,3
        weight=weights(index)
        serial=base;h=base_h;ir=base_ir;heat=base_heat;context=base_context;gf=base_gf;hf=base_hf
        correction=base_c;pressure=base_p;theta=base_theta;b=base_b;w=base_w;ch=base_ch;ext=base_ext;aux=base_aux
        trial=serial;th=h;tir=ir;theat=heat;tcorrection=999;tp=999;tcontext=context
        call serial_replay(records,good,weight)
        call require(good.and.source_calls==saved_calls,'weighted serial replay freezes external providers')
        if(weight==1) call require(serial_equal(),'explicit unit weight is bitwise native')
        gas_difference=0;held_difference=0
        do node=1,3
          gas_difference=max(gas_difference,maxval(abs(tgf(node)%x-gf(node)%x)))
          held_difference=max(held_difference,maxval(abs(thf(node)%x-hf(node)%x)))
        enddo
        if(weight<1) call require(gas_difference>1e-12_dp.and.held_difference>1e-12_dp, &
          'nonunit weight changes both gas and actual UEINT face families')
        call require(ttheta<=weight,'reported common weight bounds the accepted family')
        do rec=1,7
          do stage=0,iterations
            do node=1,3
              call require(all(replayed(rec)%external_stages(stage)%fields(node)%values== &
                records(rec)%external_stages(stage)%fields(node)%values).and. &
                all(replayed(rec)%auxiliary_stages(stage)%fields(node)%values== &
                records(rec)%auxiliary_stages(stage)%fields(node)%values),'all volume forcing stages remain frozen')
            enddo
          enddo
        enddo
        serial=trial;h=th;ir=tir;heat=theat;correction=tcorrection;pressure=tp;theta=ttheta
        gf=tgf;hf=thf;context=tcontext;b=tb;w=tw;ch=tch;ext=te;aux=ta
        call check_inventory()
        call check_weighted_held_inventory()
        ug=gas;ut=temp;uh=mh;uir=mir;uheat=mheat;uc=999;up=999;ucontext=mcontext
        call mpi_replay(mrecords,good,weight)
        call require(good.and.source_calls==saved_calls,'weighted MPI replay freezes sources on idle ranks')
        call check_mpi_ancestors(replayed,mreplayed)
        gas=ug;temp=ut;mh=uh;mir=uir;mheat=uheat;mc=uc;mp=up;mgf=ugf;mhf=uhf;msteps=usteps;mcontext=ucontext
        call require(owner_equal(),'each weighted gas/T/H/IR/heat/face/ledger is bitwise serial/MPI')
        do case_id=1,4
          bad_weight=weight
          select case(case_id)
          case(1)
            if(rank==1) bad_weight=merge(.875_dp,.25_dp,weight==1)
          case(2)
            if(rank==1) bad_weight=ieee_value(0._dp,ieee_quiet_nan)
          case(3)
            if(rank==1) bad_weight=-.25_dp
          end select
          root_reactor_calls=0;fine_refusals=0
          if(case_id==4.and.rank==1) then
            call mpi_replay(mrecords,good)
          else
            call mpi_replay(mrecords,good,bad_weight)
          endif
          call require(.not.good.and.source_calls==saved_calls.and.root_reactor_calls==0.and.fine_refusals==0, &
            'idle rank mismatched/NaN/range/absent weight refuses before chemistry or sources')
          call require(mpi_equal(),'invalid collective weight retains accepted weighted outputs')
        enddo
        refuse_chemistry=.true.;root_reactor_calls=0;fine_refusals=0
        call mpi_replay(mrecords,good,weight)
        call require(.not.good.and.refused,'weighted MPI child chemistry refusal aborts transaction')
        call MPI_Allreduce(root_reactor_calls,global_root_calls,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
        call MPI_Allreduce(fine_refusals,global_refusals,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
        call require(global_root_calls>0.and.global_refusals==1,'weighted root completes before one child refusal')
        call require(mpi_equal(),'weighted refusal retains complete accepted caller result')
        refuse_chemistry=.false.
        if(layout==1) then
          do case_id=1,2
            bad_weight=-.25_dp
            if(case_id==2) bad_weight=ieee_value(0._dp,ieee_quiet_nan)
            root_reactor_calls=0;fine_refusals=0
            call serial_replay(records,good,bad_weight)
            call require(.not.good.and.root_reactor_calls==0.and.fine_refusals==0.and.serial_equal(), &
              'invalid serial weight refuses before chemistry and retains outputs')
          enddo
          refuse_chemistry=.true.;root_reactor_calls=0;fine_refusals=0
          call serial_replay(records,good,weight)
          call require(.not.good.and.refused.and.root_reactor_calls>0.and.fine_refusals==1.and.serial_equal(), &
            'weighted serial root completes then refusing child rolls back every output')
          refuse_chemistry=.false.
        endif
        if(rank==0) print *, 'PASS_NATIVE_MOL_WEIGHTED_REPLAY',iterations,layout,weight,gas_difference,held_difference
      enddo
    enddo
  end subroutine check_weighted_window

  subroutine check_weighted_held_inventory()
    type(eb_reactive_tree_3d) :: ledger
    real(dp) :: before(1),after(1),advection(1),volume,reset,error
    integer :: node,x,y,z,rec
    logical :: ok
    ledger=initial
    do node=1,3
      ledger%nodes(node)%state=initial_h%fields(node)%values
    enddo
    call integrate_reactive_eb_tree_3d(ledger,before,ok);call require(ok,'initial independent held inventory')
    do node=1,3
      ledger%nodes(node)%state=h%fields(node)%values
    enddo
    call integrate_reactive_eb_tree_3d(ledger,after,ok);call require(ok,'final independent held inventory')
    do node=1,3
      ledger%nodes(node)%state=0
      associate(g=>geometry(node),f=>hf(node))
      do z=1,g%nz
      do y=1,g%ny
      do x=1,g%nx
        if(g%volume_fraction(x,y,z)<=0) cycle
        volume=g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz
        ledger%nodes(node)%state(1,x,y,z)=( &
          g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(1,x-1,y,z)-g%x_face_fraction(x,y,z)*f%x(1,x,y,z))+ &
          g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(1,x,y-1,z)-g%y_face_fraction(x,y,z)*f%y(1,x,y,z))+ &
          g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(1,x,y,z-1)-g%z_face_fraction(x,y,z)*f%z(1,x,y,z)))/volume
      enddo
      enddo
      enddo
      end associate
    enddo
    call integrate_reactive_eb_tree_3d(ledger,advection,ok);call require(ok,'composite weighted held face exchange')
    reset=0
    do rec=1,size(replayed)
      ! The reacting held projection already owns the entry reset and the
      ! gas-forced thermal change. These are not extra UEINT volume sources.
      reset=reset+replayed(rec)%internal_energy_projection_change+replayed(rec)%internal_energy_closure_change+ &
        replayed(rec)%held_clock_start_change+replayed(rec)%held_clock_end_change
    enddo
    error=abs(after(1)-before(1)-advection(1)-aux(nv+1)-reset)/max(1._dp,abs(before(1)))
    if(rank==0) print *, 'NATIVE_MOL_WEIGHTED_HELD_INVENTORY',iterations,error
    call require(error<3e-14_dp,'independent weighted actual UEINT face/source/reset inventory closes')
  end subroutine check_weighted_held_inventory
  subroutine check_mpi_ancestors(serial_records,mpi_records)
    type(eb_mol_held_interval_3d), intent(in) :: serial_records(:)
    type(mpi_eb_mol_interval_3d), intent(in) :: mpi_records(:)
    integer :: rec,node
    logical :: same,owned
    call require(size(serial_records)==size(mpi_records),'same physical window record count')
    do rec=1,size(serial_records)
      call require(all(serial_records(rec)%clock_start==mpi_records(rec)%clock_start).and. &
        all(serial_records(rec)%clock_interval==mpi_records(rec)%clock_interval), &
        'owner-neutral actual ancestor clocks retained bitwise')
      call require(size(mpi_records(rec)%fractions_start)==size(serial_records(rec)%fractions_start).and. &
        size(mpi_records(rec)%fractions_end)==size(serial_records(rec)%fractions_end),'native ancestor fractions shapes exact')
      call require(all(mpi_records(rec)%fractions_start==serial_records(rec)%fractions_start).and. &
        all(mpi_records(rec)%fractions_end==serial_records(rec)%fractions_end),'serial/MPI ancestor fractions exact')
      same=.true.
      do node=1,3
        owned=owners(node)==rank
        same=same.and.(allocated(mpi_records(rec)%ancestor_start%fields(node)%values).eqv.owned).and. &
          (allocated(mpi_records(rec)%ancestor_end%fields(node)%values).eqv.owned).and. &
          (allocated(mpi_records(rec)%ancestor_temperature_start%fields(node)%values).eqv.owned).and. &
          (allocated(mpi_records(rec)%ancestor_temperature_end%fields(node)%values).eqv.owned).and. &
          (allocated(mpi_records(rec)%ancestor_internal_energy_start%fields(node)%values).eqv.owned).and. &
          (allocated(mpi_records(rec)%ancestor_internal_energy_end%fields(node)%values).eqv.owned).and. &
          (allocated(mpi_records(rec)%ancestor_source_start%fields(node)%values).eqv.owned).and. &
          (allocated(mpi_records(rec)%ancestor_source_end%fields(node)%values).eqv.owned)
      enddo
      call require(same,'all ancestor fields published with full owner/idle allocation')
      same=.true.
      do node=1,3
        if(owners(node)/=rank) cycle
        same=same.and.all(mpi_records(rec)%ancestor_start%fields(node)%values== &
          serial_records(rec)%ancestor_start%nodes(node)%state).and. &
          all(mpi_records(rec)%ancestor_end%fields(node)%values==serial_records(rec)%ancestor_end%nodes(node)%state).and. &
          all(mpi_records(rec)%ancestor_temperature_start%fields(node)%values(1,:,:,:)== &
          serial_records(rec)%ancestor_start%nodes(node)%temperature).and. &
          all(mpi_records(rec)%ancestor_temperature_end%fields(node)%values(1,:,:,:)== &
          serial_records(rec)%ancestor_end%nodes(node)%temperature).and. &
          all(mpi_records(rec)%ancestor_internal_energy_start%fields(node)%values== &
          serial_records(rec)%ancestor_internal_energy_start%fields(node)%values).and. &
          all(mpi_records(rec)%ancestor_internal_energy_end%fields(node)%values== &
          serial_records(rec)%ancestor_internal_energy_end%fields(node)%values).and. &
          all(mpi_records(rec)%ancestor_source_start%fields(node)%values== &
          serial_records(rec)%ancestor_source_start%fields(node)%values).and. &
          all(mpi_records(rec)%ancestor_source_end%fields(node)%values==serial_records(rec)%ancestor_source_end%fields(node)%values)
      enddo
      call require(same,'serial/MPI full ancestor gas/T/H/source cache exact after explicit mover')
    enddo
  end subroutine check_mpi_ancestors
  subroutine check_ancestors()
    integer :: rec,ancestor,node,depth,parent
    do rec=2,size(records)
      depth=records(rec)%depth
      call require(size(records(rec)%fractions_start)==depth.and.size(records(rec)%fractions_end)==depth, &
        'actual interpolation dimensions retained')
      do ancestor=0,depth-1
        parent=rec-1
        do while(records(parent)%depth/=ancestor)
          parent=parent-1
        enddo
        do node=1,3
          if(topology%nodes(node)%level/=ancestor) cycle
          call require(all(records(rec)%ancestor_start%nodes(node)%state==records(parent)%initial%nodes(node)%state).and. &
            all(records(rec)%ancestor_end%nodes(node)%state==records(parent)%local_tree%nodes(node)%state), &
            'still-open ancestor endpoints retained before descendant average-down')
          call require(all(records(rec)%ancestor_internal_energy_start%fields(node)%values== &
            records(parent)%initial_internal_energy%fields(node)%values).and. &
            all(records(rec)%ancestor_internal_energy_end%fields(node)%values== &
            records(parent)%local_internal_energy%fields(node)%values),'actual ancestor UEINT endpoints retained')
        enddo
      enddo
    enddo
  end subroutine check_ancestors
end program test_mpi_eb_native_mol_replay_3d
