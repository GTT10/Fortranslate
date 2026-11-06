program test_mpi_eb_native_mol_partial_recovery_3d
  ! Fixed inert-N2 density-advection/source competition, with no fine donors.
  ! The original native clock refuses; one source-frozen retry must recover.
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: load_h2o2_full_transport
  use gas_transport_mod, only: gas_transport_species
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use dual_energy_mod, only: dual_energy_policy
  use amr_hierarchy_3d_mod, only: amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    build_eb_tree_extended_geometry_3d,advance_reactive_eb_tree_sdc_3d,advance_reactive_eb_tree_sdc_subcycle_3d, &
    integrate_reactive_eb_tree_3d,eb_mol_held_interval_3d,eb_native_mol_source_prefix_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_held_funding_queue_3d_mod, only: eb_held_funding_context_3d,valid_eb_held_funding_context_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_subcycle_3d,mpi_eb_mol_interval_3d
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_tree_3d) :: topology
  type(amr_patch_3d) :: empty_patches(0)
  type(eb_geometry_3d) :: fine,padded,geometry(1),extended_parent
  type(eb_geometry_3d), allocatable :: extended(:)
  type(eb_reactive_tree_3d) :: initial,reference
  type(amr_tree_state_3d) :: initial_h,reference_h,reference_ir
  type(eb_mol_held_interval_3d), allocatable :: reference_records(:)
  type(eb_integrated_flux_3d), allocatable :: reference_faces(:),reference_scalar(:)
  type(eb_held_funding_context_3d) :: reference_context
  type(eb_native_mol_source_prefix_3d) :: refused_prefix
  type(dual_energy_policy) :: policy
  real(dp), allocatable :: phi(:,:,:),padded_phi(:,:,:),u(:,:,:,:),temperature_field(:,:,:),primitive(:),cell(:)
  real(dp), allocatable :: reference_ledgers(:)
  real(dp), parameter :: duration=1.e-8_dp,sink=1.002_dp
  logical, parameter :: transport_flags(4)=[.false.,.false.,.false.,.false.]
  real(dp), parameter :: force(3)=[0._dp,0._dp,0._dp]
  character(len=16) :: boundaries(6)='outflow'
  integer :: rank,nranks,status,nv,naux,nitrogen,p,i,j,k,owner(1),layout
  integer :: source_calls=0,old_calls=0,reactor_calls=0,callback_mode=0
  logical :: good
  real(dp) :: density,rg,temperature,speed
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call check(nranks==3,'owner and two idle ranks')
  call load_h2o2_full_thermo(species,good);call check(good,'pinned native thermodynamics')
  call load_h2o2_full_mechanism(reactions,good);call check(good,'pinned native reactions')
  call load_h2o2_full_transport(transport,good);call check(good,'pinned native transport')
  species%pinned_eos_extrapolation=.true.
  nv=reactive_nvar(size(species));naux=nv+1;allocate(primitive(nv),cell(nv),reference_ledgers(5*nv+3))
  nitrogen=0
  do p=1,size(species)
    if(species(p)%name=='N2') nitrogen=p
  enddo
  call check(nitrogen>0,'native nitrogen species')
  primitive=0;primitive(5+nitrogen)=1
  rg=mixture_specific_gas_constant(species,primitive(6:),good);call check(good,'native gas constant')
  allocate(phi(0:8,0:8,0:8),padded_phi(0:32,0:32,0:32))
  do i=0,8
    phi(i,:,:)=real(i,dp)*.0005_dp-.00185_dp
  enddo
  do i=0,32
    padded_phi(i,:,:)=real(max(0,min(8,i-12)),dp)*.0005_dp-.00185_dp
  enddo
  call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[.004_dp,.004_dp,.004_dp],fine,good)
  call check(good,'native cut geometry')
  call coarsen_eb_geometry_3d(fine,2,geometry(1),good);call check(good,'root geometry from existing fixture')
  call build_level_set_eb_geometry_3d(padded_phi,spread(-.006_dp,1,3),spread(.010_dp,1,3),padded,good)
  call check(good,'padded cut geometry')
  call coarsen_eb_geometry_3d(padded,2,extended_parent,good);call check(good,'root geometry context')
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],empty_patches,topology,good)
  call check(good,'single-root topology')
  allocate(u(nv,4,4,4),temperature_field(4,4,4),source=0._dp)
  do k=1,4
  do j=1,4
  do i=1,4
    if(geometry(1)%volume_fraction(i,j,k)<=0) cycle
    density=2._dp**(3-i)
    primitive(1)=density;primitive(2)=1000;primitive(3:4)=0;primitive(5)=density*rg*1400
    call reactive_primitive_to_conserved(species,primitive,cell,temperature,speed,good)
    call check(good,'physical density-gradient gas')
    u(:,i,j,k)=cell;temperature_field(i,j,k)=temperature
  enddo
  enddo
  enddo
  call initialize_reactive_eb_tree_3d(species,topology,geometry(1),u,temperature_field,initial,good, &
    internal_energy=initial_h)
  call check(good,'physical density-gradient gas and independent H')
  ! The tree extracts node geometry, including its halo layout. Use that same
  ! geometry for MPI and funding-context identity checks, as for the gas owner.
  geometry(1)=initial%nodes(1)%geometry
  call build_eb_tree_extended_geometry_3d(topology,geometry,extended_parent,extended,good,halo_width=6)
  call check(good,'single-root spatial context')
  policy%allow_negative_energy=.false.;policy%eta2=2;policy%allow_small_energy=.false.
  call check_native_refusal()
  call check_serial(0)
  call check_serial(1)
  call check_serial(2)
  call check_serial(3)
  do layout=1,2
    owner=0;if(layout==2) owner=2
    call check_mpi_layout()
    call MPI_Barrier(MPI_COMM_WORLD,status);call check(status==MPI_SUCCESS,'owner layout complete')
    if(rank==0) print *, 'PASS_NATIVE_MOL_PARTIAL_ROOT_RECOVERY',owner
  enddo
  call MPI_Finalize(status)
contains
  subroutine reset_callbacks(mode)
    integer, intent(in) :: mode
    callback_mode=mode;source_calls=0;old_calls=0;reactor_calls=0
  end subroutine reset_callbacks

  subroutine zero_ir(ir,local)
    type(amr_tree_state_3d), intent(out) :: ir
    logical, intent(in) :: local
    allocate(ir%fields(1))
    if(local) allocate(ir%fields(1)%values(nv-4,4,4,4),source=0._dp)
  end subroutine zero_ir

  subroutine make_faces(faces,components)
    type(eb_integrated_flux_3d), allocatable, intent(out) :: faces(:)
    integer, intent(in) :: components
    allocate(faces(1))
    allocate(faces(1)%x(components,0:4,4,4),faces(1)%y(components,4,0:4,4), &
      faces(1)%z(components,4,4,0:4),source=0._dp)
  end subroutine make_faces

  subroutine check_native_refusal()
    type(eb_reactive_tree_3d) :: gas
    type(amr_tree_state_3d) :: held,ir,stock,receipt
    type(eb_integrated_flux_3d), allocatable :: faces(:),scalar(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(naux),theta,correction
    logical :: accepted,ready
    integer :: stage
    gas=initial;held=initial_h;call zero_ir(ir,.true.);call make_faces(faces,nv);call make_faces(scalar,1)
    allocate(stock%fields(1));allocate(stock%fields(1)%values(naux,4,4,4),source=0._dp);receipt=stock
    b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp;correction=7
    call reset_callbacks(0)
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,duration,2,.true.,transport_flags,force, &
      2.e-7_dp,1.e-12_dp,gas,extended,boundaries,b,w,ch,ext,theta,accepted,reaction_source=ir, &
      target_depth=0,mol_mode=.true.,integrated_flux=faces,internal_energy=held,internal_energy_flux=scalar, &
      energy_policy=policy,energy_correction=correction,auxiliary_source=common_source,auxiliary_change=aux, &
      auxiliary_has_internal_energy=.true.,source_funding_stock=stock,source_funding_impulse=receipt, &
      native_source_prefix=refused_prefix,native_source_prefix_ready=ready)
    call check(.not.accepted.and.ready,'original native unit-weight clock refuses and returns a prefix')
    call check(refused_prefix%accepted_source_stages>=1.and.refused_prefix%accepted_source_stages<3, &
      'the original failure has not observed the whole clock')
    call check(old_calls==1.and.source_calls==min(2,refused_prefix%accepted_source_stages), &
      'native failure evaluates only reached source packets')
    do stage=refused_prefix%accepted_source_stages,2
      call check(.not.allocated(refused_prefix%external_stages(stage)%fields).and. &
        .not.allocated(refused_prefix%auxiliary_stages(stage)%fields),'native future packets remain absent')
    enddo
    call check(all(gas%nodes(1)%state==initial%nodes(1)%state).and. &
      all(gas%nodes(1)%temperature==initial%nodes(1)%temperature).and. &
      all(held%fields(1)%values==initial_h%fields(1)%values).and.all(ir%fields(1)%values==0), &
      'original refusing clock preserves caller gas/T/H/I_R')
    call check(all([b,w,ch,ext,aux,theta,correction]== &
      [spread(17._dp,1,nv),spread(19._dp,1,nv),spread(23._dp,1,nv),spread(29._dp,1,nv), &
      spread(31._dp,1,naux),.7_dp,7._dp]),'original refusal preserves caller ledgers')
    call check(all(faces(1)%x==0).and.all(faces(1)%y==0).and.all(faces(1)%z==0).and. &
      all(scalar(1)%x==0).and.all(scalar(1)%y==0).and.all(scalar(1)%z==0).and. &
      all(receipt%fields(1)%values==0),'original refusal preserves faces and receipt')
    if(rank==0) print *, 'PARTIAL_ROOT_INITIAL_OBSERVED',refused_prefix%accepted_source_stages
  end subroutine check_native_refusal

  subroutine check_serial(mode)
    integer, intent(in) :: mode
    type(eb_reactive_tree_3d) :: gas
    type(amr_tree_state_3d) :: held,ir
    type(eb_held_funding_context_3d) :: context
    type(eb_mol_held_interval_3d), allocatable :: records(:),replay(:)
    type(eb_integrated_flux_3d), allocatable :: faces(:),scalar(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(naux),theta,correction
    integer(int64) :: retries,steps(1)
    logical :: accepted
    gas=initial;held=initial_h;call zero_ir(ir,.true.)
    b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp;correction=0;steps=-17
    if(mode==1) replay=reference_records
    if(mode>=2) then
      records=reference_records;faces=reference_faces;scalar=reference_scalar;context=reference_context
    endif
    call reset_callbacks(mode)
    if(mode==3) then
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,duration,2,.true., &
        transport_flags,force,2.e-7_dp,1.e-12_dp,gas,extended,boundaries,b,w,ch,ext,theta,accepted, &
        reaction_source=ir,mol_mode=.true.,redistribution='none',internal_energy=held,internal_energy_flux=scalar, &
        energy_policy=policy,energy_correction=correction,auxiliary_source=serial_source,auxiliary_change=aux, &
        auxiliary_has_internal_energy=.true.,integrated_flux=faces,held_interval_history=records,level_steps=steps, &
        source_funding_context=context,native_common_retries=retries,evaluate_reaction=refusal_reactor)
    else if(mode==1) then
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,duration,2,.true., &
        transport_flags,force,2.e-7_dp,1.e-12_dp,gas,extended,boundaries,b,w,ch,ext,theta,accepted, &
        reaction_source=ir,mol_mode=.true.,redistribution='none',internal_energy=held,internal_energy_flux=scalar, &
        energy_policy=policy,energy_correction=correction,auxiliary_source=serial_source,auxiliary_change=aux, &
        auxiliary_has_internal_energy=.true.,integrated_flux=faces,held_interval_history=records,level_steps=steps, &
        source_funding_context=context,native_common_retries=retries,native_replay_history=replay,native_spatial_weight=.5_dp)
    else
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,duration,2,.true., &
        transport_flags,force,2.e-7_dp,1.e-12_dp,gas,extended,boundaries,b,w,ch,ext,theta,accepted, &
        reaction_source=ir,mol_mode=.true.,redistribution='none',internal_energy=held,internal_energy_flux=scalar, &
        energy_policy=policy,energy_correction=correction,auxiliary_source=serial_source,auxiliary_change=aux, &
        auxiliary_has_internal_energy=.true.,integrated_flux=faces,held_interval_history=records,level_steps=steps, &
        source_funding_context=context,native_common_retries=retries)
    endif
    if(mode>=2) then
      call check(.not.accepted.and.source_calls==2.and.old_calls==1.and.retries<=1, &
        'serial caller refusal stops without another source evaluation')
      call check(all(gas%nodes(1)%state==initial%nodes(1)%state).and. &
        all(gas%nodes(1)%temperature==initial%nodes(1)%temperature).and. &
        all(held%fields(1)%values==initial_h%fields(1)%values).and.all(ir%fields(1)%values==0), &
        'serial caller refusal rolls back gas/T/H/I_R')
      call check(all([b,w,ch,ext,aux,theta,correction]== &
        [spread(17._dp,1,nv),spread(19._dp,1,nv),spread(23._dp,1,nv),spread(29._dp,1,nv), &
        spread(31._dp,1,naux),.7_dp,0._dp]).and.all(steps==-17_int64),'serial refusal rolls back physical ledgers')
      call check(same_faces(faces(1),reference_faces(1)).and.same_faces(scalar(1),reference_scalar(1)), &
        'serial refusal retains previous accepted faces')
      call check(size(records)==1.and.all(records(1)%local_tree%nodes(1)%state==reference%nodes(1)%state).and. &
        all(records(1)%auxiliary_stages(0)%fields(1)%values==reference_records(1)%auxiliary_stages(0)%fields(1)%values), &
        'serial refusal retains nonempty accepted history')
      call check_settled(context,[.true.])
      if(mode==3) call check(reactor_calls==1,'caller reactor refused a native-admitted input once')
      return
    endif
    call check(accepted.and.theta==.5_dp.and.all(steps==1_int64),'serial half-weight recovered physical clock')
    call check(retries==merge(1_int64,0_int64,mode==0),'automatic recovery differs from complete explicit replay')
    call check(source_calls==merge(2,0,mode==0).and.old_calls==merge(1,0,mode==0), &
      'source packets are observed once and frozen replay calls no source')
    call check(size(records)==1.and.records(1)%replayed,'accepted partial recovery is recorded as replayed')
    call check_settled(context,[.true.])
    call check_inventory(gas,held,records,faces,scalar,b,w,ch,ext,aux,correction)
    if(mode==0) then
      reference=gas;reference_h=held;reference_ir=ir;reference_records=records;reference_context=context
      reference_faces=faces;reference_scalar=scalar;reference_ledgers=[b,w,ch,ext,aux,theta,correction]
    else
      call check(all(gas%nodes(1)%state==reference%nodes(1)%state).and. &
        all(gas%nodes(1)%temperature==reference%nodes(1)%temperature).and. &
        all(held%fields(1)%values==reference_h%fields(1)%values).and. &
        all(ir%fields(1)%values==reference_ir%fields(1)%values),'complete frozen replay matches automatic gas/T/H/I_R')
      call check(all([b,w,ch,ext,aux,theta,correction]==reference_ledgers).and. &
        same_faces(faces(1),reference_faces(1)).and.same_faces(scalar(1),reference_scalar(1)), &
        'complete frozen replay matches automatic faces and ledgers')
    endif
  end subroutine check_serial

  subroutine check_mpi_layout()
    type(amr_tree_state_3d) :: gas,temp,held,ir
    type(eb_held_funding_context_3d) :: context,empty_context
    type(mpi_eb_mol_interval_3d), allocatable :: records(:),saved(:),replay(:)
    type(eb_integrated_flux_3d), allocatable :: faces(:),scalar(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(naux),theta,correction
    integer(int64) :: retries,steps(1)
    integer :: mode
    logical :: accepted,local
    local=owner(1)==rank
    do mode=0,3
      context=empty_context
      allocate(gas%fields(1),temp%fields(1),held%fields(1))
      if(local) then
        gas%fields(1)%values=initial%nodes(1)%state;held%fields(1)%values=initial_h%fields(1)%values
        allocate(temp%fields(1)%values(1,4,4,4));temp%fields(1)%values(1,:,:,:)=initial%nodes(1)%temperature
      endif
      call zero_ir(ir,local)
      b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp;correction=0;steps=-17
      if(mode==1) replay=saved
      if(mode>=2) then
        records=saved;context=reference_context
        if(.not.local) then
          if(allocated(context%pending(0)%fields(1)%values)) deallocate(context%pending(0)%fields(1)%values)
        endif
        allocate(faces(1),scalar(1))
        if(local) then
          faces(1)=reference_faces(1);scalar(1)=reference_scalar(1)
        endif
      endif
      call reset_callbacks(mode)
      if(mode==3) then
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owner,geometry,gas,temp,duration, &
          2,.true.,transport_flags,force,2.e-7_dp,1.e-12_dp,extended,boundaries,b,w,ch,ext,theta,MPI_COMM_WORLD,accepted, &
          reaction_source=ir,mol_mode=.true.,redistribution='none',internal_energy=held,internal_energy_flux=scalar, &
          energy_policy=policy,energy_correction=correction,auxiliary_source=partial_mpi_source,auxiliary_change=aux, &
          auxiliary_has_internal_energy=.true.,integrated_flux=faces,interval_history=records,level_steps=steps, &
          source_funding_context=context,native_common_retries=retries,evaluate_reaction=refusal_reactor)
      else if(mode==1) then
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owner,geometry,gas,temp,duration, &
          2,.true.,transport_flags,force,2.e-7_dp,1.e-12_dp,extended,boundaries,b,w,ch,ext,theta,MPI_COMM_WORLD,accepted, &
          reaction_source=ir,mol_mode=.true.,redistribution='none',internal_energy=held,internal_energy_flux=scalar, &
          energy_policy=policy,energy_correction=correction,auxiliary_source=partial_mpi_source,auxiliary_change=aux, &
          auxiliary_has_internal_energy=.true.,integrated_flux=faces,interval_history=records,level_steps=steps, &
          source_funding_context=context,native_common_retries=retries,native_replay_history=replay,native_spatial_weight=.5_dp)
      else
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owner,geometry,gas,temp,duration, &
          2,.true.,transport_flags,force,2.e-7_dp,1.e-12_dp,extended,boundaries,b,w,ch,ext,theta,MPI_COMM_WORLD,accepted, &
          reaction_source=ir,mol_mode=.true.,redistribution='none',internal_energy=held,internal_energy_flux=scalar, &
          energy_policy=policy,energy_correction=correction,auxiliary_source=partial_mpi_source,auxiliary_change=aux, &
          auxiliary_has_internal_energy=.true.,integrated_flux=faces,interval_history=records,level_steps=steps, &
          source_funding_context=context,native_common_retries=retries)
      endif
      if(mode>=2) then
        call check(.not.accepted.and.source_calls==2.and.old_calls==1.and.retries<=1, &
          'MPI owner/idle caller refusal does not repeat a source')
        call check(all([b,w,ch,ext,aux,theta,correction]== &
          [spread(17._dp,1,nv),spread(19._dp,1,nv),spread(23._dp,1,nv),spread(29._dp,1,nv), &
          spread(31._dp,1,naux),.7_dp,0._dp]).and.all(steps==-17_int64),'MPI refusal rolls back physical ledgers')
        if(local) then
          call check(all(gas%fields(1)%values==initial%nodes(1)%state).and. &
            all(temp%fields(1)%values(1,:,:,:)==initial%nodes(1)%temperature).and. &
            all(held%fields(1)%values==initial_h%fields(1)%values).and.all(ir%fields(1)%values==0), &
            'MPI refusal rolls back caller gas/T/H/I_R')
          call check(same_faces(faces(1),reference_faces(1)).and.same_faces(scalar(1),reference_scalar(1)), &
            'MPI refusal retains previous accepted faces')
          call check(size(records)==1.and.all(records(1)%local_state%fields(1)%values== &
            saved(1)%local_state%fields(1)%values).and.all(records(1)%auxiliary_stages(0)%fields(1)%values== &
            saved(1)%auxiliary_stages(0)%fields(1)%values),'MPI refusal retains nonempty accepted history')
        endif
        call check_settled(context,owner==rank)
        if(mode==3) call check(reactor_calls==merge(1,0,local),'only the owner calls the native-admitted reactor')
      else
        call check(accepted.and.theta==.5_dp.and.all(steps==1_int64),'MPI half-weight recovered physical clock')
        call check(retries==merge(1_int64,0_int64,mode==0).and.source_calls==merge(2,0,mode==0).and. &
          old_calls==merge(1,0,mode==0),'MPI automatic and frozen replay source/retry schedule')
        call check(size(records)==1.and.records(1)%replayed,'MPI accepted partial recovery is replayed')
        call check(all([b,w,ch,ext,aux,theta,correction]==reference_ledgers),'serial/MPI accepted ledgers match exactly')
        if(local) then
          call check(all(gas%fields(1)%values==reference%nodes(1)%state).and. &
            all(temp%fields(1)%values(1,:,:,:)==reference%nodes(1)%temperature).and. &
            all(held%fields(1)%values==reference_h%fields(1)%values).and. &
            all(ir%fields(1)%values==reference_ir%fields(1)%values),'serial/MPI recovered gas/T/H/I_R match exactly')
          call check(same_faces(faces(1),reference_faces(1)).and.same_faces(scalar(1),reference_scalar(1)), &
            'serial/MPI recovered gas/scalar faces match exactly')
          call check(all(records(1)%source_funding_impulse%fields(1)%values==0), &
            'single root recovery invents no off-level material')
        endif
        call check_settled(context,owner==rank)
        if(mode==0) saved=records
      endif
      if(allocated(records)) deallocate(records)
      if(allocated(replay)) deallocate(replay)
      if(allocated(faces)) deallocate(faces)
      if(allocated(scalar)) deallocate(scalar)
      deallocate(gas%fields,temp%fields,held%fields,ir%fields)
    enddo
  end subroutine check_mpi_layout

  subroutine check_settled(context,local)
    type(eb_held_funding_context_3d), intent(in) :: context
    logical, intent(in) :: local(:)
    call check(valid_eb_held_funding_context_3d(topology,geometry,nv,local,context,0), &
      'funding context has a complete settled native clock')
    call check(context%open(0)%clock==-1,'root recovery leaves no open funding family')
    if(local(1)) call check(all(context%pending(0)%fields(1)%values==0),'root recovery leaves no pending receipt')
  end subroutine check_settled

  subroutine check_inventory(gas,held,records,faces,scalar,b,w,ch,ext,aux,correction)
    type(eb_reactive_tree_3d), intent(in) :: gas
    type(amr_tree_state_3d), intent(in) :: held
    type(eb_mol_held_interval_3d), intent(in) :: records(:)
    type(eb_integrated_flux_3d), intent(in) :: faces(:),scalar(:)
    real(dp), intent(in) :: b(:),w(:),ch(:),ext(:),aux(:),correction
    real(dp) :: before(nv),after(nv),residual(nv),face_change(nv),expected(naux),hbefore,hafter,hfaces,hvolume
    real(dp) :: volume,impulse(1),gas_error,h_error
    integer :: x,y,z,stage
    logical :: accepted
    call integrate_reactive_eb_tree_3d(initial,before,accepted);call check(accepted,'entry gas inventory')
    call integrate_reactive_eb_tree_3d(gas,after,accepted);call check(accepted,'accepted gas inventory')
    face_change=0;expected=0;hbefore=0;hafter=0;hfaces=0;hvolume=0
    associate(g=>geometry(1))
    do z=1,4
    do y=1,4
    do x=1,4
      if(g%volume_fraction(x,y,z)<=0) cycle
      volume=g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz
      face_change=face_change+face_impulse(g,faces(1),x,y,z)
      impulse=face_impulse(g,scalar(1),x,y,z);hfaces=hfaces+impulse(1)
      hbefore=hbefore+volume*initial_h%fields(1)%values(1,x,y,z)
      hafter=hafter+volume*held%fields(1)%values(1,x,y,z)
      hvolume=hvolume+duration*volume*records(1)%internal_energy_spatial_rate%fields(1)%values(1,x,y,z)-impulse(1)
      if(x==4) then
        expected(:nv)=expected(:nv)-.5_dp*sink*volume*initial%nodes(1)%state(:,x,y,z)
        expected(naux)=expected(naux)-.5_dp*sink*volume*initial_h%fields(1)%values(1,x,y,z)
      endif
      call check(gas%nodes(1)%state(1,x,y,z)>0.and.all(gas%nodes(1)%state(6:,x,y,z)>=0).and. &
        held%fields(1)%values(1,x,y,z)>0,'accepted gas/species/H are admissible')
    enddo
    enddo
    enddo
    end associate
    call check(all(abs(aux-expected)<=32*epsilon(1._dp)*max(1._dp,abs(expected))), &
      'accepted gas/H auxiliary impulse is exactly one trapezoidal source')
    residual=after-before-face_change-w-ch-ext-aux(:nv);residual(5)=residual(5)-correction
    gas_error=maxval(abs(residual)/max(1._dp,abs(before)))
    h_error=abs(hafter-hbefore-hfaces-hvolume-records(1)%held_clock_projection_after- &
      records(1)%held_clock_reset_after)/max(1._dp,abs(hbefore))
    call check(gas_error<3.e-14_dp.and.h_error<3.e-14_dp,'independent gas/species and actual H inventories close')
    call check(all(abs(b-face_change)<=32*epsilon(1._dp)*max(1._dp,abs(face_change))), &
      'reported boundary inventory matches independently summed physical faces')
    do stage=0,refused_prefix%accepted_source_stages-1
      call check(all(records(1)%external_stages(stage)%fields(1)%values== &
        refused_prefix%external_stages(stage)%fields(1)%values).and. &
        all(records(1)%auxiliary_stages(stage)%fields(1)%values==refused_prefix%auxiliary_stages(stage)%fields(1)%values), &
        'accepted recovery retains every source observed before the native failure')
    enddo
    call check(all(records(1)%source_funding_impulse%fields(1)%values==0),'no unobserved funding credit is spent')
    if(rank==0) print *, 'PARTIAL_ROOT_INVENTORY',gas_error,h_error,aux(1),aux(naux)
  end subroutine check_inventory

  function face_impulse(g,f,x,y,z) result(value)
    type(eb_geometry_3d), intent(in) :: g
    type(eb_integrated_flux_3d), intent(in) :: f
    integer, intent(in) :: x,y,z
    real(dp) :: value(size(f%x,1))
    value=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(:,x-1,y,z)-g%x_face_fraction(x,y,z)*f%x(:,x,y,z))+ &
      g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(:,x,y-1,z)-g%y_face_fraction(x,y,z)*f%y(:,x,y,z))+ &
      g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(:,x,y,z-1)-g%z_face_fraction(x,y,z)*f%z(:,x,y,z))
  end function face_impulse

  logical function same_faces(a,b) result(same)
    type(eb_integrated_flux_3d), intent(in) :: a,b
    same=.false.
    if(.not.allocated(a%x).or..not.allocated(b%x).or..not.allocated(a%y).or..not.allocated(b%y).or. &
      .not.allocated(a%z).or..not.allocated(b%z)) return
    if(any(shape(a%x)/=shape(b%x)).or.any(shape(a%y)/=shape(b%y)).or.any(shape(a%z)/=shape(b%z))) return
    same=all(a%x==b%x).and.all(a%y==b%y).and.all(a%z==b%z)
  end function same_faces

  subroutine fill_source(state,interval,stage,rate)
    type(amr_tree_state_3d), intent(in) :: state
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage
    type(amr_tree_state_3d), intent(inout) :: rate
    call check(callback_mode/=1,'complete frozen replay never calls a source')
    source_calls=source_calls+1
    if(stage==0) then
      old_calls=old_calls+1;call check(old_calls==1,'observed old source is never requested twice')
    endif
    if(allocated(rate%fields(1)%values)) then
      rate%fields(1)%values=0
      if(stage==0) rate%fields(1)%values(:,4,:,:)=-sink/interval*state%fields(1)%values(:,4,:,:)
    endif
  end subroutine fill_source

  subroutine serial_source(state,interval,stage,depth,time,rate,accepted)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: accepted
    type(amr_tree_state_3d) :: packed
    allocate(packed%fields(1));packed%fields(1)%values=state%nodes(1)%state
    call fill_source(packed,interval,stage,rate)
    accepted=.not.(callback_mode==2.and.stage==1)
  end subroutine serial_source

  subroutine common_source(state,interval,stage,depth,rate,accepted)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: accepted
    call serial_source(state,interval,stage,depth,0._dp,rate,accepted)
  end subroutine common_source

  subroutine partial_mpi_source(state,temp,interval,stage,depth,time,rate,comm,accepted)
    type(amr_tree_state_3d), intent(in) :: state,temp
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: accepted
    call fill_source(state,interval,stage,rate)
    accepted=.not.(callback_mode==2.and.stage==1.and.rank==1)
  end subroutine partial_mpi_source

  subroutine refusal_reactor(species,reactions,energy,forcing,energy_forcing,interval,rtol,atol,density,temp, &
      completed_steps,accepted,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: energy,forcing(:),energy_forcing,interval,rtol,atol
    real(dp), intent(inout) :: density(:),temp
    integer, intent(out) :: completed_steps
    logical, intent(out) :: accepted
    integer, optional, intent(in) :: maximum_steps
    reactor_calls=reactor_calls+1
    call advance_forced_reactor_interval(species,reactions,energy,forcing,energy_forcing,interval,rtol,atol, &
      density,temp,completed_steps,accepted,maximum_steps=maximum_steps)
    call check(accepted,'caller refusal follows an independent native reactor admission')
    accepted=.false.;density=-1;temp=-1
  end subroutine refusal_reactor

  subroutine check(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if(condition) return
    print *, 'FAIL:',rank,trim(label)
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine check
end program test_mpi_eb_native_mol_partial_recovery_3d
