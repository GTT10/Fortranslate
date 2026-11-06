program test_mpi_eb_native_mol_partial_fine_3d
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
  use dual_energy_mod, only: dual_energy_policy,initialize_nasa7_internal_energy_patch
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    build_eb_tree_extended_geometry_3d,advance_reactive_eb_tree_sdc_3d,advance_reactive_eb_tree_sdc_subcycle_3d, &
    integrate_reactive_eb_tree_3d,eb_mol_held_interval_3d,eb_native_mol_source_prefix_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_held_funding_queue_3d_mod, only: eb_held_funding_context_3d,valid_eb_held_funding_context_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_subcycle_3d,mpi_eb_mol_interval_3d
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use amr_patch_tree_remap_3d_mod, only: covered
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_tree_3d) :: topology
  type(amr_patch_3d) :: patches(1)
  type(eb_geometry_3d) :: fine,padded,geometry(2)
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
  integer :: rank,nranks,status,nv,naux,nitrogen,p,i,j,k,owner(2),layout,mode
  integer :: source_calls=0,old_calls=0,reactor_calls=0,callback_mode=0
  logical :: good,refusal_armed=.false.
  real(dp) :: last_old_time=-1
  real(dp) :: density,rg,temperature,speed
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call check(nranks==3,'distributed owners and an idle rank')
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
  call coarsen_eb_geometry_3d(fine,2,geometry(1),good);call check(good,'root geometry')
  call initialize_amr_patch_3d(4,4,4,3,4,1,4,1,4,2,patches(1),good);call check(good,'child patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1],patches,topology,good);call check(good,'two-level topology')
  allocate(u(nv,4,4,4),temperature_field(4,4,4),source=0._dp)
  do k=1,4
  do j=1,4
  do i=1,4
    if(geometry(1)%volume_fraction(i,j,k)<=0) cycle
    ! Exact coarse volume averages of the geometric fine profile.
    select case(i)
    case(2);density=8
    case(3);density=3
    case(4);density=.75_dp
    end select
    primitive(1)=density;primitive(2)=1000;primitive(3:4)=0;primitive(5)=density*rg*1400
    call reactive_primitive_to_conserved(species,primitive,cell,temperature,speed,good)
    call check(good,'physical density-gradient gas')
    u(:,i,j,k)=cell;temperature_field(i,j,k)=temperature
  enddo
  enddo
  enddo
  call initialize_reactive_eb_tree_3d(species,topology,fine,u,temperature_field,initial,good,internal_energy=initial_h)
  call check(good,'physical two-level gas and independent H')
  do p=1,2
    geometry(p)=initial%nodes(p)%geometry
  enddo
  ! Fine x densities are 4,2,1,1/2; covered coarse cells retain their means.
  do k=1,geometry(2)%nz
  do j=1,geometry(2)%ny
  do i=1,geometry(2)%nx
    density=2._dp**(3-i)
    primitive(1)=density;primitive(2)=1000;primitive(3:4)=0;primitive(5)=density*rg*1400
    call reactive_primitive_to_conserved(species,primitive,cell,temperature,speed,good)
    call check(good,'physical fine density profile')
    initial%nodes(2)%state(:,i,j,k)=cell;initial%nodes(2)%temperature(i,j,k)=temperature
  enddo
  enddo
  enddo
  call initialize_nasa7_internal_energy_patch(species,initial%nodes(2)%state,initial%nodes(2)%temperature, &
    geometry(2)%volume_fraction,initial_h%fields(2)%values(1,:,:,:),good)
  call check(good,'independent fine thermal H')
  call build_level_set_eb_geometry_3d(padded_phi,spread(-.006_dp,1,3),spread(.010_dp,1,3),padded,good)
  call check(good,'padded cut geometry')
  call build_eb_tree_extended_geometry_3d(topology,geometry,padded,extended,good,halo_width=6)
  call check(good,'six-cell native spatial context')
  policy%allow_negative_energy=.false.;policy%eta2=2;policy%allow_small_energy=.false.
  call check_serial(0)
  call check_native_refusal()
  do mode=1,2
    call check_serial(mode)
  enddo
  do layout=1,2
    owner=[0,2];if(layout==2) owner=0
    call check_mpi_layout()
    call MPI_Barrier(MPI_COMM_WORLD,status);call check(status==MPI_SUCCESS,'owner layout complete')
    if(rank==0) print *, 'PASS_NATIVE_MOL_PARTIAL_FINE',owner
  enddo
  call MPI_Finalize(status)
contains
  subroutine reset_callbacks(mode)
    integer, intent(in) :: mode
    callback_mode=mode;source_calls=0;old_calls=0;reactor_calls=0;refusal_armed=.false.;last_old_time=-1
  end subroutine reset_callbacks

  subroutine zero_ir(ir,local)
    type(amr_tree_state_3d), intent(out) :: ir
    logical, intent(in) :: local(:)
    integer :: node,d(3)
    allocate(ir%fields(2))
    do node=1,2
      d=topology%nodes(node)%shape
      if(local(node)) allocate(ir%fields(node)%values(nv-4,d(1),d(2),d(3)),source=0._dp)
    enddo
  end subroutine zero_ir

  subroutine make_faces(faces,components)
    type(eb_integrated_flux_3d), allocatable, intent(out) :: faces(:)
    integer, intent(in) :: components
    integer :: node,d(3)
    allocate(faces(2))
    do node=1,2
      d=topology%nodes(node)%shape
      allocate(faces(node)%x(components,0:d(1),d(2),d(3)),faces(node)%y(components,d(1),0:d(2),d(3)), &
        faces(node)%z(components,d(1),d(2),0:d(3)),source=0._dp)
    enddo
  end subroutine make_faces

  subroutine check_native_refusal()
    type(eb_reactive_tree_3d) :: gas,entry
    type(amr_tree_state_3d) :: held,ir,stock,receipt,entry_h,entry_ir
    type(eb_integrated_flux_3d), allocatable :: faces(:),scalar(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(naux),theta,correction
    logical :: accepted,ready
    integer :: node,stage
    associate(anchor=>reference_records(2))
    gas=anchor%initial;held=anchor%initial_internal_energy;ir=anchor%initial_reaction_source
    entry=gas;entry_h=held;entry_ir=ir
    stock=anchor%source_funding_stock;receipt=stock
    do node=1,2
      receipt%fields(node)%values=0
    enddo
    call make_faces(faces,nv);call make_faces(scalar,1)
    b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp;correction=0
    call reset_callbacks(0)
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,anchor%interval,2,.true.,transport_flags,force, &
      2.e-7_dp,1.e-12_dp,gas,extended,boundaries,b,w,ch,ext,theta,accepted,reaction_source=ir, &
      target_depth=1,ancestor_start=anchor%ancestor_start,ancestor_end=anchor%ancestor_end, &
      fractions_start=anchor%fractions_start,fractions_end=anchor%fractions_end, &
      source_start=anchor%ancestor_source_start,source_end=anchor%ancestor_source_end, &
      ancestor_internal_energy_start=anchor%ancestor_internal_energy_start, &
      ancestor_internal_energy_end=anchor%ancestor_internal_energy_end, &
      mol_mode=.true.,integrated_flux=faces,internal_energy=held,internal_energy_flux=scalar, &
      energy_policy=policy,energy_correction=correction,auxiliary_source=common_source,auxiliary_change=aux, &
      auxiliary_has_internal_energy=.true.,source_funding_stock=stock,source_funding_impulse=receipt, &
      native_source_prefix=refused_prefix,native_source_prefix_ready=ready)
    call check(.not.accepted.and.ready,'unit-weight first fine clock refuses under its actual open ancestor')
    call check(refused_prefix%accepted_source_stages==2.and.old_calls==1.and.source_calls==2, &
      'original fine clock observes two of three source stages')
    stage=2
    call check(.not.allocated(refused_prefix%external_stages(stage)%fields).and. &
      .not.allocated(refused_prefix%auxiliary_stages(stage)%fields),'unobserved fine future packet stays absent')
    do node=1,2
      call check(all(gas%nodes(node)%state==entry%nodes(node)%state).and. &
        all(gas%nodes(node)%temperature==entry%nodes(node)%temperature).and. &
        all(held%fields(node)%values==entry_h%fields(node)%values).and. &
        all(ir%fields(node)%values==entry_ir%fields(node)%values),'refused fine clock preserves caller and ancestor fields')
      call check(all(faces(node)%x==0).and.all(faces(node)%y==0).and.all(faces(node)%z==0).and. &
        all(scalar(node)%x==0).and.all(scalar(node)%y==0).and.all(scalar(node)%z==0).and. &
        all(receipt%fields(node)%values==0),'refused fine clock preserves paired faces and receipts')
    enddo
    call check(all([b,w,ch,ext,aux,theta,correction]==sentinels()),'refused fine clock preserves numerical ledgers')
    do stage=0,1
      call check(all(anchor%external_stages(stage)%fields(2)%values== &
        refused_prefix%external_stages(stage)%fields(2)%values).and. &
        all(anchor%auxiliary_stages(stage)%fields(2)%values==refused_prefix%auxiliary_stages(stage)%fields(2)%values), &
        'accepted fine history retains both previously observed packets')
    enddo
    ! Consume the independently captured partial prefix under the same ancestor.
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,anchor%interval,2,.true.,transport_flags,force, &
      2.e-7_dp,1.e-12_dp,gas,extended,boundaries,b,w,ch,ext,theta,accepted,reaction_source=ir, &
      target_depth=1,ancestor_start=anchor%ancestor_start,ancestor_end=anchor%ancestor_end, &
      fractions_start=anchor%fractions_start,fractions_end=anchor%fractions_end, &
      source_start=anchor%ancestor_source_start,source_end=anchor%ancestor_source_end, &
      ancestor_internal_energy_start=anchor%ancestor_internal_energy_start, &
      ancestor_internal_energy_end=anchor%ancestor_internal_energy_end, &
      mol_mode=.true.,integrated_flux=faces,internal_energy=held,internal_energy_flux=scalar, &
      energy_policy=policy,energy_correction=correction,auxiliary_source=common_source,auxiliary_change=aux, &
      auxiliary_has_internal_energy=.true.,source_funding_stock=stock,source_funding_impulse=receipt, &
      native_source_replay=refused_prefix,native_spatial_weight=.5_dp)
    call check(accepted.and.source_calls==2.and.old_calls==1,'explicit fine prefix consumption repeats no source')
    call check(all(gas%nodes(2)%state==anchor%local_tree%nodes(2)%state).and. &
      all(gas%nodes(2)%temperature==anchor%local_tree%nodes(2)%temperature).and. &
      all(held%fields(2)%values==anchor%local_internal_energy%fields(2)%values).and. &
      all(ir%fields(2)%values==anchor%reaction_source%fields(2)%values), &
      'explicit fine prefix consumption matches the automatically accepted common clock')
    call check(same_faces(faces(2),anchor%spatial_flux(2)).and.same_faces(scalar(2),anchor%internal_energy_flux(2)), &
      'explicit fine prefix reconstructs automatic gas and actual-H faces')
    end associate
    if(rank==0) print *, 'PARTIAL_FINE_INITIAL_OBSERVED',refused_prefix%accepted_source_stages
  end subroutine check_native_refusal

  function sentinels() result(values)
    real(dp) :: values(5*nv+3)
    values=[spread(17._dp,1,nv),spread(19._dp,1,nv),spread(23._dp,1,nv),spread(29._dp,1,nv), &
      spread(31._dp,1,naux),.7_dp,0._dp]
  end function sentinels

  subroutine check_serial(mode)
    integer, intent(in) :: mode
    type(eb_reactive_tree_3d) :: gas
    type(amr_tree_state_3d) :: held,ir
    type(eb_held_funding_context_3d) :: context
    type(eb_mol_held_interval_3d), allocatable :: records(:)
    type(eb_integrated_flux_3d), allocatable :: faces(:),scalar(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(naux),theta,correction
    integer(int64) :: retries,steps(2)
    integer :: node,r
    logical :: accepted
    gas=initial;held=initial_h;call zero_ir(ir,[.true.,.true.])
    b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp;correction=0;steps=-17
    if(mode>0) then
      records=reference_records;faces=reference_faces;scalar=reference_scalar;context=reference_context
    endif
    call reset_callbacks(mode)
    if(mode==2) then
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,duration,2,.true., &
        transport_flags,force,2.e-7_dp,1.e-12_dp,gas,extended,boundaries,b,w,ch,ext,theta,accepted, &
        reaction_source=ir,mol_mode=.true.,redistribution='none',internal_energy=held,internal_energy_flux=scalar, &
        energy_policy=policy,energy_correction=correction,auxiliary_source=serial_source,auxiliary_change=aux, &
        auxiliary_has_internal_energy=.true.,integrated_flux=faces,held_interval_history=records,level_steps=steps, &
        source_funding_context=context,native_common_retries=retries,evaluate_reaction=refusal_reactor)
    else
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,duration,2,.true., &
        transport_flags,force,2.e-7_dp,1.e-12_dp,gas,extended,boundaries,b,w,ch,ext,theta,accepted, &
        reaction_source=ir,mol_mode=.true.,redistribution='none',internal_energy=held,internal_energy_flux=scalar, &
        energy_policy=policy,energy_correction=correction,auxiliary_source=serial_source,auxiliary_change=aux, &
        auxiliary_has_internal_energy=.true.,integrated_flux=faces,held_interval_history=records,level_steps=steps, &
        source_funding_context=context,native_common_retries=retries)
    endif
    if(mode>0) then
      if(rank==0) print *, 'PARTIAL_FINE_SERIAL_REFUSAL',mode,accepted,retries,source_calls,reactor_calls
      call check(.not.accepted.and.source_calls==6.and.old_calls==2.and.retries==1, &
        'late serial child refusal stops after fine recovery and one completed child')
      call check(all([b,w,ch,ext,aux,theta,correction]==sentinels()).and.all(steps==-17_int64), &
        'late serial refusal preserves physical and numerical ledgers')
      call check(size(records)==3,'late serial refusal retains nonempty accepted history')
      do node=1,2
        call check(all(gas%nodes(node)%state==initial%nodes(node)%state).and. &
          all(gas%nodes(node)%temperature==initial%nodes(node)%temperature).and. &
          all(held%fields(node)%values==initial_h%fields(node)%values).and.all(ir%fields(node)%values==0), &
          'late serial refusal rolls back every caller gas/T/H/I_R')
        call check(same_faces(faces(node),reference_faces(node)).and.same_faces(scalar(node),reference_scalar(node)), &
          'late serial refusal retains accepted gas and scalar faces')
        do r=1,3
          call check(all(records(r)%local_tree%nodes(node)%state==reference_records(r)%local_tree%nodes(node)%state), &
            'late serial refusal retains all accepted clock states')
        enddo
      enddo
      if(mode==2) call check(reactor_calls==1,'late native-admitted serial reactor refusal occurs once')
    else
      call check(accepted.and.theta==.5_dp.and.all(steps==[1_int64,2_int64]).and.retries==1, &
        'serial partial fine recovery completes the two child clocks')
      call check(source_calls==6.and.old_calls==2,'each serial source packet is evaluated once')
      call check(size(records)==3.and.(.not.records(1)%replayed).and. &
        records(2)%replayed.and.(.not.records(3)%replayed), &
        'accepted history replaces only the recovered first common fine clock')
      call check_inventory(gas,held,records,faces,scalar,b,w,ch,ext,aux,correction)
      reference=gas;reference_h=held;reference_ir=ir;reference_records=records;reference_context=context
      reference_faces=faces;reference_scalar=scalar;reference_ledgers=[b,w,ch,ext,aux,theta,correction]
    endif
    call check_settled(context,[.true.,.true.])
  end subroutine check_serial

  subroutine check_mpi_layout()
    type(amr_tree_state_3d) :: gas,temp,held,ir
    type(eb_held_funding_context_3d) :: context,empty_context
    type(mpi_eb_mol_interval_3d), allocatable :: records(:),saved(:)
    type(eb_integrated_flux_3d), allocatable :: faces(:),scalar(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(naux),theta,correction
    integer(int64) :: retries,steps(2)
    integer :: mode,node,r,bucket,d(3)
    logical :: accepted,local(2)
    local=owner==rank
    do mode=0,2
      context=empty_context
      allocate(gas%fields(2),temp%fields(2),held%fields(2))
      do node=1,2
        if(.not.local(node)) cycle
        d=topology%nodes(node)%shape
        gas%fields(node)%values=initial%nodes(node)%state;held%fields(node)%values=initial_h%fields(node)%values
        allocate(temp%fields(node)%values(1,d(1),d(2),d(3)))
        temp%fields(node)%values(1,:,:,:)=initial%nodes(node)%temperature
      enddo
      call zero_ir(ir,local)
      b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp;correction=0;steps=-17
      if(mode>0) then
        records=saved;context=reference_context;allocate(faces(2),scalar(2))
        do node=1,2
          if(local(node)) then
            faces(node)=reference_faces(node);scalar(node)=reference_scalar(node)
          else
            do bucket=0,1
              deallocate(context%pending(bucket)%fields(node)%values)
            enddo
          endif
        enddo
      endif
      call reset_callbacks(mode)
      if(mode==2) then
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owner,geometry,gas,temp,duration, &
          2,.true.,transport_flags,force,2.e-7_dp,1.e-12_dp,extended,boundaries,b,w,ch,ext,theta,MPI_COMM_WORLD,accepted, &
          reaction_source=ir,mol_mode=.true.,redistribution='none',internal_energy=held,internal_energy_flux=scalar, &
          energy_policy=policy,energy_correction=correction,auxiliary_source=partial_mpi_source,auxiliary_change=aux, &
          auxiliary_has_internal_energy=.true.,integrated_flux=faces,interval_history=records,level_steps=steps, &
          source_funding_context=context,native_common_retries=retries,evaluate_reaction=refusal_reactor)
      else
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owner,geometry,gas,temp,duration, &
          2,.true.,transport_flags,force,2.e-7_dp,1.e-12_dp,extended,boundaries,b,w,ch,ext,theta,MPI_COMM_WORLD,accepted, &
          reaction_source=ir,mol_mode=.true.,redistribution='none',internal_energy=held,internal_energy_flux=scalar, &
          energy_policy=policy,energy_correction=correction,auxiliary_source=partial_mpi_source,auxiliary_change=aux, &
          auxiliary_has_internal_energy=.true.,integrated_flux=faces,interval_history=records,level_steps=steps, &
          source_funding_context=context,native_common_retries=retries)
      endif
      if(mode>0) then
        if(rank==0) print *, 'PARTIAL_FINE_MPI_REFUSAL',owner,mode,accepted,retries,source_calls
        call check(.not.accepted.and.source_calls==6.and.old_calls==2.and.retries==1, &
          'late MPI owner/idle refusal follows fine recovery and the first child')
        call check(all([b,w,ch,ext,aux,theta,correction]==sentinels()).and.all(steps==-17_int64), &
          'late MPI refusal preserves caller physical and numerical ledgers')
        call check(size(records)==3,'late MPI refusal retains nonempty accepted history')
        if(mode==2) call check(reactor_calls==merge(1,0,local(2)),'only the final child owner refuses a reactor')
      else
        call check(accepted.and.theta==.5_dp.and.all(steps==[1_int64,2_int64]).and.retries==1, &
          'MPI partial fine recovery completes two child clocks')
        call check(source_calls==6.and.old_calls==2,'MPI observed source packets are never repeated')
        call check(size(records)==3.and.(.not.records(1)%replayed).and. &
        records(2)%replayed.and.(.not.records(3)%replayed), &
          'MPI history replaces only the recovered first fine clock')
        call check(all([b,w,ch,ext,aux,theta,correction]==reference_ledgers),'serial/MPI physical ledgers are exact')
      endif
      do node=1,2
        if(.not.local(node)) then
          call check(.not.allocated(gas%fields(node)%values).and..not.allocated(held%fields(node)%values), &
            'idle MPI fields remain unallocated')
          cycle
        endif
        if(mode>0) then
          call check(all(gas%fields(node)%values==initial%nodes(node)%state).and. &
            all(temp%fields(node)%values(1,:,:,:)==initial%nodes(node)%temperature).and. &
            all(held%fields(node)%values==initial_h%fields(node)%values).and.all(ir%fields(node)%values==0), &
            'late MPI refusal rolls back both root and child gas/T/H/I_R')
        else
          call check(all(gas%fields(node)%values==reference%nodes(node)%state).and. &
            all(temp%fields(node)%values(1,:,:,:)==reference%nodes(node)%temperature).and. &
            all(held%fields(node)%values==reference_h%fields(node)%values).and. &
            all(ir%fields(node)%values==reference_ir%fields(node)%values),'serial/MPI physical fields are exact')
        endif
        call check(same_faces(faces(node),reference_faces(node)).and.same_faces(scalar(node),reference_scalar(node)), &
          'accepted or retained serial/MPI gas and scalar faces are exact')
        do r=1,3
          if(mode>0) then
            if(topology%nodes(node)%level==records(r)%depth) then
              call check(all(records(r)%local_state%fields(node)%values==saved(r)%local_state%fields(node)%values), &
                'late MPI refusal retains selected-level accepted clock states')
            else
              call check(.not.allocated(records(r)%local_state%fields(node)%values), &
                'late MPI refusal retains the sparse accepted-history schema')
            endif
          else
            call check(all(records(r)%source_funding_impulse%fields(node)%values== &
              reference_records(r)%source_funding_impulse%fields(node)%values).and. &
              all(records(r)%source_funding_stock%fields(node)%values== &
              reference_records(r)%source_funding_stock%fields(node)%values),'serial/MPI native funding packets are exact')
            call check(valid_eb_held_funding_context_3d(topology,geometry,nv,local, &
              records(r)%clock_source_funding_context,records(r)%depth),'MPI history preserves owner-sparse open clocks')
          endif
        enddo
      enddo
      call check_settled(context,local)
      if(mode==0) saved=records
      deallocate(records,faces,scalar,gas%fields,temp%fields,held%fields,ir%fields)
    enddo
  end subroutine check_mpi_layout

  subroutine check_settled(context,local)
    type(eb_held_funding_context_3d), intent(in) :: context
    logical, intent(in) :: local(:)
    integer :: node,bucket
    call check(valid_eb_held_funding_context_3d(topology,geometry,nv,local,context,0), &
      'funding context remains valid and settled at the root boundary')
    call check(all(context%open%clock==-1),'partial recovery leaves no open clock')
    do bucket=0,1
      do node=1,2
        if(local(node)) call check(all(context%pending(bucket)%fields(node)%values==0), &
          'partial recovery leaves no pending receipt')
      enddo
    enddo
  end subroutine check_settled

  subroutine check_inventory(gas,held,records,faces,scalar,b,w,ch,ext,aux,correction)
    type(eb_reactive_tree_3d), intent(in) :: gas
    type(amr_tree_state_3d), intent(in) :: held
    type(eb_mol_held_interval_3d), intent(in) :: records(:)
    type(eb_integrated_flux_3d), intent(in) :: faces(:),scalar(:)
    real(dp), intent(in) :: b(:),w(:),ch(:),ext(:),aux(:),correction
    real(dp) :: before(nv),after(nv),residual(nv),face_change(nv),expected(naux),hbefore,hafter,hfaces,hvolume
    real(dp) :: volume,impulse(1),gas_error,h_error,loan_max
    integer :: node,x,y,z,stage,r
    logical :: accepted
    call integrate_reactive_eb_tree_3d(initial,before,accepted);call check(accepted,'entry composite gas inventory')
    call integrate_reactive_eb_tree_3d(gas,after,accepted);call check(accepted,'accepted composite gas inventory')
    face_change=0;expected=0;hbefore=0;hafter=0;hfaces=0;hvolume=0;loan_max=0
    do r=1,3
      call check(valid_eb_held_funding_context_3d(topology,geometry,nv,[.true.,.true.], &
        records(r)%clock_source_funding_context,records(r)%depth),'native pre-entry funding anchors remain valid')
      do node=1,2
        loan_max=max(loan_max,maxval(abs(records(r)%source_funding_impulse%fields(node)%values)))
      enddo
    enddo
    call check(records(2)%clock_source_funding_context%open(0)%clock==0.and. &
      records(3)%clock_source_funding_context%open(0)%clock==0,'both child clocks retain the open root face family')
    do node=1,2
    associate(g=>geometry(node))
    do z=1,g%nz
    do y=1,g%ny
    do x=1,g%nx
      if(g%volume_fraction(x,y,z)<=0) cycle
      call check(gas%nodes(node)%state(1,x,y,z)>0.and.all(gas%nodes(node)%state(6:,x,y,z)>=0).and. &
        held%fields(node)%values(1,x,y,z)>0,'accepted root and child gas/species/H are admissible')
      if(covered(topology,node,[x,y,z])) cycle
      volume=g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz
      face_change=face_change+face_impulse(g,faces(node),x,y,z)
      impulse=face_impulse(g,scalar(node),x,y,z);hfaces=hfaces+impulse(1)
      hbefore=hbefore+volume*initial_h%fields(node)%values(1,x,y,z)
      hafter=hafter+volume*held%fields(node)%values(1,x,y,z)
      do r=1,3
        if(records(r)%depth/=topology%nodes(node)%level) cycle
        impulse=face_impulse(g,records(r)%internal_energy_flux(node),x,y,z)
        hvolume=hvolume+records(r)%interval*volume* &
          records(r)%internal_energy_spatial_rate%fields(node)%values(1,x,y,z)-impulse(1)
      enddo
      if(node==2.and.x==g%nx) then
        do r=2,3
          expected(:nv)=expected(:nv)-.5_dp*sink*volume*records(r)%initial%nodes(node)%state(:,x,y,z)
          expected(naux)=expected(naux)-.5_dp*sink*volume*records(r)%initial_internal_energy%fields(node)%values(1,x,y,z)
        enddo
      endif
    enddo
    enddo
    enddo
    end associate
    enddo
    call check(all(abs(aux-expected)<=32*epsilon(1._dp)*max(1._dp,abs(expected))), &
      'fine physical gas/H withdrawal integrates exactly once across the forest')
    residual=after-before-face_change-w-ch-ext-aux(:nv);residual(5)=residual(5)-correction
    gas_error=maxval(abs(residual)/max(1._dp,abs(before)))
    h_error=abs(hafter-hbefore-hfaces-hvolume-records(1)%held_clock_projection_after- &
      records(1)%held_clock_reset_after)/max(1._dp,abs(hbefore))
    if(rank==0) print *, 'PARTIAL_FINE_INVENTORY',gas_error,h_error,loan_max,aux(1),aux(naux)
    call check(gas_error<3.e-14_dp.and.h_error<3.e-14_dp,'independent composite gas/species and actual H inventories close')
    call check(all(abs(b-face_change)<=32*epsilon(1._dp)*max(1._dp,abs(face_change))), &
      'reported boundary matches independently summed composite physical faces')

    call check(loan_max==0,'this fixture requires no off-level material loan')
  end subroutine check_inventory

  subroutine fill_source(state,interval,stage,depth,time,rate)
    type(amr_tree_state_3d), intent(in) :: state
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    integer :: node
    source_calls=source_calls+1
    if(depth==1.and.stage==0) then
      old_calls=old_calls+1
      call check(time>last_old_time,'each physical fine old source is observed only once')
      last_old_time=time
    endif
    refusal_armed=depth==1.and.time>0.and.stage==1
    do node=1,2
      if(.not.allocated(rate%fields(node)%values)) cycle
      rate%fields(node)%values=0
      if(node==2.and.depth==1.and.stage==0) &
        rate%fields(node)%values(:,4,:,:)=-sink/interval*state%fields(node)%values(:,4,:,:)
    enddo
  end subroutine fill_source

  subroutine serial_source(state,interval,stage,depth,time,rate,accepted)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: accepted
    type(amr_tree_state_3d) :: packed
    integer :: node
    allocate(packed%fields(2))
    do node=1,2
      packed%fields(node)%values=state%nodes(node)%state
    enddo
    call fill_source(packed,interval,stage,depth,time,rate)
    accepted=.not.(callback_mode==1.and.refusal_armed)
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
    call fill_source(state,interval,stage,depth,time,rate)
    accepted=.not.(callback_mode==1.and.refusal_armed.and.rank==1)
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
    call advance_forced_reactor_interval(species,reactions,energy,forcing,energy_forcing,interval,rtol,atol, &
      density,temp,completed_steps,accepted,maximum_steps=maximum_steps)
    if(.not.refusal_armed) return
    call check(accepted,'late caller refusal follows independent native reactor admission')
    reactor_calls=reactor_calls+1;accepted=.false.;density=-1;temp=-1
  end subroutine refusal_reactor
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

  subroutine check(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if(condition) return
    print *, 'FAIL:',rank,trim(label)
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine check

end program test_mpi_eb_native_mol_partial_fine_3d
