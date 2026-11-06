program test_mpi_eb_native_mol_future_source_3d
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: load_h2o2_full_transport
  use gas_transport_mod, only: gas_transport_species
  use mixture_thermo_mod, only: mixture_specific_gas_constant,mixture_mass_properties
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use dual_energy_mod, only: dual_energy_policy,initialize_nasa7_internal_energy_patch
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    build_eb_tree_extended_geometry_3d,advance_reactive_eb_tree_sdc_subcycle_3d,integrate_reactive_eb_tree_3d, &
    eb_mol_held_interval_3d,eb_native_mol_source_prefix_3d,advance_reactive_eb_tree_sdc_3d
  use eb_held_funding_queue_3d_mod, only: eb_held_funding_context_3d,valid_eb_held_funding_context_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_subcycle_3d,mpi_eb_mol_interval_3d, &
    advance_mpi_eb_tree_sdc_3d
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: box(1)
  type(amr_patch_tree_3d) :: topology
  type(eb_reactive_tree_3d) :: initial
  type(amr_tree_state_3d) :: initial_h
  type(eb_geometry_3d) :: finest,padded
  type(eb_geometry_3d), allocatable :: geometry(:),stage_geometry(:)
  type(dual_energy_policy) :: policy
  real(dp), allocatable :: root(:,:,:,:),root_t(:,:,:),thermal(:,:,:),phi(:,:,:),padded_phi(:,:,:)
  real(dp), allocatable :: q(:),cell(:)
  real(dp), parameter :: dt=1e-8_dp,prefix=3.5_dp,material_rate=9e7_dp
  real(dp) :: sound,t,r,source_held_density,mw,rg,cp,cv,gamma,enthalpy,energy,entropy
  logical, parameter :: recursive_future_on=.true.,recursive_material_on=.false.
  integer :: rank,nranks,status,nv,p,j,sp,iterations,layout,owners(2),dims(3)
  integer :: auxiliary_evaluations=0
  integer :: prefix_failure=0,prefix_source_calls=0,prefix_reactor_calls=0,prefix_retry_case=0
  logical :: prefix_forbid_old=.false.,selected_momentum=.false.,selected_pressure=.false.
  real(dp), allocatable :: prefix_packet(:)
  character(len=16) :: faces(6)='outflow'
  logical :: good
  logical :: replay_root,replay_receipt,forbid_sources=.false.,focused_zero
  character(len=16) :: mode,selection
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call get_command_argument(1,mode)
  call get_command_argument(2,selection);focused_zero=selection=='zero-probe'
  selected_pressure=selection=='pressure'
  selected_momentum=selection=='momentum'.or.selected_pressure
  if(selected_momentum.and.mode/='consumer'.and.mode/='retry') error stop 'Momentum option requires consumer or retry'
  replay_root=mode=='replay';replay_receipt=mode=='receipt'
  call require(nranks==3,'distributed owners and idle ranks')
  call load_h2o2_full_thermo(species,good);call require(good,'thermo')
  call load_h2o2_full_mechanism(reactions,good);call require(good,'mechanism')
  call load_h2o2_full_transport(transport,good);call require(good,'transport')
  species%pinned_eos_extrapolation=.true.
  nv=reactive_nvar(size(species));allocate(q(nv),cell(nv))
  q=0;q(1)=1;q(2:4)=[1._dp,2._dp,3._dp];q(6:)=1e-9_dp
  do sp=1,size(species)
    select case(species(sp)%name)
    case('H2');q(sp+5)=.028_dp
    case('O2');q(sp+5)=.224_dp
    case('N2');q(sp+5)=.748_dp-real(size(species)-3,dp)*1e-9_dp
    end select
  enddo
  r=mixture_specific_gas_constant(species,q(6:),good);call require(good,'initial gas constant')
  q(5)=r*1400
  call reactive_primitive_to_conserved(species,q,cell,t,sound,good);call require(good,'initial gas energy')
  allocate(root(nv,4,4,4),root_t(4,4,4),phi(0:8,0:8,0:8),padded_phi(0:32,0:32,0:32))
  do p=1,nv
    root(p,:,:,:)=cell(p)
  enddo
  root_t=t
  do j=0,8
    phi(j,:,:)=real(j,dp)*.0005_dp-.00185_dp
  enddo
  do j=0,32
    padded_phi(j,:,:)=real(max(0,min(8,j-12)),dp)*.0005_dp-.00185_dp
  enddo
  call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[.004_dp,.004_dp,.004_dp],finest,good)
  call require(good,'cut geometry')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,box(1),good);call require(good,'patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1],box,topology,good);call require(good,'tree')
  call initialize_reactive_eb_tree_3d(species,topology,finest,root,root_t,initial,good,internal_energy=initial_h)
  call require(good,'initial gas/UEINT')
  allocate(geometry(2))
  do p=1,2
    geometry(p)=initial%nodes(p)%geometry;dims=topology%nodes(p)%shape
    allocate(thermal(dims(1),dims(2),dims(3)),source=1400._dp)
    call initialize_nasa7_internal_energy_patch(species,initial%nodes(p)%state,thermal,geometry(p)%volume_fraction, &
      initial_h%fields(p)%values(1,:,:,:),good)
    call require(good,'independent held EOS');deallocate(thermal)
  enddo
  call build_level_set_eb_geometry_3d(padded_phi,spread(-.006_dp,1,3),spread(.010_dp,1,3),padded,good)
  call require(good,'padded geometry')
  call build_eb_tree_extended_geometry_3d(topology,geometry,padded,stage_geometry,good,halo_width=6)
  call require(good,'MOL six-cell source context')
  policy%allow_negative_energy=.false.;policy%eta2=2;policy%allow_small_energy=.false.
  if(mode=='prefix'.or.mode=='consumer'.or.mode=='retry') then
    call check_native_source_prefix()
    call MPI_Finalize(status)
    stop
  endif
  cell(6:)=0
  do sp=1,size(species)
    if(species(sp)%name=='N2') cell(5+sp)=1
  enddo
  q=cell;q(2:4)=cell(2:4)/cell(1)
  r=mixture_specific_gas_constant(species,q(6:),good);call require(good,'source inert gas constant')
  q(5)=cell(1)*r*1400._dp
  call reactive_primitive_to_conserved(species,q,cell,t,sound,good);call require(good,'source physical gas EOS')
  call mixture_mass_properties(species,q(6:),1400._dp,mw,rg,cp,cv,gamma,enthalpy,energy,entropy,good)
  call require(good,'source independent held EOS');source_held_density=cell(1)*energy
  call seed_future_donors()
  do iterations=1,2
  do layout=1,2
    if(focused_zero.and.(iterations/=1.or.layout/=1)) cycle
    owners=[0,2]
    if(layout==2) owners=0
    call check_recursive_native_driver()
  enddo
  enddo
  if(rank==0) print *, 'PASS: native MOL actual future fine receipt, gas/UEINT closure and owner layouts'
  call MPI_Finalize(status)
contains
  subroutine check_native_source_prefix()
    type(eb_reactive_tree_3d) :: seed,serial
    type(amr_tree_state_3d) :: seed_h,held,ir,gas,temp
    type(eb_native_mol_source_prefix_3d) :: serial_prefix,mpi_prefix,kept_serial,kept_mpi
    type(eb_native_mol_source_prefix_3d), allocatable :: optional_prefix
    type(eb_integrated_flux_3d), allocatable :: scalar(:),seed_scalar(:),gas_faces(:),seed_gas_faces(:)
    real(dp) :: primitive(nv),inert(nv),temperature,rg,sound,b(nv),w(nv),ch(nv),ext(nv),aux(nv+1),theta,correction
    integer :: nitrogen,p,x,y,z,layout,failure,d(3)
    logical :: good,ready
    nitrogen=0
    do p=1,size(species)
      if(species(p)%name=='N2') nitrogen=p
    enddo
    call prefix_require(nitrogen>0,'prefix fixture identifies full-mechanism nitrogen')
    primitive=0;primitive(1)=1;primitive(5+nitrogen)=1
    rg=mixture_specific_gas_constant(species,primitive(6:),good)
    call prefix_require(good,'prefix inert gas constant')
    primitive(5)=rg*1400._dp
    call reactive_primitive_to_conserved(species,primitive,inert,temperature,sound,good)
    call prefix_require(good,'prefix native inert gas/UEINT input')
    allocate(prefix_packet(nv+1));prefix_packet(:nv)=inert;prefix_packet(nv+1)=inert(5)
    seed=initial;seed_h=initial_h;allocate(seed_scalar(2),seed_gas_faces(2),ir%fields(2))
    do p=1,2
      d=topology%nodes(p)%shape;seed%nodes(p)%state=0;seed%nodes(p)%temperature=0;seed_h%fields(p)%values=0
      do z=1,d(3)
      do y=1,d(2)
      do x=1,d(1)
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        seed%nodes(p)%state(:,x,y,z)=inert;seed%nodes(p)%temperature(x,y,z)=temperature
        seed_h%fields(p)%values(1,x,y,z)=inert(5)
      enddo
      enddo
      enddo
      allocate(seed_gas_faces(p)%x(nv,0:d(1),d(2),d(3)),seed_gas_faces(p)%y(nv,d(1),0:d(2),d(3)), &
        seed_gas_faces(p)%z(nv,d(1),d(2),0:d(3)),source=0._dp)
      allocate(seed_scalar(p)%x(1,0:d(1),d(2),d(3)),seed_scalar(p)%y(1,d(1),0:d(2),d(3)), &
        seed_scalar(p)%z(1,d(1),d(2),0:d(3)),ir%fields(p)%values(nv-4,d(1),d(2),d(3)),source=0._dp)
    enddo
    call prefix_require(geometry(1)%volume_fraction(3,2,2)>0.and..not.covered(topology,1,[3,2,2]), &
      'prefix refusal cell is exposed and fluid')
    do layout=1,2
      owners=[0,2];if(layout==2) owners=0
      serial_prefix%depth=-9;mpi_prefix%depth=-9
      do failure=1,4
        prefix_failure=failure;prefix_source_calls=0;prefix_reactor_calls=0
        if(failure<=3) then
          serial=seed;held=seed_h;scalar=seed_scalar;gas_faces=seed_gas_faces
          do p=1,2
            if(.not.allocated(ir%fields(p)%values)) then
              d=topology%nodes(p)%shape
              allocate(ir%fields(p)%values(nv-4,d(1),d(2),d(3)))
            endif
            ir%fields(p)%values=0
          enddo
          b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp;correction=7;ready=.true.
          if(failure==3) then
          call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,2,.true., &
            [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,serial, &
            stage_geometry,faces,b,w,ch,ext,theta,good,reaction_source=ir,target_depth=0,mol_mode=.true., &
            integrated_flux=gas_faces, &
          internal_energy=held,internal_energy_flux=scalar,energy_policy=policy,energy_correction=correction, &
            auxiliary_source=prefix_serial_source,auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
            evaluate_reaction=prefix_reactor, &
            native_source_prefix=serial_prefix,native_source_prefix_ready=ready, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
          else
          call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,2,.true., &
            [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,serial, &
            stage_geometry,faces,b,w,ch,ext,theta,good,reaction_source=ir,target_depth=0,mol_mode=.true., &
            integrated_flux=gas_faces, &
          internal_energy=held,internal_energy_flux=scalar,energy_policy=policy,energy_correction=correction, &
            auxiliary_source=prefix_serial_source,auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
            native_source_prefix=serial_prefix,native_source_prefix_ready=ready, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
          endif
          call prefix_require(.not.good.and.correction==7.and.all(b==17).and.all(w==19).and.all(ch==23).and. &
            all(ext==29).and.all(aux==31).and.theta==.7_dp,'serial prefix refusal preserves accepted ledgers')
          do p=1,2
            call prefix_require(all(serial%nodes(p)%state==seed%nodes(p)%state).and. &
              all(serial%nodes(p)%temperature==seed%nodes(p)%temperature).and. &
              all(held%fields(p)%values==seed_h%fields(p)%values).and.all(ir%fields(p)%values==0), &
              'serial prefix capture preserves caller gas/T/H/IR')
            call prefix_require(all(gas_faces(p)%x==0).and.all(gas_faces(p)%y==0).and.all(gas_faces(p)%z==0), &
              'serial prefix capture preserves caller gas faces')

            call prefix_require(all(scalar(p)%x==0).and.all(scalar(p)%y==0).and.all(scalar(p)%z==0), &
              'serial prefix capture preserves caller scalar faces')
          enddo
          if(failure==1) then
            call prefix_require(ready.and.prefix_source_calls==1.and.prefix_reactor_calls==0, &
              'serial native prefix captures one observed source without a future call')
            call check_prefix(serial_prefix)
            kept_serial=serial_prefix
          else
            call prefix_require(.not.ready,'serial caller refusal is not a numerical recovery prefix')
            call prefix_require(serial_prefix%depth==kept_serial%depth.and. &
              serial_prefix%accepted_source_stages==kept_serial%accepted_source_stages.and. &
              all(serial_prefix%auxiliary_stages(0)%fields(1)%values== &
                kept_serial%auxiliary_stages(0)%fields(1)%values),'serial declined source retains old prefix')
            call prefix_require(prefix_source_calls==2,'serial declined new source retains one old source without retry')
            if(failure==3) call prefix_require(prefix_reactor_calls==1,'serial native-admitted reactor refusal is not retried')
          endif
        endif
        prefix_source_calls=0;prefix_reactor_calls=0
        allocate(gas%fields(2),temp%fields(2));held=seed_h;scalar=seed_scalar;gas_faces=seed_gas_faces
        do p=1,2
          if(allocated(ir%fields(p)%values)) deallocate(ir%fields(p)%values)
          if(owners(p)/=rank) then
            deallocate(held%fields(p)%values,scalar(p)%x,scalar(p)%y,scalar(p)%z, &
              gas_faces(p)%x,gas_faces(p)%y,gas_faces(p)%z);cycle
          endif
          d=topology%nodes(p)%shape;gas%fields(p)%values=seed%nodes(p)%state
          allocate(temp%fields(p)%values(1,d(1),d(2),d(3)),ir%fields(p)%values(nv-4,d(1),d(2),d(3)),source=0._dp)
          temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
        enddo
        if(failure<4.or.rank==0) allocate(optional_prefix,source=mpi_prefix)
        b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp;correction=7;ready=.true.
        if(failure==3) then
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,stage_geometry,dt,2,.true., &
          [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,faces, &
          b,w,ch,ext,theta,MPI_COMM_WORLD,good,reaction_source=ir,target_depth=0,mol_mode=.true., &
          integrated_flux=gas_faces, &
          internal_energy=held,internal_energy_flux=scalar,energy_policy=policy,energy_correction=correction, &
          auxiliary_source=prefix_mpi_source,auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
          evaluate_reaction=prefix_reactor, &
          native_source_prefix=optional_prefix,native_source_prefix_ready=ready, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
        else
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,stage_geometry,dt,2,.true., &
          [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,faces, &
          b,w,ch,ext,theta,MPI_COMM_WORLD,good,reaction_source=ir,target_depth=0,mol_mode=.true., &
          integrated_flux=gas_faces, &
          internal_energy=held,internal_energy_flux=scalar,energy_policy=policy,energy_correction=correction, &
          auxiliary_source=prefix_mpi_source,auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
          native_source_prefix=optional_prefix,native_source_prefix_ready=ready, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
        endif
        call prefix_require(.not.good.and.correction==7.and.all(b==17).and.all(w==19).and.all(ch==23).and. &
          all(ext==29).and.all(aux==31).and.theta==.7_dp,'MPI prefix refusal preserves owner/idle ledgers')
        do p=1,2
          if(owners(p)/=rank) cycle
          call prefix_require(all(gas%fields(p)%values==seed%nodes(p)%state).and. &
            all(temp%fields(p)%values(1,:,:,:)==seed%nodes(p)%temperature).and. &
            all(held%fields(p)%values==seed_h%fields(p)%values).and.all(ir%fields(p)%values==0), &
            'MPI prefix capture preserves owned gas/T/H/IR')
          call prefix_require(all(gas_faces(p)%x==0).and.all(gas_faces(p)%y==0).and.all(gas_faces(p)%z==0), &
            'MPI prefix capture preserves caller gas faces')

          call prefix_require(all(scalar(p)%x==0).and.all(scalar(p)%y==0).and.all(scalar(p)%z==0), &
            'MPI prefix capture preserves owned scalar faces')
        enddo
        if(failure==1) then
          call prefix_require(ready.and.prefix_source_calls==1.and.prefix_reactor_calls==0, &
            'MPI owner/idle native prefix captures one observed source')
          mpi_prefix=optional_prefix;kept_mpi=mpi_prefix;call check_prefix(mpi_prefix,.true.)
          if(rank==owners(1)) call prefix_require(all(mpi_prefix%auxiliary_stages(0)%fields(1)%values== &
            kept_serial%auxiliary_stages(0)%fields(1)%values).and. &
            all(mpi_prefix%external_stages(0)%fields(1)%values==kept_serial%external_stages(0)%fields(1)%values), &
            'serial/MPI observed native source packets match exactly')
        else
          call prefix_require(.not.ready,'MPI declined or rank-mismatched call is not a native recovery prefix')
          if(allocated(optional_prefix)) then
            call prefix_require(optional_prefix%depth==kept_mpi%depth.and. &
              optional_prefix%accepted_source_stages==kept_mpi%accepted_source_stages, &
              'MPI declined call retains accepted prefix metadata')
            if(rank==owners(1)) call prefix_require(all(optional_prefix%auxiliary_stages(0)%fields(1)%values== &
              kept_mpi%auxiliary_stages(0)%fields(1)%values),'MPI declined call retains accepted source packet')
          endif
          call prefix_require(prefix_source_calls==merge(0,2,failure==4), &
            'MPI source refusal and presence mismatch do not retry or evaluate future sources')
          if(failure==3) call prefix_require(prefix_reactor_calls==merge(1,0,rank==owners(1)), &
            'MPI native-admitted provider refusal occurs once on its owner')
        endif
        if(allocated(optional_prefix)) deallocate(optional_prefix)
        deallocate(gas%fields,temp%fields)
      enddo
      call check_prefix_success(seed,seed_h,seed_scalar,seed_gas_faces)
      if(mode=='consumer') call check_prefix_consumer(seed,seed_h,seed_scalar,seed_gas_faces,kept_serial,kept_mpi)
      if(mode=='retry') call check_recursive_prefix_retry(seed,seed_h,seed_scalar,seed_gas_faces)
      call MPI_Barrier(MPI_COMM_WORLD,status)
      call prefix_require(status==MPI_SUCCESS,'prefix layout checks finish collectively')
      if(rank==0) print *, 'PASS_NATIVE_MOL_SOURCE_PREFIX',layout
    enddo
  end subroutine check_native_source_prefix

  subroutine check_recursive_prefix_retry(seed,seed_h,seed_scalar,seed_faces)
    type(eb_reactive_tree_3d), intent(in) :: seed
    type(amr_tree_state_3d), intent(in) :: seed_h
    type(eb_integrated_flux_3d), intent(in) :: seed_scalar(:),seed_faces(:)
    type(eb_reactive_tree_3d) :: serial,serial_result
    type(amr_tree_state_3d) :: gas,temp,held,ir,serial_h,serial_ir
    type(eb_held_funding_context_3d) :: context,empty_context
    type(eb_integrated_flux_3d), allocatable :: scalar(:),gas_faces(:),serial_scalar(:),serial_faces(:)
    type(eb_mol_held_interval_3d), allocatable :: records(:)
    type(mpi_eb_mol_interval_3d), allocatable :: mpi_records(:)
    integer(int64), allocatable :: retries
    integer(int64) :: steps(2)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(nv+1),theta,correction,pressure,serial_ledgers(5*nv+4)
    integer :: family,trial,p,d(3),expected_calls,r
    logical :: good,refused
    do trial=1,8
    do family=1,2
      if(trial>=7.and.family==1) cycle
      prefix_retry_case=trial;prefix_source_calls=0;prefix_reactor_calls=0;prefix_forbid_old=.false.
      serial=seed;held=seed_h;scalar=seed_scalar;gas_faces=seed_faces;context=empty_context
      allocate(ir%fields(2))
      if(family==2) then
        deallocate(held%fields)
        call localize(seed,seed_h,gas,temp,held)
      endif
      do p=1,2
        if(family==2.and.owners(p)/=rank) then
          deallocate(scalar(p)%x,scalar(p)%y,scalar(p)%z,gas_faces(p)%x,gas_faces(p)%y,gas_faces(p)%z);cycle
        endif
        d=topology%nodes(p)%shape
        allocate(ir%fields(p)%values(nv-4,d(1),d(2),d(3)),source=0._dp)
      enddo
      b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp;correction=7;pressure=11;steps=-17_int64
      if(trial/=7.or.rank==0) allocate(retries,source=-1_int64)
      if(family==1) then
        if(trial==4) then
          call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,2,.true., &
            [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,serial, &
            stage_geometry,faces,b,w,ch,ext,theta,good,reaction_source=ir,mol_mode=.true.,redistribution='none', &
            internal_energy=held,internal_energy_flux=scalar,energy_policy=policy,energy_correction=correction, &
            auxiliary_source=retry_serial_source,auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
            integrated_flux=gas_faces,level_steps=steps,held_interval_history=records,source_funding_context=context, &
            internal_energy_pressure_work_change=pressure,native_common_retries=retries,evaluate_reaction=prefix_reactor, &
        separate_momentum=selected_momentum.neqv.(trial==8.and.rank==0), &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
        else
          call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,2,.true., &
            [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,serial, &
            stage_geometry,faces,b,w,ch,ext,theta,good,reaction_source=ir,mol_mode=.true.,redistribution='none', &
            internal_energy=held,internal_energy_flux=scalar,energy_policy=policy,energy_correction=correction, &
            auxiliary_source=retry_serial_source,auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
            integrated_flux=gas_faces,level_steps=steps,held_interval_history=records,source_funding_context=context, &
            internal_energy_pressure_work_change=pressure,native_common_retries=retries, &
        separate_momentum=selected_momentum.neqv.(trial==8.and.rank==0), &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
        endif
      else
        if(trial==4) then
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp,dt,2,.true., &
            [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,stage_geometry,faces, &
            b,w,ch,ext,theta,MPI_COMM_WORLD,good,reaction_source=ir,mol_mode=.true.,redistribution='none', &
            internal_energy=held,internal_energy_flux=scalar,energy_policy=policy,energy_correction=correction, &
            auxiliary_source=retry_mpi_source,auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
            integrated_flux=gas_faces,level_steps=steps,interval_history=mpi_records,source_funding_context=context, &
            internal_energy_pressure_work_change=pressure,native_common_retries=retries,evaluate_reaction=prefix_reactor, &
        separate_momentum=selected_momentum.neqv.(trial==8.and.rank==0), &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
        else
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,gas,temp,dt,2,.true., &
            [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,stage_geometry,faces, &
            b,w,ch,ext,theta,MPI_COMM_WORLD,good,reaction_source=ir,mol_mode=.true.,redistribution='none', &
            internal_energy=held,internal_energy_flux=scalar,energy_policy=policy,energy_correction=correction, &
            auxiliary_source=retry_mpi_source,auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
            integrated_flux=gas_faces,level_steps=steps,interval_history=mpi_records,source_funding_context=context, &
            internal_energy_pressure_work_change=pressure,native_common_retries=retries, &
        separate_momentum=selected_momentum.neqv.(trial==8.and.rank==0), &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
        endif
      endif
      refused=trial/=6
      call prefix_require(good.neqv.refused,'recursive prefix retry acceptance class')
      expected_calls=2
      if(trial==1.or.trial==2) expected_calls=1
      if(trial==6) expected_calls=6
      if(trial>=7) expected_calls=0
      call prefix_require(prefix_source_calls==expected_calls,'automatic retries evaluate each reached source only once')
      if(allocated(retries)) then
        call prefix_require(retries==merge(13_int64,0_int64,trial==1.or.trial==5), &
          'native retry diagnostics count attempts separately from accepted physical steps')
      endif
      if(refused) then
        call prefix_require(all([b,w,ch,ext,aux,theta,correction,pressure]== &
          [spread(17._dp,1,nv),spread(19._dp,1,nv),spread(23._dp,1,nv),spread(29._dp,1,nv), &
          spread(31._dp,1,nv+1),.7_dp,7._dp,11._dp]).and.all(steps==-17_int64), &
          'failed recursive retries preserve accepted ledgers and step counts')
        call prefix_require(.not.allocated(context%pending).and..not.allocated(context%open).and. &
          .not.allocated(context%topology%nodes).and..not.allocated(context%geometry), &
          'failed recursive retry publishes no funding queue or open family')
        call prefix_require(.not.allocated(records).and..not.allocated(mpi_records),'failed retry publishes no accepted history')
      else
        call prefix_require(all(steps==[1_int64,2_int64]),'successful passthrough advances actual root/two fine clocks')
        if(family==1) then
          call prefix_require(all(records%separate_momentum.eqv.selected_momentum), &
            'serial recursive history retains halo policy')
          serial_result=serial;serial_h=held;serial_ir=ir;serial_scalar=scalar;serial_faces=gas_faces
          serial_ledgers=[b,w,ch,ext,aux,theta,correction,pressure]
        else
          call prefix_require(all([b,w,ch,ext,aux,theta,correction,pressure]==serial_ledgers), &
            'native recursive passthrough ledger parity')
          call prefix_require(all(mpi_records%separate_momentum.eqv.selected_momentum), &
            'MPI recursive history retains halo policy')
          do r=1,size(mpi_records)
          do p=1,2
            call prefix_require(allocated(mpi_records(r)%initial_reaction_source%fields(p)%values).eqv. &
              (owners(p)==rank.and.topology%nodes(p)%level==mpi_records(r)%depth), &
              'prefix capture preserves the selected-level accepted replay schema')
          enddo
          enddo
        endif
      endif
      do p=1,2
        if(family==2.and.owners(p)/=rank) cycle
        if(refused) then
          call prefix_require(all(held%fields(p)%values==seed_h%fields(p)%values).and. &
            all(ir%fields(p)%values==0),'recursive retry preserves caller H/I_R on refusal')
          call prefix_require(same_prefix_faces(scalar(p),seed_scalar(p)).and. &
            same_prefix_faces(gas_faces(p),seed_faces(p)),'recursive retry preserves caller scalar/gas faces')
          if(family==1) then
            call prefix_require(all(serial%nodes(p)%state==seed%nodes(p)%state).and. &
              all(serial%nodes(p)%temperature==seed%nodes(p)%temperature),'recursive serial rollback gas/T')
          else
            call prefix_require(all(gas%fields(p)%values==seed%nodes(p)%state).and. &
              all(temp%fields(p)%values(1,:,:,:)==seed%nodes(p)%temperature),'recursive MPI rollback gas/T')
          endif
        else if(family==2) then
          call prefix_require(all(gas%fields(p)%values==serial_result%nodes(p)%state).and. &
            all(temp%fields(p)%values(1,:,:,:)==serial_result%nodes(p)%temperature).and. &
            all(held%fields(p)%values==serial_h%fields(p)%values).and.all(ir%fields(p)%values==serial_ir%fields(p)%values), &
            'native recursive passthrough gas/T/H/I_R parity')
          call prefix_require(same_prefix_faces(scalar(p),serial_scalar(p)).and. &
            same_prefix_faces(gas_faces(p),serial_faces(p)),'native recursive passthrough face parity')
        endif
      enddo
      if(allocated(retries)) deallocate(retries)
      if(allocated(records)) deallocate(records)
      if(allocated(mpi_records)) deallocate(mpi_records)
      deallocate(ir%fields)
      if(family==2) deallocate(gas%fields,temp%fields)
    enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,status)
    call prefix_require(status==MPI_SUCCESS,'recursive retry layout checks finish collectively')
    if(rank==0) print *, 'PASS_NATIVE_MOL_RECURSIVE_PREFIX_RETRY',owners
  end subroutine check_recursive_prefix_retry

  subroutine retry_serial_source(state,interval,stage,depth,time,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p
    prefix_source_calls=prefix_source_calls+1
    do p=1,size(rate%fields)
      if(.not.allocated(rate%fields(p)%values)) cycle
      rate%fields(p)%values=0
      if(topology%nodes(p)%level/=depth) cycle
      if((prefix_retry_case==1.and.stage==0).or.(prefix_retry_case==5.and.stage==1)) &
        rate%fields(p)%values=-100/interval*state%nodes(p)%state
    enddo
    good=.not.((prefix_retry_case==2.and.stage==0).or.(prefix_retry_case==3.and.stage==1))
  end subroutine retry_serial_source

  subroutine retry_mpi_source(state,temp,interval,stage,depth,time,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temp
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p
    prefix_source_calls=prefix_source_calls+1
    do p=1,size(rate%fields)
      if(.not.allocated(rate%fields(p)%values)) cycle
      rate%fields(p)%values=0
      if(topology%nodes(p)%level/=depth) cycle
      if((prefix_retry_case==1.and.stage==0).or.(prefix_retry_case==5.and.stage==1)) &
        rate%fields(p)%values=-100/interval*state%fields(p)%values
    enddo
    good=.not.(((prefix_retry_case==2.and.stage==0).or.(prefix_retry_case==3.and.stage==1)).and.rank==nranks-1)
  end subroutine retry_mpi_source

  subroutine check_prefix_consumer(seed,seed_h,seed_scalar,seed_faces,observed_serial,observed_mpi)
    type(eb_reactive_tree_3d), intent(in) :: seed
    type(amr_tree_state_3d), intent(in) :: seed_h
    type(eb_integrated_flux_3d), intent(in) :: seed_scalar(:),seed_faces(:)
    type(eb_native_mol_source_prefix_3d), intent(in) :: observed_serial,observed_mpi
    type(eb_native_mol_source_prefix_3d), allocatable :: replay
    type(eb_native_mol_source_prefix_3d) :: template,published,retained
    type(eb_reactive_tree_3d) :: serial,reference_serial
    type(amr_tree_state_3d) :: gas,temp,held,ir,reference_gas,reference_temp,reference_h,reference_ir
    type(eb_mol_held_interval_3d) :: serial_record
    type(mpi_eb_mol_interval_3d) :: mpi_record
    type(eb_integrated_flux_3d), allocatable :: scalar(:),gas_faces(:),reference_scalar(:),reference_faces(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(nv+1),theta,correction,reference_ledgers(5*nv+3),weight
    integer :: family,trial,p,d(3),stages,expected_calls,stage
    logical :: good,ready,refusal
    do family=1,2
      if(family==1) then
        template=observed_serial
      else
        template=observed_mpi
      endif
      do trial=0,13
        prefix_failure=5;prefix_source_calls=0;prefix_reactor_calls=0
        prefix_forbid_old=trial/=0
        serial=seed;held=seed_h;scalar=seed_scalar;gas_faces=seed_faces
        allocate(ir%fields(2))
        if(family==2) allocate(gas%fields(2),temp%fields(2))
        do p=1,2
          d=topology%nodes(p)%shape
          if(family==2.and.owners(p)/=rank) then
            deallocate(held%fields(p)%values,scalar(p)%x,scalar(p)%y,scalar(p)%z, &
              gas_faces(p)%x,gas_faces(p)%y,gas_faces(p)%z);cycle
          endif
          allocate(ir%fields(p)%values(nv-4,d(1),d(2),d(3)),source=0._dp)
          if(family==2) then
            gas%fields(p)%values=seed%nodes(p)%state
            allocate(temp%fields(p)%values(1,d(1),d(2),d(3)))
            temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
          endif
        enddo
        if(trial>0.and.trial/=3) then
          if(trial==12) then
            allocate(replay,source=retained)
          else
            allocate(replay,source=template)
          endif
          if(trial/=9.and.trial/=12) then
            ! Slice actual packets from the native producer control. No future
            ! packet is invented and no accepted trajectory is used as a retry.
            stages=1;if(trial==2) stages=2
            replay%accepted_source_stages=stages
            do stage=0,stages-1
              if(family==1) then
                replay%external_stages(stage)=serial_record%external_stages(stage)
                replay%auxiliary_stages(stage)=serial_record%auxiliary_stages(stage)
              else
                replay%external_stages(stage)=mpi_record%external_stages(stage)
                replay%auxiliary_stages(stage)=mpi_record%auxiliary_stages(stage)
              endif
            enddo
          endif
          if(trial==6.and.(family==1.or.rank==nranks-1)) &
            deallocate(replay%external_stages(0)%fields)
          if(trial==7.and.(family==1.or.rank==nranks-1)) &
            replay%external_stages(2)=replay%external_stages(0)
          if(trial==8.and.family==2.and.rank/=0) deallocate(replay)
          if(trial==8.and.family==1) replay%force(1)=1
          ! A packet generated with a different halo policy cannot be reused,
          ! including a mismatch visible only on one MPI owner.
          if(trial==13.and.(family==1.or.rank==nranks-1)) replay%separate_momentum=.not.selected_momentum
        endif
        if(trial==5) prefix_failure=2
        if(trial==11) prefix_failure=6
        weight=1;if(trial==3.or.trial==4) weight=.5_dp
        b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp;correction=7
        published%depth=-9;ready=.true.
        if(family==1) then
          if(trial==0) then
            call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,2,.true., &
              [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,serial, &
              stage_geometry,faces,b,w,ch,ext,theta,good,reaction_source=ir,target_depth=0,mol_mode=.true., &
              integrated_flux=gas_faces,internal_energy=held,internal_energy_flux=scalar,energy_policy=policy, &
              energy_correction=correction,auxiliary_source=prefix_serial_source,auxiliary_change=aux, &
              auxiliary_has_internal_energy=.true.,held_interval_history=serial_record, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
          else if(trial==3) then
            call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,2,.true., &
              [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,serial, &
              stage_geometry,faces,b,w,ch,ext,theta,good,reaction_source=ir,target_depth=0,mol_mode=.true., &
              integrated_flux=gas_faces,internal_energy=held,internal_energy_flux=scalar,energy_policy=policy, &
              energy_correction=correction,auxiliary_source=prefix_serial_source,auxiliary_change=aux, &
              auxiliary_has_internal_energy=.true.,replay_interval=serial_record, &
              freeze_native_external_sources=.true.,native_spatial_weight=weight, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
          else if(trial==10) then
            call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,2,.true., &
              [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,serial, &
              stage_geometry,faces,b,w,ch,ext,theta,good,reaction_source=ir,target_depth=0,mol_mode=.true., &
              integrated_flux=gas_faces,internal_energy=held,internal_energy_flux=scalar,energy_policy=policy, &
              energy_correction=correction,auxiliary_source=prefix_serial_source,auxiliary_change=aux, &
              auxiliary_has_internal_energy=.true.,native_source_replay=replay,native_spatial_weight=weight, &
              native_source_prefix=published,native_source_prefix_ready=ready,evaluate_reaction=prefix_reactor, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
          else
            call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,2,.true., &
              [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,serial, &
              stage_geometry,faces,b,w,ch,ext,theta,good,reaction_source=ir,target_depth=0,mol_mode=.true., &
              integrated_flux=gas_faces,internal_energy=held,internal_energy_flux=scalar,energy_policy=policy, &
              energy_correction=correction,auxiliary_source=prefix_serial_source,auxiliary_change=aux, &
              auxiliary_has_internal_energy=.true.,native_source_replay=replay,native_spatial_weight=weight, &
              native_source_prefix=published,native_source_prefix_ready=ready, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
          endif
        else
          if(trial==0) then
            call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
              stage_geometry,dt,2,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
              2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta,MPI_COMM_WORLD,good,reaction_source=ir,target_depth=0, &
              mol_mode=.true.,integrated_flux=gas_faces,internal_energy=held,internal_energy_flux=scalar, &
              energy_policy=policy,energy_correction=correction,auxiliary_source=prefix_mpi_source, &
              auxiliary_change=aux,auxiliary_has_internal_energy=.true.,interval_history=mpi_record, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
          else if(trial==3) then
            call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
              stage_geometry,dt,2,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
              2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta,MPI_COMM_WORLD,good,reaction_source=ir,target_depth=0, &
              mol_mode=.true.,integrated_flux=gas_faces,internal_energy=held,internal_energy_flux=scalar, &
              energy_policy=policy,energy_correction=correction,auxiliary_source=prefix_mpi_source, &
              auxiliary_change=aux,auxiliary_has_internal_energy=.true.,replay_interval=mpi_record, &
              freeze_native_external_sources=.true.,native_spatial_weight=weight, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
          else if(trial==10) then
            call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
              stage_geometry,dt,2,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
              2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta,MPI_COMM_WORLD,good,reaction_source=ir,target_depth=0, &
              mol_mode=.true.,integrated_flux=gas_faces,internal_energy=held,internal_energy_flux=scalar, &
              energy_policy=policy,energy_correction=correction,auxiliary_source=prefix_mpi_source, &
              auxiliary_change=aux,auxiliary_has_internal_energy=.true.,native_source_replay=replay, &
              native_spatial_weight=weight,native_source_prefix=published,native_source_prefix_ready=ready, &
              evaluate_reaction=prefix_reactor, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
          else
            call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp, &
              stage_geometry,dt,2,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
              2.e-7_dp,1.e-12_dp,faces,b,w,ch,ext,theta,MPI_COMM_WORLD,good,reaction_source=ir,target_depth=0, &
              mol_mode=.true.,integrated_flux=gas_faces,internal_energy=held,internal_energy_flux=scalar, &
              energy_policy=policy,energy_correction=correction,auxiliary_source=prefix_mpi_source, &
              auxiliary_change=aux,auxiliary_has_internal_energy=.true.,native_source_replay=replay, &
              native_spatial_weight=weight,native_source_prefix=published,native_source_prefix_ready=ready, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
          endif
        endif
        refusal=trial>=5
        expected_calls=0
        if(trial==0) expected_calls=2
        if(trial==1.or.trial==4.or.trial==5.or.trial==10.or.trial==11) expected_calls=1
        call prefix_require(prefix_source_calls==expected_calls,'prefix consumer evaluates only first-reached sources')
        call prefix_require(good.neqv.refusal,'prefix consumer native result classification')
        if(refusal) then
          call prefix_require(all([b,w,ch,ext,aux,theta,correction]== &
            [spread(17._dp,1,nv),spread(19._dp,1,nv),spread(23._dp,1,nv),spread(29._dp,1,nv), &
            spread(31._dp,1,nv+1),.7_dp,7._dp]),'prefix consumer refuses transactionally')
          call prefix_require(ready.eqv.(trial==9.or.trial==11.or.trial==12),'only a native refusal can extend the observed prefix')
          if(trial==11.or.trial==12) then
            call prefix_require(published%accepted_source_stages==2,'accepted new source survives native chemistry refusal')
            call prefix_require(.not.allocated(published%auxiliary_stages(2)%fields), &
              'a native refusal publishes no absent future packet')
            if(trial==11) retained=published
          endif
          if(trial/=9.and.trial<11) call prefix_require(published%depth==-9,'declined trial retains previous prefix output')
        endif
        if(trial==0.or.trial==3) then
          reference_serial=serial;reference_h=held;reference_ir=ir
          if(family==2) then
            reference_gas=gas;reference_temp=temp
          endif
          reference_scalar=scalar;reference_faces=gas_faces;reference_ledgers=[b,w,ch,ext,aux,theta,correction]
        else
          if(.not.refusal) call prefix_require(all([b,w,ch,ext,aux,theta,correction]==reference_ledgers), &
            'prefix consumer preserves native reference ledgers exactly')
          do p=1,2
            if(family==2.and.owners(p)/=rank) cycle
            if(refusal) then
              call prefix_require(all(held%fields(p)%values==seed_h%fields(p)%values).and. &
                all(ir%fields(p)%values==0),'refused prefix consumer preserves independent H and I_R')
              call prefix_require(all(gas_faces(p)%x==0).and.all(gas_faces(p)%y==0).and.all(gas_faces(p)%z==0).and. &
                all(scalar(p)%x==0).and.all(scalar(p)%y==0).and.all(scalar(p)%z==0), &
                'refused prefix consumer preserves caller gas/scalar faces')
              if(family==1) then
                call prefix_require(all(serial%nodes(p)%state==seed%nodes(p)%state).and. &
                  all(serial%nodes(p)%temperature==seed%nodes(p)%temperature),'refused serial entry unchanged')
              else
                call prefix_require(all(gas%fields(p)%values==seed%nodes(p)%state).and. &
                  all(temp%fields(p)%values(1,:,:,:)==seed%nodes(p)%temperature),'refused MPI owner entry unchanged')
              endif
            else
              call prefix_require(all(held%fields(p)%values==reference_h%fields(p)%values).and. &
                all(ir%fields(p)%values==reference_ir%fields(p)%values),'prefix native H/I_R parity')
              call prefix_require(same_prefix_faces(gas_faces(p),reference_faces(p)).and. &
                same_prefix_faces(scalar(p),reference_scalar(p)),'prefix rebuilds native gas/scalar faces exactly')
              if(family==1) then
                call prefix_require(all(serial%nodes(p)%state==reference_serial%nodes(p)%state).and. &
                  all(serial%nodes(p)%temperature==reference_serial%nodes(p)%temperature),'serial prefix native gas/T parity')
              else
                call prefix_require(all(gas%fields(p)%values==reference_gas%fields(p)%values).and. &
                  all(temp%fields(p)%values==reference_temp%fields(p)%values),'MPI prefix native gas/T parity')
              endif
            endif
          enddo
        endif
        if(allocated(replay)) deallocate(replay)
        deallocate(ir%fields)
        if(family==2) deallocate(gas%fields,temp%fields)
      enddo
    enddo
    prefix_forbid_old=.false.
    call MPI_Barrier(MPI_COMM_WORLD,status)
    call prefix_require(status==MPI_SUCCESS,'consumer layout finishes collectively')
    if(rank==0) print *, 'PASS_NATIVE_MOL_PREFIX_CONSUMER',owners
  end subroutine check_prefix_consumer

  subroutine check_prefix_success(seed,seed_h,seed_scalar,seed_gas_faces)
    type(eb_reactive_tree_3d), intent(in) :: seed
    type(amr_tree_state_3d), intent(in) :: seed_h
    type(eb_integrated_flux_3d), intent(in) :: seed_scalar(:),seed_gas_faces(:)
    type(amr_tree_state_3d) :: gas,temp,held,ir,kept_gas,kept_temp,kept_h,kept_ir
    type(eb_integrated_flux_3d), allocatable :: scalar(:),kept_scalar(:),gas_faces(:),kept_gas_faces(:)
    type(eb_native_mol_source_prefix_3d), allocatable :: prefix
    logical, allocatable :: prefix_ready
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),aux(nv+1),theta,correction,kept_ledgers(4*nv+nv+3)
    integer :: trial,p,d(3)
    logical :: good
    prefix_failure=0
    do trial=1,2
      prefix_source_calls=0;allocate(gas%fields(2),temp%fields(2),ir%fields(2))
      held=seed_h;scalar=seed_scalar;gas_faces=seed_gas_faces
      do p=1,2
        if(owners(p)/=rank) then
          deallocate(held%fields(p)%values,scalar(p)%x,scalar(p)%y,scalar(p)%z, &
              gas_faces(p)%x,gas_faces(p)%y,gas_faces(p)%z);cycle
        endif
        d=topology%nodes(p)%shape;gas%fields(p)%values=seed%nodes(p)%state
        allocate(temp%fields(p)%values(1,d(1),d(2),d(3)),ir%fields(p)%values(nv-4,d(1),d(2),d(3)),source=0._dp)
        temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
      enddo
      if(trial==2) then
        allocate(prefix,prefix_ready);prefix%depth=-9;prefix_ready=.true.
      endif
      b=17;w=19;ch=23;ext=29;aux=31;theta=.7_dp;correction=7
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,stage_geometry,dt,2,.true., &
        [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,theta,MPI_COMM_WORLD,good,reaction_source=ir,target_depth=0,mol_mode=.true., &
        integrated_flux=gas_faces, &
          internal_energy=held,internal_energy_flux=scalar,energy_policy=policy,energy_correction=correction, &
        auxiliary_source=prefix_mpi_source,auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
        native_source_prefix=prefix,native_source_prefix_ready=prefix_ready, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
      call prefix_require(good.and.prefix_source_calls==2,'successful native two-correction source schedule remains unchanged')
      if(trial==1) then
        kept_gas=gas;kept_temp=temp;kept_h=held;kept_ir=ir;kept_scalar=scalar;kept_gas_faces=gas_faces
        kept_ledgers=[b,w,ch,ext,aux,theta,correction]
      else
        call prefix_require(.not.prefix_ready.and.prefix%depth==-9,'accepted native interval is not published as a refused prefix')
        call prefix_require(all([b,w,ch,ext,aux,theta,correction]==kept_ledgers), &
          'prefix capture option preserves successful native ledgers exactly')
        do p=1,2
          if(owners(p)/=rank) cycle
          call prefix_require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
            all(temp%fields(p)%values==kept_temp%fields(p)%values).and. &
            all(held%fields(p)%values==kept_h%fields(p)%values).and.all(ir%fields(p)%values==kept_ir%fields(p)%values), &
            'prefix capture option preserves successful native gas/T/H/IR exactly')
          call prefix_require(same_prefix_faces(gas_faces(p),kept_gas_faces(p)), &
            'prefix capture option preserves successful native gas faces exactly')
          call prefix_require(same_prefix_faces(scalar(p),kept_scalar(p)), &
            'prefix capture option preserves successful native scalar faces exactly')
        enddo
        deallocate(prefix,prefix_ready)
      endif
      deallocate(gas%fields,temp%fields,ir%fields)
    enddo
  end subroutine check_prefix_success

  subroutine check_prefix(prefix,parallel)
    type(eb_native_mol_source_prefix_3d), intent(in) :: prefix
    logical, optional, intent(in) :: parallel
    integer :: stage,p
    logical :: owned
    call prefix_require(prefix%separate_momentum.eqv.selected_momentum,'prefix binds selected halo policy')
    call prefix_require(prefix%depth==0.and.prefix%iterations==2.and.prefix%accepted_source_stages==1.and. &
      prefix%interval==dt.and.all(prefix%force==0).and.prefix%chemistry.and. &
      .not.any(prefix%transport_flags).and.prefix%initial_energy_correction==7.and. &
      prefix%has_auxiliary_source.and.prefix%auxiliary_has_internal_energy, &
      'observed prefix retains native source and independent H metadata')
    call prefix_require(lbound(prefix%external_stages,1)==0.and.ubound(prefix%external_stages,1)==2.and. &
      lbound(prefix%auxiliary_stages,1)==0.and.ubound(prefix%auxiliary_stages,1)==2, &
      'observed prefix retains the actual native stage profile')
    do stage=1,2
      call prefix_require(.not.allocated(prefix%external_stages(stage)%fields).and. &
        .not.allocated(prefix%auxiliary_stages(stage)%fields),'unobserved future sources are absent')
    enddo
    do p=1,2
      owned=.true.;if(present(parallel)) owned=rank==owners(p)
      if(.not.owned.or.p/=1) cycle
      call prefix_require(all(prefix%external_stages(0)%fields(p)%values==0),'observed old native force packet is exact')
      call prefix_require(all(prefix%auxiliary_stages(0)%fields(p)%values(:nv,3,2,2)== &
        -2/dt*prefix%initial_gas%fields(p)%values(:,3,2,2)).and. &
        prefix%auxiliary_stages(0)%fields(p)%values(nv+1,3,2,2)== &
        -2/dt*prefix%initial_internal_energy%fields(p)%values(1,3,2,2), &
        'observed native gas/UEINT withdrawal is retained exactly')
    enddo
  end subroutine check_prefix

  subroutine prefix_serial_source(state,interval,stage,depth,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p
    prefix_source_calls=prefix_source_calls+1
    call prefix_require(.not.(prefix_forbid_old.and.stage==0),'observed old packet must never call its provider again')
    do p=1,size(rate%fields)
      rate%fields(p)%values=0
    enddo
    if(prefix_failure==1.and.stage==0) rate%fields(1)%values(:,3,2,2)=-2/interval*state%nodes(1)%state(:,3,2,2)
    if(prefix_failure==5) rate%fields(1)%values(:,3,2,2)=.1_dp*(stage+1)/interval*prefix_packet
    if(prefix_failure==6.and.stage==1) rate%fields(1)%values(:,3,2,2)=-4/interval*prefix_packet
    good=.not.(prefix_failure==2.and.stage==1)
  end subroutine prefix_serial_source

  subroutine prefix_mpi_source(state,temp,interval,stage,depth,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temp
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p
    prefix_source_calls=prefix_source_calls+1
    call prefix_require(.not.(prefix_forbid_old.and.stage==0),'observed old packet must never call its provider again')
    do p=1,size(rate%fields)
      if(allocated(rate%fields(p)%values)) rate%fields(p)%values=0
    enddo
    if(prefix_failure==1.and.stage==0.and.allocated(rate%fields(1)%values)) &
      rate%fields(1)%values(:,3,2,2)=-2/interval*state%fields(1)%values(:,3,2,2)
    if(prefix_failure==5.and.allocated(rate%fields(1)%values)) &
      rate%fields(1)%values(:,3,2,2)=.1_dp*(stage+1)/interval*prefix_packet
    if(prefix_failure==6.and.stage==1.and.allocated(rate%fields(1)%values)) &
      rate%fields(1)%values(:,3,2,2)=-4/interval*prefix_packet
    good=.not.(prefix_failure==2.and.stage==1.and.rank==nranks-1)
  end subroutine prefix_mpi_source

  subroutine prefix_reactor(species,reactions,energy,forcing,energy_forcing,duration,rtol,atol,density,temp, &
      steps,good,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: energy,forcing(:),energy_forcing,duration,rtol,atol
    real(dp), intent(inout) :: density(:),temp
    integer, intent(out) :: steps
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    prefix_reactor_calls=prefix_reactor_calls+1
    call advance_forced_reactor_interval(species,reactions,energy,forcing,energy_forcing,duration,rtol,atol, &
      density,temp,steps,good,maximum_steps=maximum_steps)
    call prefix_require(good,'prefix provider refuses an independently native-admitted input')
    density=-1;temp=-1;steps=0;good=.false.
  end subroutine prefix_reactor

  logical function same_prefix_faces(actual,expected) result(same)
    type(eb_integrated_flux_3d), intent(in) :: actual,expected
    same=.false.
    if(allocated(actual%x).neqv.allocated(expected%x)) return
    if(allocated(actual%y).neqv.allocated(expected%y)) return
    if(allocated(actual%z).neqv.allocated(expected%z)) return
    if(allocated(actual%x)) then
      if(any(shape(actual%x)/=shape(expected%x))) return
      if(any(actual%x/=expected%x)) return
    endif
    if(allocated(actual%y)) then
      if(any(shape(actual%y)/=shape(expected%y))) return
      if(any(actual%y/=expected%y)) return
    endif
    if(allocated(actual%z)) then
      if(any(shape(actual%z)/=shape(expected%z))) return
      if(any(actual%z/=expected%z)) return
    endif
    same=.true.
  end function same_prefix_faces

  subroutine prefix_require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    ! Owner-only checks and per-cell providers cannot enter an all-rank
    ! assertion collective. A local failure still aborts the whole MPI test.
    if(condition) return
    print *, 'FAIL_NATIVE_PREFIX:',rank,trim(label)
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine prefix_require

  subroutine seed_future_donors()
    real(dp) :: primitive(nv),conserved(nv),density,rg,temperature,sound
    integer :: node,i,j,k,n(3),nitrogen
    logical :: accepted
    nitrogen=0
    do i=1,size(species)
      if(species(i)%name=='N2') nitrogen=i+5
    enddo
    call require(nitrogen>5,'future fixture identifies inert nitrogen')
    ! Other exposed coarse cells have no nitrogen inventory. Fine donors
    ! retain the original composition at density4 and1400K, so the true CF
    ! donor path can fund the single coarse source without near-empty donors.
    do node=1,2
      n=topology%nodes(node)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(geometry(node)%volume_fraction(i,j,k)<=0.or.covered(topology,node,[i,j,k])) cycle
        if(node==1.and.all([i,j,k]==[3,2,2])) cycle
        conserved=initial%nodes(node)%state(:,i,j,k)
        primitive=conserved
        primitive(2:4)=conserved(2:4)/conserved(1)
        if(node==1) then
          primitive(nitrogen)=0
          density=sum(primitive(6:))
        else
          primitive(6:)=4*primitive(6:);density=4*conserved(1)
        endif
        primitive(1)=density;primitive(6:)=primitive(6:)/density
        rg=mixture_specific_gas_constant(species,primitive(6:),accepted)
        call require(accepted,'future seed mixture')
        primitive(5)=density*rg*1400._dp
        call reactive_primitive_to_conserved(species,primitive,conserved,temperature,sound,accepted)
        call require(accepted,'future seed consistent gas EOS')
        initial%nodes(node)%state(:,i,j,k)=conserved
        initial%nodes(node)%temperature(i,j,k)=temperature
      enddo
      enddo
      enddo
      allocate(thermal(n(1),n(2),n(3)),source=1400._dp)
      call initialize_nasa7_internal_energy_patch(species,initial%nodes(node)%state,thermal, &
        geometry(node)%volume_fraction,initial_h%fields(node)%values(1,:,:,:),accepted)
      call require(accepted,'future seed independent held EOS');deallocate(thermal)
    enddo
  end subroutine seed_future_donors

  subroutine future_source(rate,depth,stage)
    type(amr_tree_state_3d), intent(inout) :: rate
    integer, intent(in) :: depth,stage
    if(depth/=0.or.stage<=0) return
    if(.not.allocated(rate%fields(1)%values)) return
    ! Outer coarse interval2dt and native corrector half weight give a0.9
    ! pure-N2 withdrawal. MOL2 retains only the zero old auxiliary source.
    rate%fields(1)%values(:nv,3,2,2)=-material_rate*cell
    rate%fields(1)%values(nv+1,3,2,2)=-material_rate*source_held_density
  end subroutine future_source

  function scalar_face_impulse(g,f,i,j,k) result(value)
    type(eb_geometry_3d), intent(in) :: g
    type(eb_integrated_flux_3d), intent(in) :: f
    integer, intent(in) :: i,j,k
    real(dp) :: value
    value=g%dy*g%dz*(g%x_face_fraction(i-1,j,k)*f%x(1,i-1,j,k)-g%x_face_fraction(i,j,k)*f%x(1,i,j,k))+ &
      g%dx*g%dz*(g%y_face_fraction(i,j-1,k)*f%y(1,i,j-1,k)-g%y_face_fraction(i,j,k)*f%y(1,i,j,k))+ &
      g%dx*g%dy*(g%z_face_fraction(i,j,k-1)*f%z(1,i,j,k-1)-g%z_face_fraction(i,j,k)*f%z(1,i,j,k))
  end function scalar_face_impulse

  subroutine check_recursive_native_driver()
    type(eb_reactive_tree_3d) :: result,ledger
    type(eb_reactive_tree_3d) :: replayed
    type(amr_tree_state_3d) :: h,ir,ug,tg,hg,mir
    type(eb_held_funding_context_3d) :: context,owner_context
    type(eb_held_funding_context_3d) :: replay_context,mpi_replay_context
    type(eb_mol_held_interval_3d), allocatable :: records(:)
    type(mpi_eb_mol_interval_3d), allocatable :: owner_records(:)
    type(eb_integrated_flux_3d), allocatable :: gf(:),sf(:),mgf(:),msf(:)
    type(eb_integrated_flux_3d), allocatable :: rgf(:),rsf(:),rmgf(:),rmsf(:)
    type(eb_mol_held_interval_3d), allocatable :: window(:),replay_records(:)
    type(mpi_eb_mol_interval_3d), allocatable :: owner_window(:),mpi_replay_records(:)
    type(amr_tree_state_3d) :: rh,rir,rug,rtg,rhg,rmir
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),source(nv+1),mb(nv),mw(nv),mc(nv),me(nv),msource(nv+1)
    real(dp) :: h0(1),h1(1),adv(1),volume_rate(1),before_gas(nv),after_gas(nv),owner_after(nv)
    real(dp) :: correction,mcorrection,pressure,mpressure,limiter,mlimiter,error,loan_max,pending_max,volume,duration
    real(dp) :: expected_source(nv+1),source_scale
    real(dp), parameter :: weights(3)=[1._dp,.5_dp,0._dp]
    real(dp) :: rb(nv),rw(nv),rc(nv),re(nv),ra(nv+1),rm_b(nv),rm_w(nv),rm_c(nv),rm_e(nv),rm_a(nv+1)
    real(dp) :: rcorrection,rpressure,rtheta,rmcorrection,rmpressure,rmtheta,weight,replay_dt
    integer(int64) :: rsteps(2),rmsteps(2),expected_steps(2)
    integer(int64) :: steps(2),msteps(2)
    integer :: node,n(3),r,i,j,k,bucket,calls,mcalls,index,stage
    logical :: accepted,equal
    result=initial;h=initial_h;correction=prefix;pressure=0
    allocate(ir%fields(2))
    do node=1,2
      n=topology%nodes(node)%shape
      allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)),source=0._dp)
    enddo
    call integrate_reactive_eb_tree_3d(initial,before_gas,accepted)
    call require(accepted,'recursive original gas inventory')
    ledger=initial
    do node=1,2
      ledger%nodes(node)%state=initial_h%fields(node)%values
    enddo
    call integrate_reactive_eb_tree_3d(ledger,h0,accepted)
    call require(accepted,'recursive original independent H inventory')
    auxiliary_evaluations=0
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,2*dt,iterations,.true., &
      [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp,result,stage_geometry,faces, &
      b,w,ch,ext,limiter,accepted,ir,mol_mode=.true.,redistribution='none', &
      internal_energy=h,internal_energy_flux=sf,energy_policy=policy,energy_correction=correction, &
      auxiliary_source=recursive_serial_source,auxiliary_change=source,auxiliary_has_internal_energy=.true., &
      internal_energy_pressure_work_change=pressure,held_interval_history=records, &
      integrated_flux=gf,level_steps=steps,source_funding_context=context, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
    if(rank==0) print *, 'NATIVE_RECURSIVE_SERIAL_ACCEPT',recursive_material_on,iterations,accepted,auxiliary_evaluations
    call require(accepted,'actual native recursive serial gas/H/context accepts')
    calls=auxiliary_evaluations
    call require(size(records)==3.and.all(steps==[1_int64,2_int64]),'native physical root/two fine clocks retained')
    call require(all(records%has_source_funding).and.all(records%held_clock_events_closed).and. &
      all(records%held_clock_prefixes_valid),'accepted funding clocks archive paired receipts and complete scalar events')
    equal=valid_eb_held_funding_context_3d(topology,geometry,nv,[.true.,.true.],context,0)
    do bucket=0,1
      do node=1,2
        equal=equal.and.all(context%pending(bucket)%fields(node)%values==0)
      enddo
    enddo
    call require(equal,'complete native recursive interval returns fully settled queue and no open family')
    loan_max=0;pending_max=0
    do r=1,3
      call require(valid_eb_held_funding_context_3d(topology,geometry,nv,[.true.,.true.], &
        records(r)%clock_source_funding_context,records(r)%depth),'every native clock holds complete pre-entry funding anchor')
      do node=1,2
        loan_max=max(loan_max,maxval(abs(records(r)%source_funding_impulse%fields(node)%values)))
        do bucket=0,1
          pending_max=max(pending_max,maxval(abs(records(r)%clock_source_funding_context% &
            pending(bucket)%fields(node)%values)))
        enddo
      enddo
    enddo
    if(recursive_material_on.and.iterations==1) then
      call require(loan_max>0.and.pending_max>0,'real accepted offclock loan remains in pending stock of subsequent fine clock')
    endif
    if(recursive_future_on.and.iterations==1) then
      call require(records(1)%depth==0.and.records(2)%depth==1,'future donor precedes active fine clock')
      call require(maxval(abs(records(1)%source_funding_impulse%fields(2)%values))>0, &
        'real root source produces future fine donor gas/H debit')
      call require(all(records(2)%clock_source_funding_context%pending(0)%fields(2)%values== &
        records(1)%source_funding_impulse%fields(2)%values),'first fine pre-entry anchor retains actual root debit')
      call require(all(records(2)%initial%nodes(2)%state==records(2)%clock_initial%nodes(2)%state+ &
        records(1)%source_funding_impulse%fields(2)%values(:nv,:,:,:)), &
        'actual native fine old gas pays exactly the queued receipt before forming its predictor')
      call require(all(records(2)%initial_internal_energy%fields(2)%values== &
        records(2)%clock_internal_energy%fields(2)%values+ &
        records(1)%source_funding_impulse%fields(2)%values(nv+1:nv+1,:,:,:)), &
        'actual native fine old UEINT pays exactly the same queued receipt')
      call require(all(records(2)%source_funding_stock%fields(2)%values==0), &
        'paid active fine debit is absent from readonly inactive stock')
      call require(all(records(3)%clock_source_funding_context%pending(0)%fields(2)%values==0), &
        'second fine clock does not pay first fine debit twice')
      if(rank==0) print *, 'PASS_ACTUAL_FUTURE_FINE_ENTRY', &
        maxval(abs(records(1)%source_funding_impulse%fields(2)%values))
    endif
    if(rank==0) print *, 'NATIVE_RECURSIVE_LOAN_PENDING',iterations,loan_max,pending_max
    source_scale=0
    if(recursive_material_on.and.iterations==1) source_scale=1
    if(recursive_future_on.and.iterations==1) source_scale=1
    volume=sum(geometry(2)%volume_fraction)*geometry(2)%dx*geometry(2)%dy*geometry(2)%dz
    if(recursive_future_on) volume=geometry(1)%volume_fraction(3,2,2)*geometry(1)%dx*geometry(1)%dy*geometry(1)%dz
    expected_source(:nv)=-source_scale*material_rate*dt*volume*cell
    expected_source(nv+1)=-source_scale*material_rate*dt*volume*source_held_density
    call require(matching(source,expected_source),'first fine native source integrates once; later fine source is zero')
    call integrate_reactive_eb_tree_3d(result,after_gas,accepted)
    call require(accepted,'accepted native recursive gas inventory')
    after_gas=after_gas-before_gas-b-w-ch-ext-source(:nv)
    after_gas(5)=after_gas(5)-(correction-prefix)
    error=maxval(abs(after_gas)/max(1._dp,abs(before_gas)))
    if(rank==0) print *, 'NATIVE_RECURSIVE_GAS_BUDGET',iterations,error
    call require(error<3e-14_dp,'native recursive gas/species/source/receipt/CF inventory closes independently')
    ! Independent H accounting uses original/final H, complete accepted
    ! faces, accepted local spatial rates and measured projection/reset events.
    ledger=result
    do node=1,2
      ledger%nodes(node)%state=h%fields(node)%values
    enddo
    call integrate_reactive_eb_tree_3d(ledger,h1,accepted)
    call require(accepted,'accepted independent H inventory')
    do node=1,2
      n=topology%nodes(node)%shape;ledger%nodes(node)%state=0
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
        volume=geometry(node)%volume_fraction(i,j,k)*geometry(node)%dx*geometry(node)%dy*geometry(node)%dz
        ledger%nodes(node)%state(1,i,j,k)=scalar_face_impulse(geometry(node),sf(node),i,j,k)/volume
      enddo
      enddo
      enddo
    enddo
    call integrate_reactive_eb_tree_3d(ledger,adv,accepted)
    call require(accepted,'complete physical scalar-face inventory')
    do node=1,2
      ledger%nodes(node)%state=0
    enddo
    do r=1,3
      duration=records(r)%interval
      do node=1,2
        if(topology%nodes(node)%level/=records(r)%depth) cycle
        n=topology%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          volume=geometry(node)%volume_fraction(i,j,k)*geometry(node)%dx*geometry(node)%dy*geometry(node)%dz
          ledger%nodes(node)%state(1,i,j,k)=ledger%nodes(node)%state(1,i,j,k)+ &
            duration*records(r)%internal_energy_spatial_rate%fields(node)%values(1,i,j,k)- &
            scalar_face_impulse(geometry(node),records(r)%internal_energy_flux(node),i,j,k)/volume
        enddo
        enddo
        enddo
      enddo
    enddo
    call integrate_reactive_eb_tree_3d(ledger,volume_rate,accepted)
    call require(accepted,'accepted independent local pressure/volume inventory')
    error=abs(h1(1)-h0(1)-adv(1)-volume_rate(1)- &
      records(1)%held_clock_projection_after-records(1)%held_clock_reset_after)/max(1._dp,abs(h0(1)))
    if(rank==0) print *, 'NATIVE_RECURSIVE_H_BUDGET',iterations,error,h0,h1,adv,volume_rate, &
      records(1)%held_clock_projection_after,records(1)%held_clock_reset_after
    call require(error<3e-14_dp,'actual independent H/face/volume/projection/reset inventory closes')
    call localize(initial,initial_h,ug,tg,hg)
    allocate(mir%fields(2))
    do node=1,2
      if(owners(node)/=rank) cycle
      n=topology%nodes(node)%shape
      allocate(mir%fields(node)%values(nv-4,n(1),n(2),n(3)),source=0._dp)
    enddo
    mcorrection=prefix;mpressure=0;auxiliary_evaluations=0
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,ug,tg,2*dt, &
      iterations,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp, &
      stage_geometry,faces,mb,mw,mc,me,mlimiter,MPI_COMM_WORLD,accepted,mir,mol_mode=.true.,redistribution='none', &
      internal_energy=hg,internal_energy_flux=msf,energy_policy=policy,energy_correction=mcorrection, &
      auxiliary_source=recursive_mpi_source,auxiliary_change=msource,auxiliary_has_internal_energy=.true., &
      internal_energy_pressure_work_change=mpressure,interval_history=owner_records,integrated_flux=mgf, &
      level_steps=msteps,source_funding_context=owner_context, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
    if(rank==0) print *, 'NATIVE_RECURSIVE_MPI_ACCEPT',recursive_material_on,iterations,accepted,auxiliary_evaluations
    call require(accepted,'actual distributed native recursive gas/H/context accepts')
    mcalls=auxiliary_evaluations
    equal=correction==mcorrection.and.pressure==mpressure.and.calls==mcalls.and. &
      all(steps==msteps).and.size(owner_records)==3.and.all(b==mb).and.all(w==mw).and. &
      all(ch==mc).and.all(ext==me).and.all(source==msource)
    do node=1,2
      if(owners(node)/=rank) cycle
      equal=equal.and.all(result%nodes(node)%state==ug%fields(node)%values).and. &
        all(result%nodes(node)%temperature==tg%fields(node)%values(1,:,:,:)).and. &
        all(h%fields(node)%values==hg%fields(node)%values).and. &
        all(gf(node)%x==mgf(node)%x).and.all(gf(node)%y==mgf(node)%y).and.all(gf(node)%z==mgf(node)%z).and. &
        all(sf(node)%x==msf(node)%x).and.all(sf(node)%y==msf(node)%y).and.all(sf(node)%z==msf(node)%z)
      do r=1,3
        equal=equal.and.all(records(r)%source_funding_stock%fields(node)%values== &
          owner_records(r)%source_funding_stock%fields(node)%values).and. &
          all(records(r)%source_funding_impulse%fields(node)%values==owner_records(r)%source_funding_impulse%fields(node)%values)
      enddo
    enddo
    call require(equal,'complete native serial/MPI gas/T/H/receipt/stock/faces/ledger/clock agreement is exact')
    equal=valid_eb_held_funding_context_3d(topology,geometry,nv,owners==rank,owner_context,0)
    do bucket=0,1
      do node=1,2
        if(owners(node)/=rank) cycle
        equal=equal.and.all(owner_context%pending(bucket)%fields(node)%values==0)
      enddo
    enddo
    call require(equal,'distributed final queue fully settles with exact owner/idle mask')
    do r=1,3
      call require(valid_eb_held_funding_context_3d(topology,geometry,nv,owners==rank, &
        owner_records(r)%clock_source_funding_context,owner_records(r)%depth),'MPI pre-entry complete anchor remains owner sparse')
    enddo
    if(replay_root.or.replay_receipt) then
      window=records;owner_window=owner_records;expected_steps=[1_int64,2_int64];replay_dt=2*dt
      if(replay_receipt) then
        window=records(2:2);owner_window=owner_records(2:2);expected_steps=[0_int64,1_int64];replay_dt=dt
      endif
      forbid_sources=.true.
      do index=1,3
        if(focused_zero.and.index/=3) cycle
        weight=weights(index)
        ! Replace a fabricated completed tail using the actual pre-payment
        ! clock anchor, including the root-produced gas/UEINT pending debt.
        replayed=result;rh=h;rir=ir;replay_context=context;rcorrection=999;rpressure=999
        call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,replay_dt,iterations,.true., &
          [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp,replayed,stage_geometry,faces, &
          rb,rw,rc,re,rtheta,accepted,rir,mol_mode=.true.,redistribution='none', &
          internal_energy=rh,internal_energy_flux=rsf,energy_policy=policy,energy_correction=rcorrection, &
          auxiliary_source=recursive_serial_source,auxiliary_change=ra,auxiliary_has_internal_energy=.true., &
          internal_energy_pressure_work_change=rpressure,held_interval_history=replay_records,integrated_flux=rgf, &
          level_steps=rsteps,source_funding_context=replay_context,native_replay_history=window,native_spatial_weight=weight, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
        call require(accepted.and.auxiliary_evaluations==mcalls.and.all(rsteps==expected_steps), &
          'serial material replay accepts complete frozen clock without source callbacks')
        if(weight==1) then
          equal=.true.
          do node=1,2
            if(replay_root) then
              equal=equal.and.all(replayed%nodes(node)%state==result%nodes(node)%state).and. &
                all(replayed%nodes(node)%temperature==result%nodes(node)%temperature).and. &
                all(rh%fields(node)%values==h%fields(node)%values).and.all(rir%fields(node)%values==ir%fields(node)%values)
            else
              equal=equal.and.all(replayed%nodes(node)%state==records(3)%clock_initial%nodes(node)%state).and. &
                all(replayed%nodes(node)%temperature==records(3)%clock_initial%nodes(node)%temperature).and. &
                all(rh%fields(node)%values==records(3)%clock_internal_energy%fields(node)%values).and. &
                all(rir%fields(node)%values==records(3)%clock_reaction_source%fields(node)%values)
            endif
          enddo
          call require(equal,'unit material replay matches actual native next-clock gas/T/UEINT/IR bitwise')
        endif
        do r=1,size(window)
          do stage=0,iterations
            do node=1,2
              call require(all(replay_records(r)%auxiliary_stages(stage)%fields(node)%values== &
                window(r)%auxiliary_stages(stage)%fields(node)%values).and. &
                all(replay_records(r)%external_stages(stage)%fields(node)%values== &
                window(r)%external_stages(stage)%fields(node)%values),'material replay retains actual frozen source packets')
            enddo
          enddo
        enddo
        call check_replay_inventory(window(1),replayed,rh,replay_records,replay_context,rgf,rsf,rw,rc,re,ra,rcorrection)
        rug=ug;rtg=tg;rhg=hg;rmir=mir;mpi_replay_context=owner_context;rmcorrection=999;rmpressure=999
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,rug,rtg,replay_dt, &
          iterations,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp, &
          stage_geometry,faces,rm_b,rm_w,rm_c,rm_e,rmtheta,MPI_COMM_WORLD,accepted,rmir, &
          mol_mode=.true.,redistribution='none',internal_energy=rhg,internal_energy_flux=rmsf,energy_policy=policy, &
          energy_correction=rmcorrection,auxiliary_source=recursive_mpi_source,auxiliary_change=rm_a, &
          auxiliary_has_internal_energy=.true.,internal_energy_pressure_work_change=rmpressure, &
          interval_history=mpi_replay_records,integrated_flux=rmgf,level_steps=rmsteps, &
          source_funding_context=mpi_replay_context,native_replay_history=owner_window,native_spatial_weight=weight, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_pressure,paired_pressure=selected_pressure)
        call require(accepted.and.auxiliary_evaluations==mcalls.and.all(rmsteps==expected_steps), &
          'MPI material replay accepts complete frozen clock without source callbacks')
        equal=rcorrection==rmcorrection.and.rpressure==rmpressure.and.rtheta==rmtheta.and. &
          all(rb==rm_b).and.all(rw==rm_w).and.all(rc==rm_c).and.all(re==rm_e).and.all(ra==rm_a)
        do node=1,2
          if(owners(node)/=rank) cycle
          equal=equal.and.all(replayed%nodes(node)%state==rug%fields(node)%values).and. &
            all(replayed%nodes(node)%temperature==rtg%fields(node)%values(1,:,:,:)).and. &
            all(rh%fields(node)%values==rhg%fields(node)%values).and.all(rir%fields(node)%values==rmir%fields(node)%values).and. &
            all(rgf(node)%x==rmgf(node)%x).and.all(rgf(node)%y==rmgf(node)%y).and.all(rgf(node)%z==rmgf(node)%z).and. &
            all(rsf(node)%x==rmsf(node)%x).and.all(rsf(node)%y==rmsf(node)%y).and.all(rsf(node)%z==rmsf(node)%z)
          do bucket=0,1
            equal=equal.and.all(replay_context%pending(bucket)%fields(node)%values== &
              mpi_replay_context%pending(bucket)%fields(node)%values)
          enddo
          do r=1,size(window)
            equal=equal.and.all(replay_records(r)%source_funding_stock%fields(node)%values== &
              mpi_replay_records(r)%source_funding_stock%fields(node)%values).and. &
              all(replay_records(r)%source_funding_impulse%fields(node)%values== &
              mpi_replay_records(r)%source_funding_impulse%fields(node)%values)
          enddo
        enddo
        call require(equal,'material replay serial/MPI gas/T/UEINT/IR/faces/debts/budgets agree bitwise')
        if(rank==0) print *, 'PASS_NATIVE_MOL_MATERIAL_REPLAY',trim(mode),iterations,layout,weight
      enddo
      forbid_sources=.false.
    endif
    if(rank==0) print *, 'PASS_NATIVE_MOL_FUTURE_SOURCE',iterations,layout
  end subroutine check_recursive_native_driver

  subroutine check_replay_inventory(anchor,endpoint,held,history,context,gas_faces,held_faces,wall,chemical,external,aux,correction)
    type(eb_mol_held_interval_3d), intent(in) :: anchor,history(:)
    type(eb_reactive_tree_3d), intent(in) :: endpoint
    type(amr_tree_state_3d), intent(in) :: held
    type(eb_held_funding_context_3d), intent(in) :: context
    type(eb_integrated_flux_3d), intent(in) :: gas_faces(:),held_faces(:)
    real(dp), intent(in) :: wall(:),chemical(:),external(:),aux(:),correction
    type(eb_reactive_tree_3d) :: ledger
    real(dp) :: before(nv),after(nv),face_change(nv),residual(nv),pending(nv+1),hbefore(1),hafter(1)
    real(dp) :: hfaces(1),hvolume(1),volume,error,herror
    integer :: node,i,j,k,r,bucket
    logical :: accepted
    call integrate_reactive_eb_tree_3d(anchor%clock_initial,before,accepted);call require(accepted,'replay entry gas inventory')
    call integrate_reactive_eb_tree_3d(endpoint,after,accepted);call require(accepted,'replay endpoint gas inventory')
    ledger=endpoint
    do node=1,2
      ledger%nodes(node)%state=anchor%clock_internal_energy%fields(node)%values
    enddo
    call integrate_reactive_eb_tree_3d(ledger,hbefore,accepted);call require(accepted,'replay entry actual UEINT inventory')
    do node=1,2
      ledger%nodes(node)%state=held%fields(node)%values
    enddo
    call integrate_reactive_eb_tree_3d(ledger,hafter,accepted);call require(accepted,'replay endpoint actual UEINT inventory')
    ! Pending physical receipts are real conserved debt at an open clock.
    ! Include each bucket once; no hypothetical future flux is material stock.
    do bucket=0,1
      do node=1,2
        ledger%nodes(node)%state=anchor%clock_source_funding_context%pending(bucket)%fields(node)%values
      enddo
      call integrate_reactive_eb_tree_3d(ledger,pending,accepted);call require(accepted,'replay pending entry inventory')
      before=before+pending(:nv);hbefore=hbefore+pending(nv+1)
      do node=1,2
        ledger%nodes(node)%state=context%pending(bucket)%fields(node)%values
      enddo
      call integrate_reactive_eb_tree_3d(ledger,pending,accepted);call require(accepted,'replay pending endpoint inventory')
      after=after+pending(:nv);hafter=hafter+pending(nv+1)
    enddo
    do node=1,2
      allocate_dummy: block
        real(dp) :: impulse(nv)
        ledger%nodes(node)%state=endpoint%nodes(node)%state*0
        associate(g=>geometry(node),f=>gas_faces(node))
        do k=1,g%nz
        do j=1,g%ny
        do i=1,g%nx
          if(g%volume_fraction(i,j,k)<=0) cycle
          volume=g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz
          impulse=g%dy*g%dz*(g%x_face_fraction(i-1,j,k)*f%x(:,i-1,j,k)-g%x_face_fraction(i,j,k)*f%x(:,i,j,k))+ &
            g%dx*g%dz*(g%y_face_fraction(i,j-1,k)*f%y(:,i,j-1,k)-g%y_face_fraction(i,j,k)*f%y(:,i,j,k))+ &
            g%dx*g%dy*(g%z_face_fraction(i,j,k-1)*f%z(:,i,j,k-1)-g%z_face_fraction(i,j,k)*f%z(:,i,j,k))
          ledger%nodes(node)%state(:,i,j,k)=impulse/volume
        enddo
        enddo
        enddo
        end associate
      end block allocate_dummy
    enddo
    call integrate_reactive_eb_tree_3d(ledger,face_change,accepted);call require(accepted,'replay gas face inventory')
    residual=after-before-face_change-wall-chemical-external-aux(:nv)
    residual(5)=residual(5)-(correction-anchor%clock_energy_correction)
    error=maxval(abs(residual)/max(1._dp,abs(before)))
    call require(error<3e-14_dp,'material replay independent gas/species/source/debt inventory closes')
    do node=1,2
      ledger%nodes(node)%state=held%fields(node)%values*0
      do k=1,geometry(node)%nz
      do j=1,geometry(node)%ny
      do i=1,geometry(node)%nx
        if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
        volume=geometry(node)%volume_fraction(i,j,k)*geometry(node)%dx*geometry(node)%dy*geometry(node)%dz
        ledger%nodes(node)%state(1,i,j,k)=scalar_face_impulse(geometry(node),held_faces(node),i,j,k)/volume
      enddo
      enddo
      enddo
    enddo
    call integrate_reactive_eb_tree_3d(ledger,hfaces,accepted);call require(accepted,'replay scalar face inventory')
    do node=1,2
      ledger%nodes(node)%state=0
    enddo
    do r=1,size(history)
      do node=1,2
        if(topology%nodes(node)%level/=history(r)%depth) cycle
        do k=1,geometry(node)%nz
        do j=1,geometry(node)%ny
        do i=1,geometry(node)%nx
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          volume=geometry(node)%volume_fraction(i,j,k)*geometry(node)%dx*geometry(node)%dy*geometry(node)%dz
          ledger%nodes(node)%state(1,i,j,k)=ledger%nodes(node)%state(1,i,j,k)+ &
            history(r)%interval*history(r)%internal_energy_spatial_rate%fields(node)%values(1,i,j,k)- &
            scalar_face_impulse(geometry(node),history(r)%internal_energy_flux(node),i,j,k)/volume
        enddo
        enddo
        enddo
      enddo
    enddo
    call integrate_reactive_eb_tree_3d(ledger,hvolume,accepted);call require(accepted,'replay scalar volume inventory')
    herror=abs(hafter(1)-hbefore(1)-hfaces(1)-hvolume(1)- &
      history(1)%held_clock_projection_after-history(1)%held_clock_reset_after)/max(1._dp,abs(hbefore(1)))
    call require(herror<3e-14_dp,'material replay independent actual UEINT/debt inventory closes')
    if(rank==0) print *, 'NATIVE_MOL_MATERIAL_REPLAY_INVENTORY',trim(mode),iterations,error,herror
  end subroutine check_replay_inventory

  subroutine recursive_serial_source(state,interval,stage,depth,time,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: node
    do node=1,size(rate%fields)
      if(allocated(rate%fields(node)%values)) rate%fields(node)%values=0
    enddo
    good=.not.forbid_sources;auxiliary_evaluations=auxiliary_evaluations+1
    if(recursive_future_on.and.time==0) call future_source(rate,depth,stage)
  end subroutine recursive_serial_source

  subroutine recursive_mpi_source(state,temperature,interval,stage,depth,time,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: node
    do node=1,size(rate%fields)
      if(allocated(rate%fields(node)%values)) rate%fields(node)%values=0
    enddo
    good=.not.forbid_sources;auxiliary_evaluations=auxiliary_evaluations+1
    if(recursive_future_on.and.time==0) call future_source(rate,depth,stage)
  end subroutine recursive_mpi_source

  logical function matching(a,b)
    real(dp), intent(in) :: a(:),b(:)
    matching=all(abs(a-b)<=32*epsilon(1._dp)*max(abs(a),abs(b),1e-30_dp))
  end function matching

  subroutine localize(all_gas,all_h,u,t,h)
    type(eb_reactive_tree_3d), intent(in) :: all_gas
    type(amr_tree_state_3d), intent(in) :: all_h
    type(amr_tree_state_3d), intent(out) :: u,t,h
    integer :: p,d(3)
    allocate(u%fields(2),t%fields(2),h%fields(2))
    do p=1,2
      if(owners(p)/=rank) cycle
      d=topology%nodes(p)%shape
      u%fields(p)%values=all_gas%nodes(p)%state
      allocate(t%fields(p)%values(1,d(1),d(2),d(3)))
      t%fields(p)%values(1,:,:,:)=all_gas%nodes(p)%temperature
      h%fields(p)%values=all_h%fields(p)%values
    enddo
  end subroutine localize

  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    logical :: agreed
    integer :: ierr
    call MPI_Allreduce(condition,agreed,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) then
      print *, 'FAIL:',rank,trim(label)
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end subroutine require
end program test_mpi_eb_native_mol_future_source_3d
