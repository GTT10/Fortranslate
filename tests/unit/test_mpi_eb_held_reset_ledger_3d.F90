program test_mpi_eb_held_reset_ledger_3d
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: load_h2o2_full_transport
  use gas_transport_mod, only: gas_transport_species
  use mixture_thermo_mod, only: mixture_specific_gas_constant,temperature_from_internal_energy,mixture_mass_properties
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use dual_energy_mod, only: dual_energy_policy,initialize_nasa7_internal_energy_patch
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    build_eb_tree_extended_geometry_3d,advance_reactive_eb_tree_sdc_3d,integrate_reactive_eb_tree_3d, &
    eb_godunov_interval_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,compute_accepted_eb_held_clock_prefixes_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: box(1)
  type(amr_patch_tree_3d) :: topology
  type(eb_reactive_tree_3d) :: initial,serial
  type(amr_tree_state_3d) :: initial_h,serial_h,gas,temp,held
  type(eb_geometry_3d) :: finest,padded
  type(eb_geometry_3d), allocatable :: geometry(:),extended(:),mol_extended(:),stage_geometry(:)
  type(eb_integrated_flux_3d), allocatable :: serial_faces(:),mpi_faces(:)
  type(dual_energy_policy) :: policy
  real(dp), allocatable :: root(:,:,:,:),root_t(:,:,:),thermal(:,:,:),phi(:,:,:),padded_phi(:,:,:)
  real(dp), allocatable :: q(:),cell(:),before(:),after(:),local_after(:),forcing(:),boundary(:),wall(:),chemical(:),external(:)
  real(dp), allocatable :: mb(:),mw(:),mc(:),me(:)
  real(dp), allocatable :: aux(:),maux(:)
  real(dp) :: sound,t,r,serial_prefix,mpi_prefix,serial_work,mpi_work,serial_projection,mpi_projection,theta,mt,residual
  real(dp) :: serial_closure,mpi_closure,maximum_closure=0
  real(dp), parameter :: dt=1e-8_dp,prefix=3.5_dp
  real(dp), parameter :: material_rate=1.5e8_dp
  character(len=32) :: source_kind
  logical :: material_source,covered_source,corrector_source,diffusive_source,source_forbidden=.false.
  logical :: reactor_mode,caloric_provider,kinetic_reactor,builtin_reactor,reactor_audit_ok=.true.
  real(dp) :: caloric_rate,mix_mw,mix_rg,mix_cp,mix_cv,mix_gamma,mix_h,mix_e,mix_s
  integer :: caloric_denials=0,reactor_replenishments=0
  character(len=16) :: reactor_request,reset_request
  integer :: reactor_audits=0
  integer :: seed_observations=0
  real(dp) :: reactor_minimum_temperature=huge(1._dp)
  logical :: source_transport_flags(4)
  type(eb_godunov_interval_3d), allocatable :: serial_record,mpi_record
  integer :: auxiliary_evaluations=0
  integer :: rank,nranks,status,nv,p,j,sp,mode,iterations,owners(2),dims(3),halo_width,padded_count
  character(len=16) :: faces(6)='outflow',redistribution_method,requested_mode
  logical :: good,mol,same
  call get_command_argument(1,requested_mode)
  call get_command_argument(2,source_kind)
  call get_command_argument(3,reactor_request)
  call get_command_argument(4,reset_request)
  kinetic_reactor=reactor_request=='kinetic'.or.reactor_request=='builtin'
  builtin_reactor=reactor_request=='builtin'
  reactor_mode=reactor_request=='reactor'.or.kinetic_reactor
  caloric_provider=source_kind=='caloric_provider'
  diffusive_source=source_kind=='diffusive_corrector_material'
  source_transport_flags=.false.
  if(diffusive_source) source_transport_flags(:3)=.true.
  corrector_source=source_kind=='corrector_material'.or.source_kind=='covered_corrector_material'.or.diffusive_source
  corrector_source=corrector_source.or.caloric_provider
  material_source=source_kind=='material'.or.source_kind=='covered_material'.or.corrector_source
  if(caloric_provider) material_source=.false.
  covered_source=source_kind=='covered_material'.or.source_kind=='covered_corrector_material'
  redistribution_method='none'
  halo_width=6
  if(requested_mode=='registers') then
    redistribution_method='state_redist';halo_width=9
  endif
  padded_count=8+4*halo_width
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call require(nranks==3,'two owners and one idle rank')
  call check_accepted_clock_prefixes()
  owners=[0,2]
  call load_h2o2_full_thermo(species,good);call require(good,'thermo')
  call load_h2o2_full_mechanism(reactions,good);call require(good,'mechanism')
  call load_h2o2_full_transport(transport,good);call require(good,'transport')
  species%pinned_eos_extrapolation=requested_mode=='eos_source'
  nv=reactive_nvar(size(species))
  allocate(aux(nv+1),maux(nv+1))
  allocate(q(nv),cell(nv),before(nv),after(nv),local_after(nv),forcing(nv),boundary(nv),wall(nv),chemical(nv),external(nv))
  allocate(mb(nv),mw(nv),mc(nv),me(nv))
  q=0;q(1)=1;q(2:4)=[1._dp,2._dp,3._dp]
  do sp=1,size(species)
    if(species(sp)%name=='N2') q(sp+5)=1
  enddo
  if(kinetic_reactor) then
    q(6:)=1e-9_dp
    do sp=1,size(species)
      select case(species(sp)%name)
      case('H2');q(sp+5)=.028_dp
      case('O2');q(sp+5)=.224_dp
      case('N2');q(sp+5)=.748_dp-real(size(species)-3,dp)*1e-9_dp
      end select
    enddo
  endif
  if(caloric_provider) then
    call require(reactor_mode,'manufactured caloric response needs the supplied reactor')
    q=0;q(1)=10;q(6:)=1e-6_dp
    do sp=1,size(species)
      select case(species(sp)%name)
      case('H2');q(sp+5)=.028_dp
      case('O2');q(sp+5)=.224_dp
      case('N2');q(sp+5)=.748_dp-real(size(species)-3,dp)*1e-6_dp
      end select
    enddo
  endif
  r=mixture_specific_gas_constant(species,q(6:),good);call require(good,'mixture gas constant')
  q(5)=r*1000
  if(kinetic_reactor) q(5)=r*1800
  if(caloric_provider) q(5)=q(1)*r*1400
  call reactive_primitive_to_conserved(species,q,cell,t,sound,good);call require(good,'initial gas energy')
  if(caloric_provider) then
    call mixture_mass_properties(species,q(6:),200._dp,mix_mw,mix_rg,mix_cp,mix_cv,mix_gamma, &
      mix_h,mix_e,mix_s,good)
    call require(good,'manufactured caloric forced endpoint')
    caloric_rate=(q(1)*(mix_e-1e6_dp)-cell(5))/dt
  endif
  allocate(root(nv,4,4,4),root_t(4,4,4),phi(0:8,0:8,0:8),padded_phi(0:padded_count,0:padded_count,0:padded_count))
  do p=1,nv
    root(p,:,:,:)=cell(p)
  enddo
  root_t=t
  do j=0,8
    phi(:,j,:)=real(j,dp)*.0005_dp-.00125_dp
  enddo
  do j=0,padded_count
    padded_phi(:,j,:)=real(max(0,min(8,j-2*halo_width)),dp)*.0005_dp-.00125_dp
  enddo
  if(requested_mode=='eos_source') then
    do j=0,8
      phi(j,:,:)=real(j,dp)*.0005_dp-.00185_dp
    enddo
    do j=0,padded_count
      padded_phi(j,:,:)=real(max(0,min(8,j-2*halo_width)),dp)*.0005_dp-.00185_dp
    enddo
  endif
  call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[.004_dp,.004_dp,.004_dp],finest,good)
  call require(good,'cut geometry')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,box(1),good);call require(good,'patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1],box,topology,good);call require(good,'tree')
  call initialize_reactive_eb_tree_3d(species,topology,finest,root,root_t,initial,good,internal_energy=initial_h)
  call require(good,'initial gas/UEINT')
  allocate(geometry(2))
  do p=1,2
    geometry(p)=initial%nodes(p)%geometry;dims=topology%nodes(p)%shape
    allocate(thermal(dims(1),dims(2),dims(3)),source=700._dp)
    if(caloric_provider.or.kinetic_reactor) thermal=1400
    call initialize_nasa7_internal_energy_patch(species,initial%nodes(p)%state,thermal,geometry(p)%volume_fraction, &
      initial_h%fields(p)%values(1,:,:,:),good)
    call require(good,'independent700K UEINT');deallocate(thermal)
  enddo
  call build_level_set_eb_geometry_3d(padded_phi,spread(-real(halo_width,dp)*.001_dp,1,3), &
    spread(.004_dp+real(halo_width,dp)*.001_dp,1,3),padded,good)
  call require(good,'padded geometry')
  call build_eb_tree_extended_geometry_3d(topology,geometry,padded,extended,good,halo_width=halo_width)
  call require(good,'Godunov source context')
  mol_extended=extended
  if(halo_width==9) then
    do p=1,2
      dims=topology%nodes(p)%shape
      call extract_eb_geometry_3d(extended(p),[4,4,4],dims+15,mol_extended(p),good)
      call require(good,'MOL six-cell source context')
    enddo
  endif
  policy%allow_negative_energy=.false.;policy%eta2=2
  if(reset_request=='total_closure') then
    call require(requested_mode=='eos_source'.and.source_kind=='corrector_material'.and..not.reactor_mode, &
      'total-energy closure countercase uses the independently qualified source and inert policy')
    policy%eta2=0
  endif
  ! Match the retained native cooling diagnostic: signed formation energies
  ! are admissible above the EOS energy at the pinned very-small temperature.
  if(requested_mode=='eos_source') policy%allow_small_energy=.false.
  forcing=0;forcing(5)=1e8_dp
  call integrate_reactive_eb_tree_3d(initial,before,good);call require(good,'initial composite inventory')
  if(reset_request=='seed') then
    call check_seed_admission()
    call MPI_Finalize(status)
    stop
  endif
  if(requested_mode=='eos_source') call check_policy_refusal()
  do mode=0,1
    mol=mode==1
    if(mol) then
      stage_geometry=mol_extended
    else
      stage_geometry=extended
    endif
    do iterations=1,2
      if(allocated(serial_record)) deallocate(serial_record)
      if(allocated(mpi_record)) deallocate(mpi_record)
      if(material_source.and..not.mol) allocate(serial_record,mpi_record)
      reactor_audit_ok=.true.;reactor_audits=0;reactor_minimum_temperature=huge(1._dp)
      serial=initial;serial_h=initial_h;serial_prefix=prefix;serial_work=42;serial_projection=43;serial_closure=44
      if(reactor_mode.and..not.builtin_reactor) then
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,iterations,reactor_mode, &
        source_transport_flags,forcing,1e-8_dp,1e-12_dp,serial,stage_geometry,faces, &
        boundary,wall,chemical,external,theta,good, &
        redistribution=redistribution_method,mol_mode=mol, &
        internal_energy=serial_h,internal_energy_flux=serial_faces,energy_policy=policy,energy_correction=serial_prefix, &
        internal_energy_pressure_work_change=serial_work, &
        internal_energy_projection_change=serial_projection, &
        internal_energy_closure_change=serial_closure,auxiliary_source=serial_source, &
        auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
        evaluate_reaction=audit_reactor,godunov_interval_history=serial_record)
      else
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,iterations,reactor_mode, &
        source_transport_flags,forcing,1e-8_dp,1e-12_dp,serial,stage_geometry,faces, &
        boundary,wall,chemical,external,theta,good, &
        redistribution=redistribution_method,mol_mode=mol, &
        internal_energy=serial_h,internal_energy_flux=serial_faces,energy_policy=policy,energy_correction=serial_prefix, &
        internal_energy_pressure_work_change=serial_work, &
        internal_energy_projection_change=serial_projection, &
        internal_energy_closure_change=serial_closure,auxiliary_source=serial_source, &
        auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
        godunov_interval_history=serial_record)
      endif
      call require(good,'serial accepted held SDC iterate')
      if(caloric_provider) call check_caloric_admission()
      if(kinetic_reactor) call check_kinetic_changes()
      if(reactor_mode.and..not.builtin_reactor) &
        call require(reactor_audit_ok.and.reactor_audits>0,'serial actual reactor energy audited')
      if(requested_mode=='eos_source') call check_held_inventory()
      if(material_source) then
        associate(g=>geometry(2))
          call require(matching(aux(:nv),-accepted_material_scale()*material_rate*dt* &
            g%volume_fraction(4,2,2)*g%dx*g%dy*g%dz*cell), &
            'actual accepted complete material source counted once on physical leaves')
        end associate
        call require(serial%nodes(2)%state(1,4,2,2)>0.and.serial%nodes(2)%temperature(4,2,2)>0, &
          'funded fine source endpoint has physical density and temperature')
      endif
      call integrate_reactive_eb_tree_3d(serial,after,good);call require(good,'serial final inventory')
      after(5)=after(5)-(serial_prefix-prefix)
      residual=ledger_residual(after,boundary+wall+chemical+external+aux(:nv))
      if(rank==0.and.requested_mode=='eos_source') then
        print *, 'Accepted gas source ledger:',mol,iterations,residual,serial_prefix-prefix
        print *, 'Gas inventory/current/exchange:',before(5),after(5),boundary(5),wall(5),chemical(5),external(5),aux(5)
        print *, 'Cold cell H/T/E:',serial_h%fields(2)%values(1,4,2,2), &
          serial%nodes(2)%temperature(4,2,2),serial%nodes(2)%state(5,4,2,2)
      endif
      call require(residual<=3e-12_dp,'serial reset ledger counts actual accepted forcing only')
      call localize(initial,initial_h,gas,temp,held)
      mpi_prefix=prefix;mpi_work=42;mpi_projection=43;mpi_closure=44
      reactor_audit_ok=.true.;reactor_audits=0
      if(reactor_mode.and..not.builtin_reactor) then
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,stage_geometry, &
        dt,iterations,reactor_mode,source_transport_flags,forcing,1e-8_dp,1e-12_dp, &
        faces,mb,mw,mc,me,mt,MPI_COMM_WORLD,good,redistribution=redistribution_method,mol_mode=mol, &
        internal_energy=held,internal_energy_flux=mpi_faces,energy_policy=policy,energy_correction=mpi_prefix, &
        internal_energy_pressure_work_change=mpi_work, &
        internal_energy_projection_change=mpi_projection, &
        internal_energy_closure_change=mpi_closure,auxiliary_source=native_mpi_auxiliary, &
        auxiliary_change=maux,auxiliary_has_internal_energy=.true., &
        evaluate_reaction=audit_reactor,godunov_interval_history=mpi_record)
      else
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,stage_geometry, &
        dt,iterations,reactor_mode,source_transport_flags,forcing,1e-8_dp,1e-12_dp, &
        faces,mb,mw,mc,me,mt,MPI_COMM_WORLD,good,redistribution=redistribution_method,mol_mode=mol, &
        internal_energy=held,internal_energy_flux=mpi_faces,energy_policy=policy,energy_correction=mpi_prefix, &
        internal_energy_pressure_work_change=mpi_work, &
        internal_energy_projection_change=mpi_projection, &
        internal_energy_closure_change=mpi_closure,auxiliary_source=native_mpi_auxiliary, &
        auxiliary_change=maux,auxiliary_has_internal_energy=.true., &
        godunov_interval_history=mpi_record)
      endif
      call require(good,'MPI accepted held SDC iterate')
      if(reactor_mode.and..not.builtin_reactor) then
        call require(reactor_audit_ok,'MPI actual reactor energy audited')
        call require((rank==1.and.reactor_audits==0).or.(rank/=1.and.reactor_audits>0), &
          'owner-only reactor calls leave idle rank empty')
      endif
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,gas,local_after,MPI_COMM_WORLD,good)
      call require(good,'MPI final inventory')
      local_after(5)=local_after(5)-(mpi_prefix-prefix)
      residual=ledger_residual(local_after,mb+mw+mc+me+maux(:nv))
      call require(residual<=3e-12_dp,'MPI reset ledger counts actual accepted forcing only')
      same=matching(boundary,mb).and.matching(wall,mw).and.matching(chemical,mc).and.matching(external,me)
      same=same.and.abs(serial_prefix-mpi_prefix)<=16*epsilon(1._dp)*max(1._dp,abs(serial_prefix))
      same=same.and.matching([serial_work],[mpi_work]).and.matching(aux,maux)
      call require(all(ieee_is_finite([serial_projection,mpi_projection])), 'accepted local scalar projection is finite')
      call require(abs(serial_projection-mpi_projection)<=3e-12_dp*abs(before(5)), &
        'owner-neutral local scalar projection agrees with serial physical integral')
      if(allocated(serial_record)) call require(serial_projection==serial_record%internal_energy_projection_change, &
        'accepted serial Godunov history owns exact local scalar projection')
      if(allocated(mpi_record)) call require(mpi_projection==mpi_record%internal_energy_projection_change, &
        'accepted owner Godunov history owns replicated local scalar projection')
      call require(all(ieee_is_finite([serial_closure,mpi_closure])), 'accepted spatial closure scalar reset is finite')
      call require(abs(serial_closure-mpi_closure)<=3e-12_dp*abs(before(5)), &
        'owner-neutral scalar closure reset agrees with serial physical integral')
      if(allocated(serial_record)) call require(serial_closure==serial_record%internal_energy_closure_change, &
        'accepted serial Godunov history owns exact scalar closure reset')
      if(allocated(mpi_record)) call require(mpi_closure==mpi_record%internal_energy_closure_change, &
        'accepted owner Godunov history owns replicated scalar closure reset')
      if(requested_mode=='eos_source'.and..not.material_source.and..not.caloric_provider) then
        call require(serial%nodes(2)%temperature(4,2,2)>0.and. &
          serial%nodes(2)%temperature(4,2,2)<300,'native positive extrapolated fine-cell temperature')
        call require(aux(nv+1)<0,'accepted cooling source integral')
      endif
      do p=1,2
        if(owners(p)/=rank) cycle
        same=same.and.all(gas%fields(p)%values==serial%nodes(p)%state)
        same=same.and.all(temp%fields(p)%values(1,:,:,:)==serial%nodes(p)%temperature)
        same=same.and.all(held%fields(p)%values==serial_h%fields(p)%values)
      enddo
      call require(same,'serial/MPI state, actual UEINT, accepted ledgers and pressure work agree')
      if(material_source.and..not.mol) call check_funded_godunov_replay()
      maximum_closure=max(maximum_closure,abs(serial_closure))
      if(rank==0) print *, 'Accepted held reset ledger:',trim(redistribution_method),mol,iterations,residual,mpi_prefix-prefix
    enddo
  enddo
  if(rank==0) print *, 'PASS: serial/MPI accepted-iterate held reset bookkeeping'
  if(reset_request=='total_closure') call require(maximum_closure>3e-12_dp*abs(before(5)), &
    'actual scalar closure reset is resolved above the unchanged independent ledger tolerance')
  call MPI_Finalize(status)
contains
  subroutine check_accepted_clock_prefixes()
    real(dp), allocatable :: pb(:),rb(:),pa(:),ra(:)
    real(dp) :: saved_pb(3),saved_rb(3),saved_pa(3),saved_ra(3)
    logical :: valid
    call compute_accepted_eb_held_clock_prefixes_3d([0,1,1],[1._dp,2._dp,4._dp], &
      [1._dp,2._dp,4._dp],[8._dp,16._dp,32._dp],[64._dp,128._dp,256._dp],0._dp,0._dp,pb,rb,pa,ra,valid)
    call require(valid,'accepted parent/sibling clock prefix construction')
    ! Entry root: reset1, projection1/reset8. First child: reset2,
    ! projection2/reset16, exit128. Second child: reset4,
    ! projection4/reset32, exit256. Parent exit64 is last, not first.
    call require(all(pb==[0._dp,1._dp,3._dp]).and.all(rb==[1._dp,11._dp,159._dp]).and. &
      all(pa==[7._dp,3._dp,7._dp]).and.all(ra==[511._dp,155._dp,447._dp]), &
      'parent EOS reset follows descendants in immutable accepted scalar prefixes')
    call compute_accepted_eb_held_clock_prefixes_3d([1,2,1],[2._dp,4._dp,8._dp], &
      [2._dp,4._dp,8._dp],[16._dp,32._dp,64._dp],[128._dp,256._dp,512._dp],1._dp,9._dp,pb,rb,pa,ra,valid)
    call require(valid.and.all(pb==[1._dp,3._dp,7._dp]).and.all(rb==[11._dp,31._dp,455._dp]).and. &
      all(pa==[7._dp,7._dp,15._dp]).and.all(ra==[447._dp,319._dp,1031._dp]), &
      'accepted deep-subtree replacement retains the original external prefix once')
    saved_pb=pb;saved_rb=rb;saved_pa=pa;saved_ra=ra
    call compute_accepted_eb_held_clock_prefixes_3d([0,2,1],[0._dp,0._dp,0._dp], &
      [0._dp,0._dp,0._dp],[0._dp,0._dp,0._dp],[0._dp,0._dp,0._dp],0._dp,0._dp,pb,rb,pa,ra,valid)
    call require(.not.valid.and.all(pb==saved_pb).and.all(rb==saved_rb).and.all(pa==saved_pa).and. &
      all(ra==saved_ra),'invalid clock-depth jump leaves all four accepted prefix arrays unchanged')
    call compute_accepted_eb_held_clock_prefixes_3d([0],[0._dp],[huge(1._dp)],[0._dp],[0._dp], &
      huge(1._dp),0._dp,pb,rb,pa,ra,valid)
    call require(.not.valid.and.all(pb==saved_pb).and.all(rb==saved_rb).and.all(pa==saved_pa).and. &
      all(ra==saved_ra),'unrepresentable scalar clock sum refuses atomically without overflow')
    call compute_accepted_eb_held_clock_prefixes_3d([0,1],[0._dp],[0._dp,0._dp],[0._dp,0._dp], &
      [0._dp,0._dp],0._dp,0._dp,pb,rb,pa,ra,valid)
    call require(.not.valid.and.all(pb==saved_pb).and.all(rb==saved_rb).and.all(pa==saved_pa).and. &
      all(ra==saved_ra),'mismatched scalar clock packet count preserves the accepted prefix arrays')
  end subroutine check_accepted_clock_prefixes

  real(dp) function ledger_residual(current,exchange) result(relative)
    real(dp), intent(in) :: current(:),exchange(:)
    real(dp) :: scale(size(current))
    scale=max(abs(before),abs(current),1e-30_dp)
    ! The zero-velocity mixture has cancellation-sized pressure impulses.
    ! Use initial mass times thermal speed as its momentum scale; leave
    ! every original fixture and the 3e-12 balance criterion unchanged.
    if(caloric_provider) scale(2:4)=max(scale(2:4),abs(before(1))*sqrt(r*1400._dp))
    relative=maxval(abs(current-before-exchange)/scale)
  end function ledger_residual
  subroutine check_kinetic_changes()
    real(dp) :: atoms(4),changes(4),mass_scale
    integer :: sp
    call require(maxval(abs(chemical(6:)))>1e-15_dp*abs(before(1)), &
      'pinned standard kinetics produces nonzero accepted chemical species change')
    call require(abs(sum(chemical(6:)))<3e-12_dp*abs(before(1)),'accepted kinetic species conserve mass')
    changes=0
    do sp=1,size(species)
      select case(species(sp)%name)
      case('H2');atoms=[2._dp,0._dp,0._dp,0._dp]
      case('H');atoms=[1._dp,0._dp,0._dp,0._dp]
      case('O');atoms=[0._dp,1._dp,0._dp,0._dp]
      case('O2');atoms=[0._dp,2._dp,0._dp,0._dp]
      case('OH');atoms=[1._dp,1._dp,0._dp,0._dp]
      case('H2O');atoms=[2._dp,1._dp,0._dp,0._dp]
      case('HO2');atoms=[1._dp,2._dp,0._dp,0._dp]
      case('H2O2');atoms=[2._dp,2._dp,0._dp,0._dp]
      case('N2');atoms=[0._dp,0._dp,2._dp,0._dp]
      case('AR');atoms=[0._dp,0._dp,0._dp,1._dp]
      case default
        call require(.false.,'known pinned H2/O2 species for elemental conservation')
      end select
      changes=changes+atoms*chemical(sp+5)/species(sp)%molecular_weight
    enddo
    mass_scale=abs(before(1))/minval(species%molecular_weight)
    call require(maxval(abs(changes))<3e-12_dp*mass_scale,'accepted kinetic species conserve H O N Ar atoms')
    if(rank==0) print *, 'Actual kinetic chemical change/atom residual:',mol,iterations, &
      maxval(abs(chemical(6:)))/abs(before(1)),maxval(abs(changes))/mass_scale
  end subroutine check_kinetic_changes
  subroutine check_caloric_admission()
    integer :: node,i,j,k,n(3)
    call require(all(aux(:4)==0).and.all(aux(6:nv)==0),'caloric source has no material components')
    if(accepted_material_scale()>0) call require(caloric_denials>0,'provider accepts refused caloric preview')
    if(.not.mol.and.iterations==2) call require(reactor_replenishments>0,'reaction replenishes negative forced water')
    if(iterations/=1) return
    do node=1,size(serial%nodes)
      n=topology%nodes(node)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(geometry(node)%volume_fraction(i,j,k)<=0.or.covered(topology,node,[i,j,k])) cycle
        call require(abs(serial%nodes(node)%state(1,i,j,k)-initial%nodes(node)%state(1,i,j,k))< &
          3e-12_dp*q(1),'caloric provider accepts without unnecessary material loans')
      enddo
      enddo
      enddo
    enddo
  end subroutine check_caloric_admission
  subroutine check_funded_godunov_replay()
    type(eb_integrated_flux_3d), allocatable :: gas_faces(:)
    integer :: evaluations,node
    logical :: valid,same
    if(diffusive_source) then
      same=.false.
      do node=1,2
        same=same.or.any(serial_record%sources(1)%rate%fields(node)%values/=0)
      enddo
      call require(same,'diffusive corrector includes an actual nonzero new-time diffusion rate')
      same=.false.
      do node=1,2
        same=same.or.any(serial_record%sources(1)%faces(node)%x/=0).or. &
          any(serial_record%sources(1)%faces(node)%y/=0).or.any(serial_record%sources(1)%faces(node)%z/=0)
      enddo
      call require(same,'diffusive corrector includes actual nonzero spatial faces')
    endif
    serial=initial;serial_h=initial_h;serial_prefix=prefix;serial_work=42;serial_projection=43;serial_closure=44
    evaluations=auxiliary_evaluations;source_forbidden=.true.
    if(reactor_mode.and..not.builtin_reactor) then
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,iterations,reactor_mode, &
      source_transport_flags,forcing,1e-8_dp,1e-12_dp,serial,stage_geometry,faces, &
      boundary,wall,chemical,external,theta,valid,redistribution=redistribution_method,mol_mode=.false., &
      internal_energy=serial_h,internal_energy_flux=serial_faces,energy_policy=policy,energy_correction=serial_prefix, &
      internal_energy_pressure_work_change=serial_work, &
        internal_energy_projection_change=serial_projection, &
        internal_energy_closure_change=serial_closure,auxiliary_source=serial_source, &
      auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
        evaluate_reaction=audit_reactor,godunov_replay_interval=serial_record, &
      integrated_flux=gas_faces)
    else
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,iterations,reactor_mode, &
      source_transport_flags,forcing,1e-8_dp,1e-12_dp,serial,stage_geometry,faces, &
      boundary,wall,chemical,external,theta,valid,redistribution=redistribution_method,mol_mode=.false., &
      internal_energy=serial_h,internal_energy_flux=serial_faces,energy_policy=policy,energy_correction=serial_prefix, &
      internal_energy_pressure_work_change=serial_work, &
        internal_energy_projection_change=serial_projection, &
        internal_energy_closure_change=serial_closure,auxiliary_source=serial_source, &
      auxiliary_change=aux,auxiliary_has_internal_energy=.true., &
        godunov_replay_interval=serial_record, &
      integrated_flux=gas_faces)
    endif
    source_forbidden=.false.
    call require(valid.and.evaluations==auxiliary_evaluations,'serial funded replay accepted without callbacks')
    same=.true.
    do node=1,2
      same=same.and.all(serial%nodes(node)%state==serial_record%final_gas%fields(node)%values)
      same=same.and.all(serial%nodes(node)%temperature==serial_record%final_temperature%fields(node)%values(1,:,:,:))
      same=same.and.all(serial_h%fields(node)%values==serial_record%final_internal_energy%fields(node)%values)
      same=same.and.all(gas_faces(node)%x==serial_record%faces(node)%x).and. &
        all(gas_faces(node)%y==serial_record%faces(node)%y).and.all(gas_faces(node)%z==serial_record%faces(node)%z)
      same=same.and.all(serial_faces(node)%x==serial_record%scalar_faces(node)%x).and. &
        all(serial_faces(node)%y==serial_record%scalar_faces(node)%y).and. &
        all(serial_faces(node)%z==serial_record%scalar_faces(node)%z)
    enddo
    same=same.and.all(boundary==serial_record%boundary_change).and.all(wall==serial_record%wall_change).and. &
      all(chemical==serial_record%chemical_change).and.all(external==serial_record%source_change).and. &
      all(aux==serial_record%auxiliary_change)
    same=same.and.serial_prefix==serial_record%initial_energy_correction+serial_record%numerical_energy_change.and. &
      serial_work==serial_record%pressure_work_change.and.theta==serial_record%minimum_theta
    same=same.and.serial_projection==serial_record%internal_energy_projection_change
    same=same.and.serial_closure==serial_record%internal_energy_closure_change
    call require(same,'serial funded gas H temperatures faces and ledgers replay exactly')
    deallocate(gas_faces)
    call localize(initial,initial_h,gas,temp,held)
    mpi_prefix=prefix;mpi_work=42;mpi_projection=43;mpi_closure=44;evaluations=auxiliary_evaluations;source_forbidden=.true.
    if(reactor_mode.and..not.builtin_reactor) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,stage_geometry, &
      dt,iterations,reactor_mode,source_transport_flags,forcing,1e-8_dp,1e-12_dp, &
      faces,mb,mw,mc,me,mt,MPI_COMM_WORLD,valid,redistribution=redistribution_method,mol_mode=.false., &
      internal_energy=held,internal_energy_flux=mpi_faces,energy_policy=policy,energy_correction=mpi_prefix, &
      internal_energy_pressure_work_change=mpi_work, &
        internal_energy_projection_change=mpi_projection, &
        internal_energy_closure_change=mpi_closure,auxiliary_source=native_mpi_auxiliary, &
      auxiliary_change=maux,auxiliary_has_internal_energy=.true., &
        evaluate_reaction=audit_reactor,godunov_replay_interval=mpi_record, &
      integrated_flux=gas_faces)
    else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,stage_geometry, &
      dt,iterations,reactor_mode,source_transport_flags,forcing,1e-8_dp,1e-12_dp, &
      faces,mb,mw,mc,me,mt,MPI_COMM_WORLD,valid,redistribution=redistribution_method,mol_mode=.false., &
      internal_energy=held,internal_energy_flux=mpi_faces,energy_policy=policy,energy_correction=mpi_prefix, &
      internal_energy_pressure_work_change=mpi_work, &
        internal_energy_projection_change=mpi_projection, &
        internal_energy_closure_change=mpi_closure,auxiliary_source=native_mpi_auxiliary, &
      auxiliary_change=maux,auxiliary_has_internal_energy=.true., &
        godunov_replay_interval=mpi_record, &
      integrated_flux=gas_faces)
    endif
    source_forbidden=.false.
    call require(valid.and.evaluations==auxiliary_evaluations,'MPI funded replay accepted without callbacks')
    same=.true.
    do node=1,2
      if(owners(node)/=rank) cycle
      same=same.and.all(gas%fields(node)%values==mpi_record%final_gas%fields(node)%values)
      same=same.and.all(temp%fields(node)%values==mpi_record%final_temperature%fields(node)%values)
      same=same.and.all(held%fields(node)%values==mpi_record%final_internal_energy%fields(node)%values)
      same=same.and.all(gas_faces(node)%x==mpi_record%faces(node)%x).and. &
        all(gas_faces(node)%y==mpi_record%faces(node)%y).and.all(gas_faces(node)%z==mpi_record%faces(node)%z)
      same=same.and.all(mpi_faces(node)%x==mpi_record%scalar_faces(node)%x).and. &
        all(mpi_faces(node)%y==mpi_record%scalar_faces(node)%y).and.all(mpi_faces(node)%z==mpi_record%scalar_faces(node)%z)
    enddo
    same=same.and.all(mb==mpi_record%boundary_change).and.all(mw==mpi_record%wall_change).and. &
      all(mc==mpi_record%chemical_change).and.all(me==mpi_record%source_change).and.all(maux==mpi_record%auxiliary_change)
    same=same.and.mpi_prefix==mpi_record%initial_energy_correction+mpi_record%numerical_energy_change.and. &
      mpi_work==mpi_record%pressure_work_change.and.mt==mpi_record%minimum_theta
    same=same.and.mpi_projection==mpi_record%internal_energy_projection_change
    same=same.and.mpi_closure==mpi_record%internal_energy_closure_change
    ! Idle ranks participate once, after owner-local inspection, in exactly
    ! the same collective order as both owners.
    call require(same,'MPI funded gas H temperatures faces and ledgers replay exactly with idle rank')
  end subroutine check_funded_godunov_replay

  subroutine check_policy_refusal()
    integer :: trial,node
    logical :: same,valid
    do trial=1,2
      call localize(initial,initial_h,gas,temp,held)
      species%pinned_eos_extrapolation=.true.
      species%pinned_eos_temperature_min=0;species%pinned_eos_temperature_max=0
      if(rank==2) then
        if(trial==1) then
          species%pinned_eos_temperature_min=90;species%pinned_eos_temperature_max=4000
        else
          species%pinned_eos_extrapolation=.false.
        endif
      endif
      auxiliary_evaluations=0
      mb=17;mw=19;mc=23;me=29;maux=31;mpi_prefix=prefix;mpi_work=42;mpi_projection=43;mpi_closure=44
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,extended, &
        dt,1,.false.,source_transport_flags,forcing,1e-8_dp,1e-12_dp, &
        faces,mb,mw,mc,me,mt,MPI_COMM_WORLD,valid,redistribution='none', &
        internal_energy=held,internal_energy_flux=mpi_faces,energy_policy=policy,energy_correction=mpi_prefix, &
        internal_energy_pressure_work_change=mpi_work, &
        internal_energy_projection_change=mpi_projection, &
        internal_energy_closure_change=mpi_closure,auxiliary_source=native_mpi_auxiliary, &
        auxiliary_change=maux,auxiliary_has_internal_energy=.true.)
      same=.not.valid.and.auxiliary_evaluations==0.and..not.allocated(mpi_faces)
      same=same.and.mpi_prefix==prefix.and.mpi_work==42.and.all(maux==31)
      same=same.and.mpi_projection==43.and.mpi_closure==44
      same=same.and.all(mb==17).and.all(mw==19).and.all(mc==23).and.all(me==29)
      do node=1,2
        if(owners(node)/=rank) cycle
        same=same.and.all(gas%fields(node)%values==initial%nodes(node)%state)
        same=same.and.all(temp%fields(node)%values(1,:,:,:)==initial%nodes(node)%temperature)
        same=same.and.all(held%fields(node)%values==initial_h%fields(node)%values)
      enddo
      call require(same,'rank-disagreeing EOS policy refuses all state and ledger publication')
    enddo
    species%pinned_eos_extrapolation=.true.
    species%pinned_eos_temperature_min=0;species%pinned_eos_temperature_max=0
  end subroutine check_policy_refusal
  subroutine audit_reactor(table,rates,old_energy,external_species,external_energy,duration,rtol,atol, &
      rhoy,temperature,steps,ok,maximum_steps)
    type(nasa7_species), intent(in) :: table(:)
    type(elementary_reaction), intent(in) :: rates(:)
    real(dp), intent(in) :: old_energy,external_species(:),external_energy,duration,rtol,atol
    real(dp), intent(inout) :: rhoy(:),temperature
    integer, intent(out) :: steps
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    real(dp) :: endpoint(size(rhoy)),rho,energy,test_t,original_h
    logical :: valid
    steps=0;ok=.false.;reactor_audits=reactor_audits+1
    if(kinetic_reactor) then
      original_h=initial_h%fields(2)%values(1,4,2,2)
      valid=abs(old_energy-original_h)<=32*epsilon(1._dp)*max(1._dp,abs(original_h))
      if(.not.valid) then
        reactor_audit_ok=.false.;return
      endif
      energy=old_energy+duration*external_energy
      ! Let the pinned kinetic model replenish signed forced species. Its
      ! actual accepted endpoint, not a pre-reaction preview, is audited.
      call advance_forced_reactor_interval(table,rates,old_energy,external_species,external_energy,duration,rtol,atol, &
        rhoy,temperature,steps,ok,maximum_steps=maximum_steps)
      if(ok) then
        rho=sum(rhoy);valid=rho>0.and.all(rhoy>=0)
        if(valid) call temperature_from_internal_energy(table,rhoy/rho,energy/rho,temperature,test_t,valid)
        reactor_audit_ok=reactor_audit_ok.and.valid
        if(valid) reactor_minimum_temperature=min(reactor_minimum_temperature,test_t)
      endif
      reactor_audit_ok=reactor_audit_ok.and.ok
      return
    endif
    if(caloric_provider) then
      call manufacture_caloric_response(table,old_energy,external_species,external_energy,duration, &
        rhoy,temperature,steps,ok)
      reactor_audit_ok=reactor_audit_ok.and.ok
      return
    endif
    ! No collective checks here: an idle MPI rank need not enter a reactor.
    original_h=initial_h%fields(2)%values(1,4,2,2)
    valid=abs(old_energy-original_h)<=32*epsilon(1._dp)*max(1._dp,abs(original_h))
    endpoint=rhoy+duration*external_species;rho=sum(endpoint)
    valid=valid.and.rho>0.and.all(endpoint>=0)
    if(.not.valid) then
      reactor_audit_ok=.false.;return
    endif
    energy=old_energy+duration*external_energy;test_t=temperature
    call temperature_from_internal_energy(table,endpoint/rho,energy/rho,temperature,test_t,valid)
    if(.not.valid.or.test_t<=0) then
      reactor_audit_ok=.false.;return
    endif
    reactor_minimum_temperature=min(reactor_minimum_temperature,test_t)
    call advance_forced_reactor_interval(table,rates,old_energy,external_species,external_energy,duration,rtol,atol, &
      rhoy,temperature,steps,ok,maximum_steps=maximum_steps)
    reactor_audit_ok=reactor_audit_ok.and.ok
  end subroutine audit_reactor

  subroutine manufacture_caloric_response(table,old_energy,external_species,external_energy,duration, &
      rhoy,temperature,steps,ok)
    type(nasa7_species), intent(in) :: table(:)
    real(dp), intent(in) :: old_energy,external_species(:),external_energy,duration
    real(dp), intent(inout) :: rhoy(:),temperature
    integer, intent(out) :: steps
    logical, intent(out) :: ok
    real(dp) :: endpoint(size(rhoy)),rho,energy,test_t,extent,original_h,preview_t
    integer :: h2,o2,water,sp
    logical :: valid
    steps=0;ok=.false.
    original_h=initial_h%fields(2)%values(1,4,2,2)
    if(abs(old_energy-original_h)>32*epsilon(1._dp)*max(1._dp,abs(original_h))) then
      reactor_audit_ok=.false.;return
    endif
    endpoint=rhoy+duration*external_species;rho=sum(endpoint)
    h2=0;o2=0;water=0
    do sp=1,size(table)
      select case(table(sp)%name)
      case('H2');h2=sp
      case('O2');o2=sp
      case('H2O');water=sp
      end select
    enddo
    if(min(h2,o2,water)==0.or.rho<=0) then
      reactor_audit_ok=.false.;return
    endif
    ! Water is produced by this reaction and may have a negative forced
    ! endpoint before chemistry, exactly like signed native reactor forcing.
    ! No unrelated negative species is accepted or clipped.
    do sp=1,size(table)
      if(sp==water) cycle
      if(endpoint(sp)<0) then
        reactor_audit_ok=.false.;return
      endif
    enddo
    energy=old_energy+duration*external_energy
    if(all(endpoint>=0)) then
      call temperature_from_internal_energy(table,endpoint/rho,energy/rho,temperature,preview_t,valid)
      if(.not.valid) caloric_denials=caloric_denials+1
    endif
    if(endpoint(water)<0) reactor_replenishments=reactor_replenishments+1
    ! 2H2+O2 -> 2H2O, with exact pinned molecular weights. The fixed extent
    ! is a manufactured integrator response, not a chemical rate model.
    extent=.75_dp*min(endpoint(h2)/(2*table(h2)%molecular_weight),endpoint(o2)/table(o2)%molecular_weight)
    endpoint(h2)=endpoint(h2)-2*extent*table(h2)%molecular_weight
    endpoint(o2)=endpoint(o2)-extent*table(o2)%molecular_weight
    endpoint(water)=endpoint(water)+2*extent*table(water)%molecular_weight
    energy=old_energy+duration*external_energy
    call temperature_from_internal_energy(table,endpoint/rho,energy/rho,temperature,test_t,valid)
    if(.not.valid.or.test_t<=0.or.any(endpoint<0).or.abs(sum(endpoint)-rho)>3e-12_dp*rho) then
      reactor_audit_ok=.false.;return
    endif
    reactor_minimum_temperature=min(reactor_minimum_temperature,test_t)
    rhoy=endpoint;temperature=test_t;steps=1;ok=.true.
  end subroutine manufacture_caloric_response

  subroutine check_reactor_held_endpoint()
    real(dp) :: expected_source,pv(5+size(species)),tg,sound
    integer :: node,i,j,k,n(3)
    logical :: valid
    do node=1,size(serial%nodes)
      n=topology%nodes(node)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(geometry(node)%volume_fraction(i,j,k)<=0.or.covered(topology,node,[i,j,k])) cycle
        call recover_reactive_eb_cell_3d(species,geometry(node),[i,j,k],serial%nodes(node)%state(:,i,j,k), &
          serial%nodes(node)%temperature(i,j,k),pv,tg,sound,valid, &
          internal_energy_density=serial_h%fields(node)%values(1,i,j,k))
        call require(valid,'accepted independent actual H EOS after AMR synchronization')
      enddo
      enddo
      enddo
    enddo
    associate(g=>geometry(2))
      expected_source=-material_rate*initial_h%fields(2)%values(1,4,2,2)* &
        dt*g%volume_fraction(4,2,2)*g%dx*g%dy*g%dz
    end associate
    if(caloric_provider) expected_source=caloric_rate*dt*geometry(2)%volume_fraction(4,2,2)* &
      geometry(2)%dx*geometry(2)%dy*geometry(2)%dz
    expected_source=accepted_material_scale()*expected_source
    call require(matching([expected_source],aux(nv+1:)),'actual accepted independent UEINT source retained once')
    if(rank==0) print *, 'Audited actual reactor energy: ',mol,iterations,reactor_audits,reactor_minimum_temperature
  end subroutine check_reactor_held_endpoint
  subroutine check_held_inventory()
    type(eb_reactive_tree_3d) :: ledger
    real(dp) :: initial_inventory(1),final_inventory(1),advective_change(1),volume,expected_source,relative
    integer :: node,i,j,k,n(3)
    logical :: valid
    if(reactor_mode) then
      call check_reactor_held_endpoint()
      return
    endif
    ledger=initial
    do node=1,size(ledger%nodes)
      ledger%nodes(node)%state=initial_h%fields(node)%values
    enddo
    call integrate_reactive_eb_tree_3d(ledger,initial_inventory,valid)
    call require(valid,'initial physical-leaf UEINT inventory')
    do node=1,size(ledger%nodes)
      ledger%nodes(node)%state=serial_h%fields(node)%values
    enddo
    call integrate_reactive_eb_tree_3d(ledger,final_inventory,valid)
    call require(valid,'accepted physical-leaf UEINT inventory')
    do node=1,size(ledger%nodes)
      ledger%nodes(node)%state=0;n=topology%nodes(node)%shape
      associate(g=>geometry(node),f=>serial_faces(node))
      volume=g%dx*g%dy*g%dz
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(g%volume_fraction(i,j,k)<=0) cycle
        ledger%nodes(node)%state(1,i,j,k)=( &
          g%dy*g%dz*(g%x_face_fraction(i-1,j,k)*f%x(1,i-1,j,k)-g%x_face_fraction(i,j,k)*f%x(1,i,j,k))+ &
          g%dx*g%dz*(g%y_face_fraction(i,j-1,k)*f%y(1,i,j-1,k)-g%y_face_fraction(i,j,k)*f%y(1,i,j,k))+ &
          g%dx*g%dy*(g%z_face_fraction(i,j,k-1)*f%z(1,i,j,k-1)-g%z_face_fraction(i,j,k)*f%z(1,i,j,k)))/ &
          (g%volume_fraction(i,j,k)*volume)
      enddo
      enddo
      enddo
      end associate
    enddo
    call integrate_reactive_eb_tree_3d(ledger,advective_change,valid)
    call require(valid,'composite aperture-weighted accepted UEINT face exchange')
    associate(g=>geometry(2))
      expected_source=-3.1e13_dp*dt*g%volume_fraction(4,2,2)*g%dx*g%dy*g%dz
      if(material_source) expected_source=-material_rate*initial_h%fields(2)%values(1,4,2,2)* &
        dt*g%volume_fraction(4,2,2)*g%dx*g%dy*g%dz
    end associate
    call require(matching([expected_source],aux(nv+1:)),'accepted cold-cell source counted exactly once')
    relative=abs(final_inventory(1)-initial_inventory(1)-advective_change(1)-serial_work-aux(nv+1)-external(5)- &
      serial_projection-serial_closure)/ &
      max(abs(initial_inventory(1)),abs(final_inventory(1)),1e-30_dp)
    if(rank==0) print *, 'Independent accepted cold-source UEINT ledger:',mol,iterations,relative,aux(nv+1)
    call require(relative<=3e-12_dp,'independent accepted UEINT advection, pressure work and source balance')
  end subroutine check_held_inventory
  subroutine serial_source(state,interval,stage,depth,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: node
    if(reset_request=='seed'.and.mol.and.stage==0) then
      do node=1,size(state%nodes)
        if(depth>=0.and.topology%nodes(node)%level/=depth) cycle
        call require(all(state%nodes(node)%temperature==initial%nodes(node)%temperature.or. &
          geometry(node)%volume_fraction<=0),'serial admission preserves stored temperature seed')
      enddo
      seed_observations=seed_observations+1
      good=.false.;return
    endif
    call fill_source(rate,depth,good,stage)
  end subroutine serial_source
  subroutine native_mpi_auxiliary(state,temperature,interval,stage,depth,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: node
    logical :: seed_valid
    if(reset_request=='seed'.and.mol.and.stage==0) then
      seed_valid=.true.
      do node=1,size(temperature%fields)
        if(.not.allocated(temperature%fields(node)%values)) cycle
        if(depth>=0.and.topology%nodes(node)%level/=depth) cycle
        seed_valid=seed_valid.and.all(temperature%fields(node)%values(1,:,:,:)== &
          initial%nodes(node)%temperature.or.geometry(node)%volume_fraction<=0)
      enddo
      call require(seed_valid,'MPI admission preserves stored temperature seed')
      seed_observations=seed_observations+1
      good=.false.;return
    endif
    call fill_source(rate,depth,good,stage)
  end subroutine native_mpi_auxiliary
  subroutine check_seed_admission()
    integer :: node
    logical :: rollback_valid
    ! Observe the admitted old clock, then refuse the source deliberately.
    ! This isolates admission and transactional rollback from chemistry.
    mol=.true.;seed_observations=0
    serial=initial;serial_h=initial_h;serial_prefix=prefix
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,1,.true., &
      source_transport_flags,forcing,1e-8_dp,1e-12_dp,serial,mol_extended,faces, &
      boundary,wall,chemical,external,theta,good,mol_mode=.true., &
      internal_energy=serial_h,internal_energy_flux=serial_faces,energy_policy=policy,energy_correction=serial_prefix, &
      auxiliary_source=serial_source,auxiliary_change=aux,auxiliary_has_internal_energy=.true.)
    call require(.not.good.and.seed_observations>0,'serial seed observed before intentional refusal')
    do node=1,2
      call require(all(serial%nodes(node)%state==initial%nodes(node)%state).and. &
        all(serial%nodes(node)%temperature==initial%nodes(node)%temperature).and. &
        all(serial_h%fields(node)%values==initial_h%fields(node)%values),'serial seed refusal rollback')
    enddo
    call require(serial_prefix==prefix,'serial seed refusal ledger rollback')
    call localize(initial,initial_h,gas,temp,held)
    mpi_prefix=prefix;seed_observations=0
    call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,gas,temp,mol_extended, &
      dt,1,.true.,source_transport_flags,forcing,1e-8_dp,1e-12_dp,faces,mb,mw,mc,me,mt,MPI_COMM_WORLD,good, &
      mol_mode=.true.,internal_energy=held,internal_energy_flux=mpi_faces, &
      energy_policy=policy,energy_correction=mpi_prefix, &
      auxiliary_source=native_mpi_auxiliary,auxiliary_change=maux,auxiliary_has_internal_energy=.true.)
    call require(.not.good.and.seed_observations>0,'MPI seed observed before intentional refusal')
    rollback_valid=.true.
    do node=1,2
      if(owners(node)/=rank) cycle
      rollback_valid=rollback_valid.and.all(gas%fields(node)%values==initial%nodes(node)%state).and. &
        all(temp%fields(node)%values(1,:,:,:)==initial%nodes(node)%temperature).and. &
        all(held%fields(node)%values==initial_h%fields(node)%values)
    enddo
    call require(rollback_valid,'MPI seed refusal rollback')
    call require(mpi_prefix==prefix,'MPI seed refusal ledger rollback')
    if(rank==0) print *, 'PASS: serial/MPI temperature seed admission and refusal rollback'
  end subroutine check_seed_admission
  subroutine fill_source(rate,depth,good,stage)
    type(amr_tree_state_3d), intent(inout) :: rate
    integer, intent(in) :: depth,stage
    logical, intent(out) :: good
    integer :: p
    real(dp) :: rate_scale
    rate_scale=1
    if(corrector_source) then
      rate_scale=0
      if(stage>0) rate_scale=2
    endif
    call require(.not.source_forbidden,'frozen source replay performs no auxiliary callback')
    auxiliary_evaluations=auxiliary_evaluations+1
    good=.true.
    do p=1,size(rate%fields)
      if(.not.allocated(rate%fields(p)%values)) cycle
      rate%fields(p)%values=0
      if(requested_mode/='eos_source') cycle
      if(depth>=0.and.topology%nodes(p)%level/=depth) cycle
      if(p==1.and.covered_source) then
        ! This virtual coarse cell owns no composite source inventory.
        rate%fields(p)%values(:nv,2,1,1)=-rate_scale*material_rate*cell
        rate%fields(p)%values(nv+1,2,1,1)=-rate_scale*material_rate*initial_h%fields(1)%values(1,2,1,1)
      endif
      if(p/=2) cycle
      if(material_source) then
        rate%fields(p)%values(:nv,4,2,2)=-rate_scale*material_rate*cell
        rate%fields(p)%values(nv+1,4,2,2)=-rate_scale*material_rate*initial_h%fields(2)%values(1,4,2,2)
      else if(caloric_provider) then
        rate%fields(p)%values(5,4,2,2)=rate_scale*caloric_rate
        rate%fields(p)%values(nv+1,4,2,2)=rate_scale*caloric_rate
      else
        rate%fields(p)%values(nv+1,4,2,2)=-3.1e13_dp
      endif
    enddo
  end subroutine fill_source
  real(dp) function accepted_material_scale() result(scale)
    scale=1
    if(reactor_mode.and.mol.and.iterations>1) then
      ! Later pinned MOL passes retain only the old half of volume forcing.
      scale=.5_dp
      if(corrector_source) scale=0
    endif
  end function accepted_material_scale
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
end program test_mpi_eb_held_reset_ledger_3d
