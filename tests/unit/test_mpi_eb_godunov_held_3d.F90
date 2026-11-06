program test_mpi_eb_godunov_held_3d
  use iso_fortran_env, only: int64,output_unit
  use mpi_f08
  use precision_mod, only: dp
  use dual_energy_mod, only: dual_energy_policy
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use gas_transport_mod, only: gas_transport_species
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_3d,advance_mpi_eb_tree_sdc_subcycle_3d
  use transport_database_mod, only: load_h2o2_elementary_transport
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    build_eb_tree_extended_geometry_3d,compute_sdc_eb_tree_hydro_3d,compute_godunov_state_redist_eb_tree_hydro_3d, &
    eb_godunov_interval_3d,eb_godunov_source_stage_3d, &
    advance_reactive_eb_tree_sdc_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_sdc_subcycle_3d
  use mpi_eb_tree_hydro_3d_mod, only: compute_sdc_mpi_eb_tree_hydro_3d,compute_godunov_state_redist_mpi_eb_tree_hydro_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_redistribution_3d_mod, only: state_redist_register_3d,accumulate_state_redist_register_3d
  implicit none
  logical :: scalar_clock_reset_test=.false.
  integer, parameter :: n=4,np=5
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,padded
  type(eb_geometry_3d), allocatable :: geometry(:),extended(:)
  type(eb_reactive_tree_3d) :: tree,ancestor_end
  type(amr_tree_state_3d) :: held,held_end,predictor,predictor_end,state,temperature,local_held,local_predictor
  type(amr_tree_state_3d) :: end_state,end_temperature,local_held_end,local_predictor_end
  type(amr_tree_state_3d) :: gas_rate,h_rate,gas_raw,gas_old,h_raw,h_old,pwork
  type(amr_tree_state_3d) :: ref_gas,ref_h,ref_gas_raw,ref_gas_old,ref_h_raw,ref_h_old,ref_pwork
  type(eb_integrated_flux_3d), allocatable :: gas_flux(:),h_flux(:),gas_base(:),h_base(:)
  type(eb_integrated_flux_3d), allocatable :: ref_gas_flux(:),ref_h_flux(:),ref_gas_base(:),ref_h_base(:)
  type(state_redist_register_3d), allocatable :: gas_registers(:),h_registers(:),ref_gas_registers(:),ref_h_registers(:)
  real(dp), allocatable :: initial(:,:,:,:),initial_t(:,:,:),q(:),selected_weight
  real(dp), parameter :: test_weights(3)=[1._dp,0._dp,.375_dp]
  type(amr_tree_state_3d) :: native_pwork
  real(dp) :: sound,nan_value,position,inventory,boundary,max_residual
  real(dp) :: fractions(2)=[.35_dp,.5_dp]
  character(len=24) :: faces(6)='outflow'
  integer :: nv,rank,nranks,status,p,i,j,k,mode,layout,stage,selection,owners(np),refusal_mode,failure,only(3),weight_case
  character(len=32) :: argument
  logical :: good,same,changed,common_test,common_transport,recursive_test,weight_test,common_react,history_test,clock_test
  logical :: prefix_replay_test=.false.,prefix_forbid_parent=.false.,recovery_replay_test=.false.,recovery_fixture=.false.
  logical :: subtree_replay_test=.false.,routed_fixture=.false.,family_replay_test,source_replay_test,provider_must_fail=.false.
  logical :: source_stages_test=.false.,stage_sink=.false.,stage_parent_test=.false.
  logical :: packet_test=.false.,raw_refusal=.false.,raw_malformed=.false.
  integer :: stage_sink_depth=-1
  integer :: aux_calls=0,raw_calls=0,reactor_calls=0
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  source_replay_test=.false.;family_replay_test=.false.
  only=-1;common_test=.false.;common_transport=.false.;recursive_test=.false.
  weight_test=.false.;common_react=.false.;history_test=.false.;clock_test=.false.
  if(command_argument_count()==4) then
    call get_command_argument(1,argument)
    scalar_clock_reset_test=argument=='recursive_reset'
    packet_test=argument=='source_packets'
    stage_parent_test=argument=='stage_parent'
    source_stages_test=argument=='source_stages'
    recovery_replay_test=argument=='recovery_replay'
    prefix_replay_test=argument=='prefix_replay'
    subtree_replay_test=argument=='subtree_replay'.or.prefix_replay_test
    routed_fixture=argument=='routed_family_replay'.or.subtree_replay_test
    family_replay_test=argument=='family_replay'.or.routed_fixture
    source_replay_test=argument=='source_replay'.or.family_replay_test.or.recovery_replay_test.or. &
      source_stages_test.or.stage_parent_test.or.packet_test
    routed_fixture=routed_fixture.or.recovery_replay_test.or.stage_parent_test.or.packet_test
    clock_test=argument=='clocks'.or.argument=='clocks_diff'
    history_test=argument=='history'.or.argument=='history_diff'
    common_react=argument=='common_weights_react'.or.history_test
    weight_test=argument=='weights'.or.argument=='common_weights'.or.common_react
    recursive_test=argument=='recursive'.or.argument=='recursive_diff'.or.clock_test.or.scalar_clock_reset_test
    common_transport=argument=='common_diff'.or.argument=='recursive_diff'.or. &
      argument=='history_diff'.or.argument=='clocks_diff'
    common_test=argument=='common'.or.common_transport.or.recursive_test.or. &
      argument=='common_weights'.or.common_react
    call require(common_test.or.weight_test.or.source_replay_test,'supported diagnostic family')
    do i=1,3
      call get_command_argument(i+1,argument);read(argument,*) only(i)
    enddo
  endif
  if(command_argument_count()==3) then
    do i=1,3
      call get_command_argument(i,argument);read(argument,*) only(i)
    enddo
  endif
  call require(nranks==3.or.(nranks==1.and.all(only>=0)), &
    'three-rank fixture or single-rank targeted diagnostic')
  call load_h2o2_elementary_thermo(species,good)
  call require(good,'thermo')
  nv=reactive_nvar(size(species));nan_value=ieee_value(0._dp,ieee_quiet_nan)
  max_residual=0;refusal_mode=0
  allocate(initial(nv,n,n,n),initial_t(n,n,n),q(nv))
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,2,3,2,boxes(1),good)
  call require(good,'child left')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,2,3,2,boxes(2),good)
  call require(good,'child right')
  call initialize_amr_patch_3d(4,4,4,3,4,3,4,2,3,2,boxes(3),good)
  call require(good,'grandchild left')
  call initialize_amr_patch_3d(4,4,4,1,2,3,4,2,3,2,boxes(4),good)
  call require(good,'grandchild right')
  call initialize_amr_patch_tree_3d([n,n,n],[1,1,2,3],boxes,topology,good)
  call require(good,'three-level peer topology')
  do k=1,n
  do j=1,n
  do i=1,n
    q=0;q(1)=.8_dp+.001_dp*(i+j+k);q(2:4)=[4._dp,2._dp,-1._dp]-.1_dp*[i,j,k]
    q(5)=500000+100*j
    do p=1,size(species)
      if(species(p)%name=='N2') q(p+5)=1
    enddo
    call reactive_primitive_to_conserved(species,q,initial(:,i,j,k),initial_t(i,j,k),sound,good)
    call require(good,'initial gas')
  enddo
  enddo
  enddo
  do mode=0,1
    if(only(1)>=0.and.mode/=only(1)) cycle
    position=merge(-4._dp,.46_dp,mode==0)
    call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'z',position,finest,good)
    call require(good,'regular/cut finest geometry')
    call initialize_reactive_eb_tree_3d(species,topology,finest,initial,initial_t,tree,good,internal_energy=held)
    call require(good,'actual held tree initialization')
    if(allocated(geometry)) deallocate(geometry)
    allocate(geometry(np))
    do p=1,np
      geometry(p)=tree%nodes(p)%geometry
      do k=1,n
      do j=1,n
      do i=1,n
        if(geometry(p)%volume_fraction(i,j,k)<=0) cycle
        held%fields(p)%values(1,i,j,k)=held%fields(p)%values(1,i,j,k)+ &
          137*(geometry(p)%x_lower+(i-.5_dp)*geometry(p)%dx)+ &
          41*(geometry(p)%y_lower+(j-.5_dp)*geometry(p)%dy)+23*(geometry(p)%z_lower+(k-.5_dp)*geometry(p)%dz)
      enddo
      enddo
      enddo
    enddo
    call build_axis_plane_eb_geometry_3d(96,96,96,-2.5_dp,3.5_dp,-2.5_dp,3.5_dp,-2.5_dp,3.5_dp, &
      'z',position,padded,good)
    call require(good,'ten-layer finest context')
    call build_eb_tree_extended_geometry_3d(topology,geometry,padded,extended,good,halo_width=10)
    call require(good,'ten-layer patch context')
    if(allocated(predictor%fields)) deallocate(predictor%fields)
    allocate(predictor%fields(np))
    do p=1,np
      allocate(predictor%fields(p)%values(nv,n,n,n),source=0._dp)
      predictor%fields(p)%values(2,:,:,:)=.3_dp*p
      predictor%fields(p)%values(5,:,:,:)=10000._dp*p
    enddo
    ancestor_end=tree;held_end=held;predictor_end=predictor
    do p=1,np
      ancestor_end%nodes(p)%state=1.01_dp*tree%nodes(p)%state
      held_end%fields(p)%values=1.03_dp*held%fields(p)%values
      predictor_end%fields(p)%values=2*predictor%fields(p)%values
    enddo
    if(source_replay_test) then
      if(packet_test) then
        call check_source_packets()
        call check_packet_parent()
      else if(stage_parent_test) then
        call check_stage_parent()
      else if(source_stages_test) then
        call check_source_stages()
      else if(recovery_replay_test) then
        call check_recovery_replay()
      else if(subtree_replay_test) then
        call check_subtree_replay()
      else if(family_replay_test) then
        call check_family_replay()
      else
        call check_source_replay()
      endif
      cycle
    endif
    do selection=0,1
      if(only(2)>=0.and.selection/=only(2)) cycle
      if(selection==1.and.mode==0) cycle
      if(common_test) then
        do stage=0,1
          if(only(3)>=0.and.stage/=only(3)) cycle
          do weight_case=merge(3,merge(0,1,weight_test),history_test),merge(3,1,weight_test)
            if(weight_test.and.weight_case==0) then
              if(allocated(selected_weight)) deallocate(selected_weight)
            endif
            if(weight_test.and.weight_case>0) then
              if(.not.allocated(selected_weight)) allocate(selected_weight)
              selected_weight=test_weights(weight_case)
              if(rank==0) print *, 'SDC paired face weight:',selected_weight,' native reactor/reset:',common_react
            endif
            call check_common_step()
          enddo
        enddo
        cycle
      endif
      do stage=0,1
        if(only(3)>=0.and.stage/=only(3)) cycle
        do weight_case=1,merge(3,1,weight_test)
        if(weight_test) then
          if(.not.allocated(selected_weight)) allocate(selected_weight)
          selected_weight=test_weights(weight_case)
          if(rank==0) print *, 'paired face weight:',selected_weight
        endif
        if(rank==0) print *, 'Godunov held source case (geometry, selected, StateRedist):',mode,selection,stage
        call evaluate_serial()
        call require(good,'serial actual held Godunov raw/StateRedist source')
        if(weight_test) then
          if(weight_case==1) native_pwork=ref_pwork
          call require(equal_state(ref_pwork,native_pwork),'native pressure work independent of face weight')
          if(stage==0.and.selected_weight==0) then
            do p=1,np
              call require(equal_flux(ref_gas_flux(p),ref_gas_base(p)).and.equal_flux(ref_h_flux(p),ref_h_base(p)), &
                'raw zero weight retains exact gas and UEINT PCM faces, including inactive levels')
            enddo
          endif
        endif
        if(stage==0) then
          do p=1,np
            if(.not.allocated(ref_h_old%fields(p)%values)) cycle
            call require(all(pack(ref_h_old%fields(p)%values(1,7:10,7:10,7:10),geometry(p)%volume_fraction>0)== &
              pack(held%fields(p)%values(1,:,:,:),geometry(p)%volume_fraction>0)), &
              'raw old halo retains actual physical held UEINT')
            inventory=1e-8_dp*geometry(p)%dx*geometry(p)%dy*geometry(p)%dz* &
              sum(geometry(p)%volume_fraction*(ref_h%fields(p)%values(1,:,:,:)+ref_pwork%fields(p)%values(1,:,:,:)))
            boundary=geometry(p)%dy*geometry(p)%dz*sum( &
              geometry(p)%x_face_fraction(0,:,:)*ref_h_flux(p)%x(1,0,:,:)- &
              geometry(p)%x_face_fraction(n,:,:)*ref_h_flux(p)%x(1,n,:,:))+ &
              geometry(p)%dx*geometry(p)%dz*sum(geometry(p)%y_face_fraction(:,0,:)*ref_h_flux(p)%y(1,:,0,:)- &
              geometry(p)%y_face_fraction(:,n,:)*ref_h_flux(p)%y(1,:,n,:))+ &
              geometry(p)%dx*geometry(p)%dy*sum(geometry(p)%z_face_fraction(:,:,0)*ref_h_flux(p)%z(1,:,:,0)- &
              geometry(p)%z_face_fraction(:,:,n)*ref_h_flux(p)%z(1,:,:,n))
            max_residual=max(max_residual,abs(inventory-boundary))
            call require(abs(inventory-boundary)<1e-15_dp,'independent per-patch raw held advection plus pressure work')
          enddo
        endif
        do layout=1,merge(2,0,nranks==3)
          owners=[0,0,1,0,1]
          if(layout==2) owners=[1,0,1,1,0]
          call distribute()
          call evaluate_mpi()
          call require(good,'MPI actual held Godunov raw/StateRedist source')
          same=.true.
          do p=1,np
            same=same.and.(allocated(h_rate%fields(p)%values).eqv.(owners(p)==rank))
            if(owners(p)/=rank) cycle
            same=same.and.all(gas_rate%fields(p)%values==ref_gas%fields(p)%values)
            same=same.and.all(h_rate%fields(p)%values==ref_h%fields(p)%values)
            same=same.and.all(pwork%fields(p)%values==ref_pwork%fields(p)%values)
            if(selection==1.and.topology%nodes(p)%level/=1) cycle
            same=same.and.equal_flux(gas_flux(p),ref_gas_flux(p)).and.equal_flux(h_flux(p),ref_h_flux(p))
            same=same.and.equal_flux(gas_base(p),ref_gas_base(p)).and.equal_flux(h_base(p),ref_h_base(p))
            same=same.and.all(gas_raw%fields(p)%values==ref_gas_raw%fields(p)%values)
            same=same.and.all(gas_old%fields(p)%values==ref_gas_old%fields(p)%values)
            same=same.and.all(h_raw%fields(p)%values==ref_h_raw%fields(p)%values)
            same=same.and.all(h_old%fields(p)%values==ref_h_old%fields(p)%values)
            if(stage==1) same=same.and.equal_register(gas_registers(p),ref_gas_registers(p)).and. &
              equal_register(h_registers(p),ref_h_registers(p))
          enddo
          call require(same,'exact serial/MPI gas and held fields, halos, faces and registered transfers')
          if(stage==1.and.mode==1) then
            changed=.false.
            do p=1,np
              if(owners(p)/=rank.or..not.allocated(h_raw%fields(p)%values)) cycle
              changed=changed.or.any(h_rate%fields(p)%values/=h_raw%fields(p)%values(:,7:10,7:10,7:10))
            enddo
            call MPI_Allreduce(MPI_IN_PLACE,changed,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,status)
            call require(changed,'cut-cell StateRedist changes independent UEINT rate')
          endif
          if(layout==2.and.mode==1.and.selection==0.and.stage==1) call check_refusal()
          if(weight_test.and.layout==2.and.weight_case==3) call check_weight_refusal()
        enddo
        enddo
      enddo
    enddo
  enddo
  if(rank==0) then
    if(common_test) then
      print *, 'PASS: actual-held Godunov SDC, recursive/native transport:',recursive_test,common_transport, &
        ' exact fields/faces/registers, bounded ledgers, owner/idle rollback; residual=',max_residual
    else
      print *, 'PASS: actual-held Godunov raw/StateRedist source, ancestors, peers, owner/idle rollback; residual=', &
        max_residual
    endif
  endif
  call MPI_Finalize(status)
contains
  subroutine check_source_packets()
    type(eb_reactive_tree_3d) :: result,expected
    type(amr_tree_state_3d) :: result_h,expected_h,saved_gas,saved_temp,saved_h
    type(eb_godunov_interval_3d) :: record,recovery,cache,local_record,local_recovery,bad
    type(eb_integrated_flux_3d), allocatable :: f(:),h(:),lf(:),lh(:),ef(:),eh(:)
    type(state_redist_register_3d), allocatable :: gr(:),hr(:),lgr(:),lhr(:)
    real(dp) :: b(nv),w(nv),c(nv),v(nv),a(nv+1),theta,work,prefix
    integer :: node,invalid
    logical :: accepted,refused,ready,equal
    owners=rank;result=tree;result_h=held;prefix=7.5_dp;work=9.5_dp
    aux_calls=0;raw_calls=0;provider_must_fail=.false.;stage_sink=.true.;stage_sink_depth=-1
    call serial_source_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted, &
      recovery=recovery,recovery_ready=ready,refused=refused)
    call require(.not.accepted.and.ready.and..not.refused,'native density post refuses after both packets accepted')
    call require(recovery%source_prefix_only.and.recovery%accepted_source_stages==0.and. &
      recovery%sources(0)%auxiliary_accepted.and.recovery%sources(0)%provider_accepted, &
      'packet frontier retained before any complete source snapshot')
    call require(aux_calls==1.and.raw_calls==1,'both accepted declarations called once')
    equal=prefix==7.5_dp.and.work==9.5_dp.and.record%accepted_source_stages==0
    do node=1,np
      equal=equal.and.all(result%nodes(node)%state==tree%nodes(node)%state).and. &
        all(result_h%fields(node)%values==held%fields(node)%values)
    enddo
    call require(equal.and..not.allocated(f).and..not.allocated(h),'native source failure publishes no physical result')
    cache=recovery;stage_sink=.false.;provider_must_fail=.true.;aux_calls=0;raw_calls=0
    call serial_source_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted, &
      replay=cache,partial=.true.,recovery=recovery,recovery_ready=ready,refused=refused)
    call require(.not.accepted.and.ready.and..not.refused.and.aux_calls==0.and.raw_calls==0, &
      'native retry freezes accepted packets and retains its density-sum refusal')
    provider_must_fail=.false.
    do invalid=1,2
      result=tree;result_h=held;raw_refusal=invalid==1;raw_malformed=invalid==2;aux_calls=0;raw_calls=0
      call serial_source_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted, &
        recovery=bad,recovery_ready=ready,refused=refused)
      call require(.not.accepted.and..not.ready.and.aux_calls==1.and.raw_calls==1, &
        'refused or malformed provider after accepted auxiliary blocks retry-cache publication')
    enddo
    raw_refusal=.false.;raw_malformed=.false.
    bad=cache
    bad%sources(0)%provider_accepted=.false.;bad%sources(0)%provider_kind=0
    deallocate(bad%sources(0)%provider_rate%fields,bad%sources(0)%provider_halos%fields,bad%sources(0)%provider_faces)
    result=tree;result_h=held;prefix=.125_dp;aux_calls=0;raw_calls=0
    call serial_source_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted,replay=bad,partial=.true.)
    call require(accepted.and.aux_calls==2.and.raw_calls==3,'auxiliary-only frontier first-evaluates its unvisited provider')
    expected=result;expected_h=result_h;ef=f;eh=h
    if(nranks==1) return
    owners=[1,0,1,0,1];call distribute();stage_sink=.true.;prefix=7.5_dp;work=9.5_dp;aux_calls=0;raw_calls=0
    saved_gas=state;saved_temp=temperature;saved_h=local_held
    call mpi_source_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused, &
      recovery=local_recovery,recovery_ready=ready)
    call require(.not.accepted.and.ready.and..not.refused.and.aux_calls==1.and.raw_calls==1, &
      'MPI accepts both declarations before native density post refuses')
    call require(local_recovery%accepted_source_stages==0.and.local_recovery%sources(0)%auxiliary_accepted.and. &
      local_recovery%sources(0)%provider_accepted,'MPI count0 packet frontier on two owners and idle rank')
    call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
      equal_state(local_held,saved_h).and.prefix==7.5_dp.and.work==9.5_dp,'MPI source failure full rollback')
    cache=local_recovery;stage_sink=.false.;provider_must_fail=.true.;aux_calls=0;raw_calls=0
    call mpi_source_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused, &
      replay=cache,partial=.true.,recovery=local_recovery,recovery_ready=ready)
    call require(.not.accepted.and.ready.and..not.refused.and.aux_calls==0.and.raw_calls==0, &
      'MPI native replay freezes accepted declarations without suppressing native failure')
    provider_must_fail=.false.
    do invalid=1,4
      bad=cache
      if(invalid==1.and.rank==2) bad%accepted_source_stages=1
      if(invalid==2.and.rank==owners(2)) bad%sources(0)%auxiliary%fields(2)%values(nv+1,1,1,1)=nan_value
      if(invalid==3.and.rank==owners(2)) bad%sources(0)%provider_rate%fields(2)%values(2,1,1,1)=nan_value
      if(invalid==4.and.rank==2) then
        bad%sources(0)%provider_accepted=.false.;bad%sources(0)%provider_kind=0
        deallocate(bad%sources(0)%provider_rate%fields,bad%sources(0)%provider_halos%fields,bad%sources(0)%provider_faces)
      endif
      aux_calls=0;raw_calls=0
      call mpi_source_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused, &
        replay=bad,partial=.true.,recovery=local_recovery,recovery_ready=ready)
      call require(.not.accepted.and..not.ready.and.aux_calls==0.and.raw_calls==0, &
        'invalid/mismatched packet frontier refuses collectively before choosing callbacks')
      call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
        equal_state(local_held,saved_h).and.prefix==7.5_dp.and.work==9.5_dp,'bad packet cache full rollback')
    enddo
    bad=cache
    bad%sources(0)%provider_accepted=.false.;bad%sources(0)%provider_kind=0
    deallocate(bad%sources(0)%provider_rate%fields,bad%sources(0)%provider_halos%fields,bad%sources(0)%provider_faces)
    prefix=.125_dp;aux_calls=0;raw_calls=0
    call mpi_source_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused,replay=bad,partial=.true.)
    call require(accepted.and.aux_calls==2.and.raw_calls==3,'MPI auxiliary-only frontier first-evaluates raw provider')
    equal=.true.
    do node=1,np
      if(owners(node)/=rank) cycle
      equal=equal.and.all(state%fields(node)%values==expected%nodes(node)%state).and. &
        all(local_held%fields(node)%values==expected_h%fields(node)%values).and. &
        equal_flux(lf(node),ef(node)).and.equal_flux(lh(node),eh(node))
    enddo
    call require(equal.and.record_equal(local_record,record),'exact serial/MPI independent auxiliary/provider native completion')
    if(rank==0) print *, 'PASS: independent accepted packet frontier, native density refusal and first provider evaluation'
  end subroutine check_source_packets
  subroutine check_packet_parent()
    type(eb_reactive_tree_3d) :: result,fresh,expected
    type(amr_tree_state_3d) :: result_h,fresh_h,expected_h,saved_gas,saved_temp,saved_h,expected_impulse
    type(eb_godunov_interval_3d), allocatable :: reference(:),record(:),recovery(:),cache(:),local_record(:),local_cache(:),bad(:)
    type(eb_integrated_flux_3d), allocatable :: f(:),h(:),lf(:),lh(:),ef(:),eh(:)
    real(dp) :: b(nv),w(nv),c(nv),v(nv),a(nv+1),theta,work,prefix
    integer(int64) :: counts(3)
    integer :: node,r,invalid
    logical :: accepted,refused,ready,equal
    owners=rank;result=tree;result_h=held;prefix=.125_dp;stage_sink=.false.;provider_must_fail=.false.
    call serial_family_step(result,result_h,reference,b,w,c,v,a,theta,work,prefix,f,h,accepted)
    call require(accepted,'native routed reference with actual old/new endpoint packets')
    result=tree;result_h=held;prefix=7.5_dp;work=9.5_dp;counts=-19
    stage_sink=.true.;stage_sink_depth=1;aux_calls=0;raw_calls=0
    call serial_subtree_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,accepted, &
      recovery=recovery,recovery_ready=ready,refused=refused,counts=counts)
    call require(.not.accepted.and.ready.and..not.refused,'automatic parent escalation returns actual child-stage prefix')
    call require(allocated(recovery),'failed parent-stage source bank present')
    call require(size(recovery)==2,'root complete plus first child source0 are retained')
    call require(aux_calls==4.and.raw_calls==4,'accepted providers never repeat through parent/root attempts')
    call require(recovery(2)%accepted_source_stages==0.and.recovery(2)%source_prefix_only,'native child stage0 frontier')
    equal=prefix==7.5_dp.and.work==9.5_dp.and.all(counts==-19)
    do node=1,np
      equal=equal.and.all(result%nodes(node)%state==tree%nodes(node)%state).and. &
        all(result_h%fields(node)%values==held%fields(node)%values)
    enddo
    call require(equal.and..not.allocated(record),'failed automatic repair retains outer state/history/counts/ledgers')
    call require(equal_state(recovery(2)%routed_auxiliary_old_rate,reference(2)%routed_auxiliary_old_rate), &
      'child off-level old packet survives failed common step without inference')
    call require(recovery(2)%sources(0)%auxiliary_accepted.and.recovery(2)%sources(0)%provider_accepted, &
      'child native source failure retains both accepted declarations before completing stage0')
    cache=recovery;stage_sink=.false.;stage_sink_depth=-1;provider_must_fail=.true.;aux_calls=0;raw_calls=0
    call serial_subtree_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,accepted, &
      replay=cache,replay_prefix=.true.,recovery=recovery,recovery_ready=ready,refused=refused,counts=counts)
    call require(.not.accepted.and.ready.and..not.refused.and.aux_calls==0.and.raw_calls==0, &
      'recursive source failure reruns native operators without repeating accepted declarations')
    call require(recovery(2)%accepted_source_stages==0.and.prefix==7.5_dp.and.work==9.5_dp.and.all(counts==-19), &
      'repeated source failure preserves stage frontier and caller ledgers')
    provider_must_fail=.false.
    if(nranks==1) return
    owners=[1,0,1,0,1];call distribute();prefix=7.5_dp;work=9.5_dp;counts=-19
    stage_sink=.true.;stage_sink_depth=1;aux_calls=0;raw_calls=0
    saved_gas=state;saved_temp=temperature;saved_h=local_held
    call mpi_subtree_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused, &
      recovery=recovery,recovery_ready=ready,counts=counts)
    call require(.not.accepted.and.ready.and..not.refused,'MPI automatic child/parent failures preserve accepted source frontier')
    call require(size(recovery)==2.and.aux_calls==4.and.raw_calls==4,'MPI parent attempts do not repeat accepted stages')
    call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
      equal_state(local_held,saved_h).and.prefix==7.5_dp.and.work==9.5_dp.and.all(counts==-19), &
      'MPI failed automatic stage recovery retains physical state and ledgers')
    call require(recovery(2)%accepted_source_stages==0.and.recovery(2)%sources(0)%auxiliary_accepted.and. &
      recovery(2)%sources(0)%provider_accepted,'MPI native child source frontier agrees across owners and idle rank')
    local_cache=recovery;stage_sink=.false.;stage_sink_depth=-1;provider_must_fail=.true.;aux_calls=0;raw_calls=0
    call mpi_subtree_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused, &
      replay=local_cache,replay_prefix=.true.,recovery=recovery,recovery_ready=ready,counts=counts)
    call require(.not.accepted.and.ready.and..not.refused.and.aux_calls==0.and.raw_calls==0, &
      'MPI recursive native source failure freezes every previously accepted component')
    call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
      equal_state(local_held,saved_h).and.prefix==7.5_dp.and.work==9.5_dp.and.all(counts==-19), &
      'MPI repeated source failure full rollback')
    provider_must_fail=.false.
    saved_gas=state;saved_temp=temperature;saved_h=local_held
    do invalid=1,2
      bad=local_cache
      if(invalid==1.and.rank==owners(2)) bad(2)%routed_auxiliary_old_rate%fields(2)%values(nv+1,1,1,1)=nan_value
      if(invalid==2.and.rank==2) bad(2)%accepted_source_stages=1
      aux_calls=0;raw_calls=0;prefix=7.5_dp;work=9.5_dp
      call mpi_subtree_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused, &
        replay=bad,replay_prefix=.true.,counts=counts)
      call require(.not.accepted.and.aux_calls==0.and.raw_calls==0,'malformed partial/routed cache rejects before providers')
      call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
        equal_state(local_held,saved_h).and.prefix==7.5_dp.and.work==9.5_dp,'malformed partial clock full rollback')
    enddo
    if(rank==0) print *, 'PASS: independent source packets through recursive parent recovery, native refusal and MPI rollback'
  end subroutine check_packet_parent
  subroutine check_stage_parent()
    type(eb_reactive_tree_3d) :: result,fresh,expected
    type(amr_tree_state_3d) :: result_h,fresh_h,expected_h,saved_gas,saved_temp,saved_h,expected_impulse
    type(eb_godunov_interval_3d), allocatable :: reference(:),record(:),recovery(:),cache(:),local_record(:),local_cache(:),bad(:)
    type(eb_integrated_flux_3d), allocatable :: f(:),h(:),lf(:),lh(:),ef(:),eh(:)
    real(dp) :: b(nv),w(nv),c(nv),v(nv),a(nv+1),theta,work,prefix
    integer(int64) :: counts(3)
    integer :: node,r,invalid
    logical :: accepted,refused,ready,equal
    owners=rank;result=tree;result_h=held;prefix=.125_dp;stage_sink=.false.;provider_must_fail=.false.
    call serial_family_step(result,result_h,reference,b,w,c,v,a,theta,work,prefix,f,h,accepted)
    call require(accepted,'native routed reference with actual old/new endpoint packets')
    result=tree;result_h=held;prefix=7.5_dp;work=9.5_dp;counts=-19
    stage_sink=.true.;stage_sink_depth=1;aux_calls=0;raw_calls=0
    call serial_subtree_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,accepted, &
      recovery=recovery,recovery_ready=ready,refused=refused,counts=counts)
    call require(.not.accepted.and.ready.and..not.refused,'automatic parent escalation returns actual child-stage prefix')
    call require(allocated(recovery),'failed parent-stage source bank present')
    call require(size(recovery)==2,'root complete plus first child source0 are retained')
    call require(aux_calls==4.and.raw_calls==4,'accepted providers never repeat through parent/root attempts')
    call require(recovery(2)%accepted_source_stages==1.and.recovery(2)%source_prefix_only,'native child stage0 frontier')
    equal=prefix==7.5_dp.and.work==9.5_dp.and.all(counts==-19)
    do node=1,np
      equal=equal.and.all(result%nodes(node)%state==tree%nodes(node)%state).and. &
        all(result_h%fields(node)%values==held%fields(node)%values)
    enddo
    call require(equal.and..not.allocated(record),'failed automatic repair retains outer state/history/counts/ledgers')
    call require(equal_state(recovery(2)%routed_auxiliary_old_rate,reference(2)%routed_auxiliary_old_rate), &
      'child off-level old packet survives failed common step without inference')
    cache=recovery;stage_sink=.false.;stage_sink_depth=-1;fresh=tree;fresh_h=held
    do node=1,np
      fresh%nodes(node)%state=1e9_dp*fresh%nodes(node)%state
      fresh_h%fields(node)%values=1e9_dp*fresh_h%fields(node)%values
    enddo
    result=fresh;result_h=fresh_h;prefix=.125_dp;aux_calls=0;raw_calls=0
    call serial_subtree_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,accepted, &
      replay=cache,replay_prefix=.true.,counts=counts)
    call require(accepted,'native whole-family replay completes partial child source table on fresh state')
    call require(all(counts==[1_int64,2_int64,4_int64]).and.aux_calls==17.and.raw_calls==17, &
      'only seventeen previously unvisited source stages run')
    call require(size(record)==7.and.all(record%accepted_source_stages==3).and..not.any(record%source_prefix_only), &
      'complete source/clock family returns after partial replay')
    call require(equal_state(record(2)%routed_auxiliary_old_rate,cache(2)%routed_auxiliary_old_rate), &
      'old off-level packet stays frozen while new endpoint is evaluated')
    expected_impulse=cache(2)%routed_auxiliary_old_rate
    do node=1,np
      expected_impulse%fields(node)%values=record(2)%interval*expected_impulse%fields(node)%values
    enddo
    call require(equal_state(record(2)%routed_auxiliary_impulse,expected_impulse), &
      'actual partial old/new off-level packets integrate to independent constant-rate impulse')
    expected=result;expected_h=result_h;ef=f;eh=h
    if(nranks==1) return
    owners=[1,0,1,0,1];call distribute();prefix=7.5_dp;work=9.5_dp;counts=-19
    stage_sink=.true.;stage_sink_depth=1;aux_calls=0;raw_calls=0
    saved_gas=state;saved_temp=temperature;saved_h=local_held
    call mpi_subtree_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused, &
      recovery=recovery,recovery_ready=ready,counts=counts)
    call require(.not.accepted.and.ready.and..not.refused,'MPI automatic child/parent failures preserve accepted source frontier')
    call require(size(recovery)==2.and.aux_calls==4.and.raw_calls==4,'MPI parent attempts do not repeat accepted stages')
    call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
      equal_state(local_held,saved_h).and.prefix==7.5_dp.and.work==9.5_dp.and.all(counts==-19), &
      'MPI failed automatic stage recovery retains physical state and ledgers')
    local_cache=recovery;stage_sink=.false.;stage_sink_depth=-1;prefix=.125_dp
    do node=1,np
      if(owners(node)/=rank) cycle
      state%fields(node)%values=fresh%nodes(node)%state
      local_held%fields(node)%values=fresh_h%fields(node)%values
    enddo
    aux_calls=0;raw_calls=0
    call mpi_subtree_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused, &
      replay=local_cache,replay_prefix=.true.,counts=counts)
    call require(accepted.and.aux_calls==17.and.raw_calls==17.and.all(counts==[1_int64,2_int64,4_int64]), &
      'MPI partial-clock native replay evaluates only unseen source stages')
    equal=.true.
    do node=1,np
      if(owners(node)/=rank) cycle
      equal=equal.and.all(state%fields(node)%values==expected%nodes(node)%state).and. &
        all(local_held%fields(node)%values==expected_h%fields(node)%values).and. &
        equal_flux(lf(node),ef(node)).and.equal_flux(lh(node),eh(node))
    enddo
    do r=1,7
      equal=equal.and.record_equal(local_record(r),record(r))
    enddo
    call require(equal,'exact serial/MPI fields/stages/faces after child source-prefix replay')
    call require(record_field_equal(local_record(2)%routed_auxiliary_impulse,expected_impulse), &
      'owner-local routed impulse matches independently integrated old/new packet')
    saved_gas=state;saved_temp=temperature;saved_h=local_held
    do invalid=1,2
      bad=local_cache
      if(invalid==1.and.rank==owners(2)) bad(2)%routed_auxiliary_old_rate%fields(2)%values(nv+1,1,1,1)=nan_value
      if(invalid==2.and.rank==2) bad(2)%accepted_source_stages=0
      aux_calls=0;raw_calls=0;prefix=7.5_dp;work=9.5_dp
      call mpi_subtree_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused, &
        replay=bad,replay_prefix=.true.,counts=counts)
      call require(.not.accepted.and.aux_calls==0.and.raw_calls==0,'malformed partial/routed cache rejects before providers')
      call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
        equal_state(local_held,saved_h).and.prefix==7.5_dp.and.work==9.5_dp,'malformed partial clock full rollback')
    enddo
    if(rank==0) print *, 'PASS: partial child source-stage parent recovery, routed endpoints, native family and MPI rollback'
  end subroutine check_stage_parent
  subroutine check_source_stages()
    type(eb_reactive_tree_3d) :: result,fresh,expected
    type(amr_tree_state_3d) :: result_h,fresh_h,expected_h,saved_gas,saved_temp,saved_h
    type(eb_godunov_interval_3d) :: record,recovery,cache,local_record,local_recovery,bad
    type(eb_godunov_source_stage_3d) :: empty_stage
    type(eb_integrated_flux_3d), allocatable :: f(:),h(:),lf(:),lh(:),ef(:),eh(:)
    type(state_redist_register_3d), allocatable :: gr(:),hr(:),lgr(:),lhr(:)
    real(dp) :: b(nv),w(nv),c(nv),v(nv),a(nv+1),theta,work,prefix,expected_work,expected_prefix
    integer :: node,invalid
    logical :: accepted,refused,ready,equal
    owners=rank;result=tree;result_h=held;prefix=7.5_dp;work=9.5_dp
    aux_calls=0;raw_calls=0;provider_must_fail=.false.;stage_sink=.true.
    call serial_source_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted, &
      recovery=recovery,recovery_ready=ready,refused=refused)
    call require(.not.accepted.and.ready.and..not.refused,'failed native predictor yields accepted old-source prefix')
    call require(recovery%source_prefix_only.and.recovery%accepted_source_stages==1,'only stage0 is admitted')
    call require(aux_calls==1.and.raw_calls==1,'accepted stage0 providers run once before failure')
    equal=prefix==7.5_dp.and.work==9.5_dp.and.record%accepted_source_stages==0
    do node=1,np
      equal=equal.and.all(result%nodes(node)%state==tree%nodes(node)%state).and. &
        all(result_h%fields(node)%values==held%fields(node)%values)
    enddo
    call require(equal,'failed partial stage retains caller physical state, history and ledgers')
    call require(.not.allocated(f).and..not.allocated(h).and..not.allocated(gr).and..not.allocated(hr), &
      'no partial physical faces or registers publish')
    cache=recovery;fresh=tree;fresh_h=held;stage_sink=.false.
    do node=1,np
      fresh%nodes(node)%state=1e9_dp*fresh%nodes(node)%state
      fresh_h%fields(node)%values=1e9_dp*fresh_h%fields(node)%values
    enddo
    result=fresh;result_h=fresh_h;prefix=.125_dp;aux_calls=0;raw_calls=0
    call serial_source_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted, &
      replay=cache,partial=.true.,recovery=recovery,recovery_ready=ready)
    call require(accepted.and..not.ready,'partial stage replay reevaluates native physics on fresh admissible caller')
    call require(aux_calls==2.and.raw_calls==2,'only previously unvisited stages1/2 call providers')
    call require(record%accepted_source_stages==3.and..not.record%source_prefix_only,'complete native source table returns')
    expected=result;expected_h=result_h;ef=f;eh=h;expected_work=work;expected_prefix=prefix
    result=fresh;result_h=fresh_h;prefix=.125_dp;provider_must_fail=.true.;aux_calls=0;raw_calls=0
    call serial_source_step(result,result_h,recovery,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted,replay=record)
    call require(accepted.and.aux_calls==0.and.raw_calls==0,'all now accepted source stages freeze in later trial')
    equal=work==expected_work.and.prefix==expected_prefix
    do node=1,np
      equal=equal.and.all(result%nodes(node)%state==expected%nodes(node)%state).and. &
        all(result_h%fields(node)%values==expected_h%fields(node)%values).and. &
        equal_flux(f(node),ef(node)).and.equal_flux(h(node),eh(node))
    enddo
    call require(equal.and.record_equal(record,recovery),'complete replay preserves native field/stage/face result')
    provider_must_fail=.false.
    if(nranks==1) return
    owners=[1,0,1,0,1];call distribute();stage_sink=.true.;prefix=7.5_dp;work=9.5_dp;aux_calls=0;raw_calls=0
    saved_gas=state;saved_temp=temperature;saved_h=local_held
    call mpi_source_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused, &
      recovery=local_recovery,recovery_ready=ready)
    call require(.not.accepted.and.ready.and..not.refused,'MPI failed predictor returns only accepted source prefix')
    call require(local_recovery%accepted_source_stages==1.and.local_recovery%source_prefix_only,'MPI stage0 recovery identity')
    call require(aux_calls==1.and.raw_calls==1,'MPI old source accepted once on owners and idle rank')
    call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
      equal_state(local_held,saved_h).and.prefix==7.5_dp.and.work==9.5_dp,'MPI partial source failure full rollback')
    cache=local_recovery;stage_sink=.false.;prefix=.125_dp
    do node=1,np
      if(owners(node)/=rank) cycle
      state%fields(node)%values=fresh%nodes(node)%state
      local_held%fields(node)%values=fresh_h%fields(node)%values
    enddo
    aux_calls=0;raw_calls=0
    call mpi_source_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused, &
      replay=cache,partial=.true.,recovery=local_recovery,recovery_ready=ready)
    call require(accepted.and..not.ready.and.aux_calls==2.and.raw_calls==2,'MPI partial replay runs only new stages')
    equal=.true.
    do node=1,np
      if(owners(node)/=rank) cycle
      equal=equal.and.all(state%fields(node)%values==expected%nodes(node)%state).and. &
        all(local_held%fields(node)%values==expected_h%fields(node)%values).and. &
        equal_flux(lf(node),ef(node)).and.equal_flux(lh(node),eh(node))
    enddo
    call require(equal.and.record_equal(local_record,record),'exact serial/MPI fresh native partial-stage replay')
    saved_gas=state;saved_temp=temperature;saved_h=local_held
    do invalid=1,3
      bad=cache
      if(invalid==1.and.rank==2) bad%accepted_source_stages=0
      if(invalid==2.and.rank==2) then
        bad=local_record;bad%accepted_source_stages=2;bad%source_prefix_only=.true.;bad%sources(2)=empty_stage
      endif
      if(invalid==3.and.rank==owners(2)) bad%sources(0)%auxiliary%fields(2)%values(nv+1,1,1,1)=nan_value
      aux_calls=0;raw_calls=0;prefix=7.5_dp;work=9.5_dp
      call mpi_source_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused, &
        replay=bad,partial=.true.,recovery=local_recovery,recovery_ready=ready)
      call require(.not.accepted.and..not.ready.and.aux_calls==0.and.raw_calls==0,'bad/idle stage prefixes reject before callbacks')
      call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
        equal_state(local_held,saved_h).and.prefix==7.5_dp.and.work==9.5_dp,'bad partial cache full rollback')
    enddo
    provider_must_fail=.true.;aux_calls=0;raw_calls=0
    call mpi_source_step(local_record,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused, &
      replay=cache,partial=.true.,recovery=local_recovery,recovery_ready=ready)
    call require(.not.accepted.and..not.ready.and.aux_calls<=1.and.raw_calls<=1,'refused unseen stage does not publish a cache')
    call require(aux_calls+raw_calls>0,'new refused stage was reached exactly once')
    provider_must_fail=.false.
    if(rank==0) print *, 'PASS: accepted source-stage prefix, fresh native replay, first evaluation/frozen reuse and MPI rollback'
  end subroutine check_source_stages
  subroutine check_recovery_replay()
    type(eb_reactive_tree_3d) :: result,saved_result
    type(amr_tree_state_3d) :: result_h,saved_h,saved_gas,saved_temp
    type(eb_godunov_interval_3d), allocatable :: reference(:),record(:),recovery(:),cache(:),local_reference(:)
    type(eb_integrated_flux_3d), allocatable :: f(:),h(:)
    type(state_redist_register_3d), allocatable :: gr(:),hr(:),saved_gr(:),saved_hr(:)
    type(amr_tree_state_3d), allocatable :: pending(:),saved_pending(:)
    real(dp) :: b(nv),w(nv),c(nv),v(nv),a(nv+1),theta,work,prefix,saved_work,saved_prefix
    integer(int64) :: counts(3),saved_counts(3)
    integer :: pass,node,r,bucket
    logical :: accepted,refused,ready,equal
    owners=rank;result=tree;result_h=held;prefix=.125_dp;provider_must_fail=.false.;recovery_fixture=.false.
    call serial_family_step(result,result_h,reference,b,w,c,v,a,theta,work,prefix,f,h,accepted)
    call require(accepted,'native reference for failed-prefix recovery')
    record=reference;saved_result=result;saved_h=result_h
    cache=reference(5:5);prefix=7.5_dp;work=9.5_dp;counts=-19
    saved_work=work;saved_prefix=prefix;saved_counts=counts
    do pass=1,2
      aux_calls=0;raw_calls=0;provider_must_fail=pass==2;recovery_fixture=.true.
      call serial_subtree_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,accepted, &
        replay=cache,gr=gr,hr=hr,pending=pending,counts=counts,replay_prefix=pass==1, &
        recovery=recovery,recovery_ready=ready,refused=refused)
      call require(.not.accepted.and..not.refused.and.ready,'failed serial closure yields only accepted source cache')
      call require(allocated(recovery),'serial failed closure recovery table present')
      call require(size(recovery)==3,'serial recovery includes both newly evaluated child clocks')
      call require(aux_calls==merge(6,0,pass==1).and.raw_calls==merge(6,0,pass==1),'serial recovered sources never repeat')
      call require(all(recovery(2:)%replayed.eqv.(pass==2)),'serial recovery distinguishes native and replayed children')
      equal=work==saved_work.and.prefix==saved_prefix.and.all(counts==saved_counts)
      do node=1,np
        equal=equal.and.all(result%nodes(node)%state==saved_result%nodes(node)%state).and. &
          all(result%nodes(node)%temperature==saved_result%nodes(node)%temperature).and. &
          all(result_h%fields(node)%values==saved_h%fields(node)%values)
      enddo
      do r=1,7
        equal=equal.and.record_equal(record(r),reference(r))
      enddo
      call require(equal,'serial diagnostic publication retains all physical state and accepted history')
      call require(.not.allocated(gr).and..not.allocated(hr).and..not.allocated(pending),'serial failed ledgers remain absent')
      call require(all(recovery%parent_record==[1,5,5]).and.all(recovery%subtree_last==[7,6,7]), &
        'serial diagnostic recovery retains global clock metadata')
      cache=recovery
    enddo
    recovery_fixture=.false.;provider_must_fail=.false.
    if(nranks==1) return
    owners=[1,0,1,0,1];call distribute();prefix=.125_dp
    call mpi_family_step(local_reference,b,w,c,v,a,theta,work,prefix,f,h,accepted,refused)
    call require(accepted,'MPI reference for failed-prefix recovery')
    record=local_reference;saved_gas=state;saved_temp=temperature;saved_h=local_held
    cache=local_reference(5:5);prefix=7.5_dp;work=9.5_dp;counts=-19
    saved_work=work;saved_prefix=prefix;saved_counts=counts
    do pass=1,2
      aux_calls=0;raw_calls=0;provider_must_fail=pass==2;recovery_fixture=.true.
      call mpi_subtree_step(record,b,w,c,v,a,theta,work,prefix,f,h,accepted,refused, &
        replay=cache,gr=gr,hr=hr,pending=pending,counts=counts,replay_prefix=pass==1, &
        recovery=recovery,recovery_ready=ready)
      call require(.not.accepted.and..not.refused.and.ready,'MPI failed closure returns actual accepted common sources')
      call require(allocated(recovery),'MPI recovery table present on every rank')
      call require(size(recovery)==3,'MPI recovery includes both first-evaluated children')
      call require(aux_calls==merge(6,0,pass==1).and.raw_calls==merge(6,0,pass==1),'MPI recovered providers run once only')
      equal=equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and.equal_state(local_held,saved_h).and. &
        work==saved_work.and.prefix==saved_prefix.and.all(counts==saved_counts)
      do r=1,7
        equal=equal.and.record_equal(record(r),local_reference(r))
      enddo
      call require(equal,'MPI failed source recovery retains accepted physical fields/history/counts/ledgers')
      call require(.not.allocated(gr).and..not.allocated(hr).and..not.allocated(pending),'MPI failed parent/pending not published')
      call require(all(recovery%parent_record==[1,5,5]).and.all(recovery%subtree_last==[7,6,7]), &
        'MPI failed recovery retains source clock identity including idle rank')
      cache=recovery
    enddo
    ! A refused new provider cannot publish a fabricated complete source prefix.
    cache=local_reference(5:5);aux_calls=0;raw_calls=0;provider_must_fail=.true.
    call mpi_subtree_step(record,b,w,c,v,a,theta,work,prefix,f,h,accepted,refused, &
      replay=cache,gr=gr,hr=hr,pending=pending,counts=counts,replay_prefix=.true., &
      recovery=recovery,recovery_ready=ready)
    call require(.not.accepted.and..not.ready,'incomplete future source has no recovery table')
    call require(aux_calls<=1.and.raw_calls<=1.and.aux_calls+raw_calls>0,'refused unseen provider does not repeat')
    call require(size(recovery)==3,'incomplete recovery preserves previous diagnostic cache')
    recovery_fixture=.false.;provider_must_fail=.false.
    if(rank==0) print *, 'PASS: failed native closure source recovery, first evaluation then frozen reuse, complete rollback'
  end subroutine check_recovery_replay
  subroutine check_subtree_replay()
    type(eb_reactive_tree_3d) :: result
    type(amr_tree_state_3d) :: result_h,saved_gas,saved_temp,saved_h
    type(eb_godunov_interval_3d), allocatable :: reference(:),replacement(:),local_reference(:),local_result(:),bad(:)
    type(eb_godunov_interval_3d), allocatable :: saved_records(:)
    type(eb_integrated_flux_3d), allocatable :: f(:),h(:),lf(:),lh(:)
    type(state_redist_register_3d), allocatable :: gr(:),hr(:),expected_gr(:),expected_hr(:),local_gr(:),local_hr(:)
    type(state_redist_register_3d), allocatable :: saved_gr(:),saved_hr(:)
    type(state_redist_register_3d) :: part
    type(amr_tree_state_3d), allocatable :: pending(:),expected_pending(:),local_pending(:),saved_pending(:)
    real(dp) :: b(nv),w(nv),c(nv),v(nv),a(nv+1),theta,work,prefix,expected_work,saved_work,saved_prefix
    integer(int64) :: counts(3),saved_counts(3)
    integer :: node,r,invalid,bucket,expected_calls
    logical :: accepted,refused,equal,nonzero,prefix_mode
    owners=rank;result=tree;result_h=held;prefix=.125_dp;provider_must_fail=.false.
    call serial_family_step(result,result_h,reference,b,w,c,v,a,theta,work,prefix,f,h,accepted)
    call require(accepted.and.size(reference)==7,'native routed seven-clock reference')
    expected_gr=reference(5)%open_ancestor_registers;expected_hr=reference(5)%open_ancestor_scalar_registers
    expected_pending=reference(5)%clock_routed_pending
    nonzero=.false.
    do node=1,np
      nonzero=nonzero.or.any(expected_pending(0)%fields(node)%values/=0)
      if(topology%nodes(node)%level/=1) cycle
      part=reference(5)%registers(node)
      if(allocated(part%coarse_delta)) part%coarse_delta=0
      call accumulate_state_redist_register_3d(expected_gr(node),part,nv,[n,n,n]+12,accepted)
      call require(accepted,'independent native parent fine gas registration')
      part=reference(5)%scalar_registers(node)
      if(allocated(part%coarse_delta)) part%coarse_delta=0
      call accumulate_state_redist_register_3d(expected_hr(node),part,1,[n,n,n]+12,accepted)
      call require(accepted,'independent native parent fine UEINT registration')
    enddo
    call require(nonzero,'second parent group starts with prior sibling/ancestor routed source')
    do r=5,7
      do node=1,np
        bucket=min(reference(r)%depth,topology%nodes(node)%level)
        if(bucket/=0) cycle
        expected_pending(0)%fields(node)%values=expected_pending(0)%fields(node)%values+ &
          reference(r)%routed_auxiliary_impulse%fields(node)%values
      enddo
    enddo
    do bucket=1,2
      do node=1,np
        expected_pending(bucket)%fields(node)%values=0
      enddo
    enddo
    expected_work=sum(reference(5:7)%pressure_work_change)
    ! An already advanced and deliberately perturbed caller tail must be
    ! discarded in favour of the complete accepted pre-reset group clock.
    do node=1,np
      result%nodes(node)%state=1.1_dp*result%nodes(node)%state
      result_h%fields(node)%values=1.4_dp*result_h%fields(node)%values
    enddo
    expected_calls=merge(6,0,prefix_replay_test)
    prefix=7.5_dp;counts=-19;aux_calls=0;raw_calls=0;provider_must_fail=.not.prefix_replay_test
    prefix_forbid_parent=prefix_replay_test
    call serial_subtree_step(result,result_h,replacement,b,w,c,v,a,theta,work,prefix,f,h,accepted, &
      replay=reference(5:merge(5,7,prefix_replay_test)),gr=gr,hr=hr,pending=pending,counts=counts, &
      replay_prefix=prefix_replay_test)
    call require(accepted,'serial restored second parent group')
    call require(all(counts==[0_int64,1_int64,2_int64]).and.aux_calls==expected_calls.and.raw_calls==expected_calls, &
      'group clocks and first evaluations run only for unvisited prefix tail')
    equal=abs(work-expected_work)<=8*epsilon(1._dp)*max(tiny(1._dp),abs(expected_work))
    equal=equal.and.prefix==reference(5)%closed_energy_correction
    do node=1,np
      equal=equal.and.all(result%nodes(node)%state==reference(5)%closed_gas%fields(node)%values).and. &
        all(result%nodes(node)%temperature==reference(5)%closed_temperature%fields(node)%values(1,:,:,:)).and. &
        all(result_h%fields(node)%values==reference(5)%closed_internal_energy%fields(node)%values)
      equal=equal.and.equal_flux(f(node),reference(5)%closed_faces(node)).and. &
        equal_flux(h(node),reference(5)%closed_scalar_faces(node)).and. &
        equal_register(gr(node),expected_gr(node)).and.equal_register(hr(node),expected_hr(node))
    enddo
    do r=1,3
      equal=equal.and.replacement(r)%clock_energy_correction==reference(r+4)%clock_energy_correction.and. &
        replacement(r)%closed_energy_correction==reference(r+4)%closed_energy_correction.and. &
        abs(replacement(r)%clock_pressure_work_change-reference(r+4)%clock_pressure_work_change)<= &
        8*epsilon(1._dp)*max(tiny(1._dp),abs(reference(r+4)%clock_pressure_work_change)).and. &
        abs(replacement(r)%closed_pressure_work_change-reference(r+4)%closed_pressure_work_change)<= &
        8*epsilon(1._dp)*max(tiny(1._dp),abs(reference(r+4)%closed_pressure_work_change))
      equal=equal.and.record_equal(replacement(r),reference(r+4)).and. &
        replacement(r)%parent_record==reference(r+4)%parent_record.and. &
        replacement(r)%subtree_last==reference(r+4)%subtree_last.and.replacement(r)%time==reference(r+4)%time
    enddo
    do bucket=0,2
      equal=equal.and.equal_state(pending(bucket),expected_pending(bucket))
    enddo
    call require(equal,'serial group restores actual clock, hierarchy, faces, stages, parent registers and routed buckets')
    call require(replacement(1)%replayed,'accepted parent source is replayed')
    call require(all(replacement(2:)%replayed.eqv..not.prefix_replay_test),'unvisited tail runs native only once')
    prefix_forbid_parent=.false.
    if(prefix_replay_test) then
      bad=replacement;aux_calls=0;raw_calls=0;provider_must_fail=.true.
      call serial_subtree_step(result,result_h,replacement,b,w,c,v,a,theta,work,prefix,f,h,accepted, &
        replay=bad,gr=gr,hr=hr,pending=pending,counts=counts)
      call require(accepted,'complete serial prefix result can seed another native trial')
      call require(aux_calls==0.and.raw_calls==0.and.all(replacement%replayed),'new tail packets freeze in later trials')
      equal=prefix==reference(5)%closed_energy_correction
      do node=1,np
        equal=equal.and.all(result%nodes(node)%state==reference(5)%closed_gas%fields(node)%values).and. &
          all(result_h%fields(node)%values==reference(5)%closed_internal_energy%fields(node)%values).and. &
          equal_flux(f(node),reference(5)%closed_faces(node)).and.equal_flux(h(node),reference(5)%closed_scalar_faces(node))
      enddo
      do r=1,3
        equal=equal.and.record_equal(replacement(r),reference(r+4))
      enddo
      call require(equal,'frozen serial tail preserves native physical/stage/face result')
    endif
    if(nranks==1) return
    owners=[1,0,1,0,1];provider_must_fail=.false.;call distribute();prefix=.125_dp
    call mpi_family_step(local_reference,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused)
    call require(accepted,'MPI owner-local routed native reference')
    do node=1,np
      if(owners(node)/=rank) cycle
      state%fields(node)%values=1.1_dp*state%fields(node)%values
      local_held%fields(node)%values=1.4_dp*local_held%fields(node)%values
    enddo
    prefix=7.5_dp;counts=-19;provider_must_fail=.not.prefix_replay_test;aux_calls=0;raw_calls=0
    prefix_forbid_parent=prefix_replay_test
    call mpi_subtree_step(local_result,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused, &
      replay=local_reference(5:merge(5,7,prefix_replay_test)),gr=local_gr,hr=local_hr,pending=local_pending,counts=counts, &
      replay_prefix=prefix_replay_test)
    call require(accepted,'MPI restored second parent group accepted')
    call require(all(counts==[0_int64,1_int64,2_int64]).and.aux_calls==expected_calls.and.raw_calls==expected_calls, &
      'MPI native group clock/source policy')
    equal=abs(work-expected_work)<=8*epsilon(1._dp)*max(tiny(1._dp),abs(expected_work))
    equal=equal.and.abs(prefix-reference(5)%closed_energy_correction)<= &
      8*epsilon(1._dp)*max(tiny(1._dp),abs(reference(5)%closed_energy_correction))
    do node=1,np
      if(owners(node)/=rank) cycle
      equal=equal.and.all(state%fields(node)%values==result%nodes(node)%state).and. &
        all(temperature%fields(node)%values(1,:,:,:)==result%nodes(node)%temperature).and. &
        all(local_held%fields(node)%values==result_h%fields(node)%values)
      equal=equal.and.equal_flux(lf(node),f(node)).and.equal_flux(lh(node),h(node)).and. &
        equal_register(local_gr(node),expected_gr(node)).and.equal_register(local_hr(node),expected_hr(node))
    enddo
    do r=1,3
      equal=equal.and.abs(local_result(r)%clock_energy_correction-reference(r+4)%clock_energy_correction)<= &
        8*epsilon(1._dp)*max(tiny(1._dp),abs(reference(r+4)%clock_energy_correction)).and. &
        abs(local_result(r)%closed_energy_correction-reference(r+4)%closed_energy_correction)<= &
        8*epsilon(1._dp)*max(tiny(1._dp),abs(reference(r+4)%closed_energy_correction)).and. &
        abs(local_result(r)%clock_pressure_work_change-reference(r+4)%clock_pressure_work_change)<= &
        8*epsilon(1._dp)*max(tiny(1._dp),abs(reference(r+4)%clock_pressure_work_change)).and. &
        abs(local_result(r)%closed_pressure_work_change-reference(r+4)%closed_pressure_work_change)<= &
        8*epsilon(1._dp)*max(tiny(1._dp),abs(reference(r+4)%closed_pressure_work_change))
      equal=equal.and.record_equal(local_result(r),reference(r+4))
    enddo
    do bucket=0,2
      equal=equal.and.record_field_equal(local_pending(bucket),expected_pending(bucket))
    enddo
    call require(equal,'exact serial/MPI restored parent group and prior sibling ledgers')
    call require(local_result(1)%replayed,'MPI prefix parent source is frozen')
    call require(all(local_result(2:)%replayed.eqv..not.prefix_replay_test),'MPI unvisited clocks evaluate native sources')
    prefix_forbid_parent=.false.;provider_must_fail=.true.;aux_calls=0;raw_calls=0
    if(prefix_replay_test) then
      bad=local_result
      call mpi_subtree_step(local_result,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused, &
        replay=bad,gr=local_gr,hr=local_hr,pending=local_pending,counts=counts)
      call require(accepted,'complete MPI prefix result seeds a later native trial')
      call require(aux_calls==0.and.raw_calls==0.and.all(local_result%replayed),'MPI new tail packets freeze without callbacks')
      equal=.true.
      do node=1,np
        if(owners(node)/=rank) cycle
        equal=equal.and.all(state%fields(node)%values==result%nodes(node)%state).and. &
          all(local_held%fields(node)%values==result_h%fields(node)%values).and. &
          equal_flux(lf(node),f(node)).and.equal_flux(lh(node),h(node))
      enddo
      do r=1,3
        equal=equal.and.record_equal(local_result(r),reference(r+4))
      enddo
      call require(equal,'MPI frozen prefix tail preserves actual native states/stages/faces')
    endif
    saved_gas=state;saved_temp=temperature;saved_h=local_held;saved_work=work;saved_prefix=prefix;saved_counts=counts
    saved_gr=local_gr;saved_hr=local_hr;saved_pending=local_pending;saved_records=local_result
    do invalid=1,merge(7,5,prefix_replay_test)
      bad=local_reference(5:merge(5,7,prefix_replay_test))
      prefix_mode=prefix_replay_test
      if(invalid==6.and.rank==2) prefix_mode=.false.
      if(invalid==7.and.rank==2) bad=local_reference(5:6)
      if(rank==owners(2)) then
        select case(invalid)
        case(1)
          bad(1)%open_ancestor_scalar_registers(2)%scale=.5_dp
        case(2)
          bad(1)%clock_reaction_source%fields(2)%values(1,1,1,1)=nan_value
        case(3)
          bad(1)%clock_routed_pending(0)%fields(2)%values(nv+1,1,1,1)=nan_value
        case(4)
          bad(1)%open_ancestor_registers(2)%coarse_flags(1,1,1)= &
            mod(bad(1)%open_ancestor_registers(2)%coarse_flags(1,1,1)+1,3)
        end select
      endif
      if(invalid==5.and.rank==2) bad(size(bad))%time=bad(size(bad))%time+1e-10_dp
      call mpi_subtree_step(local_result,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused, &
        replay=bad,gr=local_gr,hr=local_hr,pending=local_pending,counts=counts,replay_prefix=prefix_mode)
      call require(.not.accepted.and..not.refused.and.aux_calls==0.and.raw_calls==0,'bad subtree anchor refuses before providers')
      equal=equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and.equal_state(local_held,saved_h).and. &
        work==saved_work.and.prefix==saved_prefix.and.all(counts==saved_counts)
      do node=1,np
        equal=equal.and.equal_register(local_gr(node),saved_gr(node)).and.equal_register(local_hr(node),saved_hr(node))
      enddo
      do bucket=0,2
        equal=equal.and.equal_state(local_pending(bucket),saved_pending(bucket))
      enddo
      do r=1,3
        equal=equal.and.record_equal(local_result(r),saved_records(r))
      enddo
      call require(equal,'bad subtree anchor preserves all prior state, clocks, parent ledgers and routed buckets')
    enddo
    if(prefix_replay_test) then
      aux_calls=0;raw_calls=0
      call mpi_subtree_step(local_result,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused, &
        replay=local_reference(5:5),gr=local_gr,hr=local_hr,pending=local_pending,counts=counts,replay_prefix=.true.)
      call require(.not.accepted,'new tail provider refusal rejects whole group')
      call require(aux_calls<=1.and.raw_calls<=1.and.aux_calls+raw_calls>0,'refused future source is never retried')
      equal=equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and.equal_state(local_held,saved_h).and. &
        work==saved_work.and.prefix==saved_prefix.and.all(counts==saved_counts)
      do node=1,np
        equal=equal.and.equal_register(local_gr(node),saved_gr(node)).and.equal_register(local_hr(node),saved_hr(node))
      enddo
      do bucket=0,2
        equal=equal.and.equal_state(local_pending(bucket),saved_pending(bucket))
      enddo
      do r=1,3
        equal=equal.and.record_equal(local_result(r),saved_records(r))
      enddo
      call require(equal,'future source refusal preserves accepted caller state and all published ledgers')
    endif
    if(rank==0) print *, 'PASS: restored second Godunov parent group, native stages/IR, open parent and routed sibling ledgers'
  end subroutine check_subtree_replay
  subroutine serial_subtree_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,accepted,replay, &
      declared_source,evaluate_reaction,refused,gr,hr,pending,counts,replay_prefix,recovery,recovery_ready)
    type(eb_reactive_tree_3d), intent(inout) :: result
    type(amr_tree_state_3d), intent(inout) :: result_h
    type(eb_godunov_interval_3d), allocatable, intent(inout) :: record(:)
    type(eb_godunov_interval_3d), optional, intent(in) :: replay(:)
    real(dp), intent(inout) :: b(:),w(:),c(:),v(:),a(:),theta,work,prefix
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: f(:),h(:)
    logical, optional, intent(out) :: refused
    logical, intent(out) :: accepted
    logical, optional, intent(in) :: declared_source,replay_prefix
    type(eb_godunov_interval_3d), allocatable, optional, intent(inout) :: recovery(:)
    logical, optional, intent(out) :: recovery_ready
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: gr(:),hr(:)
    type(amr_tree_state_3d), allocatable, optional, intent(inout) :: pending(:)
    integer(int64), optional, intent(inout) :: counts(:)
    real(dp) :: duration
    logical :: is_source
    type(elementary_reaction), allocatable :: reactions(:)
    type(gas_transport_species), allocatable :: transport(:)
    duration=1e-8_dp
    if(present(replay)) duration=replay(1)%interval
    is_source=.true.
    if(present(declared_source)) is_source=declared_source
    call load_h2o2_elementary_mechanism(reactions,accepted)
    call load_h2o2_elementary_transport(transport,accepted)
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,duration,2,.true., &
      [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,result,extended,faces, &
      b,w,c,v,theta,accepted,redistribution='state_redist',integrated_flux=f, &
      internal_energy=result_h,internal_energy_flux=h, &
      internal_energy_pressure_work_change=work,energy_correction=prefix,spatial_flux_weight=.375_dp, &
      auxiliary_source=serial_recursive_aux,auxiliary_change=a,auxiliary_has_internal_energy=.true., &
      evaluate_spatial_diffusion=serial_raw,spatial_diffusion_is_source=is_source, &
      godunov_interval_history=record,godunov_replay_history=replay,route_auxiliary=routed_fixture, &
      godunov_replay_prefix=replay_prefix,godunov_recovery_history=recovery,godunov_recovery_ready=recovery_ready, &
      godunov_parent_registers=gr,godunov_parent_scalar_registers=hr,godunov_routed_pending=pending,level_steps=counts, &
      evaluate_reaction=evaluate_reaction,provider_refused=refused)
  end subroutine serial_subtree_step
  subroutine mpi_subtree_step(record,b,w,c,v,a,theta,work,prefix,f,h,accepted,refused,replay,declared_source,evaluate_reaction, &
      gr,hr,pending,counts,replay_prefix,recovery,recovery_ready)
    type(eb_godunov_interval_3d), allocatable, intent(inout) :: record(:)
    type(eb_godunov_interval_3d), optional, intent(in) :: replay(:)
    real(dp), intent(inout) :: b(:),w(:),c(:),v(:),a(:),theta,work,prefix
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: f(:),h(:)
    logical, intent(out) :: accepted,refused
    logical, optional, intent(in) :: declared_source,replay_prefix
    type(eb_godunov_interval_3d), allocatable, optional, intent(inout) :: recovery(:)
    logical, optional, intent(out) :: recovery_ready
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: gr(:),hr(:)
    type(amr_tree_state_3d), allocatable, optional, intent(inout) :: pending(:)
    integer(int64), optional, intent(inout) :: counts(:)
    real(dp) :: duration
    logical :: is_source
    type(elementary_reaction), allocatable :: reactions(:)
    type(gas_transport_species), allocatable :: transport(:)
    duration=1e-8_dp
    if(present(replay)) duration=replay(1)%interval
    is_source=.true.
    if(present(declared_source)) is_source=declared_source
    call load_h2o2_elementary_mechanism(reactions,accepted)
    call load_h2o2_elementary_transport(transport,accepted)
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,state,temperature, &
      duration,2,.true.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,extended,faces, &
      b,w,c,v,theta,MPI_COMM_WORLD,accepted,redistribution='state_redist',integrated_flux=f, &
      internal_energy=local_held,internal_energy_flux=h, &
      internal_energy_pressure_work_change=work,energy_correction=prefix,spatial_flux_weight=.375_dp, &
      auxiliary_source=mpi_recursive_aux,auxiliary_change=a,auxiliary_has_internal_energy=.true., &
      evaluate_spatial_diffusion=mpi_raw,spatial_diffusion_is_source=is_source,provider_refused=refused, &
      godunov_interval_history=record,godunov_replay_history=replay,route_auxiliary=routed_fixture, &
      godunov_replay_prefix=replay_prefix,godunov_recovery_history=recovery,godunov_recovery_ready=recovery_ready, &
      godunov_parent_registers=gr,godunov_parent_scalar_registers=hr,godunov_routed_pending=pending,level_steps=counts, &
      evaluate_reaction=evaluate_reaction)
  end subroutine mpi_subtree_step
  subroutine check_family_replay()
    type(elementary_reaction), allocatable :: reactions(:)
    type(gas_transport_species), allocatable :: transport(:)
    type(eb_reactive_tree_3d) :: result,expected,fresh
    type(amr_tree_state_3d) :: result_h,expected_h,saved_gas,saved_temp,saved_h
    type(eb_godunov_interval_3d), allocatable :: original(:),reference(:),replayed(:), &
      local_original(:),local_result(:),bad(:),saved_records(:)
    type(eb_integrated_flux_3d), allocatable :: f(:),h(:),ef(:),eh(:),lf(:),lh(:)
    real(dp) :: b(nv),w(nv),c(nv),v(nv),a(nv+1),theta,work,prefix,expected_work,expected_prefix
    real(dp) :: eb(nv),ew(nv),ec(nv),ev(nv),ea(nv+1),etheta,saved_work,saved_prefix
    integer :: pass,node,owner_layout,invalid,r
    logical :: accepted,equal,refused
    call load_h2o2_elementary_mechanism(reactions,accepted);call require(accepted,'replay reaction DB')
    call load_h2o2_elementary_transport(transport,accepted);call require(accepted,'replay transport DB')
    fresh=tree
    do node=1,np
      fresh%nodes(node)%state=1.01_dp*tree%nodes(node)%state
    enddo
    do pass=1,3
      result=tree;result_h=held
      if(pass>1) then
        result=fresh
        do node=1,np
          result_h%fields(node)%values=1.03_dp*held%fields(node)%values
        enddo
      endif
      prefix=.125_dp;aux_calls=0;raw_calls=0;provider_must_fail=pass==3
      if(pass==1) then
        call serial_family_step(result,result_h,original,b,w,c,v,a,theta,work,prefix,f,h,accepted)
      else if(pass==2) then
        call serial_family_step(result,result_h,reference,b,w,c,v,a,theta,work,prefix,f,h,accepted)
      else
        call serial_family_step(result,result_h,replayed,b,w,c,v,a,theta,work,prefix,f,h,accepted,original)
      endif
      call require(accepted,'fresh native Godunov replay accepted')
      if(pass<3) then
        call require(aux_calls==21.and.raw_calls==21,'original providers run at each actual source stage')
      else
        call require(aux_calls==0.and.raw_calls==0.and.all(replayed%replayed),'cached sources never reinvoke providers')
      endif
      if(pass==2) then
        expected=result;expected_h=result_h;ef=f;eh=h
        eb=b;ew=w;ec=c;ev=v;ea=a;etheta=theta;expected_work=work;expected_prefix=prefix
        call require(work/=sum(original%pressure_work_change),'fresh native pressure work changes with fresh state')
      endif
      if(pass/=3) cycle
      equal=prefix==expected_prefix.and.work==expected_work.and.all(b==eb).and.all(w==ew).and. &
        all(c==ec).and.all(v==ev).and.all(a==ea).and.theta==etheta
      do node=1,np
        equal=equal.and.all(result%nodes(node)%state==expected%nodes(node)%state).and. &
          all(result%nodes(node)%temperature==expected%nodes(node)%temperature).and. &
          all(result_h%fields(node)%values==expected_h%fields(node)%values)
        equal=equal.and.equal_flux(f(node),ef(node)).and.equal_flux(h(node),eh(node))
      enddo
      call require(equal,'serial cached source replay equals fresh live native computation')
    enddo
    if(nranks==1) return
    do owner_layout=1,2
      owners=[0,1,2,0,1];if(owner_layout==2) owners=[1,0,1,0,1]
      provider_must_fail=.false.;aux_calls=0;raw_calls=0
      call distribute();prefix=.125_dp
      call mpi_family_step(local_original,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused)
      call require(accepted.and.aux_calls==21.and.raw_calls==21,'owner-local native provider capture')
      call distribute()
      do node=1,np
        if(owners(node)/=rank) cycle
        state%fields(node)%values=fresh%nodes(node)%state
        local_held%fields(node)%values=1.03_dp*held%fields(node)%values
      enddo
      provider_must_fail=.true.;aux_calls=0;raw_calls=0;prefix=.125_dp
      call mpi_family_step(local_result,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused,local_original)
      call require(accepted,'MPI source replay accepted before output inspection')
      call require(allocated(local_result),'MPI accepted family history allocated')
      call require(accepted.and.aux_calls==0.and.raw_calls==0.and.all(local_result%replayed),'MPI source replay')
      equal=abs(work-expected_work)<=8*epsilon(1._dp)*max(tiny(1._dp),abs(expected_work)).and. &
        abs(prefix-expected_prefix)<=8*epsilon(1._dp)*max(tiny(1._dp),abs(expected_prefix))
      do node=1,np
        if(owners(node)/=rank) cycle
        equal=equal.and.all(state%fields(node)%values==expected%nodes(node)%state).and. &
          all(temperature%fields(node)%values(1,:,:,:)==expected%nodes(node)%temperature).and. &
          all(local_held%fields(node)%values==expected_h%fields(node)%values)
        equal=equal.and.equal_flux(lf(node),ef(node)).and.equal_flux(lh(node),eh(node))
      enddo
      call require(equal,'MPI fresh replay gas/T/UEINT/faces/registers equal independent live serial result')
      do r=1,7
        call require(record_equal(local_result(r),reference(r)), &
          'all fresh native source/hydro/register/closed-clock payloads match live serial history')
      enddo
      saved_records=local_result
      saved_gas=state;saved_temp=temperature;saved_h=local_held;saved_work=work;saved_prefix=prefix
      do invalid=1,merge(4,3,routed_fixture)
        bad=local_original
        if(rank==owners(4)) then
          if(invalid==1) bad(7)%sources(2)%auxiliary%fields(4)%values(nv+1,1,1,1)=nan_value
          if(invalid==2) bad(7)%sources(2)%provider_faces(4)%x(2,0,1,1)=nan_value
        endif
        if(invalid==3.and.rank==2) bad(7)%parent_record=1
        if(invalid==4.and.rank==owners(4)) &
          bad(7)%routed_auxiliary_impulse%fields(4)%values(nv+1,1,1,1)=nan_value
        call mpi_family_step(local_result,b,w,c,v,a,theta,work,prefix,lf,lh,accepted,refused,bad)
        call require(.not.accepted.and..not.refused.and.aux_calls==0.and.raw_calls==0, &
          'late cached NaN/parent corruption refuses whole-family preflight')
        call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
          equal_state(local_held,saved_h).and.work==saved_work.and.prefix==saved_prefix, &
          'malformed cache preserves accepted state and scalar ledgers')
        do r=1,7
          call require(record_equal(local_result(r),saved_records(r)),'cache refusal preserves prior accepted clock family')
        enddo
      enddo
    enddo
    if(rank==0) print *, 'PASS: fresh Godunov full native clock replay, exact live serial/MPI state and faces; late cache rollback'
  end subroutine check_family_replay
  subroutine serial_family_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,accepted,replay, &
      declared_source,evaluate_reaction,refused)
    type(eb_reactive_tree_3d), intent(inout) :: result
    type(amr_tree_state_3d), intent(inout) :: result_h
    type(eb_godunov_interval_3d), allocatable, intent(inout) :: record(:)
    type(eb_godunov_interval_3d), optional, intent(in) :: replay(:)
    real(dp), intent(inout) :: b(:),w(:),c(:),v(:),a(:),theta,work,prefix
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: f(:),h(:)
    logical, optional, intent(out) :: refused
    logical, intent(out) :: accepted
    logical, optional, intent(in) :: declared_source
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    logical :: is_source
    type(elementary_reaction), allocatable :: reactions(:)
    type(gas_transport_species), allocatable :: transport(:)
    is_source=.true.
    if(present(declared_source)) is_source=declared_source
    call load_h2o2_elementary_mechanism(reactions,accepted)
    call load_h2o2_elementary_transport(transport,accepted)
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1e-8_dp,2,.true., &
      [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,result,extended,faces, &
      b,w,c,v,theta,accepted,redistribution='state_redist',integrated_flux=f, &
      internal_energy=result_h,internal_energy_flux=h, &
      internal_energy_pressure_work_change=work,energy_correction=prefix,spatial_flux_weight=.375_dp, &
      auxiliary_source=serial_recursive_aux,auxiliary_change=a,auxiliary_has_internal_energy=.true., &
      evaluate_spatial_diffusion=serial_raw,spatial_diffusion_is_source=is_source, &
      godunov_interval_history=record,godunov_replay_history=replay,route_auxiliary=routed_fixture, &
      evaluate_reaction=evaluate_reaction,provider_refused=refused)
  end subroutine serial_family_step
  subroutine mpi_family_step(record,b,w,c,v,a,theta,work,prefix,f,h,accepted,refused,replay,declared_source,evaluate_reaction)
    type(eb_godunov_interval_3d), allocatable, intent(inout) :: record(:)
    type(eb_godunov_interval_3d), optional, intent(in) :: replay(:)
    real(dp), intent(inout) :: b(:),w(:),c(:),v(:),a(:),theta,work,prefix
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: f(:),h(:)
    logical, intent(out) :: accepted,refused
    logical, optional, intent(in) :: declared_source
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    logical :: is_source
    type(elementary_reaction), allocatable :: reactions(:)
    type(gas_transport_species), allocatable :: transport(:)
    is_source=.true.
    if(present(declared_source)) is_source=declared_source
    call load_h2o2_elementary_mechanism(reactions,accepted)
    call load_h2o2_elementary_transport(transport,accepted)
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,state,temperature, &
      1e-8_dp,2,.true.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,extended,faces, &
      b,w,c,v,theta,MPI_COMM_WORLD,accepted,redistribution='state_redist',integrated_flux=f, &
      internal_energy=local_held,internal_energy_flux=h, &
      internal_energy_pressure_work_change=work,energy_correction=prefix,spatial_flux_weight=.375_dp, &
      auxiliary_source=mpi_recursive_aux,auxiliary_change=a,auxiliary_has_internal_energy=.true., &
      evaluate_spatial_diffusion=mpi_raw,spatial_diffusion_is_source=is_source,provider_refused=refused, &
      godunov_interval_history=record,godunov_replay_history=replay,route_auxiliary=routed_fixture, &
      evaluate_reaction=evaluate_reaction)
  end subroutine mpi_family_step
  subroutine serial_recursive_aux(gas,interval,source_stage,depth,time,rate,ok)
    type(eb_reactive_tree_3d), intent(in) :: gas
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: source_stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: ok
    call fill_aux(rate,depth,ok)
  end subroutine serial_recursive_aux
  subroutine mpi_recursive_aux(gas,temp,interval,source_stage,depth,time,rate,comm,ok)
    type(amr_tree_state_3d), intent(in) :: gas,temp
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: source_stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    call fill_aux(rate,depth,ok)
  end subroutine mpi_recursive_aux
  subroutine check_source_replay()
    type(elementary_reaction), allocatable :: reactions(:)
    type(gas_transport_species), allocatable :: transport(:)
    type(eb_reactive_tree_3d) :: result,expected,fresh
    type(amr_tree_state_3d) :: result_h,expected_h,saved_gas,saved_temp,saved_h
    type(eb_godunov_interval_3d) :: original,reference,replayed,local_original,local_result,bad
    type(eb_integrated_flux_3d), allocatable :: f(:),h(:),ef(:),eh(:),lf(:),lh(:)
    type(state_redist_register_3d), allocatable :: gr(:),hr(:),egr(:),ehr(:),lgr(:),lhr(:)
    real(dp) :: b(nv),w(nv),c(nv),v(nv),a(nv+1),theta,work,prefix,expected_work,expected_prefix
    real(dp) :: eb(nv),ew(nv),ec(nv),ev(nv),ea(nv+1),etheta,saved_work,saved_prefix
    integer :: pass,node,owner_layout,invalid
    logical :: accepted,equal,refused
    call load_h2o2_elementary_mechanism(reactions,accepted);call require(accepted,'replay reaction DB')
    call load_h2o2_elementary_transport(transport,accepted);call require(accepted,'replay transport DB')
    fresh=tree
    do node=1,np
      fresh%nodes(node)%state=1.01_dp*tree%nodes(node)%state
    enddo
    do pass=1,3
      result=tree;result_h=held
      if(pass>1) then
        result=fresh
        do node=1,np
          result_h%fields(node)%values=1.03_dp*held%fields(node)%values
        enddo
      endif
      prefix=.125_dp;aux_calls=0;raw_calls=0;provider_must_fail=pass==3
      if(pass==1) then
        call serial_source_step(result,result_h,original,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted)
      else if(pass==2) then
        call serial_source_step(result,result_h,reference,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted)
      else
        call serial_source_step(result,result_h,replayed,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted,original)
      endif
      call require(accepted,'fresh native Godunov replay accepted')
      if(pass<3) then
        call require(aux_calls==3.and.raw_calls==3,'original providers run at each actual source stage')
      else
        call require(aux_calls==0.and.raw_calls==0.and.replayed%replayed,'cached sources never reinvoke providers')
      endif
      if(pass==2) then
        expected=result;expected_h=result_h;ef=f;eh=h;egr=gr;ehr=hr
        eb=b;ew=w;ec=c;ev=v;ea=a;etheta=theta;expected_work=work;expected_prefix=prefix
        call require(work/=original%pressure_work_change,'fresh native pressure work changes with fresh state')
      endif
      if(pass/=3) cycle
      equal=prefix==expected_prefix.and.work==expected_work.and.all(b==eb).and.all(w==ew).and. &
        all(c==ec).and.all(v==ev).and.all(a==ea).and.theta==etheta
      do node=1,np
        equal=equal.and.all(result%nodes(node)%state==expected%nodes(node)%state).and. &
          all(result%nodes(node)%temperature==expected%nodes(node)%temperature).and. &
          all(result_h%fields(node)%values==expected_h%fields(node)%values)
        if(topology%nodes(node)%level/=1) cycle
        equal=equal.and.equal_flux(f(node),ef(node)).and.equal_flux(h(node),eh(node)).and. &
          equal_register(gr(node),egr(node)).and.equal_register(hr(node),ehr(node))
      enddo
      call require(equal,'serial cached source replay equals fresh live native computation')
    enddo
    result=tree;result_h=held;prefix=.125_dp;provider_must_fail=.false.
    call serial_source_step(result,result_h,original,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted, &
      declared_source=.false.)
    call require(accepted,'physical provider archive')
    provider_must_fail=.true.;aux_calls=0;raw_calls=0
    call serial_source_step(result,result_h,replayed,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted, &
      replay=original,declared_source=.false.,refused=refused)
    call require(.not.accepted.and.refused.and.aux_calls==0.and.raw_calls==1, &
      'physical provider reevaluates fresh state and refusal never retries')
    result=tree;result_h=held;prefix=.125_dp;provider_must_fail=.false.;reactor_calls=0
    call serial_source_step(result,result_h,replayed,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted, &
      evaluate_reaction=opaque_reactor,refused=refused)
    call require(.not.accepted.and.refused.and.reactor_calls==1,'native-admitted opaque Godunov reactor refuses once')
    if(nranks==1) return
    do owner_layout=1,2
      owners=[0,1,2,0,1];if(owner_layout==2) owners=[1,0,1,0,1]
      provider_must_fail=.false.;aux_calls=0;raw_calls=0
      call distribute();prefix=.125_dp
      call mpi_source_step(local_original,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused)
      call require(accepted.and.aux_calls==3.and.raw_calls==3,'owner-local native provider capture')
      call distribute()
      do node=1,np
        if(owners(node)/=rank) cycle
        state%fields(node)%values=fresh%nodes(node)%state
        local_held%fields(node)%values=1.03_dp*held%fields(node)%values
      enddo
      provider_must_fail=.true.;aux_calls=0;raw_calls=0;prefix=.125_dp
      call mpi_source_step(local_result,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused,local_original)
      call require(accepted.and.aux_calls==0.and.raw_calls==0.and.local_result%replayed,'MPI source replay')
      equal=abs(work-expected_work)<=8*epsilon(1._dp)*max(tiny(1._dp),abs(expected_work)).and. &
        abs(prefix-expected_prefix)<=8*epsilon(1._dp)*max(tiny(1._dp),abs(expected_prefix))
      do node=1,np
        if(owners(node)/=rank) cycle
        equal=equal.and.all(state%fields(node)%values==expected%nodes(node)%state).and. &
          all(temperature%fields(node)%values(1,:,:,:)==expected%nodes(node)%temperature).and. &
          all(local_held%fields(node)%values==expected_h%fields(node)%values)
        if(topology%nodes(node)%level/=1) cycle
        equal=equal.and.equal_flux(lf(node),ef(node)).and.equal_flux(lh(node),eh(node)).and. &
          equal_register(lgr(node),egr(node)).and.equal_register(lhr(node),ehr(node))
      enddo
      call require(equal,'MPI fresh replay gas/T/UEINT/faces/registers equal independent live serial result')
      saved_gas=state;saved_temp=temperature;saved_h=local_held;saved_work=work;saved_prefix=prefix
      do invalid=1,2
        bad=local_original
        if(rank==owners(2)) then
          if(invalid==1) bad%sources(2)%auxiliary%fields(2)%values(nv+1,1,1,1)=nan_value
          if(invalid==2) bad%sources(2)%provider_faces(2)%x(2,0,1,1)=nan_value
        endif
        call mpi_source_step(local_result,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused,bad)
        call require(.not.accepted.and..not.refused.and.aux_calls==0.and.raw_calls==0,'late cached NaN refuses preflight')
        call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
          equal_state(local_held,saved_h).and.work==saved_work.and.prefix==saved_prefix, &
          'malformed cache preserves accepted state and scalar ledgers')
      enddo
      provider_must_fail=.false.;call distribute();prefix=.125_dp
      call mpi_source_step(local_original,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused, &
        declared_source=.false.)
      call require(accepted,'MPI physical provider archive')
      saved_gas=state;saved_temp=temperature;saved_h=local_held;saved_work=work;saved_prefix=prefix
      provider_must_fail=.true.;aux_calls=0;raw_calls=0
      call mpi_source_step(local_result,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused, &
        replay=local_original,declared_source=.false.)
      call require(.not.accepted.and.refused.and.aux_calls==0.and.raw_calls==1, &
        'MPI physical provider reevaluates once and collectively refuses')
      call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
        equal_state(local_held,saved_h).and.work==saved_work.and.prefix==saved_prefix,'physical provider rollback')
      provider_must_fail=.false.;reactor_calls=0
      call mpi_source_step(local_result,b,w,c,v,a,theta,work,prefix,lf,lh,lgr,lhr,accepted,refused, &
        evaluate_reaction=opaque_reactor)
      call require(.not.accepted.and.refused,'MPI native-admitted opaque reactor collectively refuses')
      call require(reactor_calls<=1,'opaque reactor never retries a refused cell')
      call require(equal_state(state,saved_gas).and.equal_state(temperature,saved_temp).and. &
        equal_state(local_held,saved_h).and.work==saved_work.and.prefix==saved_prefix,'opaque reactor rollback')
    enddo
    if(rank==0) print *, 'PASS: fresh Godunov source replay, exact live serial/MPI state and faces; late cache rollback'
  end subroutine check_source_replay
  subroutine serial_source_step(result,result_h,record,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted,replay, &
      declared_source,evaluate_reaction,refused,partial,recovery,recovery_ready)
    type(eb_reactive_tree_3d), intent(inout) :: result
    type(amr_tree_state_3d), intent(inout) :: result_h
    type(eb_godunov_interval_3d), intent(inout) :: record
    type(eb_godunov_interval_3d), optional, intent(in) :: replay
    real(dp), intent(inout) :: b(:),w(:),c(:),v(:),a(:),theta,work,prefix
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: f(:),h(:)
    type(state_redist_register_3d), allocatable, intent(inout) :: gr(:),hr(:)
    logical, optional, intent(out) :: refused
    logical, intent(out) :: accepted
    logical, optional, intent(in) :: declared_source,partial
    type(eb_godunov_interval_3d), optional, intent(inout) :: recovery
    logical, optional, intent(out) :: recovery_ready
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    logical :: is_source
    type(elementary_reaction), allocatable :: reactions(:)
    type(gas_transport_species), allocatable :: transport(:)
    is_source=.true.
    if(present(declared_source)) is_source=declared_source
    call load_h2o2_elementary_mechanism(reactions,accepted)
    call load_h2o2_elementary_transport(transport,accepted)
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1e-8_dp,2,.true., &
      [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,result,extended,faces, &
      b,w,c,v,theta,accepted,redistribution='state_redist',target_depth=1,integrated_flux=f, &
      state_redist_registers=gr,internal_energy=result_h,internal_energy_flux=h,held_state_redist_registers=hr, &
      internal_energy_pressure_work_change=work,energy_correction=prefix,spatial_flux_weight=.375_dp, &
      auxiliary_source=serial_aux,auxiliary_change=a,auxiliary_has_internal_energy=.true., &
      evaluate_spatial_diffusion=serial_raw,spatial_diffusion_is_source=is_source, &
      godunov_interval_history=record,godunov_replay_interval=replay,godunov_replay_partial=partial, &
      godunov_source_recovery=recovery,godunov_source_recovery_ready=recovery_ready, &
      evaluate_reaction=evaluate_reaction,provider_refused=refused)
  end subroutine serial_source_step
  subroutine mpi_source_step(record,b,w,c,v,a,theta,work,prefix,f,h,gr,hr,accepted,refused,replay, &
      declared_source,evaluate_reaction, &
      partial,recovery,recovery_ready)
    type(eb_godunov_interval_3d), intent(inout) :: record
    type(eb_godunov_interval_3d), optional, intent(in) :: replay
    real(dp), intent(inout) :: b(:),w(:),c(:),v(:),a(:),theta,work,prefix
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: f(:),h(:)
    type(state_redist_register_3d), allocatable, intent(inout) :: gr(:),hr(:)
    logical, intent(out) :: accepted,refused
    logical, optional, intent(in) :: declared_source,partial
    type(eb_godunov_interval_3d), optional, intent(inout) :: recovery
    logical, optional, intent(out) :: recovery_ready
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    logical :: is_source
    type(elementary_reaction), allocatable :: reactions(:)
    type(gas_transport_species), allocatable :: transport(:)
    is_source=.true.
    if(present(declared_source)) is_source=declared_source
    call load_h2o2_elementary_mechanism(reactions,accepted)
    call load_h2o2_elementary_transport(transport,accepted)
    call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temperature,extended, &
      1e-8_dp,2,.true.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,faces, &
      b,w,c,v,theta,MPI_COMM_WORLD,accepted,redistribution='state_redist',target_depth=1,integrated_flux=f, &
      state_redist_registers=gr,internal_energy=local_held,internal_energy_flux=h,held_state_redist_registers=hr, &
      internal_energy_pressure_work_change=work,energy_correction=prefix,spatial_flux_weight=.375_dp, &
      auxiliary_source=mpi_aux,auxiliary_change=a,auxiliary_has_internal_energy=.true., &
      evaluate_spatial_diffusion=mpi_raw,spatial_diffusion_is_source=is_source,provider_refused=refused, &
      godunov_interval_history=record,godunov_replay_interval=replay,godunov_replay_partial=partial, &
      godunov_source_recovery=recovery,godunov_source_recovery_ready=recovery_ready, &
      evaluate_reaction=evaluate_reaction)
  end subroutine mpi_source_step
  subroutine opaque_reactor(species,reactions,energy,source_species,source_energy,duration,rtol,atol, &
      density,temperature,steps,ok,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: energy,source_species(:),source_energy,duration,rtol,atol
    real(dp), intent(inout) :: density(:),temperature
    integer, intent(out) :: steps
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    reactor_calls=reactor_calls+1;steps=0;ok=.false.
  end subroutine opaque_reactor
  subroutine fill_aux(rate,depth,ok)
    type(amr_tree_state_3d), intent(inout) :: rate
    integer, intent(in) :: depth
    logical, intent(out) :: ok
    integer :: node,index
    aux_calls=aux_calls+1;ok=.not.provider_must_fail
    if(prefix_forbid_parent.and.depth==1) ok=.false.
    if(.not.ok) return
    do node=1,np
      if(.not.allocated(rate%fields(node)%values)) cycle
      rate%fields(node)%values=0
      if(topology%nodes(node)%level/=depth.and..not.routed_fixture) cycle
      where(geometry(node)%volume_fraction>0)
        rate%fields(node)%values(5,:,:,:)=.03125_dp
        rate%fields(node)%values(nv+1,:,:,:)=.03125_dp
      endwhere
      if(stage_sink.and..not.packet_test.and.(stage_sink_depth<0.or.depth==stage_sink_depth).and. &
        topology%nodes(node)%level==depth) then
        where(geometry(node)%volume_fraction>0) rate%fields(node)%values(1,:,:,:)=-1e16_dp
        do index=1,size(species)
          if(species(index)%name/='N2') cycle
          where(geometry(node)%volume_fraction>0) rate%fields(node)%values(index+5,:,:,:)=-1e16_dp
        enddo
      endif
      if(recovery_fixture.and.depth==2.and.topology%nodes(node)%level==1) then
        where(geometry(node)%volume_fraction>0) rate%fields(node)%values(1,:,:,:)=-1e16_dp
        do index=1,size(species)
          if(species(index)%name/='N2') cycle
          where(geometry(node)%volume_fraction>0) rate%fields(node)%values(index+5,:,:,:)=-1e16_dp
        enddo
      endif
    enddo
  end subroutine fill_aux
  subroutine serial_aux(gas,interval,source_stage,depth,rate,ok)
    type(eb_reactive_tree_3d), intent(in) :: gas
    real(dp), intent(in) :: interval
    integer, intent(in) :: source_stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: ok
    call fill_aux(rate,depth,ok)
  end subroutine serial_aux
  subroutine mpi_aux(gas,temp,interval,source_stage,depth,rate,comm,ok)
    type(amr_tree_state_3d), intent(in) :: gas,temp
    real(dp), intent(in) :: interval
    integer, intent(in) :: source_stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    call fill_aux(rate,depth,ok)
  end subroutine mpi_aux
  subroutine fill_raw(rate,halos,flux,depth,limiter,ok)
    type(amr_tree_state_3d), intent(inout) :: rate,halos
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    integer, intent(in) :: depth
    real(dp), intent(inout) :: limiter
    logical, intent(out) :: ok
    integer :: node,index
    raw_calls=raw_calls+1;ok=.not.(provider_must_fail.or.raw_refusal)
    if(prefix_forbid_parent.and.depth==1) ok=.false.
    if(.not.ok) return
    if(allocated(halos%fields)) deallocate(halos%fields)
    if(allocated(flux)) deallocate(flux)
    allocate(halos%fields(np),flux(np));limiter=1
    do node=1,np
      if(.not.allocated(rate%fields(node)%values)) cycle
      rate%fields(node)%values=0
      if(topology%nodes(node)%level/=depth) cycle
      where(geometry(node)%volume_fraction>0) rate%fields(node)%values(2,:,:,:)=.003_dp
      if(packet_test.and.stage_sink.and.(stage_sink_depth<0.or.depth==stage_sink_depth)) then
        ! Each declared component is finite. Native density posting refuses
        ! their species sum before performing an overflowing addition.
        where(geometry(node)%volume_fraction>0)
          rate%fields(node)%values(6,:,:,:)=.75_dp*huge(1._dp)
          rate%fields(node)%values(7,:,:,:)=.75_dp*huge(1._dp)
        endwhere
      endif
      if(raw_malformed) rate%fields(node)%values(1,1,1,1)=nan_value
      allocate(halos%fields(node)%values(nv,n+12,n+12,n+12),source=0._dp)
      halos%fields(node)%values(2,:,:,:)=.003_dp
      halos%fields(node)%values(:,7:n+6,7:n+6,7:n+6)=rate%fields(node)%values
      allocate(flux(node)%x(nv,0:n,n,n),flux(node)%y(nv,n,0:n,n),flux(node)%z(nv,n,n,0:n))
      flux(node)%x=0;flux(node)%y=0;flux(node)%z=0
    enddo
  end subroutine fill_raw
  subroutine serial_raw(gas,padded,old_halos,interval,depth,rate,halos,flux,limiter,ok,frame)
    type(eb_reactive_tree_3d), intent(in) :: gas
    type(eb_geometry_3d), intent(in) :: padded(:)
    type(amr_tree_state_3d), intent(in) :: old_halos
    real(dp), intent(in) :: interval
    integer, intent(in) :: depth
    type(amr_tree_state_3d), intent(inout) :: rate,halos
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: limiter
    logical, intent(out) :: ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    call fill_raw(rate,halos,flux,depth,limiter,ok)
  end subroutine serial_raw
  subroutine mpi_raw(gas,temp,padded,old_halos,interval,depth,rate,halos,flux,limiter,comm,ok,frame)
    type(amr_tree_state_3d), intent(in) :: gas,temp,old_halos
    type(eb_geometry_3d), intent(in) :: padded(:)
    real(dp), intent(in) :: interval
    integer, intent(in) :: depth
    type(amr_tree_state_3d), intent(inout) :: rate,halos
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: limiter
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    call fill_raw(rate,halos,flux,depth,limiter,ok)
  end subroutine mpi_raw
  subroutine check_common_step()
    type(eb_godunov_interval_3d), allocatable :: ref_record,local_record,saved_record
    type(eb_godunov_interval_3d), allocatable :: ref_clocks(:),local_clocks(:),saved_clocks(:)
    type(elementary_reaction), allocatable :: reactions(:)
    type(gas_transport_species), allocatable :: transport(:)
    type(eb_reactive_tree_3d) :: reference,ledger
    type(amr_tree_state_3d) :: ref_held,saved_gas,saved_temp,saved_held
    type(eb_integrated_flux_3d), allocatable :: ref_faces(:),local_faces(:),ref_hfaces(:),local_hfaces(:), &
      saved_faces(:),saved_gas_faces(:)
    type(state_redist_register_3d), allocatable :: ref_regs(:),local_regs(:),ref_hregs(:),local_hregs(:), &
      saved_regs(:),saved_hregs(:)
    integer, allocatable :: depth
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),sb(nv),sw(nv),sch(nv),sext(nv),theta,stheta
    integer(int64) :: counts(3),ref_counts(3),saved_counts(3)
    real(dp) :: saved_b(nv),saved_w(nv),saved_ch(nv),saved_ext(nv),saved_theta
    real(dp) :: work,ref_work,saved_work,v0(1),v1(1),adv(1),residual,volume
    real(dp), allocatable :: ref_prefix,local_prefix
    type(dual_energy_policy), allocatable :: clock_policy
    real(dp) :: raw_residual,reset_total,projection_total
    real(dp) :: saved_prefix
    real(dp) :: clock_work
    integer :: node,x,y,z,owner_layout,invalid,r
    integer, parameter :: clock_depths(7)=[0,1,2,2,1,2,2],clock_parents(7)=[0,1,2,2,1,5,5]
    integer, parameter :: clock_last(7)=[7,4,3,4,7,6,7]
    real(dp), parameter :: clock_times(7)=1e-8_dp*[0._dp,0._dp,0._dp,.25_dp,.5_dp,.5_dp,.75_dp]
    type(eb_reactive_tree_3d), save :: native_result
    type(amr_tree_state_3d), save :: native_h
    type(eb_integrated_flux_3d), allocatable, save :: native_faces(:),native_hfaces(:)
    type(state_redist_register_3d), allocatable, save :: native_regs(:),native_hregs(:)
    real(dp), allocatable, save :: native_ledger(:)
    real(dp), save :: native_prefix=0
    logical :: accepted,equal,gas_changed,held_changed
    character(len=16) :: method
    if(common_react) then
      call load_h2o2_elementary_mechanism(reactions,accepted)
      call require(accepted,'native reaction database for reactor/reset fixture')
    else
      allocate(reactions(0))
    endif
    call load_h2o2_elementary_transport(transport,accepted);call require(accepted,'native transport database')
    method='none';if(stage==1) method='state_redist'
    call require(.not.recursive_test.or.selection==0,'recursive fixture evolves all levels')
    if(selection==1) then
      allocate(depth);depth=1
      allocate(ref_regs(np),ref_hregs(np),local_regs(np),local_hregs(np))
    endif
    if(common_react.or.scalar_clock_reset_test) then
      if(scalar_clock_reset_test) then
        allocate(clock_policy);clock_policy%eta2=0;clock_policy%allow_negative_energy=.false.
      endif
      allocate(ref_prefix,local_prefix);ref_prefix=.125_dp;local_prefix=.125_dp
    endif
    if(.not.recursive_test) allocate(ref_record,local_record)
    reference=tree;ref_held=held;ref_work=huge(1._dp)
    if(rank==0) print *, 'Godunov held step (recursive, geometry, selected, StateRedist):',recursive_test,mode,selection,stage
    if(selection==1.and.stage==1) then
      if(recursive_test) then
        call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1e-8_dp,2,.false., &
          [common_transport,common_transport,common_transport,common_transport], &
          [0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,reference,extended,faces, &
          sb,sw,sch,sext,stheta,accepted,level_steps=ref_counts,redistribution=method,integrated_flux=ref_faces, &
          internal_energy=ref_held,internal_energy_flux=ref_hfaces, &
          internal_energy_pressure_work_change=ref_work, &
            godunov_interval_history=ref_clocks, &
            energy_policy=clock_policy,energy_correction=ref_prefix)
      else
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1e-8_dp,2,common_react, &
          [common_transport,common_transport,common_transport,common_transport], &
          [0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,reference,extended,faces, &
          sb,sw,sch,sext,stheta,accepted,redistribution=method,target_depth=depth,integrated_flux=ref_faces, &
          state_redist_registers=ref_regs,internal_energy=ref_held,internal_energy_flux=ref_hfaces, &
          held_state_redist_registers=ref_hregs,internal_energy_pressure_work_change=ref_work, &
            spatial_flux_weight=selected_weight,energy_correction=ref_prefix, &
            godunov_interval_history=ref_record)
      endif
    else
      if(recursive_test) then
        call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1e-8_dp,2,.false., &
          [common_transport,common_transport,common_transport,common_transport], &
          [0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,reference,extended,faces, &
          sb,sw,sch,sext,stheta,accepted,level_steps=ref_counts,redistribution=method,integrated_flux=ref_faces, &
          internal_energy=ref_held,internal_energy_flux=ref_hfaces, &
          internal_energy_pressure_work_change=ref_work, &
            godunov_interval_history=ref_clocks, &
            energy_policy=clock_policy,energy_correction=ref_prefix)
      else
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1e-8_dp,2,common_react, &
          [common_transport,common_transport,common_transport,common_transport], &
          [0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,reference,extended,faces, &
          sb,sw,sch,sext,stheta,accepted,redistribution=method,target_depth=depth,integrated_flux=ref_faces, &
          internal_energy=ref_held,internal_energy_flux=ref_hfaces, &
          internal_energy_pressure_work_change=ref_work, &
            spatial_flux_weight=selected_weight,energy_correction=ref_prefix, &
            godunov_interval_history=ref_record)
      endif
    endif
    call require(accepted,'serial actual held Godunov time step')
    if(weight_test.and..not.history_test) then
      if(weight_case==0) then
        native_result=reference;native_h=ref_held;native_faces=ref_faces;native_hfaces=ref_hfaces
        native_ledger=[sb,sw,sch,sext,stheta,ref_work]
        if(allocated(ref_prefix)) native_prefix=ref_prefix
        if(selection==1.and.stage==1) then
          native_regs=ref_regs;native_hregs=ref_hregs
        endif
      else
        equal=.true.;gas_changed=.false.;held_changed=.false.
        do node=1,np
          if(weight_case==1) equal=equal.and.all(reference%nodes(node)%state==native_result%nodes(node)%state).and. &
            all(reference%nodes(node)%temperature==native_result%nodes(node)%temperature).and. &
            all(ref_held%fields(node)%values==native_h%fields(node)%values)
          if(weight_case==1) equal=equal.and.equal_flux(ref_faces(node),native_faces(node)).and. &
            equal_flux(ref_hfaces(node),native_hfaces(node))
          if(weight_case==1.and.selection==1.and.stage==1) equal=equal.and. &
            equal_register(ref_regs(node),native_regs(node)).and.equal_register(ref_hregs(node),native_hregs(node))
          gas_changed=gas_changed.or..not.equal_flux(ref_faces(node),native_faces(node))
          held_changed=held_changed.or..not.equal_flux(ref_hfaces(node),native_hfaces(node))
        enddo
        if(weight_case==1) then
          equal=equal.and.all([sb,sw,sch,sext,stheta,ref_work]==native_ledger)
          if(allocated(ref_prefix)) equal=equal.and.ref_prefix==native_prefix
          call require(equal,'SDC weight-one preserves exact native gas/T/UEINT, faces, registers and all ledgers')
        endif
        if(weight_case==2) call require(gas_changed.and.held_changed,'SDC zero coefficient actually changes both face families')
      endif
    endif
    if(recursive_test) call require(all(ref_counts==[1_int64,2_int64,4_int64]),'serial recursive clocks [1,2,4]')
    call require(abs(ref_work)>0,'actual nonzero Godunov pressure work ledger')
    if(recursive_test) then
      do node=1,np
        if(topology%nodes(node)%level/=2) cycle
        call require(any(ref_hfaces(node)%x/=0).or.any(ref_hfaces(node)%y/=0).or.any(ref_hfaces(node)%z/=0), &
          'complete grandchild scalar interval flux reaches outer caller')
        call require(any(ref_faces(node)%x/=0).or.any(ref_faces(node)%y/=0).or.any(ref_faces(node)%z/=0), &
          'complete grandchild gas interval flux reaches outer caller')
      enddo
    endif
    ! Independently integrate the complete native hierarchy face family,
    ! including coarse/fine masking and every descendant interval.
    if(selection==0.and..not.common_react) then
      ledger=tree
      do node=1,np
        ledger%nodes(node)%state=held%fields(node)%values
      enddo
      call integrate_reactive_eb_tree_3d(ledger,v0,accepted);call require(accepted,'old scalar composite')
      do node=1,np
        ledger%nodes(node)%state=ref_held%fields(node)%values
      enddo
      call integrate_reactive_eb_tree_3d(ledger,v1,accepted);call require(accepted,'new scalar composite')
      do node=1,np
        ledger%nodes(node)%state=0
        associate(g=>geometry(node),f=>ref_hfaces(node))
        volume=g%dx*g%dy*g%dz
        do z=1,n
        do y=1,n
        do x=1,n
          if(g%volume_fraction(x,y,z)<=0) cycle
          ledger%nodes(node)%state(1,x,y,z)=( &
            g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(1,x-1,y,z)-g%x_face_fraction(x,y,z)*f%x(1,x,y,z))+ &
            g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(1,x,y-1,z)-g%y_face_fraction(x,y,z)*f%y(1,x,y,z))+ &
            g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(1,x,y,z-1)-g%z_face_fraction(x,y,z)*f%z(1,x,y,z)))/ &
            (g%volume_fraction(x,y,z)*volume)
        enddo
        enddo
        enddo
        end associate
      enddo
      call integrate_reactive_eb_tree_3d(ledger,adv,accepted);call require(accepted,'scalar face composite')
      raw_residual=abs(v1(1)-v0(1)-adv(1)-ref_work)/max(1._dp,abs(v0(1)))
      projection_total=0;reset_total=0
      if(recursive_test) then
        call require(all(ref_clocks%held_clock_events_closed).and.all(ref_clocks%held_clock_prefixes_valid), &
          'complete accepted native clock forest has measured and valid scalar entry/exit budgets')
        projection_total=ref_clocks(1)%held_clock_projection_after
        reset_total=ref_clocks(1)%held_clock_reset_after
      endif
      residual=abs(v1(1)-v0(1)-adv(1)-ref_work-projection_total-reset_total)/max(1._dp,abs(v0(1)))
      if(scalar_clock_reset_test) then
        call require(raw_residual>3e-14_dp,'actual native three-level reset changes the independent H inventory')
        if(rank==0) print *, 'NATIVE_CLOCK_H raw, closed, projection, reset:', &
          raw_residual,residual,projection_total,reset_total
      endif
      if(rank==0) then
        print *, 'Godunov actual held step inventory residual:',residual
        print *, 'Actual held initial/final/advection/pressure-work:',v0(1),v1(1),adv(1),ref_work
        flush(output_unit)
      endif
      call MPI_Barrier(MPI_COMM_WORLD,status)
      call require(residual<3e-14_dp,'independent UEINT advection and pressure work step inventory')
      max_residual=max(max_residual,residual)
    endif
    if(nranks==1) return
    do owner_layout=1,2
      owners=[0,0,1,0,1];if(owner_layout==2) owners=[1,0,1,1,0]
      call distribute();work=huge(1._dp)
      if(allocated(local_prefix)) local_prefix=.125_dp
      if(selection==1.and.stage==1) then
        if(recursive_test) then
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,state,temperature, &
            1e-8_dp,2,.false., &
            [common_transport,common_transport,common_transport,common_transport], &
            [0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,extended,faces, &
            b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,level_steps=counts,redistribution=method,integrated_flux=local_faces, &
            internal_energy=local_held,internal_energy_flux=local_hfaces, &
            internal_energy_pressure_work_change=work, &
            godunov_interval_history=local_clocks, &
            energy_policy=clock_policy,energy_correction=local_prefix)
        else
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temperature,extended, &
            1e-8_dp,2,common_react, &
            [common_transport,common_transport,common_transport,common_transport],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,faces, &
            b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,redistribution=method,target_depth=depth,integrated_flux=local_faces, &
            state_redist_registers=local_regs,internal_energy=local_held,internal_energy_flux=local_hfaces, &
            held_state_redist_registers=local_hregs,internal_energy_pressure_work_change=work, &
            spatial_flux_weight=selected_weight,energy_correction=local_prefix, &
            godunov_interval_history=local_record)
        endif
      else
        if(recursive_test) then
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,state,temperature, &
            1e-8_dp,2,.false., &
            [common_transport,common_transport,common_transport,common_transport], &
            [0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,extended,faces, &
            b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,level_steps=counts,redistribution=method,integrated_flux=local_faces, &
            internal_energy=local_held,internal_energy_flux=local_hfaces, &
            internal_energy_pressure_work_change=work, &
            godunov_interval_history=local_clocks, &
            energy_policy=clock_policy,energy_correction=local_prefix)
        else
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temperature,extended, &
            1e-8_dp,2,common_react, &
            [common_transport,common_transport,common_transport,common_transport],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,faces, &
            b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,redistribution=method,target_depth=depth,integrated_flux=local_faces, &
            internal_energy=local_held,internal_energy_flux=local_hfaces, &
            internal_energy_pressure_work_change=work, &
            spatial_flux_weight=selected_weight,energy_correction=local_prefix, &
            godunov_interval_history=local_record)
        endif
      endif
      call require(accepted,'MPI actual held Godunov time step')
      if(recursive_test) call require(all(counts==ref_counts),'MPI recursive clocks [1,2,4]')
      ! Field and face payloads must agree exactly. Composite sums use different
      ! serial/MPI reduction ordering; bound only their few-ULP rounding.
      equal=abs(work-ref_work)<=8*epsilon(1._dp)*max(tiny(1._dp),abs(ref_work)).and. &
        all(b==sb).and.all(abs(w-sw)<=8*epsilon(1._dp)*max(tiny(1._dp),abs(sw))).and. &
        all(ch==sch).and.all(ext==sext).and.theta==stheta
      do node=1,np
        if(owners(node)/=rank) cycle
        equal=equal.and.all(state%fields(node)%values==reference%nodes(node)%state).and. &
          all(temperature%fields(node)%values(1,:,:,:)==reference%nodes(node)%temperature).and. &
          all(local_held%fields(node)%values==ref_held%fields(node)%values)
        if(selection==1.and.topology%nodes(node)%level/=1) cycle
        equal=equal.and.equal_flux(local_faces(node),ref_faces(node)).and.equal_flux(local_hfaces(node),ref_hfaces(node))
        if(selection==1.and.stage==1) equal=equal.and.equal_register(local_regs(node),ref_regs(node)).and. &
          equal_register(local_hregs(node),ref_hregs(node))
      enddo
      if(.not.equal) then
        print *, 'common diagnostic rank/layout/work/b/w/ch/ext/theta:',rank,owner_layout, &
          work-ref_work,maxval(abs(b-sb)),maxval(abs(w-sw)),maxval(abs(ch-sch)),maxval(abs(ext-sext)),theta-stheta
        do node=1,np
          if(owners(node)/=rank) cycle
          print *, 'common node/gas/T/held:',node,maxval(abs(state%fields(node)%values-reference%nodes(node)%state)), &
            maxval(abs(temperature%fields(node)%values(1,:,:,:)-reference%nodes(node)%temperature)), &
            maxval(abs(local_held%fields(node)%values-ref_held%fields(node)%values))
        enddo
      endif
      if(common_react) equal=equal.and.abs(local_prefix-ref_prefix)<= &
        8*epsilon(1._dp)*max(tiny(1._dp),abs(ref_prefix))
      call require(equal,'exact serial/MPI common fields/faces/registers and bounded ledger reductions')
      if(recursive_test) then
        call require(allocated(ref_clocks).and.allocated(local_clocks),'whole recursive clock archive publishes')
        call require(size(local_clocks)==7.and.size(ref_clocks)==7,'complete [1,2,4] preorder clock family')
        clock_work=0
        do r=1,7
          equal=record_equal(local_clocks(r),ref_clocks(r))
          equal=equal.and.local_clocks(r)%depth==clock_depths(r).and. &
            local_clocks(r)%parent_record==clock_parents(r).and.local_clocks(r)%subtree_last==clock_last(r)
          equal=equal.and.abs(local_clocks(r)%time-clock_times(r))<=32*epsilon(1._dp)*1e-8_dp.and. &
            local_clocks(r)%interval==1e-8_dp/2._dp**clock_depths(r)
          equal=equal.and.local_clocks(r)%clock_pressure_work_change==clock_work
          clock_work=clock_work+local_clocks(r)%pressure_work_change
          if(clock_depths(r)>0.and.stage==1) equal=equal.and. &
            allocated(local_clocks(r)%open_ancestor_registers).and. &
            allocated(local_clocks(r)%open_ancestor_scalar_registers)
          call require(equal,'exact native clock stage payload, parent/range/time and pressure prefix')
        enddo
        call require(clock_work==work,'native clock pressure-work sum equals accepted outer diagnostic')
        call require(all(local_clocks%held_clock_events_closed).and.all(local_clocks%held_clock_prefixes_valid), &
          'owner/idle native clock scalar forest is complete before publication')
        do r=1,7
          projection_total=local_clocks(r)%held_clock_projection_before+ &
            sum(local_clocks(r:clock_last(r))%internal_energy_projection_change)- &
            sum(local_clocks(r:clock_last(r))%initial_internal_energy_change)
          reset_total=local_clocks(r)%held_clock_reset_before+ &
            sum(local_clocks(r:clock_last(r))%internal_energy_closure_change)+ &
            sum(local_clocks(r:clock_last(r))%held_clock_end_change)
          if(r<clock_last(r)) reset_total=reset_total+sum(local_clocks(r+1:clock_last(r))%held_clock_start_change)
          call require(abs(projection_total-local_clocks(r)%held_clock_projection_after)<3e-14_dp.and. &
            abs(reset_total-local_clocks(r)%held_clock_reset_after)<3e-14_dp, &
            'actual native subtree prefix contains its accepted descendants and delayed parent reset once')
          call require(maxval(abs([local_clocks(r)%held_clock_start_change,local_clocks(r)%held_clock_end_change, &
            local_clocks(r)%held_clock_projection_before,local_clocks(r)%held_clock_reset_before, &
            local_clocks(r)%held_clock_projection_after,local_clocks(r)%held_clock_reset_after]- &
            [ref_clocks(r)%held_clock_start_change,ref_clocks(r)%held_clock_end_change, &
            ref_clocks(r)%held_clock_projection_before,ref_clocks(r)%held_clock_reset_before, &
            ref_clocks(r)%held_clock_projection_after,ref_clocks(r)%held_clock_reset_after]))<3e-14_dp, &
            'actual serial/MPI owner/idle scalar clock prefixes agree at the unchanged inventory tolerance')
        enddo

        call require(record_field_equal(local_clocks(1)%closed_gas,state).and. &
          record_field_equal(local_clocks(1)%closed_temperature,temperature).and. &
          record_field_equal(local_clocks(1)%closed_internal_energy,local_held), &
          'root closed history is the accepted synchronized physical state')
        do node=1,np
          if(owners(node)/=rank) cycle
          call require(equal_flux(local_clocks(1)%closed_faces(node),local_faces(node)).and. &
            equal_flux(local_clocks(1)%closed_scalar_faces(node),local_hfaces(node)), &
            'root closed history retains every descendant gas and UEINT face')
        enddo
      endif
      if(.not.recursive_test) then
        call require(allocated(local_record).and.allocated(ref_record),'accepted Godunov clock archives publish')
        call require(record_equal(local_record,ref_record),'exact owner-local source/hydro stage histories and bounded ledgers')
        call require(local_record%pressure_work_change==work,'archive work is the accepted native ledger')
        call require(all(local_record%boundary_change==b).and.all(local_record%wall_change==w).and. &
          all(local_record%chemical_change==ch).and.all(local_record%source_change==ext), &
          'archive retains actual accepted volume and boundary diagnostics')
        if(history_test) then
          call require(lbound(local_record%sources,1)==0.and.ubound(local_record%sources,1)==2.and. &
            size(local_record%hydro)==2,'separate source stages 0:2 and native hydro passes 1:2')
          call require(local_record%sources(2)%coarse_registers_accumulated.and. &
            local_record%hydro(2)%coarse_registers_accumulated,'archive explicitly identifies cumulative coarse chronology')
        endif
      endif
      if(owner_layout/=2) cycle
      saved_prefix=0
      if(allocated(local_prefix)) saved_prefix=local_prefix
      saved_work=work;saved_b=b;saved_w=w;saved_ch=ch;saved_ext=ext;saved_theta=theta
      saved_gas_faces=local_faces
      if(recursive_test) then
        saved_counts=counts;saved_clocks=local_clocks
      endif
      if(selection==1.and.stage==1) then
        saved_regs=local_regs;saved_hregs=local_hregs
      endif
      if(allocated(local_record)) saved_record=local_record
      saved_gas=state;saved_temp=temperature;saved_held=local_held;saved_faces=local_hfaces
      ! The idle rank omits a ledger: reject the collective transaction before
      ! replacing any prior accepted state, face family, or diagnostic.
      if(rank==2) then
        if(selection==1.and.stage==1) then
          if(recursive_test) then
            call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,state,temperature, &
              1e-8_dp,2,.false., &
              [common_transport,common_transport,common_transport,common_transport], &
              [0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,extended,faces, &
              b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,level_steps=counts,redistribution=method,integrated_flux=local_faces, &
              internal_energy=local_held,internal_energy_flux=local_hfaces, &
            energy_policy=clock_policy,energy_correction=local_prefix)
          else
            call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temperature,extended, &
              1e-8_dp,2,common_react, &
              [common_transport,common_transport,common_transport,common_transport],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,faces, &
              b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,redistribution=method,target_depth=depth,integrated_flux=local_faces, &
              state_redist_registers=local_regs,internal_energy=local_held,internal_energy_flux=local_hfaces, &
              held_state_redist_registers=local_hregs, &
            spatial_flux_weight=selected_weight,energy_correction=local_prefix, &
            godunov_interval_history=local_record)
          endif
        else
          if(recursive_test) then
            call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,state,temperature, &
              1e-8_dp,2,.false., &
              [common_transport,common_transport,common_transport,common_transport], &
              [0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,extended,faces, &
              b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,level_steps=counts,redistribution=method,integrated_flux=local_faces, &
              internal_energy=local_held,internal_energy_flux=local_hfaces, &
            energy_policy=clock_policy,energy_correction=local_prefix)
          else
            call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temperature,extended, &
              1e-8_dp,2,common_react, &
              [common_transport,common_transport,common_transport,common_transport],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,faces, &
              b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,redistribution=method,target_depth=depth,integrated_flux=local_faces, &
              internal_energy=local_held,internal_energy_flux=local_hfaces, &
            spatial_flux_weight=selected_weight,energy_correction=local_prefix, &
            godunov_interval_history=local_record)
          endif
        endif
      else
        if(selection==1.and.stage==1) then
          if(recursive_test) then
            call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,state,temperature, &
              1e-8_dp,2,.false., &
              [common_transport,common_transport,common_transport,common_transport], &
              [0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,extended,faces, &
              b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,level_steps=counts,redistribution=method,integrated_flux=local_faces, &
              internal_energy=local_held,internal_energy_flux=local_hfaces, &
              internal_energy_pressure_work_change=work, &
            godunov_interval_history=local_clocks, &
            energy_policy=clock_policy,energy_correction=local_prefix)
          else
            call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temperature,extended, &
              1e-8_dp,2,common_react, &
              [common_transport,common_transport,common_transport,common_transport],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,faces, &
              b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,redistribution=method,target_depth=depth,integrated_flux=local_faces, &
              state_redist_registers=local_regs,internal_energy=local_held,internal_energy_flux=local_hfaces, &
              held_state_redist_registers=local_hregs,internal_energy_pressure_work_change=work, &
            spatial_flux_weight=selected_weight,energy_correction=local_prefix, &
            godunov_interval_history=local_record)
          endif
        else
          if(recursive_test) then
            call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,state,temperature, &
              1e-8_dp,2,.false., &
              [common_transport,common_transport,common_transport,common_transport], &
              [0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,extended,faces, &
              b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,level_steps=counts,redistribution=method,integrated_flux=local_faces, &
              internal_energy=local_held,internal_energy_flux=local_hfaces, &
              internal_energy_pressure_work_change=work, &
            godunov_interval_history=local_clocks, &
            energy_policy=clock_policy,energy_correction=local_prefix)
          else
            call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temperature,extended, &
              1e-8_dp,2,common_react, &
              [common_transport,common_transport,common_transport,common_transport],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,faces, &
              b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,redistribution=method,target_depth=depth,integrated_flux=local_faces, &
              internal_energy=local_held,internal_energy_flux=local_hfaces, &
              internal_energy_pressure_work_change=work, &
            spatial_flux_weight=selected_weight,energy_correction=local_prefix, &
            godunov_interval_history=local_record)
          endif
        endif
      endif
      call require(.not.accepted,'idle-rank pressure-work family mismatch refuses')
      if(recursive_test) call require(all(counts==saved_counts),'recursive refusal preserves accepted clock counts')
      equal=work==saved_work.and.all(b==saved_b).and.all(w==saved_w).and.all(ch==saved_ch).and. &
        all(ext==saved_ext).and.theta==saved_theta
      do node=1,np
        if(owners(node)/=rank) cycle
        equal=equal.and.all(state%fields(node)%values==saved_gas%fields(node)%values).and. &
          all(temperature%fields(node)%values==saved_temp%fields(node)%values).and. &
          all(local_held%fields(node)%values==saved_held%fields(node)%values)
        if(selection==1.and.topology%nodes(node)%level/=1) cycle
        equal=equal.and.equal_flux(local_hfaces(node),saved_faces(node)).and. &
          equal_flux(local_faces(node),saved_gas_faces(node))
        if(selection==1.and.stage==1) equal=equal.and.equal_register(local_regs(node),saved_regs(node)).and. &
          equal_register(local_hregs(node),saved_hregs(node))
      enddo
      if(allocated(local_prefix)) equal=equal.and.local_prefix==saved_prefix
      if(allocated(saved_record)) equal=equal.and.record_equal(local_record,saved_record).and. &
        local_record%pressure_work_change==saved_record%pressure_work_change
      if(recursive_test) then
        do r=1,7
          equal=equal.and.record_equal(local_clocks(r),saved_clocks(r))
        enddo
      endif
      call require(equal,'common step collective refusal preserves prior publication')
      if(.not.weight_test.or.weight_case/=3) cycle
      do invalid=1,6
        if(rank==2) then
          select case(invalid)
          case(1)
            selected_weight=.5_dp
          case(2)
            selected_weight=nan_value
          case(3)
            selected_weight=-.1_dp
          case(4)
            selected_weight=1.1_dp
          case(5)
            deallocate(selected_weight)
          end select
        endif
        if(invalid==6.and.rank==2) then
        if(selection==1.and.stage==1) then
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temperature,extended, &
            1e-8_dp,2,common_react, &
            [common_transport,common_transport,common_transport,common_transport],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,faces, &
            b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,redistribution=method,target_depth=depth,integrated_flux=local_faces, &
            state_redist_registers=local_regs,internal_energy=local_held,internal_energy_flux=local_hfaces, &
            held_state_redist_registers=local_hregs,internal_energy_pressure_work_change=work, &
            spatial_flux_weight=selected_weight,energy_correction=local_prefix)
        else
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temperature,extended, &
            1e-8_dp,2,common_react, &
            [common_transport,common_transport,common_transport,common_transport],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,faces, &
            b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,redistribution=method,target_depth=depth,integrated_flux=local_faces, &
            internal_energy=local_held,internal_energy_flux=local_hfaces,internal_energy_pressure_work_change=work, &
            spatial_flux_weight=selected_weight,energy_correction=local_prefix)
        endif
        else
        if(selection==1.and.stage==1) then
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temperature,extended, &
            1e-8_dp,2,common_react, &
            [common_transport,common_transport,common_transport,common_transport],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,faces, &
            b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,redistribution=method,target_depth=depth,integrated_flux=local_faces, &
            state_redist_registers=local_regs,internal_energy=local_held,internal_energy_flux=local_hfaces, &
            held_state_redist_registers=local_hregs,internal_energy_pressure_work_change=work, &
            spatial_flux_weight=selected_weight,energy_correction=local_prefix, &
            godunov_interval_history=local_record)
        else
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temperature,extended, &
            1e-8_dp,2,common_react, &
            [common_transport,common_transport,common_transport,common_transport],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,faces, &
            b,w,ch,ext,theta,MPI_COMM_WORLD,accepted,redistribution=method,target_depth=depth,integrated_flux=local_faces, &
            internal_energy=local_held,internal_energy_flux=local_hfaces,internal_energy_pressure_work_change=work, &
            spatial_flux_weight=selected_weight,energy_correction=local_prefix, &
            godunov_interval_history=local_record)
        endif
        endif
        call require(.not.accepted,'SDC collective weight controls refuse before native providers and reset')
        if(.not.allocated(selected_weight)) allocate(selected_weight)
        selected_weight=test_weights(weight_case)
        equal=work==saved_work.and.all(b==saved_b).and.all(w==saved_w).and.all(ch==saved_ch).and. &
          all(ext==saved_ext).and.theta==saved_theta.and.equal_state(state,saved_gas).and. &
          equal_state(temperature,saved_temp).and.equal_state(local_held,saved_held)
        if(allocated(local_prefix)) equal=equal.and.local_prefix==saved_prefix
        do node=1,np
          equal=equal.and.equal_flux(local_hfaces(node),saved_faces(node)).and. &
            equal_flux(local_faces(node),saved_gas_faces(node))
          if(selection==1.and.stage==1) equal=equal.and.equal_register(local_regs(node),saved_regs(node)).and. &
            equal_register(local_hregs(node),saved_hregs(node))
        enddo
        equal=equal.and.record_equal(local_record,saved_record).and. &
          local_record%pressure_work_change==saved_record%pressure_work_change
        call require(equal,'SDC weight refusal preserves state, stages, faces, registers, work, ledgers and reset prefix')
      enddo
    enddo
  end subroutine check_common_step
  subroutine distribute()
    integer :: node
    state%fields=tree_state_fields();temperature%fields=held%fields
    local_held=held;local_predictor=predictor;end_state=state;end_temperature=temperature
    local_held_end=held_end;local_predictor_end=predictor_end
    do node=1,np
      state%fields(node)%values=tree%nodes(node)%state
      temperature%fields(node)%values(1,:,:,:)=tree%nodes(node)%temperature
      end_state%fields(node)%values=ancestor_end%nodes(node)%state
      end_temperature%fields(node)%values(1,:,:,:)=ancestor_end%nodes(node)%temperature
      if(owners(node)==rank) cycle
      deallocate(state%fields(node)%values,temperature%fields(node)%values,local_held%fields(node)%values, &
        local_predictor%fields(node)%values,end_state%fields(node)%values,end_temperature%fields(node)%values, &
        local_held_end%fields(node)%values,local_predictor_end%fields(node)%values)
    enddo
  end subroutine distribute
  function tree_state_fields() result(fields)
    use amr_patch_tree_3d_mod, only: amr_tree_field_3d
    type(amr_tree_field_3d), allocatable :: fields(:)
    integer :: node
    allocate(fields(np))
    do node=1,np
      fields(node)%values=tree%nodes(node)%state
    enddo
  end function tree_state_fields
  subroutine evaluate_serial()
    integer, allocatable :: depth
    type(eb_reactive_tree_3d), allocatable :: first,last
    type(amr_tree_state_3d), allocatable :: source_first,source_last,h_first,h_last
    real(dp), allocatable :: fraction(:)
    if(selection==1) then
      allocate(depth,first,last,source_first,source_last,h_first,h_last,fraction(1))
      depth=1;first=tree;last=ancestor_end;source_first=predictor;source_last=predictor_end
      h_first=held;h_last=held_end;fraction=fractions(:1)
    endif
    if(stage==0) then
      call compute_sdc_eb_tree_hydro_3d(species,tree,predictor,extended,1e-8_dp,ref_gas,ref_gas_flux,good, &
        boundary_faces=faces,baseline_flux=ref_gas_base,spatial_rate_halos=ref_gas_raw,old_state_halos=ref_gas_old, &
        target_depth=depth,ancestor_start=first,ancestor_end=last,ancestor_fraction=fraction, &
        source_start=source_first,source_end=source_last,source_fraction=fraction, &
        internal_energy_state=held,internal_energy_rate=ref_h,internal_energy_flux=ref_h_flux, &
        internal_energy_baseline_flux=ref_h_base,internal_energy_spatial_rate_halos=ref_h_raw, &
        internal_energy_old_state_halos=ref_h_old,ancestor_internal_energy_start=h_first, &
        ancestor_internal_energy_end=h_last,internal_energy_pressure_work=ref_pwork,spatial_flux_weight=selected_weight)
    else
      call compute_godunov_state_redist_eb_tree_hydro_3d(species,tree,predictor,extended,1e-8_dp,ref_gas,ref_gas_flux,good, &
        boundary_faces=faces,baseline_flux=ref_gas_base,spatial_rate_halos=ref_gas_raw,old_state_halos=ref_gas_old, &
        registers=ref_gas_registers,max_order=2,correct_forward_flux=.true.,target_depth=depth, &
        ancestor_start=first,ancestor_end=last,ancestor_fraction=fraction,source_start=source_first, &
        source_end=source_last,source_fraction=fraction,internal_energy_state=held,internal_energy_rate=ref_h, &
        internal_energy_flux=ref_h_flux,internal_energy_registers=ref_h_registers,internal_energy_baseline_flux=ref_h_base, &
        internal_energy_spatial_rate_halos=ref_h_raw,internal_energy_old_state_halos=ref_h_old, &
        ancestor_internal_energy_start=h_first,ancestor_internal_energy_end=h_last,internal_energy_pressure_work=ref_pwork, &
        spatial_flux_weight=selected_weight)
    endif
  end subroutine evaluate_serial
  subroutine evaluate_mpi()
    integer, allocatable :: depth
    type(amr_tree_state_3d), allocatable :: first,last,t_first,t_last,source_first,source_last,h_first,h_last
    real(dp), allocatable :: fraction(:)
    if(selection==1) then
      allocate(depth,first,last,t_first,t_last,source_first,source_last,h_first,h_last,fraction(2))
      depth=1;first=state;last=end_state;t_first=temperature;t_last=end_temperature
      source_first=local_predictor;source_last=local_predictor_end;h_first=local_held;h_last=local_held_end;fraction=fractions
    endif
    if(refusal_mode>0.and.rank==nranks-1) then
      if(refusal_mode==1) then
        call compute_godunov_state_redist_mpi_eb_tree_hydro_3d(species,topology,owners,geometry, &
          state,temperature,local_predictor,extended,1e-8_dp,gas_rate,gas_flux,MPI_COMM_WORLD,good, &
          boundary_faces=faces,baseline_flux=gas_base,spatial_rate_halos=gas_raw,old_state_halos=gas_old, &
          registers=gas_registers,max_order=2,correct_forward_flux=.true.)
      else
        call compute_godunov_state_redist_mpi_eb_tree_hydro_3d(species,topology,owners,geometry, &
          state,temperature,local_predictor,extended,1e-8_dp,gas_rate,gas_flux,MPI_COMM_WORLD,good, &
          boundary_faces=faces,baseline_flux=gas_base,spatial_rate_halos=gas_raw,old_state_halos=gas_old, &
          registers=gas_registers,max_order=2,correct_forward_flux=.true., &
          internal_energy_state=local_held,internal_energy_rate=h_rate,internal_energy_flux=h_flux)
      endif
      return
    endif
    if(stage==0) then
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,geometry,state,temperature,local_predictor, &
        extended,1e-8_dp,gas_rate,gas_flux,MPI_COMM_WORLD,good,boundary_faces=faces,baseline_flux=gas_base, &
        spatial_rate_halos=gas_raw,old_state_halos=gas_old,target_depth=depth,ancestor_start=first,ancestor_end=last, &
        temperature_start=t_first,temperature_end=t_last,ancestor_fraction=fraction,source_start=source_first, &
        source_end=source_last,source_fraction=fraction,internal_energy_state=local_held,internal_energy_rate=h_rate, &
        internal_energy_flux=h_flux,internal_energy_baseline_flux=h_base,internal_energy_spatial_rate_halos=h_raw, &
        internal_energy_old_state_halos=h_old,ancestor_internal_energy_start=h_first,ancestor_internal_energy_end=h_last, &
        internal_energy_pressure_work=pwork,spatial_flux_weight=selected_weight)
    else
      call compute_godunov_state_redist_mpi_eb_tree_hydro_3d(species,topology,owners,geometry,state,temperature,local_predictor, &
        extended,1e-8_dp,gas_rate,gas_flux,MPI_COMM_WORLD,good,boundary_faces=faces,baseline_flux=gas_base, &
        spatial_rate_halos=gas_raw,old_state_halos=gas_old,registers=gas_registers,max_order=2,correct_forward_flux=.true., &
        target_depth=depth,ancestor_start=first,ancestor_end=last,temperature_start=t_first,temperature_end=t_last, &
        ancestor_fraction=fraction,source_start=source_first,source_end=source_last,source_fraction=fraction, &
        internal_energy_state=local_held,internal_energy_rate=h_rate,internal_energy_flux=h_flux, &
        internal_energy_registers=h_registers,internal_energy_baseline_flux=h_base,internal_energy_spatial_rate_halos=h_raw, &
        internal_energy_old_state_halos=h_old,ancestor_internal_energy_start=h_first,ancestor_internal_energy_end=h_last, &
        internal_energy_pressure_work=pwork,spatial_flux_weight=selected_weight)
    endif
  end subroutine evaluate_mpi
  subroutine check_refusal()
    type(amr_tree_state_3d) :: kept_gas,kept_h,kept_raw,kept_old,kept_h_raw,kept_h_old,kept_work,kept_input
    type(eb_integrated_flux_3d), allocatable :: kept_gas_faces(:),kept_h_faces(:),kept_gas_base(:),kept_h_base(:)
    type(state_redist_register_3d), allocatable :: kept_gas_regs(:),kept_h_regs(:)
    integer :: node
    kept_gas=gas_rate;kept_h=h_rate;kept_gas_faces=gas_flux;kept_h_faces=h_flux
    kept_raw=gas_raw;kept_old=gas_old;kept_h_raw=h_raw;kept_h_old=h_old;kept_work=pwork
    kept_gas_base=gas_base;kept_h_base=h_base;kept_gas_regs=gas_registers;kept_h_regs=h_registers;kept_input=local_held
    do failure=1,4
      refusal_mode=0
      if(failure<=2) refusal_mode=failure
      if(failure==3.and.rank==owners(3)) local_held%fields(3)%values(1,n,n,n)=nan_value
      if(failure==4.and.rank==owners(3)) then
        deallocate(local_held%fields(3)%values)
        allocate(local_held%fields(3)%values(1,0:n-1,n,n))
        local_held%fields(3)%values=kept_input%fields(3)%values
      endif
      call evaluate_mpi()
      call require(.not.good,'idle presence/partial family, late owner NaN and zero-based scalar refusal')
      local_held=kept_input;refusal_mode=0;same=.true.
      do node=1,np
        if(owners(node)/=rank) cycle
        same=same.and.all(gas_rate%fields(node)%values==kept_gas%fields(node)%values)
        same=same.and.all(h_rate%fields(node)%values==kept_h%fields(node)%values)
        same=same.and.all(gas_raw%fields(node)%values==kept_raw%fields(node)%values)
        same=same.and.all(gas_old%fields(node)%values==kept_old%fields(node)%values)
        same=same.and.all(h_raw%fields(node)%values==kept_h_raw%fields(node)%values)
        same=same.and.all(h_old%fields(node)%values==kept_h_old%fields(node)%values)
        same=same.and.all(pwork%fields(node)%values==kept_work%fields(node)%values)
        same=same.and.equal_flux(gas_flux(node),kept_gas_faces(node)).and.equal_flux(h_flux(node),kept_h_faces(node))
        same=same.and.equal_flux(gas_base(node),kept_gas_base(node)).and.equal_flux(h_base(node),kept_h_base(node))
        same=same.and.equal_register(gas_registers(node),kept_gas_regs(node)).and. &
          equal_register(h_registers(node),kept_h_regs(node))
      enddo
      call require(same,'collective scalar refusal preserves both rates, halos, pressure work, face and register families')
    enddo
  end subroutine check_refusal
  subroutine check_weight_refusal()
    type(amr_tree_state_3d) :: kept_gas,kept_h,kept_raw,kept_old,kept_h_raw,kept_h_old,kept_work
    type(eb_integrated_flux_3d), allocatable :: kept_gas_faces(:),kept_h_faces(:),kept_gas_base(:),kept_h_base(:)
    type(state_redist_register_3d), allocatable :: kept_gas_regs(:),kept_h_regs(:)
    integer :: node,invalid
    kept_gas=gas_rate;kept_h=h_rate;kept_gas_faces=gas_flux;kept_h_faces=h_flux
    kept_raw=gas_raw;kept_old=gas_old;kept_h_raw=h_raw;kept_h_old=h_old;kept_work=pwork
    kept_gas_base=gas_base;kept_h_base=h_base
    if(allocated(gas_registers)) kept_gas_regs=gas_registers
    if(allocated(h_registers)) kept_h_regs=h_registers
    do invalid=1,5
      if(rank==nranks-1) then
        select case(invalid)
        case(1)
          selected_weight=.5_dp
        case(2)
          selected_weight=nan_value
        case(3)
          selected_weight=-.1_dp
        case(4)
          selected_weight=1.1_dp
        case(5)
          deallocate(selected_weight)
        end select
      endif
      call evaluate_mpi()
      call require(.not.good,'idle weight disagreement, NaN, bounds and presence jointly refuse')
      if(.not.allocated(selected_weight)) allocate(selected_weight)
      selected_weight=test_weights(weight_case)
      same=equal_state(gas_rate,kept_gas).and.equal_state(h_rate,kept_h).and. &
        equal_state(gas_raw,kept_raw).and.equal_state(gas_old,kept_old).and. &
        equal_state(h_raw,kept_h_raw).and.equal_state(h_old,kept_h_old).and.equal_state(pwork,kept_work)
      same=same.and.(allocated(gas_registers).eqv.allocated(kept_gas_regs)).and. &
        (allocated(h_registers).eqv.allocated(kept_h_regs))
      do node=1,np
        same=same.and.equal_flux(gas_flux(node),kept_gas_faces(node)).and.equal_flux(h_flux(node),kept_h_faces(node))
        same=same.and.equal_flux(gas_base(node),kept_gas_base(node)).and.equal_flux(h_base(node),kept_h_base(node))
        if(allocated(kept_gas_regs)) same=same.and.equal_register(gas_registers(node),kept_gas_regs(node))
        if(allocated(kept_h_regs)) same=same.and.equal_register(h_registers(node),kept_h_regs(node))
      enddo
      call require(same,'weight refusal retains all gas/held fields, raw populations, faces, work and registers')
    enddo
  end subroutine check_weight_refusal
  logical function record_field_equal(a,b) result(same)
    type(amr_tree_state_3d), intent(in) :: a,b
    integer :: node
    same=allocated(a%fields).eqv.allocated(b%fields)
    if(.not.same.or..not.allocated(a%fields)) return
    same=size(a%fields)==size(b%fields)
    if(.not.same) return
    do node=1,size(a%fields)
      if(owners(node)/=rank) then
        if(allocated(a%fields(node)%values)) same=.false.
        cycle
      endif
      if(allocated(a%fields(node)%values).neqv.allocated(b%fields(node)%values)) then
        same=.false.;return
      endif
      if(.not.allocated(a%fields(node)%values)) cycle
      same=same.and.all(a%fields(node)%values==b%fields(node)%values)
    enddo
  end function record_field_equal
  logical function record_equal(a,b) result(same)
    type(eb_godunov_interval_3d), intent(in) :: a,b
    integer :: node,stage
    same=a%depth==b%depth.and.a%iterations==b%iterations.and.a%interval==b%interval.and.a%face_weight==b%face_weight
    same=same.and.(a%routed_auxiliary.eqv.b%routed_auxiliary).and. &
      record_field_equal(a%routed_auxiliary_impulse,b%routed_auxiliary_impulse)
    same=same.and.record_field_equal(a%clock_gas,b%clock_gas).and. &
      record_field_equal(a%clock_temperature,b%clock_temperature).and. &
      record_field_equal(a%clock_reaction_source,b%clock_reaction_source).and. &
      record_field_equal(a%old_predictor,b%old_predictor).and. &
      record_field_equal(a%new_predictor,b%new_predictor).and. &
      record_field_equal(a%clock_internal_energy,b%clock_internal_energy).and. &
      record_field_equal(a%initial_gas,b%initial_gas).and.record_field_equal(a%initial_temperature,b%initial_temperature).and. &
      record_field_equal(a%initial_internal_energy,b%initial_internal_energy).and. &
      record_field_equal(a%initial_reaction_source,b%initial_reaction_source).and. &
      record_field_equal(a%final_gas,b%final_gas).and.record_field_equal(a%final_temperature,b%final_temperature).and. &
      record_field_equal(a%final_internal_energy,b%final_internal_energy).and. &
      record_field_equal(a%reaction_source,b%reaction_source)
    same=same.and.abs(a%pressure_work_change-b%pressure_work_change)<= &
      8*epsilon(1._dp)*max(tiny(1._dp),abs(b%pressure_work_change)).and. &
      abs(a%numerical_energy_change-b%numerical_energy_change)<= &
      8*epsilon(1._dp)*max(tiny(1._dp),abs(b%numerical_energy_change))
    do stage=0,a%iterations
      same=same.and.(a%sources(stage)%auxiliary_accepted.eqv.b%sources(stage)%auxiliary_accepted).and. &
        (a%sources(stage)%provider_accepted.eqv.b%sources(stage)%provider_accepted)
      same=same.and.a%sources(stage)%index==stage.and.b%sources(stage)%index==stage.and. &
        record_field_equal(a%sources(stage)%gas,b%sources(stage)%gas).and. &
        record_field_equal(a%sources(stage)%temperature,b%sources(stage)%temperature).and. &
        record_field_equal(a%sources(stage)%internal_energy,b%sources(stage)%internal_energy).and. &
        record_field_equal(a%sources(stage)%external,b%sources(stage)%external).and. &
        record_field_equal(a%sources(stage)%auxiliary,b%sources(stage)%auxiliary).and. &
        record_field_equal(a%sources(stage)%rate,b%sources(stage)%rate).and. &
        record_field_equal(a%sources(stage)%scalar_rate,b%sources(stage)%scalar_rate)
      do node=1,np
        if(owners(node)/=rank.or.(a%depth>=0.and.topology%nodes(node)%level/=a%depth)) cycle
        if(allocated(a%sources(stage)%faces)) same=same.and. &
          equal_flux(a%sources(stage)%faces(node),b%sources(stage)%faces(node))
        same=same.and.equal_flux(a%sources(stage)%scalar_faces(node),b%sources(stage)%scalar_faces(node))
        if(allocated(a%sources(stage)%registers)) same=same.and. &
          equal_register(a%sources(stage)%registers(node),b%sources(stage)%registers(node)).and. &
          equal_register(a%sources(stage)%scalar_registers(node),b%sources(stage)%scalar_registers(node))
      enddo
    enddo
    do stage=1,a%iterations
      same=same.and.a%hydro(stage)%index==stage.and.b%hydro(stage)%index==stage.and. &
        record_field_equal(a%hydro(stage)%predictor,b%hydro(stage)%predictor).and. &
        record_field_equal(a%hydro(stage)%rate,b%hydro(stage)%rate).and. &
        record_field_equal(a%hydro(stage)%scalar_rate,b%hydro(stage)%scalar_rate).and. &
        record_field_equal(a%hydro(stage)%pressure_work,b%hydro(stage)%pressure_work)
      do node=1,np
        if(owners(node)/=rank.or.(a%depth>=0.and.topology%nodes(node)%level/=a%depth)) cycle
        same=same.and.equal_flux(a%hydro(stage)%faces(node),b%hydro(stage)%faces(node)).and. &
          equal_flux(a%hydro(stage)%scalar_faces(node),b%hydro(stage)%scalar_faces(node))
        if(allocated(a%hydro(stage)%registers)) same=same.and. &
          equal_register(a%hydro(stage)%registers(node),b%hydro(stage)%registers(node)).and. &
          equal_register(a%hydro(stage)%scalar_registers(node),b%hydro(stage)%scalar_registers(node))
      enddo
    enddo
    same=same.and.record_field_equal(a%closed_gas,b%closed_gas).and. &
      record_field_equal(a%closed_temperature,b%closed_temperature).and. &
      record_field_equal(a%closed_internal_energy,b%closed_internal_energy).and. &
      record_field_equal(a%closed_reaction_source,b%closed_reaction_source).and. &
      record_field_equal(a%clock_heat,b%clock_heat)
    if(allocated(a%clock_start)) same=same.and.all(a%clock_start==b%clock_start).and. &
      all(a%clock_interval==b%clock_interval)
    do node=1,np
      if(owners(node)/=rank) cycle
      if(allocated(a%closed_faces)) same=same.and.equal_flux(a%closed_faces(node),b%closed_faces(node)).and. &
        equal_flux(a%closed_scalar_faces(node),b%closed_scalar_faces(node))
      if(allocated(a%open_ancestor_registers)) same=same.and. &
        equal_register(a%open_ancestor_registers(node),b%open_ancestor_registers(node)).and. &
        equal_register(a%open_ancestor_scalar_registers(node),b%open_ancestor_scalar_registers(node))
    enddo
    do node=1,np
      if(owners(node)/=rank.or.(a%depth>=0.and.topology%nodes(node)%level/=a%depth)) cycle
      same=same.and.equal_flux(a%faces(node),b%faces(node)).and.equal_flux(a%scalar_faces(node),b%scalar_faces(node)).and. &
        equal_flux(a%registered_faces(node),b%registered_faces(node))
      if(allocated(a%registers)) same=same.and.equal_register(a%registers(node),b%registers(node)).and. &
        equal_register(a%scalar_registers(node),b%scalar_registers(node))
    enddo
  end function record_equal
  logical function equal_state(a,b) result(same)
    type(amr_tree_state_3d), intent(in) :: a,b
    integer :: node
    same=allocated(a%fields).eqv.allocated(b%fields)
    if(.not.same.or..not.allocated(a%fields)) return
    same=size(a%fields)==size(b%fields)
    if(.not.same) return
    do node=1,size(a%fields)
      same=allocated(a%fields(node)%values).eqv.allocated(b%fields(node)%values)
      if(.not.same) return
      if(.not.allocated(a%fields(node)%values)) cycle
      same=all(lbound(a%fields(node)%values)==lbound(b%fields(node)%values)).and. &
        all(ubound(a%fields(node)%values)==ubound(b%fields(node)%values))
      if(.not.same) return
      same=all(a%fields(node)%values==b%fields(node)%values)
      if(.not.same) return
    enddo
  end function equal_state
  logical function equal_flux(a,b) result(same)
    type(eb_integrated_flux_3d), intent(in) :: a,b
    same=(allocated(a%x).eqv.allocated(b%x)).and.(allocated(a%y).eqv.allocated(b%y)).and. &
      (allocated(a%z).eqv.allocated(b%z))
    if(.not.same) return
    if(allocated(a%x)) same=same.and.all(a%x==b%x)
    if(allocated(a%y)) same=same.and.all(a%y==b%y)
    if(allocated(a%z)) same=same.and.all(a%z==b%z)
  end function equal_flux
  logical function equal_register(a,b) result(same)
    type(state_redist_register_3d), intent(in) :: a,b
    same=(allocated(a%coarse_delta).eqv.allocated(b%coarse_delta)).and. &
      (allocated(a%fine_delta).eqv.allocated(b%fine_delta))
    if(.not.same) return
    if(allocated(a%coarse_delta)) same=same.and.all(a%coarse_delta==b%coarse_delta)
    if(allocated(a%fine_delta)) same=same.and.all(a%fine_delta==b%fine_delta)
  end function equal_register
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(condition) return
    print *, 'FAIL rank ',rank,': ',label
    call MPI_Abort(MPI_COMM_WORLD,1,status)
    error stop 1
  end subroutine require
end program test_mpi_eb_godunov_held_3d
