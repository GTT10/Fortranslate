module mpi_eb_tree_full_3d_mod
  use mpi_f08
  use eb_reactive_transfer_3d_mod, only: reconstruct_eb_pressure_state_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_next_after
  use eb_held_funding_queue_3d_mod, only: eb_held_funding_context_3d, &
    valid_eb_held_funding_context_3d,move_eb_held_funding_context_3d, &
    initialize_eb_held_funding_context_3d,copy_eb_held_funding_context_3d, &
    set_eb_held_funding_clock_faces_3d,close_eb_held_funding_clock_faces_3d, &
    read_eb_held_funding_open_faces_3d,take_eb_held_funding_target_3d,take_eb_held_funding_bucket_3d, &
    read_eb_held_funding_stock_3d,queue_eb_held_funding_receipt_3d
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_momentum_source_3d,rotating_potential_3d
  use constant_volume_reactor_mod, only: forced_reactor_integrator,advance_forced_reactor_interval, &
    forced_reactor_admission,reaction_admission_not_attempted,reaction_admission_accepted, &
    reaction_admission_rejected,reaction_admission_unknown,reaction_admission_provider_refused
  use dual_energy_mod, only: dual_energy_policy,reset_nasa7_dual_energy_cell
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use reactive_3d_mod, only: advance_reactive_chemistry_3d,react_sdc_state_3d
  use reactive_eb_3d_driver_mod, only: react_sdc_eb_state_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,extract_eb_geometry_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d,fetch_mpi_amr_tree_cells_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d,average_down_mpi_eb_tree_3d, &
    validate_mpi_eb_tree_geometry_3d,average_down_mpi_eb_tree_field_3d
  use mpi_eb_tree_hydro_3d_mod, only: advance_mpi_eb_tree_hydro_3d,compute_sdc_mpi_eb_tree_hydro_3d, &
    compute_mol_mpi_eb_tree_hydro_3d,apply_mpi_eb_flux_correction_3d, &
    compute_godunov_state_redist_mpi_eb_tree_hydro_3d
  use mpi_eb_tree_transport_3d_mod, only: advance_mpi_eb_tree_transport_3d,compute_sdc_mpi_eb_tree_transport_3d, &
    compute_sdc_mpi_eb_tree_transport_halos_3d
  use mpi_eb_tree_subcycle_hydro_3d_mod, only: advance_mpi_eb_tree_subcycle_hydro_3d
  use mpi_eb_tree_subcycle_transport_3d_mod, only: advance_mpi_eb_tree_subcycle_transport_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_auxiliary_nvar,pack_eb_held_auxiliary_state_3d, &
    apply_eb_held_auxiliary_increment_3d,eb_godunov_interval_3d,move_godunov_interval_3d, &
    eb_godunov_source_stage_3d,eb_native_mol_source_prefix_3d, &
    move_godunov_source_stage_3d,extend_godunov_source_stage_3d, &
    copy_godunov_field_3d,copy_godunov_fluxes_3d,copy_godunov_registers_3d,valid_godunov_source_cache_3d, &
    valid_godunov_clock_anchor_3d,same_godunov_predictor_3d,rebuild_godunov_held_clock_prefixes_3d, &
    blend_native_mol_face_pair_3d,valid_native_mol_clock_profile_3d,valid_native_mol_source_prefix_3d
  use eb_reactive_two_level_3d_mod, only: matching_eb_metrics_3d
  use eb_reactive_reflux_3d_mod, only: compute_reactive_eb_reactor_bias_3d,admit_reactive_eb_reactor_energy_3d, &
    compute_reactive_eb_held_projection,compute_accepted_eb_held_clock_prefixes_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,valid_eb_integrated_flux_3d
  use eb_reactive_wall_flux_3d_mod, only: reactive_eb_slip_wall_source_3d
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use mpi_eb_tree_redistribution_3d_mod, only: redistribute_sdc_mpi_eb_tree_flux_rate_3d, &
    redistribute_sdc_mpi_eb_tree_level_rate_3d,redistribute_mpi_eb_tree_state_3d, &
    redistribute_sdc_mpi_eb_tree_state_rate_3d,redistribute_sdc_mpi_eb_tree_state_level_rate_3d, &
    compute_state_redist_mpi_level_stage_3d,assemble_state_redist_mpi_interface_rereflux_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d
  use mpi_eb_tree_flux_3d_mod, only: synchronize_mpi_eb_tree_flux_3d,assemble_mpi_eb_tree_reflux_3d, &
    assemble_mpi_eb_held_funding_debt_3d, &
    repair_mpi_eb_tree_interface_budget_3d,repair_mpi_eb_tree_held_interface_budget_3d, &
    synchronize_mpi_eb_tree_register_flux_3d,synchronize_mpi_eb_tree_scalar_flux_3d, &
    synchronize_mpi_eb_tree_scalar_register_flux_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use eb_reactive_redistribution_3d_mod, only: state_redist_register_3d,accumulate_state_redist_register_3d, &
    initialize_state_redist_tree_register_3d,accumulate_state_redist_coarse_register_3d, &
    accumulate_state_redist_mol_registration_stage_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: agree_mpi_constant_transport_coefficients
  use eb_reactive_redistribution_3d_mod, only: post_eb_spatial_density_rate_3d,initialize_flux_rereflux_register_3d
  implicit none
  private
  public :: initialize_mpi_eb_tree_reaction_source_3d
  public :: react_mpi_eb_tree_3d,advance_mpi_eb_tree_full_3d,react_sdc_mpi_eb_tree_3d
  public :: advance_mpi_eb_tree_sdc_3d,advance_mpi_eb_tree_sdc_subcycle_3d
  public :: reconstruct_mpi_eb_mol_interval_3d,replay_mpi_eb_mol_held_interval_3d
  public :: reconstruct_mpi_eb_mol_state_redist_interval_3d
  ! Auxiliary sources receive the old state at stage 0 and provisional new
  ! states at stages 1..iterations. Rates are per fluid volume per time.
  ! Callback-owned physical state must remain private until the driver succeeds.
  ! Opt-in auxiliary_has_internal_energy uses Nv+1 source input/rate/budget
  ! rows, with actual held UEINT in the final row; gas operators use rows1:Nv.
  public :: mpi_eb_sdc_diffusion_source,mpi_eb_sdc_spatial_diffusion_source
  abstract interface
    ! Rates per fluid volume/time; flux integrated over interval per open face area.
    subroutine mpi_eb_sdc_diffusion_source(state,temperature,interval,depth,rate,flux,theta,comm,ok, &
        ancestor_start,ancestor_end,temperature_start,temperature_end,fractions)
      import :: amr_tree_state_3d,eb_integrated_flux_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: state,temperature
      real(dp), intent(in) :: interval
      integer, intent(in) :: depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: theta
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
      type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
      real(dp), optional, intent(in) :: fractions(:)
    end subroutine

    ! Additional diffusion before StateRedist. The supplied old_halos are
    ! the actual six-layer published-energy state used by built-in operators.
    ! Return owner-local one-based valid rates and one-based spatial arrays
    ! of shape valid+12 (grow3 evaluated); valid crops match rate exactly.
    ! Same-level/periodic halo rates are refreshed from accepted valid owners.
    subroutine mpi_eb_sdc_spatial_diffusion_source(state,temperature,geometry,old_halos,interval,depth, &
        rate,spatial_rate,flux,theta,comm,ok,frame)
      import :: amr_tree_state_3d,eb_geometry_3d,eb_integrated_flux_3d,dp,MPI_Comm,rotating_frame_3d
      type(amr_tree_state_3d), intent(in) :: state,temperature,old_halos
      type(eb_geometry_3d), intent(in) :: geometry(:)
      real(dp), intent(in) :: interval
      integer, intent(in) :: depth
      type(amr_tree_state_3d), intent(inout) :: rate,spatial_rate
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: theta
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
      type(rotating_frame_3d), optional, intent(in) :: frame
    end subroutine mpi_eb_sdc_spatial_diffusion_source

    ! Called after descendants and gas synchronization, before publication.
    ! Callback state must remain private until the entire outer step succeeds.
    subroutine mpi_eb_sdc_post_level(depth,time,interval,iteration,ncycle,comm,ok)
      import :: dp,MPI_Comm
      integer, intent(in) :: depth,iteration,ncycle
      real(dp), intent(in) :: time,interval
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    ! Build the complete predictor halo for iteration 1..iterations.
    ! Required width is 3 for valid-only CTU, 6 for spatial Godunov,
    ! and 7 for corrected spatial Godunov; wider cached halos are cropped.
    ! Valid fluid cells must equal rate; callback state remains caller-owned.
    subroutine mpi_eb_sdc_predictor_builder(rate,interval,iteration,depth,halo,comm,ok)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: rate
      real(dp), intent(in) :: interval
      integer, intent(in) :: iteration,depth
      type(amr_tree_state_3d), intent(inout) :: halo
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    ! Recursive predictor donors use the midpoint of the level interval.
    subroutine mpi_eb_sdc_recursive_predictor_builder(rate,interval,iteration,depth,time,halo,comm,ok, &
        source_start,source_end,source_fraction)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: rate,source_start,source_end
      real(dp), intent(in) :: interval,time,source_fraction(:)
      integer, intent(in) :: iteration,depth
      type(amr_tree_state_3d), intent(inout) :: halo
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    subroutine mpi_eb_sdc_auxiliary_source(state,temperature,interval,stage,depth,rate,comm,ok)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: state,temperature
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    ! Time is the level interval start relative to the outer step, for all stages.
    subroutine mpi_eb_sdc_recursive_source(state,temperature,interval,stage,depth,time,rate,comm,ok)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: state,temperature
      real(dp), intent(in) :: interval,time
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    subroutine mpi_eb_sdc_recursive_halo_source(state,temperature,interval,stage,depth,time,rate,comm,ok, &
        ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: state,temperature
      type(amr_tree_state_3d), intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
      real(dp), intent(in) :: interval,time,ancestor_fraction(:)
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
  end interface
  public :: mpi_eb_sdc_post_level
  public :: mpi_eb_sdc_predictor_builder,mpi_eb_sdc_recursive_predictor_builder
  public :: mpi_eb_sdc_auxiliary_source,mpi_eb_sdc_recursive_source,mpi_eb_sdc_recursive_halo_source
  ! One accepted local MOL solve, before its coarse/fine synchronization.
  ! Fields are sparse: only nodes owned by this rank at depth are allocated.
  ! The complete forcing includes wall/volume terms absent from face fluxes.
  ! synchronized_flux is filled by the recursive driver after level closure;
  ! ancestor closure can subsequently change descendant faces and states.
  ! replayed records hold the replacement solve; old_flux keeps its captured
  ! original spatial anchor, while initial/forcing/endpoint describe the replay.
  ! Re-evaluating StateRedist needs each actual old/new source stage. Its
  ! state/slopes and register geometry cannot be reconstructed from a final
  ! mean RHS or a frozen accepted register delta.
  type, public :: mpi_eb_mol_state_redist_stage_3d
    integer :: source_stage=-1,depth=-1,max_order=0
    real(dp) :: interval=0,target_volume_fraction=.5_dp,register_scale=.5_dp,spatial_weight=1
    type(amr_tree_state_3d) :: state_halos,raw_rate_halos,spatial_rate
    type(amr_tree_state_3d) :: external_rate,auxiliary_rate
    type(amr_tree_state_3d) :: low_spatial_rate
    type(eb_integrated_flux_3d), allocatable :: faces(:),native_faces(:),low_faces(:)
    type(state_redist_register_3d), allocatable :: registers(:)
    logical :: has_internal_energy=.false.
    type(amr_tree_state_3d) :: held_state_halos,held_raw_rate_halos,held_spatial_rate,held_low_spatial_rate
    type(eb_integrated_flux_3d), allocatable :: held_faces(:),held_native_faces(:),held_low_faces(:)
    type(state_redist_register_3d), allocatable :: held_registers(:)
  end type mpi_eb_mol_state_redist_stage_3d

  type, public :: mpi_eb_mol_interval_3d
    logical :: separate_momentum=.false.,consistent_velocity=.false.,paired_pressure=.false.,connected_directions=.false.
    integer :: depth=-1
    ! Owner-local declared auxiliary stages, separate from recomputed native
    ! force work by default. Parent trials can freeze external_stages explicitly;
    ! every replacement loan/face family is rebuilt. Prepared caches retain
    ! their existing nonlinear operator/source contract.
    integer :: native_iterations=0
    ! Accepted ordinary native gas/UEINT spatial blend for this clock.
    real(dp) :: native_spatial_weight=1
    real(dp) :: native_force(3)=0,native_acceleration(3)=0
    logical :: has_native_frame=.false.
    type(rotating_frame_3d) :: native_frame
    type(amr_tree_state_3d), allocatable :: auxiliary_stages(:)
    ! Keep force/acceleration/body-work/frame packets with native source stages;
    ! parent source-frozen spatial trials use these accepted packets explicitly.
    type(amr_tree_state_3d), allocatable :: external_stages(:)
    logical :: native_providers(2)=.false.
    integer :: native_provider_halo_width=6
    type(eb_godunov_source_stage_3d), allocatable :: provider_stages(:,:)
    type(eb_geometry_3d), allocatable :: native_provider_geometry(:)
    ! True when native replay reran this local solve, including parent repair.
    logical :: replayed=.false.
    character(len=16) :: redistribution='none'
    logical :: prepared_redistribution=.false.
    ! Held reacting MOL commits only the old/final mean source. Intermediate
    ! evaluations remain in prepared_stages, but carry no accepted CF impulse.
    ! This differs from the pinned direct provisional coarse accumulation;
    ! retain an explicit contract so replay never guesses the ledger family.
    logical :: accepted_coarse_stages=.false.
    logical :: spatial_diffusion_is_source=.false.,auxiliary_has_internal_energy=.false.,routed_auxiliary=.false.
    real(dp) :: time=0,interval=0
    type(amr_tree_state_3d) :: initial,initial_temperature,forcing,mass_source
    ! Held UEINT and its spatial family are distinct from the projected gas
    ! reaction forcing. Only an accepted common interval publishes all fields.
    logical :: has_internal_energy=.false.
    type(dual_energy_policy) :: internal_energy_policy
    real(dp) :: initial_energy_correction=0,numerical_energy_correction=0
    real(dp) :: internal_energy_projection_change=0,internal_energy_closure_change=0
    real(dp) :: initial_internal_energy_change=0
    ! Accepted off-clock borrowing remains a separate gas/actual-UEINT
    ! impulse. Replay must retain its frozen stock and inactive face family.
    logical :: has_source_funding=.false.
    type(amr_tree_state_3d) :: source_funding_stock,source_funding_impulse
    ! The wrapper clock anchor precedes target receipt payment. Restore its
    ! whole pending queue and every open ancestor family with gas/H/T.
    type(eb_held_funding_context_3d) :: clock_source_funding_context
    ! Wrapper EOS impulses occur at clock entry and after all descendants.
    ! Prefixes become valid only for a completely observed accepted forest.
    real(dp) :: held_clock_start_change=0,held_clock_end_change=0
    real(dp) :: held_clock_projection_before=0,held_clock_reset_before=0
    real(dp) :: held_clock_projection_after=0,held_clock_reset_after=0
    logical :: held_clock_events_closed=.false.,held_clock_prefixes_valid=.false.
    ! A predictor repair retained in the frozen gas forcing is still charged
    ! on replay. Endpoint/closure resets are recomputed from the new solve.
    real(dp) :: forcing_energy_correction=0
    type(amr_tree_state_3d) :: clock_internal_energy
    real(dp) :: clock_energy_correction=0
    type(amr_tree_state_3d) :: initial_internal_energy,local_internal_energy
    type(amr_tree_state_3d) :: old_internal_energy_spatial_rate,internal_energy_spatial_rate
    type(eb_integrated_flux_3d), allocatable :: old_internal_energy_flux(:),internal_energy_flux(:)
    type(amr_tree_state_3d) :: local_state,local_temperature,reaction_source
    type(amr_tree_state_3d) :: initial_reaction_source
    type(amr_tree_state_3d) :: old_predictor,new_predictor
    ! Recursive replay restores the complete clock start, including ancestors
    ! modified by a completed descendant rereflux, rather than only this level.
    type(amr_tree_state_3d) :: clock_initial,clock_temperature,clock_reaction_source,clock_heat
    ! Native sub-parent replay retains the still-open ancestor trajectories;
    ! closed descendant averages cannot reconstruct their original halo inputs.
    type(amr_tree_state_3d) :: ancestor_start,ancestor_end,ancestor_temperature_start,ancestor_temperature_end
    type(amr_tree_state_3d) :: ancestor_internal_energy_start,ancestor_internal_energy_end
    type(amr_tree_state_3d) :: ancestor_source_start,ancestor_source_end
    real(dp), allocatable :: fractions_start(:),fractions_end(:)
    real(dp), allocatable :: clock_start(:),clock_interval(:)
    type(state_redist_register_3d), allocatable :: open_ancestor_registers(:),open_ancestor_held_registers(:)
    logical :: has_auxiliary_source=.false.
    ! Frozen off-level source impulse, applied at min(source,target) closure.
    ! Its covered/solid/active-source-level entries are zero.
    type(amr_tree_state_3d) :: routed_source
    ! Redistribution replay needs the operator family, not only its faces.
    ! These sparse rates are captured without subtracting near-equal fluxes.
    type(amr_tree_state_3d) :: old_volume_source,volume_source,old_spatial_rate,spatial_rate
    type(amr_tree_state_3d) :: non_face_spatial_rate
    ! Unprepared transfer impulses span every owned target level.
    type(amr_tree_state_3d) :: old_redistribution_transfer,redistribution_transfer
    ! Accepted closure spill to ancestor buckets, replaced on parent replay.
    type(amr_tree_state_3d) :: synchronization_transfer
    ! Prepared StateRedist uses its original masks, scales and halo deltas.
    type(state_redist_register_3d), allocatable :: old_redist_registers(:),redist_registers(:)
    type(state_redist_register_3d), allocatable :: old_internal_energy_registers(:),internal_energy_registers(:)
    ! Bounds 0:mol_iters retain every source stage, before coarse-register
    ! accumulation across reacting corrections. These are published only on
    ! acceptance, together with the final weighted register family above.
    type(mpi_eb_mol_state_redist_stage_3d), allocatable :: prepared_stages(:)
    type(eb_integrated_flux_3d), allocatable :: old_flux(:),spatial_flux(:),synchronized_flux(:)
    real(dp), allocatable :: chemical_change(:),wall_change(:),source_change(:),auxiliary_change(:)
  end type mpi_eb_mol_interval_3d
contains
  subroutine rebuild_mpi_mol_held_clock_prefixes_3d(records,ok)
    type(mpi_eb_mol_interval_3d), intent(inout) :: records(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: pb(:),rb(:),pa(:),ra(:),local_projection(:)
    integer :: r
    ok=.false.
    if(size(records)<1) return
    if(.not.all(records%held_clock_events_closed)) then
      ! A replaced group must first measure its new entry/exit events.
      ! Its physical result remains valid; stale diagnostic prefixes do not.
      records%held_clock_prefixes_valid=.false.
      records%held_clock_projection_before=0;records%held_clock_reset_before=0
      records%held_clock_projection_after=0;records%held_clock_reset_after=0
      ok=.true.;return
    endif
    ! Native local projection includes the entry EOS repair. Hierarchy
    ! prefixes put that repair in the start event, exactly once, so a
    ! canonical replay can retain its historical start without a stale solve.
    allocate(local_projection(size(records)))
    do r=1,size(records)
      call compute_reactive_eb_held_projection(records(r)%initial_internal_energy_change,0._dp,0._dp, &
        records(r)%internal_energy_projection_change,local_projection(r),ok)
      if(.not.ok) return
    enddo
    call compute_accepted_eb_held_clock_prefixes_3d(records%depth,records%held_clock_start_change, &
      local_projection,records%internal_energy_closure_change, &
      records%held_clock_end_change,0._dp,0._dp,pb,rb,pa,ra,ok)
    if(.not.ok) return
    do r=1,size(records)
      records(r)%held_clock_projection_before=pb(r);records(r)%held_clock_reset_before=rb(r)
      records(r)%held_clock_projection_after=pa(r);records(r)%held_clock_reset_after=ra(r)
      records(r)%held_clock_prefixes_valid=.true.
    enddo
  end subroutine rebuild_mpi_mol_held_clock_prefixes_3d

  subroutine agree_held_clock_packet(values,comm,ok)
    real(dp), intent(in) :: values(9)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    real(dp) :: lower(9),upper(9)
    integer :: status
    call agree(all(ieee_is_finite(values)),comm,ok)
    if(.not.ok) return
    call MPI_Allreduce(values,lower,9,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) then
      ok=.false.;return
    endif
    call MPI_Allreduce(values,upper,9,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    ok=status==MPI_SUCCESS.and.all(lower==upper)
  end subroutine agree_held_clock_packet

  subroutine agree_held_projection(value,comm,ok)
    real(dp), intent(in) :: value
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    real(dp) :: lower,upper
    integer :: status
    call agree(ieee_is_finite(value),comm,ok)
    if(.not.ok) return
    call MPI_Allreduce(value,lower,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) then
      ok=.false.;return
    endif
    call MPI_Allreduce(value,upper,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    ok=status==MPI_SUCCESS
    if(ok) ok=lower==upper
  end subroutine agree_held_projection

  ! Replace only the face-divergence part of an accepted interval's forcing.
  ! Keep its actual MOL wall, volume, auxiliary and redistribution contributions.
  ! The returned full owner-local tree can be passed to the reaction solver;
  ! inactive levels are zero. This does not synchronize faces or admit a state.
  subroutine reconstruct_mpi_eb_mol_interval_3d(tree,owners,geometry,history,flux,nvar,rates,comm,ok, &
      internal_energy_flux,internal_energy_rate)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(mpi_eb_mol_interval_3d), intent(in) :: history
    type(eb_integrated_flux_3d), intent(in) :: flux(:)
    type(amr_tree_state_3d), intent(inout) :: rates
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    ! These are paired spatial sources, never the reactor energy forcing.
    type(eb_integrated_flux_3d), optional, intent(in) :: internal_energy_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate
    type(amr_tree_state_3d) :: candidate,held_candidate
    integer :: held_options,policy_flags(3),policy_min(3),policy_max(3)
    real(dp) :: held_controls(15),held_min(15),held_max(15),held_change
    integer :: rank,status,node,n(3),ios,i,j,k,controls(5),low(5),high(5)
    real(dp) :: duration_min,duration_max,volume,change(nvar)
    logical :: valid,ready
    ok=.false.
    held_options=merge(1,0,history%has_internal_energy)+merge(2,0,present(internal_energy_flux))+ &
      merge(4,0,present(internal_energy_rate))
    controls=[nvar,history%depth,held_options,merge(1,0,history%has_source_funding), &
      merge(1,0,allocated(history%clock_source_funding_context%pending))]
    call MPI_Allreduce(controls,low,5,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,5,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    if(held_options/=0.and.held_options/=7) return
    valid=nvar>=1.and.history%depth>=0.and.history%depth<tree%level_count()
    valid=valid.and.ieee_is_finite(history%interval).and.history%interval>0
    valid=valid.and.size(flux)==size(tree%nodes).and.allocated(history%forcing%fields).and. &
      allocated(history%spatial_flux)
    if(allocated(history%forcing%fields)) valid=valid.and.size(history%forcing%fields)==size(tree%nodes)
    if(allocated(history%spatial_flux)) valid=valid.and.size(history%spatial_flux)==size(tree%nodes)
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(history%interval,duration_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(history%interval,duration_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.duration_min/=duration_max) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nvar,comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    valid=allocated(history%clock_source_funding_context%pending).eqv. &
      allocated(history%clock_source_funding_context%open)
    if(allocated(history%clock_source_funding_context%pending)) then
      valid=valid.and.history%has_source_funding
      valid=valid.and.valid_eb_held_funding_context_3d(tree,geometry,nvar,owners==rank, &
        history%clock_source_funding_context,history%depth)
    else
      valid=valid.and..not.allocated(history%clock_source_funding_context%topology%nodes).and. &
        .not.allocated(history%clock_source_funding_context%geometry)
    endif
    call agree(valid,comm,ready);if(.not.ready) return
    if(history%has_internal_energy) then
      ! A gas-only archive cannot describe the held operator or its reset
      ! ledger. Policy and global corrections must agree even on idle ranks.
      held_controls=[history%internal_energy_policy%eta2,history%initial_energy_correction, &
        history%numerical_energy_correction,history%forcing_energy_correction,history%internal_energy_projection_change, &
        history%internal_energy_closure_change, &
        history%initial_internal_energy_change, &
        history%held_clock_start_change,history%held_clock_end_change, &
        history%held_clock_projection_before,history%held_clock_reset_before, &
        history%held_clock_projection_after,history%held_clock_reset_after, &
        merge(1._dp,0._dp,history%held_clock_events_closed), &
        merge(1._dp,0._dp,history%held_clock_prefixes_valid)]
      valid=all(ieee_is_finite(held_controls)).and.held_controls(1)>=0
      valid=valid.and.(.not.history%held_clock_prefixes_valid.or.history%held_clock_events_closed)
      valid=valid.and.history%redistribution=='none'.and..not.history%prepared_redistribution
      valid=valid.and.size(internal_energy_flux)==size(tree%nodes).and. &
        allocated(history%internal_energy_flux).and.allocated(history%old_internal_energy_flux)
      if(allocated(history%internal_energy_flux)) valid=valid.and.size(history%internal_energy_flux)==size(tree%nodes)
      if(allocated(history%old_internal_energy_flux)) &
        valid=valid.and.size(history%old_internal_energy_flux)==size(tree%nodes)
      valid=valid.and.valid_held_interval_field(history%initial_internal_energy,tree,owners,history%depth,rank).and. &
        valid_held_interval_field(history%local_internal_energy,tree,owners,history%depth,rank).and. &
        valid_held_interval_field(history%old_internal_energy_spatial_rate,tree,owners,history%depth,rank).and. &
        valid_held_interval_field(history%internal_energy_spatial_rate,tree,owners,history%depth,rank)
      call agree(valid,comm,ready);if(.not.ready) return
      call MPI_Allreduce(held_controls,held_min,15,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(held_controls,held_max,15,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(held_min/=held_max)) return
      policy_flags=[merge(1,0,history%internal_energy_policy%allow_small_energy), &
        merge(1,0,history%internal_energy_policy%allow_negative_energy), &
        merge(1,0,history%internal_energy_policy%update_total_from_internal)]
      call MPI_Allreduce(policy_flags,policy_min,3,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(policy_flags,policy_max,3,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(policy_min/=policy_max)) return
    endif
    if(history%has_source_funding) then
      valid=history%has_internal_energy.and..not.history%prepared_redistribution
      valid=valid.and.valid_held_interval_field(history%source_funding_stock,tree,owners,-1,rank,nvar+1).and. &
        valid_held_interval_field(history%source_funding_impulse,tree,owners,-1,rank,nvar+1).and. &
        valid_held_interval_field(history%clock_initial,tree,owners,-1,rank,nvar).and. &
        valid_held_interval_field(history%clock_internal_energy,tree,owners,-1,rank)
      call agree(valid,comm,ready);if(.not.ready) return
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(tree%nodes(node)%level==history%depth) then
          valid=valid.and.all(history%source_funding_stock%fields(node)%values==0).and. &
            all(history%source_funding_impulse%fields(node)%values==0)
        endif
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(.not.covered(tree,node,[i,j,k]).and.geometry(node)%volume_fraction(i,j,k)>0) cycle
          valid=valid.and.all(history%source_funding_stock%fields(node)%values(:,i,j,k)==0).and. &
            all(history%source_funding_impulse%fields(node)%values(:,i,j,k)==0)
        enddo
        enddo
        enddo
      enddo
      call agree(valid,comm,ready);if(.not.ready) return
    else
      call agree(.not.allocated(history%source_funding_stock%fields).and. &
        .not.allocated(history%source_funding_impulse%fields),comm,ready)
      if(.not.ready) return
    endif
    valid=.true.
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      if(tree%nodes(node)%level/=history%depth) then
        if(.not.history%has_source_funding) cycle
        valid=valid_eb_integrated_flux_3d(flux(node),geometry(node),nvar).and. &
          valid_eb_integrated_flux_3d(history%spatial_flux(node),geometry(node),nvar).and. &
          valid_eb_integrated_flux_3d(internal_energy_flux(node),geometry(node),1).and. &
          valid_eb_integrated_flux_3d(history%internal_energy_flux(node),geometry(node),1)
        if(.not.valid) exit
        cycle
      endif
      if(.not.allocated(history%forcing%fields(node)%values)) then
        valid=.false.;exit
      endif
      n=tree%nodes(node)%shape
      valid=all(shape(history%forcing%fields(node)%values)==[nvar,n]).and. &
        all(lbound(history%forcing%fields(node)%values)==1)
      valid=valid.and.all(ieee_is_finite(history%forcing%fields(node)%values))
      valid=valid.and.valid_eb_integrated_flux_3d(flux(node),geometry(node),nvar).and. &
        valid_eb_integrated_flux_3d(history%spatial_flux(node),geometry(node),nvar)
      if(history%has_internal_energy) then
        valid=valid.and.valid_eb_integrated_flux_3d(internal_energy_flux(node),geometry(node),1).and. &
          valid_eb_integrated_flux_3d(history%internal_energy_flux(node),geometry(node),1).and. &
          valid_eb_integrated_flux_3d(history%old_internal_energy_flux(node),geometry(node),1)
      endif
      if(.not.valid) exit
    enddo
    call agree(valid,comm,ready)
    if(.not.ready) return
    allocate(candidate%fields(size(tree%nodes)),stat=ios)
    call agree(ios==0,comm,ready)
    if(.not.ready) return
    if(history%has_internal_energy) then
      allocate(held_candidate%fields(size(tree%nodes)),stat=ios)
      call agree(ios==0,comm,ready);if(.not.ready) return
    endif
    valid=.true.
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(candidate%fields(node)%values(nvar,n(1),n(2),n(3)),source=0._dp,stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      endif
      if(history%has_internal_energy) then
        allocate(held_candidate%fields(node)%values(1,n(1),n(2),n(3)),source=0._dp,stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        endif
      endif
      if(tree%nodes(node)%level/=history%depth) cycle
      candidate%fields(node)%values=history%forcing%fields(node)%values
      if(history%has_internal_energy) &
        held_candidate%fields(node)%values=history%internal_energy_spatial_rate%fields(node)%values
      associate(g=>geometry(node),f=>flux(node),old=>history%spatial_flux(node))
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(g%volume_fraction(i,j,k)<=0) cycle
        volume=g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz
        change=g%dy*g%dz*(g%x_face_fraction(i-1,j,k)*(f%x(:,i-1,j,k)-old%x(:,i-1,j,k))- &
          g%x_face_fraction(i,j,k)*(f%x(:,i,j,k)-old%x(:,i,j,k)))+ &
          g%dx*g%dz*(g%y_face_fraction(i,j-1,k)*(f%y(:,i,j-1,k)-old%y(:,i,j-1,k))- &
          g%y_face_fraction(i,j,k)*(f%y(:,i,j,k)-old%y(:,i,j,k)))+ &
          g%dx*g%dy*(g%z_face_fraction(i,j,k-1)*(f%z(:,i,j,k-1)-old%z(:,i,j,k-1))- &
          g%z_face_fraction(i,j,k)*(f%z(:,i,j,k)-old%z(:,i,j,k)))
        candidate%fields(node)%values(:,i,j,k)=candidate%fields(node)%values(:,i,j,k)+change/volume/history%interval
      enddo
      enddo
      enddo
      end associate
      if(history%has_internal_energy) then
        associate(g=>geometry(node),f=>internal_energy_flux(node),old=>history%internal_energy_flux(node))
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(g%volume_fraction(i,j,k)<=0) cycle
          volume=g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz
          held_change=g%dy*g%dz*(g%x_face_fraction(i-1,j,k)*(f%x(1,i-1,j,k)-old%x(1,i-1,j,k))- &
            g%x_face_fraction(i,j,k)*(f%x(1,i,j,k)-old%x(1,i,j,k)))+ &
            g%dx*g%dz*(g%y_face_fraction(i,j-1,k)*(f%y(1,i,j-1,k)-old%y(1,i,j-1,k))- &
            g%y_face_fraction(i,j,k)*(f%y(1,i,j,k)-old%y(1,i,j,k)))+ &
            g%dx*g%dy*(g%z_face_fraction(i,j,k-1)*(f%z(1,i,j,k-1)-old%z(1,i,j,k-1))- &
            g%z_face_fraction(i,j,k)*(f%z(1,i,j,k)-old%z(1,i,j,k)))
          if(held_change/=0) held_candidate%fields(node)%values(1,i,j,k)= &
            held_candidate%fields(node)%values(1,i,j,k)+held_change/volume/history%interval
        enddo
        enddo
        enddo
        end associate
        valid=valid.and.all(ieee_is_finite(held_candidate%fields(node)%values))
      endif
      valid=valid.and.all(ieee_is_finite(candidate%fields(node)%values))
    enddo
    call agree(valid,comm,ready)
    if(.not.ready) return
    if(history%has_internal_energy) call move_alloc(held_candidate%fields,internal_energy_rate%fields)
    call move_alloc(candidate%fields,rates%fields)
    ok=.true.
  end subroutine reconstruct_mpi_eb_mol_interval_3d

  ! Accepted interval fields are sparse by owner and selected source level.
  ! Their inactive entries must remain absent, including on entirely idle ranks.
  logical function valid_held_interval_field(field,tree,owners,depth,rank,components,padding) result(valid)
    type(amr_tree_state_3d), intent(in) :: field
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),depth,rank
    integer, optional, intent(in) :: components,padding
    integer :: p,ncomp,grow
    logical :: expected
    ncomp=1
    if(present(components)) ncomp=components
    grow=0
    if(present(padding)) grow=padding
    valid=.false.
    if(.not.allocated(field%fields)) return
    if(size(field%fields)/=size(tree%nodes).or.size(owners)/=size(tree%nodes)) return
    if(lbound(field%fields,1)/=1) return
    do p=1,size(tree%nodes)
      expected=owners(p)==rank.and.(depth<0.or.tree%nodes(p)%level==depth)
      if(allocated(field%fields(p)%values).neqv.expected) return
      if(.not.expected) cycle
      if(any(shape(field%fields(p)%values)/=[ncomp,tree%nodes(p)%shape+2*grow])) return
      if(any(lbound(field%fields(p)%values)/=1)) return
      if(.not.all(ieee_is_finite(field%fields(p)%values))) return
    enddo
    valid=.true.
  end function valid_held_interval_field


  logical function valid_mpi_mol_native_provider_cache_3d(tree,owners,geometry,history,nv,rank) result(valid)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nv,rank
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(mpi_eb_mol_interval_3d), intent(in) :: history
    integer :: stage,kind,p,d(3)
    logical :: expected
    valid=.false.
    if(history%native_iterations>0) then
      if(.not.allocated(history%external_stages)) return
      if(lbound(history%external_stages,1)/=0.or.ubound(history%external_stages,1)/=history%native_iterations) return
      do stage=0,history%native_iterations
        if(.not.valid_held_interval_field(history%external_stages(stage),tree,owners,-1,rank,nv)) return
      enddo
    else
      if(allocated(history%external_stages)) return
    endif
    expected=history%spatial_diffusion_is_source.and.any(history%native_providers)
    if(allocated(history%provider_stages).neqv.expected) return
    if(allocated(history%native_provider_geometry).neqv.expected) return
    if(.not.expected) then
      valid=.true.;return
    endif
    if(size(history%native_provider_geometry)/=size(tree%nodes)) return
    do p=1,size(tree%nodes)
      if(.not.history%native_provider_geometry(p)%is_valid().or..not.geometry(p)%is_valid()) return
      if(any([geometry(p)%dx,geometry(p)%dy,geometry(p)%dz, &
        geometry(p)%x_lower,geometry(p)%y_lower,geometry(p)%z_lower]/= &
        [history%native_provider_geometry(p)%dx,history%native_provider_geometry(p)%dy, &
        history%native_provider_geometry(p)%dz,history%native_provider_geometry(p)%x_lower, &
        history%native_provider_geometry(p)%y_lower,history%native_provider_geometry(p)%z_lower])) return
      if(.not.matching_eb_metrics_3d(geometry(p),history%native_provider_geometry(p))) return
    enddo
    if(any(lbound(history%provider_stages)/=[1,0]).or. &
      any(ubound(history%provider_stages)/=[2,history%native_iterations])) return
    do stage=0,history%native_iterations
      do kind=1,2
        associate(saved=>history%provider_stages(kind,stage))
        expected=history%native_providers(kind).and.stage<=1
        if(saved%provider_accepted.neqv.expected) return
        if(.not.expected) then
          if(allocated(saved%provider_rate%fields).or.allocated(saved%provider_halos%fields).or. &
            allocated(saved%provider_faces)) return
          cycle
        endif
        if(saved%index/=stage.or.saved%provider_kind/=kind) return
        if(.not.ieee_is_finite(saved%provider_limiter)) return
        if(saved%provider_limiter<0.or.saved%provider_limiter>1) return
        if(.not.valid_held_interval_field(saved%provider_rate,tree,owners,-1,rank,nv)) return
        if(.not.allocated(saved%provider_faces)) return
        if(size(saved%provider_faces)/=size(tree%nodes).or.lbound(saved%provider_faces,1)/=1) return
        if(kind==1) then
          if(.not.valid_held_interval_field(saved%provider_halos,tree,owners,history%depth,rank,nv,padding=6)) return
        else
          if(allocated(saved%provider_halos%fields)) return
        endif
        do p=1,size(tree%nodes)
          if(owners(p)/=rank.or.tree%nodes(p)%level/=history%depth) cycle
          if(.not.valid_eb_integrated_flux_3d(saved%provider_faces(p),geometry(p),nv)) return
          if(kind/=1) cycle
          d=tree%nodes(p)%shape
          if(any(saved%provider_halos%fields(p)%values(:,7:d(1)+6,7:d(2)+6,7:d(3)+6)/= &
            saved%provider_rate%fields(p)%values)) return
        enddo
        end associate
      enddo
    enddo
    valid=.true.
  end function valid_mpi_mol_native_provider_cache_3d

  ! Replay one accepted ordinary held MOL interval using an already assembled
  ! paired face family. The full clock group is restored by the recursive
  ! driver; this operator admits the local solve and all its ledger/history
  ! changes atomically. Held spatial RHS never replaces gas reaction forcing.
  subroutine replay_mpi_eb_mol_held_interval_3d(species,reactions,tree,owners,geometry,history,flux,held_flux, &
      state,temperature,internal_energy,reaction_source,energy_correction,chemical_change,rtol,atol,comm,ok, &
      evaluate_reaction,frame,heat_release,source_funding_stock,source_funding_impulse, &
      evaluate_reaction_admission,reaction_admission_status)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(mpi_eb_mol_interval_3d), intent(inout) :: history
    type(eb_integrated_flux_3d), intent(in) :: flux(:),held_flux(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature,internal_energy,reaction_source
    real(dp), intent(inout) :: energy_correction,chemical_change(:)
    real(dp), intent(in) :: rtol,atol
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(forced_reactor_admission), optional :: evaluate_reaction_admission
    integer, optional, intent(out) :: reaction_admission_status
    integer :: admission_state,admission_global,admission_ierr
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(inout) :: heat_release,source_funding_impulse
    type(amr_tree_state_3d), optional, intent(in) :: source_funding_stock
    type(mpi_eb_mol_interval_3d) :: record
    type(amr_tree_state_3d) :: old,old_t,old_held,work,work_t,work_held,work_ir,rates,held_rates,ledger,reset_delta
    type(amr_tree_state_3d) :: projection_delta
    real(dp) :: projection(1),transported
    real(dp) :: initial_scalar_reset
    type(amr_tree_state_3d), allocatable :: work_heat
    real(dp) :: initial_reset,final_reset,correction,amount(reactive_nvar(size(species))),volume
    real(dp) :: divergence(reactive_nvar(size(species))),prefix_min,prefix_max
    integer :: nv,rank,status,p,x,y,z,c,depth,ios,options,options_min,options_max,face_depth
    logical :: valid,ready
    admission_state=reaction_admission_not_attempted
    if(present(reaction_admission_status)) then
      reaction_admission_status=admission_state
      if(present(evaluate_reaction_admission)) reaction_admission_status=reaction_admission_unknown
    endif
    ok=.false.
    call run_replay()
    call MPI_Allreduce(admission_state,admission_global,1,MPI_INTEGER,MPI_MAX,comm,admission_ierr)
    if(admission_ierr/=MPI_SUCCESS) then
      ok=.false.;admission_global=reaction_admission_unknown
    endif
    if(.not.ok.and.present(evaluate_reaction_admission).and.admission_global==reaction_admission_not_attempted) &
      admission_global=reaction_admission_unknown
    if(ok.and.admission_global>=0.and.admission_global/=reaction_admission_provider_refused) &
      admission_global=reaction_admission_accepted
    if(present(reaction_admission_status)) reaction_admission_status=admission_global
  contains
    subroutine run_replay()
    ok=.false.;nv=reactive_nvar(size(species))
    options=merge(1,0,present(heat_release))+merge(2,0,present(frame))+merge(4,0,present(evaluate_reaction))+ &
      merge(8,0,present(source_funding_stock))+merge(16,0,present(source_funding_impulse))+ &
      merge(32,0,present(evaluate_reaction_admission))
    call MPI_Allreduce(options,options_min,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(options,options_max,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.options_min/=options_max) return
    call agree(.not.present(evaluate_reaction_admission).or.present(evaluate_reaction),comm,ready)
    if(.not.ready) return
    call agree(history%has_internal_energy.and.history%redistribution=='none'.and. &
      .not.history%prepared_redistribution.and. &
      size(chemical_change)==nv.and.ieee_is_finite(energy_correction),comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(energy_correction,prefix_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(energy_correction,prefix_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.prefix_min/=prefix_max) return
    call reconstruct_mpi_eb_mol_interval_3d(tree,owners,geometry,history,flux,nv,rates,comm,ready, &
      internal_energy_flux=held_flux,internal_energy_rate=held_rates)
    if(.not.ready) return
    call agree((present(source_funding_stock).eqv.present(source_funding_impulse)).and. &
      (present(source_funding_stock).eqv.history%has_source_funding),comm,ready)
    if(.not.ready) return
    if(present(source_funding_stock)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,source_funding_stock,nv+1,comm,ready)
      if(.not.ready) return
    endif
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    depth=history%depth;face_depth=depth
    if(history%has_source_funding) face_depth=-1
    valid=valid_held_interval_field(history%mass_source,tree,owners,depth,rank,nv).and. &
      valid_held_interval_field(history%initial,tree,owners,depth,rank,nv).and. &
      valid_held_interval_field(history%initial_temperature,tree,owners,depth,rank).and. &
      valid_held_interval_field(history%local_state,tree,owners,depth,rank,nv).and. &
      valid_held_interval_field(history%local_temperature,tree,owners,depth,rank).and. &
      valid_held_interval_field(history%initial_reaction_source,tree,owners,depth,rank,nv-4).and. &
      valid_held_interval_field(history%reaction_source,tree,owners,depth,rank,nv-4)
    call agree(valid,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,ready)
    if(.not.ready) return
    if(history%has_source_funding) then
      valid=.true.
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        valid=valid.and.all(source_funding_stock%fields(p)%values==history%source_funding_stock%fields(p)%values)
        if(tree%nodes(p)%level==depth) cycle
        valid=valid.and.all(state%fields(p)%values==history%clock_initial%fields(p)%values).and. &
          all(internal_energy%fields(p)%values==history%clock_internal_energy%fields(p)%values)
        valid=valid.and.all(flux(p)%x==history%spatial_flux(p)%x).and. &
          all(flux(p)%y==history%spatial_flux(p)%y).and.all(flux(p)%z==history%spatial_flux(p)%z).and. &
          all(held_flux(p)%x==history%internal_energy_flux(p)%x).and. &
          all(held_flux(p)%y==history%internal_energy_flux(p)%y).and. &
          all(held_flux(p)%z==history%internal_energy_flux(p)%z)
      enddo
      call agree(valid,comm,ready);if(.not.ready) return
    endif
    old=state;old_t=temperature;old_held=internal_energy;work_ir=reaction_source;record=history
    call reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,old,old_t,old_held,reset_delta, &
      history%internal_energy_policy,depth,initial_reset,comm,ready,frame, &
      internal_energy_change=initial_scalar_reset,reset_level=-2)
    if(.not.ready) return
    ! Freeze the canonical selected clock after its start repair. Restoring this
    ! record later retains that repair in the prefix instead of charging it twice.
    record%clock_initial=old;record%clock_temperature=old_t;record%clock_internal_energy=old_held
    record%clock_reaction_source=work_ir;record%clock_energy_correction=energy_correction+initial_reset
    if(present(heat_release)) then
      allocate(work_heat,source=heat_release,stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      record%clock_heat=work_heat
    else
      if(allocated(record%clock_heat%fields)) deallocate(record%clock_heat%fields)
    endif
    ! Direct volume sources avoid cancellation of resolved trace rates. Energy
    ! and momentum retain the wall/non-face family from paired reconstruction.
    do p=1,size(owners)
      if(owners(p)/=rank.or.tree%nodes(p)%level/=depth) cycle
      associate(g=>geometry(p),f=>flux(p))
      do z=1,g%nz
      do y=1,g%ny
      do x=1,g%nx
        if(g%volume_fraction(x,y,z)<=0) cycle
        volume=g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz
        divergence=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(:,x-1,y,z)- &
          g%x_face_fraction(x,y,z)*f%x(:,x,y,z))+ &
          g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(:,x,y-1,z)- &
          g%y_face_fraction(x,y,z)*f%y(:,x,y,z))+ &
          g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(:,x,y,z-1)- &
          g%z_face_fraction(x,y,z)*f%z(:,x,y,z))
        do c=1,nv
          if(c>=2.and.c<=5) cycle
          rates%fields(p)%values(c,x,y,z)=history%mass_source%fields(p)%values(c,x,y,z)+ &
            divergence(c)/volume/history%interval
        enddo
      enddo
      enddo
      enddo
      end associate
    enddo
    work=old;work_t=old_t;work_held=old_held
    if(present(evaluate_reaction_admission)) then
    call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old,old_t,rates, &
      history%interval,rtol,atol,work,work_t,work_ir,comm,ready,target_depth=depth, &
      evaluate_reaction=replay_admission_adapter,frame=frame,heat_release=work_heat, &
      internal_energy=work_held,old_internal_energy=old_held)
    else
    call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old,old_t,rates, &
      history%interval,rtol,atol,work,work_t,work_ir,comm,ready,target_depth=depth, &
      evaluate_reaction=evaluate_reaction,frame=frame,heat_release=work_heat, &
      internal_energy=work_held,old_internal_energy=old_held)
    endif
    if(.not.ready) return
    call reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,work,work_t,work_held,reset_delta, &
      history%internal_energy_policy,depth,final_reset,comm,ready,frame)
    if(.not.ready) return
    ledger=rates
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      ledger%fields(p)%values=0
      if(tree%nodes(p)%level/=depth) cycle
      ledger%fields(p)%values=work%fields(p)%values-(old%fields(p)%values+history%interval*rates%fields(p)%values)
      record%initial%fields(p)=old%fields(p);record%initial_temperature%fields(p)=old_t%fields(p)
      record%initial_internal_energy%fields(p)=old_held%fields(p)
      record%initial_reaction_source%fields(p)=reaction_source%fields(p)
      record%forcing%fields(p)=rates%fields(p)
      record%local_state%fields(p)=work%fields(p);record%local_temperature%fields(p)=work_t%fields(p)
      record%local_internal_energy%fields(p)=work_held%fields(p)
      record%reaction_source%fields(p)=work_ir%fields(p)
      record%internal_energy_spatial_rate%fields(p)=held_rates%fields(p)
    enddo
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,ledger,amount,comm,ready)
    if(.not.ready) return
    amount(5)=amount(5)-final_reset
    correction=energy_correction+initial_reset+history%forcing_energy_correction+final_reset
    call agree(ieee_is_finite(correction).and.all(ieee_is_finite(amount)),comm,ready)
    if(.not.ready) return
    record%replayed=.true.;record%initial_energy_correction=initial_reset
    record%numerical_energy_correction=initial_reset+history%forcing_energy_correction+final_reset
    projection_delta=old_held;valid=.true.
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      projection_delta%fields(p)%values=0
      if(tree%nodes(p)%level/=depth) cycle
      do z=1,geometry(p)%nz
      do y=1,geometry(p)%ny
      do x=1,geometry(p)%nx
        if(geometry(p)%volume_fraction(x,y,z)<=0.or.covered(tree,p,[x,y,z])) cycle
        transported=old_held%fields(p)%values(1,x,y,z)+history%interval*held_rates%fields(p)%values(1,x,y,z)
        call compute_reactive_eb_held_projection(internal_energy%fields(p)%values(1,x,y,z), &
          old_held%fields(p)%values(1,x,y,z),transported,work_held%fields(p)%values(1,x,y,z), &
          projection_delta%fields(p)%values(1,x,y,z),ready)
        valid=valid.and.ready
      enddo
      enddo
      enddo
    enddo
    call agree(valid,comm,ready)
    if(.not.ready) return
    projection=0
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,projection_delta,projection,comm,ready)
    if(.not.ready) return
    record%internal_energy_projection_change=projection(1)
    record%initial_internal_energy_change=initial_scalar_reset
    record%held_clock_start_change=initial_scalar_reset;record%held_clock_end_change=0
    ! This replacement ends at its local solve; any following hierarchy
    ! reset is owned by the surrounding clock group, not the previous record.
    record%internal_energy_closure_change=0
    record%held_clock_events_closed=.false.;record%held_clock_prefixes_valid=.false.
    record%held_clock_projection_before=0;record%held_clock_reset_before=0
    record%held_clock_projection_after=0;record%held_clock_reset_after=0
    record%chemical_change=amount
    call copy_mol_level_flux(tree,owners,rank,face_depth,flux,record%spatial_flux,comm,ready)
    if(.not.ready) return
    record%synchronized_flux=record%spatial_flux
    call copy_mol_level_flux(tree,owners,rank,face_depth,held_flux,record%internal_energy_flux,comm,ready)
    if(.not.ready) return
    call move_alloc(work%fields,state%fields);call move_alloc(work_t%fields,temperature%fields)
    call move_alloc(work_held%fields,internal_energy%fields);call move_alloc(work_ir%fields,reaction_source%fields)
    if(present(heat_release)) call move_alloc(work_heat%fields,heat_release%fields)
    if(present(source_funding_impulse)) source_funding_impulse=record%source_funding_impulse
    call move_mol_interval(record,history)
    energy_correction=correction;chemical_change=amount;ok=.true.
    end subroutine run_replay
    subroutine replay_admission_adapter(table,mechanism,energy,source_species,source_energy,duration,rtol,atol, &
        density,temperature,steps,good,maximum_steps)
      type(nasa7_species), intent(in) :: table(:)
      type(elementary_reaction), intent(in) :: mechanism(:)
      real(dp), intent(in) :: energy,source_species(:),source_energy,duration,rtol,atol
      real(dp), intent(inout) :: density(:),temperature
      integer, intent(out) :: steps
      logical, intent(out) :: good
      integer, optional, intent(in) :: maximum_steps
      logical :: refused
      call evaluate_classified_reactor(evaluate_reaction,table,mechanism,energy,source_species,source_energy, &
        duration,rtol,atol,density,temperature,steps,good,refused,maximum_steps, &
        evaluate_admission=evaluate_reaction_admission,admission_status=admission_state)
    end subroutine replay_admission_adapter
  end subroutine replay_mpi_eb_mol_held_interval_3d

  ! Re-evaluate every supplied prepared MOL stage through native StateRedist.
  ! The caller may replace state/raw-RHS halos when rebuilding a clock group;
  ! accepted external volume forcing remains frozen in the interval archive.
  ! Intermediate new stages contribute coarse deltas, while only the final
  ! new stage contributes the fine delta and the accepted mean spatial RHS.
  ! No chemistry, face synchronization, callback or caller publication occurs
  ! here. Rates/registers are transactional, including collective rejection.
  subroutine reconstruct_mpi_eb_mol_state_redist_interval_3d(tree,owners,geometry,extended_geometry,history, &
      stages,nvar,rates,registers,comm,ok,boundary_faces,internal_energy_rate,internal_energy_registers, &
      register_provisional_coarse)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    type(mpi_eb_mol_interval_3d), intent(in) :: history
    type(mpi_eb_mol_state_redist_stage_3d), intent(in) :: stages(0:)
    type(amr_tree_state_3d), intent(inout) :: rates
    type(state_redist_register_3d), allocatable, intent(inout) :: registers(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), intent(in) :: boundary_faces(6)
    ! Explicit true returns the pinned chronological all-coarse/old-final-fine
    ! registration. RHS reconstruction and the archive's state-update policy
    ! are unchanged. Absent retains the existing archive policy.
    logical, optional, intent(in) :: register_provisional_coarse
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: internal_energy_registers(:)
    type(amr_tree_state_3d) :: spatial,old_spatial,accepted,held_spatial,held_old_spatial
    type(state_redist_register_3d), allocatable :: held_current(:),held_old(:),held_new(:)
    type(state_redist_register_3d), allocatable :: old_registers(:),new_registers(:),current_registers(:)
    integer :: rank,status,allocation_status,p,n(3),stage,steps,first_depth,last_depth,depth,held_options
    type(state_redist_register_3d), allocatable :: registered(:),held_registered(:)
    integer :: controls(8),low(8),high(8),registration_mode
    logical :: accumulate_provisional
    real(dp) :: weight_min,weight_max
    logical :: valid,ready
    ok=.false.;steps=ubound(stages,1)
    registration_mode=-1
    if(present(register_provisional_coarse)) registration_mode=merge(1,0,register_provisional_coarse)
    accumulate_provisional=registration_mode<0.and..not.history%accepted_coarse_stages
    held_options=merge(1,0,history%has_internal_energy)+merge(2,0,present(internal_energy_rate))+ &
      merge(4,0,present(internal_energy_registers))
    controls=[nvar,history%depth,steps,held_options,merge(1,0,history%accepted_coarse_stages), &
      merge(1,0,history%spatial_diffusion_is_source),merge(1,0,history%auxiliary_has_internal_energy),registration_mode]
    call MPI_Allreduce(controls,low,8,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,8,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    valid=nvar>=6.and.steps>=1.and.history%redistribution=='state_redist'.and.history%prepared_redistribution
    valid=valid.and.history%native_iterations==0.and. &
      .not.allocated(history%auxiliary_stages).and..not.allocated(history%provider_stages).and. &
      .not.allocated(history%external_stages).and. &
      .not.allocated(history%native_provider_geometry).and..not.any(history%native_providers)
    valid=valid.and..not.history%has_source_funding.and. &
      .not.allocated(history%source_funding_stock%fields).and. &
      .not.allocated(history%source_funding_impulse%fields)
    valid=valid.and.(held_options==0.or.held_options==7)
    valid=valid.and.(.not.history%accepted_coarse_stages.or.history%has_internal_energy)
    valid=valid.and.(.not.history%auxiliary_has_internal_energy.or. &
      (history%has_internal_energy.and.history%has_auxiliary_source))
    valid=valid.and.history%depth>=-1.and.history%depth<tree%level_count()
    valid=valid.and.ieee_is_finite(history%interval).and.history%interval>0
    valid=valid.and.size(owners)==size(tree%nodes).and.size(geometry)==size(tree%nodes).and. &
      size(extended_geometry)==size(tree%nodes).and.allocated(history%volume_source%fields)
    if(allocated(history%volume_source%fields)) valid=valid.and.size(history%volume_source%fields)==size(tree%nodes)
    valid=valid.and.allocated(history%forcing%fields).and.allocated(history%old_spatial_rate%fields).and. &
      allocated(history%spatial_rate%fields)
    if(allocated(history%forcing%fields)) valid=valid.and.size(history%forcing%fields)==size(tree%nodes)
    if(allocated(history%old_spatial_rate%fields)) valid=valid.and.size(history%old_spatial_rate%fields)==size(tree%nodes)
    if(allocated(history%spatial_rate%fields)) valid=valid.and.size(history%spatial_rate%fields)==size(tree%nodes)
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    allocate(spatial%fields(size(owners)),stat=allocation_status)
    call agree(allocation_status==0,comm,ready)
    if(.not.ready) return
    if(history%has_internal_energy) allocate(held_spatial%fields(size(owners)),stat=allocation_status)
    call agree(allocation_status==0,comm,ready)
    if(.not.ready) return
    valid=.true.
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape
      allocate(spatial%fields(p)%values(nvar,n(1),n(2),n(3)),source=0._dp,stat=allocation_status)
      if(allocation_status/=0) then
        valid=.false.;exit
      endif
      if(history%has_internal_energy) then
        allocate(held_spatial%fields(p)%values(1,n(1),n(2),n(3)),source=0._dp,stat=allocation_status)
        if(allocation_status/=0) then
          valid=.false.;exit
        endif
      endif
      if(history%depth>=0.and.tree%nodes(p)%level/=history%depth) cycle
      if(.not.allocated(history%volume_source%fields(p)%values)) then
        valid=.false.;exit
      endif
      if(.not.allocated(history%forcing%fields(p)%values).or. &
        .not.allocated(history%old_spatial_rate%fields(p)%values).or. &
        .not.allocated(history%spatial_rate%fields(p)%values)) then
        valid=.false.;exit
      endif
      valid=all(shape(history%volume_source%fields(p)%values)==[nvar,n]).and. &
        all(ieee_is_finite(history%volume_source%fields(p)%values))
      valid=valid.and.all(shape(history%forcing%fields(p)%values)==[nvar,n]).and. &
        all(shape(history%old_spatial_rate%fields(p)%values)==[nvar,n]).and. &
        all(shape(history%spatial_rate%fields(p)%values)==[nvar,n])
      valid=valid.and.all(ieee_is_finite(history%forcing%fields(p)%values)).and. &
        all(ieee_is_finite(history%old_spatial_rate%fields(p)%values)).and. &
        all(ieee_is_finite(history%spatial_rate%fields(p)%values))
      if(.not.valid) exit
    enddo
    call agree(valid,comm,ready)
    if(.not.ready) return
    if(history%has_internal_energy) then
      valid=valid_held_interval_field(history%old_internal_energy_spatial_rate,tree,owners,history%depth,rank).and. &
        valid_held_interval_field(history%internal_energy_spatial_rate,tree,owners,history%depth,rank)
      call agree(valid,comm,ready)
      if(.not.ready) return
    endif
    first_depth=max(0,history%depth);last_depth=first_depth
    if(history%depth<0) last_depth=tree%level_count()-1
    do stage=0,steps
      valid=stages(stage)%source_stage==stage.and.stages(stage)%depth==history%depth.and. &
        stages(stage)%interval==history%interval.and.stages(stage)%register_scale==.5_dp
      valid=valid.and.(stages(stage)%has_internal_energy.eqv.history%has_internal_energy)
      valid=valid.and.stages(stage)%max_order==stages(0)%max_order.and. &
        stages(stage)%target_volume_fraction==stages(0)%target_volume_fraction
      valid=valid.and.ieee_is_finite(stages(stage)%spatial_weight).and.stages(stage)%spatial_weight>=0.and. &
        stages(stage)%spatial_weight<=1.and.stages(stage)%spatial_weight==stages(0)%spatial_weight
      call agree(valid,comm,ready)
      if(.not.ready) return
      ! Every rank must select the same family before entering its collectives.
      ! In particular, a unit-weight rank must not skip a mixed-family exchange.
      call MPI_Allreduce(stages(stage)%spatial_weight,weight_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(stages(stage)%spatial_weight,weight_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.weight_min/=weight_max) return
      if(allocated(current_registers)) deallocate(current_registers)
      do depth=first_depth,last_depth
        call compute_state_redist_mpi_level_stage_3d(tree,owners,extended_geometry,stages(stage)%state_halos, &
          stages(stage)%raw_rate_halos,nvar,depth,history%interval,spatial,current_registers,comm,ready,boundary_faces, &
          stages(stage)%target_volume_fraction,stages(stage)%max_order,register_scale=stages(stage)%register_scale)
        if(.not.ready) return
      enddo
      valid=.true.
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        if(history%depth>=0.and.tree%nodes(p)%level/=history%depth) cycle
        call post_eb_spatial_density_rate_3d(geometry(p),spatial%fields(p)%values,valid)
        if(.not.valid) exit
      enddo
      call agree(valid,comm,ready)
      if(.not.ready) return
      if(stages(stage)%spatial_weight<1) then
        call blend_mol_state_redist_stage_3d(tree,owners,geometry,history%depth,nvar,stages(stage)%spatial_weight, &
          stages(stage)%low_spatial_rate,spatial,current_registers,comm,ready)
        if(.not.ready) return
      endif
      if(registration_mode==1) then
        if(stage==0) then
          allocate(registered(size(owners)),stat=allocation_status)
          call agree(allocation_status==0,comm,ready);if(.not.ready) return
        endif
        valid=.true.
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          if(history%depth>=0.and.tree%nodes(p)%level/=history%depth) cycle
          n=tree%nodes(p)%shape
          call accumulate_state_redist_mol_registration_stage_3d(registered(p),current_registers(p), &
            stage,steps,nvar,n+12,ready)
          valid=valid.and.ready
        enddo
        call agree(valid,comm,ready);if(.not.ready) return
      endif
      if(history%has_internal_energy) then
        valid=valid_held_interval_field(stages(stage)%held_state_halos,tree,owners,history%depth,rank,padding=6).and. &
          valid_held_interval_field(stages(stage)%held_raw_rate_halos,tree,owners,history%depth,rank,padding=6)
        if(stages(stage)%spatial_weight<1) valid=valid.and. &
          valid_held_interval_field(stages(stage)%held_low_spatial_rate,tree,owners,history%depth,rank)
        call agree(valid,comm,ready)
        if(.not.ready) return
        if(allocated(held_current)) deallocate(held_current)
        do depth=first_depth,last_depth
          call compute_state_redist_mpi_level_stage_3d(tree,owners,extended_geometry,stages(stage)%held_state_halos, &
            stages(stage)%held_raw_rate_halos,1,depth,history%interval,held_spatial,held_current,comm,ready,boundary_faces, &
            stages(stage)%target_volume_fraction,stages(stage)%max_order,register_scale=stages(stage)%register_scale)
          if(.not.ready) return
        enddo
        if(stages(stage)%spatial_weight<1) then
          call blend_mol_state_redist_stage_3d(tree,owners,geometry,history%depth,1,stages(stage)%spatial_weight, &
            stages(stage)%held_low_spatial_rate,held_spatial,held_current,comm,ready)
          if(.not.ready) return
        endif
        if(registration_mode==1) then
          if(stage==0) then
            allocate(held_registered(size(owners)),stat=allocation_status)
            call agree(allocation_status==0,comm,ready);if(.not.ready) return
          endif
          valid=.true.
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            if(history%depth>=0.and.tree%nodes(p)%level/=history%depth) cycle
            n=tree%nodes(p)%shape
            call accumulate_state_redist_mol_registration_stage_3d(held_registered(p),held_current(p), &
              stage,steps,1,n+12,ready)
            valid=valid.and.ready
          enddo
          call agree(valid,comm,ready);if(.not.ready) return
        endif
        if(history%auxiliary_has_internal_energy) then
          valid=history%has_auxiliary_source.and.allocated(stages(stage)%auxiliary_rate%fields)
          if(valid) valid=size(stages(stage)%auxiliary_rate%fields)==size(owners).and. &
            lbound(stages(stage)%auxiliary_rate%fields,1)==1
          call agree(valid,comm,ready);if(.not.ready) return
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            if(history%depth>=0.and.tree%nodes(p)%level/=history%depth) cycle
            n=tree%nodes(p)%shape
            if(.not.allocated(stages(stage)%auxiliary_rate%fields(p)%values)) then
              valid=.false.;exit
            endif
            if(any(shape(stages(stage)%auxiliary_rate%fields(p)%values)/=[nvar+1,n]).or. &
              any(lbound(stages(stage)%auxiliary_rate%fields(p)%values)/=1)) then
              valid=.false.;exit
            endif
            if(.not.all(ieee_is_finite(stages(stage)%auxiliary_rate%fields(p)%values))) then
              valid=.false.;exit
            endif
          enddo
          call agree(valid,comm,ready);if(.not.ready) return
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            if(history%depth>=0.and.tree%nodes(p)%level/=history%depth) cycle
            if(any(abs(held_spatial%fields(p)%values)>huge(1._dp)- &
              abs(stages(stage)%auxiliary_rate%fields(p)%values(nvar+1:nvar+1,:,:,:)))) then
              valid=.false.;exit
            endif
            held_spatial%fields(p)%values=held_spatial%fields(p)%values+ &
              stages(stage)%auxiliary_rate%fields(p)%values(nvar+1:nvar+1,:,:,:)
            valid=all(ieee_is_finite(held_spatial%fields(p)%values));if(.not.valid) exit
          enddo
          call agree(valid,comm,ready);if(.not.ready) return
        endif
        if(stage==0) then
          held_old_spatial=held_spatial
          call move_alloc(held_current,held_old)
        else
          if(stage>1.and.accumulate_provisional) then
            valid=.true.
            do p=1,size(owners)
              if(owners(p)/=rank) cycle
              if(history%depth>=0.and.tree%nodes(p)%level/=history%depth) cycle
              n=tree%nodes(p)%shape
              call accumulate_state_redist_coarse_register_3d(held_current(p),held_new(p),1,n+12,ready)
              valid=valid.and.ready
            enddo
            call agree(valid,comm,ready)
            if(.not.ready) return
          endif
          call move_alloc(held_current,held_new)
        endif
      endif
      if(stage==0) then
        allocate(old_spatial%fields(size(owners)),stat=allocation_status)
        call agree(allocation_status==0,comm,ready)
        if(.not.ready) return
        valid=.true.
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          allocate(old_spatial%fields(p)%values,source=spatial%fields(p)%values,stat=allocation_status)
          if(allocation_status/=0) then
            valid=.false.;exit
          endif
        enddo
        call agree(valid,comm,ready)
        if(.not.ready) return
        call move_alloc(current_registers,old_registers)
      else
        if(stage>1.and.accumulate_provisional) then
          valid=.true.
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            if(history%depth>=0.and.tree%nodes(p)%level/=history%depth) cycle
            n=tree%nodes(p)%shape
            call accumulate_state_redist_coarse_register_3d(current_registers(p),new_registers(p),nvar,n+12,ready)
            if(.not.ready) then
              valid=.false.;exit
            endif
          enddo
          call agree(valid,comm,ready)
          if(.not.ready) return
        endif
        call move_alloc(current_registers,new_registers)
      endif
    enddo
    call move_alloc(spatial%fields,accepted%fields)
    valid=.true.
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      if(history%depth>=0.and.tree%nodes(p)%level/=history%depth) cycle
      ! Preserve the accepted anchor exactly when its actual native spatial
      ! operators are unchanged. A changed operator is composed directly with
      ! frozen volume forcing, avoiding an unfunded subtraction residual.
      if(all(old_spatial%fields(p)%values==history%old_spatial_rate%fields(p)%values).and. &
        all(.5_dp*(old_spatial%fields(p)%values+accepted%fields(p)%values)==history%spatial_rate%fields(p)%values)) then
        accepted%fields(p)%values=history%forcing%fields(p)%values
      else
        accepted%fields(p)%values=history%volume_source%fields(p)%values+ &
          .5_dp*(old_spatial%fields(p)%values+accepted%fields(p)%values)
      endif
      valid=all(ieee_is_finite(accepted%fields(p)%values))
      if(.not.valid) exit
      n=tree%nodes(p)%shape
      if(history%has_internal_energy) then
        held_spatial%fields(p)%values=.5_dp*(held_old_spatial%fields(p)%values+held_spatial%fields(p)%values)
        valid=all(ieee_is_finite(held_spatial%fields(p)%values))
        if(.not.valid) exit
        if(registration_mode/=1) then
        call accumulate_state_redist_register_3d(held_old(p),held_new(p),1,n+12,ready)
        if(.not.ready) then
          valid=.false.;exit
        endif
        endif
      endif
      if(registration_mode/=1) then
      call accumulate_state_redist_register_3d(old_registers(p),new_registers(p),nvar,n+12,ready)
      if(.not.ready) then
        valid=.false.;exit
      endif
      endif
    enddo
    call agree(valid,comm,ready)
    if(.not.ready) return
    if(registration_mode==1) then
      call move_alloc(registered,old_registers)
      if(history%has_internal_energy) call move_alloc(held_registered,held_old)
    endif
    if(history%has_internal_energy) then
      call move_alloc(held_spatial%fields,internal_energy_rate%fields)
      call move_alloc(held_old,internal_energy_registers)
    endif
    call move_alloc(accepted%fields,rates%fields)
    call move_alloc(old_registers,registers)
    ok=.true.
  end subroutine reconstruct_mpi_eb_mol_state_redist_interval_3d

  subroutine evaluate_classified_reactor(evaluator,species,reactions,energy,source_species,source_energy, &
      duration,rtol,atol,density,temperature,steps,ok,refused,maximum_steps,evaluate_admission,admission_status)
    procedure(forced_reactor_integrator) :: evaluator
    procedure(forced_reactor_admission), optional :: evaluate_admission
    integer, optional, intent(inout) :: admission_status
    integer :: result
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: energy,source_species(:),source_energy,duration,rtol,atol
    real(dp), intent(inout) :: density(:),temperature
    integer, intent(out) :: steps
    logical, intent(out) :: ok,refused
    integer, optional, intent(in) :: maximum_steps
    real(dp) :: original(size(density)),original_t,probe(size(density)),probe_t
    integer :: probe_steps
    logical :: admitted
    original=density;original_t=temperature;refused=.false.
    if(present(evaluate_admission)) then
      if(present(admission_status)) then
        if(admission_status==reaction_admission_provider_refused) then
          ok=.false.;steps=0;refused=.true.;return
        endif
      endif
      call evaluate_admission(species,reactions,energy,source_species,source_energy,duration,rtol,atol, &
        original,original_t,result)
      if(result<reaction_admission_accepted.or.result>reaction_admission_provider_refused) &
        result=reaction_admission_unknown
      if(present(admission_status)) admission_status=max(admission_status,result)
      if(result/=reaction_admission_accepted) then
        ok=.false.;steps=0;refused=result==reaction_admission_provider_refused;return
      endif
    endif
    call evaluator(species,reactions,energy,source_species,source_energy,duration,rtol,atol, &
      density,temperature,steps,ok,maximum_steps)
    if(ok) return
    if(present(evaluate_admission)) then
      refused=.true.
      if(present(admission_status)) admission_status=reaction_admission_provider_refused
      return
    endif
    ! A solver refusal of a native-admissible cell aborts the transaction.
    ! Probe the identical cell privately, retaining every failed callback output
    ! and never retrying that evaluator on the same rejected input.
    probe=original;probe_t=original_t
    call advance_forced_reactor_interval(species,reactions,energy,source_species,source_energy,duration,rtol,atol, &
      probe,probe_t,probe_steps,admitted,maximum_steps=maximum_steps)
    refused=admitted
  end subroutine evaluate_classified_reactor

  subroutine blend_mol_state_redist_stage_3d(tree,owners,geometry,depth,nvar,weight,low,high,registers,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),depth,nvar
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(in) :: weight
    type(amr_tree_state_3d), intent(in) :: low
    type(amr_tree_state_3d), intent(inout) :: high
    type(state_redist_register_3d), intent(inout) :: registers(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: rank,status,p,n(3)
    real(dp) :: weight_min,weight_max
    logical :: valid,ready
    ok=.false.
    valid=ieee_is_finite(weight).and.weight>=0.and.weight<=1
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(weight,weight_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(weight,weight_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.weight_min/=weight_max) return
    valid=allocated(low%fields).and.allocated(high%fields)
    if(valid) valid=size(low%fields)==size(owners).and.size(high%fields)==size(owners).and.size(registers)==size(owners)
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      if(depth>=0.and.tree%nodes(p)%level/=depth) cycle
      n=tree%nodes(p)%shape
      if(.not.allocated(low%fields(p)%values).or..not.allocated(high%fields(p)%values)) then
        valid=.false.;exit
      endif
      valid=all(shape(low%fields(p)%values)==[nvar,n]).and.all(shape(high%fields(p)%values)==[nvar,n])
      valid=valid.and.all(ieee_is_finite(low%fields(p)%values)).and.all(ieee_is_finite(high%fields(p)%values))
      valid=valid.and.registers(p)%is_valid(nvar,n+12)
      if(.not.valid) exit
    enddo
    call agree(valid,comm,ready)
    if(.not.ready) return
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      if(depth>=0.and.tree%nodes(p)%level/=depth) cycle
      if(weight==0) then
        high%fields(p)%values=low%fields(p)%values
      elseif(weight<1) then
        high%fields(p)%values=(1-weight)*low%fields(p)%values+weight*high%fields(p)%values
      endif
      if(weight<1) then
        if(allocated(registers(p)%coarse_delta)) registers(p)%coarse_delta=weight*registers(p)%coarse_delta
        if(allocated(registers(p)%fine_delta)) registers(p)%fine_delta=weight*registers(p)%fine_delta
      endif
      if(nvar>=6) call post_eb_spatial_density_rate_3d(geometry(p),high%fields(p)%values,valid)
      if(.not.valid) exit
    enddo
    call agree(valid,comm,ok)
  end subroutine blend_mol_state_redist_stage_3d

  subroutine move_mol_interval(source,destination)
    type(mpi_eb_mol_interval_3d), intent(inout) :: source,destination
    destination%depth=source%depth;destination%time=source%time;destination%interval=source%interval
    destination%separate_momentum=source%separate_momentum
    destination%consistent_velocity=source%consistent_velocity
    destination%paired_pressure=source%paired_pressure
    destination%connected_directions=source%connected_directions
    destination%native_iterations=source%native_iterations
    destination%native_spatial_weight=source%native_spatial_weight
    destination%native_force=source%native_force;destination%native_acceleration=source%native_acceleration
    destination%has_native_frame=source%has_native_frame;destination%native_frame=source%native_frame
    call move_alloc(source%auxiliary_stages,destination%auxiliary_stages)
    call move_alloc(source%external_stages,destination%external_stages)
    destination%native_providers=source%native_providers
    destination%native_provider_halo_width=source%native_provider_halo_width
    call move_alloc(source%provider_stages,destination%provider_stages)
    call move_alloc(source%native_provider_geometry,destination%native_provider_geometry)
    destination%replayed=source%replayed
    destination%redistribution=source%redistribution
    destination%prepared_redistribution=source%prepared_redistribution
    destination%accepted_coarse_stages=source%accepted_coarse_stages
    destination%spatial_diffusion_is_source=source%spatial_diffusion_is_source
    destination%auxiliary_has_internal_energy=source%auxiliary_has_internal_energy
    destination%routed_auxiliary=source%routed_auxiliary
    destination%has_auxiliary_source=source%has_auxiliary_source
    destination%has_internal_energy=source%has_internal_energy
    destination%has_source_funding=source%has_source_funding
    call move_eb_held_funding_context_3d(source%clock_source_funding_context, &
      destination%clock_source_funding_context)
    call move_alloc(source%source_funding_stock%fields,destination%source_funding_stock%fields)
    call move_alloc(source%source_funding_impulse%fields,destination%source_funding_impulse%fields)
    destination%internal_energy_policy=source%internal_energy_policy
    destination%initial_energy_correction=source%initial_energy_correction
    destination%numerical_energy_correction=source%numerical_energy_correction
    destination%internal_energy_projection_change=source%internal_energy_projection_change
    destination%internal_energy_closure_change=source%internal_energy_closure_change
    destination%initial_internal_energy_change=source%initial_internal_energy_change
    destination%held_clock_start_change=source%held_clock_start_change
    destination%held_clock_end_change=source%held_clock_end_change
    destination%held_clock_projection_before=source%held_clock_projection_before
    destination%held_clock_reset_before=source%held_clock_reset_before
    destination%held_clock_projection_after=source%held_clock_projection_after
    destination%held_clock_reset_after=source%held_clock_reset_after
    destination%held_clock_events_closed=source%held_clock_events_closed
    destination%held_clock_prefixes_valid=source%held_clock_prefixes_valid
    destination%forcing_energy_correction=source%forcing_energy_correction
    destination%clock_energy_correction=source%clock_energy_correction
    call move_alloc(source%clock_internal_energy%fields,destination%clock_internal_energy%fields)
    call move_alloc(source%initial_internal_energy%fields,destination%initial_internal_energy%fields)
    call move_alloc(source%local_internal_energy%fields,destination%local_internal_energy%fields)
    call move_alloc(source%old_internal_energy_spatial_rate%fields,destination%old_internal_energy_spatial_rate%fields)
    call move_alloc(source%internal_energy_spatial_rate%fields,destination%internal_energy_spatial_rate%fields)
    call move_alloc(source%old_internal_energy_flux,destination%old_internal_energy_flux)
    call move_alloc(source%internal_energy_flux,destination%internal_energy_flux)
    call move_alloc(source%initial%fields,destination%initial%fields)
    call move_alloc(source%initial_temperature%fields,destination%initial_temperature%fields)
    call move_alloc(source%forcing%fields,destination%forcing%fields)
    call move_alloc(source%mass_source%fields,destination%mass_source%fields)
    call move_alloc(source%local_state%fields,destination%local_state%fields)
    call move_alloc(source%local_temperature%fields,destination%local_temperature%fields)
    call move_alloc(source%reaction_source%fields,destination%reaction_source%fields)
    call move_alloc(source%initial_reaction_source%fields,destination%initial_reaction_source%fields)
    call move_alloc(source%routed_source%fields,destination%routed_source%fields)
    call move_alloc(source%old_volume_source%fields,destination%old_volume_source%fields)
    call move_alloc(source%volume_source%fields,destination%volume_source%fields)
    call move_alloc(source%old_spatial_rate%fields,destination%old_spatial_rate%fields)
    call move_alloc(source%spatial_rate%fields,destination%spatial_rate%fields)
    call move_alloc(source%non_face_spatial_rate%fields,destination%non_face_spatial_rate%fields)
    call move_alloc(source%old_redistribution_transfer%fields,destination%old_redistribution_transfer%fields)
    call move_alloc(source%redistribution_transfer%fields,destination%redistribution_transfer%fields)
    call move_alloc(source%synchronization_transfer%fields,destination%synchronization_transfer%fields)
    call move_alloc(source%prepared_stages,destination%prepared_stages)
    call move_alloc(source%clock_initial%fields,destination%clock_initial%fields)
    call move_alloc(source%clock_temperature%fields,destination%clock_temperature%fields)
    call move_alloc(source%clock_reaction_source%fields,destination%clock_reaction_source%fields)
    call move_alloc(source%clock_heat%fields,destination%clock_heat%fields)
    call move_alloc(source%ancestor_start%fields,destination%ancestor_start%fields)
    call move_alloc(source%ancestor_end%fields,destination%ancestor_end%fields)
    call move_alloc(source%ancestor_temperature_start%fields,destination%ancestor_temperature_start%fields)
    call move_alloc(source%ancestor_temperature_end%fields,destination%ancestor_temperature_end%fields)
    call move_alloc(source%ancestor_internal_energy_start%fields,destination%ancestor_internal_energy_start%fields)
    call move_alloc(source%ancestor_internal_energy_end%fields,destination%ancestor_internal_energy_end%fields)
    call move_alloc(source%ancestor_source_start%fields,destination%ancestor_source_start%fields)
    call move_alloc(source%ancestor_source_end%fields,destination%ancestor_source_end%fields)
    call move_alloc(source%fractions_start,destination%fractions_start)
    call move_alloc(source%fractions_end,destination%fractions_end)
    call move_alloc(source%clock_start,destination%clock_start)
    call move_alloc(source%clock_interval,destination%clock_interval)
    call move_alloc(source%open_ancestor_registers,destination%open_ancestor_registers)
    call move_alloc(source%open_ancestor_held_registers,destination%open_ancestor_held_registers)
    call move_alloc(source%old_redist_registers,destination%old_redist_registers)
    call move_alloc(source%redist_registers,destination%redist_registers)
    call move_alloc(source%old_internal_energy_registers,destination%old_internal_energy_registers)
    call move_alloc(source%internal_energy_registers,destination%internal_energy_registers)
    call move_alloc(source%old_flux,destination%old_flux)
    call move_alloc(source%spatial_flux,destination%spatial_flux)
    call move_alloc(source%synchronized_flux,destination%synchronized_flux)
    call move_alloc(source%chemical_change,destination%chemical_change)
    call move_alloc(source%wall_change,destination%wall_change)
    call move_alloc(source%source_change,destination%source_change)
    call move_alloc(source%auxiliary_change,destination%auxiliary_change)
    call move_alloc(source%old_predictor%fields,destination%old_predictor%fields)
    call move_alloc(source%new_predictor%fields,destination%new_predictor%fields)
  end subroutine move_mol_interval

  subroutine copy_mol_level_flux(tree,owners,rank,depth,source,destination,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),rank,depth
    type(eb_integrated_flux_3d), intent(in) :: source(:)
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: destination(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(eb_integrated_flux_3d), allocatable :: candidate(:)
    integer :: node,ios
    logical :: valid
    allocate(candidate(size(tree%nodes)),stat=ios)
    call agree(ios==0,comm,ok)
    if(.not.ok) return
    valid=.true.
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
      allocate(candidate(node)%x,source=source(node)%x,stat=ios)
      if(ios==0) allocate(candidate(node)%y,source=source(node)%y,stat=ios)
      if(ios==0) allocate(candidate(node)%z,source=source(node)%z,stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      endif
    enddo
    call agree(valid,comm,ok)
    if(ok) call move_alloc(candidate,destination)
  end subroutine copy_mol_level_flux

  subroutine agree(condition,comm,ok)
    logical, intent(in) :: condition
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: status
    call MPI_Allreduce(condition,ok,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if (status/=MPI_SUCCESS) ok=.false.
  end subroutine

  ! Owner-local counterpart of the serial fresh MOL history initializer.
  subroutine initialize_mpi_eb_tree_reaction_source_3d(species,reactions,tree,owners,geometry,state,temperature, &
      dt,subcycling,rtol,atol,reaction_source,comm,ok,evaluate_reaction,internal_energy,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: subcycling
    type(amr_tree_state_3d), intent(inout) :: reaction_source
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d) :: candidate,temps,rates,sources
    type(amr_tree_state_3d), allocatable :: held
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: node,depth,nv,rank,status,mode,low,high
    real(dp) :: level_dt
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species));mode=merge(1,0,subcycling)
    call MPI_Allreduce(mode,low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(mode,high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low/=high) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,valid)
    if(.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,valid)
    if(.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,reaction_source,nv-4,comm,valid)
    if(.not.valid) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    candidate=state;temps=temperature;rates=state;sources=reaction_source
    if(present(internal_energy)) then
      allocate(held);held=internal_energy
    endif
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      rates%fields(node)%values=0;sources%fields(node)%values=0
    enddo
    level_dt=dt
    do depth=0,tree%level_count()-1
      call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,state,temperature,rates,level_dt, &
        rtol,atol,candidate,temps,sources,comm,valid,target_depth=depth,evaluate_reaction=evaluate_reaction, &
        internal_energy=held,old_internal_energy=internal_energy,frame=frame)
      if(.not.valid) return
      if(depth==0) then
        call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,sources,nv-4,comm,valid)
        if(.not.valid) return
      endif
      if(subcycling.and.depth<tree%level_count()-1) level_dt=level_dt/real(ratios(depth+2),dp)
    enddo
    call move_alloc(sources%fields,reaction_source%fields);ok=.true.
  end subroutine initialize_mpi_eb_tree_reaction_source_3d

  subroutine react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old_state,old_temperature, &
      nonreacting_rates,dt,rtol,atol,state,temperature,reaction_source,comm,ok,target_depth, &
      evaluate_reaction,frame,heat_release,internal_energy,energy_policy,energy_correction,old_internal_energy)
    type(amr_tree_state_3d), optional, intent(inout) :: heat_release,internal_energy
    type(amr_tree_state_3d), optional, intent(in) :: old_internal_energy
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    real(dp), optional, intent(inout) :: energy_correction
    type(amr_tree_state_3d), allocatable :: heat,held
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: old_state,old_temperature,nonreacting_rates
    real(dp), intent(in) :: dt,rtol,atol
    type(amr_tree_state_3d), intent(inout) :: state,temperature,reaction_source
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d) :: candidate,temps,sources
    real(dp) :: controls(3),low(3),high(3),stamp
    real(dp) :: frame_values(4),frame_min(4),frame_max(4)
    real(dp) :: policy_value,policy_min,policy_max,initial_correction,correction_min,correction_max
    real(dp) :: local_correction,node_correction,total_correction
    integer :: policy_flags(3),policy_flags_min(3),policy_flags_max(3)
    integer :: nv,depth,counts(4),min_counts(4),max_counts(4),rank,status,node,n(3)
    integer :: supplied,minimum_supplied,maximum_supplied,axis_min,axis_max
    logical :: valid,ready
    real(dp), allocatable :: node_heat(:,:,:),node_old_held(:,:,:)
    integer :: heat_supplied,heat_min,heat_max,heat_status
    ok=.false.;nv=reactive_nvar(size(species));depth=-1;controls=[dt,rtol,atol]
    heat_supplied=merge(1,0,present(heat_release))
    call MPI_Allreduce(heat_supplied,heat_min,1,MPI_INTEGER,MPI_MIN,comm,heat_status)
    if(heat_status/=MPI_SUCCESS) return
    call MPI_Allreduce(heat_supplied,heat_max,1,MPI_INTEGER,MPI_MAX,comm,heat_status)
    if(heat_status/=MPI_SUCCESS.or.heat_min/=heat_max) return
    heat_supplied=merge(1,0,present(internal_energy))
    call MPI_Allreduce(heat_supplied,heat_min,1,MPI_INTEGER,MPI_MIN,comm,heat_status)
    if(heat_status/=MPI_SUCCESS) return
    call MPI_Allreduce(heat_supplied,heat_max,1,MPI_INTEGER,MPI_MAX,comm,heat_status)
    if(heat_status/=MPI_SUCCESS.or.heat_min/=heat_max) return
    heat_supplied=merge(1,0,present(old_internal_energy))
    call MPI_Allreduce(heat_supplied,heat_min,1,MPI_INTEGER,MPI_MIN,comm,heat_status)
    if(heat_status/=MPI_SUCCESS) return
    call MPI_Allreduce(heat_supplied,heat_max,1,MPI_INTEGER,MPI_MAX,comm,heat_status)
    if(heat_status/=MPI_SUCCESS.or.heat_min/=heat_max) return
    heat_supplied=merge(1,0,present(energy_correction))
    call MPI_Allreduce(heat_supplied,heat_min,1,MPI_INTEGER,MPI_MIN,comm,heat_status)
    if(heat_status/=MPI_SUCCESS) return
    call MPI_Allreduce(heat_supplied,heat_max,1,MPI_INTEGER,MPI_MAX,comm,heat_status)
    if(heat_status/=MPI_SUCCESS.or.heat_min/=heat_max) return
    heat_supplied=merge(1,0,present(energy_policy))
    call MPI_Allreduce(heat_supplied,heat_min,1,MPI_INTEGER,MPI_MIN,comm,heat_status)
    if(heat_status/=MPI_SUCCESS) return
    call MPI_Allreduce(heat_supplied,heat_max,1,MPI_INTEGER,MPI_MAX,comm,heat_status)
    if(heat_status/=MPI_SUCCESS.or.heat_min/=heat_max) return
    valid=(.not.present(energy_correction).or.present(internal_energy)).and. &
      (.not.present(energy_policy).or.present(energy_correction)).and. &
      (.not.present(old_internal_energy).or.present(internal_energy))
    if(present(energy_correction)) valid=valid.and.ieee_is_finite(energy_correction)
    call agree(valid,comm,ready)
    if(.not.ready) return
    initial_correction=0;local_correction=0
    if(present(energy_correction)) then
      initial_correction=energy_correction
      call MPI_Allreduce(initial_correction,correction_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(initial_correction,correction_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.correction_min/=correction_max) return
    endif
    if(present(energy_policy)) then
      policy_flags=[merge(1,0,energy_policy%allow_small_energy), &
        merge(1,0,energy_policy%allow_negative_energy),merge(1,0,energy_policy%update_total_from_internal)]
      policy_value=energy_policy%eta2
      call agree(ieee_is_finite(policy_value),comm,ready)
      if(.not.ready) return
      call MPI_Allreduce(policy_flags,policy_flags_min,3,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(policy_flags,policy_flags_max,3,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(policy_flags_min/=policy_flags_max)) return
      call MPI_Allreduce(policy_value,policy_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(policy_value,policy_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.policy_min/=policy_max) return
    endif
    supplied=merge(1,0,present(frame))
    call MPI_Allreduce(supplied,minimum_supplied,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(supplied,maximum_supplied,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.minimum_supplied/=maximum_supplied) return
    if(present(frame)) then
      call agree(frame%is_valid(),comm,ready)
      if(.not.ready) return
      call MPI_Allreduce(frame%axis,axis_min,1,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame%axis,axis_max,1,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS) return
      frame_values=[frame%angular_speed,frame%axis_point]
      call MPI_Allreduce(frame_values,frame_min,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame_values,frame_max,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.axis_min/=axis_max.or.any(frame_min/=frame_max)) return
    end if
    valid=all(ieee_is_finite(controls)).and.min(dt,rtol,atol)>0.and.size(reactions)>0
    if(present(target_depth)) then
      depth=target_depth
      valid=valid.and.depth>=0
    end if
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(controls,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    counts=[size(species),size(reactions),depth,merge(1,0,present(evaluate_reaction))]
    call MPI_Allreduce(counts,min_counts,4,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(counts,max_counts,4,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(min_counts/=max_counts)) return
    call validate_mpi_amr_tree_local_3d(tree,owners,old_state,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,old_temperature,1,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,nonreacting_rates,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,reaction_source,nv-4,comm,ready)
    if(.not.ready) return
    call agree(depth<tree%level_count().and.size(geometry)==size(tree%nodes),comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    candidate=old_state;temps=old_temperature;sources=reaction_source
    if(depth>=0) then
      call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ready)
      if(.not.ready) return
      candidate=state;temps=temperature
    end if
    if(present(internal_energy)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,valid)
      if(.not.valid) return
      allocate(held);held=internal_energy
    endif
    if(present(old_internal_energy)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,old_internal_energy,1,comm,valid)
      if(.not.valid) return
    endif
    valid=.true.
    if(present(heat_release)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,heat_release,1,comm,valid)
      if(.not.valid) return
      allocate(heat);heat=heat_release
    endif
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
      n=tree%nodes(node)%shape
      if(.not.geometry(node)%is_valid()) then
        valid=.false.;exit
      end if
      if(any([geometry(node)%nx,geometry(node)%ny,geometry(node)%nz]/=n)) then
        valid=.false.;exit
      end if
      stamp=tree_work_start()
      if(present(old_internal_energy)) node_old_held=old_internal_energy%fields(node)%values(1,:,:,:)
      if(present(heat_release)) allocate(node_heat(n(1),n(2),n(3)),source=0._dp)
      if(present(energy_correction)) then
        node_correction=0
        call react_sdc_eb_state_3d(species,reactions,old_state%fields(node)%values, &
          old_temperature%fields(node)%values(1,:,:,:),geometry(node),nonreacting_rates%fields(node)%values, &
          dt,rtol,atol,candidate%fields(node)%values,temps%fields(node)%values(1,:,:,:), &
          sources%fields(node)%values,valid, &
          evaluate_reaction=evaluate_reaction,frame=frame,heat_release=node_heat, &
          internal_energy_density=held%fields(node)%values(1,:,:,:), &
          energy_policy=energy_policy,energy_correction=node_correction,old_internal_energy_density=node_old_held)
        if(valid) local_correction=local_correction+node_correction
      else if(present(internal_energy)) then
        call react_sdc_eb_state_3d(species,reactions,old_state%fields(node)%values, &
          old_temperature%fields(node)%values(1,:,:,:),geometry(node),nonreacting_rates%fields(node)%values, &
          dt,rtol,atol,candidate%fields(node)%values,temps%fields(node)%values(1,:,:,:), &
          sources%fields(node)%values,valid, &
          evaluate_reaction=evaluate_reaction,frame=frame,heat_release=node_heat, &
          internal_energy_density=held%fields(node)%values(1,:,:,:),old_internal_energy_density=node_old_held)
      else
        call react_sdc_eb_state_3d(species,reactions,old_state%fields(node)%values, &
          old_temperature%fields(node)%values(1,:,:,:),geometry(node),nonreacting_rates%fields(node)%values, &
          dt,rtol,atol,candidate%fields(node)%values,temps%fields(node)%values(1,:,:,:), &
          sources%fields(node)%values,valid, &
          evaluate_reaction=evaluate_reaction,frame=frame,heat_release=node_heat, &
          old_internal_energy_density=node_old_held)
      end if
      call tree_work_stop(node,stamp)
      if(.not.valid) exit
      if(present(heat_release)) then
        heat%fields(node)%values(1,:,:,:)=node_heat
        deallocate(node_heat)
      endif
    end do
    call agree(valid,comm,ready)
    if(.not.ready) return
    if(present(energy_correction)) then
      call agree(ieee_is_finite(local_correction),comm,ready)
      if(.not.ready) return
      call MPI_Allreduce(local_correction,total_correction,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
      if(status/=MPI_SUCCESS) return
      call agree(ieee_is_finite(initial_correction+total_correction),comm,ready)
      if(.not.ready) return
      energy_correction=initial_correction+total_correction
    endif
    if(present(heat_release)) heat_release=heat
    if(present(internal_energy)) internal_energy=held
    state=candidate;temperature=temps;reaction_source=sources;ok=.true.
  end subroutine react_sdc_mpi_eb_tree_3d

  subroutine react_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,dt,rtol,atol,state,temperature, &
      chemical_change,comm,ok,integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(in) :: dt,rtol,atol
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(inout) :: chemical_change(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    type(amr_tree_state_3d) :: candidate,ct
    real(dp) :: work_stamp
    real(dp) :: before(size(chemical_change)),after(size(chemical_change)),parameters(3),low(3),high(3)
    logical, allocatable :: active(:,:,:)
    integer :: nv,rank,status,node,child,n(3),lo(3),hi(3),ios
    character(len=32) :: method,reference
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species));parameters=[dt,rtol,atol]
    valid=all(ieee_is_finite(parameters))
    if (valid) valid=all(parameters>0._dp)
    valid=valid.and.size(chemical_change)==nv.and.size(reactions)>0
    method='@default'
    if (present(integrator)) then
      valid=valid.and.len_trim(integrator)<=len(method)
      method=integrator
    end if
    call agree(valid,comm,ready)
    if (.not.ready) return
    call MPI_Allreduce(parameters,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(parameters,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.any(low/=high)) return
    reference=method
    call MPI_Bcast(reference,len(reference),MPI_CHARACTER,0,comm,status)
    call agree(status==MPI_SUCCESS.and.method==reference,comm,ready)
    if (.not.ready) return
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,state,before,comm,ready)
    if (.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    candidate=state;ct=temperature
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(active(n(1),n(2),n(3)),stat=ios)
      valid=ios==0
      if (.not.valid) exit
      active=geometry(node)%volume_fraction>0._dp
      do child=node+1,size(owners)
        if (tree%nodes(child)%parent/=node) cycle
        associate(p=>tree%nodes(child)%patch)
        lo=[p%coarse_i_lower,p%coarse_j_lower,p%coarse_k_lower]
        hi=[p%coarse_i_upper,p%coarse_j_upper,p%coarse_k_upper]
        end associate
        active(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      work_stamp=tree_work_start()
      if (any(active)) call advance_reactive_chemistry_3d(species,reactions,candidate%fields(node)%values, &
        ct%fields(node)%values(1,:,:,:),n(1),n(2),n(3),dt,rtol,atol,valid,active,integrator)
      call tree_work_stop(node,work_stamp)
      deallocate(active)
      if (.not.valid) exit
    end do
    call agree(valid,comm,ready)
    if (.not.ready) return
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,ready)
    if (.not.ready) return
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,candidate,after,comm,ready)
    if (.not.ready) return
    call agree(all(ieee_is_finite(after-before)),comm,ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(ct%fields,temperature%fields)
    chemical_change=after-before;ok=.true.
  end subroutine

  ! Transactional R/2 T/2 H T/2 R/2, common or recursive H/T, owner-local fields.
  subroutine advance_mpi_eb_tree_full_3d(species,reactions,transport,tree,owners,geometry,solver,dt, &
      chemistry_enabled,rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion, &
      state,temperature,boundary_change,wall_change,chemical_change,minimum_theta,comm,ok,integrator,subcycling, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: chemistry_enabled,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),minimum_theta
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    real(dp), optional, intent(in) :: difmag
    real(dp) :: coefficient,coefficient_low,coefficient_high
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    logical, optional, intent(in) :: subcycling
    logical :: recursive_steps
    type(amr_tree_state_3d) :: candidate,ct
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: exchange(size(boundary_change)),transport_sum(size(boundary_change)),theta,stage_theta
    integer :: flags(3),low(3),high(3),status,nv,half
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    valid=ieee_is_finite(dt).and.ieee_is_finite(coefficient).and.coefficient>=0
    if (valid) valid=dt>0._dp
    if(present(wall_temperature)) then
      valid=valid.and.present(boundary_faces)
      valid=valid.and.all(ieee_is_finite(wall_temperature)).and.all(wall_temperature>0)
    end if
    if(present(boundary_faces)) then
      if(any(boundary_faces=='isothermal_wall')) then
        valid=valid.and.present(wall_temperature).and.transport_enabled.and.conduction
      end if
    end if
    call agree(valid.and.size(boundary_change)==nv.and.size(wall_change)==nv.and.size(chemical_change)==nv,comm,ready)
    if (.not.ready) return
    call MPI_Allreduce(coefficient,coefficient_low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(coefficient,coefficient_high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.coefficient_low/=coefficient_high) return
    recursive_steps=.false.
    if (present(subcycling)) recursive_steps=subcycling
    flags=merge(1,0,[chemistry_enabled,transport_enabled,recursive_steps])
    call MPI_Allreduce(flags,low,3,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,3,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.any(low/=high)) return
    candidate=state;ct=temperature
    chemical=0._dp;transport_sum=0._dp;theta=1._dp
    do half=1,2
      if (half==1.and.chemistry_enabled) then
        call reaction_half(valid)
        if (.not.valid) return
      end if
      if (transport_enabled) then
        if (recursive_steps) then
          call advance_mpi_eb_tree_subcycle_transport_3d(species,transport,tree,owners,geometry,0.5_dp*dt, &
            viscosity,conduction,diffusion,barodiffusion,candidate,ct,exchange,stage_theta,comm,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,wall_temperature=wall_temperature)
        else
          call advance_mpi_eb_tree_transport_3d(species,transport,tree,owners,geometry,0.5_dp*dt, &
            viscosity,conduction,diffusion,barodiffusion,candidate,ct,exchange,stage_theta,comm,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,wall_temperature=wall_temperature)
        end if
        if (.not.valid) return
        transport_sum=transport_sum+exchange;theta=min(theta,stage_theta)
      end if
      if (half==1) then
        if (recursive_steps) then
          call advance_mpi_eb_tree_subcycle_hydro_3d(species,tree,owners,geometry,solver,dt,candidate,ct, &
            boundary,wall,comm,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state, &
            difmag=coefficient,extended_geometry=extended_geometry)
        else
          call advance_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,solver,dt,candidate,ct,boundary,wall,comm,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state, &
            difmag=coefficient,extended_geometry=extended_geometry)
        end if
        if (.not.valid) return
      else if (chemistry_enabled) then
        call reaction_half(valid)
        if (.not.valid) return
      end if
    end do
    boundary=boundary+transport_sum
    call agree(all(ieee_is_finite(boundary)).and.all(ieee_is_finite(chemical)),comm,valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(ct%fields,temperature%fields)
    boundary_change=boundary;wall_change=wall;chemical_change=chemical;minimum_theta=theta;ok=.true.
  contains
    subroutine reaction_half(good)
      logical, intent(out) :: good
      call react_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,0.5_dp*dt,rtol,atol, &
        candidate,ct,exchange,comm,good,integrator)
      if (good) chemical=chemical+exchange
    end subroutine
  end subroutine
  ! Common-step or selected-level SDC; stages and publication are collective.
  ! Selected levels export fluxes without initial restriction or final reflux.
  subroutine advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry, &
      state,temperature,extended_geometry,dt,iterations,chemistry,transport_flags,force,rtol,atol, &
      boundary_faces,boundary_change,wall_change,chemical_change,source_change,minimum_theta,comm,ok, &
      reaction_source,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
      target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start,fractions_end, &
      source_start,source_end,integrated_flux,old_predictor,new_predictor,redistribution,redistribution_transfer, &
      target_volume_fraction,auxiliary_source,auxiliary_change,gas_acceleration,predictor_builder,evaluate_diffusion, &
      evaluate_reaction,embedded_wall_temperature,embedded_noslip,mol_mode,mol_lax_friedrichs,frame, &
      state_redist_max_order,state_redist_registers,evaluate_spatial_diffusion,spatial_diffusion_halo_width, &
      constant_coefficients,heat_release,registered_flux,internal_energy,internal_energy_flux, &
      ancestor_internal_energy_start,ancestor_internal_energy_end,held_state_redist_registers, &
      accepted_forcing,accepted_mass_source,old_spatial_flux,interval_history,replay_interval,replay_spatial_weight, &
      provider_refused,energy_policy,energy_correction, &
      spatial_diffusion_is_source,auxiliary_has_internal_energy,internal_energy_pressure_work_change,spatial_flux_weight, &
      godunov_interval_history,godunov_replay_interval,godunov_replay_partial, &
      godunov_source_recovery,godunov_source_recovery_ready,source_funding_stock,source_funding_impulse, &
      internal_energy_projection_change,internal_energy_closure_change,freeze_native_external_sources,native_spatial_weight, &
      native_source_prefix,native_source_prefix_ready,native_source_replay, &
      separate_momentum,consistent_velocity,paired_pressure,connected_directions, &
      evaluate_reaction_admission,reaction_admission_status)
    procedure(forced_reactor_admission), optional :: evaluate_reaction_admission
    integer, optional, intent(out) :: reaction_admission_status
    integer :: admission_state,admission_global,admission_low,admission_high,admission_ierr,child_admission
    integer :: probe_admission,admission_cache_count,admission_cache_next
    integer(int64), allocatable :: admission_cache_packets(:,:)
    integer :: admission_cache_results(16)
    logical, optional, intent(in) :: separate_momentum,consistent_velocity,paired_pressure,connected_directions
    logical :: selected_momentum,selected_velocity,selected_pressure,selected_connected
    integer :: momentum_min,momentum_max
    type(eb_native_mol_source_prefix_3d), optional, intent(in) :: native_source_replay
    type(eb_native_mol_source_prefix_3d), optional, intent(inout) :: native_source_prefix
    logical, optional, intent(out) :: native_source_prefix_ready
    logical :: native_prefix_agreed,native_source_capture,native_prefix_scope
    integer :: native_prefix_observed
    type(amr_tree_state_3d) :: native_prefix_initial_reaction_source
    logical, optional, intent(in) :: freeze_native_external_sources
    real(dp), optional, intent(in) :: native_spatial_weight
    real(dp) :: native_mol_weight,native_mol_weight_min,native_mol_weight_max
    type(eb_integrated_flux_3d), allocatable :: native_mol_low_faces(:),native_mol_low_held_faces(:)
    logical :: native_external_frozen
    ! Selected clocks retain off-level gas/actual-UEINT inventory as a
    ! packed integrated impulse. Existing pending stock is read-only and is
    ! charged in donor admission, never added to the transported loan itself.
    type(amr_tree_state_3d), optional, intent(in) :: source_funding_stock
    type(amr_tree_state_3d), optional, intent(inout) :: source_funding_impulse
    type(amr_tree_state_3d) :: funding_receipt,funding_raw_h,godunov_predictor_forcing
    type(amr_tree_state_3d) :: godunov_predictor_hydro,godunov_predictor_held,godunov_predictor_work
    type(eb_integrated_flux_3d), allocatable :: godunov_predictor_faces(:),godunov_predictor_scalar_faces(:)
    logical :: native_held_funding,funding_corrector_stage
    integer :: funding_flags,funding_low,funding_high
    type(eb_godunov_interval_3d), optional, intent(in) :: godunov_replay_interval
    logical, optional, intent(in) :: godunov_replay_partial
    type(eb_godunov_interval_3d), optional, intent(inout) :: godunov_source_recovery
    logical, optional, intent(out) :: godunov_source_recovery_ready
    logical :: godunov_partial,godunov_stage_capture,godunov_recovery_agreed
    logical :: godunov_recovery_forbidden,godunov_source_attempt_incomplete,godunov_packet_pending
    integer :: godunov_source_stage
    type(eb_godunov_interval_3d), optional, intent(inout) :: godunov_interval_history
    type(eb_godunov_interval_3d), allocatable :: godunov_history
    integer :: history_status
    ! Native Godunov gas/UEINT face family; external volume forcing and pdivu
    ! are reevaluated at each actual SDC iterate, before StateRedist.
    real(dp), optional, intent(in) :: spatial_flux_weight
    ! Signed time-integrated -p div(u), before redistribution, on composite leaves.
    real(dp), optional, intent(inout) :: internal_energy_pressure_work_change
    ! Local scalar reassignment after the final accepted SDC iterate, before
    ! hierarchy closure. Picard iterates replace this impulse; rejected trials
    ! publish none. Initial scalar resets are included, gas resets are separate.
    real(dp), optional, intent(inout) :: internal_energy_projection_change
    ! Scalar reset after accepted common-clock spatial closure, separate from
    ! the last local SDC projection and from gas numerical-energy correction.
    real(dp), optional, intent(inout) :: internal_energy_closure_change
    real(dp) :: pressure_work_change,pressure_work_inventory(1),projection_inventory(1)
    real(dp) :: closure_inventory(1),scalar_reset_change,initial_scalar_reset
    type(amr_tree_state_3d) :: projection_delta,projection_transport
    type(amr_tree_state_3d) :: held_hydro,held_hydro_work
    type(eb_integrated_flux_3d), allocatable :: held_hydro_flux(:),held_hydro_baseline(:)
    type(state_redist_register_3d), allocatable :: held_hydro_registers(:),iteration_held_hydro_registers(:), &
      prior_held_hydro_registers(:)
    logical, optional, intent(in) :: auxiliary_has_internal_energy
    logical :: held_auxiliary
    integer :: naux,held_aux_low,held_aux_high,auxiliary_alloc_status
    logical, optional, intent(in) :: spatial_diffusion_is_source
    logical :: external_diffusion_source
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    real(dp), optional, intent(inout) :: energy_correction
    type(dual_energy_policy) :: held_policy
    type(amr_tree_state_3d) :: reset_delta,predictor_reset
    real(dp) :: initial_reset,predictor_reset_change,final_reset,post_reset,reset_change,numerical_change,accepted_predictor_reset
    real(dp) :: held_controls(2),held_controls_low(2),held_controls_high(2)
    integer :: held_policy_flags(3),held_policy_min(3),held_policy_max(3)
    ! Optional MOL interval history for conservative reaction/flux reconstruction.
    ! Rates are those of the final local reaction solve, before any hierarchy
    ! flux synchronization. mass_source contains only density/species volume
    ! forcing. old_spatial_flux retains f0 before any stage registration weights.
    ! All three publish with the accepted state, never on a rejected transaction.
    type(mpi_eb_mol_interval_3d), optional, intent(inout) :: interval_history
    type(mpi_eb_mol_interval_3d) :: pending_history
    ! Recompute native prepared MOL stages from the supplied fresh gas state.
    ! External volume families and the initial I_R predictor are frozen from
    ! the accepted archive; auxiliary providers are never called again.
    type(mpi_eb_mol_interval_3d), optional, intent(in) :: replay_interval
    ! One owner-neutral weight for all conservative spatial rates, faces and
    ! prepared register deltas. An absent weight retains the native operator.
    real(dp), optional, intent(in) :: replay_spatial_weight
    logical, optional, intent(out) :: provider_refused
    logical :: reactor_refused,reactor_refused_local,prepared_reaction_done
    procedure(forced_reactor_integrator), pointer :: prepared_integrator
    real(dp) :: spatial_weight,spatial_weight_min,spatial_weight_max
    real(dp) :: native_weight,native_weight_min,native_weight_max
    type(eb_integrated_flux_3d), allocatable :: prepared_low_faces(:),prepared_low_held_faces(:)
    type(amr_tree_state_3d), optional, intent(inout) :: accepted_forcing,accepted_mass_source
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: old_spatial_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: heat_release
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_internal_energy_start,ancestor_internal_energy_end
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: internal_energy_flux(:)
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: held_state_redist_registers(:)
    ! Time-integrated spatial face flux with pinned stage registration weights,
    ! before peer/coarse-fine synchronization. Separate from the physical flux
    ! used by the conservative state update and boundary inventory.
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: registered_flux(:)
    type(eb_integrated_flux_3d), allocatable :: registration(:),physical_flux(:)
    type(amr_tree_state_3d), allocatable :: heat
    real(dp), optional, intent(in) :: constant_coefficients(:)
    logical :: constant_agreed
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(rotating_frame_3d), optional, intent(in) :: frame
    procedure(mpi_eb_sdc_diffusion_source), optional :: evaluate_diffusion
    procedure(mpi_eb_sdc_spatial_diffusion_source), optional :: evaluate_spatial_diffusion
    integer, optional, intent(in) :: spatial_diffusion_halo_width
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),iterations
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), intent(in) :: dt,force(3),rtol,atol
    logical, intent(in) :: chemistry,transport_flags(4)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),minimum_theta
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: fractions_start(:),fractions_end(:)
    type(amr_tree_state_3d), optional, intent(in) :: source_start,source_end
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: integrated_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: old_predictor,new_predictor
    type(amr_tree_state_3d), optional, intent(inout) :: redistribution_transfer
    type(amr_tree_state_3d) :: nonface0,nonface1
    type(amr_tree_state_3d) :: transfer0,transfer1,hydro_transfer,final_transfer
    ! Accepted FCT means are authoritative. Re-encoding them as a new stage
    ! and averaging again can erase a resolved trace under a large old RHS.
    type(amr_tree_state_3d) :: admitted_mean_rate,admitted_mean_transfer
    type(eb_integrated_flux_3d), allocatable :: admitted_mean_faces(:),admitted_held_faces(:)
    type(amr_tree_state_3d) :: admitted_held_rate
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(in) :: target_volume_fraction
    procedure(mpi_eb_sdc_auxiliary_source), optional :: auxiliary_source
    real(dp), optional, intent(inout) :: auxiliary_change(:)
    procedure(mpi_eb_sdc_predictor_builder), optional :: predictor_builder
    type(amr_tree_state_3d), allocatable :: predictor_halo
    type(amr_tree_state_3d) :: a0,a1
    real(dp), allocatable :: auxiliary(:),auxiliary_inventory0(:),auxiliary_inventory1(:)
    real(dp), optional, intent(in) :: gas_acceleration(3)
    real(dp) :: accel(3)
    integer :: redist_order,redist_low,redist_high
    real(dp) :: redist_target
    character(len=16) :: redistribution_method
    type(amr_tree_state_3d), allocatable :: redist_halos
    type(eb_geometry_3d), allocatable :: redist_geometry(:)
    type(amr_tree_state_3d) :: redistributed,before_sync,sync_rate
    type(amr_tree_state_3d) :: op,np
    real(dp), allocatable :: midpoint(:),fraction_low(:),fraction_high(:)
    integer :: selected_depth,selection(6),minselection(6),maxselection(6),levels
    type(MPI_Comm), intent(in) :: comm
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: state_redist_registers(:)
    type(state_redist_register_3d), allocatable :: stage_registers0(:),stage_registers1(:),final_registers(:)
    logical :: prepared_redist,prepared_godunov
    logical :: capture_prepared
    integer :: godunov_padding
    type(state_redist_register_3d), allocatable :: hydro_registers(:),iteration_registers(:),prior_hydro_registers(:)
    type(eb_geometry_3d), allocatable :: register_geometry(:)
    integer :: prepared_flag,prepared_min,prepared_max
    integer, optional, intent(in) :: state_redist_max_order
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: wall_temperature(6),embedded_wall_temperature
    type(eb_geometry_3d), allocatable :: heat_geometry(:),mol_geometry(:)
    logical, optional, intent(in) :: embedded_noslip,mol_mode,mol_lax_friedrichs
    logical :: noslip,wall_active,mol
    logical :: mol_source_funded,source_funded
    real(dp), allocatable :: heat_stage,heat_frequency,wall_impulse(:)
    real(dp) :: impulse_old(3),impulse_new(3)
    real(dp) :: heat_old,heat_new,mol_wall_old(size(wall_change)),mol_wall_new(size(wall_change))
    type(amr_tree_state_3d) :: old,ot,candidate,ct,e0,e1,d0,d1,rates,hydro,ir,ledger
    type(amr_tree_state_3d) :: held_old,held_candidate,held_rate0,held_rate1,auxiliary_packed0,auxiliary_packed1
    type(eb_integrated_flux_3d), allocatable :: f0(:),f1(:),hf(:),baseline_hf(:)
    type(eb_integrated_flux_3d), allocatable :: held_flux0(:),held_flux1(:),held_flux(:)
    type(state_redist_register_3d), allocatable :: held_registers0(:),held_registers1(:),held_registers(:)
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: external(size(source_change)),v0(size(source_change)),v1(size(source_change)),theta,stage_theta
    real(dp) :: controls(16),low(16),high(16),recovered,speed
    real(dp) :: frame_values(4),frame_min(4),frame_max(4)
    real(dp), allocatable :: primitive(:),wall_rate(:,:,:,:)
    integer :: spatial_callback_width
    integer :: counts(30),mincounts(30),maxcounts(30),nv,rank,status,node,n(3),it,i,j,k,c,ios,sync_depth
    integer :: supplied,minimum_supplied,maximum_supplied,axis_min,axis_max
    logical :: valid,ready,diffusive,has_diffusion
    logical :: unsafe_forward,correct_forward,correction_periodic(3)
    real(dp) :: forward_state(reactive_nvar(size(species)))
    integer :: heat_supplied,heat_min,heat_max,heat_status
    integer :: held_supplied,held_min,held_max
    if(present(provider_refused)) provider_refused=.false.
    admission_cache_count=0;admission_cache_next=1
    admission_state=reaction_admission_not_attempted
    if(present(reaction_admission_status)) then
      reaction_admission_status=admission_state
      if(present(evaluate_reaction_admission)) reaction_admission_status=reaction_admission_unknown
    endif
    ok=.false.
    call MPI_Allreduce(merge(1,0,present(evaluate_reaction_admission)),admission_low,1,MPI_INTEGER,MPI_MIN,comm,admission_ierr)
    if(admission_ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,present(evaluate_reaction_admission)),admission_high,1,MPI_INTEGER,MPI_MAX,comm,admission_ierr)
    if(admission_ierr/=MPI_SUCCESS.or.admission_low/=admission_high) then
      if(present(reaction_admission_status)) reaction_admission_status=reaction_admission_unknown
      return
    endif
    call agree(.not.present(evaluate_reaction_admission).or. &
      (chemistry.and.present(evaluate_reaction)),comm,ready)
    if(.not.ready) then
      if(present(reaction_admission_status)) reaction_admission_status=reaction_admission_unknown
      return
    endif
    godunov_partial=.false.
    if(present(godunov_replay_partial)) godunov_partial=godunov_replay_partial
    godunov_stage_capture=present(godunov_interval_history).or.present(godunov_source_recovery)
    godunov_recovery_agreed=.false.;godunov_recovery_forbidden=.false.;godunov_source_attempt_incomplete=.false.
    godunov_packet_pending=.false.
    if(present(godunov_source_recovery_ready)) godunov_source_recovery_ready=.false.
    native_prefix_agreed=.false.;native_prefix_observed=0
    native_source_capture=present(interval_history).or.present(native_source_prefix)
    if(present(native_source_prefix_ready)) native_source_prefix_ready=.false.
    ok=.false.
    call run_common()
    call MPI_Allreduce(admission_state,admission_global,1,MPI_INTEGER,MPI_MAX,comm,admission_ierr)
    if(admission_ierr/=MPI_SUCCESS) then
      ok=.false.;admission_global=reaction_admission_unknown
    endif
    admission_state=admission_global
    if(.not.ok.and.present(evaluate_reaction_admission).and.admission_state==reaction_admission_not_attempted) &
      admission_state=reaction_admission_unknown
    if(ok.and.admission_state/=reaction_admission_not_attempted.and. &
      admission_state/=reaction_admission_provider_refused) admission_state=reaction_admission_accepted
    if(present(reaction_admission_status)) reaction_admission_status=admission_state
    if(admission_state==reaction_admission_provider_refused) then
      ok=.false.
      if(present(provider_refused)) provider_refused=.true.
    endif

    if(.not.ok.and.godunov_recovery_agreed.and.admission_state/=reaction_admission_provider_refused) &
      call publish_source_prefix()
    if(.not.ok.and.native_prefix_agreed.and.admission_state/=reaction_admission_provider_refused) &
      call publish_native_source_prefix()
  contains
    subroutine run_common()
    held_auxiliary=.false.
    if(present(auxiliary_has_internal_energy)) held_auxiliary=auxiliary_has_internal_energy
    external_diffusion_source=.false.
    if(present(spatial_diffusion_is_source)) external_diffusion_source=spatial_diffusion_is_source
    ok=.false.;nv=reactive_nvar(size(species))
    naux=eb_auxiliary_nvar(nv,auxiliary_has_internal_energy);theta=1
    allocate(auxiliary(naux) ,auxiliary_inventory0(naux),auxiliary_inventory1(naux),stat=auxiliary_alloc_status)
    call agree(auxiliary_alloc_status==0,comm,ready);if(.not.ready) return
    initial_reset=0;predictor_reset_change=0;final_reset=0;post_reset=0;accepted_predictor_reset=0
    initial_scalar_reset=0
    mol_source_funded=.false.;source_funded=.false.
    held_policy=dual_energy_policy()
    if(present(energy_policy)) held_policy=energy_policy
    call MPI_Allreduce(merge(1,0,held_auxiliary),held_aux_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,held_auxiliary),held_aux_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.held_aux_low/=held_aux_high) return
    call agree(.not.held_auxiliary.or.(present(internal_energy).and.present(auxiliary_source)),comm,ready)
    if(.not.ready) return
    if(present(provider_refused)) provider_refused=.false.
    heat_supplied=merge(1,0,present(heat_release))
    call MPI_Allreduce(heat_supplied,heat_min,1,MPI_INTEGER,MPI_MIN,comm,heat_status)
    if(heat_status/=MPI_SUCCESS) return
    call MPI_Allreduce(heat_supplied,heat_max,1,MPI_INTEGER,MPI_MAX,comm,heat_status)
    if(heat_status/=MPI_SUCCESS.or.heat_min/=heat_max) return
    held_supplied=merge(1,0,present(internal_energy))+2*merge(1,0,present(internal_energy_flux))+ &
      4*merge(1,0,present(ancestor_internal_energy_start))+8*merge(1,0,present(ancestor_internal_energy_end))+ &
      16*merge(1,0,present(energy_correction))+32*merge(1,0,present(energy_policy))+ &
      64*merge(1,0,present(internal_energy_pressure_work_change))+ &
      128*merge(1,0,present(internal_energy_projection_change))+ &
      256*merge(1,0,present(internal_energy_closure_change))
    call MPI_Allreduce(held_supplied,held_min,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(held_supplied,held_max,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.held_min/=held_max) return
    if(iand(held_supplied,15)/=0.and.iand(held_supplied,15)/=3.and.iand(held_supplied,15)/=15) return
    valid=(.not.present(energy_policy).or.present(energy_correction)).and. &
      (.not.present(energy_correction).or.present(internal_energy))
    if(present(energy_correction)) valid=valid.and.ieee_is_finite(energy_correction).and.ieee_is_finite(held_policy%eta2)
    call agree(valid,comm,ready);if(.not.ready) return
    if(present(energy_correction)) then
      held_controls=[energy_correction,held_policy%eta2]
      held_policy_flags=[merge(1,0,held_policy%allow_small_energy), &
        merge(1,0,held_policy%allow_negative_energy),merge(1,0,held_policy%update_total_from_internal)]
      call MPI_Allreduce(held_controls,held_controls_low,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(held_controls,held_controls_high,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(held_controls_low/=held_controls_high)) return
      call MPI_Allreduce(held_policy_flags,held_policy_min,3,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(held_policy_flags,held_policy_max,3,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(held_policy_min/=held_policy_max)) return
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
    endif
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,constant_agreed,species=species)
    if(.not.constant_agreed) return
    spatial_callback_width=6
    if(present(spatial_diffusion_halo_width)) spatial_callback_width=spatial_diffusion_halo_width
    prepared_redist=present(state_redist_registers);prepared_flag=merge(1,0,prepared_redist)
    call MPI_Allreduce(prepared_flag,prepared_min,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(prepared_flag,prepared_max,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.prepared_min/=prepared_max) return
    redist_order=0
    if(present(state_redist_max_order)) redist_order=state_redist_max_order
    call agree(redist_order==0.or.redist_order==2.or.redist_order==4,comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(redist_order,redist_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(redist_order,redist_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.redist_low/=redist_high) return
    supplied=merge(1,0,present(frame))
    call MPI_Allreduce(supplied,minimum_supplied,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(supplied,maximum_supplied,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.minimum_supplied/=maximum_supplied) return
    if(present(frame)) then
      call agree(frame%is_valid(),comm,ready)
      if(.not.ready) return
      call MPI_Allreduce(frame%axis,axis_min,1,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame%axis,axis_max,1,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS) return
      frame_values=[frame%angular_speed,frame%axis_point]
      call MPI_Allreduce(frame_values,frame_min,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame_values,frame_max,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.axis_min/=axis_max.or.any(frame_min/=frame_max)) return
    end if
    has_diffusion=present(evaluate_diffusion)
    accel=0
    if(present(gas_acceleration)) accel=gas_acceleration
    redistribution_method='none'
    if(present(redistribution)) redistribution_method=redistribution
    redist_target=.5_dp
    if(present(target_volume_fraction)) redist_target=target_volume_fraction
    call agree(redistribution_method=='none'.or.redistribution_method=='flux_redist'.or. &
      redistribution_method=='state_redist',comm,ready)
    if(.not.ready) return
    controls(:11)=[dt,force,rtol,atol,merge(1._dp,0._dp,redistribution_method=='flux_redist')+ &
      merge(2._dp,0._dp,redistribution_method=='state_redist'),redist_target,accel]
    mol=.false.
    if(present(mol_mode)) mol=mol_mode
    selected_momentum=.false.;selected_velocity=.false.;selected_pressure=.false.;selected_connected=.false.
    if(present(separate_momentum)) selected_momentum=separate_momentum
    if(present(consistent_velocity)) selected_velocity=consistent_velocity
    if(present(paired_pressure)) selected_pressure=paired_pressure
    if(present(connected_directions)) selected_connected=connected_directions
    call MPI_Allreduce(merge(1,0,selected_momentum)+merge(2,0,selected_velocity)+merge(4,0, &
      selected_pressure)+merge(8,0,selected_connected), &
      momentum_min,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,selected_momentum)+merge(2,0,selected_velocity)+merge(4,0, &
      selected_pressure)+merge(8,0,selected_connected), &
      momentum_max,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.momentum_min/=momentum_max) return
    call agree((.not.selected_pressure.or.selected_velocity).and. &
      (.not.selected_pressure.or..not.present(frame)),comm,ready)
    if(.not.ready) return
    call agree((.not.selected_velocity.or.selected_momentum).and. &
      (.not.selected_connected.or.selected_momentum),comm,ready)
    if(.not.ready) return
    call agree(.not.selected_momentum.or.(mol.and.present(internal_energy).and. &
      redistribution_method=='none'.and..not.present(state_redist_registers)),comm,ready)
    if(.not.ready) return
    controls(15)=merge(1._dp,0._dp,mol)
    controls(16)=0
    if(present(mol_lax_friedrichs)) controls(16)=merge(1._dp,0._dp,mol_lax_friedrichs)
    mol_wall_old=0;mol_wall_new=0
    noslip=.false.
    if(present(embedded_noslip)) noslip=embedded_noslip
    wall_active=noslip.or.present(embedded_wall_temperature)
    controls(14)=merge(1._dp,0._dp,noslip)
    controls(12:13)=0
    if(present(embedded_wall_temperature)) controls(12:13)=[1._dp,embedded_wall_temperature]
    heat_old=0;heat_new=0;impulse_old=0;impulse_new=0
    valid=all(ieee_is_finite(controls)).and.min(dt,rtol,atol)>0.and.iterations>=1
    valid=valid.and.(controls(16)==0.or.mol)
    valid=valid.and.redist_target>0.and.redist_target<=1.and.all(abs(accel)<=1e6_dp)
    valid=valid.and.all([size(boundary_change),size(wall_change),size(chemical_change),size(source_change)]==nv)
    valid=valid.and..not.(transport_flags(4).and..not.transport_flags(3))
    valid=valid.and.(present(auxiliary_source).eqv.present(auxiliary_change))
    if(present(auxiliary_change)) valid=valid.and.size(auxiliary_change)==naux
    valid=valid.and.(chemistry.or..not.present(evaluate_reaction))
    valid=valid.and.(.not.external_diffusion_source.or.present(evaluate_spatial_diffusion).or.has_diffusion)
    if(present(embedded_wall_temperature)) valid=valid.and.transport_flags(2).and.embedded_wall_temperature>0
    if(noslip) valid=valid.and.transport_flags(1)
    ! Particle predictor callbacks still require their own MOL staging.
    if(mol) valid=valid.and..not.present(predictor_builder)
    if(present(frame)) valid=valid.and.mol
    if(present(accepted_forcing).or.present(accepted_mass_source).or.present(old_spatial_flux)) &
      valid=valid.and.mol
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(controls,low,16,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,16,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    native_external_frozen=.false.
    if(present(freeze_native_external_sources)) native_external_frozen=freeze_native_external_sources
    counts=[iterations,merge(1,0,chemistry),merge(1,0,transport_flags),merge(1,0,present(reaction_source)), &
      merge(1,0,present(auxiliary_source)),merge(1,0,present(predictor_builder)),merge(1,0,present(evaluate_diffusion)), &
      merge(1,0,present(evaluate_reaction)),merge(1,0,present(evaluate_spatial_diffusion)), &
      spatial_callback_width,merge(1,0,present(interval_history)),merge(1,0,present(replay_interval)), &
      merge(1,0,present(replay_spatial_weight)),merge(1,0,present(provider_refused)), &
      merge(1,0,external_diffusion_source),merge(1,0,present(spatial_flux_weight)),merge(1,0,present(godunov_interval_history)), &
      merge(1,0,present(godunov_replay_interval)),merge(1,0,godunov_partial), &
      merge(1,0,present(godunov_source_recovery)),merge(1,0,present(godunov_source_recovery_ready)), &
      merge(1,0,present(freeze_native_external_sources)),merge(1,0,native_external_frozen), &
      merge(1,0,present(native_spatial_weight)), &
      merge(1,0,present(native_source_prefix)),merge(1,0,present(native_source_prefix_ready)), &
      merge(1,0,present(native_source_replay))]
    call MPI_Allreduce(counts,mincounts,30,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(counts,maxcounts,30,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(mincounts/=maxcounts)) return
    call agree(.not.native_external_frozen.or.(mol.and.chemistry.and.present(replay_interval).and. &
      present(internal_energy).and.present(energy_correction)),comm,ready)
    if(.not.ready) return
    native_mol_weight=1
    if(present(native_spatial_weight)) native_mol_weight=native_spatial_weight
    call agree(ieee_is_finite(native_mol_weight),comm,ready);if(.not.ready) return
    call agree(native_mol_weight>=0.and.native_mol_weight<=1.and. &
      (.not.present(native_spatial_weight).or.native_external_frozen.or.present(native_source_replay)),comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(native_mol_weight,native_mol_weight_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(native_mol_weight,native_mol_weight_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.native_mol_weight_min/=native_mol_weight_max) return
    call agree(present(native_source_prefix).eqv.present(native_source_prefix_ready),comm,ready)
    if(.not.ready) return
    native_prefix_agreed=present(native_source_prefix)
    native_prefix_scope=native_prefix_agreed.or.present(native_source_replay)
    if(native_prefix_scope) then
      valid=mol.and.chemistry.and.present(target_depth).and.redistribution_method=='none'.and. &
        present(internal_energy).and.present(energy_correction).and. &
        .not.present(evaluate_spatial_diffusion).and..not.present(evaluate_diffusion).and. &
        .not.present(predictor_builder).and..not.present(godunov_source_recovery)
      call agree(valid,comm,ready);if(.not.ready) return
    endif
    call agree(present(godunov_source_recovery).eqv.present(godunov_source_recovery_ready),comm,ready)
    if(.not.ready) return
    call agree(.not.godunov_partial.or.present(godunov_replay_interval),comm,ready)
    if(.not.ready) return
    godunov_recovery_agreed=present(godunov_source_recovery)
    native_weight=1
    if(present(spatial_flux_weight)) native_weight=spatial_flux_weight
    call agree(ieee_is_finite(native_weight),comm,ready)
    if(.not.ready) return
    call agree(native_weight>=0.and.native_weight<=1.and. &
      (.not.present(spatial_flux_weight).or.(.not.mol.and..not.present(replay_spatial_weight).and. &
      .not.present(replay_interval))),comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(native_weight,native_weight_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(native_weight,native_weight_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.native_weight_min/=native_weight_max) return
    spatial_weight=1
    if(present(replay_spatial_weight)) then
      spatial_weight=replay_spatial_weight
      call agree(ieee_is_finite(spatial_weight).and.spatial_weight>=0.and.spatial_weight<=1.and. &
        present(replay_interval),comm,ready)
      if(.not.ready) return
      call MPI_Allreduce(spatial_weight,spatial_weight_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(spatial_weight,spatial_weight_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.spatial_weight_min/=spatial_weight_max) return
    endif
    selected_depth=-1
    if(present(target_depth)) selected_depth=target_depth
    if(native_prefix_scope) then
      call agree(selected_depth>=0,comm,ready);if(.not.ready) return
    endif
    funding_flags=merge(1,0,present(source_funding_stock))+2*merge(1,0,present(source_funding_impulse))
    call MPI_Allreduce(funding_flags,funding_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(funding_flags,funding_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.funding_low/=funding_high) return
    call agree(present(source_funding_stock).eqv.present(source_funding_impulse),comm,ready)
    if(.not.ready) return
    native_held_funding=selected_depth<0.or.funding_flags==3
    if(funding_flags==3) then
      call agree(present(internal_energy),comm,ready)
      if(.not.ready) return
      call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,source_funding_stock,nv+1,comm,ready)
      if(.not.ready) return
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        if(selected_depth<0.or.tree%nodes(node)%level==selected_depth) &
          valid=valid.and.all(source_funding_stock%fields(node)%values==0)
        ! Pending loans belong to exposed fluid inventory. Virtual covered
        ! cells and solid EB cells cannot own a debt that admission later
        ! clears without settlement.
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(covered(tree,node,[i,j,k]).or.geometry(node)%volume_fraction(i,j,k)<=0) then
            valid=valid.and.all(source_funding_stock%fields(node)%values(:,i,j,k)==0)
          endif
        enddo
        enddo
        enddo
      enddo
      call agree(valid,comm,ready);if(.not.ready) return
    endif
    selection=[selected_depth,sum(merge(1,0,[present(ancestor_start),present(ancestor_end), &
      present(temperature_start),present(temperature_end),present(fractions_start),present(fractions_end)])), &
      sum(merge(1,0,[present(source_start),present(source_end)])),merge(1,0,present(integrated_flux)), &
      2*merge(1,0,present(old_predictor))+merge(1,0,present(new_predictor)), &
      merge(1,0,present(redistribution_transfer))]
    if(present(interval_history)) then
      call agree(mol.and.chemistry.and.selected_depth>=0,comm,ready)
      if(.not.ready) return
    endif
    valid=selected_depth>=-1.and.spatial_callback_width>=6.and.spatial_callback_width<=huge(1)/4
    if(spatial_callback_width>6) valid=valid.and.present(evaluate_spatial_diffusion)
    valid=valid.and..not.(has_diffusion.and.present(evaluate_spatial_diffusion))
    if(present(evaluate_spatial_diffusion)) then
      valid=valid.and.(redistribution_method=='state_redist'.or.(mol.and.redistribution_method=='none'))
      if(selected_depth>=0.and.redistribution_method=='state_redist') &
        valid=valid.and.present(state_redist_registers)
    endif
    if(present(target_depth)) valid=valid.and.selected_depth>=0.and.present(integrated_flux)
    valid=valid.and.(selection(2)==0.or.selection(2)==6).and.(selection(3)==0.or.selection(3)==2)
    if(selected_depth<0.and.mol.and.redistribution_method=='state_redist'.and..not.has_diffusion) &
      prepared_redist=.true.
    if(selected_depth<0.and..not.mol.and.redistribution_method=='state_redist'.and..not.has_diffusion) then
      if(size(extended_geometry)==size(tree%nodes)) then
        prepared_godunov=.true.
        do node=1,size(tree%nodes)
          n=tree%nodes(node)%shape
          if(.not.(all([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]==n+18).or. &
            all([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]==n+20))) &
            prepared_godunov=.false.
        enddo
        prepared_redist=prepared_redist.or.prepared_godunov
      endif
    endif
    prepared_godunov=prepared_redist.and..not.mol
    if(present(source_funding_stock).and.prepared_redist) return
    godunov_padding=9
    if(prepared_godunov.and.size(extended_geometry)>0) then
      godunov_padding=(extended_geometry(1)%nx-tree%nodes(1)%shape(1))/2
    endif
    if(prepared_redist) then
      valid=valid.and.redistribution_method=='state_redist'
      valid=valid.and..not.present(redistribution_transfer).and..not.has_diffusion
    else if(selected_depth>=0.and.redistribution_method/='none') then
      valid=valid.and.present(redistribution_transfer)
    else
      valid=valid.and..not.present(redistribution_transfer)
    end if
    if(present(evaluate_spatial_diffusion)) valid=valid.and.(prepared_redist.or.(mol.and.redistribution_method=='none'))
    if(selection(2)>0) valid=valid.and.selected_depth>=0
    if(selection(3)>0) valid=valid.and.selection(2)==6
    valid=valid.and.(.not.present(internal_energy_pressure_work_change).or.present(internal_energy))
    valid=valid.and.(.not.present(internal_energy_projection_change).or.present(internal_energy))
    valid=valid.and.(.not.present(internal_energy_closure_change).or.present(internal_energy))
    if(present(spatial_flux_weight)) valid=valid.and. &
      (redistribution_method=='none'.or.prepared_godunov)
    if(godunov_stage_capture) valid=valid.and..not.mol.and.present(internal_energy).and. &
      (redistribution_method=='none'.or.prepared_godunov).and..not.present(replay_interval).and. &
      .not.present(interval_history).and.(.not.present(ancestor_start).or. &
      (present(ancestor_internal_energy_start).and.present(ancestor_internal_energy_end)))
    if(present(internal_energy)) valid=valid.and. &
      (.not.chemistry.or.((redistribution_method=='none'.or.prepared_redist).and.present(energy_correction))).and. &
      (.not.has_diffusion.or.redistribution_method=='none').and. &
      (.not.present(replay_interval).or.prepared_redist.or.(mol.and.redistribution_method=='none')).and. &
      (mol.or.redistribution_method=='none'.or.prepared_godunov).and. &
      (redistribution_method=='none'.or.redistribution_method=='state_redist').and. &
      .not.(selected_depth>=0.and.redistribution_method=='state_redist'.and. &
      .not.present(held_state_redist_registers)).and. &
      (.not.present(auxiliary_source).or.chemistry.or.held_auxiliary)
    if(present(held_state_redist_registers)) valid=valid.and.present(internal_energy).and.selected_depth>=0.and.prepared_redist
    if(present(ancestor_internal_energy_start)) valid=valid.and.selected_depth>=0.and.selection(2)==6
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(selection,minselection,6,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(selection,maxselection,6,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(minselection/=maxselection)) return
    if(present(godunov_replay_interval)) then
      call agree(.not.mol.and.present(internal_energy).and..not.present(interval_history).and. &
        .not.present(replay_interval).and.(redistribution_method=='none'.or.prepared_godunov),comm,ready)
      if(.not.ready) return
      call validate_godunov_replay(ready)
      if(.not.ready) return
    endif
    old=state;ot=temperature
    if(present(internal_energy)) held_old=internal_energy
    if(present(energy_correction)) then
      call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,old,nv,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,ot,1,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,ready)
      if(.not.ready) return
      held_old=internal_energy
      call reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,old,ot,held_old,reset_delta, &
        held_policy,selected_depth,initial_reset,comm,ready,frame, &
        internal_energy_change=initial_scalar_reset,reset_level=merge(-2,-1,mol.and.chemistry.and.redistribution_method=='none'))
      if(.not.ready) return
    endif
    if(selected_depth<0.and..not.(mol.and.chemistry.and.redistribution_method=='none'.and.present(energy_correction))) then
      if(present(internal_energy)) then
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,old,ot,comm,ready,frame,held_old)
      else
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,old,ot,comm,ready,frame)
      endif
      if(.not.ready) return
    else
      call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,old,nv,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,ot,1,comm,ready)
      if(.not.ready) return
      if(selected_depth>=tree%level_count()) return
    end if
    if(present(internal_energy)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,ready)
      if(.not.ready) return
      if(present(ancestor_internal_energy_start)) then
        call validate_mpi_amr_tree_local_3d(tree,owners,ancestor_internal_energy_start,1,comm,ready)
        if(.not.ready) return
        call validate_mpi_amr_tree_local_3d(tree,owners,ancestor_internal_energy_end,1,comm,ready)
        if(.not.ready) return
      endif
      held_candidate=held_old
    endif
    if(prepared_redist) then
      call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,ready, &
        padding=merge(godunov_padding,6,prepared_godunov))
      if(.not.ready) return
      allocate(register_geometry(size(owners)),stat=ios)
      call agree(ios==0,comm,ready);if(.not.ready) return
      valid=.true.
      do node=1,size(owners)
        n=tree%nodes(node)%shape
        if(prepared_godunov) then
          call extract_eb_geometry_3d(extended_geometry(node),[1,1,1]*(godunov_padding-5), &
            n+godunov_padding+6,register_geometry(node),valid)
          if(.not.valid) exit
        else
          register_geometry(node)=extended_geometry(node)
        endif
      enddo
      call agree(valid,comm,ready);if(.not.ready) return
    endif
    if(mol) then
      call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,ready,padding=6)
      if(.not.ready) return
      allocate(mol_geometry(size(owners)),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      do node=1,size(owners)
        n=tree%nodes(node)%shape
        call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],n+9,mol_geometry(node),valid)
        if(.not.valid) exit
      end do
      call agree(valid,comm,ready)
      if(.not.ready) return
    end if
    if(wall_active) then
      call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,ready, &
        padding=merge(godunov_padding,6,prepared_godunov))
      if(.not.ready) return
      if(noslip) then
        allocate(wall_impulse(3),stat=ios)
        call agree(ios==0,comm,ready)
        if(.not.ready) return
      end if
      allocate(heat_geometry(size(owners)),heat_stage,heat_frequency,stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      do node=1,size(owners)
        n=tree%nodes(node)%shape
        if(prepared_godunov) then
          call extract_eb_geometry_3d(register_geometry(node),[5,5,5],n+8,heat_geometry(node),valid)
        else
          call extract_eb_geometry_3d(extended_geometry(node),[5,5,5],n+8,heat_geometry(node),valid)
        endif
        if(.not.valid) exit
      end do
      call agree(valid,comm,ready)
      if(.not.ready) return
    end if
    if(selection(2)==6) then
      levels=tree%level_count()-1
      call agree(size(fractions_start)==levels.and.size(fractions_end)==levels,comm,ready)
      if(.not.ready) return
      call agree(all(ieee_is_finite(fractions_start)).and.all(ieee_is_finite(fractions_end)),comm,ready)
      if(.not.ready) return
      call agree(all(fractions_start>=0).and.all(fractions_end<=1).and.all(fractions_end>=fractions_start),comm,ready)
      if(.not.ready) return
      allocate(fraction_low(2*levels),fraction_high(2*levels),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      call MPI_Allreduce([fractions_start,fractions_end],fraction_low,2*levels,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce([fractions_start,fractions_end],fraction_high,2*levels,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(fraction_low/=fraction_high)) return
      if(selection(3)==2) midpoint=.5_dp*(fractions_start+fractions_end)
    end if
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    allocate(ir%fields(size(owners)),primitive(reactive_nprim(size(species))),stat=ios)
    call agree(ios==0,comm,ready)
    if(.not.ready) return
    if(present(heat_release)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,heat_release,1,comm,valid)
      if(.not.valid) return
      allocate(heat);heat=heat_release
    endif
    if(present(heat_release).and..not.chemistry) then
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(selected_depth>=0.and.tree%nodes(node)%level/=selected_depth) cycle
        heat%fields(node)%values=0
      enddo
    endif
    if(present(reaction_source)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,reaction_source,nv-4,comm,ready)
      if(.not.ready) return
      ir=reaction_source
    else
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)),source=0._dp,stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        end if
      end do
      call agree(valid,comm,ready)
      if(.not.ready) return
    end if
    do node=1,size(owners)
      if(owners(node)/=rank.or..not.active(node)) cycle
      ! The pinned disabled-chemistry path resets I_R before predicting.
      if(.not.chemistry) ir%fields(node)%values=0
      do c=1,nv-4
        where(geometry(node)%volume_fraction<=0) ir%fields(node)%values(c,:,:,:)=0
      end do
      if(.not.all(ieee_is_finite(ir%fields(node)%values))) valid=.false.
    end do
    call agree(valid,comm,ready)
    if(.not.ready) return
    if(present(replay_interval)) then
      if(prepared_redist) then
        call validate_prepared_replay(ready)
      else
        call validate_native_mol_replay(ready)
      endif
      if(.not.ready) return
      do node=1,size(owners)
        if(owners(node)/=rank.or..not.active(node)) cycle
        ir%fields(node)%values=replay_interval%initial_reaction_source%fields(node)%values
      enddo
    endif
    if(present(native_source_replay)) then
      call validate_native_source_replay(ready)
      if(.not.ready) return
    endif
    capture_prepared=prepared_redist.and.mol.and.(chemistry.or.present(interval_history))
    nullify(prepared_integrator)
    if(present(evaluate_reaction)) prepared_integrator=>prepared_reactor_adapter
    reactor_refused=.false.;reactor_refused_local=.false.
    if(present(provider_refused)) provider_refused=.false.
    if(capture_prepared.or.native_source_capture) then
      call copy_field(ir,pending_history%initial_reaction_source,ready)
      if(.not.ready) return
    endif
    if(native_prefix_agreed) then
      ! Prefix reuse binds the whole owned entry, while an accepted interval
      ! retains the existing selected-level reaction-source contract.
      call copy_field(ir,native_prefix_initial_reaction_source,ready,all_levels=.true.)
      if(.not.ready) return
    endif
    if(native_source_capture.and.mol.and.chemistry.and.present(internal_energy).and. &
      selected_depth>=0.and.redistribution_method=='none') then
      pending_history%separate_momentum=selected_momentum
      pending_history%consistent_velocity=selected_velocity
      pending_history%paired_pressure=selected_pressure
      pending_history%connected_directions=selected_connected
      pending_history%native_iterations=iterations;pending_history%native_force=force
      allocate(pending_history%external_stages(0:iterations),stat=ios)
      call agree(ios==0,comm,ready);if(.not.ready) return
      pending_history%native_acceleration=accel;pending_history%has_native_frame=present(frame)
      pending_history%native_providers=[present(evaluate_spatial_diffusion),present(evaluate_diffusion)]
      pending_history%native_provider_halo_width=spatial_callback_width
      if(external_diffusion_source.and.any(pending_history%native_providers)) then
        allocate(pending_history%provider_stages(2,0:iterations),stat=ios)
        call agree(ios==0,comm,ready);if(.not.ready) return
        allocate(pending_history%native_provider_geometry,source=geometry,stat=ios)
        call agree(ios==0,comm,ready);if(.not.ready) return
      endif
      if(present(frame)) pending_history%native_frame=frame
      if(present(auxiliary_source)) then
        allocate(pending_history%auxiliary_stages(0:iterations),stat=ios)
        call agree(ios==0,comm,ready);if(.not.ready) return
      endif
    endif
    if(godunov_stage_capture) then
      call begin_godunov_history(ready)
      if(.not.ready) return
    endif
    candidate=old;ct=ot;diffusive=any(transport_flags(:3))
    call sources(old,ot,e0,d0,a0,0,f0,transfer0,ready,fractions_start)
    if(.not.ready) then
      godunov_source_attempt_incomplete=.true.;return
    endif
    if(native_prefix_agreed) native_prefix_observed=max(native_prefix_observed,1)
    if(godunov_stage_capture) then
      call capture_godunov_source(0,old,ot,e0,d0,a0,f0,ready)
      if(.not.ready) then
        godunov_source_attempt_incomplete=.true.;return
      endif
    endif
    e1=e0;d1=d0;heat_new=heat_old;impulse_new=impulse_old;mol_wall_new=mol_wall_old
    if(prepared_redist) stage_registers1=stage_registers0
    if(present(auxiliary_source)) a1=a0
    call predictor_field(e0,d0,a0,op,ready)
    if(.not.ready) return
    if(mol.and.present(internal_energy)) then
      held_hydro=held_old;held_hydro_flux=held_flux0
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        held_hydro%fields(node)%values=0
        if(allocated(held_hydro_flux(node)%x)) held_hydro_flux(node)%x=0
        if(allocated(held_hydro_flux(node)%y)) held_hydro_flux(node)%y=0
        if(allocated(held_hydro_flux(node)%z)) held_hydro_flux(node)%z=0
      enddo
    endif
    if(mol) then
      hydro=d0;hf=f0
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        hydro%fields(node)%values=0
        if(.not.active(node)) cycle
        if(allocated(hf(node)%x)) hf(node)%x=0
        if(allocated(hf(node)%y)) hf(node)%y=0
        if(allocated(hf(node)%z)) hf(node)%z=0
      end do
      if(present(redistribution_transfer)) hydro_transfer=hydro
    end if
    do it=1,merge(1,iterations,mol.and..not.chemistry)
      if(mol.and.chemistry.and.it>1.and.mol_source_funded) then
        ! A material loan belongs to its actual native mean. Later reacting
        ! MOL passes rebuild the new spatial RHS and omit its prior new-time
        ! volume source. Retain that provisional state as the new RHS input,
        ! but replace the previous numerical loan/face family from scratch.
        ! In common, nonprepared MOL these hydro fields contain only funding.
        do node=1,size(owners)
        if(owners(node)/=rank) cycle
          hydro%fields(node)%values=0;held_hydro%fields(node)%values=0
          if(allocated(hf(node)%x)) hf(node)%x=0
        if(allocated(hf(node)%y)) hf(node)%y=0
        if(allocated(hf(node)%z)) hf(node)%z=0
          if(allocated(held_hydro_flux(node)%x)) held_hydro_flux(node)%x=0
        if(allocated(held_hydro_flux(node)%y)) held_hydro_flux(node)%y=0
        if(allocated(held_hydro_flux(node)%z)) held_hydro_flux(node)%z=0
        enddo
        mol_source_funded=.false.;source_funded=.false.
      endif
      ! Only the accepted SDC iterate contributes its predictor reset to
      ! the published forcing ledger. Earlier iterates are provisional.
      accepted_predictor_reset=0
      if(allocated(admitted_mean_rate%fields)) deallocate(admitted_mean_rate%fields)
      if(allocated(admitted_mean_transfer%fields)) deallocate(admitted_mean_transfer%fields)
      if(allocated(admitted_mean_faces)) deallocate(admitted_mean_faces)
      if(allocated(admitted_held_faces)) deallocate(admitted_held_faces)
      if(allocated(admitted_held_rate%fields)) deallocate(admitted_held_rate%fields)
      ! Later reacting MOL iterations re-evaluate the previous reacted state;
      ! only the first iteration forms the Euler predictor.
      if(.not.mol.or.it==1) then
      call combined_rates(.true.)
      if(present(predictor_builder)) then
        if(allocated(predictor_halo)) deallocate(predictor_halo)
        allocate(predictor_halo,stat=ios)
        call agree(ios==0,comm,ready)
        if(.not.ready) return
        call predictor_builder(rates,dt,it,selected_depth,predictor_halo,comm,valid)
        call agree(valid,comm,ready)
        if(.not.ready) return
      end if
      if(.not.mol) then
      if(prepared_godunov) then
        if(it>1) prior_hydro_registers=hydro_registers
        if(it>1.and.present(internal_energy)) prior_held_hydro_registers=held_hydro_registers
        if(present(internal_energy)) then
      call compute_godunov_state_redist_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,old,ot,rates,extended_geometry, &
          dt,hydro,hf,comm,ready,difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start, &
          source_start,source_end,midpoint,predictor_halo,registers=iteration_registers, &
          max_order=redist_order,target_volume_fraction=redist_target, &
          internal_energy_state=held_old,internal_energy_rate=held_hydro, &
          internal_energy_flux=held_hydro_flux,internal_energy_pressure_work=held_hydro_work, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end, &
          internal_energy_registers=iteration_held_hydro_registers, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_godunov_state_redist_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,old,ot,rates,extended_geometry, &
          dt,hydro,hf,comm,ready,difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start, &
          source_start,source_end,midpoint,predictor_halo,registers=iteration_registers, &
          max_order=redist_order,target_volume_fraction=redist_target, &
        spatial_flux_weight=spatial_flux_weight)
      endif
        if(.not.ready) return
        valid=.true.
        if(it>1) then
          do node=1,size(owners)
            if(owners(node)/=rank.or..not.active(node)) cycle
            n=tree%nodes(node)%shape
            call accumulate_state_redist_coarse_register_3d(iteration_registers(node),prior_hydro_registers(node),nv,n+12,ready)
            valid=valid.and.ready
            if(present(internal_energy)) then
              call accumulate_state_redist_coarse_register_3d(iteration_held_hydro_registers(node), &
                prior_held_hydro_registers(node),1,n+12,ready)
              valid=valid.and.ready
            endif
          enddo
        endif
        call agree(valid,comm,ready);if(.not.ready) return
        call move_alloc(iteration_registers,hydro_registers)
        if(present(internal_energy)) call move_alloc(iteration_held_hydro_registers,held_hydro_registers)
      else
      if(present(internal_energy)) then
      call compute_sdc_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,old,ot,rates,extended_geometry, &
        dt,hydro,hf,comm,ready,difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start, &
        source_start,source_end,midpoint,predictor_halo, &
          internal_energy_state=held_old,internal_energy_rate=held_hydro, &
          internal_energy_flux=held_hydro_flux,internal_energy_pressure_work=held_hydro_work, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_sdc_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,old,ot,rates,extended_geometry, &
        dt,hydro,hf,comm,ready,difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start, &
        source_start,source_end,midpoint,predictor_halo, &
        spatial_flux_weight=spatial_flux_weight)
      endif
      endif
      if(.not.ready) return
      ! The full forward predictor includes diffusion, auxiliary sources and
      ! the current reaction estimate. Only endangered predictors need the
      ! extra normal baseline and same-level donor exchange.
      unsafe_forward=.false.;valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank.or..not.active(node)) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          forward_state=old%fields(node)%values(:,i,j,k)+ &
            dt*(rates%fields(node)%values(:,i,j,k)+hydro%fields(node)%values(:,i,j,k))
          valid=valid.and.all(ieee_is_finite(forward_state))
          unsafe_forward=unsafe_forward.or.forward_state(1)<=0.or.any(forward_state(6:)<0)
        enddo
        enddo
        enddo
      enddo
      call agree(valid,comm,ready);if(.not.ready) return
      call MPI_Allreduce(unsafe_forward,correct_forward,1,MPI_LOGICAL,MPI_LOR,comm,status)
      if(status/=MPI_SUCCESS) return
      if(correct_forward.and.present(internal_energy).and.native_held_funding.and..not.prepared_redist) then
        call fund_held_native_stage(ready)
        if(ready) correct_forward=.false.
      endif
      if(correct_forward) then
        if(prepared_godunov) then
        if(present(internal_energy)) then
      call compute_godunov_state_redist_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,old,ot,rates,extended_geometry, &
          dt,hydro,hf,comm,ready,difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start, &
          source_start,source_end,midpoint,predictor_halo,registers=iteration_registers, &
          max_order=redist_order,target_volume_fraction=redist_target,correct_forward_flux=.true., &
          internal_energy_state=held_old,internal_energy_rate=held_hydro, &
          internal_energy_flux=held_hydro_flux,internal_energy_pressure_work=held_hydro_work, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end, &
          internal_energy_registers=iteration_held_hydro_registers, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_godunov_state_redist_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,old,ot,rates,extended_geometry, &
          dt,hydro,hf,comm,ready,difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start, &
          source_start,source_end,midpoint,predictor_halo,registers=iteration_registers, &
          max_order=redist_order,target_volume_fraction=redist_target,correct_forward_flux=.true., &
        spatial_flux_weight=spatial_flux_weight)
      endif
        if(.not.ready) return
        valid=.true.
        if(it>1) then
          do node=1,size(owners)
            if(owners(node)/=rank.or..not.active(node)) cycle
            n=tree%nodes(node)%shape
            call accumulate_state_redist_coarse_register_3d(iteration_registers(node),prior_hydro_registers(node),nv,n+12,ready)
            valid=valid.and.ready
            if(present(internal_energy)) then
              call accumulate_state_redist_coarse_register_3d(iteration_held_hydro_registers(node), &
                prior_held_hydro_registers(node),1,n+12,ready)
              valid=valid.and.ready
            endif
          enddo
        endif
        call agree(valid,comm,ready);if(.not.ready) return
        call move_alloc(iteration_registers,hydro_registers)
        if(present(internal_energy)) call move_alloc(iteration_held_hydro_registers,held_hydro_registers)
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank.or..not.active(node)) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          forward_state=old%fields(node)%values(:,i,j,k)+ &
            dt*(rates%fields(node)%values(:,i,j,k)+hydro%fields(node)%values(:,i,j,k))
          valid=valid.and.all(ieee_is_finite(forward_state))
          valid=valid.and.forward_state(1)>0.and.all(forward_state(6:)>=0)
        enddo
        enddo
        enddo
      enddo
      call agree(valid,comm,ready);if(.not.ready) return
        else
        if(present(internal_energy)) then
      call compute_sdc_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,old,ot,rates,extended_geometry, &
          dt,hydro,hf,comm,ready,difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start, &
          source_start,source_end,midpoint,predictor_halo,baseline_hf, &
          internal_energy_state=held_old,internal_energy_rate=held_hydro, &
          internal_energy_flux=held_hydro_flux,internal_energy_pressure_work=held_hydro_work, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end,internal_energy_baseline_flux=held_hydro_baseline, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_sdc_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,old,ot,rates,extended_geometry, &
          dt,hydro,hf,comm,ready,difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start, &
          source_start,source_end,midpoint,predictor_halo,baseline_hf, &
        spatial_flux_weight=spatial_flux_weight)
      endif
        if(.not.ready) return
        correction_periodic= &
          boundary_faces([1,3,5])=='periodic'.and.boundary_faces([2,4,6])=='periodic'
        if(present(internal_energy)) then
        call apply_mpi_eb_flux_correction_3d(tree,owners,geometry,old,rates,baseline_hf,dt,hydro,hf,comm,ready, &
          correction_periodic,target_depth, &
          internal_energy_baseline=held_hydro_baseline,internal_energy_rate=held_hydro,internal_energy_flux=held_hydro_flux)
        else
        call apply_mpi_eb_flux_correction_3d(tree,owners,geometry,old,rates,baseline_hf,dt,hydro,hf,comm,ready, &
          correction_periodic,target_depth)
        endif
        if(.not.ready) return
        endif
      endif
      if(redistribution_method=='state_redist'.and..not.prepared_godunov) then
        if(redist_order>0) then
          call redist_boundary_samples(old,ot,ready)
          if(.not.ready) return
        endif
        if(selected_depth>=0) then
          call redistribute_sdc_mpi_eb_tree_state_level_rate_3d(tree,owners,geometry,old,hydro,nv,selected_depth,dt, &
            redistributed,hydro_transfer,comm,ready,redist_target,boundary_faces, &
      max_order=redist_order,boundary_halos=redist_halos,boundary_geometry=redist_geometry)
        else
          call redistribute_sdc_mpi_eb_tree_state_rate_3d(tree,owners,geometry,old,hydro,nv,dt,redistributed, &
            comm,ready,redist_target,boundary_faces, &
      max_order=redist_order,boundary_halos=redist_halos,boundary_geometry=redist_geometry)
        end if
        if(.not.ready) return
        hydro=redistributed
      else if(redistribution_method=='flux_redist') then
        if(selected_depth>=0) then
          call redistribute_sdc_mpi_eb_tree_level_rate_3d(tree,owners,geometry,hydro,nv,selected_depth,dt, &
            redistributed,hydro_transfer,comm,ready,boundary_faces)
        else
          call redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree,owners,geometry,hydro,nv,redistributed, &
            comm,ready,boundary_faces)
        end if
        if(.not.ready) return
        hydro=redistributed
      end if
      ! Publish no source family unless every owner completes native EB postprocessing.
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank.or..not.active(node)) cycle
        call post_eb_spatial_density_rate_3d(geometry(node),hydro%fields(node)%values,valid)
        if(.not.valid) exit
      enddo
      call agree(valid,comm,ready);if(.not.ready) return
      end if
      if(mol.and.present(internal_energy).and.native_held_funding.and..not.prepared_redist) then
        unsafe_forward=.false.
        do node=1,size(owners)
          if(owners(node)/=rank.or..not.active(node)) cycle
          n=tree%nodes(node)%shape
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
            forward_state=old%fields(node)%values(:,i,j,k)+ &
              dt*(rates%fields(node)%values(:,i,j,k)+hydro%fields(node)%values(:,i,j,k))
            unsafe_forward=unsafe_forward.or.any(.not.ieee_is_finite(forward_state)).or. &
              forward_state(1)<=0.or.any(forward_state(6:)<0)
          enddo
          enddo
          enddo
        enddo
        call MPI_Allreduce(unsafe_forward,correct_forward,1,MPI_LOGICAL,MPI_LOR,comm,status)
        if(status/=MPI_SUCCESS) return
        if(correct_forward) then
          call fund_held_native_stage(ready)
          if(.not.ready) return
        endif
      endif
      if(.not.mol.and.present(internal_energy)) then
        godunov_predictor_forcing=rates;godunov_predictor_hydro=hydro
        godunov_predictor_held=held_hydro;godunov_predictor_work=held_hydro_work
        godunov_predictor_faces=hf;godunov_predictor_scalar_faces=held_hydro_flux
      endif
      if(godunov_stage_capture) then
        call capture_godunov_hydro(it,ready)
        if(.not.ready) return
      endif
      do node=1,size(owners)
        if(owners(node)/=rank.or..not.active(node)) cycle
        candidate%fields(node)%values=old%fields(node)%values+dt*(rates%fields(node)%values+hydro%fields(node)%values)
      end do
      if(present(energy_correction)) then
        call held_update(held_old,held_rate0,held_rate0,.5_dp,held_candidate,ready)
        if(.not.ready) return
      endif
      call recover(ready,.false.,1)
      if(.not.ready) return
      end if
      if(present(internal_energy).and..not.present(energy_correction)) then
        call held_update(held_old,held_rate0,held_rate0,.5_dp,held_candidate,ready)
        if(.not.ready) return
      endif
      call sources(candidate,ct,e1,d1,a1,it,f1,transfer1,ready,fractions_end)
      if(.not.ready) then
        godunov_source_attempt_incomplete=.true.;return
      endif
      if(native_prefix_agreed) native_prefix_observed=max(native_prefix_observed,it+1)
      if(godunov_stage_capture) then
        call capture_godunov_source(it,candidate,ct,e1,d1,a1,f1,ready)
        if(.not.ready) then
          godunov_source_attempt_incomplete=.true.;return
        endif
      endif
      if(mol.and.it>1) then
        ! Pinned do_mol_advance rebuilds only the spatial new-time RHS on
        ! later reacting passes. molSrc_old retains old external forcing;
        ! the first predictor's new external forcing is not added again.
        do node=1,size(owners)
          if(owners(node)/=rank.or..not.active(node)) cycle
          e1%fields(node)%values=0
        enddo
      endif
      call combined_rates(.false.)
      prepared_reaction_done=.false.
      if(mol.and.chemistry.and.prepared_redist.and..not.present(replay_interval)) then
        call limit_prepared_mol_pass(it,ready,prepared_reaction_done)
        if(.not.ready) return
        call combined_rates(.false.)
      endif
      if(mol.and.chemistry.and..not.prepared_redist) then
        if(present(internal_energy).and.native_held_funding) then
          block
          type(amr_tree_state_3d) :: unreset_rates
          ! Fund the actual mean before a flux-only safeguard can refuse a
          ! volume-source deficit. Its temporary predictor reset is the same
          ! forcing applied below. Restore rates by assignment so an admitted
          ! unfunded family keeps its original represented arithmetic.
          unreset_rates=rates
          if(present(energy_correction).and.it==1.and..not.prepared_redist) then
            do node=1,size(owners)
            if(owners(node)/=rank.or..not.active(node)) cycle
              rates%fields(node)%values(5,:,:,:)=rates%fields(node)%values(5,:,:,:)+ &
                .5_dp*predictor_reset%fields(node)%values(5,:,:,:)/dt
            enddo
          endif
          call admit_held_native_corrector(ready)
          rates=unreset_rates
          if(.not.ready) return
          end block
        endif
        call limit_mol_pass(ready)
        if(.not.ready) return
        call combined_rates(.false.)
      endif
      if(present(energy_correction).and.it==1.and..not.prepared_redist) then
        accepted_predictor_reset=.5_dp*predictor_reset_change
        do node=1,size(owners)
          if(owners(node)/=rank.or..not.active(node)) cycle
          rates%fields(node)%values(5,:,:,:)=rates%fields(node)%values(5,:,:,:)+ &
            .5_dp*predictor_reset%fields(node)%values(5,:,:,:)/dt
        enddo
      endif
      if(present(internal_energy).and.native_held_funding.and..not.prepared_redist) then
        call admit_held_native_corrector(ready)
        if(.not.ready) return
      endif
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        rates%fields(node)%values=rates%fields(node)%values+hydro%fields(node)%values
      end do
      if(chemistry) then
        if(prepared_redist.and.mol) then
          if(.not.prepared_reaction_done) then
            call run_prepared_reaction(ready)
            if(.not.ready.and..not.reactor_refused.and..not.present(replay_interval).and.spatial_weight>0) then
              call limit_prepared_mol_pass(it,ready,prepared_reaction_done,force_search=.true.)
              if(.not.ready) return
              call combined_rates(.false.)
              if(.not.prepared_reaction_done) call run_prepared_reaction(ready)
            endif
            if(.not.ready) return
          endif
        else
        if(present(internal_energy)) then
        if(native_held_funding.and..not.prepared_redist) then
          call prepare_covered_reactor_forcing(ready)
          if(.not.ready) return
        endif
        reactor_refused_local=.false.
        call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old,ot,rates,dt,rtol,atol, &
          candidate,ct,ir,comm,ready,target_depth, &
          evaluate_reaction=godunov_reactor_adapter,frame=frame,heat_release=heat, &
          internal_energy=held_candidate,old_internal_energy=held_old)
        call MPI_Allreduce(reactor_refused_local,reactor_refused,1,MPI_LOGICAL,MPI_LOR,comm,status)
        if(status/=MPI_SUCCESS) return
        godunov_recovery_forbidden=godunov_recovery_forbidden.or.reactor_refused
        if(present(provider_refused)) provider_refused=reactor_refused
        else
        if(present(evaluate_reaction_admission)) then
          call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old,ot,rates,dt,rtol,atol, &
            candidate,ct,ir,comm,ready,target_depth, &
            evaluate_reaction=godunov_reactor_adapter,frame=frame,heat_release=heat)
        else
        call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old,ot,rates,dt,rtol,atol, &
          candidate,ct,ir,comm,ready,target_depth, &
          evaluate_reaction=evaluate_reaction,frame=frame,heat_release=heat)
        endif
        endif
        if(.not.ready) return
        if(present(energy_correction)) then
          call recover(ready,.true.,2)
          if(.not.ready) return
        endif
        endif
      else
        do node=1,size(owners)
          if(owners(node)/=rank.or..not.active(node)) cycle
          candidate%fields(node)%values=old%fields(node)%values+dt*rates%fields(node)%values
        end do
        if(present(energy_correction)) then
          call held_update(held_old,held_rate0,held_rate1,.5_dp,held_candidate,ready)
          if(.not.ready) return
        endif
        call recover(ready,.true.,2)
        if(.not.ready) return
      end if
      if(present(internal_energy).and..not.present(energy_correction)) then
        call held_update(held_old,held_rate0,held_rate1,.5_dp,held_candidate,ready)
        if(.not.ready) return
      endif
    end do
    projection_inventory=0;closure_inventory=0
    if(present(internal_energy)) then
      call measure_held_projection(ready)
      if(.not.ready) return
    endif
    ledger=old
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      ledger%fields(node)%values=0
      if(.not.active(node)) cycle
      ledger%fields(node)%values=candidate%fields(node)%values-(old%fields(node)%values+dt*rates%fields(node)%values)
    end do
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,ledger,chemical,comm,ready)
    if(.not.ready) return
    if(present(energy_correction)) chemical(5)=chemical(5)-final_reset
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,e0,v0,comm,ready)
    if(.not.ready) return
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,e1,v1,comm,ready)
    if(.not.ready) return
    external=.5_dp*dt*(v0+v1)
    if(present(auxiliary_source)) then
      if(held_auxiliary) then
        call integrate_mpi_eb_tree_3d(tree,owners,geometry,auxiliary_packed0,auxiliary_inventory0,comm,ready)
        if(.not.ready) return
        call integrate_mpi_eb_tree_3d(tree,owners,geometry,auxiliary_packed1,auxiliary_inventory1,comm,ready)
        if(.not.ready) return
        auxiliary=.5_dp*dt*(auxiliary_inventory0+auxiliary_inventory1)
      else
      call integrate_mpi_eb_tree_3d(tree,owners,geometry,a0,v0,comm,ready)
      if(.not.ready) return
      call integrate_mpi_eb_tree_3d(tree,owners,geometry,a1,v1,comm,ready)
      if(.not.ready) return
      auxiliary=.5_dp*dt*(v0+v1)
      endif
      call agree(all(ieee_is_finite(auxiliary)),comm,ready)
      if(.not.ready) return
    end if
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      if(.not.active(node)) cycle
      n=tree%nodes(node)%shape
      allocate(wall_rate(nv,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      call wall_source(old,ot,node,wall_rate,valid,0)
      if(.not.valid) exit
      ledger%fields(node)%values=dt*wall_rate
      deallocate(wall_rate)
      if(diffusive.or.has_diffusion.or.mol) then
        ! The accepted trapezoidal state uses both endpoint sources, including
        ! the single nonreacting MOL pass. Its flux ledger must use both too.
        if(allocated(admitted_mean_faces)) then
          hf(node)%x=hf(node)%x+admitted_mean_faces(node)%x
          hf(node)%y=hf(node)%y+admitted_mean_faces(node)%y
          hf(node)%z=hf(node)%z+admitted_mean_faces(node)%z
        else
          hf(node)%x=hf(node)%x+.5_dp*(f0(node)%x+f1(node)%x)
          hf(node)%y=hf(node)%y+.5_dp*(f0(node)%y+f1(node)%y)
          hf(node)%z=hf(node)%z+.5_dp*(f0(node)%z+f1(node)%z)
        endif
      end if
    end do
    call agree(valid,comm,ready)
    if(.not.ready) return
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,ledger,wall,comm,ready)
    if(.not.ready) return
    if(mol) wall=.5_dp*(mol_wall_old+mol_wall_new)
    wall(5)=wall(5)+.5_dp*(heat_old+heat_new)
    wall(2:4)=wall(2:4)+.5_dp*(impulse_old+impulse_new)
    call agree(all(ieee_is_finite(wall)),comm,ready)
    if(.not.ready) return
    if(present(registered_flux).or.godunov_stage_capture.or.(mol.and..not.chemistry.and.iterations>1)) then
      registration=hf
      if(mol.and..not.chemistry.and.iterations>1) then
        ! No later reacting pass exists to register the final new-time half.
        ! Assemble directly from f0, avoiding subtraction cancellation.
        do node=1,size(registration)
          if(owners(node)/=rank.or..not.active(node)) cycle
          registration(node)%x=.5_dp*f0(node)%x
          registration(node)%y=.5_dp*f0(node)%y
          registration(node)%z=.5_dp*f0(node)%z
        enddo
      endif
    endif
    if(mol.and..not.chemistry.and.iterations>1) then
      ! The reference never registers the provisional new half when chemistry
      ! is disabled and mol_iters>1. A selected level returns that register;
      ! a simultaneous hierarchy update synchronizes its full physical flux.
      physical_flux=hf
      if(selected_depth>=0) hf=registration
    endif
    if(present(internal_energy)) then
      call ensure_selected_scalar_faces(ready)
      if(.not.ready) return
      held_flux=held_flux0
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank.or.(.not.active(node).and..not.native_held_funding)) cycle
        if(allocated(admitted_held_faces)) then
          held_flux(node)=admitted_held_faces(node)
        else
        held_flux(node)%x=.5_dp*(held_flux0(node)%x+held_flux1(node)%x)
        held_flux(node)%y=.5_dp*(held_flux0(node)%y+held_flux1(node)%y)
        held_flux(node)%z=.5_dp*(held_flux0(node)%z+held_flux1(node)%z)
        endif
        if(.not.mol.or.mol_source_funded) then
          held_flux(node)%x=held_flux(node)%x+held_hydro_flux(node)%x
          held_flux(node)%y=held_flux(node)%y+held_hydro_flux(node)%y
          held_flux(node)%z=held_flux(node)%z+held_hydro_flux(node)%z
        endif
        valid=valid.and.all(ieee_is_finite(held_flux(node)%x)).and. &
          all(ieee_is_finite(held_flux(node)%y)).and.all(ieee_is_finite(held_flux(node)%z))
      enddo
      call agree(valid,comm,ready)
      if(.not.ready) return
      if(prepared_redist) then
        if(prepared_godunov) then
          do node=1,size(held_registers0)
            if(allocated(held_registers0(node)%fine_delta)) &
              held_registers0(node)%fine_delta=.5_dp*held_registers0(node)%fine_delta
            if(allocated(held_registers1(node)%fine_delta)) &
              held_registers1(node)%fine_delta=.5_dp*held_registers1(node)%fine_delta
          enddo
        endif
        held_registers=held_registers0;valid=.true.
        if(selected_depth>=0.and.mol.and..not.chemistry.and.iterations>1) then
          do node=1,size(held_registers1)
            if(owners(node)/=rank) cycle
            if(allocated(held_registers1(node)%fine_delta)) held_registers1(node)%fine_delta=0
          enddo
        endif
        do node=1,size(owners)
          if(owners(node)/=rank.or..not.active(node)) cycle
          n=tree%nodes(node)%shape
          call accumulate_state_redist_register_3d(held_registers(node),held_registers1(node),1,n+12,ready)
          valid=valid.and.ready
          if(prepared_godunov) then
            call accumulate_state_redist_register_3d(held_registers(node),held_hydro_registers(node),1,n+12,ready)
            valid=valid.and.ready
          endif
        enddo
        call agree(valid,comm,ready)
        if(.not.ready) return
      endif
    endif
    if(prepared_redist) then
      final_registers=stage_registers0;valid=.true.
      ! Selected-level registration follows the pinned zero new-time factor;
      ! a simultaneous hierarchy update retains both physical half steps.
      if(selected_depth>=0.and.mol.and..not.chemistry.and.iterations>1) then
        do node=1,size(stage_registers1)
          if(allocated(stage_registers1(node)%fine_delta)) stage_registers1(node)%fine_delta=0
        enddo
      endif
      if(prepared_godunov) then
        do node=1,size(owners)
          if(owners(node)/=rank) cycle
          if(allocated(final_registers(node)%fine_delta)) final_registers(node)%fine_delta= &
            .5_dp*final_registers(node)%fine_delta
          if(allocated(stage_registers1(node)%fine_delta)) stage_registers1(node)%fine_delta= &
            .5_dp*stage_registers1(node)%fine_delta
        enddo
      endif
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        if(.not.allocated(final_registers(node)%coarse_delta).and. &
          .not.allocated(final_registers(node)%fine_delta)) cycle
        n=tree%nodes(node)%shape
        ! Match the accepted endpoint source in both coarse and fine ledgers.
        call accumulate_state_redist_register_3d(final_registers(node),stage_registers1(node),nv,n+12,ready)
        valid=valid.and.ready
        if(prepared_godunov) then
          call accumulate_state_redist_register_3d(final_registers(node),hydro_registers(node),nv,n+12,ready)
          valid=valid.and.ready
        endif
      enddo
      call agree(valid,comm,ready)
      if(.not.ready) return
    endif
    if(selected_depth<0.and.present(internal_energy).and.present(energy_correction)) then
      if(tree%level_count()==1) then
        call close_held_prepared_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,held_candidate,hf,held_flux, &
          held_policy,post_reset,0,comm,ready,boundary_faces,frame, &
          register_geometry=extended_geometry,internal_energy_change=closure_inventory(1), &
          reset_current_level=mol.and.chemistry.and.redistribution_method=='none')
        if(.not.ready) return
      endif
      do sync_depth=tree%level_count()-2,0,-1
        if(prepared_redist) then
          call close_held_prepared_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,held_candidate,hf,held_flux, &
            held_policy,post_reset,sync_depth,comm,ready,boundary_faces,frame,register_geometry, &
            final_registers,held_registers, &
          internal_energy_change=closure_inventory(1))
        else
          call close_held_prepared_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,held_candidate,hf,held_flux, &
            held_policy,post_reset,sync_depth,comm,ready,boundary_faces,frame, &
          register_geometry=extended_geometry,internal_energy_change=closure_inventory(1), &
          reset_current_level=mol.and.chemistry.and.redistribution_method=='none')
        endif
        if(.not.ready) return
      enddo
      valid=.true.
      if(prepared_redist) then
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        if(allocated(final_registers(node)%coarse_delta)) then
          if(any(final_registers(node)%coarse_delta/=0)) valid=.false.
        endif
        if(allocated(final_registers(node)%fine_delta)) then
          if(any(final_registers(node)%fine_delta/=0)) valid=.false.
        endif
        if(allocated(held_registers(node)%coarse_delta)) then
          if(any(held_registers(node)%coarse_delta/=0)) valid=.false.
        endif
        if(allocated(held_registers(node)%fine_delta)) then
          if(any(held_registers(node)%fine_delta/=0)) valid=.false.
        endif
      enddo
      endif
      call agree(valid,comm,ready)
      if(.not.ready) return
    else if(selected_depth<0.and.prepared_redist) then
      if(tree%level_count()==1) then
        call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,hf,comm,ready, &
          0,0,boundary_faces=boundary_faces,frame=frame)
        if(.not.ready) return
      endif
      do sync_depth=tree%level_count()-2,0,-1
        call synchronize_mpi_eb_tree_register_flux_3d(species,tree,owners,geometry,register_geometry, &
          candidate,ct,hf,final_registers,sync_depth,comm,ready,boundary_faces,frame)
        if(.not.ready) return
      enddo
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        if(allocated(final_registers(node)%coarse_delta)) then
          if(any(final_registers(node)%coarse_delta/=0)) valid=.false.
        endif
        if(allocated(final_registers(node)%fine_delta)) then
          if(any(final_registers(node)%fine_delta/=0)) valid=.false.
        endif
      enddo
      call agree(valid,comm,ready)
      if(.not.ready) return
    else
    if(selected_depth<0) then
      if(redistribution_method/='none') before_sync=candidate
      call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,hf,comm,ready, &
        boundary_faces=boundary_faces,frame=frame)
      if(.not.ready) then
        ! At a common endpoint, finalize adjacent levels deepest first. Their
        ! canonical faces become the baseline of the next synchronization, so
        ! an already accepted fine correction is never applied twice.
        do sync_depth=tree%level_count()-2,0,-1
          call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,hf,comm,ready, &
            sync_depth,sync_depth+1,boundary_faces=boundary_faces,frame=frame)
          if(.not.ready) then
            call repair_mpi_eb_tree_interface_budget_3d(species,tree,owners,geometry,candidate,ct,hf,comm,ready, &
              sync_depth,boundary_faces,frame)
            if(.not.ready) return
            call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,hf,comm,ready, &
              sync_depth,sync_depth+1,boundary_faces=boundary_faces,frame=frame)
          endif
          if(.not.ready) return
        enddo
        call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,hf,comm,ready, &
          boundary_faces=boundary_faces,frame=frame)
      endif
      if(.not.ready) return
      if(redistribution_method/='none') then
        sync_rate=hydro
        do node=1,size(owners)
          if(owners(node)/=rank) cycle
          sync_rate%fields(node)%values=(candidate%fields(node)%values-before_sync%fields(node)%values)/dt
        end do
        if(redistribution_method=='state_redist') then
          call redistribute_mpi_eb_tree_state_3d(tree,owners,geometry,sync_rate,nv,redistributed, &
            comm,ready,redist_target,boundary_faces, &
      max_order=redist_order)
        else
          call redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree,owners,geometry,sync_rate,nv,redistributed, &
            comm,ready,boundary_faces)
        end if
        if(.not.ready) return
        do node=1,size(owners)
          if(owners(node)/=rank) cycle
          candidate%fields(node)%values=before_sync%fields(node)%values+dt*redistributed%fields(node)%values
        end do
        call recover(ready,.true.)
        if(.not.ready) return
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,ready,frame)
        if(.not.ready) return
      end if
    end if
    endif
    if(present(internal_energy)) then
      if(selected_depth<0.and..not.present(energy_correction)) then
        if(prepared_redist) then
          if(tree%level_count()==1) then
            call synchronize_mpi_eb_tree_scalar_flux_3d(species,tree,owners,geometry,held_candidate,held_flux,comm,ready, &
              0,0,boundary_faces)
            if(.not.ready) return
          endif
          do sync_depth=tree%level_count()-2,0,-1
            call synchronize_mpi_eb_tree_scalar_register_flux_3d(species,tree,owners,geometry,register_geometry, &
              held_candidate,held_flux,held_registers,sync_depth,comm,ready,boundary_faces)
            if(.not.ready) return
          enddo
          valid=.true.
          do node=1,size(held_registers)
            if(owners(node)/=rank) cycle
            if(allocated(held_registers(node)%coarse_delta)) then
              if(any(held_registers(node)%coarse_delta/=0)) valid=.false.
            endif
            if(allocated(held_registers(node)%fine_delta)) then
              if(any(held_registers(node)%fine_delta/=0)) valid=.false.
            endif
          enddo
          call agree(valid,comm,ready)
          if(.not.ready) return
        else
          call synchronize_mpi_eb_tree_scalar_flux_3d(species,tree,owners,geometry,held_candidate,held_flux,comm,ready, &
            boundary_faces=boundary_faces)
          if(.not.ready) return
        endif
      endif
    endif
    if(present(energy_correction)) then
      call recover(ready,.true.,3)
      if(.not.ready) return
      if(selected_depth<0.and..not.(mol.and.chemistry.and.redistribution_method=='none')) then
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,ready,frame,held_candidate)
        if(.not.ready) return
      endif
    endif
    if(allocated(physical_flux).and.selected_depth>=0) call move_alloc(physical_flux,hf)
    boundary=0
    if(rank==owners(1).and.selected_depth<=0) then
      n=tree%nodes(1)%shape
      associate(g=>geometry(1),f=>hf(1))
      do c=1,nv
        boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
          g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
          g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
          g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
          g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
          g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
      end do
      end associate
    end if
    call MPI_Bcast(boundary,nv,MPI_DOUBLE_PRECISION,owners(1),comm,status)
    call agree(status==MPI_SUCCESS.and.all(ieee_is_finite(boundary)).and.all(ieee_is_finite(external)),comm,ready)
    if(.not.ready) return
    call predictor_field(e1,d1,a1,np,ready)
    if(.not.ready) return
    if(present(redistribution_transfer)) then
      final_transfer=hydro_transfer
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        if(allocated(admitted_mean_transfer%fields)) then
          final_transfer%fields(node)%values=hydro_transfer%fields(node)%values+ &
            admitted_mean_transfer%fields(node)%values
        else
          final_transfer%fields(node)%values=hydro_transfer%fields(node)%values+ &
            .5_dp*(transfer0%fields(node)%values+transfer1%fields(node)%values)
        endif
        valid=valid.and.all(ieee_is_finite(final_transfer%fields(node)%values))
      end do
      call agree(valid,comm,ready)
      if(.not.ready) return
    end if
    if(present(heat_release).and.selected_depth<0) then
      call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,heat,1,comm,ready)
      if(.not.ready) return
    endif
    if(present(reaction_source)) then
      do node=1,size(owners)
        if(owners(node)/=rank.or..not.active(node)) cycle
        if(.not.chemistry) cycle
        do c=1,nv-4
          where(geometry(node)%volume_fraction<=0) &
            ir%fields(node)%values(c,:,:,:)=reaction_source%fields(node)%values(c,:,:,:)
        end do
      end do
      if(selected_depth<0) then
        call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,ir,nv-4,comm,ready)
        if(.not.ready) return
      endif
    end if
    if(present(accepted_mass_source).or.present(interval_history)) then
      ledger=rates
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        ledger%fields(node)%values=0
        if(.not.active(node)) cycle
        ledger%fields(node)%values=.5_dp*(e0%fields(node)%values+e1%fields(node)%values)
        if(present(auxiliary_source)) ledger%fields(node)%values=ledger%fields(node)%values+ &
          .5_dp*(a0%fields(node)%values+a1%fields(node)%values)
        ledger%fields(node)%values(2:5,:,:,:)=0
      enddo
    endif
    if(present(energy_correction)) then
      numerical_change=initial_reset+final_reset+post_reset+accepted_predictor_reset
      call agree(ieee_is_finite(energy_correction+numerical_change),comm,ready)
      if(.not.ready) return
    endif
    if(present(source_funding_impulse)) then
      funding_receipt=source_funding_stock
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        funding_receipt%fields(node)%values=0
        if(active(node)) cycle
        funding_receipt%fields(node)%values(:nv,:,:,:)=dt*hydro%fields(node)%values
        funding_receipt%fields(node)%values(nv+1,:,:,:)=dt*held_hydro%fields(node)%values(1,:,:,:)
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(covered(tree,node,[i,j,k])) funding_receipt%fields(node)%values(:,i,j,k)=0
        enddo
        enddo
        enddo
        valid=valid.and.all(ieee_is_finite(funding_receipt%fields(node)%values))
      enddo
      call agree(valid,comm,ready);if(.not.ready) return
    endif
    if(present(interval_history)) then
      call capture_mol_interval(ready)
      if(.not.ready) return

    endif
    pressure_work_change=0
    if((present(internal_energy_pressure_work_change).or.godunov_stage_capture).and..not.mol) then
      call integrate_mpi_eb_tree_3d(tree,owners,geometry,held_hydro_work,pressure_work_inventory,comm,ready)
      if(.not.ready) return
      pressure_work_change=-dt*pressure_work_inventory(1)
      call agree(ieee_is_finite(pressure_work_change),comm,ready);if(.not.ready) return
    endif
    if(godunov_stage_capture) then
      call finish_godunov_history(ready)
      if(.not.ready) return
      if(present(godunov_interval_history)) &
        call move_godunov_interval_3d(godunov_history,godunov_interval_history)
    endif
    if(present(interval_history)) call move_mol_interval(pending_history,interval_history)
    if(present(source_funding_impulse)) call move_alloc(funding_receipt%fields,source_funding_impulse%fields)
    if(present(reaction_source)) reaction_source=ir
    if(present(heat_release)) heat_release=heat
    if(present(auxiliary_change)) auxiliary_change=auxiliary
    if(present(redistribution_transfer)) redistribution_transfer=final_transfer
    if(present(state_redist_registers)) call move_alloc(final_registers,state_redist_registers)
    if(present(held_state_redist_registers)) call move_alloc(held_registers,held_state_redist_registers)
    if(present(integrated_flux)) integrated_flux=hf
    if(present(registered_flux)) call move_alloc(registration,registered_flux)
    if(present(energy_correction)) energy_correction=energy_correction+numerical_change
    if(present(internal_energy_pressure_work_change)) internal_energy_pressure_work_change=pressure_work_change
    if(present(internal_energy_projection_change)) internal_energy_projection_change=projection_inventory(1)
    if(present(internal_energy_closure_change)) internal_energy_closure_change=closure_inventory(1)
    if(present(internal_energy)) then
      call move_alloc(held_candidate%fields,internal_energy%fields)
      call move_alloc(held_flux,internal_energy_flux)
    endif
    if(present(accepted_forcing)) accepted_forcing=rates
    if(present(old_spatial_flux)) old_spatial_flux=f0
    if(present(accepted_mass_source)) accepted_mass_source=ledger
    if(present(old_predictor)) old_predictor=op
    if(present(new_predictor)) new_predictor=np
    call move_alloc(candidate%fields,state%fields);call move_alloc(ct%fields,temperature%fields)
    boundary_change=boundary;wall_change=wall;chemical_change=chemical;source_change=external
    minimum_theta=theta;ok=.true.
    end subroutine run_common
    logical function cached_godunov_auxiliary(stage) result(cached)
      integer, intent(in) :: stage
      cached=.false.
      if(present(godunov_replay_interval)) then
        if(stage>=0.and.stage<=godunov_replay_interval%iterations) &
          cached=godunov_replay_interval%sources(stage)%auxiliary_accepted
      endif
    end function cached_godunov_auxiliary
    logical function cached_godunov_provider(stage) result(cached)
      integer, intent(in) :: stage
      cached=.false.
      if(present(godunov_replay_interval)) then
        if(stage>=0.and.stage<=godunov_replay_interval%iterations) &
          cached=godunov_replay_interval%sources(stage)%provider_accepted
      endif
    end function cached_godunov_provider
    subroutine validate_native_source_replay(good)
      logical, intent(out) :: good
      logical :: valid
      integer :: minimum,maximum
      ! Profile agreement precedes callback/spatial work, including idle ranks.
      valid=.not.present(replay_interval).and..not.native_external_frozen
      call agree(valid,comm,good);if(.not.good) return
      valid=valid_native_mol_source_prefix_3d(native_source_replay,tree,geometry,selected_depth,dt,iterations, &
        chemistry,transport_flags,force,accel,held_policy,old,ot,held_old,ir,owners==rank, &
        energy_correction+initial_reset,present(auxiliary_source),held_auxiliary,frame,source_funding_stock, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected)
      call agree(valid,comm,good);if(.not.good) return
      call MPI_Allreduce(native_source_replay%accepted_source_stages,minimum,1,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) then
        good=.false.;return
      endif
      call MPI_Allreduce(native_source_replay%accepted_source_stages,maximum,1,MPI_INTEGER,MPI_MAX,comm,status)
      good=status==MPI_SUCCESS.and.minimum==maximum
    end subroutine validate_native_source_replay

    logical function cached_native_source(stage) result(cached)
      integer, intent(in) :: stage
      cached=.false.
      if(present(native_source_replay)) cached=stage<native_source_replay%accepted_source_stages
    end function cached_native_source

    subroutine publish_native_source_prefix()
      type(eb_native_mol_source_prefix_3d) :: result
      logical :: ready,local_ok
      integer :: stage
      local_ok=native_prefix_observed>0.and..not.godunov_recovery_forbidden.and..not.godunov_packet_pending.and. &
        allocated(pending_history%external_stages)
      call agree(local_ok,comm,ready);if(.not.ready) return
      result%depth=selected_depth;result%iterations=iterations
      result%accepted_source_stages=native_prefix_observed;result%interval=dt;result%force=force
      result%separate_momentum=selected_momentum
      result%consistent_velocity=selected_velocity
      result%paired_pressure=selected_pressure
      result%connected_directions=selected_connected
      result%chemistry=chemistry;result%transport_flags=transport_flags
      result%has_auxiliary_source=present(auxiliary_source);result%auxiliary_has_internal_energy=held_auxiliary
      result%has_frame=present(frame);if(present(frame)) result%frame=frame
      result%policy=held_policy;result%topology=tree;result%acceleration=accel
      result%geometry=geometry
      result%initial_energy_correction=energy_correction+initial_reset
      call copy_field(old,result%initial_gas,ready,all_levels=.true.)
      if(.not.ready) return
      call copy_field(ot,result%initial_temperature,ready,all_levels=.true.)
      if(.not.ready) return
      call copy_field(held_old,result%initial_internal_energy,ready,all_levels=.true.)
      if(.not.ready) return
      result%initial_reaction_source=native_prefix_initial_reaction_source
      result%has_source_funding=present(source_funding_stock)
      if(present(source_funding_stock)) result%source_funding_stock=source_funding_stock
      allocate(result%external_stages,source=pending_history%external_stages,stat=ios)
      call agree(ios==0,comm,ready);if(.not.ready) return
      if(allocated(pending_history%auxiliary_stages)) then
        allocate(result%auxiliary_stages,source=pending_history%auxiliary_stages,stat=ios)
        call agree(ios==0,comm,ready);if(.not.ready) return
      endif
      if(present(native_source_replay)) then
        result%accepted_source_stages=max(result%accepted_source_stages,native_source_replay%accepted_source_stages)
        do stage=0,native_source_replay%accepted_source_stages-1
          result%external_stages(stage)=native_source_replay%external_stages(stage)
          if(result%has_auxiliary_source) result%auxiliary_stages(stage)=native_source_replay%auxiliary_stages(stage)
        enddo
      endif
      do stage=result%accepted_source_stages,iterations
        if(allocated(result%external_stages(stage)%fields)) deallocate(result%external_stages(stage)%fields)
        if(allocated(result%auxiliary_stages)) then
          if(allocated(result%auxiliary_stages(stage)%fields)) deallocate(result%auxiliary_stages(stage)%fields)
        endif
      enddo
      ! A declined or malformed declared packet is not a numerical retry.
      native_source_prefix=result;native_source_prefix_ready=.true.
    end subroutine publish_native_source_prefix

    subroutine publish_source_prefix()
      logical :: local_ok,ready
      local_ok=allocated(godunov_history)
      call agree(local_ok,comm,ready);if(.not.ready) return
      local_ok=.not.godunov_recovery_forbidden.and..not.godunov_packet_pending
      call agree(local_ok,comm,ready);if(.not.ready) return
      local_ok=godunov_history%accepted_source_stages>0
      if(godunov_history%accepted_source_stages<=iterations) then
        associate(front=>godunov_history%sources(godunov_history%accepted_source_stages))
        local_ok=local_ok.or.front%auxiliary_accepted.or.front%provider_accepted
        end associate
      endif
      call agree(local_ok,comm,ready);if(.not.ready) return
      godunov_history%source_prefix_only=.true.
      local_ok=valid_godunov_source_cache_3d(godunov_history,tree,nv,naux,rank,owners,allow_partial=.true.)
      call agree(local_ok,comm,ready);if(.not.ready) return
      call move_godunov_interval_3d(godunov_history,godunov_source_recovery)
      godunov_source_recovery_ready=.true.
    end subroutine publish_source_prefix
    subroutine godunov_reactor_adapter(species,reactions,energy,source_species,source_energy,duration,rtol,atol, &
        density,temperature,steps,ok,maximum_steps)
      type(nasa7_species), intent(in) :: species(:)
      type(elementary_reaction), intent(in) :: reactions(:)
      real(dp), intent(in) :: energy,source_species(:),source_energy,duration,rtol,atol
      real(dp), intent(inout) :: density(:),temperature
      integer, intent(out) :: steps
      logical, intent(out) :: ok
      integer, optional, intent(in) :: maximum_steps
      logical :: refused
      if(present(evaluate_reaction)) then
        call evaluate_classified_reactor(evaluate_reaction,species,reactions,energy,source_species,source_energy, &
          duration,rtol,atol,density,temperature,steps,ok,refused,maximum_steps, &
          evaluate_admission=evaluate_reaction_admission,admission_status=admission_state)
        reactor_refused_local=reactor_refused_local.or.refused
      else
        call advance_forced_reactor_interval(species,reactions,energy,source_species,source_energy,duration,rtol,atol, &
          density,temperature,steps,ok,maximum_steps=maximum_steps)
      endif
    end subroutine godunov_reactor_adapter
    subroutine capture_auxiliary_packet(stage,packet,good)
      integer, intent(in) :: stage
      type(amr_tree_state_3d), intent(in) :: packet
      logical, intent(out) :: good
      logical :: local_ok
      call copy_godunov_field_3d(packet,godunov_history%sources(stage)%auxiliary,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      godunov_history%sources(stage)%auxiliary_accepted=.true.
      godunov_history%sources(stage)%index=stage
    end subroutine capture_auxiliary_packet
    subroutine capture_provider_packet(stage,kind,rate,halos,faces,limiter,good)
      integer, intent(in) :: stage,kind
      type(amr_tree_state_3d), intent(in) :: rate,halos
      type(eb_integrated_flux_3d), allocatable, intent(in) :: faces(:)
      real(dp), intent(in) :: limiter
      logical, intent(out) :: good
      logical :: local_ok
      associate(saved=>godunov_history%sources(stage))
      call copy_godunov_field_3d(rate,saved%provider_rate,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      if(kind==1) then
        call copy_godunov_field_3d(halos,saved%provider_halos,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      endif
      call copy_godunov_fluxes_3d(faces,saved%provider_faces,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      saved%provider_kind=kind;saved%provider_limiter=limiter
      saved%provider_accepted=.true.;saved%index=stage
      end associate
      good=.true.
    end subroutine capture_provider_packet
    subroutine validate_godunov_replay(good)
      logical, intent(out) :: good
      logical :: local_ok
      integer :: stage,p,d(3),expected,header(2),lower(2),upper(2)
      good=.false.
      associate(saved=>godunov_replay_interval)
      local_ok=valid_godunov_source_cache_3d(saved,tree,nv,naux,rank,owners,allow_partial=godunov_partial)
      call agree(local_ok,comm,good);if(.not.good) return

      call agree_held_projection(saved%internal_energy_projection_change,comm,good)
      if(.not.good) return
      call agree_held_projection(saved%internal_energy_closure_change,comm,good)
      if(.not.good) return
      call agree_held_clock_packet([saved%initial_internal_energy_change, &
        saved%held_clock_start_change,saved%held_clock_end_change, &
        saved%held_clock_projection_before,saved%held_clock_reset_before, &
        saved%held_clock_projection_after,saved%held_clock_reset_after, &
        merge(1._dp,0._dp,saved%held_clock_events_closed), &
        merge(1._dp,0._dp,saved%held_clock_prefixes_valid)],comm,good)
      if(.not.good) return
      header=[saved%accepted_source_stages,merge(1,0,saved%source_prefix_only)]
      call MPI_Allreduce(header,lower,2,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) then
        good=.false.;return
      endif
      call MPI_Allreduce(header,upper,2,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(lower/=upper)) then
        good=.false.;return
      endif

      local_ok=allocated(saved%topology%nodes).and.allocated(saved%sources)
      call agree(local_ok,comm,good);if(.not.good) return
      local_ok=size(saved%topology%nodes)==size(owners).and.lbound(saved%sources,1)==0.and. &
        ubound(saved%sources,1)==iterations.and.saved%depth==selected_depth.and.saved%iterations==iterations.and. &
        ieee_is_finite(saved%interval).and.all(ieee_is_finite(saved%force))
      call agree(local_ok,comm,good);if(.not.good) return
      local_ok=saved%interval==dt.and.all(saved%force==force).and.(saved%chemistry.eqv.chemistry).and. &
        all(saved%transport_flags.eqv.transport_flags).and.(saved%prepared_redistribution.eqv.prepared_redist).and. &
        (saved%has_auxiliary.eqv.present(auxiliary_source)).and. &
        (saved%has_spatial_provider.eqv.present(evaluate_spatial_diffusion)).and. &
        (saved%spatial_diffusion_is_source.eqv.external_diffusion_source).and. &
        (saved%auxiliary_has_internal_energy.eqv.held_auxiliary)
      local_ok=local_ok.and.(saved%has_valid_diffusion_provider.eqv.present(evaluate_diffusion))
      call agree(local_ok,comm,good);if(.not.good) return
      do stage=0,iterations
        header=[merge(1,0,saved%sources(stage)%auxiliary_accepted),merge(1,0,saved%sources(stage)%provider_accepted)]
        call MPI_Allreduce(header,lower,2,MPI_INTEGER,MPI_MIN,comm,status)
        if(status/=MPI_SUCCESS) then
          good=.false.;return
        endif
        call MPI_Allreduce(header,upper,2,MPI_INTEGER,MPI_MAX,comm,status)
        if(status/=MPI_SUCCESS.or.any(lower/=upper)) then
          good=.false.;return
        endif
      enddo
      if(.not.saved%source_prefix_only) then
        do stage=1,iterations
          header=[merge(1,0,saved%hydro(stage)%source_funded), &
            merge(1,0,saved%hydro(stage)%coarse_registers_accumulated)]
          call MPI_Allreduce(header,lower,2,MPI_INTEGER,MPI_MIN,comm,status)
          if(status/=MPI_SUCCESS) then
            good=.false.;return
          endif
          call MPI_Allreduce(header,upper,2,MPI_INTEGER,MPI_MAX,comm,status)
          if(status/=MPI_SUCCESS.or.any(lower/=upper)) then
            good=.false.;return
          endif
        enddo
      endif
      local_ok=saved%has_source_funding.eqv.present(source_funding_stock)
      call agree(local_ok,comm,good);if(.not.good) return
      if(present(source_funding_stock)) then
        local_ok=.true.
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          local_ok=local_ok.and.all(saved%source_funding_stock%fields(p)%values== &
            source_funding_stock%fields(p)%values)
        enddo
        call agree(local_ok,comm,good);if(.not.good) return
      endif
      do p=1,size(owners)
        local_ok=saved%topology%nodes(p)%parent==tree%nodes(p)%parent.and. &
          saved%topology%nodes(p)%level==tree%nodes(p)%level.and. &
          all(saved%topology%nodes(p)%shape==tree%nodes(p)%shape)
      call agree(local_ok,comm,good);if(.not.good) return
      enddo
      do stage=0,iterations
        if(.not.(saved%sources(stage)%auxiliary_accepted.or.saved%sources(stage)%provider_accepted).and. &
          stage>=saved%accepted_source_stages) cycle
        local_ok=saved%sources(stage)%index==stage
      call agree(local_ok,comm,good);if(.not.good) return
        if(saved%sources(stage)%auxiliary_accepted) then
          call valid_cached_field(saved%sources(stage)%auxiliary,naux,0,good)
          if(.not.good) return
        endif
        expected=0
        if(present(evaluate_spatial_diffusion)) then
          expected=1
          if(external_diffusion_source.and.stage>1.and..not.prepared_redist) expected=0
        endif
        if(present(evaluate_diffusion)) expected=2
        if(.not.saved%sources(stage)%provider_accepted) expected=0
        local_ok=saved%sources(stage)%provider_kind==expected
      call agree(local_ok,comm,good);if(.not.good) return
        if(expected==0) cycle
        call valid_cached_field(saved%sources(stage)%provider_rate,nv,0,good)
        if(.not.good) return
        if(expected==1) then
          call valid_cached_field(saved%sources(stage)%provider_halos,nv,6,good)
          if(.not.good) return
        endif
        local_ok=allocated(saved%sources(stage)%provider_faces).and.ieee_is_finite(saved%sources(stage)%provider_limiter)
      call agree(local_ok,comm,good);if(.not.good) return
        local_ok=size(saved%sources(stage)%provider_faces)==size(owners).and. &
          saved%sources(stage)%provider_limiter>=0.and.saved%sources(stage)%provider_limiter<=1
      call agree(local_ok,comm,good);if(.not.good) return
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          if(.not.active(p)) cycle
          local_ok=valid_eb_integrated_flux_3d(saved%sources(stage)%provider_faces(p),geometry(p),nv)
          if(.not.local_ok) exit
          if(expected==1) then
            d=tree%nodes(p)%shape
            local_ok=all(saved%sources(stage)%provider_halos%fields(p)%values(:,7:d(1)+6,7:d(2)+6,7:d(3)+6)== &
              saved%sources(stage)%provider_rate%fields(p)%values)
            if(.not.local_ok) exit
          endif
        enddo
      call agree(local_ok,comm,good);if(.not.good) return
      enddo
      end associate
      good=.true.
    end subroutine validate_godunov_replay
    subroutine valid_cached_field(field,components,padding,good)
      type(amr_tree_state_3d), intent(in) :: field
      integer, intent(in) :: components,padding
      logical, intent(out) :: good
      logical :: local_ok
      integer :: p,d(3)
      local_ok=allocated(field%fields)
      call agree(local_ok,comm,good);if(.not.good) return
      local_ok=size(field%fields)==size(owners).and.lbound(field%fields,1)==1
      call agree(local_ok,comm,good);if(.not.good) return
      do p=1,size(owners)
        if(owners(p)/=rank) then
          if(allocated(field%fields(p)%values)) local_ok=.false.
          cycle
        endif
        if(padding>0.and..not.active(p)) then
          if(allocated(field%fields(p)%values)) local_ok=.false.
          cycle
        endif
        if(.not.allocated(field%fields(p)%values)) then
          local_ok=.false.;exit
        endif
        d=tree%nodes(p)%shape
        local_ok=all(shape(field%fields(p)%values)==[components,d+2*padding]).and. &
          all(lbound(field%fields(p)%values)==1)
        if(.not.local_ok) exit
        local_ok=all(ieee_is_finite(field%fields(p)%values))
        if(.not.local_ok) exit
      enddo
      call agree(local_ok,comm,good);if(.not.good) return
      good=.true.
    end subroutine valid_cached_field
    subroutine history_copy_flux(source,destination,good)
      type(eb_integrated_flux_3d), allocatable, intent(in) :: source(:)
      type(eb_integrated_flux_3d), allocatable, intent(out) :: destination(:)
      logical, intent(out) :: good
      integer :: p,allocation_status
      logical :: valid
      good=.true.
      if(.not.allocated(source)) return
      allocate(destination(size(source)),stat=allocation_status)
      call agree(allocation_status==0,comm,good);if(.not.good) return
      valid=.true.
      do p=1,size(source)
        if(allocated(source(p)%x)) then
          allocate(destination(p)%x,source=source(p)%x,stat=allocation_status)
          if(allocation_status/=0) then
            valid=.false.;exit
          endif
        endif
        if(allocated(source(p)%y)) then
          allocate(destination(p)%y,source=source(p)%y,stat=allocation_status)
          if(allocation_status/=0) then
            valid=.false.;exit
          endif
        endif
        if(allocated(source(p)%z)) then
          allocate(destination(p)%z,source=source(p)%z,stat=allocation_status)
          if(allocation_status/=0) then
            valid=.false.;exit
          endif
        endif
      enddo
      call agree(valid,comm,good)
    end subroutine history_copy_flux
    subroutine begin_godunov_history(good)
      logical, intent(out) :: good
      allocate(godunov_history,stat=history_status)
      call agree(history_status==0,comm,good);if(.not.good) return
      allocate(godunov_history%sources(0:iterations),godunov_history%hydro(iterations),stat=history_status)
      call agree(history_status==0,comm,good);if(.not.good) return
      godunov_history%depth=selected_depth;godunov_history%interval=dt;godunov_history%iterations=iterations
      godunov_history%force=force;godunov_history%chemistry=chemistry;godunov_history%transport_flags=transport_flags
      godunov_history%prepared_redistribution=prepared_redist;godunov_history%policy=held_policy
      if(present(spatial_flux_weight)) godunov_history%face_weight=spatial_flux_weight
      godunov_history%replayed=present(godunov_replay_interval)
      godunov_history%has_auxiliary=present(auxiliary_source)
      godunov_history%has_source_funding=present(source_funding_stock)
      if(present(source_funding_stock)) godunov_history%source_funding_stock=source_funding_stock
      godunov_history%has_spatial_provider=present(evaluate_spatial_diffusion)
      godunov_history%has_valid_diffusion_provider=present(evaluate_diffusion)
      godunov_history%has_predictor_provider=present(predictor_builder)
      godunov_history%spatial_diffusion_is_source=external_diffusion_source
      godunov_history%auxiliary_has_internal_energy=held_auxiliary
      if(present(energy_correction)) godunov_history%initial_energy_correction=energy_correction
      godunov_history%topology=tree
      call copy_field(state,godunov_history%clock_gas,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(temperature,godunov_history%clock_temperature,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(internal_energy,godunov_history%clock_internal_energy,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(old,godunov_history%initial_gas,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(ot,godunov_history%initial_temperature,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(held_old,godunov_history%initial_internal_energy,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(ir,godunov_history%initial_reaction_source,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(ir,godunov_history%clock_reaction_source,good,all_levels=.true.)
      if(.not.good) return
      if(present(ancestor_start)) then
      call copy_field(ancestor_start,godunov_history%ancestor_start,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(ancestor_end,godunov_history%ancestor_end,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(temperature_start,godunov_history%ancestor_temperature_start,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(temperature_end,godunov_history%ancestor_temperature_end,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(ancestor_internal_energy_start,godunov_history%ancestor_internal_energy_start,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(ancestor_internal_energy_end,godunov_history%ancestor_internal_energy_end,good,all_levels=.true.)
      if(.not.good) return
        godunov_history%fractions_start=fractions_start;godunov_history%fractions_end=fractions_end
      endif
      if(present(source_start)) then
      call copy_field(source_start,godunov_history%source_start,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(source_end,godunov_history%source_end,good,all_levels=.true.)
      if(.not.good) return
      endif
      good=.true.
    end subroutine begin_godunov_history
    subroutine capture_godunov_source(stage,gas,temp,external,diffusion,auxiliary,faces,good)
      integer, intent(in) :: stage
      type(amr_tree_state_3d), intent(in) :: gas
      type(amr_tree_state_3d), intent(in) :: temp
      type(amr_tree_state_3d), intent(in) :: external,diffusion,auxiliary
      type(eb_integrated_flux_3d), allocatable, intent(in) :: faces(:)
      logical, intent(out) :: good
      godunov_history%sources(stage)%index=stage
      godunov_history%sources(stage)%coarse_registers_accumulated=prepared_redist.and.stage>1
      associate(saved=>godunov_history%sources(stage))
      call copy_field(gas,saved%gas,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(temp,saved%temperature,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(external,saved%external,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(diffusion,saved%rate,good,all_levels=.true.)
      if(.not.good) return
      call history_copy_flux(faces,saved%faces,good)
      if(.not.good) return
      if(present(auxiliary_source)) then
        if(stage==0) then
      call copy_field(auxiliary_packed0,saved%auxiliary,good,all_levels=.true.)
      if(.not.good) return
        else
      call copy_field(auxiliary_packed1,saved%auxiliary,good,all_levels=.true.)
      if(.not.good) return
        endif
      endif
      if(stage==0) then
      call copy_field(held_old,saved%internal_energy,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(held_rate0,saved%scalar_rate,good,all_levels=.true.)
      if(.not.good) return
      call history_copy_flux(held_flux0,saved%scalar_faces,good)
      if(.not.good) return
        if(prepared_redist) then
      call copy_registers(stage_registers0,saved%registers,good)
      if(.not.good) return
      call copy_registers(held_registers0,saved%scalar_registers,good)
      if(.not.good) return
        endif
      else
      call copy_field(held_candidate,saved%internal_energy,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(held_rate1,saved%scalar_rate,good,all_levels=.true.)
      if(.not.good) return
      call history_copy_flux(held_flux1,saved%scalar_faces,good)
      if(.not.good) return
        if(prepared_redist) then
      call copy_registers(stage_registers1,saved%registers,good)
      if(.not.good) return
      call copy_registers(held_registers1,saved%scalar_registers,good)
      if(.not.good) return
        endif
      endif
      end associate
      godunov_history%accepted_source_stages=stage+1
      good=.true.
    end subroutine capture_godunov_source
    subroutine capture_godunov_hydro(stage,good,retain_predictor)
      integer, intent(in) :: stage
      logical, intent(out) :: good
      logical, optional, intent(in) :: retain_predictor
      logical :: keep_forcing
      keep_forcing=.false.
      if(present(retain_predictor)) keep_forcing=retain_predictor
      godunov_history%hydro(stage)%index=stage
      godunov_history%hydro(stage)%coarse_registers_accumulated=prepared_redist.and.stage>1
      associate(saved=>godunov_history%hydro(stage))
      saved%source_funded=source_funded
      if(.not.keep_forcing) then
      call copy_field(rates,saved%predictor,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(hydro,saved%predictor_rate,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(held_hydro,saved%predictor_scalar_rate,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(held_hydro_work,saved%predictor_pressure_work,good,all_levels=.true.)
      if(.not.good) return
      call history_copy_flux(hf,saved%predictor_faces,good)
      if(.not.good) return
      call history_copy_flux(held_hydro_flux,saved%predictor_scalar_faces,good)
      if(.not.good) return
      endif
      call copy_field(hydro,saved%rate,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(held_hydro,saved%scalar_rate,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(held_hydro_work,saved%pressure_work,good,all_levels=.true.)
      if(.not.good) return
      call history_copy_flux(hf,saved%faces,good)
      if(.not.good) return
      call history_copy_flux(held_hydro_flux,saved%scalar_faces,good)
      if(.not.good) return
      if(prepared_redist) then
      call copy_registers(hydro_registers,saved%registers,good)
      if(.not.good) return
      call copy_registers(held_hydro_registers,saved%scalar_registers,good)
      if(.not.good) return
      endif
      end associate
      good=.true.
    end subroutine capture_godunov_hydro
    subroutine finish_godunov_history(good)
      logical, intent(out) :: good
      call copy_field(candidate,godunov_history%final_gas,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(ct,godunov_history%final_temperature,good,all_levels=.true.)
      if(.not.good) return
      if(present(source_funding_impulse)) then
        call copy_field(funding_receipt,godunov_history%source_funding_impulse,good,all_levels=.true.)
        if(.not.good) return
      endif
      call copy_field(held_candidate,godunov_history%final_internal_energy,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(ir,godunov_history%reaction_source,good,all_levels=.true.)
      if(.not.good) return
      call history_copy_flux(hf,godunov_history%faces,good)
      if(.not.good) return
      call history_copy_flux(held_flux,godunov_history%scalar_faces,good)
      if(.not.good) return
      call history_copy_flux(registration,godunov_history%registered_faces,good)
      if(.not.good) return
      if(prepared_redist) then
      call copy_registers(final_registers,godunov_history%registers,good)
      if(.not.good) return
      call copy_registers(held_registers,godunov_history%scalar_registers,good)
      if(.not.good) return
      endif
      allocate(godunov_history%boundary_change(size(boundary)),godunov_history%wall_change(size(wall)), &
        godunov_history%chemical_change(size(chemical)),godunov_history%source_change(size(external)), &
        godunov_history%auxiliary_change(naux),stat=history_status)
      call agree(history_status==0,comm,good);if(.not.good) return
      godunov_history%boundary_change=boundary;godunov_history%wall_change=wall
      godunov_history%chemical_change=chemical;godunov_history%source_change=external
      godunov_history%auxiliary_change=0
      if(present(auxiliary_source)) godunov_history%auxiliary_change=auxiliary
      godunov_history%numerical_energy_change=0
      if(present(energy_correction)) godunov_history%numerical_energy_change=numerical_change
      godunov_history%internal_energy_projection_change=projection_inventory(1)
      godunov_history%internal_energy_closure_change=closure_inventory(1)
      godunov_history%initial_internal_energy_change=initial_scalar_reset
      call copy_field(op,godunov_history%old_predictor,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(np,godunov_history%new_predictor,good,all_levels=.true.)
      if(.not.good) return
      godunov_history%pressure_work_change=pressure_work_change
      godunov_history%minimum_theta=theta
      good=.true.
    end subroutine finish_godunov_history
    subroutine prepared_reactor_adapter(species,reactions,energy,source_species,source_energy,duration,rtol,atol, &
        density,temperature,steps,ok,maximum_steps)
      type(nasa7_species), intent(in) :: species(:)
      type(elementary_reaction), intent(in) :: reactions(:)
      real(dp), intent(in) :: energy,source_species(:),source_energy,duration,rtol,atol
      real(dp), intent(inout) :: density(:),temperature
      integer, intent(out) :: steps
      logical, intent(out) :: ok
      integer, optional, intent(in) :: maximum_steps
      logical :: refused
      call evaluate_classified_reactor(evaluate_reaction,species,reactions,energy,source_species,source_energy, &
        duration,rtol,atol,density,temperature,steps,ok,refused,maximum_steps, &
          evaluate_admission=evaluate_reaction_admission,admission_status=admission_state)
      reactor_refused_local=reactor_refused_local.or.refused
    end subroutine prepared_reactor_adapter

    subroutine run_prepared_reaction(good)
      logical, intent(out) :: good
      logical :: supplied
      reactor_refused_local=.false.
      if(present(internal_energy)) then
        if(associated(prepared_integrator)) then
          call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old,ot,rates,dt,rtol,atol, &
            candidate,ct,ir,comm,good,target_depth,evaluate_reaction=prepared_integrator,frame=frame,heat_release=heat, &
            internal_energy=held_candidate,old_internal_energy=held_old)
        else
          call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old,ot,rates,dt,rtol,atol, &
            candidate,ct,ir,comm,good,target_depth,frame=frame,heat_release=heat, &
            internal_energy=held_candidate,old_internal_energy=held_old)
        endif
      else
      if(associated(prepared_integrator)) then
        call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old,ot,rates,dt,rtol,atol, &
          candidate,ct,ir,comm,good,target_depth,evaluate_reaction=prepared_integrator,frame=frame,heat_release=heat)
      else
        call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old,ot,rates,dt,rtol,atol, &
          candidate,ct,ir,comm,good,target_depth,frame=frame,heat_release=heat)
      endif
      endif
      if(good.and.present(energy_correction)) call recover(good,.true.,2)
      call agree(.not.reactor_refused_local,comm,supplied)
      reactor_refused=.not.supplied
      if(reactor_refused) then
        godunov_recovery_forbidden=.true.
        if(present(provider_refused)) provider_refused=.true.
        good=.false.
      endif
    end subroutine run_prepared_reaction

    subroutine capture_mol_interval(good)
      logical, intent(out) :: good
      integer :: allocation_status,p,face_depth
      logical :: rate_good
      face_depth=selected_depth
      if(present(source_funding_stock)) face_depth=-1
      pending_history%native_spatial_weight=native_mol_weight
      pending_history%replayed=present(replay_interval).or.present(native_source_replay)
      pending_history%redistribution=redistribution_method
      pending_history%prepared_redistribution=prepared_redist
      pending_history%accepted_coarse_stages=prepared_redist.and.mol.and.chemistry.and.present(internal_energy)
      pending_history%spatial_diffusion_is_source=external_diffusion_source
      pending_history%auxiliary_has_internal_energy=held_auxiliary
      pending_history%has_auxiliary_source=present(auxiliary_source)
      pending_history%depth=selected_depth;pending_history%interval=dt
      pending_history%has_internal_energy=present(internal_energy)
      pending_history%has_source_funding=present(source_funding_stock)
      if(pending_history%native_iterations>0.and.present(ancestor_start)) then
        call copy_field(ancestor_start,pending_history%ancestor_start,good,all_levels=.true.)
        if(.not.good) return
        call copy_field(ancestor_end,pending_history%ancestor_end,good,all_levels=.true.)
        if(.not.good) return
        call copy_field(temperature_start,pending_history%ancestor_temperature_start,good,all_levels=.true.)
        if(.not.good) return
        call copy_field(temperature_end,pending_history%ancestor_temperature_end,good,all_levels=.true.)
        if(.not.good) return
        if(present(ancestor_internal_energy_start)) then
          call copy_field(ancestor_internal_energy_start,pending_history%ancestor_internal_energy_start,good,all_levels=.true.)
          if(.not.good) return
          call copy_field(ancestor_internal_energy_end,pending_history%ancestor_internal_energy_end,good,all_levels=.true.)
          if(.not.good) return
        endif
        if(present(source_start)) then
          call copy_field(source_start,pending_history%ancestor_source_start,good,all_levels=.true.)
          if(.not.good) return
          call copy_field(source_end,pending_history%ancestor_source_end,good,all_levels=.true.)
          if(.not.good) return
        endif
        ! The recursive caller supplies a full level vector. Archive only
        ! this clock's active ancestors, including the empty root prefix.
        allocate(pending_history%fractions_start,source=fractions_start(:selected_depth),stat=allocation_status)
        call agree(allocation_status==0,comm,good);if(.not.good) return
        allocate(pending_history%fractions_end,source=fractions_end(:selected_depth),stat=allocation_status)
        call agree(allocation_status==0,comm,good);if(.not.good) return
      endif
      if(present(source_funding_stock)) then
        call copy_field(source_funding_stock,pending_history%source_funding_stock,good,all_levels=.true.)
        if(.not.good) return
        call copy_field(funding_receipt,pending_history%source_funding_impulse,good,all_levels=.true.)
        if(.not.good) return
        call copy_field(old,pending_history%clock_initial,good,all_levels=.true.)
        if(.not.good) return
        call copy_field(held_old,pending_history%clock_internal_energy,good,all_levels=.true.)
        if(.not.good) return
      endif
      if(present(internal_energy)) then
        pending_history%internal_energy_policy=held_policy
        pending_history%initial_energy_correction=initial_reset
        pending_history%numerical_energy_correction=numerical_change
        pending_history%forcing_energy_correction=accepted_predictor_reset
        pending_history%internal_energy_projection_change=projection_inventory(1)
      pending_history%internal_energy_closure_change=closure_inventory(1)
      pending_history%initial_internal_energy_change=initial_scalar_reset
        call copy_field(held_old,pending_history%initial_internal_energy,good)
        if(.not.good) return
        call copy_field(held_candidate,pending_history%local_internal_energy,good)
        if(.not.good) return
        call copy_field(held_rate0,pending_history%old_internal_energy_spatial_rate,good)
        if(.not.good) return
        call copy_field(held_rate0,pending_history%internal_energy_spatial_rate,good)
        if(.not.good) return
        rate_good=.true.
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          if(allocated(admitted_held_rate%fields)) then
            pending_history%internal_energy_spatial_rate%fields(p)%values=admitted_held_rate%fields(p)%values
          else
            pending_history%internal_energy_spatial_rate%fields(p)%values=.5_dp* &
              (held_rate0%fields(p)%values+held_rate1%fields(p)%values)
          endif
          ! held_update advances the admitted mean plus the native loan
          ! rate. Archive that same complete spatial operator with its paired
          ! funded faces; the inactive donor debit stays a separate receipt.
          if(mol_source_funded) then
            call guarded_spatial_add(pending_history%internal_energy_spatial_rate%fields(p)%values, &
              held_hydro%fields(p)%values,rate_good)
            if(.not.rate_good) exit
          endif
        enddo
        call agree(rate_good,comm,good);if(.not.good) return
        call copy_mol_level_flux(tree,owners,rank,selected_depth,held_flux0, &
          pending_history%old_internal_energy_flux,comm,good)
        if(.not.good) return
        call copy_mol_level_flux(tree,owners,rank,face_depth,held_flux, &
          pending_history%internal_energy_flux,comm,good)
        if(.not.good) return
      endif
      allocate(pending_history%chemical_change,source=chemical,stat=allocation_status)
      call agree(allocation_status==0,comm,good)
      if(.not.good) return
      allocate(pending_history%wall_change,source=wall,stat=allocation_status)
      call agree(allocation_status==0,comm,ready)
      if(.not.ready) return
      allocate(pending_history%source_change,source=external,stat=allocation_status)
      call agree(allocation_status==0,comm,ready)
      if(.not.ready) return
      allocate(pending_history%auxiliary_change(naux),stat=allocation_status)
      call agree(allocation_status==0,comm,ready)
      if(.not.ready) return
      pending_history%auxiliary_change=0
      if(present(auxiliary_source)) pending_history%auxiliary_change=auxiliary
      call agree(allocation_status==0,comm,good)
      if(.not.good) return
      call copy_field(old,pending_history%initial,good)
      if(.not.good) return
      call copy_field(ot,pending_history%initial_temperature,good)
      if(.not.good) return
      call copy_field(rates,pending_history%forcing,good)
      if(.not.good) return
      call copy_field(ledger,pending_history%mass_source,good)
      if(.not.good) return
      call copy_field(candidate,pending_history%local_state,good)
      if(.not.good) return
      call copy_field(ct,pending_history%local_temperature,good)
      if(.not.good) return
      call copy_field(ir,pending_history%reaction_source,good)
      if(.not.good) return
      call copy_mol_level_flux(tree,owners,rank,selected_depth,f0,pending_history%old_flux,comm,good)
      if(.not.good) return
      call copy_mol_level_flux(tree,owners,rank,face_depth,hf,pending_history%spatial_flux,comm,good)
      if(.not.good.or.redistribution_method=='none') return
      call copy_field(d0,pending_history%old_spatial_rate,good)
      if(.not.good) return
      call copy_field(d0,pending_history%spatial_rate,good)
      if(.not.good) return
      call copy_field(e0,pending_history%old_volume_source,good)
      if(.not.good) return
      call copy_field(e0,pending_history%volume_source,good)
      if(.not.good) return
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        if(allocated(admitted_mean_rate%fields)) then
          pending_history%spatial_rate%fields(p)%values=admitted_mean_rate%fields(p)%values
        else
          pending_history%spatial_rate%fields(p)%values=.5_dp*(d0%fields(p)%values+d1%fields(p)%values)
        endif
        if(present(auxiliary_source)) pending_history%old_volume_source%fields(p)%values= &
          e0%fields(p)%values+a0%fields(p)%values
        pending_history%volume_source%fields(p)%values=.5_dp*(e0%fields(p)%values+e1%fields(p)%values)
        if(present(energy_correction).and.(.not.chemistry.or.iterations==1)) &
          pending_history%volume_source%fields(p)%values(5,:,:,:)= &
            pending_history%volume_source%fields(p)%values(5,:,:,:)+.5_dp*predictor_reset%fields(p)%values(5,:,:,:)/dt
        if(present(auxiliary_source)) pending_history%volume_source%fields(p)%values= &
          pending_history%volume_source%fields(p)%values+.5_dp*(a0%fields(p)%values+a1%fields(p)%values)
      enddo
      if(redistribution_method=='flux_redist') then
        call copy_field(nonface0,pending_history%non_face_spatial_rate,good)
        if(.not.good) return
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          pending_history%non_face_spatial_rate%fields(p)%values=.5_dp* &
            (nonface0%fields(p)%values+nonface1%fields(p)%values)
        enddo
      endif
      if(prepared_redist) then
        call copy_registers(stage_registers0,pending_history%old_redist_registers,good)
        if(.not.good) return
        call copy_registers(final_registers,pending_history%redist_registers,good)
        if(.not.good) return
        if(present(internal_energy)) then
          call copy_registers(held_registers0,pending_history%old_internal_energy_registers,good)
          if(.not.good) return
          call copy_registers(held_registers,pending_history%internal_energy_registers,good)
          if(.not.good) return
        endif
      elseif(present(redistribution_transfer)) then
        call copy_field(transfer0,pending_history%old_redistribution_transfer,good,all_levels=.true.)
        if(.not.good) return
        call copy_field(final_transfer,pending_history%redistribution_transfer,good,all_levels=.true.)
      endif
    end subroutine capture_mol_interval
    subroutine capture_native_mol_provider(stage,kind,rate,halos,faces,limiter,good)
      integer, intent(in) :: stage,kind
      type(amr_tree_state_3d), intent(in) :: rate,halos
      type(eb_integrated_flux_3d), allocatable, intent(in) :: faces(:)
      real(dp), intent(in) :: limiter
      logical, intent(out) :: good
      logical :: local_ok
      associate(saved=>pending_history%provider_stages(kind,stage))
      call copy_godunov_field_3d(rate,saved%provider_rate,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      if(kind==1) then
        call copy_godunov_field_3d(halos,saved%provider_halos,local_ok)
        call agree(local_ok,comm,good);if(.not.good) return
      endif
      call copy_godunov_fluxes_3d(faces,saved%provider_faces,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      saved%index=stage;saved%provider_kind=kind;saved%provider_limiter=limiter;saved%provider_accepted=.true.
      end associate
    end subroutine capture_native_mol_provider

    subroutine validate_native_mol_replay(good)
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: checked,checked_held
      integer :: stage
      logical :: valid
      valid=all(ieee_is_finite(replay_interval%native_force)).and. &
        all(ieee_is_finite(replay_interval%native_acceleration)).and. &
        ieee_is_finite(replay_interval%internal_energy_policy%eta2).and. &
        ieee_is_finite(replay_interval%interval).and.allocated(replay_interval%spatial_flux).and. &
        allocated(replay_interval%internal_energy_flux)
      call agree(valid,comm,good);if(.not.good) return
      valid=mol.and.chemistry.and.present(internal_energy).and.present(energy_correction).and. &
        selected_depth>=0.and.redistribution_method=='none'.and..not.present(replay_spatial_weight)
      valid=valid.and.(replay_interval%separate_momentum.eqv.selected_momentum)
      valid=valid.and.(replay_interval%consistent_velocity.eqv.selected_velocity)
      valid=valid.and.(replay_interval%paired_pressure.eqv.selected_pressure)
      valid=valid.and.(replay_interval%connected_directions.eqv.selected_connected)
      valid=valid.and.replay_interval%has_internal_energy.and..not.replay_interval%prepared_redistribution.and. &
        replay_interval%redistribution=='none'.and.replay_interval%native_iterations==iterations.and. &
        replay_interval%depth==selected_depth.and.replay_interval%interval==dt.and. &
        all(replay_interval%native_force==force).and.all(replay_interval%native_acceleration==accel).and. &
        (replay_interval%has_native_frame.eqv.present(frame)).and. &
        (replay_interval%has_auxiliary_source.eqv.present(auxiliary_source)).and. &
        (replay_interval%spatial_diffusion_is_source.eqv.external_diffusion_source).and. &
        (replay_interval%auxiliary_has_internal_energy.eqv.held_auxiliary).and. &
        all(replay_interval%native_providers.eqv.[present(evaluate_spatial_diffusion),present(evaluate_diffusion)]).and. &
        replay_interval%native_provider_halo_width==spatial_callback_width
      if(present(frame)) valid=valid.and.replay_interval%native_frame%is_valid().and. &
        replay_interval%native_frame%axis==frame%axis.and. &
        replay_interval%native_frame%angular_speed==frame%angular_speed.and. &
        all(replay_interval%native_frame%axis_point==frame%axis_point)
      valid=valid.and.replay_interval%internal_energy_policy%eta2==held_policy%eta2.and. &
        (replay_interval%internal_energy_policy%allow_small_energy.eqv.held_policy%allow_small_energy).and. &
        (replay_interval%internal_energy_policy%allow_negative_energy.eqv.held_policy%allow_negative_energy).and. &
        (replay_interval%internal_energy_policy%update_total_from_internal.eqv.held_policy%update_total_from_internal)
      valid=valid.and.(allocated(replay_interval%auxiliary_stages).eqv.present(auxiliary_source))
      call agree(valid,comm,good);if(.not.good) return
      if(present(auxiliary_source)) then
        valid=lbound(replay_interval%auxiliary_stages,1)==0.and. &
          ubound(replay_interval%auxiliary_stages,1)==iterations
        call agree(valid,comm,good);if(.not.good) return
        do stage=0,iterations
          valid=valid_held_interval_field(replay_interval%auxiliary_stages(stage),tree,owners,selected_depth,rank,naux)
          call agree(valid,comm,good);if(.not.good) return
        enddo
      endif
      valid=valid_held_interval_field(replay_interval%initial_reaction_source,tree,owners,selected_depth,rank,nv-4)
      call agree(valid,comm,good);if(.not.good) return
      valid=valid_mpi_mol_native_provider_cache_3d(tree,owners,geometry,replay_interval,nv,rank)
      call agree(valid,comm,good);if(.not.good) return
      call reconstruct_mpi_eb_mol_interval_3d(tree,owners,geometry,replay_interval,replay_interval%spatial_flux, &
        nv,checked,comm,good,internal_energy_flux=replay_interval%internal_energy_flux,internal_energy_rate=checked_held)
    end subroutine validate_native_mol_replay

    subroutine validate_prepared_replay(good)
      logical, intent(out) :: good
      integer :: stage,p,d(3)
      logical :: valid
      call agree_held_projection(replay_interval%internal_energy_projection_change,comm,good)
      if(.not.good) return
      call agree_held_projection(replay_interval%internal_energy_closure_change,comm,good)
      if(.not.good) return
      call agree_held_clock_packet([replay_interval%initial_internal_energy_change, &
        replay_interval%held_clock_start_change,replay_interval%held_clock_end_change, &
        replay_interval%held_clock_projection_before,replay_interval%held_clock_reset_before, &
        replay_interval%held_clock_projection_after,replay_interval%held_clock_reset_after, &
        merge(1._dp,0._dp,replay_interval%held_clock_events_closed), &
        merge(1._dp,0._dp,replay_interval%held_clock_prefixes_valid)],comm,good)
      if(.not.good) return
      valid=mol.and.chemistry.and.prepared_redist.and..not.prepared_godunov.and.selected_depth>=0
      valid=valid.and.replay_interval%redistribution=='state_redist'.and.replay_interval%prepared_redistribution
      valid=valid.and..not.replay_interval%separate_momentum.and.replay_interval%native_iterations==0.and. &
        .not.allocated(replay_interval%auxiliary_stages).and..not.allocated(replay_interval%provider_stages).and. &
        .not.allocated(replay_interval%external_stages).and..not.native_external_frozen.and. &
        .not.allocated(replay_interval%native_provider_geometry).and..not.any(replay_interval%native_providers)
      valid=valid.and..not.replay_interval%has_source_funding.and. &
        .not.allocated(replay_interval%source_funding_stock%fields).and. &
        .not.allocated(replay_interval%source_funding_impulse%fields)
      valid=valid.and.replay_interval%depth==selected_depth.and.replay_interval%interval==dt
      valid=valid.and.(replay_interval%has_auxiliary_source.eqv.present(auxiliary_source))
      valid=valid.and.(replay_interval%has_internal_energy.eqv.present(internal_energy))
      valid=valid.and.(replay_interval%accepted_coarse_stages.eqv.present(internal_energy)).and. &
        (replay_interval%spatial_diffusion_is_source.eqv.external_diffusion_source).and. &
        (replay_interval%auxiliary_has_internal_energy.eqv.held_auxiliary)
      if(present(internal_energy)) then
        valid=valid.and.replay_interval%internal_energy_policy%eta2==held_policy%eta2.and. &
          (replay_interval%internal_energy_policy%allow_small_energy.eqv.held_policy%allow_small_energy).and. &
          (replay_interval%internal_energy_policy%allow_negative_energy.eqv.held_policy%allow_negative_energy).and. &
          (replay_interval%internal_energy_policy%update_total_from_internal.eqv.held_policy%update_total_from_internal)
        valid=valid.and.valid_held_interval_field(replay_interval%initial_internal_energy,tree,owners,selected_depth,rank)
      endif
      valid=valid.and.allocated(replay_interval%prepared_stages).and.allocated(replay_interval%initial_reaction_source%fields)
      call agree(valid,comm,good)
      if(.not.good) return
      valid=lbound(replay_interval%prepared_stages,1)==0.and.ubound(replay_interval%prepared_stages,1)==iterations
      valid=valid.and.size(replay_interval%initial_reaction_source%fields)==size(owners)
      call agree(valid,comm,good)
      if(.not.good) return
      do stage=0,iterations
        associate(saved=>replay_interval%prepared_stages(stage))
        valid=saved%source_stage==stage.and.saved%depth==selected_depth.and.saved%interval==dt.and. &
          saved%max_order==redist_order.and.saved%target_volume_fraction==redist_target.and.saved%register_scale==.5_dp
        valid=valid.and.(saved%has_internal_energy.eqv.present(internal_energy))
        if(present(internal_energy)) then
          valid=valid.and.valid_held_interval_field(saved%held_state_halos,tree,owners,selected_depth,rank,padding=6).and. &
            valid_held_interval_field(saved%held_raw_rate_halos,tree,owners,selected_depth,rank,padding=6)
        endif
        valid=valid.and.allocated(saved%external_rate%fields).and.allocated(saved%auxiliary_rate%fields)
        call agree(valid,comm,good)
        if(.not.good) return
        valid=size(saved%external_rate%fields)==size(owners).and.size(saved%auxiliary_rate%fields)==size(owners)
        call agree(valid,comm,good)
        if(.not.good) return
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          d=tree%nodes(p)%shape
          if(.not.allocated(saved%external_rate%fields(p)%values).or. &
            .not.allocated(saved%auxiliary_rate%fields(p)%values).or. &
            .not.allocated(replay_interval%initial_reaction_source%fields(p)%values)) then
            valid=.false.;exit
          endif
          valid=all(shape(saved%external_rate%fields(p)%values)==[nv,d]).and. &
            all(shape(saved%auxiliary_rate%fields(p)%values)==[naux,d]).and. &
            all(shape(replay_interval%initial_reaction_source%fields(p)%values)==[nv-4,d])
          valid=valid.and.all(ieee_is_finite(saved%external_rate%fields(p)%values)).and. &
            all(ieee_is_finite(saved%auxiliary_rate%fields(p)%values)).and. &
            all(ieee_is_finite(replay_interval%initial_reaction_source%fields(p)%values))
          if(.not.valid) exit
        enddo
        call agree(valid,comm,good)
        if(.not.good) return
        end associate
      enddo
    end subroutine validate_prepared_replay

    subroutine capture_prepared_volume_stage(stage,external,auxiliary,good)
      integer, intent(in) :: stage
      type(amr_tree_state_3d), intent(in) :: external,auxiliary
      logical, intent(out) :: good
      integer :: p
      call copy_field(external,pending_history%prepared_stages(stage)%external_rate,good)
      if(.not.good) return
      if(present(auxiliary_source)) then
        call copy_field(auxiliary,pending_history%prepared_stages(stage)%auxiliary_rate,good)
      else
        call copy_field(external,pending_history%prepared_stages(stage)%auxiliary_rate,good)
        if(.not.good) return
        do p=1,size(owners)
          if(owners(p)==rank.and.active(p)) pending_history%prepared_stages(stage)%auxiliary_rate%fields(p)%values=0
        enddo
      endif
    end subroutine capture_prepared_volume_stage

    subroutine capture_prepared_mol_stage(stage,state_halos,raw_halos,spatial,faces,registers,good,low,native_faces)
      integer, intent(in) :: stage
      type(amr_tree_state_3d), intent(in) :: state_halos,raw_halos,spatial
      type(eb_integrated_flux_3d), intent(in) :: faces(:)
      type(state_redist_register_3d), allocatable, intent(in) :: registers(:)
      logical, intent(out) :: good
      integer :: allocation_status
      type(amr_tree_state_3d), optional, intent(in) :: low
      type(eb_integrated_flux_3d), optional, intent(in) :: native_faces(:)
      good=.false.
      if(.not.allocated(pending_history%prepared_stages)) then
        allocate(pending_history%prepared_stages(0:iterations),stat=allocation_status)
        call agree(allocation_status==0,comm,good)
        if(.not.good) return
      endif
      call agree(stage>=0.and.stage<=iterations,comm,good)
      if(.not.good) return
      associate(record=>pending_history%prepared_stages(stage))
        record%source_stage=stage;record%depth=selected_depth
        record%max_order=redist_order;record%interval=dt
        record%target_volume_fraction=redist_target;record%register_scale=.5_dp
        record%spatial_weight=spatial_weight
        if(spatial_weight<1) then
          call agree(present(low).and.present(native_faces),comm,good)
          if(.not.good) return
          call copy_field(low,record%low_spatial_rate,good)
          if(.not.good) return
          call copy_mol_level_flux(tree,owners,rank,selected_depth,native_faces,record%native_faces,comm,good)
          if(.not.good) return
          call copy_mol_level_flux(tree,owners,rank,selected_depth,prepared_low_faces,record%low_faces,comm,good)
          if(.not.good) return
        endif
        call copy_field(state_halos,record%state_halos,good)
        if(.not.good) return
        call copy_field(raw_halos,record%raw_rate_halos,good)
        if(.not.good) return
        call copy_field(spatial,record%spatial_rate,good)
        if(.not.good) return
        call copy_mol_level_flux(tree,owners,rank,selected_depth,faces,record%faces,comm,good)
        if(.not.good) return
        call copy_registers(registers,record%registers,good)
      end associate
    end subroutine capture_prepared_mol_stage

    subroutine copy_field(source,destination,good,all_levels,unit_bounds)
      type(amr_tree_state_3d), intent(in) :: source
      type(amr_tree_state_3d), intent(out) :: destination
      logical, intent(out) :: good
      logical, optional, intent(in) :: all_levels,unit_bounds
      integer :: p,allocation_status,dims(4)
      logical :: valid,copy_all,normalize
      copy_all=.false.
      if(present(all_levels)) copy_all=all_levels
      normalize=.false.
      if(present(unit_bounds)) normalize=unit_bounds
      allocate(destination%fields(size(tree%nodes)),stat=allocation_status)
      call agree(allocation_status==0,comm,good)
      if(.not.good) return
      valid=.true.
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        if(.not.copy_all.and..not.active(p)) cycle
        if(normalize) then
          dims=shape(source%fields(p)%values)
          allocate(destination%fields(p)%values(dims(1),dims(2),dims(3),dims(4)), &
            source=source%fields(p)%values,stat=allocation_status)
        else
          allocate(destination%fields(p)%values,source=source%fields(p)%values,stat=allocation_status)
        endif
        if(allocation_status/=0) then
          valid=.false.;exit
        endif
      enddo
      call agree(valid,comm,good)
    end subroutine copy_field
    subroutine copy_registers(source,destination,good)
      type(state_redist_register_3d), allocatable, intent(in) :: source(:)
      type(state_redist_register_3d), allocatable, intent(out) :: destination(:)
      logical, intent(out) :: good
      integer :: p,allocation_status
      logical :: valid
      call agree(allocated(source),comm,good)
      if(.not.good) return
      call agree(size(source)==size(owners),comm,good)
      if(.not.good) return
      allocate(destination(size(owners)),stat=allocation_status)
      call agree(allocation_status==0,comm,good)
      if(.not.good) return
      valid=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        destination(p)%scale=source(p)%scale
        if(allocated(source(p)%coarse_flags)) then
          allocate(destination(p)%coarse_flags,source=source(p)%coarse_flags,stat=allocation_status)
          if(allocation_status/=0) then
            valid=.false.;exit
          endif
        endif
        if(allocated(source(p)%fine_ghost)) then
          allocate(destination(p)%fine_ghost,source=source(p)%fine_ghost,stat=allocation_status)
          if(allocation_status/=0) then
            valid=.false.;exit
          endif
        endif
        if(allocated(source(p)%valid_cells)) then
          allocate(destination(p)%valid_cells,source=source(p)%valid_cells,stat=allocation_status)
          if(allocation_status/=0) then
            valid=.false.;exit
          endif
        endif
        if(allocated(source(p)%coarse_delta)) then
          allocate(destination(p)%coarse_delta,source=source(p)%coarse_delta,stat=allocation_status)
          if(allocation_status/=0) then
            valid=.false.;exit
          endif
        endif
        if(allocated(source(p)%fine_delta)) then
          allocate(destination(p)%fine_delta,source=source(p)%fine_delta,stat=allocation_status)
          if(allocation_status/=0) then
            valid=.false.;exit
          endif
        endif
      enddo
      call agree(valid,comm,good)
    end subroutine copy_registers
    subroutine redist_boundary_samples(gas,temp,good)
      type(amr_tree_state_3d), intent(in) :: gas,temp
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: metrics,gt,vf
      integer :: p,dims(3),allocation_status
      logical :: local_ok,all_ok
      good=.false.
      if(allocated(redist_halos)) deallocate(redist_halos)
      if(allocated(redist_geometry)) deallocate(redist_geometry)
      allocate(redist_halos,redist_geometry(size(owners)),metrics%fields(size(owners)),stat=allocation_status)
      call agree(allocation_status==0,comm,all_ok)
      if(.not.all_ok) return
      local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        dims=tree%nodes(p)%shape
        call extract_eb_geometry_3d(extended_geometry(p),[5,5,5],dims+8,redist_geometry(p),local_ok)
        if(.not.local_ok) exit
        allocate(metrics%fields(p)%values(1,-1:dims(1)+2,-1:dims(2)+2,-1:dims(3)+2),stat=allocation_status)
        local_ok=allocation_status==0
        if(.not.local_ok) exit
        metrics%fields(p)%values(1,:,:,:)=redist_geometry(p)%volume_fraction
      enddo
      call agree(local_ok,comm,all_ok)
      if(.not.all_ok) return
      call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,gas,temp,redist_halos,gt,comm,all_ok, &
        inflow_state=inflow_state,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature, &
        inflow_pressure=inflow_pressure,boundary_faces=boundary_faces,ghost_volume_fraction=vf, &
        halo_width=2,target_volume_fraction=metrics,restrict_snapshot=.false.,frame=frame)
      good=all_ok
    end subroutine redist_boundary_samples
    logical function active(p)
      integer, intent(in) :: p
      active=selected_depth<0.or.tree%nodes(p)%level==selected_depth
    end function
    subroutine initialize_held_godunov_stage(gas,temp,stage,good,fraction)
      type(amr_tree_state_3d), intent(in) :: gas,temp
      integer, intent(in) :: stage
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fraction(:)
      type(amr_tree_state_3d) :: stored,rate,old_scalar,zero_source,old_gas,old_temp,volumes,metrics
      type(eb_integrated_flux_3d), allocatable :: faces(:)
      type(state_redist_register_3d), allocatable :: registers(:)
      real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),vf(:,:,:),ghost_held(:,:,:)
      integer :: p,d(3),level,first,last,allocation_status
      logical :: ready_stage,local_ok
      ! Pinned diffusion has no direct UEINT flux/source. Its ML
      ! redistribution still applies to old UEINT when diffusion is active.
      good=.false.
      if(stage==0) then
        stored=held_old
      else
        stored=held_candidate
      endif
      rate=stored
      allocate(faces(size(owners)),stat=allocation_status)
      call agree(allocation_status==0,comm,ready_stage);if(.not.ready_stage) return
      local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        rate%fields(p)%values=0
        if(.not.active(p)) cycle
        d=tree%nodes(p)%shape
        allocate(faces(p)%x(1,0:d(1),d(2),d(3)),faces(p)%y(1,d(1),0:d(2),d(3)), &
          faces(p)%z(1,d(1),d(2),0:d(3)),stat=allocation_status)
        if(allocation_status/=0) then
          local_ok=.false.;exit
        endif
        faces(p)%x=0;faces(p)%y=0;faces(p)%z=0
      enddo
      call agree(local_ok,comm,ready_stage);if(.not.ready_stage) return
      call initialize_inactive_faces(faces,1,ready_stage)
      if(.not.ready_stage) return
      if(prepared_godunov) then
        allocate(registers(size(owners)),stat=allocation_status)
        call agree(allocation_status==0,comm,ready_stage);if(.not.ready_stage) return
        if(diffusive.or.present(evaluate_spatial_diffusion)) then
                    allocate(metrics%fields(size(owners)),old_scalar%fields(size(owners)), &
            zero_source%fields(size(owners)),stat=allocation_status)
          call agree(allocation_status==0,comm,ready_stage);if(.not.ready_stage) return
          local_ok=.true.
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            d=tree%nodes(p)%shape
            allocate(metrics%fields(p)%values(1,1-6:d(1)+6,1-6:d(2)+6,1-6:d(3)+6),stat=allocation_status)
            if(allocation_status/=0) then
              local_ok=.false.;exit
            endif
            metrics%fields(p)%values(1,:,:,:)=register_geometry(p)%volume_fraction
          enddo
          call agree(local_ok,comm,ready_stage);if(.not.ready_stage) return
          call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,gas,temp,old_gas,old_temp,comm,ready_stage, &
            boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
            inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=volumes, &
            halo_width=6,target_volume_fraction=metrics,restrict_snapshot=.false., &
            ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
            temperature_end=temperature_end,ancestor_fraction=fraction,internal_energy=stored, &
            ghost_internal_energy=old_scalar,ancestor_internal_energy_start=ancestor_internal_energy_start, &
            ancestor_internal_energy_end=ancestor_internal_energy_end)
          if(.not.ready_stage) return
          local_ok=.true.
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            if(.not.active(p)) then
              deallocate(old_scalar%fields(p)%values);cycle
            endif
            d=tree%nodes(p)%shape
            allocate(ghost_held(d(1)+12,d(2)+12,d(3)+12), &
              zero_source%fields(p)%values(1,d(1)+12,d(2)+12,d(3)+12),source=0._dp,stat=allocation_status)
            if(allocation_status/=0) then
              local_ok=.false.;exit
            endif
            ghost_held=old_scalar%fields(p)%values(1,:,:,:)
            deallocate(old_scalar%fields(p)%values)
            allocate(old_scalar%fields(p)%values(1,d(1)+12,d(2)+12,d(3)+12),stat=allocation_status)
            if(allocation_status/=0) then
              local_ok=.false.;exit
            endif
            old_scalar%fields(p)%values(1,:,:,:)=ghost_held;deallocate(ghost_held)
          enddo
          call agree(local_ok,comm,ready_stage);if(.not.ready_stage) return
          first=max(0,selected_depth);last=first
          if(selected_depth<0) last=tree%level_count()-1
          do level=first,last
            call compute_state_redist_mpi_level_stage_3d(tree,owners,register_geometry,old_scalar,zero_source,1, &
              level,dt,rate,registers,comm,ready_stage,boundary_faces,redist_target,redist_order,register_scale=1._dp)
            if(.not.ready_stage) return
          enddo
        else
          local_ok=.true.
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            if(.not.active(p)) cycle
            call initialize_state_redist_tree_register_3d(tree,p,1,6,registers(p),ready_stage,boundary_faces,1._dp)
            local_ok=local_ok.and.ready_stage
          enddo
          call agree(local_ok,comm,ready_stage);if(.not.ready_stage) return
        endif
      endif
      if(stage==0) then
        held_rate0=rate;call move_alloc(faces,held_flux0)
        if(prepared_godunov) call move_alloc(registers,held_registers0)
      else
        if(prepared_godunov) then
          if(stage>1) then
            local_ok=.true.
            do p=1,size(owners)
              if(owners(p)/=rank) cycle
              if(.not.active(p)) cycle
              d=tree%nodes(p)%shape
              call accumulate_state_redist_coarse_register_3d(registers(p),held_registers1(p),1,d+12,ready_stage)
              local_ok=local_ok.and.ready_stage
            enddo
            call agree(local_ok,comm,ready_stage);if(.not.ready_stage) return
          endif
          call move_alloc(registers,held_registers1)
        endif
        held_rate1=rate;call move_alloc(faces,held_flux1)
      endif
      good=.true.
    end subroutine initialize_held_godunov_stage
    subroutine blend_native_mol_stage(rate,faces,scalar_rate,scalar_faces,good)
      type(amr_tree_state_3d), intent(inout) :: rate,scalar_rate
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: faces(:),scalar_faces(:)
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: ignored,ignored_held,old_halos
      integer :: p
      logical :: local_ok
      good=.true.
      if(native_mol_weight==1) return
      if(.not.allocated(native_mol_low_faces)) then
        call compute_mol_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,extended_geometry,old,ot,dt,1, &
          ignored,native_mol_low_faces,comm,good,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          lax_friedrichs=.true.,frame=frame,internal_energy_rate=ignored_held, &
          internal_energy_flux=native_mol_low_held_faces,internal_energy_state=held_old, &
          internal_energy_old_state_halos=old_halos, &
          target_depth=target_depth,ancestor_start=ancestor_start,ancestor_end=ancestor_end, &
          temperature_start=temperature_start,temperature_end=temperature_end,ancestor_fraction=fractions_start, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected)
        if(.not.good) return
      endif
      local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        call blend_native_mol_face_pair_3d(geometry(p),dt,native_mol_weight, &
          native_mol_low_faces(p),native_mol_low_held_faces(p),rate%fields(p)%values,faces(p), &
          scalar_rate%fields(p)%values,scalar_faces(p),local_ok, &
          mass_species_from_faces=.not.present(evaluate_spatial_diffusion).and..not.present(evaluate_diffusion))
        if(.not.local_ok) exit
      enddo
      call agree(local_ok,comm,good)
      if(good) theta=min(theta,native_mol_weight)
    end subroutine blend_native_mol_stage

    subroutine held_mol_source(gas,temp,stored,rate,scalar_flux,scalar_registers,good,fraction, &
        gas_rate,gas_flux,gas_spatial_halos,gas_old_halos,scalar_state_halos,scalar_raw_halos)
      type(amr_tree_state_3d), intent(in) :: gas,temp,stored
      type(amr_tree_state_3d), intent(inout) :: rate
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: scalar_flux(:)
      type(state_redist_register_3d), allocatable, intent(inout) :: scalar_registers(:)
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fraction(:)
      type(amr_tree_state_3d), intent(inout) :: gas_rate
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: gas_flux(:)
      type(amr_tree_state_3d), optional, intent(inout) :: gas_spatial_halos,gas_old_halos
      type(amr_tree_state_3d), optional, intent(inout) :: scalar_state_halos,scalar_raw_halos
      type(amr_tree_state_3d) :: old_halos,spatial_halos
      integer :: depth,first_depth,last_depth
      call compute_mol_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,extended_geometry,gas,temp,dt,2, &
        gas_rate,gas_flux,comm,good,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        lax_friedrichs=mol_lax_friedrichs,frame=frame,internal_energy_rate=rate, &
        internal_energy_flux=scalar_flux,internal_energy_state=stored,internal_energy_old_state_halos=old_halos, &
        internal_energy_spatial_rate_halos=spatial_halos, &
        spatial_rate_halos=gas_spatial_halos,old_state_halos=gas_old_halos, &
        target_depth=target_depth,ancestor_start=ancestor_start,ancestor_end=ancestor_end, &
        temperature_start=temperature_start,temperature_end=temperature_end,ancestor_fraction=fraction, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected)
      if(.not.good) return
      if(present(scalar_state_halos)) scalar_state_halos=old_halos
      if(present(scalar_raw_halos)) scalar_raw_halos=spatial_halos
      call initialize_inactive_faces(scalar_flux,1,good)
      if(.not.good) return
      if(.not.prepared_redist) return
      first_depth=max(0,selected_depth);last_depth=first_depth
      if(selected_depth<0) last_depth=tree%level_count()-1
      do depth=first_depth,last_depth
        call compute_state_redist_mpi_level_stage_3d(tree,owners,extended_geometry,old_halos,spatial_halos,1, &
          depth,dt,rate,scalar_registers,comm,good,boundary_faces,redist_target,redist_order,register_scale=.5_dp)
        if(.not.good) return
      enddo
    end subroutine held_mol_source
    ! Admit the actual nonreacting mean source before its EOS reset. A
    ! corrector can first exhaust stock at its new-time source evaluation.
    subroutine restore_cached_funding_hydro(good)
      logical, intent(out) :: good
      logical :: matches,local_ok
      integer :: p
      good=.true.
      if(mol.or.prepared_redist.or.godunov_partial) return
      if(.not.present(godunov_replay_interval)) return
      if(present(spatial_flux_weight)) then
        if(spatial_flux_weight/=1) return
      endif
      associate(saved=>godunov_replay_interval)
      if(saved%source_prefix_only.or.saved%face_weight/=1) return
      ! A raw replay of the same clock must retain the accepted conservative
      ! borrowing family. A supplied reactor must not receive its unfunded
      ! precursor merely because it does not run the built-in preview. Fresh
      ! clocks, changed predictors and weighted trials still rebuild funding.
      matches=same_godunov_predictor_3d(saved%hydro(it),godunov_predictor_forcing, &
        godunov_predictor_hydro,godunov_predictor_held,godunov_predictor_work, &
        godunov_predictor_faces,godunov_predictor_scalar_faces)
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        matches=matches.and.all(old%fields(p)%values==saved%initial_gas%fields(p)%values)
        matches=matches.and.all(ot%fields(p)%values(1,:,:,:)== &
          saved%initial_temperature%fields(p)%values(1,:,:,:))
        matches=matches.and.all(held_old%fields(p)%values==saved%initial_internal_energy%fields(p)%values)
      enddo
      call agree(matches,comm,local_ok)
      if(.not.local_ok) return
      call copy_field(saved%hydro(it)%rate,hydro,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(saved%hydro(it)%scalar_rate,held_hydro,good,all_levels=.true.)
      if(.not.good) return
      call copy_field(saved%hydro(it)%pressure_work,held_hydro_work,good,all_levels=.true.)
      if(.not.good) return
      call copy_godunov_fluxes_3d(saved%hydro(it)%faces,hf,local_ok)
      call agree(local_ok,comm,good)
      if(.not.good) return
      call copy_godunov_fluxes_3d(saved%hydro(it)%scalar_faces,held_hydro_flux,local_ok)
      call agree(local_ok,comm,good)
      if(.not.good) return
      source_funded=saved%hydro(it)%source_funded
      if(godunov_stage_capture) call capture_godunov_hydro(it,good,retain_predictor=.true.)
      end associate
    end subroutine restore_cached_funding_hydro

    subroutine admit_held_native_corrector(good)
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: trial_h
      real(dp) :: post(reactive_nvar(size(species))),qp(reactive_nprim(size(species))),tg,sound,bias
      integer :: p,d(3),x,y,z
      logical :: unsafe,cell_ok,global_unsafe
      good=.false.;unsafe=.false.
      call restore_cached_funding_hydro(good)
      if(.not.good) return
      call held_update(held_old,held_rate0,held_rate1,.5_dp,trial_h,cell_ok)
      if(.not.cell_ok) return
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          ! Virtual coarse sources own no physical inventory. Their reacting
          ! forcing is reconstructed from the admitted fine stage below.
          if(chemistry.and.covered(tree,p,[x,y,z])) cycle
          post=old%fields(p)%values(:,x,y,z)+ &
            dt*(rates%fields(p)%values(:,x,y,z)+hydro%fields(p)%values(:,x,y,z))
          cell_ok=all(ieee_is_finite(post)).and.post(1)>0
          if(.not.chemistry) cell_ok=cell_ok.and.all(post(6:)>=0)
          ! Signed species forcing can become physical through reaction.
          ! A caloric preview requires positive composition, while the private
          ! built-in trial retains signed forcing. Chemistry cannot fix rho.
          if(cell_ok.and.chemistry.and..not.present(evaluate_reaction)) then
            bias=0
            call compute_reactive_eb_reactor_bias_3d(geometry(p),[x,y,z],old%fields(p)%values(:,x,y,z), &
              held_old%fields(p)%values(1,x,y,z),bias,cell_ok,frame)
            if(cell_ok) then
              if(all(post(6:)>=0)) then
              call admit_reactive_eb_reactor_energy_3d(species,geometry(p),[x,y,z],post,bias, &
                ot%fields(p)%values(1,x,y,z),cell_ok,frame)
              else
                ! Signed forced species need the actual built-in trial too:
                ! chemistry may replenish them, but an unrecoverable deficit
                ! must request the same conservative material funding.
                cell_ok=.false.
              endif
              ! Composition can recover through reaction even when its forced
              ! caloric preview refuses. Only a refused private built-in trial
              ! asks for funding. A caller-owned reactor is attempted once in
              ! the ordinary native path, never retried or preempted here.
              if(.not.cell_ok) call try_default_reactor_cell(p,[x,y,z],cell_ok)
            endif
          endif
          if(cell_ok.and.chemistry.and.present(evaluate_reaction_admission)) &
            call try_default_reactor_cell(p,[x,y,z],cell_ok)
          ! Previously funded material retains its independently conserved H.
          ! Ordinary unfunded thermal resets keep the pinned policy path.
          if(cell_ok.and.source_funded.and..not.chemistry) then
            call recover_reactive_eb_cell_3d(species,geometry(p),[x,y,z],post, &
              ot%fields(p)%values(1,x,y,z),qp,tg,sound,cell_ok,frame,trial_h%fields(p)%values(1,x,y,z))
          endif
          unsafe=unsafe.or..not.cell_ok
        enddo
        enddo
        enddo
      enddo
      call MPI_Allreduce(admission_state,admission_global,1,MPI_INTEGER,MPI_MAX,comm,admission_ierr)
      if(admission_ierr/=MPI_SUCCESS) return
      admission_state=admission_global
      if(admission_state==reaction_admission_provider_refused) then
        good=.false.;reactor_refused=.true.;godunov_recovery_forbidden=.true.;return
      endif
      call MPI_Allreduce(unsafe,global_unsafe,1,MPI_LOGICAL,MPI_LOR,comm,status)
      if(status/=MPI_SUCCESS) return
      unsafe=global_unsafe
      good=.true.
      if(.not.unsafe) return
      call fund_held_native_stage(good,corrector=.true.)
      if(.not.good) return
      ! New-time funding changes the accepted hydro/face family, while the
      ! original Godunov predictor forcing remains its actual native input.
      if(godunov_stage_capture) call capture_godunov_hydro(it,good,retain_predictor=.true.)
    end subroutine admit_held_native_corrector

    subroutine prepare_covered_reactor_forcing(good)
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: endpoint
      integer :: p,d(3),x,y,z
      logical :: reconstruct
      good=.false.;reconstruct=source_funded
      endpoint=rates
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        endpoint%fields(p)%values=old%fields(p)%values+dt*rates%fields(p)%values
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(.not.covered(tree,p,[x,y,z])) cycle
          ! A restored funded hydro family has no mutable funding flag. A
          ! deficient virtual density still selects the same reconstruction.
          reconstruct=reconstruct.or.endpoint%fields(p)%values(1,x,y,z)<=0
        enddo
        enddo
        enddo
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,reconstruct,1,MPI_LOGICAL,MPI_LOR,comm,status)
      if(status/=MPI_SUCCESS) return
      good=.true.
      if(.not.reconstruct) return
      ! Restrict the forced gas vector, not a caloric preview: reaction may
      ! make its composition admissible. The actual old UEINT clock remains
      ! the reactor's authority, and physical leaf rates are never rewritten.
      call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,endpoint,nv,comm,good)
      if(.not.good) return
      do p=1,size(owners)
        if(.not.active(p)) cycle
        if(owners(p)/=rank) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(.not.covered(tree,p,[x,y,z])) cycle
          rates%fields(p)%values(:,x,y,z)=(endpoint%fields(p)%values(:,x,y,z)- &
            old%fields(p)%values(:,x,y,z))/dt
        enddo
        enddo
        enddo
      enddo
    end subroutine prepare_covered_reactor_forcing

    ! A private built-in trial preserves chemistry's actual admission. It
    ! publishes no state/history and invokes no caller-owned reactor provider.
    subroutine try_default_reactor_cell(node,cell,good,prospective_forcing)
      integer, intent(in) :: node,cell(3)
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: prospective_forcing(:)
      real(dp) :: initial(nv,1,1,1),forcing(nv,1,1,1),trial(nv,1,1,1),source(nv-4,1,1,1)
      real(dp) :: initial_t(1,1,1),trial_t(1,1,1),old_h(1,1,1),trial_h(1,1,1),position(3)
      real(dp) :: old_minus_projection,new_minus_projection,energy_change
      logical :: finite_projection
      real(dp), allocatable :: potential(:,:,:)
      initial(:,1,1,1)=old%fields(node)%values(:,cell(1),cell(2),cell(3))
      forcing(:,1,1,1)=rates%fields(node)%values(:,cell(1),cell(2),cell(3))+ &
        hydro%fields(node)%values(:,cell(1),cell(2),cell(3))
      good=.false.
      if(present(prospective_forcing)) then
        if(size(prospective_forcing)/=nv) return
        if(.not.all(ieee_is_finite(prospective_forcing))) return
        forcing(:,1,1,1)=prospective_forcing
      endif
      ! The native cell kernel squares momentum before division. Do not let
      ! an inconclusive extreme-value preview create a floating-point trap.
      if(maxval(abs(initial(2:4,1,1,1)))>sqrt(huge(1._dp)/3)) return
      trial=initial+dt*forcing
      if(.not.all(ieee_is_finite(trial))) return
      if(maxval(abs(trial(2:4,1,1,1)))>sqrt(huge(1._dp)/3)) return
      old_minus_projection=0;new_minus_projection=0
      call compute_reactive_eb_reactor_bias_3d(geometry(node),cell,initial(:,1,1,1),0._dp, &
        old_minus_projection,finite_projection,frame)
      if(.not.finite_projection) return
      call compute_reactive_eb_reactor_bias_3d(geometry(node),cell,trial(:,1,1,1),0._dp, &
        new_minus_projection,finite_projection,frame)
      if(.not.finite_projection) return
      if(new_minus_projection<0) then
        if(old_minus_projection>huge(1._dp)+new_minus_projection) return
      else if(new_minus_projection>0) then
        if(old_minus_projection< -huge(1._dp)+new_minus_projection) return
      endif
      energy_change=old_minus_projection-new_minus_projection
      if(dt<1) then
        if(abs(energy_change)>huge(1._dp)*dt) return
      endif
      initial_t=ot%fields(node)%values(1,cell(1),cell(2),cell(3));trial_t=initial_t
      old_h=held_old%fields(node)%values(1,cell(1),cell(2),cell(3));trial_h=old_h
      if(present(frame)) then
        allocate(potential(1,1,1))
        associate(g=>geometry(node))
        position=[g%x_lower,g%y_lower,g%z_lower]+(real(cell,dp)-.5_dp)*[g%dx,g%dy,g%dz]
        end associate
        call rotating_potential_3d(frame,position,potential(1,1,1),good)
        if(.not.good) return
      endif
      if(present(evaluate_reaction_admission)) then
        if(admission_state==reaction_admission_provider_refused) return
        probe_admission=reaction_admission_unknown
        call react_sdc_state_3d(species,reactions,initial,initial_t,forcing,dt,rtol,atol,trial,trial_t,source,good, &
          evaluate_reaction=admission_probe_adapter,rotating_potential=potential, &
          internal_energy_density=trial_h,old_internal_energy_density=old_h)
        good=probe_admission==reaction_admission_accepted
        admission_state=max(admission_state,probe_admission)
      else
      call react_sdc_state_3d(species,reactions,initial,initial_t,forcing,dt,rtol,atol,trial,trial_t,source,good, &
        rotating_potential=potential,internal_energy_density=trial_h,old_internal_energy_density=old_h)
      endif
    end subroutine try_default_reactor_cell

    subroutine admission_probe_adapter(table,mechanism,energy,source_species,source_energy,duration,rtol,atol, &
        density,temperature,steps,ok,maximum_steps)
      type(nasa7_species), intent(in) :: table(:)
      type(elementary_reaction), intent(in) :: mechanism(:)
      real(dp), intent(in) :: energy,source_species(:),source_energy,duration,rtol,atol
      real(dp), intent(inout) :: density(:),temperature
      integer, intent(out) :: steps
      logical, intent(out) :: ok
      integer, optional, intent(in) :: maximum_steps
      integer(int64) :: packet(2*size(density)+6)
      integer :: entry,allocation_status
      logical :: cached
      ! Memoize only immutable private proposals, never the final actual
      ! reactor recheck. Exact bits include temperature and both tolerances;
      ! mechanism arrays are immutable for this owning public call.
      packet=transfer([energy,source_species,source_energy,duration,rtol,atol,density,temperature],packet)
      cached=.false.
      do entry=1,admission_cache_count
        if(.not.all(packet==admission_cache_packets(:,entry))) cycle
        probe_admission=admission_cache_results(entry);cached=.true.;exit
      enddo
      if(.not.cached) then
        call evaluate_reaction_admission(table,mechanism,energy,source_species,source_energy,duration,rtol,atol, &
          density,temperature,probe_admission)
        if(probe_admission<reaction_admission_accepted.or.probe_admission>reaction_admission_provider_refused) &
          probe_admission=reaction_admission_unknown
        if(.not.allocated(admission_cache_packets)) then
          allocate(admission_cache_packets(size(packet),16),stat=allocation_status)
        endif
        if(allocated(admission_cache_packets)) then
          admission_cache_packets(:,admission_cache_next)=packet
          admission_cache_results(admission_cache_next)=probe_admission
          admission_cache_count=min(16,admission_cache_count+1)
          admission_cache_next=mod(admission_cache_next,16)+1
        endif
      endif
      ! Packing only: never unpack a fictitious integrated state.
      steps=0;ok=.false.
    end subroutine admission_probe_adapter

    ! Native selected hydro leaves inactive scalar faces unallocated.
    ! A donor transaction can touch those clocks, so represent their absent
    ! native operator with zero faces before the paired inventory/flux check.
    subroutine ensure_selected_scalar_faces(good)
      logical, intent(out) :: good
      call initialize_inactive_faces(hf,nv,good)
      if(.not.good) return
      call initialize_inactive_faces(held_hydro_flux,1,good)
      if(.not.good) return
      ! Mean-flux admission owns selected clocks only. Represent absent
      ! inactive owner faces by zero before copying the admitted family;
      ! existing donor faces and every nonowner packet remain untouched.
      if(allocated(admitted_held_faces)) &
        call initialize_inactive_faces(admitted_held_faces,1,good)
    end subroutine ensure_selected_scalar_faces

    subroutine initialize_inactive_faces(faces,components,good)
      type(eb_integrated_flux_3d), intent(inout) :: faces(:)
      integer, intent(in) :: components
      logical, intent(out) :: good
      integer :: p,d(3),allocation_status
      logical :: local_ok
      good=.true.
      if(selected_depth<0.or..not.native_held_funding) return
      local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        if(active(p)) cycle
        d=tree%nodes(p)%shape
        if(.not.allocated(faces(p)%x)) then
          allocate(faces(p)%x(components,0:d(1),d(2),d(3)),source=0._dp,stat=allocation_status)
          local_ok=local_ok.and.allocation_status==0
        endif
        if(.not.allocated(faces(p)%y)) then
          allocate(faces(p)%y(components,d(1),0:d(2),d(3)),source=0._dp,stat=allocation_status)
          local_ok=local_ok.and.allocation_status==0
        endif
        if(.not.allocated(faces(p)%z)) then
          allocate(faces(p)%z(components,d(1),d(2),0:d(3)),source=0._dp,stat=allocation_status)
          local_ok=local_ok.and.allocation_status==0
        endif
      enddo
      call agree(local_ok,comm,good)
    end subroutine initialize_inactive_faces

    subroutine admit_native_funding_cell(node,cell,pre,require_post,good)
      integer, intent(in) :: node,cell(3)
      real(dp), intent(in) :: pre(:)
      logical, intent(in) :: require_post
      logical, intent(out) :: good
      real(dp) :: delivered(reactive_nvar(size(species))),next_rate(reactive_nvar(size(species)))
      real(dp) :: next_h,stored,bias,pv(reactive_nprim(size(species))),temp,sound
      integer :: x,y,z
      logical :: cell_ok
      cell_ok=.true.
      x=cell(1);y=cell(2);z=cell(3)
      if(owners(node)==rank.and.require_post.and..not.covered(tree,node,cell)) then
        ! Preserve the native rate/impulse association exactly. A packet
        ! accepted by pre+debt can fail composition closure after cancellation
        ! in this expression, so validate it while the proposal is private.
        next_rate=hydro%fields(node)%values(:,x,y,z)+(pre(:nv)-old%fields(node)%values(:,x,y,z))/dt
        next_h=held_hydro%fields(node)%values(1,x,y,z)+ &
          (pre(nv+1)-held_old%fields(node)%values(1,x,y,z))/dt
        cell_ok=all(ieee_is_finite(next_rate)).and.ieee_is_finite(next_h)
        if(cell_ok) then
          delivered=old%fields(node)%values(:,x,y,z)+dt*(rates%fields(node)%values(:,x,y,z)+next_rate)
          stored=funding_raw_h%fields(node)%values(1,x,y,z)+ &
            dt*(next_h-held_hydro%fields(node)%values(1,x,y,z))
          if(present(source_funding_stock).and..not.active(node)) then
            delivered=delivered+source_funding_stock%fields(node)%values(:nv,x,y,z)
            stored=stored+source_funding_stock%fields(node)%values(nv+1,x,y,z)
          endif
          call recover_reactive_eb_cell_3d(species,geometry(node),cell,delivered, &
            ot%fields(node)%values(1,x,y,z),pv,temp,sound,cell_ok,frame,stored)
          if(cell_ok.and.chemistry) then
            call compute_reactive_eb_reactor_bias_3d(geometry(node),cell,old%fields(node)%values(:,x,y,z), &
              held_old%fields(node)%values(1,x,y,z),bias,cell_ok,frame)
            if(cell_ok) call admit_reactive_eb_reactor_energy_3d(species,geometry(node),cell,delivered, &
              bias,ot%fields(node)%values(1,x,y,z),cell_ok,frame)
            ! A static gas/H endpoint can leave an active donor unable to
            ! complete its actual forced reaction. Admit the same native rate
            ! with a private built-in trial for the accepted corrector family.
            ! A predictor supplies only an unreacted spatial/source estimate.
            ! Inactive clocks carry deferred receipts; caller-owned reactors
            ! retain their ordinary single attempt and are never previewed.
            if(cell_ok.and.funding_corrector_stage.and.active(node).and. &
              (.not.present(evaluate_reaction).or.present(evaluate_reaction_admission))) &
              call try_default_reactor_cell(node,cell,cell_ok,rates%fields(node)%values(:,x,y,z)+next_rate)
          endif
        endif
      endif
      call MPI_Allreduce(admission_state,admission_global,1,MPI_INTEGER,MPI_MAX,comm,admission_ierr)
      if(admission_ierr/=MPI_SUCCESS) then
        good=.false.;return
      endif
      admission_state=admission_global
      if(admission_state==reaction_admission_provider_refused) then
        cell_ok=.false.;reactor_refused=.true.;godunov_recovery_forbidden=.true.
      endif
      call agree(cell_ok,comm,good)
    end subroutine admit_native_funding_cell

    subroutine fund_held_native_stage(good,corrector)
      logical, intent(out) :: good
      logical, optional, intent(in) :: corrector
      type(amr_tree_state_3d), allocatable :: reactor_bias
      logical :: final_stage
      type(amr_tree_state_3d) :: funded,temps,gas_debt,held_debt,funded_h,raw_h,next_hydro,next_held
      type(eb_integrated_flux_3d), allocatable :: gas_faces(:),scalar_faces(:),gas_stage_faces(:),held_stage_faces(:)
      real(dp) :: post(reactive_nvar(size(species))),stored,pv(reactive_nprim(size(species))),temp,sound
      integer :: p,d(3),x,y,z,allocation_status
      logical :: cell_ok,all_ok
      final_stage=.false.
      if(present(corrector)) final_stage=corrector
      call ensure_selected_scalar_faces(good)
      if(.not.good) return
      good=.false.;funded=old;temps=ot;funded_h=held_old;gas_faces=hf;scalar_faces=held_hydro_flux
      if(mol.or.diffusive) then
        gas_stage_faces=f0;held_stage_faces=held_flux0
        if(final_stage) then
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            gas_stage_faces(p)%x=.5_dp*(f0(p)%x+f1(p)%x)
            gas_stage_faces(p)%y=.5_dp*(f0(p)%y+f1(p)%y)
            gas_stage_faces(p)%z=.5_dp*(f0(p)%z+f1(p)%z)
            held_stage_faces(p)%x=.5_dp*(held_flux0(p)%x+held_flux1(p)%x)
            held_stage_faces(p)%y=.5_dp*(held_flux0(p)%y+held_flux1(p)%y)
            held_stage_faces(p)%z=.5_dp*(held_flux0(p)%z+held_flux1(p)%z)
          enddo
        endif
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          gas_faces(p)%x=gas_faces(p)%x+gas_stage_faces(p)%x
          gas_faces(p)%y=gas_faces(p)%y+gas_stage_faces(p)%y
          gas_faces(p)%z=gas_faces(p)%z+gas_stage_faces(p)%z
          scalar_faces(p)%x=scalar_faces(p)%x+held_stage_faces(p)%x
          scalar_faces(p)%y=scalar_faces(p)%y+held_stage_faces(p)%y
          scalar_faces(p)%z=scalar_faces(p)%z+held_stage_faces(p)%z
        enddo
      endif
      if(final_stage) then
        call held_update(held_old,held_rate0,held_rate1,.5_dp,raw_h,cell_ok)
      else
      call held_update(held_old,held_rate0,held_rate0,.5_dp,raw_h,cell_ok)
      endif
      if(.not.cell_ok) return
      if(selected_depth>=0) then
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          if(active(p)) cycle
          raw_h%fields(p)%values=held_old%fields(p)%values+dt*held_hydro%fields(p)%values
        enddo
      endif
      allocate(gas_debt%fields(size(owners)),held_debt%fields(size(owners)))
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        gas_debt%fields(p)%values=dt*(rates%fields(p)%values+hydro%fields(p)%values)
        held_debt%fields(p)%values=raw_h%fields(p)%values-held_old%fields(p)%values
        if(present(source_funding_stock).and..not.active(p)) then
          gas_debt%fields(p)%values=gas_debt%fields(p)%values+source_funding_stock%fields(p)%values(:nv,:,:,:)
          held_debt%fields(p)%values(1,:,:,:)=held_debt%fields(p)%values(1,:,:,:)+ &
            source_funding_stock%fields(p)%values(nv+1,:,:,:)
        endif
        ! Covered operator values are virtual; physical fine sources already
        ! carry their own debt. Re-reflux must not charge those values twice.
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(.not.covered(tree,p,[x,y,z])) cycle
          gas_debt%fields(p)%values(:,x,y,z)=0
          held_debt%fields(p)%values(:,x,y,z)=0
        enddo
        enddo
        enddo
      enddo
      if(chemistry) then
        allocate(reactor_bias,stat=allocation_status)
        call agree(allocation_status==0,comm,good)
        if(.not.good) return
        all_ok=.true.
        reactor_bias=held_old
        do p=1,size(owners)
        if(owners(p)/=rank) cycle
          reactor_bias%fields(p)%values=0
          d=tree%nodes(p)%shape
          do z=1,d(3)
          do y=1,d(2)
          do x=1,d(1)
            if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
            call compute_reactive_eb_reactor_bias_3d(geometry(p),[x,y,z],old%fields(p)%values(:,x,y,z), &
              held_old%fields(p)%values(1,x,y,z),reactor_bias%fields(p)%values(1,x,y,z),cell_ok,frame)
            all_ok=all_ok.and.cell_ok
          enddo
          enddo
          enddo
        enddo
      call agree(all_ok,comm,good)
      if(.not.good) return
      endif
      funding_raw_h=raw_h
      funding_corrector_stage=final_stage
      call repair_mpi_eb_tree_held_interface_budget_3d(species,tree,owners,geometry,funded,temps,funded_h, &
        gas_faces,scalar_faces,comm,cell_ok,0,boundary_faces,frame,gas_debt,held_debt, &
        reactor_energy_bias=reactor_bias,native_projection=admit_native_funding_cell, &
        check_native_endpoints=present(evaluate_reaction_admission).and.final_stage)
      if(.not.cell_ok) return
      next_hydro=hydro;next_held=held_hydro;all_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        next_hydro%fields(p)%values=hydro%fields(p)%values+ &
          (funded%fields(p)%values-old%fields(p)%values)/dt
        next_held%fields(p)%values=held_hydro%fields(p)%values+ &
          (funded_h%fields(p)%values-held_old%fields(p)%values)/dt
        all_ok=all_ok.and.all(ieee_is_finite(next_hydro%fields(p)%values)).and. &
          all(ieee_is_finite(next_held%fields(p)%values))
        if(.not.all_ok) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0.or. &
            covered(tree,p,[x,y,z])) cycle
          post=old%fields(p)%values(:,x,y,z)+dt*(rates%fields(p)%values(:,x,y,z)+ &
            next_hydro%fields(p)%values(:,x,y,z))
          stored=raw_h%fields(p)%values(1,x,y,z)+ &
            dt*(next_held%fields(p)%values(1,x,y,z)-held_hydro%fields(p)%values(1,x,y,z))
          if(present(source_funding_stock).and..not.active(p)) then
            post=post+source_funding_stock%fields(p)%values(:nv,x,y,z)
            stored=stored+source_funding_stock%fields(p)%values(nv+1,x,y,z)
          endif
          call recover_reactive_eb_cell_3d(species,geometry(p),[x,y,z],post, &
            ot%fields(p)%values(1,x,y,z),pv,temp,sound,cell_ok,frame,stored)
          if(cell_ok.and.chemistry) call admit_reactive_eb_reactor_energy_3d(species,geometry(p),[x,y,z],post, &
            reactor_bias%fields(p)%values(1,x,y,z),ot%fields(p)%values(1,x,y,z),cell_ok,frame)
          if(cell_ok.and.chemistry.and.final_stage.and.active(p).and.present(evaluate_reaction_admission)) &
            call try_default_reactor_cell(p,[x,y,z],cell_ok, &
              rates%fields(p)%values(:,x,y,z)+next_hydro%fields(p)%values(:,x,y,z))
          all_ok=all_ok.and.cell_ok
        enddo
        enddo
        enddo
      enddo
      call agree(all_ok,comm,good)
      if(.not.good) return
      if(mol.or.diffusive) then
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          gas_faces(p)%x=gas_faces(p)%x-gas_stage_faces(p)%x
          gas_faces(p)%y=gas_faces(p)%y-gas_stage_faces(p)%y
          gas_faces(p)%z=gas_faces(p)%z-gas_stage_faces(p)%z
          scalar_faces(p)%x=scalar_faces(p)%x-held_stage_faces(p)%x
          scalar_faces(p)%y=scalar_faces(p)%y-held_stage_faces(p)%y
          scalar_faces(p)%z=scalar_faces(p)%z-held_stage_faces(p)%z
        enddo
        if(mol) mol_source_funded=.true.
      endif
      source_funded=.true.
      hydro=next_hydro;held_hydro=next_held;hf=gas_faces;held_hydro_flux=scalar_faces
    end subroutine fund_held_native_stage

    subroutine measure_held_projection(good)
      logical, intent(out) :: good
      integer :: p,x,y,z,d(3)
      logical :: finite,cell_ok
      call held_update(held_old,held_rate0,held_rate1,.5_dp,projection_transport,good)
      if(.not.good) return
      projection_delta=held_old;finite=.true.
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        projection_delta%fields(p)%values=0
        if(.not.active(p)) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          if(covered(tree,p,[x,y,z])) cycle
          call compute_reactive_eb_held_projection(internal_energy%fields(p)%values(1,x,y,z), &
            held_old%fields(p)%values(1,x,y,z),projection_transport%fields(p)%values(1,x,y,z), &
            held_candidate%fields(p)%values(1,x,y,z),projection_delta%fields(p)%values(1,x,y,z),cell_ok)
          finite=finite.and.cell_ok
        enddo
        enddo
        enddo
      enddo
      call agree(finite,comm,good)
      if(.not.good) return
      call integrate_mpi_eb_tree_3d(tree,owners,geometry,projection_delta,projection_inventory,comm,good)
    end subroutine measure_held_projection

    subroutine held_update(initial,first,second,weight,output,good)
      type(amr_tree_state_3d), intent(in) :: initial,first,second
      real(dp), intent(in) :: weight
      type(amr_tree_state_3d), intent(inout) :: output
      logical, intent(out) :: good
      integer :: p
      logical :: finite
      output=initial;finite=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        if(allocated(admitted_held_rate%fields)) then
          output%fields(p)%values=initial%fields(p)%values+dt*admitted_held_rate%fields(p)%values
        else
        output%fields(p)%values=initial%fields(p)%values+weight*dt* &
          (first%fields(p)%values+second%fields(p)%values)
        endif
        if(.not.mol.or.mol_source_funded) &
          output%fields(p)%values=output%fields(p)%values+dt*held_hydro%fields(p)%values
        finite=finite.and.all(ieee_is_finite(output%fields(p)%values))
      enddo
      call agree(finite,comm,good)
    end subroutine held_update
    subroutine predictor_field(ext,diff,aux,output,good)
      type(amr_tree_state_3d), intent(in) :: ext,diff,aux
      type(amr_tree_state_3d), intent(out) :: output
      logical, intent(out) :: good
      logical :: finite
      integer :: p
      output=ext;finite=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        output%fields(p)%values=ext%fields(p)%values+diff%fields(p)%values
        if(present(auxiliary_source)) output%fields(p)%values=output%fields(p)%values+aux%fields(p)%values
        output%fields(p)%values(6:,:,:,:)=output%fields(p)%values(6:,:,:,:)+ir%fields(p)%values(:nv-5,:,:,:)
        output%fields(p)%values(5,:,:,:)=output%fields(p)%values(5,:,:,:)+ir%fields(p)%values(nv-4,:,:,:)
        if(.not.all(ieee_is_finite(output%fields(p)%values))) finite=.false.
      end do
      call agree(finite,comm,good)
    end subroutine
    ! A first MOL pass can withdraw an absent species through reconstructed
    ! new-time transport. Later Picard passes can transport mass introduced
    ! by a new auxiliary family that is no longer part of the accepted RHS.
    ! Limit the combined spatial flux against a first-order old-state
    ! baseline and this pass's complete volume-source budget. This safeguard
    ! is additional to pinned MOL and leaves admissible passes unchanged.
    subroutine limit_mol_pass(good)
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: fixed_rate,new_rate,low_rate
      type(eb_integrated_flux_3d), allocatable :: baseline(:),new_flux(:)
      type(amr_tree_state_3d) :: low_held_rate,new_held_rate,old_held_halos
      type(eb_integrated_flux_3d), allocatable :: low_held_faces(:),new_held_faces(:)
      real(dp) :: endpoint(nv)
      integer :: p,d(3),x,y,z,ios
      logical :: safe,local_safe,periodic(3)
      good=.false.;local_safe=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          if(source_funded.and.covered(tree,p,[x,y,z])) cycle
          if(mol_source_funded) then
            ! Funding admits this represented native state. A separate
            ! addition of opposing rates can manufacture a negative trace
            ! and discard an already admissible conservative loan family.
            endpoint=old%fields(p)%values(:,x,y,z)+dt* &
              (rates%fields(p)%values(:,x,y,z)+hydro%fields(p)%values(:,x,y,z))
          else
            endpoint=old%fields(p)%values(:,x,y,z)+dt*rates%fields(p)%values(:,x,y,z)
          endif
          if(redistribution_method/='flux_redist') &
            endpoint(6:)=endpoint(6:)+dt*ir%fields(p)%values(:nv-5,x,y,z)
          local_safe=local_safe.and.all(ieee_is_finite(endpoint)).and.endpoint(1)>0.and.all(endpoint(6:)>=0)
        enddo
        enddo
        enddo
      enddo
      if(redistribution_method=='flux_redist'.and.present(redistribution_transfer)) then
        do p=1,size(owners)
          if(owners(p)/=rank.or.active(p)) cycle
          d=tree%nodes(p)%shape
          do z=1,d(3)
          do y=1,d(2)
          do x=1,d(1)
            if(geometry(p)%volume_fraction(x,y,z)<=0.or.covered(tree,p,[x,y,z])) cycle
            endpoint=old%fields(p)%values(:,x,y,z)+ &
              .5_dp*(transfer0%fields(p)%values(:,x,y,z)+transfer1%fields(p)%values(:,x,y,z))
            local_safe=local_safe.and.all(ieee_is_finite(endpoint)).and.endpoint(1)>0.and.all(endpoint(6:)>=0)
          enddo
          enddo
          enddo
        enddo
      endif
      call agree(local_safe,comm,safe)
      if(safe) then
        good=.true.;return
      endif
      ! Preserve the caller's nonunit forest family after paired admission.
      ! A rejected mean cannot choose an independent local blend weight.
      if(native_mol_weight/=1) return
      if(redistribution_method=='flux_redist') then
        call limit_mol_flux_redist_pass(good)
        return
      endif
      fixed_rate=e0;new_rate=d1
      if(present(internal_energy)) then
      call compute_mol_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,extended_geometry,old,ot, &
        dt,1,low_rate,baseline,comm,safe,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start,mol_lax_friedrichs,frame, &
        internal_energy_rate=low_held_rate,internal_energy_flux=low_held_faces,internal_energy_state=held_old, &
        internal_energy_old_state_halos=old_held_halos, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected)
      else
      call compute_mol_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,mol_geometry,old,ot, &
        dt,1,low_rate,baseline,comm,safe,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start,mol_lax_friedrichs,frame)
      endif
      if(.not.safe) return
      allocate(new_flux(size(owners)),stat=ios)
      call agree(ios==0,comm,safe);if(.not.safe) return
      if(present(internal_energy)) then
        new_held_rate=held_rate1
        allocate(new_held_faces(size(owners)),stat=ios)
        call agree(ios==0,comm,safe);if(.not.safe) return
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          if(.not.active(p)) cycle
          new_held_rate%fields(p)%values=.5_dp*(held_rate0%fields(p)%values+held_rate1%fields(p)%values)
          new_held_faces(p)=held_flux1(p)
          new_held_faces(p)%x=.5_dp*(held_flux0(p)%x+held_flux1(p)%x)
          new_held_faces(p)%y=.5_dp*(held_flux0(p)%y+held_flux1(p)%y)
          new_held_faces(p)%z=.5_dp*(held_flux0(p)%z+held_flux1(p)%z)
        enddo
      endif
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        fixed_rate%fields(p)%values=.5_dp*(e0%fields(p)%values+e1%fields(p)%values)
        new_rate%fields(p)%values=.5_dp*(d0%fields(p)%values+d1%fields(p)%values)
        if(present(auxiliary_source)) fixed_rate%fields(p)%values=fixed_rate%fields(p)%values+ &
          .5_dp*(a0%fields(p)%values+a1%fields(p)%values)
        if(mol_source_funded) fixed_rate%fields(p)%values=fixed_rate%fields(p)%values+hydro%fields(p)%values
        fixed_rate%fields(p)%values(6:,:,:,:)=fixed_rate%fields(p)%values(6:,:,:,:)+ &
          ir%fields(p)%values(:nv-5,:,:,:)
        fixed_rate%fields(p)%values(5,:,:,:)=fixed_rate%fields(p)%values(5,:,:,:)+ &
          ir%fields(p)%values(nv-4,:,:,:)
        new_flux(p)=f1(p)
        new_flux(p)%x=.5_dp*(f0(p)%x+f1(p)%x)
        new_flux(p)%y=.5_dp*(f0(p)%y+f1(p)%y)
        new_flux(p)%z=.5_dp*(f0(p)%z+f1(p)%z)
      enddo
      ! I_R predicts chemistry; it is not a committed volume source. A
      ! preceding provisional inflow can leave a negative chemical estimate
      ! that the physical old-state anchor cannot fund. In that case use the
      ! actual accepted external family as the low-order budget and let the
      ! native reactor recompute chemistry after conservative flux limiting.
      local_safe=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          endpoint=old%fields(p)%values(:,x,y,z)+dt*(low_rate%fields(p)%values(:,x,y,z)+fixed_rate%fields(p)%values(:,x,y,z))
          local_safe=local_safe.and.all(ieee_is_finite(endpoint)).and.endpoint(1)>0.and.all(endpoint(6:)>=0)
        enddo
        enddo
        enddo
      enddo
      call agree(local_safe,comm,safe)
      if(.not.safe) then
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          fixed_rate%fields(p)%values=.5_dp*(e0%fields(p)%values+e1%fields(p)%values)
          if(present(auxiliary_source)) fixed_rate%fields(p)%values=fixed_rate%fields(p)%values+ &
            .5_dp*(a0%fields(p)%values+a1%fields(p)%values)
          if(mol_source_funded) fixed_rate%fields(p)%values=fixed_rate%fields(p)%values+hydro%fields(p)%values
        enddo
      endif
      periodic=[boundary_faces(1)=='periodic',boundary_faces(3)=='periodic',boundary_faces(5)=='periodic']
      if(present(internal_energy)) then
      call apply_mpi_eb_flux_correction_3d(tree,owners,geometry,old,fixed_rate,baseline,dt,new_rate,new_flux, &
        comm,safe,periodic,target_depth=selected_depth, &
        internal_energy_baseline=low_held_faces,internal_energy_rate=new_held_rate,internal_energy_flux=new_held_faces)
        if(.not.safe) return
        call move_alloc(new_held_faces,admitted_held_faces)
        call copy_field(new_held_rate,admitted_held_rate,safe)
        if(.not.safe) return
      else
      call apply_mpi_eb_flux_correction_3d(tree,owners,geometry,old,fixed_rate,baseline,dt,new_rate,new_flux, &
        comm,safe,periodic,target_depth=selected_depth)
      endif
      if(.not.safe) return
      ! The correction admits the mean itself. Retain that mean rather than
      ! recovering it from an encoded endpoint after trace-budget cancellation.
      call copy_field(new_rate,admitted_mean_rate,safe)
      if(.not.safe) return
      call copy_mol_level_flux(tree,owners,rank,selected_depth,new_flux,admitted_mean_faces,comm,safe)
      if(.not.safe) return
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        new_flux(p)%x=2*new_flux(p)%x-f0(p)%x
        new_flux(p)%y=2*new_flux(p)%y-f0(p)%y
        new_flux(p)%z=2*new_flux(p)%z-f0(p)%z
        d1%fields(p)%values=2*new_rate%fields(p)%values-d0%fields(p)%values
        f1(p)=new_flux(p)
      enddo
      good=.true.
    end subroutine limit_mol_pass
    ! Face-only FCT cannot retain an unfunded redistribution residual.
    ! Blend its complete linear operator and off-level impulse with the same
    ! old spatial anchor; all conserved components use one admitted weight.
    subroutine capture_nonface_rate(source,faces,stage,good)
      type(amr_tree_state_3d), intent(in) :: source
      type(eb_integrated_flux_3d), intent(in) :: faces(:)
      integer, intent(in) :: stage
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: residual
      real(dp) :: divergence(nv),volume
      integer :: p,x,y,z,c
      call copy_field(source,residual,good,all_levels=.true.)
      if(.not.good) return
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        associate(g=>geometry(p),f=>faces(p))
        do z=1,g%nz
        do y=1,g%ny
        do x=1,g%nx
          if(g%volume_fraction(x,y,z)<=0) cycle
          volume=g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz
          divergence=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(:,x-1,y,z)- &
            g%x_face_fraction(x,y,z)*f%x(:,x,y,z))+ &
            g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(:,x,y-1,z)- &
            g%y_face_fraction(x,y,z)*f%y(:,x,y,z))+ &
            g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(:,x,y,z-1)- &
            g%z_face_fraction(x,y,z)*f%z(:,x,y,z))
          do c=1,nv
            if(.not.has_diffusion.and.(c==1.or.c>=6)) then
              residual%fields(p)%values(c,x,y,z)=0
            else
              residual%fields(p)%values(c,x,y,z)=source%fields(p)%values(c,x,y,z)-divergence(c)/volume/dt
            endif
          enddo
        enddo
        enddo
        enddo
        end associate
      enddo
      if(stage==0) then
        call move_alloc(residual%fields,nonface0%fields)
      else
        call move_alloc(residual%fields,nonface1%fields)
      endif
    end subroutine capture_nonface_rate

    subroutine reaction_mass_pools(pools,good)
      integer, intent(out) :: pools(:)
      logical, intent(out) :: good
      integer :: r,c,first,other
      good=.false.
      pools=[(c,c=1,size(pools))]
      do r=1,size(reactions)
        if(.not.allocated(reactions(r)%reactant_stoich).or..not.allocated(reactions(r)%product_stoich)) return
        if(size(reactions(r)%reactant_stoich)/=size(pools).or.size(reactions(r)%product_stoich)/=size(pools)) return
        first=0
        do c=1,size(pools)
          if(reactions(r)%reactant_stoich(c)==reactions(r)%product_stoich(c)) cycle
          if(first==0) then
            first=pools(c)
          else
            other=pools(c)
            where(pools==other) pools=first
          endif
        enddo
      enddo
      good=.true.
    end subroutine reaction_mass_pools

    subroutine limit_prepared_mol_pass(last_stage,good,reacted,force_search)
      integer, intent(in) :: last_stage
      logical, intent(out) :: good,reacted
      logical, optional, intent(in) :: force_search
      type(mpi_eb_mol_state_redist_stage_3d), allocatable :: rebuilt(:)
      type(amr_tree_state_3d) :: native,low_rate,best,best_t,best_ir,held_native,held_low,best_held
      type(amr_tree_state_3d), allocatable :: best_heat
      type(state_redist_register_3d), allocatable :: registers(:),scalar_registers(:)
      real(dp) :: low(nv),high(nv),fixed(nv),blend,weight,original_weight,lower,upper,middle,best_final_reset
      integer :: pools(nv-5)
      integer :: p,x,y,z,c,stage,depth,first_depth,last_depth,d(3),attempt
      logical :: valid,safe,unsafe,correct,forced
      good=.false.;reacted=.false.;valid=.true.;unsafe=.false.
      forced=.false.
      if(present(force_search)) forced=force_search
      original_weight=spatial_weight
      call reaction_mass_pools(pools,valid)
      call agree(valid,comm,safe)
      if(.not.safe) return
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          high=old%fields(p)%values(:,x,y,z)+dt*rates%fields(p)%values(:,x,y,z)
          valid=valid.and.all(ieee_is_finite(high))
          unsafe=unsafe.or.high(1)<=0
          ! Reactions conserve every connected changing-species mass pool.
          ! Inert/unchanged third bodies stay separate. Individual chemical
          ! intermediates may be supplied by the native reaction, whereas a
          ! negative complete pool cannot be repaired by that reaction.
          do c=1,nv-5
            if(pools(c)/=c) cycle
            unsafe=unsafe.or.sum(high(6:),mask=pools==c)<0
          enddo
        enddo
        enddo
        enddo
      enddo
      call agree(valid,comm,safe)
      if(.not.safe) return
      call MPI_Allreduce(unsafe,correct,1,MPI_LOGICAL,MPI_LOR,comm,status)
      if(status/=MPI_SUCCESS) return
      if(.not.correct.and..not.forced) then
        good=.true.;return
      endif
      ! Rebuild each actual nonlinear native stage, including intermediate
      ! coarse contributions. A changed weight never scales a frozen final
      ! register residual or extrapolates complete source ghosts.
      allocate(rebuilt(0:last_stage),stat=ios)
      call agree(ios==0,comm,safe)
      if(.not.safe) return
      rebuilt=pending_history%prepared_stages(0:last_stage)
      native=d0
      first_depth=max(0,selected_depth);last_depth=first_depth
      if(selected_depth<0) last_depth=tree%level_count()-1
      do stage=0,last_stage
        if(allocated(registers)) deallocate(registers)
        do depth=first_depth,last_depth
          call compute_state_redist_mpi_level_stage_3d(tree,owners,extended_geometry,rebuilt(stage)%state_halos, &
            rebuilt(stage)%raw_rate_halos,nv,depth,dt,native,registers,comm,safe,boundary_faces, &
            redist_target,redist_order,register_scale=.5_dp)
          if(.not.safe) return
        enddo
        valid=.true.
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          call post_eb_spatial_density_rate_3d(geometry(p),native%fields(p)%values,valid)
          if(.not.valid) exit
        enddo
        call agree(valid,comm,safe)
        if(.not.safe) return
        if(rebuilt(stage)%spatial_weight==1) rebuilt(stage)%native_faces=rebuilt(stage)%faces
        call build_prepared_low_rate(rebuilt(stage)%raw_rate_halos,native,rebuilt(stage)%native_faces,low_rate,safe)
        if(.not.safe) return
        call copy_field(native,rebuilt(stage)%spatial_rate,safe)
        if(.not.safe) return
        call copy_field(low_rate,rebuilt(stage)%low_spatial_rate,safe)
        if(.not.safe) return
        call copy_mol_level_flux(tree,owners,rank,selected_depth,prepared_low_faces,rebuilt(stage)%low_faces,comm,safe)
        if(.not.safe) return
        call copy_registers(registers,rebuilt(stage)%registers,safe)
        if(.not.safe) return
        if(present(internal_energy)) then
          held_native=held_rate0
          if(allocated(scalar_registers)) deallocate(scalar_registers)
          do depth=first_depth,last_depth
            call compute_state_redist_mpi_level_stage_3d(tree,owners,extended_geometry,rebuilt(stage)%held_state_halos, &
              rebuilt(stage)%held_raw_rate_halos,1,depth,dt,held_native,scalar_registers,comm,safe,boundary_faces, &
              redist_target,redist_order,register_scale=.5_dp)
            if(.not.safe) return
          enddo
          if(rebuilt(stage)%spatial_weight==1) rebuilt(stage)%held_native_faces=rebuilt(stage)%held_faces
          call build_prepared_low_held_rate(rebuilt(stage)%held_raw_rate_halos,held_native, &
            rebuilt(stage)%held_native_faces,held_low,safe)
          if(.not.safe) return
          call copy_field(held_native,rebuilt(stage)%held_spatial_rate,safe)
          if(.not.safe) return
          call copy_field(held_low,rebuilt(stage)%held_low_spatial_rate,safe)
          if(.not.safe) return
          call copy_mol_level_flux(tree,owners,rank,selected_depth,prepared_low_held_faces, &
            rebuilt(stage)%held_low_faces,comm,safe)
          if(.not.safe) return
          call copy_registers(scalar_registers,rebuilt(stage)%held_registers,safe)
          if(.not.safe) return
        endif
      enddo
      valid=.true.;blend=spatial_weight
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          fixed=.5_dp*(e0%fields(p)%values(:,x,y,z)+e1%fields(p)%values(:,x,y,z))
          if(present(auxiliary_source)) fixed=fixed+.5_dp* &
            (a0%fields(p)%values(:,x,y,z)+a1%fields(p)%values(:,x,y,z))
          low=old%fields(p)%values(:,x,y,z)+dt*(fixed+.5_dp* &
            (rebuilt(0)%low_spatial_rate%fields(p)%values(:,x,y,z)+ &
             rebuilt(last_stage)%low_spatial_rate%fields(p)%values(:,x,y,z)))
          high=old%fields(p)%values(:,x,y,z)+dt*(fixed+.5_dp* &
            (rebuilt(0)%spatial_rate%fields(p)%values(:,x,y,z)+ &
             rebuilt(last_stage)%spatial_rate%fields(p)%values(:,x,y,z)))
          valid=valid.and.all(ieee_is_finite(low)).and.all(ieee_is_finite(high))
          valid=valid.and.low(1)>0.and.all(low(6:)>=0)
          if(.not.valid) cycle
          do c=1,nv
            if(c>=2.and.c<=5) cycle
            if(high(c)<0.or.(c==1.and.high(c)==0)) then
              weight=low(c)/(low(c)-high(c))
              if(weight>0) weight=weight*(1-64*epsilon(1._dp))
              blend=min(blend,weight)
            endif
          enddo
        enddo
        enddo
        enddo
      enddo
      call agree(valid,comm,safe)
      if(.not.safe) return
      call MPI_Allreduce(MPI_IN_PLACE,blend,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call apply_prepared_mol_family(rebuilt,blend,safe)
      if(.not.safe) return
      if(blend==0) then
        theta=min(theta,blend);good=.true.;return
      endif
      if(.not.(forced.and.blend==original_weight)) then
        call combined_rates(.false.)
        call run_prepared_reaction(safe)
        if(reactor_refused) return
        if(safe) then
          theta=min(theta,blend);reacted=.true.;good=.true.;return
        endif
      endif
      ! The native chemical gate rejected this conservative source. Admit
      ! the complete low family first, then retain the largest bracketed weight.
      call apply_prepared_mol_family(rebuilt,0._dp,safe)
      if(.not.safe) return
      call combined_rates(.false.)
      call run_prepared_reaction(safe)
      if(.not.safe) return
      best=candidate;best_t=ct;best_ir=ir
      if(present(internal_energy)) then
        best_held=held_candidate;best_final_reset=final_reset
      endif
      if(present(heat_release)) then
        allocate(best_heat,source=heat,stat=ios)
        call agree(ios==0,comm,safe)
        if(.not.safe) return
      endif
      lower=0;upper=blend
      do attempt=1,12
        middle=.5_dp*(lower+upper)
        call apply_prepared_mol_family(rebuilt,middle,safe)
        if(.not.safe) return
        call combined_rates(.false.)
        call run_prepared_reaction(safe)
        if(reactor_refused) return
        if(safe) then
          lower=middle;best=candidate;best_t=ct;best_ir=ir
          if(present(internal_energy)) then
            best_held=held_candidate;best_final_reset=final_reset
          endif
          if(present(heat_release)) best_heat=heat
        else
          upper=middle
        endif
      enddo
      call apply_prepared_mol_family(rebuilt,lower,safe)
      if(.not.safe) return
      candidate=best;ct=best_t;ir=best_ir
      if(present(internal_energy)) then
        held_candidate=best_held;final_reset=best_final_reset
      endif
      if(present(heat_release)) heat=best_heat
      theta=min(theta,lower);reacted=.true.;good=.true.
    end subroutine limit_prepared_mol_pass

    subroutine apply_prepared_mol_family(original,blend,good)
      type(mpi_eb_mol_state_redist_stage_3d), intent(in) :: original(0:)
      real(dp), intent(in) :: blend
      logical, intent(out) :: good
      type(mpi_eb_mol_state_redist_stage_3d), allocatable :: rebuilt(:)
      type(state_redist_register_3d), allocatable :: registers(:),accumulated(:),scalar_registers(:),held_accumulated(:)
      integer :: last_stage,stage,p,d(3),allocation_status
      logical :: valid,safe
      good=.false.;last_stage=ubound(original,1)
      allocate(rebuilt(0:last_stage),stat=allocation_status)
      call agree(allocation_status==0,comm,safe)
      if(.not.safe) return
      rebuilt=original
      do stage=0,last_stage
        ! Only conserved spatial families change; accepted volume sources and
        ! the native chemical solver retain their original meanings and gates.
        call blend_mol_state_redist_stage_3d(tree,owners,geometry,selected_depth,nv,blend, &
          rebuilt(stage)%low_spatial_rate,rebuilt(stage)%spatial_rate,rebuilt(stage)%registers,comm,safe)
        if(.not.safe) return
        if(present(internal_energy)) then
          call blend_mol_state_redist_stage_3d(tree,owners,geometry,selected_depth,1,blend, &
            rebuilt(stage)%held_low_spatial_rate,rebuilt(stage)%held_spatial_rate,rebuilt(stage)%held_registers,comm,safe)
          if(.not.safe) return
          do p=1,size(owners)
            if(owners(p)/=rank.or..not.active(p)) cycle
            if(blend==0) then
              rebuilt(stage)%held_faces(p)=rebuilt(stage)%held_low_faces(p)
            else
              rebuilt(stage)%held_faces(p)%x=(1-blend)*rebuilt(stage)%held_low_faces(p)%x+ &
                blend*rebuilt(stage)%held_native_faces(p)%x
              rebuilt(stage)%held_faces(p)%y=(1-blend)*rebuilt(stage)%held_low_faces(p)%y+ &
                blend*rebuilt(stage)%held_native_faces(p)%y
              rebuilt(stage)%held_faces(p)%z=(1-blend)*rebuilt(stage)%held_low_faces(p)%z+ &
                blend*rebuilt(stage)%held_native_faces(p)%z
            endif
          enddo
          if(stage==0) then
            held_registers0=rebuilt(stage)%held_registers
          else
            scalar_registers=rebuilt(stage)%held_registers
            if(stage>1.and..not.(chemistry.and.present(internal_energy))) then
              valid=.true.
              do p=1,size(owners)
                if(owners(p)/=rank.or..not.active(p)) cycle
                d=tree%nodes(p)%shape
                call accumulate_state_redist_coarse_register_3d(scalar_registers(p),held_accumulated(p),1,d+12,safe)
                valid=valid.and.safe
              enddo
              call agree(valid,comm,safe)
              if(.not.safe) return
            endif
            call move_alloc(scalar_registers,held_accumulated)
          endif
        endif
        rebuilt(stage)%spatial_weight=blend
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          if(blend==0) then
            rebuilt(stage)%faces(p)=rebuilt(stage)%low_faces(p)
          else
            rebuilt(stage)%faces(p)%x=(1-blend)*rebuilt(stage)%low_faces(p)%x+blend*rebuilt(stage)%native_faces(p)%x
            rebuilt(stage)%faces(p)%y=(1-blend)*rebuilt(stage)%low_faces(p)%y+blend*rebuilt(stage)%native_faces(p)%y
            rebuilt(stage)%faces(p)%z=(1-blend)*rebuilt(stage)%low_faces(p)%z+blend*rebuilt(stage)%native_faces(p)%z
          endif
        enddo
        if(stage==0) then
          stage_registers0=rebuilt(stage)%registers
        else
          registers=rebuilt(stage)%registers
          if(stage>1.and..not.(chemistry.and.present(internal_energy))) then
            valid=.true.
            do p=1,size(owners)
              if(owners(p)/=rank.or..not.active(p)) cycle
              d=tree%nodes(p)%shape
              call accumulate_state_redist_coarse_register_3d(registers(p),accumulated(p),nv,d+12,safe)
              if(.not.safe) then
                valid=.false.;exit
              endif
            enddo
            call agree(valid,comm,safe)
            if(.not.safe) return
          endif
          call move_alloc(registers,accumulated)
        endif
      enddo
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        d0%fields(p)%values=rebuilt(0)%spatial_rate%fields(p)%values
        d1%fields(p)%values=rebuilt(last_stage)%spatial_rate%fields(p)%values
        f0(p)=rebuilt(0)%faces(p);f1(p)=rebuilt(last_stage)%faces(p)
        if(present(internal_energy)) then
          held_rate0%fields(p)%values=rebuilt(0)%held_spatial_rate%fields(p)%values
          held_rate1%fields(p)%values=rebuilt(last_stage)%held_spatial_rate%fields(p)%values
          held_flux0(p)=rebuilt(0)%held_faces(p);held_flux1(p)=rebuilt(last_stage)%held_faces(p)
        endif
      enddo
      if(held_auxiliary) then
        call add_held_auxiliary_rate(held_rate0,rebuilt(0)%auxiliary_rate,safe)
        if(.not.safe) return
        call add_held_auxiliary_rate(held_rate1,rebuilt(last_stage)%auxiliary_rate,safe)
        if(.not.safe) return
      endif
      if(present(internal_energy)) call move_alloc(held_accumulated,held_registers1)
      call move_alloc(accumulated,stage_registers1)
      pending_history%prepared_stages(0:last_stage)=rebuilt
      spatial_weight=blend;good=.true.
    end subroutine apply_prepared_mol_family

    subroutine limit_mol_flux_redist_pass(good)
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: fixed_rate,low_rate,low_transfer,raw_rate
      type(eb_integrated_flux_3d), allocatable :: low_flux(:)
      real(dp) :: low(nv),high(nv),blend,weight,divergence(nv),volume
      integer :: p,d(3),x,y,z,c,retreat
      logical :: valid,safe
      good=.false.;fixed_rate=e0
      call compute_mol_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,mol_geometry,old,ot, &
        dt,1,raw_rate,low_flux,comm,safe,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start,.true.,frame)
      if(.not.safe) return
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        associate(g=>geometry(p),f=>low_flux(p))
        do z=1,g%nz
        do y=1,g%ny
        do x=1,g%nx
          if(g%volume_fraction(x,y,z)<=0) cycle
          volume=g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz
          divergence=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(:,x-1,y,z)- &
            g%x_face_fraction(x,y,z)*f%x(:,x,y,z))+ &
            g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(:,x,y-1,z)- &
            g%y_face_fraction(x,y,z)*f%y(:,x,y,z))+ &
            g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(:,x,y,z-1)- &
            g%z_face_fraction(x,y,z)*f%z(:,x,y,z))
          raw_rate%fields(p)%values(:,x,y,z)=divergence/volume/dt+ &
            .5_dp*(nonface0%fields(p)%values(:,x,y,z)+nonface1%fields(p)%values(:,x,y,z))
        enddo
        enddo
        enddo
        end associate
      enddo
      low_rate=raw_rate
      if(present(redistribution_transfer)) then
        low_transfer=raw_rate
        do p=1,size(owners)
          if(owners(p)==rank) low_transfer%fields(p)%values=0
        enddo
      endif
      valid=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        call post_eb_spatial_density_rate_3d(geometry(p),low_rate%fields(p)%values,valid)
        if(.not.valid) exit
      enddo
      call agree(valid,comm,safe)
      if(.not.safe) return
      ! Only accepted external mass can fund transport. I_R predicts a
      ! future chemical increment and is not an available species inventory.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        fixed_rate%fields(p)%values=.5_dp*(e0%fields(p)%values+e1%fields(p)%values)
        if(present(auxiliary_source)) fixed_rate%fields(p)%values=fixed_rate%fields(p)%values+ &
          .5_dp*(a0%fields(p)%values+a1%fields(p)%values)
      enddo
      valid=.true.;blend=1
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          low=old%fields(p)%values(:,x,y,z)+ &
            dt*(fixed_rate%fields(p)%values(:,x,y,z)+low_rate%fields(p)%values(:,x,y,z))
          high=old%fields(p)%values(:,x,y,z)+dt*(fixed_rate%fields(p)%values(:,x,y,z)+ &
            .5_dp*(d0%fields(p)%values(:,x,y,z)+d1%fields(p)%values(:,x,y,z)))
          valid=valid.and.all(ieee_is_finite(low)).and.all(ieee_is_finite(high))
          valid=valid.and.low(1)>0.and.all(low(6:)>=0)
          if(.not.valid) cycle
          do c=1,nv
            if(c>=2.and.c<=5) cycle
            if(high(c)<0.or.(c==1.and.high(c)==0)) then
              weight=low(c)/(low(c)-high(c))
              if(weight>0) weight=weight*(1-64*epsilon(1._dp))
              blend=min(blend,weight)
            endif
          enddo
        enddo
        enddo
        enddo
      enddo
      if(present(redistribution_transfer)) then
        do p=1,size(owners)
          if(owners(p)/=rank.or.active(p)) cycle
          d=tree%nodes(p)%shape
          do z=1,d(3)
          do y=1,d(2)
          do x=1,d(1)
            if(geometry(p)%volume_fraction(x,y,z)<=0.or.covered(tree,p,[x,y,z])) cycle
            low=old%fields(p)%values(:,x,y,z)
            high=low+.5_dp*(transfer0%fields(p)%values(:,x,y,z)+transfer1%fields(p)%values(:,x,y,z))
            valid=valid.and.all(ieee_is_finite(high)).and.low(1)>0.and.all(low(6:)>=0)
            if(.not.valid) cycle
            do c=1,nv
              if(c>=2.and.c<=5) cycle
              if(high(c)<0.or.(c==1.and.high(c)==0)) then
                weight=low(c)/(low(c)-high(c))
                if(weight>0) weight=weight*(1-64*epsilon(1._dp))
                blend=min(blend,weight)
              endif
            enddo
          enddo
          enddo
          enddo
        enddo
      endif
      call agree(valid,comm,safe)
      if(.not.safe) return
      call MPI_Allreduce(blend,weight,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      blend=weight
      call copy_field(low_rate,admitted_mean_rate,safe)
      if(.not.safe) return
      call copy_mol_level_flux(tree,owners,rank,selected_depth,low_flux,admitted_mean_faces,comm,safe)
      if(.not.safe) return
      if(present(redistribution_transfer)) then
        call copy_field(low_transfer,admitted_mean_transfer,safe,all_levels=.true.)
        if(.not.safe) return
      endif
      ! Validate the same arithmetic that the reactor will consume. Encoding
      ! the selected mean as 2*mean-old_stage and averaging it back can erase
      ! the last trace-species inventory. Keep the mean authoritative, and
      ! retreat the entire conservative family together if rounding exhausts
      ! a budget. The already admitted low family remains the final anchor.
      do retreat=0,65
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          if(present(redistribution_transfer)) then
            if(blend==0) then
              admitted_mean_transfer%fields(p)%values=low_transfer%fields(p)%values
            else
              admitted_mean_transfer%fields(p)%values=(1-blend)*low_transfer%fields(p)%values+ &
                .5_dp*blend*(transfer0%fields(p)%values+transfer1%fields(p)%values)
            endif
          endif
          if(.not.active(p)) cycle
          if(blend==0) then
            admitted_mean_rate%fields(p)%values=low_rate%fields(p)%values
            admitted_mean_faces(p)%x=low_flux(p)%x
            admitted_mean_faces(p)%y=low_flux(p)%y
            admitted_mean_faces(p)%z=low_flux(p)%z
          else
            admitted_mean_rate%fields(p)%values=(1-blend)*low_rate%fields(p)%values+ &
              .5_dp*blend*(d0%fields(p)%values+d1%fields(p)%values)
            admitted_mean_faces(p)%x=(1-blend)*low_flux(p)%x+.5_dp*blend*(f0(p)%x+f1(p)%x)
            admitted_mean_faces(p)%y=(1-blend)*low_flux(p)%y+.5_dp*blend*(f0(p)%y+f1(p)%y)
            admitted_mean_faces(p)%z=(1-blend)*low_flux(p)%z+.5_dp*blend*(f0(p)%z+f1(p)%z)
          endif
        enddo
        valid=.true.
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          if(.not.active(p).and..not.present(redistribution_transfer)) cycle
          d=tree%nodes(p)%shape
          do z=1,d(3)
          do y=1,d(2)
          do x=1,d(1)
            if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
            if(active(p)) then
              high=old%fields(p)%values(:,x,y,z)+dt*(fixed_rate%fields(p)%values(:,x,y,z)+ &
                admitted_mean_rate%fields(p)%values(:,x,y,z))
            else
              if(covered(tree,p,[x,y,z])) cycle
              high=old%fields(p)%values(:,x,y,z)+admitted_mean_transfer%fields(p)%values(:,x,y,z)
            endif
            valid=valid.and.all(ieee_is_finite(high)).and.high(1)>0.and.all(high(6:)>=0)
          enddo
          enddo
          enddo
        enddo
        call agree(valid,comm,safe)
        if(safe) exit
        if(blend==0) return
        if(retreat<64) then
          blend=ieee_next_after(blend,0._dp)
        else
          blend=0
        endif
      enddo
      if(.not.safe) return
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        if(present(redistribution_transfer)) then
          transfer1%fields(p)%values=2*admitted_mean_transfer%fields(p)%values-transfer0%fields(p)%values
        endif
        if(.not.active(p)) cycle
        d1%fields(p)%values=2*admitted_mean_rate%fields(p)%values-d0%fields(p)%values
        f1(p)%x=2*admitted_mean_faces(p)%x-f0(p)%x
        f1(p)%y=2*admitted_mean_faces(p)%y-f0(p)%y
        f1(p)%z=2*admitted_mean_faces(p)%z-f0(p)%z
      enddo
      theta=min(theta,blend);good=.true.
    end subroutine limit_mol_flux_redist_pass

    subroutine combined_rates(with_reaction)
      logical, intent(in) :: with_reaction
      integer :: p
      rates=e0
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        if(allocated(admitted_mean_rate%fields)) then
          rates%fields(p)%values=.5_dp*(e0%fields(p)%values+e1%fields(p)%values)
          if(present(auxiliary_source)) rates%fields(p)%values=rates%fields(p)%values+ &
            .5_dp*(a0%fields(p)%values+a1%fields(p)%values)
          rates%fields(p)%values=rates%fields(p)%values+admitted_mean_rate%fields(p)%values
        else
          rates%fields(p)%values=.5_dp*(e0%fields(p)%values+e1%fields(p)%values+d0%fields(p)%values+d1%fields(p)%values)
          if(present(auxiliary_source)) rates%fields(p)%values=rates%fields(p)%values+ &
            .5_dp*(a0%fields(p)%values+a1%fields(p)%values)
        endif
        if(.not.with_reaction.and.present(energy_correction).and.prepared_redist.and.it==1) &
          rates%fields(p)%values(5,:,:,:)=rates%fields(p)%values(5,:,:,:)+.5_dp*predictor_reset%fields(p)%values(5,:,:,:)/dt
        if(.not.with_reaction) cycle
        rates%fields(p)%values(6:,:,:,:)=rates%fields(p)%values(6:,:,:,:)+ir%fields(p)%values(:nv-5,:,:,:)
        rates%fields(p)%values(5,:,:,:)=rates%fields(p)%values(5,:,:,:)+ir%fields(p)%values(nv-4,:,:,:)
      end do
      if(.not.with_reaction.and.present(energy_correction).and.prepared_redist.and.it==1) &
        accepted_predictor_reset=.5_dp*predictor_reset_change
    end subroutine
    subroutine add_held_auxiliary_stage(stage,good)
      integer, intent(in) :: stage
      logical, intent(out) :: good
      if(stage==0) then
        call add_held_auxiliary_rate(held_rate0,auxiliary_packed0,good)
      else
        call add_held_auxiliary_rate(held_rate1,auxiliary_packed1,good)
      endif
    end subroutine add_held_auxiliary_stage

    subroutine add_held_auxiliary_rate(rate,packet,good)
      type(amr_tree_state_3d), intent(inout) :: rate
      type(amr_tree_state_3d), intent(in) :: packet
      logical, intent(out) :: good
      integer :: p
      logical :: local_ok
      good=.true.;if(.not.held_auxiliary) return
      local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        call guarded_spatial_add(rate%fields(p)%values,packet%fields(p)%values(naux:naux,:,:,:),local_ok)
        if(.not.local_ok) exit
      enddo
      call agree(local_ok,comm,good)
    end subroutine add_held_auxiliary_rate

    subroutine sources(gas,temp,external_rate,transport_rate,auxiliary_rate,stage,flux,exchange,good,fractions)
      type(amr_tree_state_3d) :: mol_transport_source,mol_transport_old
      real(dp), allocatable :: mol_body_exchange(:),mol_body_frequency
      type(amr_tree_state_3d), allocatable :: transport_held
      type(amr_tree_state_3d), intent(in) :: gas,temp
      type(amr_tree_state_3d), intent(out) :: external_rate,transport_rate,auxiliary_rate,exchange
      integer, intent(in) :: stage
      type(eb_integrated_flux_3d), allocatable, intent(out) :: flux(:)
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fractions(:)
      integer :: p,d(3),x,y,z,alloc_status
      logical :: local_ok,all_ok,rotation_ok
      type(amr_tree_state_3d) :: packet,packed_state,source_held
      type(amr_tree_state_3d) :: extra
      type(eb_integrated_flux_3d), allocatable :: extra_flux(:)
      real(dp) :: extra_theta,global_theta,position(3),rotation_source(3)
      local_ok=.true.;good=.false.
      godunov_source_stage=stage
      allocate(external_rate%fields(size(owners)),stat=alloc_status)
      call agree(alloc_status==0,comm,all_ok)
      if(.not.all_ok) return
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        d=tree%nodes(p)%shape
        allocate(external_rate%fields(p)%values(nv,d(1),d(2),d(3)),source=0._dp,stat=alloc_status)
        if(alloc_status/=0) then
          local_ok=.false.;exit
        end if
        if(.not.active(p)) cycle
        if(cached_native_source(stage)) then
          external_rate%fields(p)%values=native_source_replay%external_stages(stage)%fields(p)%values
          cycle
        endif
        if(present(replay_interval).and.(prepared_redist.or.native_external_frozen)) then
          if(prepared_redist) then
            external_rate%fields(p)%values=replay_interval%prepared_stages(stage)%external_rate%fields(p)%values
          else
            external_rate%fields(p)%values=replay_interval%external_stages(stage)%fields(p)%values
          endif
          cycle
        endif
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          if(.not.all(ieee_is_finite(gas%fields(p)%values(:,x,y,z)))) then
            local_ok=.false.;cycle
          end if
          if(gas%fields(p)%values(1,x,y,z)<=0) then
            local_ok=.false.;cycle
          end if
          ! force is force per volume; accel is acceleration evaluated at this stage.
          external_rate%fields(p)%values(2:4,x,y,z)=force+gas%fields(p)%values(1,x,y,z)*accel
          external_rate%fields(p)%values(5,x,y,z)=dot_product( &
            external_rate%fields(p)%values(2:4,x,y,z),gas%fields(p)%values(2:4,x,y,z))/ &
            gas%fields(p)%values(1,x,y,z)
          if(present(frame)) then
            associate(g=>geometry(p))
            position=[g%x_lower,g%y_lower,g%z_lower]+(real([x,y,z],dp)-.5_dp)*[g%dx,g%dy,g%dz]
            end associate
            call rotating_momentum_source_3d(frame,position,gas%fields(p)%values(1,x,y,z), &
              gas%fields(p)%values(2:4,x,y,z),rotation_source,rotation_ok)
            if(.not.rotation_ok) then
              local_ok=.false.;cycle
            end if
            external_rate%fields(p)%values(2:4,x,y,z)=external_rate%fields(p)%values(2:4,x,y,z)+rotation_source
          end if
        end do
        end do
        end do
      end do
      call agree(local_ok,comm,all_ok)
      if(.not.all_ok) return
      if(allocated(pending_history%external_stages)) then
        ! Cache the full owner-local source field, including inactive zero packets.
        call copy_field(external_rate,pending_history%external_stages(stage),all_ok,all_levels=.true.)
        if(.not.all_ok) return
      endif
      transport_rate=external_rate
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        transport_rate%fields(p)%values=0
      end do
      if(present(auxiliary_source)) then
        auxiliary_rate=transport_rate;packet=transport_rate
        if(held_auxiliary) then
          local_ok=.true.
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            d=tree%nodes(p)%shape
            deallocate(packet%fields(p)%values)
            allocate(packet%fields(p)%values(naux,d(1),d(2),d(3)),source=0._dp,stat=alloc_status)
            if(alloc_status/=0) then
              local_ok=.false.;exit
            endif
          enddo
          call agree(local_ok,comm,all_ok);if(.not.all_ok) return
        endif
        if(cached_native_source(stage)) then
          do p=1,size(owners)
            if(owners(p)/=rank.or..not.active(p)) cycle
            packet%fields(p)%values=native_source_replay%auxiliary_stages(stage)%fields(p)%values
          enddo
        else if(cached_godunov_auxiliary(stage)) then
          call copy_godunov_field_3d(godunov_replay_interval%sources(stage)%auxiliary,packet,local_ok)
          call agree(local_ok,comm,all_ok);if(.not.all_ok) return
        else if(present(replay_interval)) then
          do p=1,size(owners)
            if(owners(p)/=rank.or..not.active(p)) cycle
            if(prepared_redist) then
              packet%fields(p)%values=replay_interval%prepared_stages(stage)%auxiliary_rate%fields(p)%values
            else
              packet%fields(p)%values=replay_interval%auxiliary_stages(stage)%fields(p)%values
            endif
          enddo
        else if(.not.(mol.and.stage>1)) then
          if(held_auxiliary) then
            if(stage==0) then
              source_held=held_old
            else
              source_held=held_candidate
            endif
            call pack_eb_held_auxiliary_state_3d(gas,source_held,packed_state,local_ok)
            call agree(local_ok,comm,all_ok);if(.not.all_ok) return
            godunov_packet_pending=.true.
            call auxiliary_source(packed_state,temp,dt,stage,selected_depth,packet,comm,local_ok)
          else
            godunov_packet_pending=.true.
            call auxiliary_source(gas,temp,dt,stage,selected_depth,packet,comm,local_ok)
          endif
          call agree(local_ok,comm,all_ok);if(.not.all_ok) return
        endif
        local_ok=allocated(packet%fields)
        if(local_ok) local_ok=size(packet%fields)==size(owners).and.lbound(packet%fields,1)==1
        if(local_ok) then
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            if(.not.allocated(packet%fields(p)%values)) then
              local_ok=.false.;exit
            endif
            if(any(lbound(packet%fields(p)%values)/=1)) then
              local_ok=.false.;exit
            endif
          enddo
        endif
        call agree(local_ok,comm,all_ok);if(.not.all_ok) return
        call validate_mpi_amr_tree_local_3d(tree,owners,packet,naux,comm,all_ok)
        if(.not.all_ok) return
        local_ok=.true.
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          if(.not.active(p)) then
            packet%fields(p)%values=0;cycle
          endif
          d=tree%nodes(p)%shape
          do z=1,d(3)
          do y=1,d(2)
          do x=1,d(1)
            if(geometry(p)%volume_fraction(x,y,z)<=0) packet%fields(p)%values(:,x,y,z)=0
          enddo
          enddo
          enddo
          local_ok=local_ok.and.all(ieee_is_finite(packet%fields(p)%values))
          auxiliary_rate%fields(p)%values=packet%fields(p)%values(:nv,:,:,:)
        enddo
        call agree(local_ok,comm,all_ok);if(.not.all_ok) return
        if(allocated(godunov_history)) then
          call capture_auxiliary_packet(stage,packet,all_ok)
          if(.not.all_ok) return
        endif
        godunov_packet_pending=.false.
        if(allocated(pending_history%auxiliary_stages)) then
          call copy_field(packet,pending_history%auxiliary_stages(stage),all_ok)
          if(.not.all_ok) return
        endif
        if(stage==0) then
          auxiliary_packed0=packet
        else
          auxiliary_packed1=packet
        endif
      endif
      if(native_prefix_agreed) native_prefix_observed=max(native_prefix_observed,stage+1)
      if(prepared_redist) then
        if(prepared_godunov) then
          call prepared_diffusion_sources(gas,temp,stage,transport_rate,flux,stage_theta,all_ok,fractions)
        else
          call prepared_mol_sources(gas,temp,stage,transport_rate,flux,stage_theta,all_ok,fractions)
        endif
        if(.not.all_ok) return
        if(capture_prepared) then
          call capture_prepared_volume_stage(stage,external_rate,packet,all_ok)
          if(.not.all_ok) return
        endif
        if(held_auxiliary) then
          call add_held_auxiliary_stage(stage,all_ok)
          if(.not.all_ok) return
        endif
        theta=min(theta,stage_theta)
        good=.true.
        return
      endif
      if(selected_depth>=0.and.redistribution_method/='none') exchange=transport_rate
      if(diffusive) then
        if(present(internal_energy)) then
          allocate(transport_held)
          if(stage==0) then
            transport_held=held_old
          else
            transport_held=held_candidate
          endif
        endif
        ! CTU also needs both tangential stencils across patch faces.
        ! Retain its original coarse interpolation and energy convention.
        if(.not.mol.or.present(internal_energy)) then
          if(wall_active) then
            allocate(mol_body_exchange(nv),mol_body_frequency,stat=alloc_status)
            call agree(alloc_status==0,comm,all_ok)
            if(.not.all_ok) return
          endif
          call compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,tree,owners,geometry,extended_geometry, &
            gas,temp,dt,transport_flags,transport_rate,flux,mol_transport_source,mol_transport_old,stage_theta,comm,all_ok, &
            boundary_faces,inflow_state,inflow_temperature,wall_temperature,target_depth,ancestor_start,ancestor_end, &
            temperature_start,temperature_end,fractions,frame,embedded_wall_temperature,embedded_noslip, &
            mol_body_exchange,mol_body_frequency,inflow_mach,inflow_pressure, &
            constant_coefficients=constant_coefficients,internal_energy=transport_held, &
            ancestor_internal_energy_start=ancestor_internal_energy_start, &
            ancestor_internal_energy_end=ancestor_internal_energy_end,coarse_linear=mol.and.present(internal_energy), &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
          primitive_from_total_energy=mol.and..not.present(frame).and..not.selected_momentum)
          if(.not.all_ok) return
          if(wall_active) then
            heat_stage=mol_body_exchange(5);heat_frequency=mol_body_frequency
            if(noslip) wall_impulse=mol_body_exchange(2:4)
          endif
        else
        call compute_sdc_mpi_eb_tree_transport_3d(species,transport,tree,owners,geometry,dt, &
          transport_flags(1),transport_flags(2),transport_flags(3),transport_flags(4),gas,temp, &
          transport_rate,flux,stage_theta,comm,all_ok,inflow_mach,inflow_temperature,inflow_pressure, &
          boundary_faces,inflow_state,wall_temperature, &
          target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions, &
          embedded_wall_temperature=embedded_wall_temperature,extended_geometry=heat_geometry, &
          embedded_wall_energy=heat_stage,embedded_frequency=heat_frequency, &
          embedded_noslip=embedded_noslip,embedded_wall_impulse=wall_impulse,frame=frame, &
          constant_coefficients=constant_coefficients,internal_energy=transport_held, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end)
        endif
        if(wall_active.and.all_ok) then
          if(stage==0) then
            heat_old=heat_stage
            if(noslip) impulse_old=wall_impulse
          else
            heat_new=heat_stage
            if(noslip) impulse_new=wall_impulse
          end if
        end if
        if(.not.all_ok) return
      else
        stage_theta=1
      end if
      if(mol) then
        if(present(internal_energy)) then
          if(stage==0) then
            call held_mol_source(gas,temp,held_old,held_rate0,held_flux0,held_registers0,all_ok,fractions, &
              extra,extra_flux)
          else
            call held_mol_source(gas,temp,held_candidate,held_rate1,held_flux1,held_registers1,all_ok,fractions, &
              extra,extra_flux)
          endif
        else
        call compute_mol_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,mol_geometry,gas,temp,dt,2, &
          extra,extra_flux,comm,all_ok,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions,mol_lax_friedrichs,frame)
        endif
        if(.not.all_ok) return
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          transport_rate%fields(p)%values=transport_rate%fields(p)%values+extra%fields(p)%values
          if(diffusive) then
            flux(p)%x=flux(p)%x+extra_flux(p)%x
            flux(p)%y=flux(p)%y+extra_flux(p)%y
            flux(p)%z=flux(p)%z+extra_flux(p)%z
          end if
        end do
        if(.not.diffusive) call move_alloc(extra_flux,flux)
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          extra%fields(p)%values=0
          if(.not.active(p)) cycle
          call wall_source(gas,temp,p,extra%fields(p)%values,local_ok,stage)
          if(.not.local_ok) exit
        end do
        call agree(local_ok,comm,all_ok)
        if(.not.all_ok) return
        if(stage==0) then
          call integrate_mpi_eb_tree_3d(tree,owners,geometry,extra,mol_wall_old,comm,all_ok)
          mol_wall_old=dt*mol_wall_old
        else
          call integrate_mpi_eb_tree_3d(tree,owners,geometry,extra,mol_wall_new,comm,all_ok)
          mol_wall_new=dt*mol_wall_new
        end if
        if(.not.all_ok) return
      end if
      if(has_diffusion.and..not.(external_diffusion_source.and.stage>1)) then
        extra=transport_rate;extra_theta=1
        do p=1,size(owners)
          if(owners(p)==rank) extra%fields(p)%values=0
        enddo
        ! The additional operator owns a separate face family; MOL scratch
        ! faces must not escape as incoming callback storage.
        if(allocated(extra_flux)) deallocate(extra_flux)
        if(cached_godunov_provider(godunov_source_stage).and.external_diffusion_source) then
          call copy_godunov_field_3d(godunov_replay_interval%sources(stage)%provider_rate,extra,local_ok)
          call agree(local_ok,comm,all_ok);if(.not.all_ok) return
          call copy_godunov_fluxes_3d(godunov_replay_interval%sources(stage)%provider_faces,extra_flux,local_ok)
          call agree(local_ok,comm,all_ok);if(.not.all_ok) return
          extra_theta=godunov_replay_interval%sources(stage)%provider_limiter
        else if(present(replay_interval).and.external_diffusion_source.and..not.prepared_redist) then
          associate(saved=>replay_interval%provider_stages(2,stage))
          call copy_godunov_field_3d(saved%provider_rate,extra,local_ok)
          call agree(local_ok,comm,all_ok);if(.not.all_ok) return
          call copy_godunov_fluxes_3d(saved%provider_faces,extra_flux,local_ok)
          call agree(local_ok,comm,all_ok);if(.not.all_ok) return
          extra_theta=saved%provider_limiter
          end associate
        else
        godunov_packet_pending=.true.
        call evaluate_diffusion(gas,temp,dt,selected_depth,extra,extra_flux,extra_theta,comm,local_ok, &
          ancestor_start,ancestor_end,temperature_start,temperature_end,fractions)
        call agree(local_ok,comm,all_ok)
        if(.not.all_ok) then
          godunov_recovery_forbidden=.true.
        if(present(provider_refused)) provider_refused=.true.
          return
        endif
        endif
        local_ok=allocated(extra%fields)
        if(local_ok) local_ok=size(extra%fields)==size(owners).and.lbound(extra%fields,1)==1
        if(local_ok) then
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            if(.not.allocated(extra%fields(p)%values)) then
              local_ok=.false.;exit
            endif
            if(any(lbound(extra%fields(p)%values)/=1)) then
              local_ok=.false.;exit
            endif
          enddo
        endif
        call agree(local_ok,comm,all_ok);if(.not.all_ok) return
        call validate_mpi_amr_tree_local_3d(tree,owners,extra,nv,comm,all_ok)
        if(.not.all_ok) return
        local_ok=allocated(extra_flux).and.ieee_is_finite(extra_theta)
        if(local_ok) local_ok=size(extra_flux)==size(owners).and.lbound(extra_flux,1)==1.and. &
          extra_theta>=0.and.extra_theta<=1
        call agree(local_ok,comm,all_ok)
        if(.not.all_ok) return
        local_ok=.true.
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          local_ok=local_ok.and.valid_eb_integrated_flux_3d(extra_flux(p),geometry(p),nv)
        enddo
        call agree(local_ok,comm,all_ok);if(.not.all_ok) return
        if(allocated(godunov_history)) then
          call capture_provider_packet(stage,2,extra,extra,extra_flux,extra_theta,all_ok)
          if(.not.all_ok) return
        endif
        if(allocated(pending_history%provider_stages)) then
          call capture_native_mol_provider(stage,2,extra,extra,extra_flux,extra_theta,all_ok)
          if(.not.all_ok) return
        endif
        godunov_packet_pending=.false.
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          local_ok=valid_eb_integrated_flux_3d(extra_flux(p),geometry(p),nv)
          local_ok=local_ok.and.all(ieee_is_finite(extra%fields(p)%values))
          if(.not.local_ok) exit
          call guarded_spatial_add(transport_rate%fields(p)%values,extra%fields(p)%values,local_ok)
          if(.not.local_ok) exit
          if(allocated(flux)) then
            call guarded_spatial_add(flux(p)%x,extra_flux(p)%x,local_ok)
            if(.not.local_ok) exit
            call guarded_spatial_add(flux(p)%y,extra_flux(p)%y,local_ok)
            if(.not.local_ok) exit
            call guarded_spatial_add(flux(p)%z,extra_flux(p)%z,local_ok)
            if(.not.local_ok) exit
          end if
        end do
        call agree(local_ok,comm,all_ok)
        if(.not.all_ok) return
        call MPI_Allreduce(extra_theta,global_theta,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
        if(status/=MPI_SUCCESS) return
        stage_theta=min(stage_theta,global_theta)
        if(.not.allocated(flux)) call move_alloc(extra_flux,flux)
      end if
      if(present(evaluate_spatial_diffusion).and..not.(external_diffusion_source.and.stage>1)) then
        if(present(internal_energy).and..not.allocated(transport_held)) then
          allocate(transport_held,stat=alloc_status)
          call agree(alloc_status==0,comm,all_ok);if(.not.all_ok) return
          if(stage==0) then
            transport_held=held_old
          else
            transport_held=held_candidate
          endif
        endif
        call ordinary_spatial_diffusion(gas,temp,transport_rate,flux,stage_theta,all_ok,fractions,transport_held)
        if(.not.all_ok) return
      endif
      if(present(internal_energy).and..not.mol) then
        call initialize_held_godunov_stage(gas,temp,stage,good,fractions)
        if(.not.good) return
      endif
      if(held_auxiliary) then
        call add_held_auxiliary_stage(stage,all_ok)
        if(.not.all_ok) return
      endif
      if(mol.and.chemistry.and.redistribution_method=='flux_redist') then
        call capture_nonface_rate(transport_rate,flux,stage,all_ok)
        if(.not.all_ok) return
      endif
      if(diffusive.or.has_diffusion.or.present(evaluate_spatial_diffusion).or.mol) then
        if(redistribution_method=='state_redist') then
          if(redist_order>0) then
            call redist_boundary_samples(gas,temp,all_ok)
            if(.not.all_ok) return
          endif
          if(selected_depth>=0) then
            call redistribute_sdc_mpi_eb_tree_state_level_rate_3d(tree,owners,geometry,gas,transport_rate,nv,selected_depth,dt, &
              redistributed,exchange,comm,all_ok,redist_target,boundary_faces, &
      max_order=redist_order,boundary_halos=redist_halos,boundary_geometry=redist_geometry)
          else
            call redistribute_sdc_mpi_eb_tree_state_rate_3d(tree,owners,geometry,gas,transport_rate,nv,dt,redistributed, &
              comm,all_ok,redist_target,boundary_faces, &
      max_order=redist_order,boundary_halos=redist_halos,boundary_geometry=redist_geometry)
          end if
          if(.not.all_ok) return
          transport_rate=redistributed
        else if(redistribution_method=='flux_redist') then
          if(selected_depth>=0) then
            call redistribute_sdc_mpi_eb_tree_level_rate_3d(tree,owners,geometry,transport_rate,nv,selected_depth,dt, &
              redistributed,exchange,comm,all_ok,boundary_faces)
          else
            call redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree,owners,geometry,transport_rate,nv,redistributed, &
              comm,all_ok,boundary_faces)
          end if
          if(.not.all_ok) return
          transport_rate=redistributed
        end if
        local_ok=.true.
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          call post_eb_spatial_density_rate_3d(geometry(p),transport_rate%fields(p)%values,local_ok)
          if(.not.local_ok) exit
        enddo
        call agree(local_ok,comm,all_ok);if(.not.all_ok) return
        theta=min(theta,stage_theta)
      end if
      if(native_mol_weight/=1) then
        if(stage==0) then
          call blend_native_mol_stage(transport_rate,flux,held_rate0,held_flux0,good)
        else
          call blend_native_mol_stage(transport_rate,flux,held_rate1,held_flux1,good)
        endif
        if(.not.good) return
      endif
      good=.true.
      if(allocated(flux)) call initialize_inactive_faces(flux,nv,good)
    end subroutine
    subroutine ordinary_spatial_diffusion(gas,temp,output,flux,limiter,good,fractions,stored)
      type(amr_tree_state_3d), intent(in) :: gas,temp
      type(amr_tree_state_3d), optional, intent(in) :: stored
      type(amr_tree_state_3d), intent(inout) :: output
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: limiter
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fractions(:)
      type(eb_geometry_3d), allocatable :: padded(:)
      type(amr_tree_state_3d) :: source,old_halo,old_temp,volume,metrics
      type(amr_tree_state_3d), allocatable :: scalar_old
      integer :: p,d(3),padding,allocation_status
      logical :: local_ok
      good=.false.;local_ok=.true.
      allocate(padded(size(owners)),metrics%fields(size(owners)),source%fields(size(owners)),stat=allocation_status)
      call agree(allocation_status==0,comm,good);if(.not.good) return
      do p=1,size(owners)
        d=tree%nodes(p)%shape; padding=(extended_geometry(p)%nx-d(1))/2
        call extract_eb_geometry_3d(extended_geometry(p),[padding-5,padding-5,padding-5], &
          d+padding+6,padded(p),local_ok)
        if(.not.local_ok) exit
        if(owners(p)/=rank) cycle
        allocate(metrics%fields(p)%values(1,-5:d(1)+6,-5:d(2)+6,-5:d(3)+6),stat=allocation_status)
        if(allocation_status/=0) then
          local_ok=.false.;exit
        endif
        metrics%fields(p)%values(1,:,:,:)=padded(p)%volume_fraction
        if(.not.active(p)) cycle
        allocate(source%fields(p)%values(nv,d(1)+12,d(2)+12,d(3)+12),source=0._dp,stat=allocation_status)
        if(allocation_status/=0) then
          local_ok=.false.;exit
        endif
      enddo
      call agree(local_ok,comm,good);if(.not.good) return
      if(present(stored)) then
        allocate(scalar_old,stat=allocation_status)
        call agree(allocation_status==0,comm,good);if(.not.good) return
      endif
      ! Gas, temperature and UEINT donors share the actual stage/ancestor clock.
      call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,gas,temp,old_halo,old_temp,comm,good, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=volume, &
        halo_width=6,target_volume_fraction=metrics,restrict_snapshot=.false.,ancestor_start=ancestor_start, &
        ancestor_end=ancestor_end,temperature_start=temperature_start,temperature_end=temperature_end, &
        ancestor_fraction=fractions,frame=frame,internal_energy=stored,ghost_internal_energy=scalar_old, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end)
      if(.not.good) return
      call add_spatial_diffusion_callback(gas,temp,padded,old_halo,output,source,flux,limiter,good,fractions,stored)
    end subroutine ordinary_spatial_diffusion

    subroutine guarded_spatial_add(destination,increment,good)
      real(dp), intent(inout) :: destination(:,:,:,:)
      real(dp), intent(in) :: increment(:,:,:,:)
      logical, intent(out) :: good
      real(dp) :: old,delta
      integer :: c,x,y,z
      good=.false.
      if(any(shape(destination)/=shape(increment))) return
      if(.not.all(ieee_is_finite(destination)).or..not.all(ieee_is_finite(increment))) return
      do z=1,size(destination,4)
      do y=1,size(destination,3)
      do x=1,size(destination,2)
      do c=1,size(destination,1)
        old=destination(c,x,y,z);delta=increment(c,x,y,z)
        if(delta>0) then
          if(old>huge(1._dp)-delta) return
        else if(delta<0) then
          if(old< -huge(1._dp)-delta) return
        endif
      enddo
      enddo
      enddo
      enddo
      destination=destination+increment;good=.true.
    end subroutine guarded_spatial_add
    subroutine refresh_callback_halos(callback_rate,callback_halos,success)
      type(amr_tree_state_3d), intent(in) :: callback_rate
      type(amr_tree_state_3d), intent(inout) :: callback_halos
      logical, intent(out) :: success
      integer(int64), allocatable :: origin(:,:),domain(:,:)
      integer, allocatable :: ratios(:),requests(:,:),destinations(:,:)
      real(dp), allocatable :: values(:,:)
      integer(int64) :: cell(3),total
      integer :: p,q,level,i,j,k,index(3),a,m(3),stage,r,allocation_status
      character(len=32) :: faces(6)
      logical :: inside,good
      faces=boundary_faces
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
      call agree(good,comm,success)
      if(.not.success) return
      total=0
      do stage=1,2
        r=0;good=.true.
        do p=1,size(owners)
          if(.not.allocated(callback_halos%fields(p)%values)) cycle
          level=tree%nodes(p)%level;m=tree%nodes(p)%shape
          do k=4,m(3)+9
          do j=4,m(2)+9
          do i=4,m(1)+9
            cell=origin(:,p)+int([i,j,k]-7,int64);inside=.true.
            do a=1,3
              if(faces(2*a-1)=='periodic') then
                cell(a)=modulo(cell(a),domain(a,level+1))
              else if(cell(a)<0.or.cell(a)>=domain(a,level+1)) then
                inside=.false.
              endif
            enddo
            if(.not.inside) cycle
            do q=1,size(owners)
              if(tree%nodes(q)%level/=level) cycle
              if(any(cell<origin(:,q)).or.any(cell>=origin(:,q)+int(tree%nodes(q)%shape,int64))) cycle
              index=int(cell-origin(:,q))+1
              if(owners(q)==rank) then
                if(stage==2) callback_halos%fields(p)%values(:,i,j,k)= &
                  callback_rate%fields(q)%values(:,index(1),index(2),index(3))
              else if(stage==1) then
                total=total+1_int64
              else
                r=r+1;requests(:,r)=[q,index];destinations(:,r)=[p,i,j,k]
              endif
              exit
            enddo
          enddo
          enddo
          enddo
        enddo
        if(stage==1) then
          call agree(total<=int(huge(1)/max(4,nv),int64),comm,success)
          if(.not.success) return
          allocate(requests(4,int(total)),destinations(4,int(total)),stat=allocation_status)
          call agree(allocation_status==0,comm,success)
          if(.not.success) return
        endif
      enddo
      call fetch_mpi_amr_tree_cells_3d(tree,owners,callback_rate,nv,requests,values,comm,success)
      if(.not.success) return
      do r=1,size(requests,2)
        p=destinations(1,r);index=destinations(2:4,r)
        callback_halos%fields(p)%values(:,index(1),index(2),index(3))=values(:,r)
      enddo
    end subroutine refresh_callback_halos
    subroutine prepared_diffusion_sources(gas,temp,stage,output,flux,limiter,good,fractions)
      type(amr_tree_state_3d), intent(in) :: gas,temp
      integer, intent(in) :: stage
      type(amr_tree_state_3d), intent(inout) :: output
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: limiter
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fractions(:)
      type(amr_tree_state_3d), allocatable :: transport_held
      type(amr_tree_state_3d) :: source,old_halo
      type(state_redist_register_3d), allocatable :: registers(:)
      real(dp), allocatable :: body_exchange(:),frequency
      integer :: p,m(3),allocation_status,depth,first_depth,last_depth
      logical :: local_ok,all_ok
      good=.false.
      if(present(internal_energy)) then
        allocate(transport_held)
        if(stage==0) then
          transport_held=held_old
        else
          transport_held=held_candidate
        endif
      endif
      if(wall_active) then
        allocate(body_exchange(nv),frequency,stat=allocation_status)
        call agree(allocation_status==0,comm,all_ok);if(.not.all_ok) return
      endif
      call compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,tree,owners,geometry,register_geometry, &
        gas,temp,dt,transport_flags,output,flux,source,old_halo,limiter,comm,all_ok, &
        boundary_faces,inflow_state,inflow_temperature,wall_temperature,target_depth,ancestor_start,ancestor_end, &
        temperature_start,temperature_end,fractions,frame,embedded_wall_temperature,embedded_noslip, &
        body_exchange,frequency,inflow_mach,inflow_pressure, &
        constant_coefficients=constant_coefficients,internal_energy=transport_held, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end)
      if(.not.all_ok) return
      call add_spatial_diffusion_callback(gas,temp,register_geometry,old_halo,output,source,flux,limiter, &
        all_ok,fractions,transport_held)
      if(.not.all_ok) return
      if(diffusive.or.present(evaluate_spatial_diffusion)) then
        first_depth=max(0,selected_depth);last_depth=first_depth
        if(selected_depth<0) last_depth=tree%level_count()-1
        do depth=first_depth,last_depth
          call compute_state_redist_mpi_level_stage_3d(tree,owners,register_geometry,old_halo,source,nv, &
            depth,dt,output,registers,comm,all_ok,boundary_faces,redist_target,redist_order,register_scale=1._dp)
          if(.not.all_ok) return
        enddo
      else
        allocate(registers(size(owners)),stat=allocation_status)
        call agree(allocation_status==0,comm,all_ok);if(.not.all_ok) return
        local_ok=.true.
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          call initialize_state_redist_tree_register_3d(tree,p,nv,6,registers(p),local_ok,boundary_faces,1._dp)
          if(.not.local_ok) exit
        enddo
        call agree(local_ok,comm,all_ok)
        if(.not.all_ok) return
      endif
      if(wall_active) then
        if(stage==0) then
          heat_old=body_exchange(5);impulse_old=body_exchange(2:4)
        else
          heat_new=body_exchange(5);impulse_new=body_exchange(2:4)
        endif
      endif
      ! Post-redistribution density consistency precedes gas/EOS stepping.
      local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        call post_eb_spatial_density_rate_3d(geometry(p),output%fields(p)%values,local_ok)
        if(.not.local_ok) exit
      enddo
      call agree(local_ok,comm,all_ok);if(.not.all_ok) return
      if(present(internal_energy)) then
        call initialize_held_godunov_stage(gas,temp,stage,good,fractions)
        if(.not.good) return
      endif
      if(stage==0) then
        call move_alloc(registers,stage_registers0)
      else
        local_ok=.true.
        if(stage>1) then
          do p=1,size(owners)
            if(owners(p)/=rank.or..not.active(p)) cycle
            m=tree%nodes(p)%shape
            call accumulate_state_redist_coarse_register_3d(registers(p),stage_registers1(p),nv,m+12,local_ok)
            if(.not.local_ok) exit
          enddo
        endif
        call agree(local_ok,comm,all_ok);if(.not.all_ok) return
        call move_alloc(registers,stage_registers1)
      endif
      good=.true.
    end subroutine prepared_diffusion_sources
    subroutine add_spatial_diffusion_callback(gas,temp,padded,old_halo,output,source,flux,limiter,good,fractions,stored)
      type(amr_tree_state_3d), intent(in) :: gas,temp,old_halo
      type(amr_tree_state_3d), optional, intent(in) :: stored
      type(amr_tree_state_3d), allocatable :: scalar_callback_old
      type(eb_geometry_3d), intent(in) :: padded(:)
      type(amr_tree_state_3d), intent(inout) :: output,source
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: limiter
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fractions(:)
      type(amr_tree_state_3d) :: callback_rate,callback_halos,callback_old,callback_temp,callback_volume,metrics
      type(eb_integrated_flux_3d), allocatable :: callback_flux(:)
      real(dp) :: callback_theta,global_theta
      integer :: p,m(3),allocation_status,width
      logical :: local_ok
      good=.true.
      if(.not.present(evaluate_spatial_diffusion)) return
        callback_old=old_halo
        if(spatial_callback_width>6) then
          ! Extra LES/filter donors use the same stage snapshots and ancestor
          ! interpolation as the built-in six-layer state. No RHS extrapolation.
          width=spatial_callback_width;local_ok=.true.
          do p=1,size(owners)
            if(any(geometry(p)%volume_fraction/=1).or.any(padded(p)%volume_fraction/=1)) local_ok=.false.
            if(any(tree%nodes(p)%shape>huge(1)-2*width)) local_ok=.false.
          enddo
          call agree(local_ok,comm,good);if(.not.good) return
          allocate(metrics%fields(size(owners)),stat=allocation_status)
          call agree(allocation_status==0,comm,good);if(.not.good) return
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            m=tree%nodes(p)%shape
            allocate(metrics%fields(p)%values(1,1-width:m(1)+width,1-width:m(2)+width, &
              1-width:m(3)+width),source=1._dp,stat=allocation_status)
            if(allocation_status/=0) then
              local_ok=.false.;exit
            endif
          enddo
          call agree(local_ok,comm,good);if(.not.good) return
          if(present(stored)) then
            allocate(scalar_callback_old,stat=allocation_status)
            call agree(allocation_status==0,comm,good);if(.not.good) return
          endif
          call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,gas,temp,callback_old,callback_temp,comm,good, &
            boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
            inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=callback_volume, &
            halo_width=width,target_volume_fraction=metrics,restrict_snapshot=.false.,ancestor_start=ancestor_start, &
            ancestor_end=ancestor_end,temperature_start=temperature_start,temperature_end=temperature_end, &
            ancestor_fraction=fractions,frame=frame,internal_energy=stored, &
            ghost_internal_energy=scalar_callback_old, &
            ancestor_internal_energy_start=ancestor_internal_energy_start, &
            ancestor_internal_energy_end=ancestor_internal_energy_end)
          if(.not.good) return
          do p=1,size(owners)
            if(owners(p)/=rank.or..not.active(p)) cycle
            m=tree%nodes(p)%shape
            local_ok=all(callback_old%fields(p)%values(:,-5:m(1)+6,-5:m(2)+6,-5:m(3)+6)== &
              old_halo%fields(p)%values)
            if(.not.local_ok) exit
          enddo
          call agree(local_ok,comm,good);if(.not.good) return
        endif
        callback_rate=output;callback_theta=1
        do p=1,size(owners)
          if(owners(p)==rank) callback_rate%fields(p)%values=0
        enddo
      if(cached_godunov_provider(godunov_source_stage).and.external_diffusion_source) then
        call copy_godunov_field_3d(godunov_replay_interval%sources(godunov_source_stage)%provider_rate,callback_rate,local_ok)
        call agree(local_ok,comm,good);if(.not.good) return
        call copy_godunov_field_3d(godunov_replay_interval%sources(godunov_source_stage)%provider_halos,callback_halos,local_ok)
        call agree(local_ok,comm,good);if(.not.good) return
        call copy_godunov_fluxes_3d(godunov_replay_interval%sources(godunov_source_stage)%provider_faces,callback_flux,local_ok)
        call agree(local_ok,comm,good);if(.not.good) return
        callback_theta=godunov_replay_interval%sources(godunov_source_stage)%provider_limiter
      else if(present(replay_interval).and.external_diffusion_source.and..not.prepared_redist) then
        associate(saved=>replay_interval%provider_stages(1,godunov_source_stage))
        call copy_godunov_field_3d(saved%provider_rate,callback_rate,local_ok)
        call agree(local_ok,comm,good);if(.not.good) return
        call copy_godunov_field_3d(saved%provider_halos,callback_halos,local_ok)
        call agree(local_ok,comm,good);if(.not.good) return
        call copy_godunov_fluxes_3d(saved%provider_faces,callback_flux,local_ok)
        call agree(local_ok,comm,good);if(.not.good) return
        callback_theta=saved%provider_limiter
        end associate
      else
        godunov_packet_pending=.true.
        call evaluate_spatial_diffusion(gas,temp,padded,callback_old,dt,selected_depth, &
          callback_rate,callback_halos,callback_flux,callback_theta,comm,local_ok,frame)
        call agree(local_ok,comm,good)
        if(.not.good) then
        godunov_recovery_forbidden=.true.
        if(present(provider_refused)) provider_refused=.true.
        return
      endif
      endif
        ! Tree indexing and peer fetches use one-based valid source arrays.
        ! Reject malformed callback containers before validators or slices
        ! can address them; collective failure retains the outer transaction.
        local_ok=allocated(callback_rate%fields)
        if(local_ok) local_ok=size(callback_rate%fields)==size(owners).and.lbound(callback_rate%fields,1)==1
        if(local_ok) then
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            if(.not.allocated(callback_rate%fields(p)%values)) then
              local_ok=.false.;exit
            endif
            if(any(lbound(callback_rate%fields(p)%values)/=1)) then
              local_ok=.false.;exit
            endif
          enddo
        endif
        call agree(local_ok,comm,good);if(.not.good) return
        call validate_mpi_amr_tree_local_3d(tree,owners,callback_rate,nv,comm,good)
        if(.not.good) return
        local_ok=allocated(callback_halos%fields).and.allocated(callback_flux).and.ieee_is_finite(callback_theta)
        if(local_ok) local_ok=size(callback_halos%fields)==size(owners).and.size(callback_flux)==size(owners).and. &
          callback_theta>=0.and.callback_theta<=1.and.lbound(callback_halos%fields,1)==1.and.lbound(callback_flux,1)==1
        call agree(local_ok,comm,good)
        if(.not.good) return
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) then
            if(allocated(callback_halos%fields(p)%values)) local_ok=.false.
            cycle
          endif
          m=tree%nodes(p)%shape
          if(.not.allocated(callback_halos%fields(p)%values)) then
            local_ok=.false.;exit
          endif
          if(any(shape(callback_halos%fields(p)%values)/=[nv,m+12]).or. &
            any(lbound(callback_halos%fields(p)%values)/=1)) then
            local_ok=.false.;exit
          endif
          local_ok=all(ieee_is_finite(callback_halos%fields(p)%values)).and. &
            all(ieee_is_finite(callback_rate%fields(p)%values)).and. &
            valid_eb_integrated_flux_3d(callback_flux(p),geometry(p),nv)
          if(.not.local_ok) exit
          local_ok=all(callback_halos%fields(p)%values(:,7:m(1)+6,7:m(2)+6,7:m(3)+6)== &
            callback_rate%fields(p)%values)
          if(.not.local_ok) exit
        enddo
        call agree(local_ok,comm,good)
        if(.not.good) return
        call refresh_callback_halos(callback_rate,callback_halos,good)
        if(.not.good) return
        call MPI_Allreduce(callback_theta,global_theta,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
        if(status/=MPI_SUCCESS) then
          good=.false.;return
        endif
        if(allocated(godunov_history)) then
          call capture_provider_packet(godunov_source_stage,1,callback_rate,callback_halos,callback_flux,global_theta,good)
          if(.not.good) return
        endif
        if(allocated(pending_history%provider_stages)) then
          call capture_native_mol_provider(godunov_source_stage,1,callback_rate,callback_halos, &
            callback_flux,global_theta,good)
          if(.not.good) return
        endif
        godunov_packet_pending=.false.
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          call guarded_spatial_add(output%fields(p)%values,callback_rate%fields(p)%values,local_ok)
          if(.not.local_ok) exit
          call guarded_spatial_add(source%fields(p)%values,callback_halos%fields(p)%values,local_ok)
          if(.not.local_ok) exit
          call guarded_spatial_add(flux(p)%x,callback_flux(p)%x,local_ok)
          if(.not.local_ok) exit
          call guarded_spatial_add(flux(p)%y,callback_flux(p)%y,local_ok)
          if(.not.local_ok) exit
          call guarded_spatial_add(flux(p)%z,callback_flux(p)%z,local_ok)
          if(.not.local_ok) exit
        enddo
        call agree(local_ok,comm,good)
        if(.not.good) return
        limiter=min(limiter,global_theta)
    end subroutine add_spatial_diffusion_callback
    subroutine build_prepared_low_rate(raw_halos,valid_rate,native_faces,low,good)
      type(amr_tree_state_3d), intent(in) :: raw_halos,valid_rate
      type(eb_integrated_flux_3d), intent(in) :: native_faces(:)
      type(amr_tree_state_3d), intent(out) :: low
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: ignored,ignored_halos,ignored_held_rate,ignored_held_halos
      real(dp) :: old_divergence(nv),new_divergence(nv),nonface(nv),volume
      integer :: p,x,y,z
      logical :: valid
      good=.false.
      if(.not.allocated(prepared_low_faces)) then
        if(present(internal_energy)) then
          call compute_mol_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,extended_geometry,old,ot, &
            dt,1,ignored,prepared_low_faces,comm,good,boundary_faces,inflow_state,inflow_mach, &
            inflow_temperature,inflow_pressure,target_depth,ancestor_start,ancestor_end,temperature_start, &
            temperature_end,fractions_start,.true.,frame,old_state_halos=ignored_halos, &
            internal_energy_rate=ignored_held_rate,internal_energy_flux=prepared_low_held_faces, &
            internal_energy_state=held_old,internal_energy_old_state_halos=ignored_held_halos, &
            ancestor_internal_energy_start=ancestor_internal_energy_start, &
            ancestor_internal_energy_end=ancestor_internal_energy_end, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected)
        else
        call compute_mol_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,extended_geometry,old,ot, &
          dt,1,ignored,prepared_low_faces,comm,good,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions_start,.true.,frame, &
          old_state_halos=ignored_halos)
        endif
        if(.not.good) return
      endif
      call copy_field(valid_rate,low,good)
      if(.not.good) return
      valid=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        associate(g=>geometry(p),old_face=>native_faces(p),new_face=>prepared_low_faces(p))
        do z=1,g%nz
        do y=1,g%ny
        do x=1,g%nx
          low%fields(p)%values(:,x,y,z)=0
          if(g%volume_fraction(x,y,z)<=0) cycle
          volume=g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz
          old_divergence=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*old_face%x(:,x-1,y,z)- &
            g%x_face_fraction(x,y,z)*old_face%x(:,x,y,z))+ &
            g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*old_face%y(:,x,y-1,z)- &
            g%y_face_fraction(x,y,z)*old_face%y(:,x,y,z))+ &
            g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*old_face%z(:,x,y,z-1)- &
            g%z_face_fraction(x,y,z)*old_face%z(:,x,y,z))
          new_divergence=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*new_face%x(:,x-1,y,z)- &
            g%x_face_fraction(x,y,z)*new_face%x(:,x,y,z))+ &
            g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*new_face%y(:,x,y-1,z)- &
            g%y_face_fraction(x,y,z)*new_face%y(:,x,y,z))+ &
            g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*new_face%z(:,x,y,z-1)- &
            g%z_face_fraction(x,y,z)*new_face%z(:,x,y,z))
          nonface=raw_halos%fields(p)%values(:,x+6,y+6,z+6)-old_divergence/volume/dt
          ! Built-in closed EB spatial operators have no non-face species or
          ! mass source. Do not invent trace credit by subtracting their fluxes.
          if(.not.present(evaluate_spatial_diffusion)) then
            nonface(1)=0;nonface(6:)=0
          endif
          low%fields(p)%values(:,x,y,z)=new_divergence/volume/dt+nonface
        enddo
        enddo
        enddo
        end associate
        call post_eb_spatial_density_rate_3d(geometry(p),low%fields(p)%values,valid)
        if(.not.valid) exit
      enddo
      call agree(valid,comm,good)
    end subroutine build_prepared_low_rate

    subroutine accept_held_prepared_source(stage,old_halos,raw_halos,rate,faces,registers,good)
      integer, intent(in) :: stage
      type(amr_tree_state_3d), intent(in) :: old_halos,raw_halos
      type(amr_tree_state_3d), intent(inout) :: rate
      type(eb_integrated_flux_3d), intent(inout) :: faces(:)
      type(state_redist_register_3d), allocatable, intent(inout) :: registers(:)
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: low
      type(eb_integrated_flux_3d), allocatable :: native_faces(:)
      integer :: p
      good=.true.
      if(spatial_weight<1) then
        call copy_mol_level_flux(tree,owners,rank,selected_depth,faces,native_faces,comm,good)
        if(.not.good) return
        call build_prepared_low_held_rate(raw_halos,rate,faces,low,good)
        if(.not.good) return
        call blend_mol_state_redist_stage_3d(tree,owners,geometry,selected_depth,1,spatial_weight,low,rate,registers,comm,good)
        if(.not.good) return
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          if(spatial_weight==0) then
            faces(p)=prepared_low_held_faces(p)
          else
            faces(p)%x=(1-spatial_weight)*prepared_low_held_faces(p)%x+spatial_weight*faces(p)%x
            faces(p)%y=(1-spatial_weight)*prepared_low_held_faces(p)%y+spatial_weight*faces(p)%y
            faces(p)%z=(1-spatial_weight)*prepared_low_held_faces(p)%z+spatial_weight*faces(p)%z
          endif
        enddo
      endif
      if(capture_prepared) then
        if(spatial_weight<1) then
          call capture_held_prepared_stage(stage,old_halos,raw_halos,rate,faces,registers,good,low,native_faces)
        else
          call capture_held_prepared_stage(stage,old_halos,raw_halos,rate,faces,registers,good)
        endif
      endif
    end subroutine accept_held_prepared_source

    subroutine build_prepared_low_held_rate(raw_halos,valid_rate,native_faces,low,good)
      type(amr_tree_state_3d), intent(in) :: raw_halos,valid_rate
      type(eb_integrated_flux_3d), intent(in) :: native_faces(:)
      type(amr_tree_state_3d), intent(out) :: low
      logical, intent(out) :: good
      real(dp) :: old_divergence,new_divergence,volume
      integer :: p,x,y,z
      logical :: valid
      call agree(allocated(prepared_low_held_faces),comm,good)
      if(.not.good) return
      call copy_field(valid_rate,low,good)
      if(.not.good) return
      valid=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        associate(g=>geometry(p),a=>native_faces(p),b=>prepared_low_held_faces(p))
        do z=1,g%nz
        do y=1,g%ny
        do x=1,g%nx
          low%fields(p)%values(:,x,y,z)=0
          if(g%volume_fraction(x,y,z)<=0) cycle
          volume=g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz
          old_divergence=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*a%x(1,x-1,y,z)- &
            g%x_face_fraction(x,y,z)*a%x(1,x,y,z))+ &
            g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*a%y(1,x,y-1,z)-g%y_face_fraction(x,y,z)*a%y(1,x,y,z))+ &
            g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*a%z(1,x,y,z-1)-g%z_face_fraction(x,y,z)*a%z(1,x,y,z))
          new_divergence=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*b%x(1,x-1,y,z)- &
            g%x_face_fraction(x,y,z)*b%x(1,x,y,z))+ &
            g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*b%y(1,x,y-1,z)-g%y_face_fraction(x,y,z)*b%y(1,x,y,z))+ &
            g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*b%z(1,x,y,z-1)-g%z_face_fraction(x,y,z)*b%z(1,x,y,z))
          low%fields(p)%values(1,x,y,z)=raw_halos%fields(p)%values(1,x+6,y+6,z+6)- &
            old_divergence/volume/dt+new_divergence/volume/dt
        enddo
        enddo
        enddo
        valid=valid.and.all(ieee_is_finite(low%fields(p)%values))
        end associate
      enddo
      call agree(valid,comm,good)
    end subroutine build_prepared_low_held_rate

    subroutine capture_held_prepared_stage(stage,old_halos,raw_halos,rate,faces,registers,good,low,native_faces)
      integer, intent(in) :: stage
      type(amr_tree_state_3d), intent(in) :: old_halos,raw_halos,rate
      type(eb_integrated_flux_3d), intent(in) :: faces(:)
      type(state_redist_register_3d), allocatable, intent(in) :: registers(:)
      logical, intent(out) :: good
      type(amr_tree_state_3d), optional, intent(in) :: low
      type(eb_integrated_flux_3d), optional, intent(in) :: native_faces(:)
      associate(record=>pending_history%prepared_stages(stage))
      record%has_internal_energy=.true.
      call copy_field(old_halos,record%held_state_halos,good,unit_bounds=.true.)
      if(.not.good) return
      call copy_field(raw_halos,record%held_raw_rate_halos,good,unit_bounds=.true.)
      if(.not.good) return
      call copy_field(rate,record%held_spatial_rate,good)
      if(.not.good) return
      call copy_mol_level_flux(tree,owners,rank,selected_depth,faces,record%held_faces,comm,good)
      if(.not.good) return
      call copy_registers(registers,record%held_registers,good)
      if(.not.good) return
      if(spatial_weight<1) then
        call agree(present(low).and.present(native_faces),comm,good)
        if(.not.good) return
        call copy_field(low,record%held_low_spatial_rate,good)
        if(.not.good) return
        call copy_mol_level_flux(tree,owners,rank,selected_depth,native_faces,record%held_native_faces,comm,good)
        if(.not.good) return
        call copy_mol_level_flux(tree,owners,rank,selected_depth,prepared_low_held_faces,record%held_low_faces,comm,good)
      endif
      end associate
    end subroutine capture_held_prepared_stage

    subroutine prepared_mol_sources(gas,temp,stage,output,flux,limiter,good,fractions)
      type(amr_tree_state_3d), allocatable :: transport_held
      type(amr_tree_state_3d), intent(in) :: gas,temp
      integer, intent(in) :: stage
      type(amr_tree_state_3d), intent(inout) :: output
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: limiter
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fractions(:)
      type(amr_tree_state_3d) :: source,old_halo,mol_source,mol_old,extra,low
      type(amr_tree_state_3d) :: scalar_old,scalar_raw
      type(state_redist_register_3d), allocatable :: prior_held_registers(:)
      type(eb_integrated_flux_3d), allocatable :: mol_flux(:),native_faces(:)
      type(state_redist_register_3d), allocatable :: registers(:)
      real(dp), allocatable :: body_exchange(:),frequency
      integer :: p,allocation_status,depth,first_depth,last_depth,m(3)
      logical :: local_ok
      good=.false.
      if(wall_active) then
        allocate(body_exchange(nv),frequency,stat=allocation_status)
        call agree(allocation_status==0,comm,good)
        if(.not.good) return
      endif
      if(present(internal_energy)) then
        allocate(transport_held)
        if(stage==0) then
          transport_held=held_old
        else
          transport_held=held_candidate
        endif
      endif
      call compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,tree,owners,geometry,extended_geometry, &
        gas,temp,dt,transport_flags,output,flux,source,old_halo,limiter,comm,good, &
        boundary_faces,inflow_state,inflow_temperature,wall_temperature,target_depth,ancestor_start,ancestor_end, &
        temperature_start,temperature_end,fractions,frame,embedded_wall_temperature,embedded_noslip, &
        body_exchange,frequency,inflow_mach,inflow_pressure, &
        constant_coefficients=constant_coefficients,internal_energy=transport_held, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end,coarse_linear=present(internal_energy), &
          primitive_from_total_energy=mol.and..not.present(frame).and..not.selected_momentum)
      if(.not.good) return
      if(wall_active) then
        if(stage==0) then
          heat_old=body_exchange(5);impulse_old=body_exchange(2:4)
        else
          heat_new=body_exchange(5);impulse_new=body_exchange(2:4)
        endif
      endif
      if(present(internal_energy)) then
        if(stage>1.and..not.chemistry) prior_held_registers=held_registers1
        if(stage==0) then
          call held_mol_source(gas,temp,held_old,held_rate0,held_flux0,held_registers0,good,fractions, &
            extra,mol_flux,mol_source,mol_old,scalar_old,scalar_raw)
        else
          call held_mol_source(gas,temp,held_candidate,held_rate1,held_flux1,held_registers1,good,fractions, &
            extra,mol_flux,mol_source,mol_old,scalar_old,scalar_raw)
        endif
      else
      call compute_mol_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,extended_geometry,gas,temp,dt,2, &
        extra,mol_flux,comm,good,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions,mol_lax_friedrichs,frame, &
        spatial_rate_halos=mol_source,old_state_halos=mol_old)
      endif
      if(.not.good) return
      local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        local_ok=local_ok.and.all(old_halo%fields(p)%values==mol_old%fields(p)%values)
        if(.not.local_ok) exit
        call guarded_spatial_add(output%fields(p)%values,extra%fields(p)%values,local_ok)
        if(.not.local_ok) exit
        call guarded_spatial_add(source%fields(p)%values,mol_source%fields(p)%values,local_ok)
        if(.not.local_ok) exit
        call guarded_spatial_add(flux(p)%x,mol_flux(p)%x,local_ok)
        if(.not.local_ok) exit
        call guarded_spatial_add(flux(p)%y,mol_flux(p)%y,local_ok)
        if(.not.local_ok) exit
        call guarded_spatial_add(flux(p)%z,mol_flux(p)%z,local_ok)
        if(.not.local_ok) exit
      enddo
      call agree(local_ok,comm,good)
      if(.not.good) return
      ! Match the serial prepared operator: molecular plus MOL, then SGS.
      ! The order is also used by archived/replayed raw rates and face fluxes.
      if(.not.(external_diffusion_source.and.stage>1)) then
        call add_spatial_diffusion_callback(gas,temp,extended_geometry,old_halo,output,source,flux,limiter,good,fractions, &
          transport_held)
        if(.not.good) return
      endif
      if(spatial_weight<1) then
        call build_prepared_low_rate(source,output,flux,low,good)
        if(.not.good) return
        if(capture_prepared) then
          call copy_mol_level_flux(tree,owners,rank,selected_depth,flux,native_faces,comm,good)
          if(.not.good) return
        endif
      endif
      first_depth=max(0,selected_depth);last_depth=first_depth
      if(selected_depth<0) last_depth=tree%level_count()-1
      do depth=first_depth,last_depth
        call compute_state_redist_mpi_level_stage_3d(tree,owners,extended_geometry,old_halo,source,nv, &
          depth,dt,output,registers,comm,good,boundary_faces,redist_target,redist_order,register_scale=.5_dp)
        if(.not.good) return
      enddo
      ! Post-redistribution density consistency precedes gas/EOS stepping.
      local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or..not.active(p)) cycle
        call post_eb_spatial_density_rate_3d(geometry(p),output%fields(p)%values,local_ok)
        if(.not.local_ok) exit
      enddo
      call agree(local_ok,comm,good);if(.not.good) return
      if(spatial_weight<1) then
        call blend_mol_state_redist_stage_3d(tree,owners,geometry,selected_depth,nv,spatial_weight, &
          low,output,registers,comm,good)
        if(.not.good) return
        do p=1,size(owners)
          if(owners(p)/=rank.or..not.active(p)) cycle
          if(spatial_weight==0) then
            flux(p)=prepared_low_faces(p)
          else
            flux(p)%x=(1-spatial_weight)*prepared_low_faces(p)%x+spatial_weight*flux(p)%x
            flux(p)%y=(1-spatial_weight)*prepared_low_faces(p)%y+spatial_weight*flux(p)%y
            flux(p)%z=(1-spatial_weight)*prepared_low_faces(p)%z+spatial_weight*flux(p)%z
          endif
        enddo
        theta=min(theta,spatial_weight)
      endif
      if(capture_prepared) then
        if(spatial_weight<1) then
          call capture_prepared_mol_stage(stage,old_halo,source,output,flux,registers,good,low,native_faces)
        else
          call capture_prepared_mol_stage(stage,old_halo,source,output,flux,registers,good)
        endif
        if(.not.good) return
      endif
      if(present(internal_energy)) then
        if(stage==0) then
          call accept_held_prepared_source(stage,scalar_old,scalar_raw,held_rate0,held_flux0,held_registers0,good)
        else
          call accept_held_prepared_source(stage,scalar_old,scalar_raw,held_rate1,held_flux1,held_registers1,good)
        endif
        if(.not.good) return
        if(stage>1.and..not.(chemistry.and.present(internal_energy))) then
          local_ok=.true.
          do p=1,size(owners)
            if(owners(p)/=rank.or..not.active(p)) cycle
            m=tree%nodes(p)%shape
            call accumulate_state_redist_coarse_register_3d(held_registers1(p),prior_held_registers(p),1,m+12,good)
            local_ok=local_ok.and.good
          enddo
          call agree(local_ok,comm,good)
          if(.not.good) return
        endif
      endif
      ! Provisional reacting passes do not occur in the final old/new mean.
      ! Their raw native ledgers were archived above. Publish only the accepted
      ! stages in the coupled gas/held family; the legacy gas policy retains
      ! direct coarse accumulation until separately qualified and migrated.
      if(stage==0) then
        call move_alloc(registers,stage_registers0)
      else
        if(stage>1.and..not.(chemistry.and.present(internal_energy))) then
          local_ok=.true.
          do p=1,size(owners)
            if(owners(p)/=rank.or..not.active(p)) cycle
            m=tree%nodes(p)%shape
            call accumulate_state_redist_coarse_register_3d(registers(p),stage_registers1(p),nv,m+12,good)
            if(.not.good) then
              local_ok=.false.;exit
            endif
          enddo
          call agree(local_ok,comm,good)
          if(.not.good) return
        endif
        call move_alloc(registers,stage_registers1)
      endif
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        extra%fields(p)%values=0
        if(.not.active(p)) cycle
        call wall_source(gas,temp,p,extra%fields(p)%values,local_ok,stage)
        if(.not.local_ok) exit
      enddo
      call agree(local_ok,comm,good)
      if(.not.good) return
      if(stage==0) then
        call integrate_mpi_eb_tree_3d(tree,owners,geometry,extra,mol_wall_old,comm,good)
        mol_wall_old=dt*mol_wall_old
      else
        call integrate_mpi_eb_tree_3d(tree,owners,geometry,extra,mol_wall_new,comm,good)
        mol_wall_new=dt*mol_wall_new
      endif
    end subroutine prepared_mol_sources
    subroutine wall_source(gas,temp,p,output,good,stage)
      type(amr_tree_state_3d), intent(in) :: gas,temp
      integer, intent(in) :: p,stage
      real(dp), intent(out) :: output(:,:,:,:)
      logical, intent(out) :: good
      real(dp), allocatable :: potential(:,:,:),stored(:,:,:)
      real(dp) :: position(3)
      integer :: x,y,z
      good=.false.;output=0
      associate(g=>geometry(p))
      if(present(frame)) then
        allocate(potential(g%nx,g%ny,g%nz),source=0._dp)
        do z=1,g%nz
        do y=1,g%ny
        do x=1,g%nx
          if(g%volume_fraction(x,y,z)<=0) cycle
          position=[g%x_lower,g%y_lower,g%z_lower]+(real([x,y,z],dp)-.5_dp)*[g%dx,g%dy,g%dz]
          call rotating_potential_3d(frame,position,potential(x,y,z),good)
          if(.not.good) return
        end do
        end do
        end do
      end if
      ! Match the actual gas divergence: ordinary MOL and Godunov use total
      ! energy. Rotating/experimental MOL retains its held-pressure policy.
      if(present(internal_energy).and.mol.and.(present(frame).or.selected_momentum)) then
        if(stage==0) then
          stored=held_old%fields(p)%values(1,:,:,:)
        else
          stored=held_candidate%fields(p)%values(1,:,:,:)
        endif
      endif
      call reactive_eb_slip_wall_source_3d(species,gas%fields(p)%values,temp%fields(p)%values(1,:,:,:), &
        g,output,good,potential,stored)
      end associate
    end subroutine wall_source

    subroutine recover(good,final_state,reset_phase)
      logical, intent(out) :: good
      logical, intent(in) :: final_state
      integer, optional, intent(in) :: reset_phase
      logical :: local_ok,cell_ok
      if(source_funded.and.selected_depth<0) then
        if(.not.present(energy_correction).and.present(reset_phase)) then
          if(reset_phase==1) then
            call held_update(held_old,held_rate0,held_rate0,.5_dp,held_candidate,good)
            if(.not.good) return
          else if(reset_phase==2.and..not.chemistry) then
            call held_update(held_old,held_rate0,held_rate1,.5_dp,held_candidate,good)
            if(.not.good) return
          endif
        endif
        ! These coarse values own no composite inventory and must represent
        ! the funded fine stage before they supply subsequent native sources.
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,good,frame,held_candidate)
        if(.not.good) return
      endif
      if(present(reset_phase).and.present(energy_correction)) then
        call reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,candidate,ct,held_candidate,reset_delta, &
          held_policy,selected_depth,reset_change,comm,good,frame, &
          internal_energy_change=scalar_reset_change)
        if(.not.good) return
        select case(reset_phase)
        case(1)
          predictor_reset=reset_delta;predictor_reset_change=reset_change
        case(2)
          final_reset=reset_change
        case(3)
          post_reset=post_reset+reset_change
          closure_inventory(1)=closure_inventory(1)+scalar_reset_change
          call agree(ieee_is_finite(closure_inventory(1)),comm,good)
          if(.not.good) return
        end select
        return
      endif
      local_ok=.true.
      do node=1,size(owners)
        if(owners(node)/=rank.or..not.active(node)) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          ! Primitive recovery tolerates roundoff negatives for evaluation;
          ! an accepted conservative state must also be writable as a restart.
          if(final_state.and.any(candidate%fields(node)%values(6:,i,j,k)<0)) then
            local_ok=.false.;cycle
          end if
          ! Native construct_Snew starts each nonreacting CTU update from
          ! S_old, including UTEMP. The pinned inverse can retain its seed.
          ! Limit this ownership correction to the total-energy pinned path.
          if(.not.mol.and..not.chemistry.and..not.present(internal_energy).and. &
             species(1)%pinned_eos_extrapolation.and.present(reset_phase)) then
            if(reset_phase==1.or.reset_phase==2) ct%fields(node)%values(1,i,j,k)=ot%fields(node)%values(1,i,j,k)
          endif
          call recover_reactive_eb_cell_3d(species,geometry(node),[i,j,k],candidate%fields(node)%values(:,i,j,k), &
            ct%fields(node)%values(1,i,j,k),primitive,recovered,speed,cell_ok,frame)
          if(.not.cell_ok) then
            local_ok=.false.;cycle
          end if
          ct%fields(node)%values(1,i,j,k)=recovered
        end do
        end do
        end do
      end do
      call agree(local_ok,comm,good)
    end subroutine
  end subroutine advance_mpi_eb_tree_sdc_3d
  subroutine reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,gas,temp,stored,delta, &
      policy,depth,change,comm,good,frame,leaf_only,internal_energy_change,reset_level,covered_only)
    ! -1 (default): reset normally; -2: validate on entry, repairing only
    ! failed held EOS states; >=0: reset that level and validate the others.
    ! Ordinary reacting MOL retains the fine-level total/held difference
    ! left by parent reflux until that level advances. A valid held EOS does
    ! not authorize overwriting conserved energy at another clock's boundary.
    integer, optional, intent(in) :: reset_level
    logical :: validate_first
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),depth
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(dual_energy_policy), intent(in) :: policy
    type(MPI_Comm), intent(in) :: comm
    type(rotating_frame_3d), optional, intent(in) :: frame
    logical, optional, intent(in) :: leaf_only,covered_only
    type(amr_tree_state_3d), intent(inout) :: gas,temp,stored
    type(amr_tree_state_3d), intent(out) :: delta
    real(dp), intent(out) :: change
    real(dp), optional, intent(inout) :: internal_energy_change
    logical, intent(out) :: good
    real(dp) :: cell_change,phi,position(3),volume,pv(reactive_nprim(size(species))),recovered_t,sound,local_change
    real(dp) :: held_before,held_cell_change,local_held_change,held_change
    integer :: p,d(3),x,y,z,rank,status,ios
    logical :: local_ok,cell_ok
    change=0;local_change=0;local_held_change=0
    good=.false.
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    allocate(delta%fields(size(tree%nodes)),stat=ios)
    call agree(ios==0,comm,good)
    if(.not.good) return
    local_ok=.true.
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      allocate(delta%fields(p)%values,source=gas%fields(p)%values,stat=ios)
      if(ios/=0) then
        local_ok=.false.;exit
      endif
    enddo
    call agree(local_ok,comm,good)
    if(.not.good) return
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      delta%fields(p)%values=0
      if(depth>=0.and.tree%nodes(p)%level/=depth) cycle
      d=tree%nodes(p)%shape
      associate(g=>geometry(p))
      volume=g%dx*g%dy*g%dz
      do z=1,d(3)
      do y=1,d(2)
      do x=1,d(1)
        if(g%volume_fraction(x,y,z)<=0) cycle
        if(present(leaf_only)) then
          if(leaf_only.and.covered(tree,p,[x,y,z])) cycle
        endif
        if(present(covered_only)) then
          if(covered_only.and..not.covered(tree,p,[x,y,z])) cycle
        endif
        phi=0
        if(present(frame)) then
          position=[g%x_lower,g%y_lower,g%z_lower]+(real([x,y,z],dp)-.5_dp)*[g%dx,g%dy,g%dz]
          call rotating_potential_3d(frame,position,phi,cell_ok)
          if(.not.cell_ok) then
            local_ok=.false.;cycle
          endif
        endif
        validate_first=.false.
        if(present(reset_level)) validate_first=reset_level==-2.or. &
          (reset_level>=0.and.tree%nodes(p)%level/=reset_level)
        if(validate_first) then
          call recover_reactive_eb_cell_3d(species,g,[x,y,z],gas%fields(p)%values(:,x,y,z), &
            temp%fields(p)%values(1,x,y,z),pv,recovered_t,sound,cell_ok,frame,stored%fields(p)%values(1,x,y,z))
          if(cell_ok) then
            ! Match serial admission: retain the native clock's EOS seed.
            if(reset_level/=-2) temp%fields(p)%values(1,x,y,z)=recovered_t
            cycle
          endif
          ! Fall through to the existing policy repair and ledger accounting.
        endif
        cell_change=0
        held_before=0
        if(present(internal_energy_change)) held_before=stored%fields(p)%values(1,x,y,z)
        call reset_nasa7_dual_energy_cell(species,gas%fields(p)%values(:,x,y,z), &
          stored%fields(p)%values(1,x,y,z),temp%fields(p)%values(1,x,y,z),policy,cell_change,cell_ok,phi)
        if(.not.cell_ok) then
          local_ok=.false.;cycle
        endif
        if(present(internal_energy_change)) then
          call compute_reactive_eb_held_projection(held_before,held_before,held_before, &
            stored%fields(p)%values(1,x,y,z),held_cell_change,cell_ok)
          if(.not.cell_ok) then
            local_ok=.false.;cycle
          endif
          if(.not.covered(tree,p,[x,y,z])) &
            local_held_change=local_held_change+volume*g%volume_fraction(x,y,z)*held_cell_change
        endif
        delta%fields(p)%values(5,x,y,z)=cell_change
        if(.not.covered(tree,p,[x,y,z])) local_change=local_change+volume*g%volume_fraction(x,y,z)*cell_change
        call recover_reactive_eb_cell_3d(species,g,[x,y,z],gas%fields(p)%values(:,x,y,z), &
          temp%fields(p)%values(1,x,y,z),pv,recovered_t,sound,cell_ok,frame,stored%fields(p)%values(1,x,y,z))
        if(.not.cell_ok) then
          local_ok=.false.;cycle
        endif
        temp%fields(p)%values(1,x,y,z)=recovered_t
      enddo
      enddo
      enddo
      end associate
    enddo
    call agree(local_ok.and.ieee_is_finite(local_change).and.ieee_is_finite(local_held_change),comm,good)
    if(.not.good) return
    call MPI_Allreduce(local_change,change,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    if(status/=MPI_SUCCESS) then
      good=.false.;return
    endif
    call agree(ieee_is_finite(change),comm,good)
    if(.not.good) return
    if(present(internal_energy_change)) then
      call MPI_Allreduce(local_held_change,held_change,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
      if(status/=MPI_SUCCESS) then
        good=.false.;return
      endif
      call agree(ieee_is_finite(held_change),comm,good)
      if(.not.good) return
      internal_energy_change=held_change
    endif
  end subroutine reset_held_mpi_eb_tree_family_3d

  ! Shared paired closure for common and recursive held prepared clocks.
  subroutine close_held_prepared_mpi_eb_tree_3d(species,tree,owners,geometry,gas,temp,stored,faces,scalar_faces, &
      policy,prefix,depth,comm,good,boundary_faces,frame,register_geometry,registers,scalar_registers, &
      internal_energy_change,reset_current_level)
    logical, optional, intent(in) :: reset_current_level
    integer :: selected_reset_level
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(eb_geometry_3d), optional, intent(in) :: register_geometry(:)
    type(dual_energy_policy), intent(in) :: policy
    type(MPI_Comm), intent(in) :: comm
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), intent(inout) :: gas,temp,stored
    type(eb_integrated_flux_3d), intent(inout) :: faces(:),scalar_faces(:)
    real(dp), intent(inout) :: prefix
    real(dp), optional, intent(inout) :: internal_energy_change
    integer, intent(in) :: depth
    logical, intent(out) :: good
    type(state_redist_register_3d), optional, intent(inout) :: registers(:),scalar_registers(:)
    type(amr_tree_state_3d) :: candidate,temps,held,increment,correction,spill,reset_delta,scalar_correction,scalar_spill
    type(eb_integrated_flux_3d), allocatable :: effective(:),scalar_effective(:),repair_faces(:),scalar_input(:)
    type(state_redist_register_3d), allocatable :: scalar_work(:)
    type(amr_tree_state_3d) :: held_before,joint_correction,joint_spill
    type(state_redist_register_3d), allocatable :: flux_registers(:)
    real(dp), allocatable :: joint_source(:,:,:,:)
    integer :: h
    logical :: flux_rereflux
    real(dp) :: reset_change,trial_prefix,value,delta
    real(dp) :: held_reset_change,held_prefix
    integer :: first_depth,last_depth,scalar_first_depth,p,n(3),x,y,z,c,nv,rank,status,closure_attempt
    logical :: valid,agreed,has_spill,has_scalar_spill
    good=.false.;nv=reactive_nvar(size(species))
    selected_reset_level=-1
    if(present(reset_current_level)) then
      if(reset_current_level) selected_reset_level=depth
    endif
    call agree((present(registers).eqv.present(scalar_registers)).and. &
      (.not.present(registers).or.present(register_geometry)).and.ieee_is_finite(prefix),comm,agreed)
    if(.not.agreed) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    candidate=gas;temps=temp;held=stored;repair_faces=faces;scalar_input=scalar_faces
    first_depth=depth;last_depth=min(depth+1,tree%level_count()-1)
    if(present(internal_energy_change)) then
      call agree_held_projection(internal_energy_change,comm,valid)
      if(.not.valid) return
    endif
    flux_rereflux=present(register_geometry).and..not.present(registers).and.depth+1<tree%level_count()
    if(flux_rereflux) allocate(flux_registers(size(tree%nodes)))
    held_reset_change=0
    do closure_attempt=1,2
    first_depth=depth
    scalar_effective=scalar_input
    held_before=held
    if(present(registers)) then
      call assemble_state_redist_mpi_interface_rereflux_3d(tree,owners,register_geometry,registers,nv, &
        depth,correction,comm,valid,boundary_faces,spill,has_spill)
      if(.not.valid) return
      if(has_spill) then
        first_depth=depth-1
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          if(tree%nodes(p)%level==first_depth) correction%fields(p)=spill%fields(p)
        enddo
      endif
      call assemble_state_redist_mpi_interface_rereflux_3d(tree,owners,register_geometry,scalar_registers,1, &
        depth,scalar_correction,comm,valid,boundary_faces,scalar_spill,has_scalar_spill)
      if(.not.valid) return
      scalar_first_depth=depth
      if(has_scalar_spill) then
        scalar_first_depth=depth-1
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          if(tree%nodes(p)%level==scalar_first_depth) scalar_correction%fields(p)=scalar_spill%fields(p)
        enddo
      endif
      first_depth=min(first_depth,scalar_first_depth)
      call assemble_mpi_eb_tree_reflux_3d(species,tree,owners,geometry,candidate,temps,repair_faces,increment, &
        effective,comm,valid,depth,last_depth,boundary_faces,frame,correction,first_depth)
      if(.not.valid) return
      scalar_work=scalar_registers
      call synchronize_mpi_eb_tree_scalar_register_flux_3d(species,tree,owners,geometry,register_geometry, &
        held,scalar_effective,scalar_work,depth,comm,valid,boundary_faces)
    else
      call assemble_mpi_eb_tree_reflux_3d(species,tree,owners,geometry,candidate,temps,repair_faces,increment, &
        effective,comm,valid,depth,last_depth,boundary_faces,frame)
      if(.not.valid) return
      call synchronize_mpi_eb_tree_scalar_flux_3d(species,tree,owners,geometry,held,scalar_effective, &
        comm,valid,minimum_level=depth,maximum_level=last_depth,boundary_faces=boundary_faces)
    endif
    if(.not.valid) return
    if(flux_rereflux) then
      valid=.true.
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        n=tree%nodes(p)%shape
        allocate(joint_source(nv+1,n(1),n(2),n(3)))
        joint_source(:nv,:,:,:)=increment%fields(p)%values
        joint_source(nv+1:nv+1,:,:,:)=held%fields(p)%values-held_before%fields(p)%values
        call initialize_flux_rereflux_register_3d(tree,p,register_geometry(p),joint_source, &
          flux_registers(p),agreed,boundary_faces)
        valid=valid.and.agreed
        deallocate(joint_source)
      enddo
      call agree(valid,comm,agreed)
      if(.not.agreed) return
      call assemble_state_redist_mpi_interface_rereflux_3d(tree,owners,register_geometry,flux_registers,nv+1, &
        depth,joint_correction,comm,valid,boundary_faces,joint_spill,has_spill)
      if(.not.valid) return
      first_depth=depth
      if(has_spill) first_depth=depth-1
      correction=increment;scalar_correction=held_before
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        n=tree%nodes(p)%shape
        if(tree%nodes(p)%level==depth) then
          h=(register_geometry(p)%nx-n(1))/2
          joint_correction%fields(p)%values=joint_correction%fields(p)%values- &
            flux_registers(p)%coarse_delta(:,h+1:h+n(1),h+1:h+n(2),h+1:h+n(3))
        else if(has_spill.and.tree%nodes(p)%level==first_depth) then
          joint_correction%fields(p)=joint_spill%fields(p)
        endif
        correction%fields(p)%values=joint_correction%fields(p)%values(:nv,:,:,:)
        scalar_correction%fields(p)%values=joint_correction%fields(p)%values(nv+1:nv+1,:,:,:)
      enddo
      ! The existing owner-local route handles fine-covered and spilled
      ! recipients. Publish gas and UEINT only after the paired EOS admission.
      call assemble_mpi_eb_tree_reflux_3d(species,tree,owners,geometry,candidate,temps,repair_faces,increment, &
        effective,comm,valid,depth,last_depth,boundary_faces,frame,correction,first_depth)
      if(.not.valid) return
      held=held_before;scalar_effective=scalar_input
      call synchronize_mpi_eb_tree_scalar_flux_3d(species,tree,owners,geometry,held,scalar_effective, &
        comm,valid,minimum_level=depth,maximum_level=last_depth,boundary_faces=boundary_faces, &
        rereflux_delta=scalar_correction,rereflux_first_level=first_depth)
      if(.not.valid) return
    endif
    valid=.true.
    do p=1,size(tree%nodes)
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        do c=1,nv
          value=candidate%fields(p)%values(c,x,y,z);delta=increment%fields(p)%values(c,x,y,z)
          if(.not.ieee_is_finite(value).or..not.ieee_is_finite(delta)) then
            valid=.false.;cycle
          endif
          if(delta>0) then
            if(value>huge(1._dp)-delta) then
              valid=.false.;cycle
            endif
          else if(delta<0) then
            if(value< -huge(1._dp)-delta) then
              valid=.false.;cycle
            endif
          endif
          candidate%fields(p)%values(c,x,y,z)=value+delta
        enddo
        if(any(candidate%fields(p)%values(6:,x,y,z)<0)) valid=.false.
      enddo
      enddo
      enddo
    enddo
    call agree(valid,comm,agreed)
    if(agreed) then
      call reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,candidate,temps,held,reset_delta, &
        policy,-1,reset_change,comm,valid,frame,leaf_only=.true., &
        internal_energy_change=held_reset_change,reset_level=selected_reset_level)
      agreed=valid
    endif
    if(agreed) exit
    if(closure_attempt==2) return
    ! A represented species or actual-UEINT EOS deficit needs conservative funding.
    ! Restart from the private pre-reflux family and repair gas and actual
    ! UEINT together; no particle/source callbacks or scalar projection occur.
    candidate=gas;temps=temp;held=stored;repair_faces=faces;scalar_input=scalar_faces
    if(present(registers).or.flux_rereflux) then
      call repair_mpi_eb_tree_held_interface_budget_3d(species,tree,owners,geometry,candidate,temps,held, &
        repair_faces,scalar_input,comm,valid,depth,boundary_faces,frame,correction,scalar_correction, &
        rereflux_first_level=first_depth)
    else
    call repair_mpi_eb_tree_held_interface_budget_3d(species,tree,owners,geometry,candidate,temps,held, &
      repair_faces,scalar_input,comm,valid,depth,boundary_faces,frame)
    endif
    if(.not.valid) return
    enddo
    ! Only the accepted private reset enters the prefix. Failed reset trials
    ! were discarded with the pre-reflux family before conservative repair.
    held_prefix=held_reset_change
    if(present(internal_energy_change)) held_prefix=internal_energy_change+held_reset_change
    call agree(ieee_is_finite(held_prefix),comm,valid)
    if(.not.valid) return
    trial_prefix=prefix+reset_change
    call agree(ieee_is_finite(trial_prefix),comm,valid)
    if(.not.valid) return
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,temps,comm,valid,frame,held)
    if(.not.valid) return
    if(selected_reset_level>=0) then
      ! Covered parent cells supply future interpolation but carry no
      ! composite inventory. Repair them only after the closing average-down.
      call reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,candidate,temps,held,reset_delta, &
        policy,selected_reset_level,reset_change,comm,valid,frame,covered_only=.true.)
      if(.not.valid.or.reset_change/=0) return
    endif
    call move_alloc(candidate%fields,gas%fields);call move_alloc(temps%fields,temp%fields)
    call move_alloc(held%fields,stored%fields)
    faces=effective;scalar_faces=scalar_effective;prefix=trial_prefix
    if(present(registers)) then
      scalar_registers=scalar_work
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        if(tree%nodes(p)%level==depth) then
          if(allocated(registers(p)%coarse_delta)) registers(p)%coarse_delta=0
        else if(tree%nodes(p)%level==depth+1) then
          if(allocated(registers(p)%fine_delta)) registers(p)%fine_delta=0
        endif
      enddo
    endif
    if(present(internal_energy_change)) internal_energy_change=held_prefix
    good=.true.
  end subroutine close_held_prepared_mpi_eb_tree_3d

  recursive subroutine advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,geometry,state,temperature, &
      dt,iterations,chemistry, &
      transport_flags,force,rtol,atol,extended_geometry,boundary_faces, &
      boundary_change,wall_change,chemical_change,source_change,minimum_theta,comm,ok,reaction_source, &
      difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature,level_steps,redistribution, &
      target_volume_fraction,auxiliary_source,auxiliary_change,route_auxiliary,gas_acceleration,halo_source, &
      predictor_builder,post_level,evaluate_diffusion, &
      evaluate_reaction,embedded_wall_temperature,embedded_noslip,mol_mode,mol_lax_friedrichs,frame, &
      state_redist_max_order,evaluate_spatial_diffusion,spatial_diffusion_halo_width, &
      constant_coefficients,heat_release,registered_flux,integrated_flux,internal_energy,internal_energy_flux, &
      interval_history,energy_policy,energy_correction, &
      spatial_diffusion_is_source,auxiliary_has_internal_energy,internal_energy_pressure_work_change, &
      godunov_interval_history,godunov_replay_history,spatial_flux_weight,provider_refused, &
      godunov_parent_registers,godunov_parent_scalar_registers,godunov_routed_pending,godunov_replay_prefix, &
      godunov_recovery_history,godunov_recovery_ready,source_funding_context,native_replay_history, &
      native_spatial_weight,native_common_retries,native_replay_prefix,native_replay_preserve_weights, &
      separate_momentum,consistent_velocity,paired_pressure,connected_directions, &
      evaluate_reaction_admission,reaction_admission_status)
    procedure(forced_reactor_admission), optional :: evaluate_reaction_admission
    integer, optional, intent(out) :: reaction_admission_status
    integer :: admission_state,admission_global,admission_low,admission_high,admission_ierr,child_admission
    logical, optional, intent(in) :: separate_momentum,consistent_velocity,paired_pressure,connected_directions
    logical :: selected_momentum,selected_velocity,selected_pressure,selected_connected
    ! Failure diagnostics; these do not commit a physical step or ledger.
    integer(int64), optional, intent(out) :: native_common_retries
    ! The full native receipt/open-face context is private with gas/H until
    ! outer success. Non-regular prepared funding still requires its native
    ! joint StateRedist stage path; this connects ordinary held MOL clocks.
    ! Complete source-frozen native MOL window replay; all fields/receipts and
    ! ledgers stay private until the last descendant and closure accept.
    type(mpi_eb_mol_interval_3d), optional, intent(in) :: native_replay_history(:)
    real(dp), optional, intent(in) :: native_spatial_weight
    real(dp) :: native_mol_weight,native_mol_weight_low,native_mol_weight_high
    ! A prefix reuses only the supplied complete clocks; later clocks evaluate
    ! their sources natively. Publication still waits for the entire window.
    logical, optional, intent(in) :: native_replay_prefix
    logical, optional, intent(in) :: native_replay_preserve_weights
    logical :: replaying_native,native_prefix,preserve_native_weights
    integer :: native_cursor
    integer :: native_replay_depth
    real(dp) :: native_replay_time
    type(eb_held_funding_context_3d), optional, intent(inout) :: source_funding_context
    type(eb_held_funding_context_3d), allocatable :: funding_context_work
    logical, allocatable :: funding_local(:)
    integer :: funding_status
    integer :: funding_controls(2),funding_low(2),funding_high(2)
    ! Replay actual native clocks from fresh caller state, freezing declared external packets.
    type(eb_godunov_interval_3d), optional, intent(in) :: godunov_replay_history(:)
    real(dp), optional, intent(in) :: spatial_flux_weight
    ! A prefix contains only accepted source tables. Later clocks run natively
    ! for the first time; no absent source packet is synthesized.
    logical, optional, intent(in) :: godunov_replay_prefix
    ! Diagnostic-only channel for a failed transaction: complete accepted
    ! common source tables, without publishing any failed physical result.
    type(eb_godunov_interval_3d), allocatable, optional, intent(inout) :: godunov_recovery_history(:)
    logical, optional, intent(out) :: godunov_recovery_ready
    logical :: godunov_prefix
    integer :: godunov_replay_count
    integer :: godunov_replay_cursor
    ! Nonroot replay restores the accepted clock and returns the still-open
    ! parent ledgers and routed buckets only with the complete new subtree.
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: godunov_parent_registers(:)
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: godunov_parent_scalar_registers(:)
    type(amr_tree_state_3d), allocatable, optional, intent(inout) :: godunov_routed_pending(:)
    type(state_redist_register_3d), allocatable :: godunov_parent_work(:),godunov_parent_scalar_work(:)
    type(amr_tree_state_3d), allocatable :: godunov_pending_work(:)
    type(state_redist_register_3d) :: godunov_expected_register
    integer :: godunov_replay_depth,godunov_replay_base,godunov_parent
    real(dp) :: godunov_replay_time,godunov_pressure_offset,godunov_minimum_offset
    integer :: godunov_header(3),godunov_header_low(3),godunov_header_high(3)
    real(dp) :: godunov_clock_values(5),godunov_clock_low(5),godunov_clock_high(5)

    real(dp) :: godunov_weight,godunov_weight_low,godunov_weight_high
    logical, optional, intent(out) :: provider_refused
    type(eb_godunov_interval_3d), allocatable, optional, intent(inout) :: godunov_interval_history(:)
    type(eb_godunov_interval_3d), allocatable :: godunov_records(:)
    integer, allocatable :: godunov_active(:)
    integer :: godunov_cursor,godunov_capacity,godunov_level_calls,godunov_status
    logical :: godunov_capture,godunov_mol
    logical :: godunov_automatic_repair,godunov_provider_rejected,godunov_source_exhausted
    type :: godunov_root_trial
      type(amr_tree_state_3d) :: gas,temperature
      type(amr_tree_state_3d) :: stored,reaction
      type(state_redist_register_3d), allocatable :: parent_registers(:),parent_scalar_registers(:)
      type(amr_tree_state_3d), allocatable :: routed_pending(:)
      type(amr_tree_state_3d), allocatable :: heat
      type(eb_godunov_interval_3d), allocatable :: records(:),recovery(:)
      logical :: recovery_ready=.false.
      type(eb_integrated_flux_3d), allocatable :: faces(:),registered_faces(:),scalar_faces(:)
      real(dp), allocatable :: boundary(:),wall(:),chemical(:),external(:),auxiliary(:)
      integer(int64), allocatable :: counts(:)
      real(dp) :: prefix=0,pressure_work=0,limiter=1
      logical :: refused=.false.
    end type godunov_root_trial

    ! Composite signed pressure work over accepted native level intervals.
    ! Publish only with the outer tree/held transaction.
    real(dp), optional, intent(inout) :: internal_energy_pressure_work_change
    real(dp) :: pressure_work_change,initial_clock_reset
    logical, optional, intent(in) :: auxiliary_has_internal_energy
    logical :: held_auxiliary
    integer :: naux,held_aux_low,held_aux_high,auxiliary_alloc_status
    logical, optional, intent(in) :: spatial_diffusion_is_source
    logical :: external_diffusion_source
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    real(dp), optional, intent(inout) :: energy_correction
    type(dual_energy_policy) :: held_policy
    real(dp), allocatable :: correction_work
    real(dp) :: reset_change,held_controls(2),held_controls_low(2),held_controls_high(2)
    integer :: held_policy_flags(3),held_policy_min(3),held_policy_max(3)
    type(amr_tree_state_3d) :: reset_difference
    type(amr_tree_state_3d), optional, intent(inout) :: heat_release
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: internal_energy_flux(:)
    ! Preorder interval history is published only with the accepted outer step.
    type(mpi_eb_mol_interval_3d), allocatable, optional, intent(inout) :: interval_history(:)
    type(mpi_eb_mol_interval_3d), allocatable :: history(:)
    integer :: history_count,history_size,history_level_count
    ! Time-integrated spatial face flux with pinned stage registration weights,
    ! before peer/coarse-fine synchronization. Separate from the physical flux
    ! used by the conservative state update and boundary inventory.
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: registered_flux(:)
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: integrated_flux(:)
    type(eb_integrated_flux_3d), allocatable :: registration(:)
    type(amr_tree_state_3d), allocatable :: heat
    real(dp), optional, intent(in) :: constant_coefficients(:)
    logical :: constant_agreed
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(rotating_frame_3d), optional, intent(in) :: frame
    procedure(mpi_eb_sdc_diffusion_source), optional :: evaluate_diffusion
    procedure(mpi_eb_sdc_spatial_diffusion_source), optional :: evaluate_spatial_diffusion
    integer, optional, intent(in) :: spatial_diffusion_halo_width
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: dt,force(3),rtol,atol
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry,transport_flags(4)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(MPI_Comm), intent(in) :: comm
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),minimum_theta
    integer, optional, intent(in) :: state_redist_max_order
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: wall_temperature(6),gas_acceleration(3),embedded_wall_temperature
    logical, optional, intent(in) :: embedded_noslip,mol_mode,mol_lax_friedrichs
    integer(int64), optional, intent(inout) :: level_steps(:)
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(in) :: target_volume_fraction
    procedure(mpi_eb_sdc_recursive_source), optional :: auxiliary_source
    procedure(mpi_eb_sdc_recursive_halo_source), optional :: halo_source
    procedure(mpi_eb_sdc_recursive_predictor_builder), optional :: predictor_builder
    procedure(mpi_eb_sdc_post_level), optional :: post_level
    logical :: has_auxiliary
    real(dp), optional, intent(inout) :: auxiliary_change(:)
    real(dp), allocatable :: auxiliary(:)
    real(dp) :: source_time
    logical, optional, intent(in) :: route_auxiliary
    logical :: routing,old_auxiliary_only
    type(amr_tree_state_3d) :: off_old,off_new
    type(amr_tree_state_3d), allocatable :: auxiliary_pending(:)
    integer :: redist_order,redist_low,redist_high
    real(dp) :: redist_target
    character(len=16) :: method
    integer :: method_code,method_low,method_high
    type :: interval_register_set
      type(state_redist_register_3d), allocatable :: fields(:)
    end type interval_register_set
    type(interval_register_set), allocatable :: interface_registers(:),held_interface_registers(:)
    logical :: prepared_redist,prepared_godunov
    integer :: godunov_padding
    type(eb_geometry_3d), allocatable :: register_geometry(:)
    type(amr_tree_state_3d), allocatable :: pending(:)
    type(amr_tree_state_3d) :: candidate,ct,starts,ends,ts,te
    type(amr_tree_state_3d) :: held_candidate,held_starts,held_ends
    type(amr_tree_state_3d), allocatable :: ir
    type(amr_tree_state_3d) :: source0,source1
    type(eb_integrated_flux_3d), allocatable :: total(:),synchronization_flux(:),held_total(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:)
    integer, allocatable :: ratios(:)
    ! Only completed physics closures may request a larger clock-group replay.
    ! Provider/callback failures leave these slots empty and roll back directly.
    integer, allocatable :: deferred_record(:),deferred_iteration(:),deferred_cycle(:)
    real(dp), allocatable :: clock_start(:),clock_dt(:)
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: external(size(source_change)),theta
    integer :: nv,node,n(3),levels,c,rank,status,ios
    logical :: valid,ready,replay_enabled,collect_history
    type :: replay_flux_family
      type(eb_integrated_flux_3d), allocatable :: faces(:),held_faces(:)
    end type replay_flux_family
    type :: mol_replay_result
      type(amr_tree_state_3d) :: state,temperature,reaction_source,internal_energy
      type(eb_held_funding_context_3d), allocatable :: funding
      real(dp) :: energy_correction=0
      type(amr_tree_state_3d) :: starts,ends,ts,te,source0,source1,held_starts,held_ends
      type(amr_tree_state_3d), allocatable :: redistribution_pending(:)
      type(amr_tree_state_3d), allocatable :: heat
      type(eb_integrated_flux_3d), allocatable :: flux(:),held_flux(:)
      type(mpi_eb_mol_interval_3d), allocatable :: records(:)
      real(dp), allocatable :: chemical(:),wall(:),source(:),auxiliary(:)
      type(interval_register_set), allocatable :: registers(:),held_registers(:)
    end type mol_replay_result
    type :: prepared_replay_work
      type(amr_tree_state_3d) :: gas,temp,reaction,starts,ends,ts,te,source0,source1,held,held_starts,held_ends
      type(amr_tree_state_3d), allocatable :: heat,routed(:)
      type(interval_register_set), allocatable :: registers(:),held_registers(:)
      type(mpi_eb_mol_interval_3d), allocatable :: records(:)
      real(dp), allocatable :: clock_start(:),clock_dt(:)
      integer :: cursor=0,first=0,last=0,depth=0
      real(dp) :: weight=1,time=0,energy_correction=0
    end type prepared_replay_work
    type(prepared_replay_work), pointer :: prepared_context=>null()
    logical :: prepared_provider_rejected=.false.
    integer :: heat_supplied,heat_min,heat_max,heat_status
    if(present(provider_refused)) provider_refused=.false.
    admission_state=reaction_admission_not_attempted
    if(present(reaction_admission_status)) then
      reaction_admission_status=admission_state
      if(present(evaluate_reaction_admission)) reaction_admission_status=reaction_admission_unknown
    endif
    ok=.false.
    call MPI_Allreduce(merge(1,0,present(evaluate_reaction_admission)),admission_low,1,MPI_INTEGER,MPI_MIN,comm,admission_ierr)
    if(admission_ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,present(evaluate_reaction_admission)),admission_high,1,MPI_INTEGER,MPI_MAX,comm,admission_ierr)
    if(admission_ierr/=MPI_SUCCESS.or.admission_low/=admission_high) then
      if(present(reaction_admission_status)) reaction_admission_status=reaction_admission_unknown
      return
    endif
    call agree(.not.present(evaluate_reaction_admission).or. &
      (chemistry.and.present(evaluate_reaction)),comm,ready)
    if(.not.ready) then
      if(present(reaction_admission_status)) reaction_admission_status=reaction_admission_unknown
      return
    endif
    call run_recursive()
    call MPI_Allreduce(admission_state,admission_global,1,MPI_INTEGER,MPI_MAX,comm,admission_ierr)
    if(admission_ierr/=MPI_SUCCESS) then
      ok=.false.;admission_global=reaction_admission_unknown
    endif
    admission_state=admission_global
    if(.not.ok.and.present(evaluate_reaction_admission).and.admission_state==reaction_admission_not_attempted) &
      admission_state=reaction_admission_unknown
    if(ok.and.admission_state/=reaction_admission_not_attempted.and. &
      admission_state/=reaction_admission_provider_refused) admission_state=reaction_admission_accepted
    if(present(reaction_admission_status)) reaction_admission_status=admission_state
    if(admission_state==reaction_admission_provider_refused) then
      ok=.false.
      if(present(provider_refused)) provider_refused=.true.
    endif
  contains
    subroutine run_recursive()
    if(present(native_common_retries)) native_common_retries=0_int64
    if(present(godunov_recovery_ready)) godunov_recovery_ready=.false.
    held_auxiliary=.false.
    if(present(auxiliary_has_internal_energy)) held_auxiliary=auxiliary_has_internal_energy
    external_diffusion_source=.false.
    if(present(spatial_diffusion_is_source)) external_diffusion_source=spatial_diffusion_is_source
    ok=.false.;nv=reactive_nvar(size(species))
    initial_clock_reset=0
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return

    godunov_capture=present(godunov_interval_history).or.present(godunov_recovery_history)
    godunov_cursor=0;godunov_replay_cursor=0
    godunov_replay_depth=0;godunov_replay_base=0;godunov_parent=0
    godunov_replay_time=0;godunov_pressure_offset=0;godunov_minimum_offset=1
    godunov_provider_rejected=.false.;godunov_source_exhausted=.false.
    godunov_prefix=.false.
    if(present(godunov_replay_prefix)) godunov_prefix=godunov_replay_prefix
    call MPI_Allreduce(merge(1,0,godunov_prefix),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,godunov_prefix),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call agree(.not.godunov_prefix.or.present(godunov_replay_history),comm,ready)
    if(.not.ready) return
    call agree(present(godunov_recovery_history).eqv.present(godunov_recovery_ready),comm,ready)
    if(.not.ready) return
    godunov_weight=1
    if(present(spatial_flux_weight)) godunov_weight=spatial_flux_weight
    if(present(provider_refused)) provider_refused=.false.
    replaying_native=present(native_replay_history);native_cursor=0
    native_prefix=.false.
    if(present(native_replay_prefix)) native_prefix=native_replay_prefix
    preserve_native_weights=.false.
    if(present(native_replay_preserve_weights)) preserve_native_weights=native_replay_preserve_weights
    call MPI_Allreduce(merge(1,0,preserve_native_weights),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,preserve_native_weights),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call agree(.not.preserve_native_weights.or.(replaying_native.and..not.present(native_spatial_weight)),comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(merge(1,0,native_prefix),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,native_prefix),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call agree(.not.native_prefix.or.replaying_native,comm,ready);if(.not.ready) return
    native_replay_depth=0;native_replay_time=0
    godunov_mol=.false.
    if(present(mol_mode)) godunov_mol=mol_mode
    selected_momentum=.false.;selected_velocity=.false.;selected_pressure=.false.;selected_connected=.false.
    if(present(separate_momentum)) selected_momentum=separate_momentum
    if(present(consistent_velocity)) selected_velocity=consistent_velocity
    if(present(paired_pressure)) selected_pressure=paired_pressure
    if(present(connected_directions)) selected_connected=connected_directions
    call MPI_Allreduce(merge(1,0,selected_momentum)+merge(2,0,selected_velocity)+merge(4,0, &
      selected_pressure)+merge(8,0,selected_connected),method_low, &
      1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,selected_momentum)+merge(2,0,selected_velocity)+merge(4,0, &
      selected_pressure)+merge(8,0,selected_connected),method_high, &
      1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call agree((.not.selected_pressure.or.selected_velocity).and. &
      (.not.selected_pressure.or..not.present(frame)),comm,ready)
    if(.not.ready) return
    call agree((.not.selected_velocity.or.selected_momentum).and. &
      (.not.selected_connected.or.selected_momentum),comm,ready)
    if(.not.ready) return
    naux=eb_auxiliary_nvar(nv,auxiliary_has_internal_energy)
    allocate(auxiliary(naux),stat=auxiliary_alloc_status)
    call agree(auxiliary_alloc_status==0,comm,ready);if(.not.ready) return
    call MPI_Allreduce(merge(1,0,held_auxiliary),held_aux_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,held_auxiliary),held_aux_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.held_aux_low/=held_aux_high) return
    call agree(.not.held_auxiliary.or.(present(internal_energy).and. &
      (present(auxiliary_source).or.present(halo_source))),comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(merge(1,0,external_diffusion_source),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,external_diffusion_source),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call agree(.not.external_diffusion_source.or.present(evaluate_spatial_diffusion).or. &
      present(evaluate_diffusion),comm,ready)
    if(.not.ready) return
    held_policy=dual_energy_policy()
    if(present(energy_policy)) held_policy=energy_policy
    heat_supplied=merge(1,0,present(heat_release))
    call MPI_Allreduce(heat_supplied,heat_min,1,MPI_INTEGER,MPI_MIN,comm,heat_status)
    if(heat_status/=MPI_SUCCESS) return
    call MPI_Allreduce(heat_supplied,heat_max,1,MPI_INTEGER,MPI_MAX,comm,heat_status)
    if(heat_status/=MPI_SUCCESS.or.heat_min/=heat_max) return
    call MPI_Allreduce(merge(1,0,present(internal_energy))+2*merge(1,0,present(internal_energy_flux))+ &
      16*merge(1,0,present(energy_correction))+32*merge(1,0,present(energy_policy))+ &
      64*merge(1,0,present(internal_energy_pressure_work_change))+128*merge(1,0,present(godunov_interval_history))+ &
      256*merge(1,0,present(godunov_replay_history))+512*merge(1,0,present(spatial_flux_weight))+ &
      1024*merge(1,0,present(provider_refused))+2048*merge(1,0,present(godunov_parent_registers))+ &
      4096*merge(1,0,present(godunov_parent_scalar_registers))+8192*merge(1,0,present(godunov_routed_pending))+ &
      16384*merge(1,0,present(godunov_recovery_history))+32768*merge(1,0,present(godunov_recovery_ready))+ &
      65536*merge(1,0,replaying_native)+131072*merge(1,0,present(native_spatial_weight))+ &
      262144*merge(1,0,present(native_common_retries)), &
      method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,present(internal_energy))+2*merge(1,0,present(internal_energy_flux))+ &
      16*merge(1,0,present(energy_correction))+32*merge(1,0,present(energy_policy))+ &
      64*merge(1,0,present(internal_energy_pressure_work_change))+128*merge(1,0,present(godunov_interval_history))+ &
      256*merge(1,0,present(godunov_replay_history))+512*merge(1,0,present(spatial_flux_weight))+ &
      1024*merge(1,0,present(provider_refused))+2048*merge(1,0,present(godunov_parent_registers))+ &
      4096*merge(1,0,present(godunov_parent_scalar_registers))+8192*merge(1,0,present(godunov_routed_pending))+ &
      16384*merge(1,0,present(godunov_recovery_history))+32768*merge(1,0,present(godunov_recovery_ready))+ &
      65536*merge(1,0,replaying_native)+131072*merge(1,0,present(native_spatial_weight))+ &
      262144*merge(1,0,present(native_common_retries)), &
      method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    if(iand(method_low,3)/=0.and.iand(method_low,3)/=3) return
    native_mol_weight=1
    if(present(native_spatial_weight)) native_mol_weight=native_spatial_weight
    call agree(ieee_is_finite(native_mol_weight),comm,ready);if(.not.ready) return
    call agree(native_mol_weight>=0.and.native_mol_weight<=1.and. &
      (.not.present(native_spatial_weight).or.replaying_native),comm,ready);if(.not.ready) return
    call MPI_Allreduce(native_mol_weight,native_mol_weight_low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(native_mol_weight,native_mol_weight_high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.native_mol_weight_low/=native_mol_weight_high) return
    call agree(ieee_is_finite(godunov_weight),comm,ready);if(.not.ready) return
    call agree(godunov_weight>=0.and.godunov_weight<=1,comm,ready);if(.not.ready) return
    call MPI_Allreduce(godunov_weight,godunov_weight_low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(godunov_weight,godunov_weight_high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.godunov_weight_low/=godunov_weight_high) return
    valid=(.not.present(internal_energy_pressure_work_change).or.present(internal_energy)).and. &
      (.not.present(energy_policy).or.present(energy_correction)).and. &
      (.not.present(energy_correction).or.present(internal_energy))
    if(present(energy_correction)) valid=valid.and.ieee_is_finite(energy_correction).and.ieee_is_finite(held_policy%eta2)
    call agree(valid,comm,ready);if(.not.ready) return
    if(present(energy_correction)) then
      held_controls=[energy_correction,held_policy%eta2]
      held_policy_flags=[merge(1,0,held_policy%allow_small_energy), &
        merge(1,0,held_policy%allow_negative_energy),merge(1,0,held_policy%update_total_from_internal)]
      call MPI_Allreduce(held_controls,held_controls_low,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(held_controls,held_controls_high,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(held_controls_low/=held_controls_high)) return
      call MPI_Allreduce(held_policy_flags,held_policy_min,3,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(held_policy_flags,held_policy_max,3,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(held_policy_min/=held_policy_max)) return
      allocate(correction_work,source=energy_correction,stat=ios)
      call agree(ios==0,comm,ready);if(.not.ready) return
    endif
    ! Registration allocation/overflow acceptance is collective in recursion.
    call MPI_Allreduce(merge(1,0,present(registered_flux)),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,present(registered_flux)),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call MPI_Allreduce(merge(1,0,present(interval_history)),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,present(interval_history)),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,constant_agreed,species=species)
    if(.not.constant_agreed) return
    redist_order=0
    if(present(state_redist_max_order)) redist_order=state_redist_max_order
    call agree(redist_order==0.or.redist_order==2.or.redist_order==4,comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(redist_order,redist_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(redist_order,redist_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.redist_low/=redist_high) return
    call MPI_Allreduce(merge(1,0,present(evaluate_spatial_diffusion)),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,present(evaluate_spatial_diffusion)),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    has_auxiliary=present(auxiliary_source).or.present(halo_source)
    routing=.false.
    if(present(route_auxiliary)) routing=route_auxiliary
    old_auxiliary_only=.false.
    if(present(mol_mode)) old_auxiliary_only=mol_mode.and.chemistry.and.iterations>1
    valid=ieee_is_finite(dt).and..not.(present(auxiliary_source).and.present(halo_source))
    if(valid) valid=dt>0
    valid=valid.and.all([size(boundary_change),size(wall_change),size(chemical_change),size(source_change)]==nv)
    valid=valid.and.(.not.routing.or.has_auxiliary)
    valid=valid.and.(has_auxiliary.eqv.present(auxiliary_change))
    if(present(auxiliary_change)) valid=valid.and.size(auxiliary_change)==naux
    valid=valid.and.(chemistry.or..not.present(evaluate_reaction))
    if(present(internal_energy)) then
      valid=valid.and. &
        (.not.has_auxiliary.or.(held_auxiliary.or.(chemistry.and..not.routing)))
    endif
    if(present(interval_history)) then
      valid=valid.and.chemistry
      if(present(mol_mode)) then
        valid=valid.and.mol_mode
      else
        valid=.false.
      endif
    endif
    call agree(valid,comm,ready)
    if(.not.ready) return
    call MPI_Allreduce(merge(1,0,has_auxiliary)+merge(2,0,present(halo_source)),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,has_auxiliary)+merge(2,0,present(halo_source)),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call MPI_Allreduce(merge(1,0,present(predictor_builder)),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,present(predictor_builder)),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call MPI_Allreduce(merge(1,0,present(post_level)),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,present(post_level)),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call MPI_Allreduce(merge(1,0,routing),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,routing),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call MPI_Allreduce(merge(1,0,present(evaluate_reaction)),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,present(evaluate_reaction)),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    method='none'
    if(present(redistribution)) method=redistribution
    call MPI_Allreduce(merge(1,0,selected_momentum),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,selected_momentum),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call agree((.not.selected_pressure.or.selected_velocity).and. &
      (.not.selected_pressure.or..not.present(frame)),comm,ready)
    if(.not.ready) return
    call agree((.not.selected_velocity.or.selected_momentum).and. &
      (.not.selected_connected.or.selected_momentum),comm,ready)
    if(.not.ready) return
    call agree(.not.selected_momentum.or.(godunov_mol.and.present(internal_energy).and.method=='none'),comm,ready)
    if(.not.ready) return
    call agree(method=='none'.or.method=='flux_redist'.or.method=='state_redist',comm,ready)
    if(.not.ready) return
    if(present(internal_energy)) then
      call agree((method=='none'.or.method=='state_redist').and. &
        (.not.chemistry.or.present(energy_correction)),comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,ready)
      if(.not.ready) return
    endif
    method_code=merge(1,0,method=='flux_redist')+merge(2,0,method=='state_redist')
    call MPI_Allreduce(method_code,method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(method_code,method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    prepared_redist=.false.
    if(present(mol_mode)) prepared_redist=mol_mode.and.method=='state_redist'
    ! Grow9 or grow10 geometry selects the composed Godunov stage; grow6 callers retain
    ! their existing route until application geometry/callback migration.
    prepared_godunov=.false.
    if(method=='state_redist'.and..not.prepared_redist) then
      if(size(extended_geometry)==size(tree%nodes)) then
        prepared_godunov=.true.
        do node=1,size(tree%nodes)
          n=tree%nodes(node)%shape
          if(.not.(all([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]==n+18).or. &
            all([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]==n+20))) &
            prepared_godunov=.false.
        enddo
      endif
    endif
    prepared_redist=prepared_redist.or.prepared_godunov
    call agree(.not.(present(internal_energy).and.method=='state_redist'.and..not.prepared_redist),comm,ready)
    if(.not.ready) return
    godunov_padding=9
    if(prepared_godunov.and.size(extended_geometry)>0) then
      godunov_padding=(extended_geometry(1)%nx-tree%nodes(1)%shape(1))/2
    endif
    call MPI_Allreduce(merge(1,0,prepared_redist),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,prepared_redist),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    redist_target=.5_dp
    if(present(target_volume_fraction)) redist_target=target_volume_fraction
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,comm,ready)
    if(.not.ready) return
    levels=size(ratios)
    if(present(godunov_replay_history)) then
      call initialize_godunov_replay(valid)
      if(.not.valid) return
    else
      valid=.not.present(godunov_parent_registers).and..not.present(godunov_parent_scalar_registers).and. &
        .not.present(godunov_routed_pending)
      call agree(valid,comm,ready);if(.not.ready) return
    endif
    if(replaying_native) then
      call validate_native_window_history(ready)
      if(.not.ready) return
    endif
    candidate=state;ct=temperature
    if(present(internal_energy)) held_candidate=internal_energy
    if(allocated(correction_work)) then
      call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,candidate,nv,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,ct,1,comm,ready)
      if(.not.ready) return
      held_candidate=internal_energy
      if(godunov_replay_depth==0.and..not.replaying_native) then
      call reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,candidate,ct,held_candidate,reset_difference, &
        held_policy,-1,reset_change,comm,ready,frame, &
        internal_energy_change=initial_clock_reset,reset_level=merge(-2,-1,godunov_mol.and.chemistry.and.method=='none'))
      if(.not.ready) return
      correction_work=correction_work+reset_change
      endif
    endif
    if(godunov_replay_depth==0) then
    if(.not.(godunov_mol.and.chemistry.and.method=='none'.and.allocated(correction_work))) then
    if(present(internal_energy)) then
      call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,ready,frame,held_candidate)
    else
      call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,ready,frame)
    endif
    if(.not.ready) return
    endif
    endif
    if(replaying_native) then
      candidate=native_replay_history(1)%clock_initial;ct=native_replay_history(1)%clock_temperature
      held_candidate=native_replay_history(1)%clock_internal_energy
      correction_work=native_replay_history(1)%clock_energy_correction
    endif
    replay_enabled=.false.
    if(present(mol_mode)) replay_enabled=mol_mode.and.chemistry.and. &
      (method=='none'.or.method=='flux_redist'.or.(method=='state_redist'.and.prepared_redist))
    call MPI_Allreduce(merge(1,0,replay_enabled),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,replay_enabled),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    call agree(.not.godunov_capture.or.(.not.godunov_mol.and.present(internal_energy).and. &
      .not.present(interval_history).and.(method=='none'.or.prepared_godunov)),comm,ready)
    if(.not.ready) return
    if(present(godunov_replay_history).or.present(spatial_flux_weight)) then
      valid=.not.godunov_mol.and.present(internal_energy).and.(method=='none'.or.prepared_godunov)
      valid=valid.and..not.present(interval_history).and..not.present(post_level)
      call agree(valid,comm,ready);if(.not.ready) return
    endif
    if(present(godunov_replay_history)) then
      call validate_godunov_family(valid)
      if(.not.valid) return
    endif
    godunov_automatic_repair=.not.godunov_mol.and.chemistry.and.present(internal_energy).and. &
      .not.present(godunov_replay_history).and.(method=='none'.or.prepared_godunov)
    godunov_automatic_repair=godunov_automatic_repair.and..not.present(post_level).and..not.present(interval_history)
    call MPI_Allreduce(merge(1,0,godunov_automatic_repair),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,godunov_automatic_repair),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.method_low/=method_high) return
    godunov_capture=godunov_capture.or.godunov_automatic_repair
    if(godunov_capture) then
      godunov_capacity=1;godunov_level_calls=1
      do c=godunov_replay_depth+2,levels
        if(godunov_level_calls>huge(godunov_capacity)/ratios(c)) return
        godunov_level_calls=godunov_level_calls*ratios(c)
        if(godunov_capacity>huge(godunov_capacity)-godunov_level_calls) return
        godunov_capacity=godunov_capacity+godunov_level_calls
      enddo
      allocate(godunov_records(godunov_capacity),godunov_active(0:levels-1),stat=godunov_status)
      call agree(godunov_status==0,comm,ready);if(.not.ready) return
      godunov_active=0
      if(godunov_replay_depth>0) godunov_active(godunov_replay_depth-1)=godunov_parent
    endif
    collect_history=present(interval_history).or.replay_enabled.or.replaying_native
    replay_enabled=replay_enabled.and..not.replaying_native
    history_count=0
    if(collect_history) then
      history_size=1;history_level_count=1;valid=.true.
      do c=native_replay_depth+2,levels
        if(ratios(c)<=0) then
          valid=.false.;exit
        endif
        if(history_level_count>huge(history_size)/ratios(c)) then
          valid=.false.;exit
        endif
        history_level_count=history_level_count*ratios(c)
        if(history_size>huge(history_size)-history_level_count) then
          valid=.false.;exit
        endif
        history_size=history_size+history_level_count
      enddo
      call agree(valid,comm,ready)
      if(.not.ready) return
      allocate(history(history_size),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
    endif
    if(present(level_steps)) then
      call agree(size(level_steps)==levels,comm,ready)
    else
      call agree(.true.,comm,ready)
    end if
    if(.not.ready) return
    if(godunov_replay_depth>0) then
      associate(anchor=>godunov_replay_history(1))
      call copy_godunov_field_3d(anchor%clock_gas,candidate,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      call copy_godunov_field_3d(anchor%clock_temperature,ct,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      call copy_godunov_field_3d(anchor%clock_internal_energy,held_candidate,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      if(allocated(correction_work)) correction_work=anchor%clock_energy_correction
      end associate
    endif
    funding_controls=[merge(1,0,present(source_funding_context)),0]
    if(present(source_funding_context)) &
      funding_controls(2)=merge(1,0,allocated(source_funding_context%pending))
    call MPI_Allreduce(funding_controls,funding_low,2,MPI_INTEGER,MPI_MIN,comm,status)
    call agree(status==MPI_SUCCESS,comm,ready);if(.not.ready) return
    call MPI_Allreduce(funding_controls,funding_high,2,MPI_INTEGER,MPI_MAX,comm,status)
    call agree(status==MPI_SUCCESS,comm,ready);if(.not.ready) return
    if(any(funding_low/=funding_high)) return
    if(present(source_funding_context)) then
      if(.not.present(mol_mode)) return
      if(.not.mol_mode.or..not.present(internal_energy).or..not.present(energy_correction)) return
      if(method/='none') return
      if(present(godunov_replay_history)) return
      allocate(funding_context_work,funding_local(size(tree%nodes)),stat=funding_status)
      call agree(funding_status==0,comm,ready);if(.not.ready) return
      funding_local=owners==rank
      if(replaying_native) then
        call copy_eb_held_funding_context_3d(tree,geometry,nv,funding_local, &
          native_replay_history(1)%clock_source_funding_context,funding_context_work,valid)
      else if(allocated(source_funding_context%pending)) then
        call copy_eb_held_funding_context_3d(tree,geometry,nv,funding_local,source_funding_context, &
          funding_context_work,valid)
      else
        valid=.not.allocated(source_funding_context%open).and. &
          .not.allocated(source_funding_context%geometry).and..not.allocated(source_funding_context%topology%nodes)
        call agree(valid,comm,ready);if(.not.ready) return
        call initialize_eb_held_funding_context_3d(tree,geometry,nv,funding_local, &
          funding_context_work,valid)
      endif
      call agree(valid,comm,ready);if(.not.ready) return
      valid=valid_eb_held_funding_context_3d(tree,geometry,nv,funding_local,funding_context_work,native_replay_depth)
      call agree(valid,comm,ready);if(.not.ready) return
    endif
    starts=candidate;ends=candidate;ts=ct;te=ct
    if(present(internal_energy)) then
      held_starts=held_candidate;held_ends=held_candidate
    endif
    allocate(clock_start(levels),clock_dt(levels),counts(levels),source0%fields(size(tree%nodes)),stat=ios)
    call agree(ios==0,comm,ready)
    if(.not.ready) return
    counts=0;clock_start=0;clock_dt=dt
    allocate(deferred_record(0:levels-1),deferred_iteration(0:levels-1),deferred_cycle(0:levels-1),stat=ios)
    call agree(ios==0,comm,ready)
    if(.not.ready) return
    deferred_record=0;deferred_iteration=0;deferred_cycle=0
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(source0%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp,stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
    end do
    call agree(valid,comm,ready)
    if(.not.ready) return
    source1=source0
    if(replaying_native) then
      associate(anchor=>native_replay_history(1))
      clock_start=anchor%clock_start;clock_dt=anchor%clock_interval
      starts=anchor%ancestor_start;ends=anchor%ancestor_end
      ts=anchor%ancestor_temperature_start;te=anchor%ancestor_temperature_end
      held_starts=anchor%ancestor_internal_energy_start;held_ends=anchor%ancestor_internal_energy_end
      source0=anchor%ancestor_source_start;source1=anchor%ancestor_source_end
      end associate
    endif
    if(godunov_replay_depth>0) then
      associate(anchor=>godunov_replay_history(1))
      clock_start=anchor%clock_start;clock_dt=anchor%clock_interval
      call copy_godunov_field_3d(anchor%ancestor_start,starts,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      call copy_godunov_field_3d(anchor%ancestor_end,ends,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      call copy_godunov_field_3d(anchor%ancestor_temperature_start,ts,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      call copy_godunov_field_3d(anchor%ancestor_temperature_end,te,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      call copy_godunov_field_3d(anchor%ancestor_internal_energy_start,held_starts,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      call copy_godunov_field_3d(anchor%ancestor_internal_energy_end,held_ends,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      call copy_godunov_field_3d(anchor%source_start,source0,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      call copy_godunov_field_3d(anchor%source_end,source1,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      end associate
    endif

    if(routing) then
      allocate(auxiliary_pending(0:levels-1))
      do c=0,levels-1
        auxiliary_pending(c)=source0
        if(held_auxiliary) then
          status=0
          do node=1,size(tree%nodes)
            if(owners(node)/=rank) cycle
            n=tree%nodes(node)%shape
            deallocate(auxiliary_pending(c)%fields(node)%values)
            allocate(auxiliary_pending(c)%fields(node)%values(naux,n(1),n(2),n(3)),source=0._dp,stat=status)
            if(status/=0) exit
          enddo
          call agree(status==0,comm,ready);if(.not.ready) return
        endif
      end do
    end if
    if(godunov_replay_depth>0.and.routing) then
      do c=0,levels-1
        call copy_godunov_field_3d(godunov_replay_history(1)%clock_routed_pending(c),auxiliary_pending(c),valid)
        call agree(valid,comm,ready);if(.not.ready) return
      enddo
    endif
    if(method/='none'.and..not.prepared_redist) then
      allocate(pending(0:levels-1),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      do c=0,levels-1
        pending(c)=source0
      end do
    end if
    if(prepared_redist) then
      allocate(register_geometry(size(tree%nodes)),stat=ios)
      call agree(ios==0,comm,ready)
      if(.not.ready) return
      call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,ready, &
        padding=merge(godunov_padding,6,prepared_godunov))
      if(.not.ready) return
      valid=.true.
      do node=1,size(tree%nodes)
        n=tree%nodes(node)%shape
        if(prepared_godunov) then
          call extract_eb_geometry_3d(extended_geometry(node),[1,1,1]*(godunov_padding-5), &
            n+godunov_padding+6,register_geometry(node),valid)
          if(.not.valid) exit
        else
          register_geometry(node)=extended_geometry(node)
        endif
      enddo
      call agree(valid,comm,ready)
      if(.not.ready) return
      allocate(interface_registers(0:levels-2))
      if(present(internal_energy)) allocate(held_interface_registers(0:levels-2))
      valid=.true.
      do c=0,levels-2
        allocate(interface_registers(c)%fields(size(tree%nodes)))
        if(present(internal_energy)) allocate(held_interface_registers(c)%fields(size(tree%nodes)))
        do node=1,size(tree%nodes)
          if(owners(node)/=rank) cycle
          call initialize_state_redist_tree_register_3d(tree,node,nv,6,interface_registers(c)%fields(node), &
            ready,boundary_faces,scale=merge(1._dp,.5_dp,prepared_godunov))
          if(.not.ready) then
            valid=.false.;exit
          endif
          if(present(internal_energy)) then
            call initialize_state_redist_tree_register_3d(tree,node,1,6,held_interface_registers(c)%fields(node), &
              ready,boundary_faces,scale=merge(1._dp,.5_dp,prepared_godunov))
            if(.not.ready) then
              valid=.false.;exit
            endif
          endif
        enddo
        if(.not.valid) exit
      enddo
      call agree(valid,comm,ready)
      if(.not.ready) return
    endif
    if(godunov_replay_depth>0.and.prepared_redist) then
      associate(anchor=>godunov_replay_history(1))
      valid=.true.
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        godunov_expected_register=interface_registers(godunov_replay_depth-1)%fields(node)
        call accumulate_state_redist_register_3d(godunov_expected_register,anchor%open_ancestor_registers(node), &
          nv,n+12,valid)
        if(.not.valid) exit
        godunov_expected_register=held_interface_registers(godunov_replay_depth-1)%fields(node)
        call accumulate_state_redist_register_3d(godunov_expected_register,anchor%open_ancestor_scalar_registers(node), &
          1,n+12,valid)
        if(.not.valid) exit
      enddo
      call agree(valid,comm,ready);if(.not.ready) return
      call copy_godunov_registers_3d(anchor%open_ancestor_registers,interface_registers(godunov_replay_depth-1)%fields,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      call copy_godunov_registers_3d(anchor%open_ancestor_scalar_registers, &
        held_interface_registers(godunov_replay_depth-1)%fields,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      end associate
    endif
    if(present(heat_release)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,heat_release,1,comm,valid)
      if(.not.valid) return
      allocate(heat);heat=heat_release
    endif
    if(chemistry.or.present(reaction_source).or.godunov_capture.or.present(godunov_replay_history)) then
      allocate(ir)
      if(present(reaction_source)) then
        ir=reaction_source
      else
        allocate(ir%fields(size(tree%nodes)))
        do node=1,size(tree%nodes)
          if(owners(node)/=rank) cycle
          n=tree%nodes(node)%shape
          allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)),source=0._dp,stat=ios)
          if(ios/=0) then
            valid=.false.;exit
          end if
        end do
      end if
    end if
    call agree(valid,comm,ready)
    if(.not.ready) return
    if(godunov_replay_depth>0) then
      call copy_godunov_field_3d(godunov_replay_history(1)%clock_reaction_source,ir,valid)
      call agree(valid,comm,ready);if(.not.ready) return
      if(present(heat_release)) then
        call copy_godunov_field_3d(godunov_replay_history(1)%clock_heat,heat,valid)
        call agree(valid,comm,ready);if(.not.ready) return
      endif
    endif
    if(replaying_native) then
      ir=native_replay_history(1)%clock_reaction_source
      if(present(heat_release)) heat=native_replay_history(1)%clock_heat
    endif
    auxiliary=0;source_time=0
    wall=0;chemical=0;external=0;theta=godunov_minimum_offset;pressure_work_change=0
    if(replaying_native) then
      call advance_level(native_replay_depth,native_replay_time,dt,1,1,synchronization_flux,total,held_total,valid)
    else
      call advance_level(godunov_replay_depth,godunov_replay_time,dt,1,1,synchronization_flux,total,held_total,valid)
    endif
    if(.not.valid) then
      if(godunov_automatic_repair.and..not.godunov_provider_rejected.and..not.godunov_source_exhausted) &
        call repair_godunov_root()
      if(.not.ok.and.present(godunov_recovery_history)) call publish_godunov_recovery()
      return
    endif
    if(allocated(correction_work)) then
      if(.not.(godunov_mol.and.chemistry.and.method=='none')) then
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,ready,frame,held_candidate)
        if(.not.ready) return
      endif
      call agree(ieee_is_finite(correction_work),comm,ready)
      if(.not.ready) return
    endif
    if(prepared_redist) then
      valid=.true.
      do c=godunov_replay_depth,levels-2
        do node=1,size(tree%nodes)
          if(owners(node)/=rank) cycle
          if(allocated(interface_registers(c)%fields(node)%coarse_delta)) then
            if(any(interface_registers(c)%fields(node)%coarse_delta/=0)) valid=.false.
          endif
          if(allocated(interface_registers(c)%fields(node)%fine_delta)) then
            if(any(interface_registers(c)%fields(node)%fine_delta/=0)) valid=.false.
          endif
          if(present(internal_energy)) then
            if(allocated(held_interface_registers(c)%fields(node)%coarse_delta)) then
              if(any(held_interface_registers(c)%fields(node)%coarse_delta/=0)) valid=.false.
            endif
            if(allocated(held_interface_registers(c)%fields(node)%fine_delta)) then
              if(any(held_interface_registers(c)%fields(node)%fine_delta/=0)) valid=.false.
            endif
          endif
        enddo
      enddo
      call agree(valid,comm,ready)
      if(.not.ready) return
    endif
    boundary=0
    if(rank==owners(1)) then
    n=tree%nodes(1)%shape
    associate(g=>geometry(1),f=>total(1))
    do c=1,nv
      boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
        g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
        g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
        g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
        g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
        g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
    end do
    end associate
    end if
    call MPI_Bcast(boundary,nv,MPI_DOUBLE_PRECISION,owners(1),comm,status)
    call agree(status==MPI_SUCCESS.and.all(ieee_is_finite(boundary)).and.all(ieee_is_finite(wall)).and. &
      all(ieee_is_finite(chemical)).and.all(ieee_is_finite(external)),comm,ready)
    if(.not.ready) return
    if(replaying_native) then
      call agree(native_cursor==size(native_replay_history),comm,ready);if(.not.ready) return
    endif
    if(present(godunov_replay_history)) then
      call agree(godunov_replay_cursor==godunov_replay_count,comm,ready)
      if(.not.ready) return
    endif
    if(godunov_capture) then
      call rebuild_godunov_held_clock_prefixes_3d(godunov_records,valid)
      call agree(valid,comm,ready);if(.not.ready) return
    endif
    if(collect_history.and.present(internal_energy)) then
      call rebuild_mpi_mol_held_clock_prefixes_3d(history,valid)
      call agree(valid,comm,ready);if(.not.ready) return
    endif
    call prepare_godunov_parent_outputs(valid)
    if(.not.valid) return
    if(present(interval_history)) then
      call agree(history_count==size(history),comm,ready)
      if(.not.ready) return
      call move_alloc(history,interval_history)
    endif
    if(godunov_capture) then
      call agree(godunov_cursor==size(godunov_records),comm,ready);if(.not.ready) return
      if(present(godunov_interval_history)) call move_alloc(godunov_records,godunov_interval_history)
    endif
    if(present(godunov_parent_registers)) call move_alloc(godunov_parent_work,godunov_parent_registers)
    if(present(godunov_parent_scalar_registers)) &
      call move_alloc(godunov_parent_scalar_work,godunov_parent_scalar_registers)
    if(present(godunov_routed_pending)) call move_alloc(godunov_pending_work,godunov_routed_pending)
    if(present(energy_correction)) energy_correction=correction_work
    if(present(internal_energy_pressure_work_change)) internal_energy_pressure_work_change=pressure_work_change
    call move_alloc(candidate%fields,state%fields);call move_alloc(ct%fields,temperature%fields)
    boundary_change=boundary;wall_change=wall;chemical_change=chemical;source_change=external
    if(present(auxiliary_change)) auxiliary_change=auxiliary
    minimum_theta=theta
    if(present(heat_release)) heat_release=heat
    if(present(reaction_source)) reaction_source=ir
    if(present(level_steps)) level_steps=counts
    if(present(registered_flux)) call move_alloc(registration,registered_flux)
    if(present(integrated_flux)) call move_alloc(total,integrated_flux)
    if(present(internal_energy)) then
      internal_energy=held_candidate
      call move_alloc(held_total,internal_energy_flux)
    endif
    if(present(source_funding_context)) &
      call move_eb_held_funding_context_3d(funding_context_work,source_funding_context)
    ok=.true.
    end subroutine run_recursive
    subroutine validate_native_window_history(good)
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: checked,checked_held
      integer :: r,lev,c,stage,count_min,count_max,seen(0:levels-1),header(4),lo(4),hi(4)
      real(dp) :: cursor(0:levels-1),parent_start(0:levels-1),parent_dt(0:levels-1)
      real(dp) :: expected_dt,tolerance,expected_start,expected_end,clock(3),clock_lo(3),clock_hi(3)
      real(dp) :: ending,profile(2*levels),profile_lo(2*levels),profile_hi(2*levels)
      integer :: p
      logical :: valid
      valid=godunov_mol.and.chemistry.and.present(internal_energy).and.present(energy_correction).and. &
        method=='none'.and..not.routing.and..not.godunov_capture.and..not.present(godunov_replay_history).and. &
        .not.present(spatial_flux_weight).and..not.present(halo_source).and..not.present(predictor_builder).and. &
        .not.present(post_level)
      call agree(valid,comm,good);if(.not.good) return
      call MPI_Allreduce(size(native_replay_history),count_min,1,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) then
        good=.false.;return
      endif
      call MPI_Allreduce(size(native_replay_history),count_max,1,MPI_INTEGER,MPI_MAX,comm,status)
      good=status==MPI_SUCCESS.and.count_min==count_max.and.count_min>0
      if(.not.good) return
      native_replay_depth=native_replay_history(1)%depth;native_replay_time=native_replay_history(1)%time
      valid=valid_native_mol_clock_profile_3d(native_replay_history(1)%clock_start, &
        native_replay_history(1)%clock_interval,native_replay_depth,native_replay_time,dt,levels)
      call agree(valid,comm,good);if(.not.good) return
      valid=native_replay_history(1)%interval==dt.and. &
        (.not.(native_replay_depth==0).or.native_replay_time==0).and. &
        (native_replay_depth==0.or.present(source_funding_context))
      call agree(valid,comm,good);if(.not.good) return
      ending=native_replay_time+dt
      call agree(ieee_is_finite(ending),comm,good);if(.not.good) return
      call agree(ending>native_replay_time,comm,good);if(.not.good) return
      cursor=native_replay_time;seen=0
      parent_start=native_replay_history(1)%clock_start;parent_dt=native_replay_history(1)%clock_interval
      tolerance=32*epsilon(1._dp)*max(abs(ending),tiny(1._dp))
      do r=1,count_min
        associate(saved=>native_replay_history(r))
        call agree(ieee_is_finite(saved%native_spatial_weight),comm,good);if(.not.good) return
        call agree(saved%native_spatial_weight>=0.and.saved%native_spatial_weight<=1,comm,good)
        if(.not.good) return
        call MPI_Allreduce(saved%native_spatial_weight,native_mol_weight_low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
        if(status/=MPI_SUCCESS) then
          good=.false.;return
        endif
        call MPI_Allreduce(saved%native_spatial_weight,native_mol_weight_high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
        good=status==MPI_SUCCESS.and.native_mol_weight_low==native_mol_weight_high
        if(.not.good) return
        header=[saved%depth,saved%native_iterations,merge(1,0,saved%has_auxiliary_source), &
          merge(1,0,saved%has_source_funding)]
        call MPI_Allreduce(header,lo,4,MPI_INTEGER,MPI_MIN,comm,status)
        if(status/=MPI_SUCCESS) then
          good=.false.;return
        endif
        call MPI_Allreduce(header,hi,4,MPI_INTEGER,MPI_MAX,comm,status)
        good=status==MPI_SUCCESS.and.all(lo==hi)
        if(.not.good) return
        lev=saved%depth
        valid=(saved%separate_momentum.eqv.selected_momentum).and. &
          (saved%consistent_velocity.eqv.selected_velocity).and. &
          (saved%paired_pressure.eqv.selected_pressure).and. &
          (saved%connected_directions.eqv.selected_connected).and. &
          lev>=native_replay_depth.and.lev<levels.and.saved%native_iterations==iterations.and.saved%has_internal_energy.and. &
          .not.saved%prepared_redistribution.and..not.saved%routed_auxiliary.and.saved%redistribution=='none'.and. &
          (saved%has_auxiliary_source.eqv.has_auxiliary).and. &
          (saved%has_source_funding.eqv.present(source_funding_context)).and. &
          (saved%spatial_diffusion_is_source.eqv.external_diffusion_source).and. &
          (saved%auxiliary_has_internal_energy.eqv.held_auxiliary).and. &
          all(saved%native_providers.eqv.[present(evaluate_spatial_diffusion),present(evaluate_diffusion)])
        valid=valid.and.(.not.any(saved%native_providers).or.external_diffusion_source)
        valid=valid.and.all(ieee_is_finite(saved%native_force)).and.all(saved%native_force==force).and. &
          all(ieee_is_finite(saved%native_acceleration)).and.(saved%has_native_frame.eqv.present(frame))
        if(present(gas_acceleration)) then
          valid=valid.and.all(saved%native_acceleration==gas_acceleration)
        else
          valid=valid.and.all(saved%native_acceleration==0)
        endif
        if(present(frame)) valid=valid.and.saved%native_frame%is_valid().and. &
          saved%native_frame%axis==frame%axis.and.saved%native_frame%angular_speed==frame%angular_speed.and. &
          all(saved%native_frame%axis_point==frame%axis_point)
        valid=valid.and.saved%internal_energy_policy%eta2==held_policy%eta2.and. &
          (saved%internal_energy_policy%allow_small_energy.eqv.held_policy%allow_small_energy).and. &
          (saved%internal_energy_policy%allow_negative_energy.eqv.held_policy%allow_negative_energy).and. &
          (saved%internal_energy_policy%update_total_from_internal.eqv.held_policy%update_total_from_internal)
        valid=valid.and.allocated(saved%spatial_flux).and.allocated(saved%internal_energy_flux).and. &
          allocated(saved%fractions_start).and.allocated(saved%fractions_end)
        valid=valid.and.(allocated(saved%auxiliary_stages).eqv.has_auxiliary)
        call agree(valid,comm,good);if(.not.good) return
        valid=valid_native_mol_clock_profile_3d(saved%clock_start,saved%clock_interval, &
          lev,saved%time,saved%interval,levels)
        call agree(valid,comm,good);if(.not.good) return
        profile=[saved%clock_start,saved%clock_interval]
        call MPI_Allreduce(profile,profile_lo,2*levels,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
        if(status/=MPI_SUCCESS) then
          good=.false.;return
        endif
        call MPI_Allreduce(profile,profile_hi,2*levels,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
        good=status==MPI_SUCCESS.and.all(profile_lo==profile_hi)
        if(.not.good) return
        if(r==1) then
          expected_dt=dt
          do c=native_replay_depth-1,0,-1
            call agree(expected_dt<=huge(1._dp)/real(ratios(c+2),dp),comm,good);if(.not.good) return
            expected_dt=expected_dt*real(ratios(c+2),dp)
            call agree(abs(parent_dt(c)-expected_dt)<=tolerance,comm,good);if(.not.good) return
          enddo
        endif
        if(has_auxiliary) then
          valid=lbound(saved%auxiliary_stages,1)==0.and.ubound(saved%auxiliary_stages,1)==iterations
          call agree(valid,comm,good);if(.not.good) return
          do stage=0,iterations
            valid=valid_held_interval_field(saved%auxiliary_stages(stage),tree,owners,lev,rank,naux)
            call agree(valid,comm,good);if(.not.good) return
          enddo
        endif
        clock=[saved%time,saved%interval,saved%clock_energy_correction]
        call agree(all(ieee_is_finite(clock)),comm,good);if(.not.good) return
        call MPI_Allreduce(clock,clock_lo,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
        if(status/=MPI_SUCCESS) then
          good=.false.;return
        endif
        call MPI_Allreduce(clock,clock_hi,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
        good=status==MPI_SUCCESS.and.all(clock_lo==clock_hi)
        if(.not.good) return
        valid=abs(saved%time-cursor(lev))<=tolerance.and.size(saved%fractions_start)==lev.and. &
          size(saved%fractions_end)==lev
        if(r>1) valid=valid.and.lev<=native_replay_history(r-1)%depth+1
        call agree(valid,comm,good);if(.not.good) return
        expected_dt=dt
        do c=native_replay_depth+2,lev+1
          expected_dt=expected_dt/real(ratios(c),dp)
        enddo
        valid=abs(saved%interval-expected_dt)<=tolerance
        do c=0,lev-1
          valid=valid.and.saved%clock_start(c+1)==parent_start(c).and.saved%clock_interval(c+1)==parent_dt(c)
          expected_start=(saved%time-parent_start(c))/parent_dt(c)
          expected_end=(saved%time+saved%interval-parent_start(c))/parent_dt(c)
          valid=valid.and.ieee_is_finite(saved%fractions_start(c+1)).and.ieee_is_finite(saved%fractions_end(c+1)).and. &
            abs(saved%fractions_start(c+1)-expected_start)<=32*epsilon(1._dp).and. &
            abs(saved%fractions_end(c+1)-expected_end)<=32*epsilon(1._dp)
        enddo
        call agree(valid,comm,good);if(.not.good) return
        parent_start(lev)=saved%time;parent_dt(lev)=saved%interval
        cursor(lev)=saved%time+saved%interval;seen(lev)=seen(lev)+1
        call agree(cursor(lev)<=ending+tolerance,comm,good);if(.not.good) return
        valid=valid_held_interval_field(saved%ancestor_start,tree,owners,-1,rank,nv).and. &
          valid_held_interval_field(saved%ancestor_end,tree,owners,-1,rank,nv).and. &
          valid_held_interval_field(saved%ancestor_temperature_start,tree,owners,-1,rank,1).and. &
          valid_held_interval_field(saved%ancestor_temperature_end,tree,owners,-1,rank,1).and. &
          valid_held_interval_field(saved%ancestor_internal_energy_start,tree,owners,-1,rank,1).and. &
          valid_held_interval_field(saved%ancestor_internal_energy_end,tree,owners,-1,rank,1).and. &
          valid_held_interval_field(saved%ancestor_source_start,tree,owners,-1,rank,nv).and. &
          valid_held_interval_field(saved%ancestor_source_end,tree,owners,-1,rank,nv)
        valid=valid.and.valid_mpi_mol_native_provider_cache_3d(tree,owners,geometry,saved,nv,rank)
        call agree(valid,comm,good);if(.not.good) return
        do p=1,size(owners)
          if(owners(p)/=rank.or.tree%nodes(p)%level>=native_replay_depth) cycle
          associate(anchor=>native_replay_history(1))
          valid=valid.and.all(saved%ancestor_start%fields(p)%values==anchor%ancestor_start%fields(p)%values).and. &
            all(saved%ancestor_end%fields(p)%values==anchor%ancestor_end%fields(p)%values).and. &
            all(saved%ancestor_temperature_start%fields(p)%values==anchor%ancestor_temperature_start%fields(p)%values).and. &
            all(saved%ancestor_temperature_end%fields(p)%values==anchor%ancestor_temperature_end%fields(p)%values).and. &
            all(saved%ancestor_internal_energy_start%fields(p)%values== &
              anchor%ancestor_internal_energy_start%fields(p)%values).and. &
            all(saved%ancestor_internal_energy_end%fields(p)%values== &
              anchor%ancestor_internal_energy_end%fields(p)%values).and. &
            all(saved%ancestor_source_start%fields(p)%values==anchor%ancestor_source_start%fields(p)%values).and. &
            all(saved%ancestor_source_end%fields(p)%values==anchor%ancestor_source_end%fields(p)%values)
          end associate
        enddo
        call agree(valid,comm,good);if(.not.good) return
        call reconstruct_mpi_eb_mol_interval_3d(tree,owners,geometry,saved,saved%spatial_flux,nv,checked,comm,good, &
          internal_energy_flux=saved%internal_energy_flux,internal_energy_rate=checked_held)
        if(.not.good) return
        end associate
      enddo
      if(.not.native_prefix) then
        call agree(all(seen(native_replay_depth:)>0).and. &
          all(abs(cursor(native_replay_depth:)-ending)<=tolerance),comm,good);if(.not.good) return
      endif
      associate(anchor=>native_replay_history(1))
      valid=valid_held_interval_field(anchor%clock_initial,tree,owners,-1,rank,nv).and. &
        valid_held_interval_field(anchor%clock_temperature,tree,owners,-1,rank,1).and. &
        valid_held_interval_field(anchor%clock_internal_energy,tree,owners,-1,rank,1).and. &
        valid_held_interval_field(anchor%clock_reaction_source,tree,owners,-1,rank,nv-4)
      if(present(heat_release)) valid=valid.and.valid_held_interval_field(anchor%clock_heat,tree,owners,-1,rank,1)
      call agree(valid,comm,good)
      end associate
    end subroutine validate_native_window_history

    subroutine apply_native_funding_receipt(receipt,good)
      type(amr_tree_state_3d), intent(in) :: receipt
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: gas_work,held_work
      logical :: local_ok
      call copy_godunov_field_3d(candidate,gas_work,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      call copy_godunov_field_3d(held_candidate,held_work,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      call apply_eb_held_auxiliary_increment_3d(gas_work,held_work,receipt,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      call move_alloc(gas_work%fields,candidate%fields)
      call move_alloc(held_work%fields,held_candidate%fields)
    end subroutine apply_native_funding_receipt

    subroutine initialize_godunov_replay(good)
      logical, intent(out) :: good
      logical :: local_ok
      integer :: lev,calls
      local_ok=size(godunov_replay_history)>0.and.lbound(godunov_replay_history,1)==1
      call agree(local_ok,comm,good);if(.not.good) return
      associate(anchor=>godunov_replay_history(1))
      local_ok=anchor%depth>=0.and.anchor%depth<levels.and. &
        anchor%subtree_last>=size(godunov_replay_history).and.ieee_is_finite(anchor%time).and. &
        ieee_is_finite(anchor%interval)
      call agree(local_ok,comm,good);if(.not.good) return
      godunov_replay_depth=anchor%depth;godunov_replay_time=anchor%time
      godunov_replay_count=1;calls=1
      do lev=godunov_replay_depth+2,levels
        local_ok=calls<=huge(calls)/ratios(lev)
      call agree(local_ok,comm,good);if(.not.good) return
        calls=calls*ratios(lev)
        local_ok=godunov_replay_count<=huge(calls)-calls
      call agree(local_ok,comm,good);if(.not.good) return
        godunov_replay_count=godunov_replay_count+calls
      enddo
      local_ok=size(godunov_replay_history)<=godunov_replay_count
      if(.not.godunov_prefix) local_ok=size(godunov_replay_history)==godunov_replay_count
      call agree(local_ok,comm,good);if(.not.good) return
      godunov_replay_base=anchor%subtree_last-godunov_replay_count;godunov_parent=anchor%parent_record
      local_ok=anchor%interval==dt.and.godunov_replay_time>=0
      if(godunov_replay_depth==0) then
        local_ok=local_ok.and.godunov_replay_base==0.and.godunov_parent==0.and.godunov_replay_time==0
      else
        local_ok=local_ok.and.godunov_replay_base>=1.and.godunov_parent>=1.and.godunov_parent<=godunov_replay_base
      endif
      call agree(local_ok,comm,good);if(.not.good) return
      if(godunov_replay_depth>0) then
        local_ok=all(ieee_is_finite([anchor%clock_pressure_work_change,anchor%clock_energy_correction, &
      anchor%clock_minimum_theta]))
      call agree(local_ok,comm,good);if(.not.good) return
      endif
      call MPI_Allreduce(size(godunov_replay_history),method_low,1,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) then
        good=.false.;return
      endif
      call MPI_Allreduce(size(godunov_replay_history),method_high,1,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.method_low/=method_high) then
        good=.false.;return
      endif
      godunov_header=[godunov_replay_depth,godunov_replay_base,godunov_parent]
      godunov_clock_values=[godunov_replay_time,anchor%interval,0._dp,0._dp,1._dp]
      if(godunov_replay_depth>0) &
        godunov_clock_values(3:5)=[anchor%clock_pressure_work_change,anchor%clock_energy_correction,anchor%clock_minimum_theta]
      call MPI_Allreduce(godunov_header,godunov_header_low,3,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) then
        good=.false.;return
      endif
      call MPI_Allreduce(godunov_header,godunov_header_high,3,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(godunov_header_low/=godunov_header_high)) then
        good=.false.;return
      endif
      call MPI_Allreduce(godunov_clock_values,godunov_clock_low,5,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) then
        good=.false.;return
      endif
      call MPI_Allreduce(godunov_clock_values,godunov_clock_high,5,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(godunov_clock_low/=godunov_clock_high)) then
        good=.false.;return
      endif
      if(godunov_replay_depth>0) then
        local_ok=valid_godunov_clock_anchor_3d(anchor,tree,nv,naux,levels,rank,owners)
      call agree(local_ok,comm,good);if(.not.good) return
        local_ok=.true.
        if(prepared_redist) local_ok=present(godunov_parent_registers).and.present(godunov_parent_scalar_registers)
        if(routing) local_ok=local_ok.and.present(godunov_routed_pending)
        if(present(heat_release)) local_ok=local_ok.and.allocated(anchor%clock_heat%fields)
      call agree(local_ok,comm,good);if(.not.good) return
        local_ok=ieee_is_finite(held_policy%eta2)
      call agree(local_ok,comm,good);if(.not.good) return
        local_ok=anchor%policy%eta2==held_policy%eta2.and. &
          (anchor%policy%allow_small_energy.eqv.held_policy%allow_small_energy).and. &
          (anchor%policy%allow_negative_energy.eqv.held_policy%allow_negative_energy).and. &
          (anchor%policy%update_total_from_internal.eqv.held_policy%update_total_from_internal)
      call agree(local_ok,comm,good);if(.not.good) return
        godunov_pressure_offset=anchor%clock_pressure_work_change
        godunov_minimum_offset=anchor%clock_minimum_theta
      endif
      end associate
      good=.true.
    end subroutine initialize_godunov_replay
    subroutine prepare_godunov_parent_outputs(good)
      logical, intent(out) :: good
      logical :: local_ok
      integer :: bucket,p,allocation_status
      good=.true.
      if(godunov_replay_depth>0.and.prepared_redist) then
        call copy_godunov_registers_3d(interface_registers(godunov_replay_depth-1)%fields,godunov_parent_work,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
        call copy_godunov_registers_3d(held_interface_registers(godunov_replay_depth-1)%fields, &
          godunov_parent_scalar_work,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      endif
      if(present(godunov_routed_pending).and.routing) then
        allocate(godunov_pending_work(0:levels-1),stat=allocation_status)
        local_ok=allocation_status==0
      call agree(local_ok,comm,good);if(.not.good) return
        do bucket=0,levels-1
          local_ok=.true.
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            local_ok=local_ok.and.all(ieee_is_finite(auxiliary_pending(bucket)%fields(p)%values))
          enddo
      call agree(local_ok,comm,good);if(.not.good) return
          call copy_godunov_field_3d(auxiliary_pending(bucket),godunov_pending_work(bucket),local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
        enddo
      endif
      good=.true.
    end subroutine prepare_godunov_parent_outputs
    subroutine seal_godunov_group_prefix(first,last,good)
      integer, intent(in) :: first,last
      logical, intent(out) :: good
      integer :: cursor
      logical :: local_ok
      local_ok=first>=1.and.last>=first.and.last==godunov_cursor
      call agree(local_ok,comm,good);if(.not.good) return
      cursor=first-1
      call seal_godunov_prefix_clock(godunov_records(first)%depth,godunov_records(first)%parent_record, &
        godunov_records(first)%time,godunov_records(first)%interval,cursor,last,good)
      if(.not.good) return
      local_ok=last<=cursor
      call agree(local_ok,comm,good);if(.not.good) return
    end subroutine seal_godunov_group_prefix
    recursive subroutine seal_godunov_prefix_clock(depth,parent,time,duration,cursor,last,good)
      integer, intent(in) :: depth,parent,last
      real(dp), intent(in) :: time,duration
      integer, intent(inout) :: cursor
      logical, intent(out) :: good
      integer :: here,step
      real(dp) :: elapsed,child_dt,tolerance
      logical :: local_ok
      local_ok=cursor<size(godunov_records).and.depth>=0.and.depth<levels
      call agree(local_ok,comm,good);if(.not.good) return
      cursor=cursor+1;here=cursor
      if(here<=last) then
        associate(saved=>godunov_records(here))
        local_ok=valid_godunov_source_cache_3d(saved,tree,nv,naux,rank,owners,allow_partial=.true.)
      call agree(local_ok,comm,good);if(.not.good) return
        local_ok=saved%accepted_source_stages==iterations+1.or.here==last
      call agree(local_ok,comm,good);if(.not.good) return
        tolerance=32*epsilon(1._dp)*max(abs(time+duration),tiny(1._dp))
        local_ok=saved%depth==depth.and.saved%parent_record==parent.and. &
          ieee_is_finite(saved%time).and.ieee_is_finite(saved%interval)
      call agree(local_ok,comm,good);if(.not.good) return
        local_ok=abs(saved%time-time)<=tolerance.and.abs(saved%interval-duration)<=tolerance
      call agree(local_ok,comm,good);if(.not.good) return
        end associate
      endif
      if(depth+1<levels) then
        elapsed=0
        do step=1,ratios(depth+2)
          child_dt=duration/real(ratios(depth+2),dp)
          if(step==ratios(depth+2)) child_dt=duration-elapsed
          call seal_godunov_prefix_clock(depth+1,here+godunov_replay_base,time+elapsed,child_dt,cursor,last,good)
          if(.not.good) return
          elapsed=elapsed+child_dt
        enddo
      endif
      if(here<=last) godunov_records(here)%subtree_last=cursor+godunov_replay_base
      good=.true.
    end subroutine seal_godunov_prefix_clock
    subroutine publish_godunov_recovery()
      type(eb_godunov_interval_3d), allocatable :: work(:)
      integer :: r,allocation_status
      logical :: local_ok,ready
      ! Only failure after complete accepted common clocks is recoverable here.
      ! A refused provider or missing stage cannot be turned into a source cache.
      if(godunov_provider_rejected.or.godunov_source_exhausted) return
      local_ok=allocated(godunov_records).and.godunov_cursor>0
      call agree(local_ok,comm,ready);if(.not.ready) return
      local_ok=godunov_cursor<=size(godunov_records)
      call agree(local_ok,comm,ready);if(.not.ready) return
      call seal_godunov_group_prefix(1,godunov_cursor,ready)
      if(.not.ready) return
      allocate(work(godunov_cursor),stat=allocation_status)
      local_ok=allocation_status==0
      call agree(local_ok,comm,ready);if(.not.ready) return
      do r=1,godunov_cursor
        call move_godunov_interval_3d(godunov_records(r),work(r))
      enddo
      call move_alloc(work,godunov_recovery_history)
      godunov_recovery_ready=.true.
    end subroutine publish_godunov_recovery
    subroutine retain_godunov_recovery(first,last,result,good)
      integer, intent(in) :: first,last
      type(godunov_root_trial), intent(inout) :: result
      logical, intent(out) :: good
      integer :: r,known,recovered_last,q,old_count,new_count
      logical :: local_ok
      good=.false.
      if(.not.result%recovery_ready) return
      local_ok=allocated(result%recovery)
      call agree(local_ok,comm,good);if(.not.good) return
      known=last-first+1
      recovered_last=first-1+size(result%recovery)
      local_ok=recovered_last<=size(godunov_records)
      call agree(local_ok,comm,good);if(.not.good) return
      ! Keep original accepted packets and their original incoming clock.
      ! Add only genuinely new accepted tables from the quarantined trial.
      do r=1,size(result%recovery)
        local_ok=valid_godunov_source_cache_3d(result%recovery(r),tree,nv,naux,rank,owners,allow_partial=.true.)
      call agree(local_ok,comm,good);if(.not.good) return
      enddo
      do r=1,min(known,size(result%recovery))
        associate(original=>godunov_records(first+r-1),recovered=>result%recovery(r))
        old_count=original%accepted_source_stages;new_count=recovered%accepted_source_stages
        if(.not.original%source_prefix_only.or.new_count<old_count) cycle
        ! A fresh native trial may stop before reaching a previously accepted
        ! provider. Its absence here does not revoke the original declaration.
        do q=old_count,min(new_count,iterations)
          if(q==new_count.and..not.(recovered%sources(q)%auxiliary_accepted.or. &
            recovered%sources(q)%provider_accepted)) cycle
          call extend_godunov_source_stage_3d(recovered%sources(q),original%sources(q))
        enddo
        if(routing.and.allocated(recovered%routed_auxiliary_new_rate%fields)) &
          call move_alloc(recovered%routed_auxiliary_new_rate%fields,original%routed_auxiliary_new_rate%fields)
        original%accepted_source_stages=new_count
        end associate
      enddo
      do r=known+1,size(result%recovery)
        call move_godunov_interval_3d(result%recovery(r),godunov_records(first+r-1))
      enddo
      godunov_cursor=max(godunov_cursor,recovered_last)
      good=.true.
    end subroutine retain_godunov_recovery
    subroutine repair_godunov_parent(first,last,depth,faces,physical_faces,scalar_faces,good)
      integer, intent(in) :: first,last,depth
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: faces(:),physical_faces(:),scalar_faces(:)
      logical, intent(out) :: good
      type(godunov_root_trial) :: best,trial
      type(eb_integrated_flux_3d), allocatable :: replacement_physical(:),replacement_registration(:)
      real(dp) :: lower,upper,weight,prefix
      real(dp) :: replacement_wall(nv),replacement_chemical(nv),replacement_external(nv),replacement_auxiliary(naux)
      real(dp) :: replacement_pressure_work
      integer :: r,p,ending,lev,latest,attempt,c,replacement_last
      logical :: local_ok,valid
      good=.false.
      if(godunov_provider_rejected.or.godunov_source_exhausted.or. &
        first<1.or.last<first.or.last>godunov_cursor) return
      ! Accepted prefix packets are frozen. Unvisited clocks are evaluated
      ! once by the PCM trial, whose complete accepted archive seeds all later
      ! weights. A refused or incomplete native source clock cannot be replayed.
      call seal_godunov_group_prefix(first,last,good)
      if(.not.good) return
      call build_godunov_group_trial(first,last,0._dp,best,valid)
      if(best%refused) then
        godunov_provider_rejected=.true.
        if(present(provider_refused)) provider_refused=.true.
        good=.false.;return
      endif
      good=.false.
      if(.not.valid) then
        call retain_godunov_recovery(first,last,best,local_ok)
        if(.not.local_ok.and.(last<godunov_records(first)%subtree_last.or. &
          any(godunov_records(first:last)%source_prefix_only)).and. &
          (has_auxiliary.or.external_diffusion_source)) godunov_source_exhausted=.true.
        return
      endif
      lower=0;upper=godunov_weight
      call build_godunov_group_trial(first,last,upper,trial,valid,cache=best%records)
      if(trial%refused) then
        godunov_provider_rejected=.true.
        if(present(provider_refused)) provider_refused=.true.
        return
      endif
      if(valid) then
        best=trial;lower=upper
      else
        do attempt=1,12
          weight=.5_dp*(lower+upper)
          call build_godunov_group_trial(first,last,weight,trial,valid,cache=best%records)
          if(trial%refused) then
            godunov_provider_rejected=.true.
            if(present(provider_refused)) provider_refused=.true.
            return
          endif
          if(valid) then
            best=trial;lower=weight
          else
            upper=weight
          endif
        enddo
      endif
      replacement_last=first-1+size(best%records)
      local_ok=replacement_last<=size(godunov_records)
      call agree(local_ok,comm,good);if(.not.good) return
      replacement_wall=0;replacement_chemical=0;replacement_external=0;replacement_auxiliary=0
      replacement_pressure_work=0
      do r=1,replacement_last
        if(r>=first.and.r<=replacement_last) then
          replacement_wall=replacement_wall+best%records(r-first+1)%wall_change
          replacement_chemical=replacement_chemical+best%records(r-first+1)%chemical_change
          replacement_external=replacement_external+best%records(r-first+1)%source_change
          if(has_auxiliary) replacement_auxiliary=replacement_auxiliary+best%records(r-first+1)%auxiliary_change
          replacement_pressure_work=replacement_pressure_work+best%records(r-first+1)%pressure_work_change
        else
          replacement_wall=replacement_wall+godunov_records(r)%wall_change
          replacement_chemical=replacement_chemical+godunov_records(r)%chemical_change
          replacement_external=replacement_external+godunov_records(r)%source_change
          if(has_auxiliary) replacement_auxiliary=replacement_auxiliary+godunov_records(r)%auxiliary_change
          replacement_pressure_work=replacement_pressure_work+godunov_records(r)%pressure_work_change
        endif
      enddo
      local_ok=all(ieee_is_finite(replacement_wall)).and.all(ieee_is_finite(replacement_chemical)).and. &
        all(ieee_is_finite(replacement_external)).and.all(ieee_is_finite(replacement_auxiliary)).and. &
        ieee_is_finite(replacement_pressure_work)
      call agree(local_ok,comm,good);if(.not.good) return
      ! Prepare physical and outer registration families before adoption.
      call copy_godunov_fluxes_3d(best%faces,replacement_physical,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      if(present(registered_flux)) then
        call copy_godunov_fluxes_3d(registration,replacement_registration,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          replacement_registration(p)%x=0;replacement_registration(p)%y=0;replacement_registration(p)%z=0
        enddo
        do r=1,replacement_last
          if(r>=first) then
            lev=best%records(r-first+1)%depth
          else
            lev=godunov_records(r)%depth
          endif
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            if(tree%nodes(p)%level/=lev) cycle
            if(r>=first.and.r<=replacement_last) then
              replacement_registration(p)%x=replacement_registration(p)%x+best%records(r-first+1)%registered_faces(p)%x
              replacement_registration(p)%y=replacement_registration(p)%y+best%records(r-first+1)%registered_faces(p)%y
              replacement_registration(p)%z=replacement_registration(p)%z+best%records(r-first+1)%registered_faces(p)%z
            else
              replacement_registration(p)%x=replacement_registration(p)%x+godunov_records(r)%registered_faces(p)%x
              replacement_registration(p)%y=replacement_registration(p)%y+godunov_records(r)%registered_faces(p)%y
              replacement_registration(p)%z=replacement_registration(p)%z+godunov_records(r)%registered_faces(p)%z
            endif
          enddo
        enddo
        local_ok=.true.
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          local_ok=local_ok.and.all(ieee_is_finite(replacement_registration(p)%x)).and. &
            all(ieee_is_finite(replacement_registration(p)%y)).and.all(ieee_is_finite(replacement_registration(p)%z))
        enddo
      call agree(local_ok,comm,good);if(.not.good) return
      endif
      call move_alloc(best%gas%fields,candidate%fields)
      call move_alloc(best%temperature%fields,ct%fields)
      call move_alloc(best%stored%fields,held_candidate%fields)
      call move_alloc(best%reaction%fields,ir%fields)
      if(present(heat_release)) call move_alloc(best%heat,heat)
      correction_work=best%prefix
      if(depth>0.and.prepared_redist) then
        call move_alloc(best%parent_registers,interface_registers(depth-1)%fields)
        call move_alloc(best%parent_scalar_registers,held_interface_registers(depth-1)%fields)
      endif
      if(prepared_redist) then
        do c=depth,levels-2
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            if(allocated(interface_registers(c)%fields(p)%coarse_delta)) interface_registers(c)%fields(p)%coarse_delta=0
            if(allocated(interface_registers(c)%fields(p)%fine_delta)) interface_registers(c)%fields(p)%fine_delta=0
            if(allocated(held_interface_registers(c)%fields(p)%coarse_delta)) &
              held_interface_registers(c)%fields(p)%coarse_delta=0
            if(allocated(held_interface_registers(c)%fields(p)%fine_delta)) &
              held_interface_registers(c)%fields(p)%fine_delta=0
          enddo
        enddo
      endif
      if(routing) call move_alloc(best%routed_pending,auxiliary_pending)
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        lev=tree%nodes(p)%level
        if(lev<depth) cycle
        latest=0
        do r=1,size(best%records)
          if(best%records(r)%depth==lev) latest=r
        enddo
        if(latest==0) cycle
        associate(saved=>best%records(latest))
        starts%fields(p)%values=saved%initial_gas%fields(p)%values
        ends%fields(p)%values=saved%final_gas%fields(p)%values
        ts%fields(p)%values=saved%initial_temperature%fields(p)%values
        te%fields(p)%values=saved%final_temperature%fields(p)%values
        held_starts%fields(p)%values=saved%initial_internal_energy%fields(p)%values
        held_ends%fields(p)%values=saved%final_internal_energy%fields(p)%values
        source0%fields(p)%values=saved%old_predictor%fields(p)%values
        source1%fields(p)%values=saved%new_predictor%fields(p)%values
        end associate
      enddo
      do r=1,size(best%records)
        lev=best%records(r)%depth
        clock_start(lev+1)=best%records(r)%time;clock_dt(lev+1)=best%records(r)%interval
      enddo
      prefix=godunov_records(first)%clock_minimum_theta
      do r=first,replacement_last
        call move_godunov_interval_3d(best%records(r-first+1),godunov_records(r))
      enddo
      godunov_cursor=replacement_last
      counts=0
      do r=1,godunov_cursor
        lev=godunov_records(r)%depth
        counts(lev+1)=counts(lev+1)+1
      enddo
      wall=replacement_wall;chemical=replacement_chemical;external=replacement_external
      auxiliary=replacement_auxiliary;pressure_work_change=0
      do r=1,replacement_last
        godunov_records(r)%clock_pressure_work_change=pressure_work_change
        pressure_work_change=pressure_work_change+godunov_records(r)%pressure_work_change
      enddo
      do r=first,replacement_last
        ending=godunov_records(r)%subtree_last
        godunov_records(r)%closed_pressure_work_change=godunov_records(ending)%clock_pressure_work_change+ &
          godunov_records(ending)%pressure_work_change
      enddo
      theta=min(prefix,best%limiter,lower)
      if(present(registered_flux)) call move_alloc(replacement_registration,registration)
      call move_alloc(best%faces,faces);call move_alloc(replacement_physical,physical_faces)
      call move_alloc(best%scalar_faces,scalar_faces)
      good=.true.
    end subroutine repair_godunov_parent
    subroutine repair_godunov_root()
      type(godunov_root_trial) :: best,trial
      real(dp) :: lower,upper,blend
      integer :: r,last,next,attempt
      logical :: valid,local_ok
      ! Only complete native selected clocks may be replayed. A local source,
      ! reaction or allocation refusal has no accepted packet for that clock.
      ! Never fill a missing clock or invoke its external providers again.
      local_ok=godunov_cursor==size(godunov_records)
      call agree(local_ok,comm,valid);if(.not.valid) return
      do r=1,size(godunov_records)
        local_ok=valid_godunov_source_cache_3d(godunov_records(r),tree,nv,naux,rank,owners)
      call agree(local_ok,comm,valid);if(.not.valid) return
        last=size(godunov_records)
        do next=r+1,size(godunov_records)
          if(godunov_records(next)%depth<=godunov_records(r)%depth) then
            last=next-1;exit
          endif
        enddo
        godunov_records(r)%subtree_last=last
      enddo
      call build_godunov_group_trial(1,godunov_cursor,0._dp,best,valid)
      if(best%refused) then
        godunov_provider_rejected=.true.
        if(present(provider_refused)) provider_refused=.true.
        return
      endif
      if(.not.valid) return
      lower=0;upper=godunov_weight
      call build_godunov_group_trial(1,godunov_cursor,upper,trial,valid)
      if(trial%refused) then
        godunov_provider_rejected=.true.
        if(present(provider_refused)) provider_refused=.true.
        return
      endif
      if(valid) then
        best=trial;lower=upper
      else
        do attempt=1,12
          blend=.5_dp*(lower+upper)
          call build_godunov_group_trial(1,godunov_cursor,blend,trial,valid)
          if(trial%refused) then
            godunov_provider_rejected=.true.
            if(present(provider_refused)) provider_refused=.true.
            return
          endif
          if(valid) then
            best=trial;lower=blend
          else
            upper=blend
          endif
        enddo
      endif
      ! A successful complete recursive native driver has closed both gas and
      ! held register families and passed all unchanged reaction/EOS gates.
      ! Replace every state, face, clock and numerical/pressure ledger together.
      call move_alloc(best%gas%fields,state%fields)
      call move_alloc(best%temperature%fields,temperature%fields)
      call move_alloc(best%stored%fields,internal_energy%fields)
      call move_alloc(best%scalar_faces,internal_energy_flux)
      if(present(integrated_flux)) call move_alloc(best%faces,integrated_flux)
      if(present(registered_flux)) call move_alloc(best%registered_faces,registered_flux)
      if(present(godunov_interval_history)) call move_alloc(best%records,godunov_interval_history)
      if(present(reaction_source)) reaction_source=best%reaction
      if(present(heat_release)) heat_release=best%heat
      if(present(level_steps)) level_steps=best%counts
      energy_correction=best%prefix
      if(present(internal_energy_pressure_work_change)) internal_energy_pressure_work_change=best%pressure_work
      boundary_change=best%boundary;wall_change=best%wall;chemical_change=best%chemical;source_change=best%external
      if(present(auxiliary_change)) auxiliary_change=best%auxiliary
      minimum_theta=min(best%limiter,lower)
      ok=.true.
    end subroutine repair_godunov_root
    subroutine build_godunov_group_trial(first,last,weight,result,good,cache)
      integer, intent(in) :: first,last
      real(dp), intent(in) :: weight
      type(godunov_root_trial), intent(out) :: result
      logical, intent(out) :: good
      type(eb_godunov_interval_3d), optional, intent(in) :: cache(:)
      type(eb_godunov_interval_3d), allocatable :: frozen(:)
      logical :: prefix_replay
      integer :: allocation_status
      logical :: local_ok
      good=.false.
      if(present(cache)) then
        allocate(frozen,source=cache,stat=allocation_status)
      else
        allocate(frozen,source=godunov_records(first:last),stat=allocation_status)
      endif
      local_ok=allocation_status==0
      call agree(local_ok,comm,good);if(.not.good) return
      prefix_replay=size(frozen)<frozen(1)%subtree_last-first+1.or.any(frozen%source_prefix_only)
      result%gas=state;result%temperature=temperature
      result%stored=internal_energy;result%prefix=energy_correction
      result%reaction=godunov_records(first)%clock_reaction_source
      if(present(reaction_source)) result%reaction=reaction_source
      allocate(result%boundary(nv),result%wall(nv),result%chemical(nv),result%external(nv), &
        result%counts(levels),stat=allocation_status)
      local_ok=allocation_status==0
      call agree(local_ok,comm,good);if(.not.good) return
      if(has_auxiliary) then
        allocate(result%auxiliary(naux),stat=allocation_status)
        local_ok=allocation_status==0
      call agree(local_ok,comm,good);if(.not.good) return
      endif
      if(present(heat_release)) then
        allocate(result%heat,source=heat_release,stat=allocation_status)
        local_ok=allocation_status==0
      call agree(local_ok,comm,good);if(.not.good) return
      endif
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,geometry,result%gas, &
        result%temperature,godunov_records(first)%interval,iterations,chemistry,transport_flags,force,rtol,atol, &
        extended_geometry,boundary_faces, &
        result%boundary,result%wall,result%chemical,result%external,result%limiter,comm,good, &
        reaction_source=result%reaction,difmag=difmag,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,wall_temperature=wall_temperature, &
        level_steps=result%counts,redistribution=method,target_volume_fraction=redist_target, &
        auxiliary_source=auxiliary_source,auxiliary_change=result%auxiliary,route_auxiliary=route_auxiliary, &
        evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip,frame=frame, &
        state_redist_max_order=redist_order,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width,constant_coefficients=constant_coefficients, &
        heat_release=result%heat,registered_flux=result%registered_faces,integrated_flux=result%faces, &
        internal_energy=result%stored,internal_energy_flux=result%scalar_faces,energy_policy=energy_policy, &
        energy_correction=result%prefix,internal_energy_pressure_work_change=result%pressure_work, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        godunov_interval_history=result%records,godunov_replay_history=frozen,spatial_flux_weight=weight, &
        godunov_replay_prefix=prefix_replay, &
        godunov_recovery_history=result%recovery,godunov_recovery_ready=result%recovery_ready, &
        godunov_parent_registers=result%parent_registers,godunov_parent_scalar_registers=result%parent_scalar_registers, &
        godunov_routed_pending=result%routed_pending,provider_refused=result%refused, &
        gas_acceleration=gas_acceleration,halo_source=halo_source,predictor_builder=predictor_builder, &
        evaluate_diffusion=evaluate_diffusion, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
    end subroutine build_godunov_group_trial
    subroutine validate_godunov_family(good)
      logical, intent(out) :: good
      integer :: cursor
      logical :: local_ok
      cursor=0
      local_ok=lbound(godunov_replay_history,1)==1.and.size(godunov_replay_history)>0
      call agree(local_ok,comm,good);if(.not.good) return
      call validate_godunov_clock(godunov_replay_depth,godunov_parent,godunov_replay_time,dt,cursor,good)
      if(.not.good) return
      local_ok=cursor==godunov_replay_count
      call agree(local_ok,comm,good);if(.not.good) return
    end subroutine validate_godunov_family
    recursive subroutine validate_godunov_clock(depth,parent,time,duration,cursor,good)
      integer, intent(in) :: depth,parent
      real(dp), intent(in) :: time,duration
      integer, intent(inout) :: cursor
      logical, intent(out) :: good
      integer :: here,step,stage,header(4),lower(4),upper(4)
      real(dp) :: elapsed,child_dt,tolerance
      logical :: local_ok
      local_ok=cursor<godunov_replay_count
      call agree(local_ok,comm,good);if(.not.good) return
      cursor=cursor+1;here=cursor
      if(here<=size(godunov_replay_history)) then
      associate(saved=>godunov_replay_history(here))
      local_ok=ieee_is_finite(saved%time).and.ieee_is_finite(saved%interval)
      call agree(local_ok,comm,good);if(.not.good) return
      tolerance=32*epsilon(1._dp)*max(abs(time+duration),tiny(1._dp))
      local_ok=saved%depth==depth.and.saved%parent_record==parent.and. &
        abs(saved%time-time)<=tolerance.and.abs(saved%interval-duration)<=tolerance.and. &
        saved%iterations==iterations.and.(saved%chemistry.eqv.chemistry).and. &
        all(saved%transport_flags.eqv.transport_flags).and.all(ieee_is_finite(saved%force))
      call agree(local_ok,comm,good);if(.not.good) return
      local_ok=all(saved%force==force).and.(saved%prepared_redistribution.eqv.prepared_redist).and. &
        (saved%has_auxiliary.eqv.has_auxiliary).and.(saved%has_spatial_provider.eqv.present(evaluate_spatial_diffusion)).and. &
        (saved%spatial_diffusion_is_source.eqv.external_diffusion_source).and. &
        (saved%auxiliary_has_internal_energy.eqv.held_auxiliary).and.(saved%routed_auxiliary.eqv.routing)
      call agree(local_ok,comm,good);if(.not.good) return
      local_ok=(saved%has_valid_diffusion_provider.eqv.present(evaluate_diffusion))
      call agree(local_ok,comm,good);if(.not.good) return
      local_ok=saved%accepted_source_stages==iterations+1.or.here==size(godunov_replay_history)
      call agree(local_ok,comm,good);if(.not.good) return
      local_ok=valid_godunov_source_cache_3d(saved,tree,nv,naux,rank,owners,allow_partial=godunov_prefix)
      call agree(local_ok,comm,good);if(.not.good) return
      header=[saved%accepted_source_stages,merge(1,0,saved%source_prefix_only), &
        merge(1,0,saved%has_source_funding),merge(1,0,allocated(saved%clock_source_funding_context%pending))]
      call MPI_Allreduce(header,lower,4,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) then
        good=.false.;return
      endif
      call MPI_Allreduce(header,upper,4,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(lower/=upper)) then
        good=.false.;return
      endif
      if(allocated(saved%clock_source_funding_context%pending)) then
        local_ok=valid_eb_held_funding_context_3d(tree,geometry,nv,owners==rank, &
          saved%clock_source_funding_context,depth)
        call agree(local_ok,comm,good);if(.not.good) return
      endif
      do stage=0,iterations
        header=[merge(1,0,saved%sources(stage)%auxiliary_accepted), &
          merge(1,0,saved%sources(stage)%provider_accepted),0,0]
        call MPI_Allreduce(header,lower,4,MPI_INTEGER,MPI_MIN,comm,status)
        if(status/=MPI_SUCCESS) then
          good=.false.;return
        endif
        call MPI_Allreduce(header,upper,4,MPI_INTEGER,MPI_MAX,comm,status)
        if(status/=MPI_SUCCESS.or.any(lower/=upper)) then
          good=.false.;return
        endif
      enddo
      end associate
      endif
      if(depth+1<levels) then
        elapsed=0
        do step=1,ratios(depth+2)
          child_dt=duration/real(ratios(depth+2),dp)
          if(step==ratios(depth+2)) child_dt=duration-elapsed
          call validate_godunov_clock(depth+1,here+godunov_replay_base,time+elapsed,child_dt,cursor,good)
          if(.not.good) return
          elapsed=elapsed+child_dt
        enddo
      endif
      if(here<=size(godunov_replay_history)) then
      local_ok=godunov_replay_history(here)%subtree_last==cursor+godunov_replay_base
      call agree(local_ok,comm,good);if(.not.good) return
      endif
      good=.true.
    end subroutine validate_godunov_clock
    subroutine capture_godunov_clock(level,saved,good)
      integer, intent(in) :: level
      type(eb_godunov_interval_3d), intent(inout) :: saved
      logical, intent(out) :: good
      logical :: local_ok
      integer :: bucket
      call copy_godunov_field_3d(candidate,saved%clock_gas,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      call copy_godunov_field_3d(ct,saved%clock_temperature,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      call copy_godunov_field_3d(ir,saved%clock_reaction_source,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      call copy_godunov_field_3d(held_candidate,saved%clock_internal_energy,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      if(present(heat_release)) then
      call copy_godunov_field_3d(heat,saved%clock_heat,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      endif
      saved%clock_energy_correction=0
      if(allocated(correction_work)) saved%clock_energy_correction=correction_work
      saved%clock_pressure_work_change=godunov_pressure_offset+pressure_work_change
      saved%clock_minimum_theta=theta
      saved%routed_auxiliary=routing
      if(routing) then
        allocate(saved%clock_routed_pending(0:levels-1),stat=godunov_status)
        local_ok=godunov_status==0
      call agree(local_ok,comm,good);if(.not.good) return
        do bucket=0,levels-1
          call copy_godunov_field_3d(auxiliary_pending(bucket),saved%clock_routed_pending(bucket),local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
        enddo
      endif
      allocate(saved%clock_start(levels),saved%clock_interval(levels),stat=godunov_status)
      local_ok=godunov_status==0
      call agree(local_ok,comm,good);if(.not.good) return
      saved%clock_start=clock_start;saved%clock_interval=clock_dt
      if(prepared_redist.and.level>0) then
      call copy_godunov_registers_3d(interface_registers(level-1)%fields,saved%open_ancestor_registers,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      call copy_godunov_registers_3d(held_interface_registers(level-1)%fields,saved%open_ancestor_scalar_registers,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      endif
      good=.true.
    end subroutine capture_godunov_clock
    subroutine close_godunov_clock(saved,faces,scalar_faces,good)
      type(eb_godunov_interval_3d), intent(inout) :: saved
      type(eb_integrated_flux_3d), allocatable, intent(in) :: faces(:),scalar_faces(:)
      logical, intent(out) :: good
      logical :: local_ok
      call copy_godunov_field_3d(candidate,saved%closed_gas,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      call copy_godunov_field_3d(ct,saved%closed_temperature,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      call copy_godunov_field_3d(held_candidate,saved%closed_internal_energy,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      if(allocated(ir)) then
      call copy_godunov_field_3d(ir,saved%closed_reaction_source,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      endif
      call copy_godunov_fluxes_3d(faces,saved%closed_faces,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      call copy_godunov_fluxes_3d(scalar_faces,saved%closed_scalar_faces,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
      saved%subtree_last=godunov_cursor+godunov_replay_base;saved%closed_energy_correction=0
      if(allocated(correction_work)) saved%closed_energy_correction=correction_work
      saved%closed_pressure_work_change=godunov_pressure_offset+pressure_work_change
      good=.true.
    end subroutine close_godunov_clock
    subroutine predictor_adapter(rate,interval,iteration,depth,halo,comm,ok)
      type(amr_tree_state_3d), intent(in) :: rate
      real(dp), intent(in) :: interval
      integer, intent(in) :: iteration,depth
      type(amr_tree_state_3d), intent(inout) :: halo
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
      real(dp) :: first(levels-1),last(levels-1)
      first=0;last=0
      first(:depth)=max(0._dp,min(1._dp,(source_time-clock_start(:depth))/clock_dt(:depth)))
      last(:depth)=max(first(:depth),min(1._dp,(source_time+interval-clock_start(:depth))/clock_dt(:depth)))
      call predictor_builder(rate,interval,iteration,depth,source_time,halo,comm,ok, &
        source0,source1,.5_dp*(first+last))
    end subroutine predictor_adapter
    subroutine source_adapter(state,temperature,interval,stage,depth,rate,comm,ok)
      type(amr_tree_state_3d), intent(in) :: state,temperature
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
      type(amr_tree_state_3d) :: off,packed_starts,packed_ends
      integer :: p,x,y,z,d(3)
      logical :: valid,agreed
      real(dp) :: fraction(levels-1),sample_time
      if(present(halo_source)) then
        fraction=0;sample_time=source_time
        if(stage>0) sample_time=sample_time+interval
        fraction(:depth)=max(0._dp,min(1._dp,(sample_time-clock_start(:depth))/clock_dt(:depth)))
        if(held_auxiliary) then
          call pack_eb_held_auxiliary_state_3d(starts,held_starts,packed_starts,valid)
          call agree(valid,comm,agreed);if(.not.agreed) then
            ok=.false.;return
          endif
          call pack_eb_held_auxiliary_state_3d(ends,held_ends,packed_ends,valid)
          call agree(valid,comm,agreed);if(.not.agreed) then
            ok=.false.;return
          endif
          call halo_source(state,temperature,interval,stage,depth,source_time,rate,comm,ok, &
            packed_starts,packed_ends,ts,te,fraction)
        else
        call halo_source(state,temperature,interval,stage,depth,source_time,rate,comm,ok, &
          starts,ends,ts,te,fraction)
        endif
      else
        call auxiliary_source(state,temperature,interval,stage,depth,source_time,rate,comm,ok)
      end if
      if(.not.routing) return
      call agree(ok,comm,agreed)
      ok=.false.
      if(.not.agreed) return
      valid=allocated(rate%fields)
      if(valid) valid=size(rate%fields)==size(owners).and.lbound(rate%fields,1)==1
      if(valid) then
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          if(.not.allocated(rate%fields(p)%values)) then
            valid=.false.;exit
          endif
          if(any(lbound(rate%fields(p)%values)/=1)) then
            valid=.false.;exit
          endif
        enddo
      endif
      call agree(valid,comm,agreed);if(.not.agreed) return
      call validate_mpi_amr_tree_local_3d(tree,owners,rate,naux,comm,agreed)
      if(.not.agreed) return
      off=rate;valid=.true.
      do p=1,size(tree%nodes)
        if(allocated(off%fields(p)%values)) off%fields(p)%values=0
      end do
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        if(tree%nodes(p)%level==depth) cycle
        d=tree%nodes(p)%shape
        off%fields(p)%values=0
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0.or.covered(tree,p,[x,y,z])) cycle
          off%fields(p)%values(:,x,y,z)=rate%fields(p)%values(:,x,y,z)
        end do
        end do
        end do
        valid=valid.and.all(ieee_is_finite(off%fields(p)%values))
      end do
      call agree(valid,comm,agreed)
      if(.not.agreed) return
      if(stage==0) then
        off_old=off
      else
        off_new=off
      end if
      ok=.true.

    end subroutine source_adapter
    subroutine allocate_replay_flux(faces,good,components)
      type(eb_integrated_flux_3d), allocatable, intent(out) :: faces(:)
      logical, intent(out) :: good
      integer, optional, intent(in) :: components
      integer :: p,n(3),allocation_status,ncomp
      logical :: valid
      ncomp=nv
      if(present(components)) ncomp=components
      allocate(faces(size(tree%nodes)),stat=allocation_status)
      call agree(allocation_status==0,comm,good)
      if(.not.good) return
      valid=.true.
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        n=tree%nodes(p)%shape
        allocate(faces(p)%x(ncomp,0:n(1),n(2),n(3)),faces(p)%y(ncomp,n(1),0:n(2),n(3)), &
          faces(p)%z(ncomp,n(1),n(2),0:n(3)),source=0._dp,stat=allocation_status)
        if(allocation_status/=0) then
          valid=.false.;exit
        endif
      enddo
      call agree(valid,comm,good)
    end subroutine allocate_replay_flux

    subroutine replay_parent_interval(first,last,depth,accepted_flux,good,accepted_held_flux)
      integer, intent(in) :: first,last,depth
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: accepted_flux(:)
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: accepted_held_flux(:)
      logical, intent(out) :: good
      type(mol_replay_result) :: best,trial
      type(eb_integrated_flux_3d), allocatable :: trial_registration(:)
      real(dp) :: lower,upper,blend
      real(dp) :: replacement_chemical(nv),replacement_wall(nv),replacement_source(nv),replacement_auxiliary(naux)
      integer :: attempt,r,p,lev
      logical :: valid
      good=.false.
      prepared_provider_rejected=.false.
      call build_mol_replay_trial(first,last,depth,0._dp,best,valid)
      call agree(.not.prepared_provider_rejected,comm,good)
      prepared_provider_rejected=.not.good
      if(prepared_provider_rejected) return
      good=.false.
      if(.not.valid) return
      lower=0;upper=1
      call build_mol_replay_trial(first,last,depth,1._dp,trial,valid)
      call agree(.not.prepared_provider_rejected,comm,good)
      prepared_provider_rejected=.not.good
      if(prepared_provider_rejected) return
      good=.false.
      if(valid) then
        best=trial;lower=1
      else
        ! All conserved flux components share one interval blend. Each trial
        ! still passes the unchanged reaction/EOS gates on every owner rank.
        ! These twelve deterministic refinements select an admissible blend;
        ! they do not change any state/reaction acceptance tolerance.
        do attempt=1,12
          blend=.5_dp*(lower+upper)
          call build_mol_replay_trial(first,last,depth,blend,trial,valid)
          call agree(.not.prepared_provider_rejected,comm,good)
          prepared_provider_rejected=.not.good
          if(prepared_provider_rejected) return
          good=.false.
          if(valid) then
            best=trial;lower=blend
          else
            upper=blend
          endif
        enddo
      endif
      if(present(registered_flux)) then
        call allocate_replay_flux(trial_registration,valid)
        if(.not.valid) return
        do r=1,history_count
          lev=history(r)%depth
          do p=1,size(tree%nodes)
            if(owners(p)/=rank.or.tree%nodes(p)%level/=lev) cycle
            if(r>=first.and.r<=last) then
              trial_registration(p)%x=trial_registration(p)%x+best%records(r-first+1)%spatial_flux(p)%x
              trial_registration(p)%y=trial_registration(p)%y+best%records(r-first+1)%spatial_flux(p)%y
              trial_registration(p)%z=trial_registration(p)%z+best%records(r-first+1)%spatial_flux(p)%z
            else
              trial_registration(p)%x=trial_registration(p)%x+history(r)%spatial_flux(p)%x
              trial_registration(p)%y=trial_registration(p)%y+history(r)%spatial_flux(p)%y
              trial_registration(p)%z=trial_registration(p)%z+history(r)%spatial_flux(p)%z
            endif
          enddo
        enddo
      endif
      if(allocated(best%funding)) then
        ! Recursive replay returns this window's inventories. Earlier
        ! accepted siblings remain physical and must be counted once.
        replacement_chemical=0;replacement_wall=0;replacement_source=0;replacement_auxiliary=0
        do r=1,history_count
          if(r>=first.and.r<=last) then
            replacement_chemical=replacement_chemical+best%records(r-first+1)%chemical_change
            replacement_wall=replacement_wall+best%records(r-first+1)%wall_change
            replacement_source=replacement_source+best%records(r-first+1)%source_change
            replacement_auxiliary=replacement_auxiliary+best%records(r-first+1)%auxiliary_change
          else
            replacement_chemical=replacement_chemical+history(r)%chemical_change
            replacement_wall=replacement_wall+history(r)%wall_change
            replacement_source=replacement_source+history(r)%source_change
            replacement_auxiliary=replacement_auxiliary+history(r)%auxiliary_change
          endif
        enddo
        call agree(all(ieee_is_finite(replacement_chemical)).and.all(ieee_is_finite(replacement_wall)).and. &
          all(ieee_is_finite(replacement_source)).and.all(ieee_is_finite(replacement_auxiliary)),comm,valid)
        if(.not.valid) return
        best%chemical=replacement_chemical
      endif
      call move_alloc(best%state%fields,candidate%fields)
      call move_alloc(best%temperature%fields,ct%fields)
      call move_alloc(best%reaction_source%fields,ir%fields)
      if(present(internal_energy)) then
        call move_alloc(best%internal_energy%fields,held_candidate%fields)
        correction_work=best%energy_correction
      endif
      if(method=='flux_redist') call move_alloc(best%redistribution_pending,pending)
      if(routing) then
        do lev=depth+1,levels-1
          do p=1,size(tree%nodes)
            if(owners(p)==rank) auxiliary_pending(lev)%fields(p)%values=0
          enddo
        enddo
      endif
      if(present(heat_release)) call move_alloc(best%heat,heat)
      chemical=best%chemical
      if(allocated(best%funding)) then
        call move_eb_held_funding_context_3d(best%funding,funding_context_work)
        wall=replacement_wall;external=replacement_source;auxiliary=replacement_auxiliary
      endif
      if(prepared_redist) then
        wall=best%wall
        call move_alloc(best%registers,interface_registers)
        if(present(internal_energy)) then
          call move_alloc(best%held_registers,held_interface_registers)
          call move_alloc(best%held_starts%fields,held_starts%fields)
          call move_alloc(best%held_ends%fields,held_ends%fields)
        endif
        call move_alloc(best%starts%fields,starts%fields)
        call move_alloc(best%ends%fields,ends%fields)
        call move_alloc(best%ts%fields,ts%fields)
        call move_alloc(best%te%fields,te%fields)
        call move_alloc(best%source0%fields,source0%fields)
        call move_alloc(best%source1%fields,source1%fields)
      endif
      if(present(registered_flux)) call move_alloc(trial_registration,registration)
      do r=first,last
        best%records(r-first+1)%replayed=.true.
        call move_mol_interval(best%records(r-first+1),history(r))
      enddo
      if(.not.prepared_redist) then
        do p=1,size(tree%nodes)
          if(owners(p)/=rank.or.tree%nodes(p)%level<depth) cycle
          if(allocated(best%funding)) then
            do r=first,last
              if(tree%nodes(p)%level/=history(r)%depth) cycle
              starts%fields(p)=history(r)%initial%fields(p);ends%fields(p)=history(r)%local_state%fields(p)
              ts%fields(p)=history(r)%initial_temperature%fields(p);te%fields(p)=history(r)%local_temperature%fields(p)
              held_starts%fields(p)=history(r)%initial_internal_energy%fields(p)
              held_ends%fields(p)=history(r)%local_internal_energy%fields(p)
              source0%fields(p)=history(r)%old_predictor%fields(p)
              source1%fields(p)=history(r)%new_predictor%fields(p)
            enddo
          else
            ends%fields(p)=candidate%fields(p);te%fields(p)=ct%fields(p)
            if(present(internal_energy)) held_ends%fields(p)=held_candidate%fields(p)
          endif
        enddo
      endif
      call move_alloc(best%flux,accepted_flux)
      call move_alloc(best%held_flux,accepted_held_flux)
      theta=min(theta,lower)
      good=.true.
    end subroutine replay_parent_interval

    subroutine build_prepared_replay_trial(first,last,depth,blend,result,good)
      integer, intent(in) :: first,last,depth
      real(dp), intent(in) :: blend
      type(mol_replay_result), intent(out) :: result
      logical, intent(out) :: good
      type(prepared_replay_work), target :: work
      type(eb_integrated_flux_3d), allocatable :: faces(:),scalar_faces(:)
      real(dp) :: chemical_trial(nv),wall_trial(nv)
      integer :: p,c,r,allocation_status
      logical :: valid,agreed
      good=.false.
      valid=allocated(history(first)%clock_initial%fields).and. &
        allocated(history(first)%clock_temperature%fields).and. &
        allocated(history(first)%clock_reaction_source%fields)
      if(depth>0) valid=valid.and.allocated(history(first)%open_ancestor_registers)
      if(present(internal_energy)) then
        valid=valid.and.history(first)%has_internal_energy.and. &
          allocated(history(first)%clock_internal_energy%fields).and. &
          ieee_is_finite(history(first)%clock_energy_correction)
        if(depth>0) valid=valid.and.allocated(history(first)%open_ancestor_held_registers)
      endif
      if(present(heat_release)) valid=valid.and.allocated(history(first)%clock_heat%fields)
      call agree(valid,comm,agreed)
      valid=agreed
      if(.not.valid) return
      work%gas=history(first)%clock_initial;work%temp=history(first)%clock_temperature
      work%reaction=history(first)%clock_reaction_source
      work%starts=starts;work%ends=ends;work%ts=ts;work%te=te
      work%source0=source0;work%source1=source1
      work%clock_start=clock_start;work%clock_dt=clock_dt
      work%registers=interface_registers
      if(present(internal_energy)) then
        work%held=history(first)%clock_internal_energy
        work%energy_correction=history(first)%clock_energy_correction
        work%held_starts=held_starts;work%held_ends=held_ends
        work%held_registers=held_interface_registers
        if(depth>0) work%held_registers(depth-1)%fields=history(first)%open_ancestor_held_registers
      endif
      ! Completed sibling clocks stay consumed. The first record supplies the
      ! exact still-open parent ledger before this group's fine contribution.
      if(depth>0) work%registers(depth-1)%fields=history(first)%open_ancestor_registers
      do c=depth,levels-2
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          if(allocated(work%registers(c)%fields(p)%coarse_delta)) work%registers(c)%fields(p)%coarse_delta=0
          if(allocated(work%registers(c)%fields(p)%fine_delta)) work%registers(c)%fields(p)%fine_delta=0
          if(present(internal_energy)) then
            if(allocated(work%held_registers(c)%fields(p)%coarse_delta)) work%held_registers(c)%fields(p)%coarse_delta=0
            if(allocated(work%held_registers(c)%fields(p)%fine_delta)) work%held_registers(c)%fields(p)%fine_delta=0
          endif
        enddo
      enddo
      if(present(heat_release)) then
        allocate(work%heat,stat=allocation_status)
        call agree(allocation_status==0,comm,valid)
        if(.not.valid) return
        work%heat=history(first)%clock_heat
      endif
      if(routing) then
        allocate(work%routed(depth:levels-1),stat=allocation_status)
        call agree(allocation_status==0,comm,valid)
        if(.not.valid) return
        do c=depth,levels-1
          work%routed(c)=auxiliary_pending(c)
          do p=1,size(tree%nodes)
            if(owners(p)==rank) work%routed(c)%fields(p)%values=0
          enddo
        enddo
      endif
      work%first=first;work%last=last;work%cursor=first-1;work%depth=depth;work%weight=blend
      work%records=history(first:last)
      call advance_prepared_replay_level(work,depth,history(first)%time,history(first)%interval,faces,scalar_faces,valid)
      nullify(prepared_context)
      if(.not.valid) return
      call agree(work%cursor==last,comm,valid)
      if(.not.valid) return
      chemical_trial=0;wall_trial=0
      do r=1,history_count
        if(r>=first.and.r<=last) then
          chemical_trial=chemical_trial+work%records(r-first+1)%chemical_change
          wall_trial=wall_trial+work%records(r-first+1)%wall_change
        else
          chemical_trial=chemical_trial+history(r)%chemical_change
          wall_trial=wall_trial+history(r)%wall_change
        endif
      enddo
      call agree(all(ieee_is_finite(chemical_trial)).and.all(ieee_is_finite(wall_trial)),comm,valid)
      if(.not.valid) return
      allocate(result%chemical,source=chemical_trial,stat=allocation_status)
      call agree(allocation_status==0,comm,valid)
      if(.not.valid) return
      allocate(result%wall,source=wall_trial,stat=allocation_status)
      call agree(allocation_status==0,comm,valid)
      if(.not.valid) return
      call move_alloc(work%gas%fields,result%state%fields)
      call move_alloc(work%temp%fields,result%temperature%fields)
      call move_alloc(work%reaction%fields,result%reaction_source%fields)
      call move_alloc(work%starts%fields,result%starts%fields)
      call move_alloc(work%ends%fields,result%ends%fields)
      call move_alloc(work%ts%fields,result%ts%fields)
      call move_alloc(work%te%fields,result%te%fields)
      call move_alloc(work%source0%fields,result%source0%fields)
      call move_alloc(work%source1%fields,result%source1%fields)
      call move_alloc(work%heat,result%heat)
      call move_alloc(work%records,result%records)
      call move_alloc(work%registers,result%registers)
      if(present(internal_energy)) then
        call agree(ieee_is_finite(work%energy_correction),comm,valid)
        if(.not.valid) return
        call move_alloc(work%held%fields,result%internal_energy%fields)
        call move_alloc(work%held_starts%fields,result%held_starts%fields)
        call move_alloc(work%held_ends%fields,result%held_ends%fields)
        call move_alloc(work%held_registers,result%held_registers)
        call move_alloc(scalar_faces,result%held_flux)
        result%energy_correction=work%energy_correction
      endif
      call move_alloc(faces,result%flux)
      good=.true.
    end subroutine build_prepared_replay_trial

    subroutine close_prepared_level(gas,temp,faces,scalar_faces,depth,good,paired,scalar_change)
      type(amr_tree_state_3d), intent(inout) :: gas,temp
      type(eb_integrated_flux_3d), intent(inout) :: faces(:)
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: scalar_faces(:)
      integer, intent(in) :: depth
      logical, intent(out) :: good,paired
      real(dp), optional, intent(inout) :: scalar_change
      paired=.false.
      if(present(internal_energy).and.allocated(correction_work)) then
        if(depth+1<levels) then
          call synchronize_prepared_held_interval(gas,temp,held_candidate,faces,scalar_faces,correction_work,depth,good, &
            interface_registers(depth)%fields,held_interface_registers(depth)%fields,scalar_change)
        else
          call synchronize_prepared_held_interval(gas,temp,held_candidate,faces,scalar_faces,correction_work,depth,good, &
            scalar_change=scalar_change)
        endif
        paired=good
      else if(depth+1<levels) then
        call synchronize_mpi_eb_tree_register_flux_3d(species,tree,owners,geometry,register_geometry, &
          gas,temp,faces,interface_registers(depth)%fields,depth,comm,good,boundary_faces,frame)
      else
        call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,gas,temp,faces,comm,good, &
          minimum_level=depth,maximum_level=depth,boundary_faces=boundary_faces,frame=frame)
      endif
    end subroutine close_prepared_level

    ! A held interval closes its gas and UEINT families as one transaction.
    ! Assemble signed gas increments before EOS admission; the held endpoint
    ! repair is then based on the matching scalar rereflux. Gas-only donor
    ! repair would change only one family and is deliberately not used here.
    subroutine synchronize_prepared_held_interval(gas,temp,stored,faces,scalar_faces,prefix,depth,good, &
        registers,scalar_registers,scalar_change)
      type(amr_tree_state_3d), intent(inout) :: gas,temp,stored
      type(eb_integrated_flux_3d), intent(inout) :: faces(:),scalar_faces(:)
      real(dp), intent(inout) :: prefix
      real(dp), optional, intent(inout) :: scalar_change
      integer, intent(in) :: depth
      logical, intent(out) :: good
      type(state_redist_register_3d), optional, intent(inout) :: registers(:),scalar_registers(:)
      call agree(present(registers).eqv.present(scalar_registers),comm,good)
      if(.not.good) return
      if(present(registers)) then
        call close_held_prepared_mpi_eb_tree_3d(species,tree,owners,geometry,gas,temp,stored,faces,scalar_faces, &
          held_policy,prefix,depth,comm,good,boundary_faces,frame,register_geometry,registers,scalar_registers, &
          internal_energy_change=scalar_change)
      else
        call close_held_prepared_mpi_eb_tree_3d(species,tree,owners,geometry,gas,temp,stored,faces,scalar_faces, &
          held_policy,prefix,depth,comm,good,boundary_faces,frame,register_geometry=extended_geometry, &
          internal_energy_change=scalar_change)
      endif
    end subroutine synchronize_prepared_held_interval

    subroutine prepared_predictor_adapter(rate,interval,iteration,depth,halo,comm,ok)
      type(amr_tree_state_3d), intent(in) :: rate
      real(dp), intent(in) :: interval
      integer, intent(in) :: iteration,depth
      type(amr_tree_state_3d), intent(inout) :: halo
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
      real(dp) :: first(levels-1),last(levels-1)
      first=0;last=0
      associate(w=>prepared_context)
      first(:depth)=max(0._dp,min(1._dp,(w%time-w%clock_start(:depth))/w%clock_dt(:depth)))
      last(:depth)=max(first(:depth),min(1._dp,(w%time+interval-w%clock_start(:depth))/w%clock_dt(:depth)))
      call predictor_builder(rate,interval,iteration,depth,w%time,halo,comm,ok, &
        w%source0,w%source1,.5_dp*(first+last))
      if(.not.ok) prepared_provider_rejected=.true.
      end associate
    end subroutine prepared_predictor_adapter

    subroutine prepared_diffusion_adapter(state,temperature,geometry,old_halos,interval,depth, &
        rate,spatial_rate,flux,theta,comm,ok,frame)
      type(amr_tree_state_3d), intent(in) :: state,temperature,old_halos
      type(eb_geometry_3d), intent(in) :: geometry(:)
      real(dp), intent(in) :: interval
      integer, intent(in) :: depth
      type(amr_tree_state_3d), intent(inout) :: rate,spatial_rate
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: theta
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
      type(rotating_frame_3d), optional, intent(in) :: frame
      call evaluate_spatial_diffusion(state,temperature,geometry,old_halos,interval,depth, &
        rate,spatial_rate,flux,theta,comm,ok,frame)
      if(.not.ok) prepared_provider_rejected=.true.
    end subroutine prepared_diffusion_adapter

    subroutine evaluate_prepared_replay_step(work,level,interval,fstart,fend,frozen,record,op,np, &
        step_flux,registers,scalar_state,scalar_starts,scalar_ends,scalar_correction,step_scalar_flux,scalar_registers, &
        aux,b,w,ch,ext,limiter,valid,refused,frozen_source,replay_diffusion)
      type(prepared_replay_work), intent(inout) :: work
      integer, intent(in) :: level
      real(dp), intent(in) :: interval,fstart(:),fend(:)
      type(mpi_eb_mol_interval_3d), intent(in) :: frozen
      type(mpi_eb_mol_interval_3d), intent(inout) :: record
      type(amr_tree_state_3d), intent(inout) :: op,np
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: step_flux(:),step_scalar_flux(:)
      type(state_redist_register_3d), allocatable, intent(inout) :: registers(:),scalar_registers(:)
      type(amr_tree_state_3d), allocatable, intent(inout) :: scalar_state,scalar_starts,scalar_ends
      real(dp), allocatable, intent(inout) :: scalar_correction,aux(:)
      real(dp), intent(inout) :: b(:),w(:),ch(:),ext(:),limiter
      logical, intent(out) :: valid,refused
      procedure(mpi_eb_sdc_auxiliary_source), optional :: frozen_source
      procedure(mpi_eb_sdc_spatial_diffusion_source), optional :: replay_diffusion
      if(present(internal_energy)) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,work%gas,work%temp, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,valid,work%reaction,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,work%starts,work%ends,work%ts,work%te,fstart,fend,work%source0,work%source1,step_flux,op,np, &
        redistribution=method,target_volume_fraction=redist_target,auxiliary_source=frozen_source,auxiliary_change=aux, &
        gas_acceleration=gas_acceleration,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame,state_redist_max_order=redist_order, &
        state_redist_registers=registers,evaluate_spatial_diffusion=replay_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width,constant_coefficients=constant_coefficients, &
        heat_release=work%heat,interval_history=record,replay_interval=frozen,replay_spatial_weight=work%weight, &
        provider_refused=refused,internal_energy=scalar_state,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=scalar_starts,ancestor_internal_energy_end=scalar_ends, &
        held_state_redist_registers=scalar_registers,energy_policy=energy_policy,energy_correction=scalar_correction, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,work%gas,work%temp, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,valid,work%reaction,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,work%starts,work%ends,work%ts,work%te,fstart,fend,work%source0,work%source1,step_flux,op,np, &
        redistribution=method,target_volume_fraction=redist_target,auxiliary_source=frozen_source,auxiliary_change=aux, &
        gas_acceleration=gas_acceleration,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame,state_redist_max_order=redist_order, &
        state_redist_registers=registers,evaluate_spatial_diffusion=replay_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width,constant_coefficients=constant_coefficients, &
        heat_release=work%heat,interval_history=record,replay_interval=frozen,replay_spatial_weight=work%weight, &
        provider_refused=refused, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      endif
    end subroutine evaluate_prepared_replay_step

    recursive subroutine advance_prepared_replay_level(work,level,time,interval,flux,scalar_flux,good)
      type(prepared_replay_work), target, intent(inout) :: work
      integer, intent(in) :: level
      real(dp), intent(in) :: time,interval
      type(eb_integrated_flux_3d), allocatable, intent(out) :: flux(:),scalar_flux(:)
      logical, intent(out) :: good
      type(mpi_eb_mol_interval_3d) :: frozen,record,clock_record
      type(amr_tree_state_3d) :: op,np,closing,closing_t,closing_held,reset_delta
      type(eb_integrated_flux_3d), allocatable :: step_flux(:),child_flux(:),closing_flux(:)
      type(eb_integrated_flux_3d), allocatable :: step_scalar_flux(:),child_scalar_flux(:),closing_scalar_flux(:)
      type(state_redist_register_3d), allocatable :: registers(:),closing_registers(:)
      type(state_redist_register_3d), allocatable :: scalar_registers(:),closing_scalar_registers(:)
      type(amr_tree_state_3d), allocatable :: scalar_state,scalar_starts,scalar_ends
      real(dp), allocatable :: scalar_correction
      real(dp) :: reset_change,closing_correction,clock_start_reset,clock_end_reset,source_reset,accepted_start
      type(state_redist_register_3d) :: part
      real(dp), allocatable :: aux(:)
      real(dp) :: fstart(levels-1),fend(levels-1),b(nv),w(nv),ch(nv),ext(nv),limiter,elapsed,sub_dt,tolerance
      integer :: r,index,p,n(3),step,target,bucket
      logical :: valid,agreed,refused
      good=.false.;clock_start_reset=0;clock_end_reset=0;source_reset=0
      r=work%cursor+1
      if(r>work%last) return
      tolerance=32*epsilon(1._dp)*max(abs(time+interval),tiny(1._dp))
      call agree(history(r)%depth==level.and.abs(history(r)%time-time)<=tolerance.and. &
        abs(history(r)%interval-interval)<=tolerance,comm,valid)
      if(.not.valid) return
      work%cursor=r;index=r-work%first+1
      ! Match the normal driver's canonical reset before capturing a complete
      ! clock. Retained descendant-tail resets never enter the replacement.
      if(present(internal_energy)) then
        call reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,work%gas,work%temp,work%held,reset_delta, &
          held_policy,level,reset_change,comm,valid,frame,internal_energy_change=clock_start_reset, &
          reset_level=merge(-2,-1,godunov_mol.and.chemistry.and.method=='none'))
        if(.not.valid) return
        work%energy_correction=work%energy_correction+reset_change
      endif
      work%clock_start(level+1)=time;work%clock_dt(level+1)=interval
      fstart=0;fend=0
      fstart(:level)=max(0._dp,min(1._dp,(time-work%clock_start(:level))/work%clock_dt(:level)))
      fend(:level)=max(fstart(:level),min(1._dp,(time+interval-work%clock_start(:level))/work%clock_dt(:level)))
      clock_record%clock_initial=work%gas;clock_record%clock_temperature=work%temp
      clock_record%clock_reaction_source=work%reaction
      if(present(heat_release)) clock_record%clock_heat=work%heat
      if(level>0) clock_record%open_ancestor_registers=work%registers(level-1)%fields
      if(present(internal_energy)) then
        clock_record%clock_internal_energy=work%held
        clock_record%clock_energy_correction=work%energy_correction
        if(level>0) clock_record%open_ancestor_held_registers=work%held_registers(level-1)%fields
        allocate(scalar_state,scalar_starts,scalar_ends,scalar_correction,stat=target)
        call agree(target==0,comm,valid)
        if(.not.valid) return
        scalar_state=work%held;scalar_starts=work%held_starts;scalar_ends=work%held_ends
        scalar_correction=work%energy_correction
      endif
      frozen=history(r)
      do p=1,size(tree%nodes)
        if(owners(p)/=rank.or.tree%nodes(p)%level/=level) cycle
        work%starts%fields(p)=work%gas%fields(p);work%ts%fields(p)=work%temp%fields(p)
        frozen%initial_reaction_source%fields(p)=work%reaction%fields(p)
        if(present(internal_energy)) then
          work%held_starts%fields(p)=work%held%fields(p)
          scalar_starts%fields(p)=work%held%fields(p)
        endif
      enddo
      if(has_auxiliary) allocate(aux(naux))
      work%time=time;prepared_context=>work
      if(present(predictor_builder)) return ! Native prepared MOL forbids a custom predictor.
      ! A null procedure pointer is not an absent optional callback in the
      ! checked compiler path. Supply only the actual active provider families.
      if(has_auxiliary) then
        if(present(evaluate_spatial_diffusion)) then
          call evaluate_prepared_replay_step(work,level,interval,fstart,fend,frozen,record,op,np,step_flux,registers, &
          scalar_state,scalar_starts,scalar_ends,scalar_correction,step_scalar_flux,scalar_registers, &
          aux,b,w,ch,ext,limiter,valid,refused, &
            frozen_source=source_adapter,replay_diffusion=prepared_diffusion_adapter)
        else
          call evaluate_prepared_replay_step(work,level,interval,fstart,fend,frozen,record,op,np,step_flux,registers, &
          scalar_state,scalar_starts,scalar_ends,scalar_correction,step_scalar_flux,scalar_registers, &
          aux,b,w,ch,ext,limiter,valid,refused,frozen_source=source_adapter)
        endif
      else
        if(present(evaluate_spatial_diffusion)) then
          call evaluate_prepared_replay_step(work,level,interval,fstart,fend,frozen,record,op,np,step_flux,registers, &
          scalar_state,scalar_starts,scalar_ends,scalar_correction,step_scalar_flux,scalar_registers, &
          aux,b,w,ch,ext,limiter,valid,refused, &
            replay_diffusion=prepared_diffusion_adapter)
        else
          call evaluate_prepared_replay_step(work,level,interval,fstart,fend,frozen,record,op,np,step_flux,registers, &
          scalar_state,scalar_starts,scalar_ends,scalar_correction,step_scalar_flux,scalar_registers, &
          aux,b,w,ch,ext,limiter,valid,refused)
        endif
      endif
      if(refused) prepared_provider_rejected=.true.
      if(.not.valid) return
      if(present(internal_energy)) then
        call move_alloc(scalar_state%fields,work%held%fields)
        work%energy_correction=scalar_correction
        record%old_internal_energy_flux=history(r)%old_internal_energy_flux
        record%clock_internal_energy=clock_record%clock_internal_energy
        record%clock_energy_correction=clock_record%clock_energy_correction
        call compute_reactive_eb_held_projection(0._dp,clock_start_reset,0._dp, &
          record%initial_internal_energy_change,accepted_start,valid)
        if(.not.valid) return
        if(r==work%first) then
          call compute_reactive_eb_held_projection(0._dp,history(r)%held_clock_start_change,0._dp, &
            accepted_start,source_reset,valid)
          if(.not.valid) return
          accepted_start=source_reset
        endif
        record%held_clock_start_change=accepted_start;record%held_clock_end_change=0
        record%held_clock_events_closed=.false.;record%held_clock_prefixes_valid=.false.
        record%held_clock_projection_before=0;record%held_clock_reset_before=0
        record%held_clock_projection_after=0;record%held_clock_reset_after=0
        call move_alloc(clock_record%open_ancestor_held_registers,record%open_ancestor_held_registers)
      endif
      record%time=time;record%old_flux=history(r)%old_flux
      record%clock_initial=clock_record%clock_initial;record%clock_temperature=clock_record%clock_temperature
      record%clock_reaction_source=clock_record%clock_reaction_source;record%clock_heat=clock_record%clock_heat
      call move_alloc(clock_record%open_ancestor_registers,record%open_ancestor_registers)
      if(routing) then
        record%routed_source=history(r)%routed_source
        record%routed_auxiliary=.true.
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          target=tree%nodes(p)%level;bucket=min(level,target)
          if(bucket<work%depth) cycle
          work%routed(bucket)%fields(p)%values=work%routed(bucket)%fields(p)%values+ &
            record%routed_source%fields(p)%values
        enddo
      endif
      call move_mol_interval(record,work%records(index))
      call allocate_replay_flux(flux,valid)
      if(.not.valid) return
      if(present(internal_energy)) then
        call allocate_replay_flux(scalar_flux,valid,1)
        if(.not.valid) return
      endif
      do p=1,size(tree%nodes)
        if(owners(p)/=rank.or.tree%nodes(p)%level/=level) cycle
        n=tree%nodes(p)%shape
        flux(p)=step_flux(p)
        work%ends%fields(p)=work%gas%fields(p);work%te%fields(p)=work%temp%fields(p)
        if(present(internal_energy)) then
          scalar_flux(p)=step_scalar_flux(p)
          work%held_ends%fields(p)=work%held%fields(p)
        endif
        work%source0%fields(p)=op%fields(p);work%source1%fields(p)=np%fields(p)
        if(level+1<levels) then
          part=registers(p)
          if(allocated(part%fine_delta)) part%fine_delta=0
          call accumulate_state_redist_register_3d(work%registers(level)%fields(p),part,nv,n+12,valid)
          if(.not.valid) exit
          if(present(internal_energy)) then
            part=scalar_registers(p)
            if(allocated(part%fine_delta)) part%fine_delta=0
            call accumulate_state_redist_register_3d(work%held_registers(level)%fields(p),part,1,n+12,valid)
            if(.not.valid) exit
          endif
        endif
        if(level>0) then
          part=registers(p)
          if(allocated(part%coarse_delta)) part%coarse_delta=0
          call accumulate_state_redist_register_3d(work%registers(level-1)%fields(p),part,nv,n+12,valid)
          if(.not.valid) exit
          if(present(internal_energy)) then
            part=scalar_registers(p)
            if(allocated(part%coarse_delta)) part%coarse_delta=0
            call accumulate_state_redist_register_3d(work%held_registers(level-1)%fields(p),part,1,n+12,valid)
            if(.not.valid) exit
          endif
        endif
      enddo
      call agree(valid,comm,agreed)
      valid=agreed
      if(.not.valid) return
      if(level+1<levels) then
        elapsed=0
        do step=1,ratios(level+2)
          sub_dt=interval/real(ratios(level+2),dp)
          if(step==ratios(level+2)) sub_dt=interval-elapsed
          call advance_prepared_replay_level(work,level+1,time+elapsed,sub_dt,child_flux,child_scalar_flux,valid)
          if(.not.valid) return
          do p=1,size(tree%nodes)
            ! Match the full descendant face family of the native recursion.
            if(owners(p)/=rank.or.tree%nodes(p)%level<=level) cycle
            flux(p)%x=flux(p)%x+child_flux(p)%x
            flux(p)%y=flux(p)%y+child_flux(p)%y
            flux(p)%z=flux(p)%z+child_flux(p)%z
            if(present(internal_energy)) then
              scalar_flux(p)%x=scalar_flux(p)%x+child_scalar_flux(p)%x
              scalar_flux(p)%y=scalar_flux(p)%y+child_scalar_flux(p)%y
              scalar_flux(p)%z=scalar_flux(p)%z+child_scalar_flux(p)%z
            endif
          enddo
          elapsed=elapsed+sub_dt
        enddo
      endif
      if(level==work%depth) then
        ! Admit the root closure privately. Its caller still owns publication
        ! and applies the retained interface/routed ledgers exactly once.
        closing=work%gas;closing_t=work%temp;closing_flux=flux
        if(present(internal_energy)) then
          closing_held=work%held;closing_correction=work%energy_correction;closing_scalar_flux=scalar_flux
          if(level+1<levels) then
            closing_registers=work%registers(level)%fields
            closing_scalar_registers=work%held_registers(level)%fields
            call synchronize_prepared_held_interval(closing,closing_t,closing_held,closing_flux,closing_scalar_flux, &
              closing_correction,level,valid,closing_registers,closing_scalar_registers)
          else
            call synchronize_prepared_held_interval(closing,closing_t,closing_held,closing_flux,closing_scalar_flux, &
              closing_correction,level,valid)
          endif
        else if(level+1<levels) then
          closing_registers=work%registers(level)%fields
          call synchronize_mpi_eb_tree_register_flux_3d(species,tree,owners,geometry,register_geometry, &
            closing,closing_t,closing_flux,closing_registers,level,comm,valid,boundary_faces,frame)
        else
          call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,closing,closing_t,closing_flux,comm,valid, &
            minimum_level=level,maximum_level=level,boundary_faces=boundary_faces,frame=frame)
        endif
        if(.not.valid) return
        if(routing) then
          if(held_auxiliary) then
            call apply_routed_held_source(closing,closing_t,closing_held,closing_correction,work%routed(level),valid)
          else
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            closing%fields(p)%values=closing%fields(p)%values+work%routed(level)%fields(p)%values
          enddo
          call recover_routed_replay(closing,closing_t,valid)
          endif
          if(.not.valid) return
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            valid=valid.and.all(work%routed(level)%fields(p)%values==auxiliary_pending(level)%fields(p)%values)
          enddo
          call agree(valid,comm,agreed)
      valid=agreed
          if(.not.valid) return
        endif
      else
        if(present(internal_energy)) then
          if(level+1<levels) then
            call synchronize_prepared_held_interval(work%gas,work%temp,work%held,flux,scalar_flux, &
              work%energy_correction,level,valid,work%registers(level)%fields,work%held_registers(level)%fields, &
              scalar_change=clock_end_reset)
          else
            call synchronize_prepared_held_interval(work%gas,work%temp,work%held,flux,scalar_flux, &
              work%energy_correction,level,valid,scalar_change=clock_end_reset)
          endif
        else if(level+1<levels) then
          call synchronize_mpi_eb_tree_register_flux_3d(species,tree,owners,geometry,register_geometry, &
            work%gas,work%temp,flux,work%registers(level)%fields,level,comm,valid,boundary_faces,frame)
        else
          call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,work%gas,work%temp,flux,comm,valid, &
            minimum_level=level,maximum_level=level,boundary_faces=boundary_faces,frame=frame)
        endif
        if(.not.valid) return
        if(routing) then
          if(held_auxiliary) then
            call apply_routed_held_source(work%gas,work%temp,work%held,work%energy_correction,work%routed(level),valid, &
              scalar_change=source_reset)
            if(.not.valid) return
            call compute_reactive_eb_held_projection(0._dp,clock_end_reset,0._dp,source_reset,accepted_start,valid)
            if(.not.valid) return
            clock_end_reset=accepted_start
            do p=1,size(tree%nodes)
              if(owners(p)==rank) work%routed(level)%fields(p)%values=0
            enddo
          else
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            work%gas%fields(p)%values=work%gas%fields(p)%values+work%routed(level)%fields(p)%values
            work%routed(level)%fields(p)%values=0
          enddo
          call recover_routed_replay(work%gas,work%temp,valid)
          if(.not.valid) return
          call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,work%gas,work%temp,comm,valid,frame)
          endif
          if(.not.valid) return
        endif
      endif
      if(present(internal_energy).and.level/=work%depth) then
        work%records(index)%held_clock_end_change=clock_end_reset
        work%records(index)%held_clock_events_closed=.true.
      endif
      ! Root trial admission used private copies. Its actual exit is recorded
      ! by the surrounding parent closure only after this trial is adopted.
      if(level+1<levels) then
        call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,work%reaction,nv-4,comm,valid,minimum_depth=level)
        if(.not.valid) return
        if(present(heat_release)) then
          call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,work%heat,1,comm,valid,minimum_depth=level)
          if(.not.valid) return
        endif
      endif
      call copy_mol_level_flux(tree,owners,rank,level,flux, &
        work%records(index)%synchronized_flux,comm,valid)
      if(.not.valid) return
      good=.true.
    end subroutine advance_prepared_replay_level

    subroutine reconstruct_flux_redist_replay(record,canonical,rates,spatial,transfer,good)
      integer, intent(in) :: record
      type(eb_integrated_flux_3d), intent(in) :: canonical(:)
      type(amr_tree_state_3d), intent(out) :: rates,spatial,transfer
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: difference,redistributed,exchange
      real(dp) :: impulse(nv),volume
      integer :: p,x,y,z,lev
      logical :: valid
      good=.false.;lev=history(record)%depth
      difference=source0
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        difference%fields(p)%values=0
        if(tree%nodes(p)%level/=lev) cycle
        associate(g=>geometry(p),f=>canonical(p))
        do z=1,g%nz
        do y=1,g%ny
        do x=1,g%nx
          if(g%volume_fraction(x,y,z)<=0) cycle
          volume=g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz
          impulse=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(:,x-1,y,z)- &
            g%x_face_fraction(x,y,z)*f%x(:,x,y,z))+ &
            g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(:,x,y-1,z)- &
            g%y_face_fraction(x,y,z)*f%y(:,x,y,z))+ &
            g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(:,x,y,z-1)- &
            g%z_face_fraction(x,y,z)*f%z(:,x,y,z))
          difference%fields(p)%values(:,x,y,z)=impulse/volume/history(record)%interval+ &
            history(record)%non_face_spatial_rate%fields(p)%values(:,x,y,z)
        enddo
        enddo
        enddo
        end associate
      enddo
      call redistribute_sdc_mpi_eb_tree_level_rate_3d(tree,owners,geometry,difference,nv,lev, &
        history(record)%interval,redistributed,exchange,comm,valid,boundary_faces)
      if(.not.valid) return
      rates=source0;spatial=history(record)%spatial_rate;transfer=exchange;valid=.true.
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        rates%fields(p)%values=0
        valid=valid.and.all(ieee_is_finite(transfer%fields(p)%values))
        if(.not.valid) exit
        if(tree%nodes(p)%level/=lev) cycle
        spatial%fields(p)%values=redistributed%fields(p)%values
        call post_eb_spatial_density_rate_3d(geometry(p),spatial%fields(p)%values,valid)
        if(.not.valid) exit
        rates%fields(p)%values=history(record)%volume_source%fields(p)%values+spatial%fields(p)%values
        valid=valid.and.all(ieee_is_finite(rates%fields(p)%values))
        if(.not.valid) exit
      enddo
      call agree(valid,comm,good)
    end subroutine reconstruct_flux_redist_replay

    subroutine replay_reactor_adapter(species,reactions,energy,source_species,source_energy,duration,rtol,atol, &
        density,temperature,steps,ok,maximum_steps)
      type(nasa7_species), intent(in) :: species(:)
      type(elementary_reaction), intent(in) :: reactions(:)
      real(dp), intent(in) :: energy,source_species(:),source_energy,duration,rtol,atol
      real(dp), intent(inout) :: density(:),temperature
      integer, intent(out) :: steps
      logical, intent(out) :: ok
      integer, optional, intent(in) :: maximum_steps
      logical :: refused
      call evaluate_classified_reactor(evaluate_reaction,species,reactions,energy,source_species,source_energy, &
        duration,rtol,atol,density,temperature,steps,ok,refused,maximum_steps, &
          evaluate_admission=evaluate_reaction_admission,admission_status=admission_state)
      prepared_provider_rejected=prepared_provider_rejected.or.refused
    end subroutine replay_reactor_adapter

    subroutine build_mol_replay_trial(first,last,depth,blend,result,good)
      integer, intent(in) :: first,last,depth
      real(dp), intent(in) :: blend
      type(mol_replay_result), intent(out) :: result
      logical, intent(out) :: good
      type(replay_flux_family), allocatable :: family(:)
      type(mpi_eb_mol_interval_3d), allocatable :: records(:)
      type(eb_integrated_flux_3d), allocatable :: raw(:),canonical(:),result_flux(:)
      type(eb_integrated_flux_3d), allocatable :: held_raw(:),held_result_flux(:)
      type(amr_tree_state_3d), allocatable :: work_held
      type(amr_tree_state_3d) :: held_assembly,reset_delta
      real(dp) :: correction_trial,closure_reset,scalar_reset,accepted_start
      integer :: closing_record,entry
      type(amr_tree_state_3d) :: work,work_t,work_ir,old,old_t,rates,ledger,assembly_delta
      type(amr_tree_state_3d), allocatable :: work_heat
      type(amr_tree_state_3d), allocatable :: routed_pending(:),redist_pending(:)
      type(amr_tree_state_3d) :: replay_spatial,replay_transfer,closing_state,closing_temperature
      real(dp) :: chemical_trial(nv),amount(nv),scale,time_tolerance
      real(dp) :: volume,divergence(nv)
      integer :: r,q,p,lev,reference(0:levels-1),allocation_status,x,y,z,c,target,bucket,next_depth,closure
      logical :: valid,agreed,seen(0:levels-1)
      logical :: refused
      real(dp) :: trial_boundary(nv),trial_theta
      good=.false.
      ! Routed impulses are frozen volume sources, not gas face families.
      ! Check the complete paired layout before any replay queue is indexed.
      do r=first,last
        valid=(history(r)%auxiliary_has_internal_energy.eqv.held_auxiliary).and. &
          (history(r)%routed_auxiliary.eqv.routing).and.(history(r)%has_auxiliary_source.eqv.has_auxiliary)
        call agree(valid,comm,agreed);if(.not.agreed) return
        if(routing) then
          call validate_mpi_amr_tree_local_3d(tree,owners,history(r)%routed_source,naux,comm,valid)
          if(.not.valid) return
        endif
      enddo
      if(prepared_redist) then
        call build_prepared_replay_trial(first,last,depth,blend,result,good)
        return
      endif
      if(allocated(funding_context_work)) then
        ! Only a complete observed window can enter native replay. It restores
        ! the pre-payment anchor and closes each pending/open clock privately.
        result%state=candidate;result%temperature=ct;result%reaction_source=ir
        result%internal_energy=held_candidate;result%energy_correction=correction_work
        result%records=history(first:last)
        allocate(result%funding,source=funding_context_work)
        allocate(result%chemical(nv),result%wall(nv),result%source(nv))
        if(has_auxiliary) allocate(result%auxiliary(naux))
        if(present(heat_release)) allocate(result%heat,source=heat)
        result%chemical=0;result%wall=0;result%source=0
        if(has_auxiliary) result%auxiliary=0
        trial_boundary=0;trial_theta=1;refused=.false.
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,geometry, &
          result%state,result%temperature,history(first)%interval,iterations,chemistry,transport_flags,force,rtol,atol, &
          extended_geometry,boundary_faces,trial_boundary,result%wall,result%chemical,result%source,trial_theta, &
          comm,good,result%reaction_source,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          wall_temperature,redistribution=method,target_volume_fraction=redist_target, &
          auxiliary_source=auxiliary_source,auxiliary_change=result%auxiliary,route_auxiliary=route_auxiliary, &
          gas_acceleration=gas_acceleration,halo_source=halo_source,predictor_builder=predictor_builder, &
          evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
          embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
          mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame,state_redist_max_order=redist_order, &
          evaluate_spatial_diffusion=evaluate_spatial_diffusion,spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
          constant_coefficients=constant_coefficients,heat_release=result%heat,integrated_flux=result%flux, &
          internal_energy=result%internal_energy,internal_energy_flux=result%held_flux, &
          interval_history=result%records,energy_policy=energy_policy,energy_correction=result%energy_correction, &
          spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
          provider_refused=refused,source_funding_context=result%funding, &
          native_replay_history=history(first:last),native_spatial_weight=blend, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
        prepared_provider_rejected=prepared_provider_rejected.or.refused
        return
      endif
      allocate(family(first:last),stat=allocation_status)
      call agree(allocation_status==0,comm,valid)
      if(.not.valid) return
      records=history(first:last)
      if(present(internal_energy)) then
        records%held_clock_start_change=0;records%held_clock_end_change=0
        records%held_clock_events_closed=.false.;records%held_clock_prefixes_valid=.false.
      endif
      reference=0
      do r=first,last
        lev=history(r)%depth
        if(reference(lev)==0) reference(lev)=r
      enddo
      ! A common old spatial family across the parent interval is the low-order
      ! anchor. Later fine f0 values can already contain the fictional parent
      ! endpoint whose trace species the failed reflux was trying to remove.
      do r=last,first,-1
        lev=history(r)%depth;q=reference(lev)
        scale=history(r)%interval/history(q)%interval
        call allocate_replay_flux(raw,valid)
        if(.not.valid) return
        if(present(internal_energy)) then
          call allocate_replay_flux(held_raw,valid,1)
          if(.not.valid) return
        endif
        do p=1,size(tree%nodes)
          if(owners(p)/=rank.or.tree%nodes(p)%level/=lev) cycle
          if(blend==0) then
            raw(p)%x=scale*history(q)%old_flux(p)%x
            raw(p)%y=scale*history(q)%old_flux(p)%y
            raw(p)%z=scale*history(q)%old_flux(p)%z
            if(present(internal_energy)) then
              held_raw(p)%x=scale*history(q)%old_internal_energy_flux(p)%x
              held_raw(p)%y=scale*history(q)%old_internal_energy_flux(p)%y
              held_raw(p)%z=scale*history(q)%old_internal_energy_flux(p)%z
            endif
          else
            raw(p)%x=(1-blend)*scale*history(q)%old_flux(p)%x+blend*history(r)%spatial_flux(p)%x
            raw(p)%y=(1-blend)*scale*history(q)%old_flux(p)%y+blend*history(r)%spatial_flux(p)%y
            raw(p)%z=(1-blend)*scale*history(q)%old_flux(p)%z+blend*history(r)%spatial_flux(p)%z
            if(present(internal_energy)) then
              held_raw(p)%x=(1-blend)*scale*history(q)%old_internal_energy_flux(p)%x+ &
                blend*history(r)%internal_energy_flux(p)%x
              held_raw(p)%y=(1-blend)*scale*history(q)%old_internal_energy_flux(p)%y+ &
                blend*history(r)%internal_energy_flux(p)%y
              held_raw(p)%z=(1-blend)*scale*history(q)%old_internal_energy_flux(p)%z+ &
                blend*history(r)%internal_energy_flux(p)%z
            endif
          endif
        enddo
        do q=r+1,last
          if(history(q)%depth/=lev+1) cycle
          time_tolerance=32*epsilon(1._dp)*max(abs(history(r)%time+history(r)%interval),tiny(1._dp))
          if(history(q)%time<history(r)%time-time_tolerance.or. &
            history(q)%time+history(q)%interval>history(r)%time+history(r)%interval+time_tolerance) cycle
          do p=1,size(tree%nodes)
            if(owners(p)/=rank.or.tree%nodes(p)%level/=lev+1) cycle
            raw(p)%x=raw(p)%x+family(q)%faces(p)%x
            raw(p)%y=raw(p)%y+family(q)%faces(p)%y
            raw(p)%z=raw(p)%z+family(q)%faces(p)%z
            if(present(internal_energy)) then
              held_raw(p)%x=held_raw(p)%x+family(q)%held_faces(p)%x
              held_raw(p)%y=held_raw(p)%y+family(q)%held_faces(p)%y
              held_raw(p)%z=held_raw(p)%z+family(q)%held_faces(p)%z
            endif
          enddo
        enddo
        call assemble_mpi_eb_tree_reflux_3d(species,tree,owners,geometry,candidate,ct,raw,assembly_delta, &
          canonical,comm,valid,lev,min(lev+1,levels-1),boundary_faces,frame)
        if(.not.valid) return
        call copy_mol_level_flux(tree,owners,rank,lev,canonical,family(r)%faces,comm,valid)
        if(.not.valid) return
        if(present(internal_energy)) then
          ! Use the scalar's own linear peer reconciliation and fine-face
          ! restriction. The scratch field is discarded: replay applies the
          ! resulting divergence once through the reconstructed spatial rate.
          held_assembly=held_candidate
          do p=1,size(tree%nodes)
            if(owners(p)==rank) held_assembly%fields(p)%values=0
          enddo
          call synchronize_mpi_eb_tree_scalar_flux_3d(species,tree,owners,geometry,held_assembly,held_raw, &
            comm,valid,minimum_level=lev,maximum_level=min(lev+1,levels-1),boundary_faces=boundary_faces)
          if(.not.valid) return
          call copy_mol_level_flux(tree,owners,rank,lev,held_raw,family(r)%held_faces,comm,valid)
          if(.not.valid) return
        endif
      enddo
      work=candidate;work_t=ct;work_ir=ir
      if(present(internal_energy)) then
        ! Restore the entire clock, including ancestor changes, and the numeric
        ! prefix. Previously accepted descendants must not remain charged.
        call validate_mpi_amr_tree_local_3d(tree,owners,history(first)%clock_initial,nv,comm,valid)
        if(.not.valid) return
        call validate_mpi_amr_tree_local_3d(tree,owners,history(first)%clock_temperature,1,comm,valid)
        if(.not.valid) return
        call validate_mpi_amr_tree_local_3d(tree,owners,history(first)%clock_reaction_source,nv-4,comm,valid)
        if(.not.valid) return
        call validate_mpi_amr_tree_local_3d(tree,owners,history(first)%clock_internal_energy,1,comm,valid)
        if(.not.valid) return
        call agree(ieee_is_finite(history(first)%clock_energy_correction),comm,valid)
        if(.not.valid) return
        allocate(work_held,stat=allocation_status)
        call agree(allocation_status==0,comm,valid)
        if(.not.valid) return
        work=history(first)%clock_initial;work_t=history(first)%clock_temperature
        work_ir=history(first)%clock_reaction_source;work_held=history(first)%clock_internal_energy
        correction_trial=history(first)%clock_energy_correction
      endif
      if(present(heat_release)) then
        allocate(work_heat,stat=allocation_status)
        call agree(allocation_status==0,comm,valid)
        if(.not.valid) return
        work_heat=heat
        if(present(internal_energy)) then
          call validate_mpi_amr_tree_local_3d(tree,owners,history(first)%clock_heat,1,comm,valid)
          if(.not.valid) return
          work_heat=history(first)%clock_heat
        endif
      endif
      if(routing) then
        allocate(routed_pending(depth:levels-1),stat=allocation_status)
        call agree(allocation_status==0,comm,valid)
        if(.not.valid) return
        do closure=depth,levels-1
          routed_pending(closure)=candidate
          do p=1,size(tree%nodes)
            if(owners(p)==rank) routed_pending(closure)%fields(p)%values=0
          enddo
        enddo
      endif
      if(method=='flux_redist') then
        allocate(redist_pending(0:levels-1),stat=allocation_status)
        call agree(allocation_status==0,comm,valid)
        if(.not.valid) return
        do closure=0,levels-1
          redist_pending(closure)=source0
          do p=1,size(tree%nodes)
            if(owners(p)==rank) redist_pending(closure)%fields(p)%values=0
          enddo
        enddo
        ! Rebuild still-open ancestor buckets from their current clock windows.
        ! A subtract/add replacement can erase a resolved trace impulse under
        ! a large previous total. Completed sibling windows must stay consumed.
        do r=1,history_count
          if(r>=first.and.r<=last) cycle
          lev=history(r)%depth
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            target=tree%nodes(p)%level;bucket=min(lev,target)
            if(bucket>=depth) cycle
            time_tolerance=32*epsilon(1._dp)*max(abs(clock_start(bucket+1)+clock_dt(bucket+1)),tiny(1._dp))
            if(history(r)%time<clock_start(bucket+1)-time_tolerance.or. &
              history(r)%time+history(r)%interval>clock_start(bucket+1)+clock_dt(bucket+1)+time_tolerance) cycle
            redist_pending(bucket)%fields(p)%values=redist_pending(bucket)%fields(p)%values+ &
              history(r)%redistribution_transfer%fields(p)%values
            if(allocated(history(r)%synchronization_transfer%fields)) &
              redist_pending(bucket)%fields(p)%values=redist_pending(bucket)%fields(p)%values+ &
                history(r)%synchronization_transfer%fields(p)%values
          enddo
        enddo
      endif
      seen=.false.
      if(.not.present(internal_energy)) then
      do r=first,last
        lev=history(r)%depth
        if(seen(lev)) cycle
        seen(lev)=.true.
        do p=1,size(tree%nodes)
          if(owners(p)/=rank.or.tree%nodes(p)%level/=lev) cycle
          work%fields(p)=history(r)%initial%fields(p)
          work_t%fields(p)=history(r)%initial_temperature%fields(p)
        enddo
      enddo
      endif
      do r=first,last
        lev=history(r)%depth
        if(routing) then
          call agree(allocated(history(r)%routed_source%fields),comm,valid)
          if(.not.valid) return
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            target=tree%nodes(p)%level;bucket=min(lev,target)
            if(bucket<depth) cycle
            routed_pending(bucket)%fields(p)%values=routed_pending(bucket)%fields(p)%values+ &
              history(r)%routed_source%fields(p)%values
          enddo
        endif
        if(present(internal_energy)) then
          if(present(evaluate_reaction)) then
            call replay_mpi_eb_mol_held_interval_3d(species,reactions,tree,owners,geometry,records(r-first+1), &
              family(r)%faces,family(r)%held_faces,work,work_t,work_held,work_ir,correction_trial,amount, &
              rtol,atol,comm,valid,evaluate_reaction=replay_reactor_adapter,frame=frame,heat_release=work_heat)
          else
            call replay_mpi_eb_mol_held_interval_3d(species,reactions,tree,owners,geometry,records(r-first+1), &
              family(r)%faces,family(r)%held_faces,work,work_t,work_held,work_ir,correction_trial,amount, &
              rtol,atol,comm,valid,frame=frame,heat_release=work_heat)
          endif
          call agree(.not.prepared_provider_rejected,comm,agreed)
          prepared_provider_rejected=.not.agreed
          if(prepared_provider_rejected.or..not.valid) return
          entry=r-first+1;accepted_start=records(entry)%initial_internal_energy_change
          if(r==first) then
            call compute_reactive_eb_held_projection(0._dp,history(first)%held_clock_start_change,0._dp, &
              accepted_start,scalar_reset,valid)
            if(.not.valid) return
            accepted_start=scalar_reset
          endif
          records(entry)%held_clock_start_change=accepted_start
        else
        if(method=='flux_redist') then
          call reconstruct_flux_redist_replay(r,family(r)%faces, &
            rates,replay_spatial,replay_transfer,valid)
          if(.not.valid) return
          records(r-first+1)%spatial_rate=replay_spatial
          records(r-first+1)%redistribution_transfer=replay_transfer
          records(r-first+1)%synchronization_transfer=source0
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            records(r-first+1)%synchronization_transfer%fields(p)%values=0
            target=tree%nodes(p)%level;bucket=min(lev,target)
            redist_pending(bucket)%fields(p)%values=redist_pending(bucket)%fields(p)%values+ &
              replay_transfer%fields(p)%values
          enddo
        else
        call reconstruct_mpi_eb_mol_interval_3d(tree,owners,geometry,history(r),family(r)%faces,nv,rates,comm,valid)
        if(.not.valid) return
        ! With no redistribution, density/species volume forcing is available
        ! directly. Use it rather than subtracting nearly equal trace fluxes.
        ! Momentum/energy retain the complete wall and volume residual above.
        do p=1,size(tree%nodes)
          if(owners(p)/=rank.or.tree%nodes(p)%level/=lev) cycle
          associate(g=>geometry(p),f=>family(r)%faces(p))
          do z=1,g%nz
          do y=1,g%ny
          do x=1,g%nx
            if(g%volume_fraction(x,y,z)<=0) cycle
            volume=g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz
            divergence=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(:,x-1,y,z)- &
              g%x_face_fraction(x,y,z)*f%x(:,x,y,z))+ &
              g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(:,x,y-1,z)- &
              g%y_face_fraction(x,y,z)*f%y(:,x,y,z))+ &
              g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(:,x,y,z-1)- &
              g%z_face_fraction(x,y,z)*f%z(:,x,y,z))
            do c=1,nv
              if(c>=2.and.c<=5) cycle
              rates%fields(p)%values(c,x,y,z)=history(r)%mass_source%fields(p)%values(c,x,y,z)+ &
                divergence(c)/volume/history(r)%interval
            enddo
          enddo
          enddo
          enddo
          end associate
        enddo
        endif
        old=work;old_t=work_t
        do p=1,size(owners)
          if(owners(p)/=rank.or.tree%nodes(p)%level/=lev) cycle
          records(r-first+1)%initial_reaction_source%fields(p)=work_ir%fields(p)
        enddo
        if(present(evaluate_reaction)) then
          call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old,old_t,rates, &
            history(r)%interval,rtol,atol,work,work_t,work_ir,comm,valid,target_depth=lev, &
            evaluate_reaction=replay_reactor_adapter,frame=frame,heat_release=work_heat)
        else
          call react_sdc_mpi_eb_tree_3d(species,reactions,tree,owners,geometry,old,old_t,rates, &
            history(r)%interval,rtol,atol,work,work_t,work_ir,comm,valid,target_depth=lev, &
            frame=frame,heat_release=work_heat)
        endif
        call agree(.not.prepared_provider_rejected,comm,agreed)
        prepared_provider_rejected=.not.agreed
        if(prepared_provider_rejected) return
        if(.not.valid) return
        ledger=rates
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          ledger%fields(p)%values=0
          if(tree%nodes(p)%level/=lev) cycle
          ledger%fields(p)%values=work%fields(p)%values- &
            (old%fields(p)%values+history(r)%interval*rates%fields(p)%values)
          records(r-first+1)%initial%fields(p)=old%fields(p)
          records(r-first+1)%initial_temperature%fields(p)=old_t%fields(p)
          records(r-first+1)%forcing%fields(p)=rates%fields(p)
          records(r-first+1)%local_state%fields(p)=work%fields(p)
          records(r-first+1)%local_temperature%fields(p)=work_t%fields(p)
          records(r-first+1)%reaction_source%fields(p)=work_ir%fields(p)
        enddo
        call integrate_mpi_eb_tree_3d(tree,owners,geometry,ledger,amount,comm,valid)
        if(.not.valid) return
        records(r-first+1)%chemical_change=amount
        records(r-first+1)%spatial_flux=family(r)%faces
        records(r-first+1)%synchronized_flux=family(r)%faces
        endif ! ordinary gas-only or paired held interval
        if(present(internal_energy)) then
          call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,work,work_t,comm,valid,frame,work_held)
        else
          call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,work,work_t,comm,valid,frame)
        endif
        if(.not.valid) return
        if(present(internal_energy)) then
          next_depth=depth
          if(r<last) next_depth=history(r+1)%depth
          do closure=lev,max(next_depth,depth+1),-1
            closing_record=0
            do q=r-first+1,1,-1
              if(records(q)%depth/=closure) cycle
              closing_record=q;exit
            enddo
            if(closing_record==0) return
            scalar_reset=0
            call reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,work,work_t,work_held,reset_delta, &
              held_policy,closure,closure_reset,comm,valid,frame,leaf_only=.true., &
              internal_energy_change=scalar_reset)
            if(.not.valid) return
            correction_trial=correction_trial+closure_reset
            call compute_reactive_eb_held_projection(0._dp,records(closing_record)%held_clock_end_change,0._dp, &
              scalar_reset,accepted_start,valid)
            if(.not.valid) return
            records(closing_record)%held_clock_end_change=accepted_start
            records(closing_record)%held_clock_events_closed=.true.
          enddo
        endif
        if(method=='flux_redist') then
          next_depth=depth
          if(r<last) next_depth=history(r+1)%depth
          do closure=lev,max(next_depth,depth+1),-1
            do p=1,size(tree%nodes)
              if(owners(p)/=rank) cycle
              work%fields(p)%values=work%fields(p)%values+redist_pending(closure)%fields(p)%values
              redist_pending(closure)%fields(p)%values=0
            enddo
            call recover_routed_replay(work,work_t,valid)
            if(.not.valid) return
            call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,work,work_t,comm,valid,frame)
            if(.not.valid) return
          enddo
        endif
        if(routing) then
          next_depth=depth
          if(r<last) next_depth=history(r+1)%depth
          do closure=lev,max(next_depth,depth+1),-1
            if(held_auxiliary) then
              closing_record=0
              do q=r-first+1,1,-1
                if(records(q)%depth/=closure) cycle
                closing_record=q;exit
              enddo
              if(closing_record==0) return
              scalar_reset=0
              call apply_routed_held_source(work,work_t,work_held,correction_trial,routed_pending(closure),valid, &
                scalar_change=scalar_reset)
              if(.not.valid) return
              call compute_reactive_eb_held_projection(0._dp,records(closing_record)%held_clock_end_change,0._dp, &
                scalar_reset,accepted_start,valid)
              if(.not.valid) return
              records(closing_record)%held_clock_end_change=accepted_start
              do p=1,size(tree%nodes)
                if(owners(p)==rank) routed_pending(closure)%fields(p)%values=0
              enddo
            else
            do p=1,size(tree%nodes)
              if(owners(p)/=rank) cycle
              work%fields(p)%values=work%fields(p)%values+routed_pending(closure)%fields(p)%values
              routed_pending(closure)%fields(p)%values=0
            enddo
            call recover_routed_replay(work,work_t,valid)
            if(.not.valid) return
            call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,work,work_t,comm,valid,frame)
            endif
            if(.not.valid) return
          enddo
        endif
      enddo
      if(routing) then
        ! The failed parent has not applied this bucket yet. Leave it to the
        ! original post-sync closure exactly once; interior buckets replayed
        ! above had already been consumed before the original rejection.
        valid=.true.
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          valid=valid.and.all(routed_pending(depth)%fields(p)%values==auxiliary_pending(depth)%fields(p)%values)
        enddo
        call agree(valid,comm,agreed)
        if(.not.agreed) return
      endif
      if(method=='flux_redist') then
        ! Admit the pending failed-parent closure too, while retaining it in
        ! the result for the original caller to apply once after acceptance.
        closing_state=work;closing_temperature=work_t
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          closing_state%fields(p)%values=closing_state%fields(p)%values+redist_pending(depth)%fields(p)%values
        enddo
        call recover_routed_replay(closing_state,closing_temperature,valid)
        if(.not.valid) return
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,closing_state,closing_temperature,comm,valid,frame)
        if(.not.valid) return
      endif
      call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,work_ir,nv-4,comm,valid,minimum_depth=depth)
      if(.not.valid) return
      if(present(heat_release)) then
        call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,work_heat,1,comm,valid,minimum_depth=depth)
        if(.not.valid) return
      endif
      call allocate_replay_flux(result_flux,valid)
      if(.not.valid) return
      if(present(internal_energy)) then
        call allocate_replay_flux(held_result_flux,valid,1)
        if(.not.valid) return
      endif
      do r=first,last
        lev=history(r)%depth
        ! Parent closure uses adjacent levels, but the published physical
        ! family must retain every accepted descendant interval.
        do p=1,size(tree%nodes)
          if(owners(p)/=rank.or.tree%nodes(p)%level/=lev) cycle
          result_flux(p)%x=result_flux(p)%x+family(r)%faces(p)%x
          result_flux(p)%y=result_flux(p)%y+family(r)%faces(p)%y
          result_flux(p)%z=result_flux(p)%z+family(r)%faces(p)%z
          if(present(internal_energy)) then
            held_result_flux(p)%x=held_result_flux(p)%x+family(r)%held_faces(p)%x
            held_result_flux(p)%y=held_result_flux(p)%y+family(r)%held_faces(p)%y
            held_result_flux(p)%z=held_result_flux(p)%z+family(r)%held_faces(p)%z
          endif
        enddo
      enddo
      ! Rebuild the interval ledger directly; subtracting the previous total
      ! can lose resolved trace increments when the bulk inventory is large.
      chemical_trial=0
      do r=1,history_count
        if(r>=first.and.r<=last) then
          chemical_trial=chemical_trial+records(r-first+1)%chemical_change
        else
          chemical_trial=chemical_trial+history(r)%chemical_change
        endif
      enddo
      call agree(all(ieee_is_finite(chemical_trial)),comm,valid)
      if(.not.valid) return
      allocate(result%chemical,source=chemical_trial,stat=allocation_status)
      call agree(allocation_status==0,comm,valid)
      if(.not.valid) return
      if(present(internal_energy)) then
        call agree(ieee_is_finite(correction_trial),comm,valid)
        if(.not.valid) return
        call move_alloc(work_held%fields,result%internal_energy%fields)
        result%energy_correction=correction_trial
        call move_alloc(held_result_flux,result%held_flux)
      endif
      if(method=='flux_redist') call move_alloc(redist_pending,result%redistribution_pending)
      call move_alloc(work%fields,result%state%fields)
      call move_alloc(work_t%fields,result%temperature%fields)
      call move_alloc(work_ir%fields,result%reaction_source%fields)
      call move_alloc(work_heat,result%heat)
      call move_alloc(records,result%records)
      call move_alloc(result_flux,result%flux)
      good=.true.
    end subroutine build_mol_replay_trial

    subroutine store_godunov_step(saved,clock,index,time,good)
      type(eb_godunov_interval_3d), intent(inout) :: saved,clock
      integer, intent(in) :: index
      real(dp), intent(in) :: time
      logical, intent(out) :: good
      logical :: local_ok
      good=.false.
      if(routing) then
        if(allocated(off_old%fields)) then
          call copy_godunov_field_3d(off_old,saved%routed_auxiliary_old_rate,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
        endif
        if(allocated(off_new%fields)) then
          call copy_godunov_field_3d(off_new,saved%routed_auxiliary_new_rate,local_ok)
      call agree(local_ok,comm,good);if(.not.good) return
        endif
      endif
        call move_alloc(clock%clock_gas%fields,saved%clock_gas%fields)
        call move_alloc(clock%clock_temperature%fields,saved%clock_temperature%fields)
        call move_alloc(clock%clock_internal_energy%fields,saved%clock_internal_energy%fields)
        call move_alloc(clock%clock_heat%fields,saved%clock_heat%fields)
        call move_alloc(clock%clock_reaction_source%fields,saved%clock_reaction_source%fields)
        call move_alloc(clock%open_ancestor_registers,saved%open_ancestor_registers)
        call move_alloc(clock%open_ancestor_scalar_registers,saved%open_ancestor_scalar_registers)
        call move_alloc(clock%clock_start,saved%clock_start)
        call move_alloc(clock%clock_interval,saved%clock_interval)
        call move_alloc(clock%clock_routed_pending,saved%clock_routed_pending)
        saved%routed_auxiliary=clock%routed_auxiliary
        saved%clock_energy_correction=clock%clock_energy_correction
        saved%clock_pressure_work_change=clock%clock_pressure_work_change
        call compute_reactive_eb_held_projection(0._dp,clock%held_clock_start_change,0._dp, &
          saved%initial_internal_energy_change,saved%held_clock_start_change,local_ok)
        call agree(local_ok,comm,good);if(.not.good) return
        saved%clock_minimum_theta=clock%clock_minimum_theta
        saved%parent_record=clock%parent_record;saved%time=time
        call move_godunov_interval_3d(saved,godunov_records(index))
      good=.true.
    end subroutine store_godunov_step
    subroutine complete_interval_faces(f,g,components,good)
      type(eb_integrated_flux_3d), intent(inout) :: f
      type(eb_geometry_3d), intent(in) :: g
      integer, intent(in) :: components
      logical, intent(out) :: good
      integer :: n(3),allocation_status
      good=.false.;n=[g%nx,g%ny,g%nz]
      if(.not.allocated(f%x)) then
        allocate(f%x(components,0:n(1),n(2),n(3)),source=0._dp,stat=allocation_status)
        if(allocation_status/=0) return
      endif
      if(.not.allocated(f%y)) then
        allocate(f%y(components,n(1),0:n(2),n(3)),source=0._dp,stat=allocation_status)
        if(allocation_status/=0) return
      endif
      if(.not.allocated(f%z)) then
        allocate(f%z(components,n(1),n(2),0:n(3)),source=0._dp,stat=allocation_status)
        if(allocation_status/=0) return
      endif
      good=valid_eb_integrated_flux_3d(f,g,components)
    end subroutine complete_interval_faces

    subroutine prepare_clock_pressure(depth,good)
      integer, intent(in) :: depth
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: pressure
      real(dp) :: q(reactive_nprim(size(species))),t,cs,delta,target
      integer :: p,x,y,z,m(3)
      logical :: admitted,local_ok,ready
      good=.false.
      if(depth==levels-1) then
        good=.true.;return
      endif
      pressure=held_candidate;local_ok=.true.
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        pressure%fields(p)%values=0
        if(tree%nodes(p)%level<depth) cycle
        m=tree%nodes(p)%shape
        do z=1,m(3)
        do y=1,m(2)
        do x=1,m(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          call reactive_conserved_to_primitive(species,candidate%fields(p)%values(:,x,y,z), &
            ct%fields(p)%values(1,x,y,z),q,t,cs,admitted, &
            internal_energy_density=held_candidate%fields(p)%values(1,x,y,z))
          if(.not.admitted) then
            local_ok=.false.;cycle
          endif
          pressure%fields(p)%values(1,x,y,z)=q(5)
        enddo
        enddo
        enddo
      enddo
      call agree(local_ok,comm,ready)
      if(.not.ready) return
      call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,pressure,1,comm,ready,minimum_depth=depth)
      if(.not.ready) return
      do p=1,size(owners)
        if(owners(p)/=rank.or.tree%nodes(p)%level/=depth) cycle
        m=tree%nodes(p)%shape
        do z=1,m(3)
        do y=1,m(2)
        do x=1,m(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          if(.not.covered(tree,p,[x,y,z])) cycle
          call reactive_conserved_to_primitive(species,candidate%fields(p)%values(:,x,y,z), &
            ct%fields(p)%values(1,x,y,z),q,t,cs,admitted, &
            internal_energy_density=held_candidate%fields(p)%values(1,x,y,z))
          if(.not.admitted) then
            local_ok=.false.;cycle
          endif
          target=pressure%fields(p)%values(1,x,y,z)
          if(abs(q(5)-target)<=64*epsilon(1._dp)*max(abs(q(5)),abs(target))) cycle
          call reconstruct_eb_pressure_state_3d(species,candidate%fields(p)%values(:,x,y,z), &
            held_candidate%fields(p)%values(1,x,y,z),ct%fields(p)%values(1,x,y,z),target,delta,admitted)
          local_ok=local_ok.and.admitted
        enddo
        enddo
        enddo
      enddo
      call agree(local_ok,comm,good)
    end subroutine prepare_clock_pressure

    recursive subroutine advance_level(level,time,interval,iteration,ncycle,flux,physical_flux,scalar_flux,good)
      integer, intent(in) :: level,iteration,ncycle
      real(dp), intent(in) :: time,interval
      type(eb_integrated_flux_3d), allocatable, intent(out) :: flux(:),physical_flux(:)
      type(eb_integrated_flux_3d), allocatable, intent(out) :: scalar_flux(:)
      logical, intent(out) :: good
      type(eb_held_funding_context_3d) :: funding_next,funding_clock_context
      type(amr_tree_state_3d), allocatable :: funding_stock,funding_receipt
      type(amr_tree_state_3d), allocatable :: funding_queue_next(:)
      type(amr_tree_state_3d) :: funding_payment,funding_debt
      type(eb_integrated_flux_3d), allocatable :: funding_faces(:),funding_scalar_faces(:)
      type(eb_integrated_flux_3d), allocatable :: fine_flux(:),fine_physical_flux(:)
      type(eb_integrated_flux_3d), allocatable :: fine_scalar_flux(:),step_scalar_flux(:)
      type(eb_integrated_flux_3d), allocatable :: sync_flux(:)
      type(eb_integrated_flux_3d), allocatable :: registered_replacement(:)
      type(eb_integrated_flux_3d), allocatable :: step_registration(:)
      type(amr_tree_state_3d) :: op,np,correction,redistributed,before_sync,before_t,registered_increment
      type(mpi_eb_mol_interval_3d), allocatable :: step_history,step_native_replay
      type(mpi_eb_mol_interval_3d) :: clock_record
      integer :: history_index
      type(amr_tree_state_3d) :: off_transfer
      real(dp) :: off_inventory(naux)
      real(dp), allocatable :: step_aux(:)
      real(dp), allocatable :: step_pressure_work
      type(eb_godunov_interval_3d), allocatable :: step_godunov,step_replay_godunov,step_source_recovery
      type(eb_godunov_interval_3d) :: godunov_clock
      integer :: godunov_index
      logical, allocatable :: source_prefix_ready
      logical :: partial_source_replay,complete_routed_replay
      type(eb_native_mol_source_prefix_3d), allocatable :: native_prefix_input,native_prefix_output
      logical, allocatable :: native_prefix_ready
      real(dp), allocatable :: common_native_weight
      logical :: native_retry_enabled,step_replaying_native
      integer :: native_attempt

      type(amr_tree_state_3d), allocatable :: transfer,before_pending(:)
      type(state_redist_register_3d), allocatable :: stage_registers(:),held_stage_registers(:)
      procedure(mpi_eb_sdc_auxiliary_source), pointer :: held_auxiliary_callback
      type(state_redist_register_3d) :: register_part
      real(dp) :: fstart(levels-1),fend(levels-1),b(nv),w(nv),ch(nv),ext(nv),limiter,elapsed,sub_dt
      integer :: p,d(3),step,allocation_status,target,bucket
      real(dp) :: clock_start_reset,clock_end_reset,scalar_reset_change
      logical :: accepted,local_ok,old_only,replayed_closure,paired_closure,step_refused
      good=.false.;replayed_closure=.false.;paired_closure=.false.;step_refused=.false.
      if(selected_pressure) then
        call prepare_clock_pressure(level,accepted)
        if(.not.accepted) return
      endif
      partial_source_replay=.false.;complete_routed_replay=.false.
      clock_start_reset=0;clock_end_reset=0;scalar_reset_change=0
      step_replaying_native=.false.
      if(replaying_native) then
        if(native_cursor>=size(native_replay_history)) then
          if(.not.native_prefix) return
        else
          step_replaying_native=.true.
        endif
      endif
      if(step_replaying_native) then
        native_cursor=native_cursor+1
        ! A prefix must follow the actual depth-first clock schedule. Reject
        ! a missing/reordered child before using its frozen source packets.
        associate(saved=>native_replay_history(native_cursor))
        call agree(saved%depth==level.and.saved%time==time.and.saved%interval==interval,comm,accepted)
        if(.not.accepted) return
        end associate
        allocate(step_native_replay,source=native_replay_history(native_cursor),stat=allocation_status)
        call agree(allocation_status==0,comm,accepted);if(.not.accepted) return
        do p=1,size(tree%nodes)
          if(owners(p)==rank.and.tree%nodes(p)%level==level) &
            step_native_replay%initial_reaction_source%fields(p)=ir%fields(p)
        enddo
      endif
      old_only=.false.
      if(present(mol_mode)) old_only=mol_mode.and..not.chemistry.and.iterations>1
      if(present(godunov_replay_history)) then
        if(godunov_replay_cursor>=godunov_replay_count) return
        godunov_replay_cursor=godunov_replay_cursor+1
        if(godunov_replay_cursor<=size(godunov_replay_history)) then
        allocate(step_replay_godunov,source=godunov_replay_history(godunov_replay_cursor),stat=godunov_status)
        call agree(godunov_status==0,comm,accepted);if(.not.accepted) return
        endif
      endif
      if(routing) then
        if(allocated(off_old%fields)) deallocate(off_old%fields)
        if(allocated(off_new%fields)) deallocate(off_new%fields)
      endif
      if(allocated(step_replay_godunov)) then
        partial_source_replay=step_replay_godunov%source_prefix_only
        complete_routed_replay=.not.partial_source_replay
        if(routing.and.allocated(step_replay_godunov%routed_auxiliary_old_rate%fields)) then
          call copy_godunov_field_3d(step_replay_godunov%routed_auxiliary_old_rate,off_old,local_ok)
      call agree(local_ok,comm,accepted);if(.not.accepted) return
        endif
        if(routing.and.allocated(step_replay_godunov%routed_auxiliary_new_rate%fields)) then
          call copy_godunov_field_3d(step_replay_godunov%routed_auxiliary_new_rate,off_new,local_ok)
      call agree(local_ok,comm,accepted);if(.not.accepted) return
        endif
      endif
      if(counts(level+1)==huge(1_int64)) return
      if(godunov_capture) then
        if(godunov_cursor>=size(godunov_records)) return
        godunov_cursor=godunov_cursor+1;godunov_index=godunov_cursor
        godunov_active(level)=godunov_index+godunov_replay_base
        allocate(step_godunov,step_source_recovery,source_prefix_ready,stat=godunov_status)
        call agree(godunov_status==0,comm,accepted);if(.not.accepted) return
        godunov_clock%parent_record=0
        if(level>0) godunov_clock%parent_record=godunov_active(level-1)
        call capture_godunov_clock(level,godunov_clock,accepted)
        if(.not.accepted) return
      endif
      if(allocated(funding_context_work)) then
        call copy_eb_held_funding_context_3d(tree,geometry,nv,funding_local,funding_context_work, &
          funding_clock_context,accepted)
        call agree(accepted,comm,local_ok);if(.not.local_ok) return
        if(collect_history) then
          clock_record%clock_initial=candidate;clock_record%clock_temperature=ct
          clock_record%clock_internal_energy=held_candidate;clock_record%clock_reaction_source=ir
          clock_record%clock_energy_correction=correction_work
          if(present(heat_release)) clock_record%clock_heat=heat
          clock_record%clock_source_funding_context=funding_clock_context
        endif
        call take_eb_held_funding_target_3d(tree,geometry,nv,funding_local,level, &
          funding_context_work%pending,funding_payment,funding_queue_next,accepted)
        call agree(accepted,comm,local_ok);if(.not.local_ok) return
        call apply_native_funding_receipt(funding_payment,accepted);if(.not.accepted) return
        call move_alloc(funding_queue_next,funding_context_work%pending)
      endif
      if(allocated(correction_work)) then
        call reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,candidate,ct,held_candidate,reset_difference, &
          held_policy,level,reset_change,comm,accepted,frame, &
          internal_energy_change=clock_start_reset,reset_level=merge(-2,-1,godunov_mol.and.chemistry.and.method=='none'))
        if(.not.accepted) return
        correction_work=correction_work+reset_change
      endif
      if(godunov_capture) then
        godunov_clock%held_clock_start_change=clock_start_reset
        if(godunov_index==1) godunov_clock%held_clock_start_change=clock_start_reset+initial_clock_reset
      endif
      clock_start(level+1)=time;clock_dt(level+1)=interval
      if(collect_history.and.present(internal_energy)) then
        allocate(clock_record%clock_start(levels),clock_record%clock_interval(levels),stat=allocation_status)
        call agree(allocation_status==0,comm,accepted);if(.not.accepted) return
        clock_record%clock_start=clock_start;clock_record%clock_interval=clock_dt
      endif
      fstart=0;fend=0
      fstart(:level)=max(0._dp,min(1._dp,(time-clock_start(:level))/clock_dt(:level)))
      fend(:level)=max(fstart(:level),min(1._dp,(time+interval-clock_start(:level))/clock_dt(:level)))
      do p=1,size(tree%nodes)
        if(owners(p)/=rank.or.tree%nodes(p)%level/=level) cycle
        starts%fields(p)%values=candidate%fields(p)%values
        ts%fields(p)%values=ct%fields(p)%values
        if(present(internal_energy)) held_starts%fields(p)=held_candidate%fields(p)
      end do
      if(method/='none'.and..not.prepared_redist) then
        allocate(transfer,stat=allocation_status)
        call agree(allocation_status==0,comm,accepted)
        if(.not.accepted) return
      end if
      if(has_auxiliary) allocate(step_aux(naux))
      if(collect_history) then
        allocate(step_history,stat=allocation_status)
        call agree(allocation_status==0,comm,accepted)
        if(.not.accepted) return
      endif
      if(collect_history.and.present(internal_energy)) then
        clock_record%held_clock_start_change=clock_start_reset
        if(history_count==0) then
          clock_record%held_clock_start_change=clock_start_reset+initial_clock_reset
          if(replaying_native) then
            call compute_reactive_eb_held_projection(0._dp,native_replay_history(1)%held_clock_start_change,0._dp, &
              clock_record%held_clock_start_change,scalar_reset_change,accepted)
            call agree(accepted,comm,local_ok);if(.not.local_ok) return
            clock_record%held_clock_start_change=scalar_reset_change
          endif
        endif
      endif
      if((replay_enabled.or.replaying_native).and.(prepared_redist.or.present(internal_energy))) then
        if(.not.allocated(funding_context_work)) then
        clock_record%clock_initial=candidate;clock_record%clock_temperature=ct
        clock_record%clock_reaction_source=ir
        if(present(heat_release)) clock_record%clock_heat=heat
        endif
        if(prepared_redist.and.level>0) clock_record%open_ancestor_registers=interface_registers(level-1)%fields
        if(present(internal_energy)) then
          if(.not.allocated(funding_context_work)) then
          clock_record%clock_internal_energy=held_candidate
          clock_record%clock_energy_correction=correction_work
          endif
          if(prepared_redist.and.level>0) &
            clock_record%open_ancestor_held_registers=held_interface_registers(level-1)%fields
        endif
      endif
      if(present(internal_energy_pressure_work_change).or.godunov_capture) then
        allocate(step_pressure_work,source=0._dp,stat=allocation_status)
        call agree(allocation_status==0,comm,accepted);if(.not.accepted) return
      endif
      if(allocated(funding_context_work)) then
        call read_eb_held_funding_open_faces_3d(tree,geometry,nv,funding_local,level, &
          funding_context_work,funding_faces,funding_scalar_faces,accepted)
        call agree(accepted,comm,local_ok);if(.not.local_ok) return
        call assemble_mpi_eb_held_funding_debt_3d(species,tree,owners,geometry,candidate,ct,held_candidate, &
          funding_faces,funding_scalar_faces,level,funding_debt,comm,accepted,boundary_faces,frame)
        if(.not.accepted) return
        allocate(funding_stock,funding_receipt,stat=funding_status)
        call agree(funding_status==0,comm,accepted);if(.not.accepted) return
        call read_eb_held_funding_stock_3d(tree,geometry,nv,funding_local,level, &
          funding_context_work%pending,funding_stock,accepted,open_debt=funding_debt)
        call agree(accepted,comm,local_ok);if(.not.local_ok) return
      endif
      source_time=time
      ! The gas predictor callback carries Nv rates only. UEINT prediction
      ! remains native and uses its independent scalar halos and registers.
      ! A level post callback follows paired gas/UEINT synchronization and
      ! owns only private external state until outer success.
      nullify(held_auxiliary_callback)
      if(has_auxiliary) held_auxiliary_callback=>source_adapter
      native_retry_enabled=(replay_enabled.or.(native_prefix.and..not.step_replaying_native )).and. &
        present(internal_energy).and.present(energy_correction).and. &
        method=='none'.and..not.routing.and..not.present(evaluate_spatial_diffusion).and. &
        .not.present(evaluate_diffusion).and..not.present(predictor_builder).and..not.present(post_level)
      if(native_retry_enabled) then
        allocate(native_prefix_output,native_prefix_ready)
        native_prefix_ready=.false.
      endif
      if(step_replaying_native) then
        if(present(native_spatial_weight)) then
          allocate(common_native_weight,source=native_spatial_weight)
        else if(preserve_native_weights) then
          allocate(common_native_weight,source=native_replay_history(native_cursor)%native_spatial_weight)
        endif
      endif
      do native_attempt=1,14
      if(present(internal_energy)) then
      if(prepared_redist) then
      if(has_auxiliary) then
      if(present(predictor_builder)) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,redistribution=method, &
        gas_acceleration=gas_acceleration,evaluate_reaction=evaluate_reaction, &
        predictor_builder=predictor_adapter, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
        evaluate_diffusion=evaluate_diffusion,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        internal_energy=held_candidate,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=held_starts,ancestor_internal_energy_end=held_ends, &
        target_volume_fraction=redist_target,state_redist_max_order=redist_order, &
        state_redist_registers=stage_registers,held_state_redist_registers=held_stage_registers, &
        interval_history=step_history,energy_policy=energy_policy,energy_correction=correction_work, &
        auxiliary_source=held_auxiliary_callback,auxiliary_change=step_aux, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        internal_energy_pressure_work_change=step_pressure_work, &
        godunov_interval_history=step_godunov,godunov_replay_interval=step_replay_godunov, &
        godunov_replay_partial=partial_source_replay,godunov_source_recovery=step_source_recovery, &
        godunov_source_recovery_ready=source_prefix_ready, &
        spatial_flux_weight=spatial_flux_weight,provider_refused=step_refused, &
        source_funding_stock=funding_stock,source_funding_impulse=funding_receipt, &
        replay_interval=step_native_replay,freeze_native_external_sources=step_replaying_native, &
        native_spatial_weight=common_native_weight, &
        native_source_prefix=native_prefix_output,native_source_prefix_ready=native_prefix_ready, &
        native_source_replay=native_prefix_input, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,redistribution=method, &
        gas_acceleration=gas_acceleration,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
        evaluate_diffusion=evaluate_diffusion,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        internal_energy=held_candidate,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=held_starts,ancestor_internal_energy_end=held_ends, &
        target_volume_fraction=redist_target,state_redist_max_order=redist_order, &
        state_redist_registers=stage_registers,held_state_redist_registers=held_stage_registers, &
        interval_history=step_history,energy_policy=energy_policy,energy_correction=correction_work, &
        auxiliary_source=held_auxiliary_callback,auxiliary_change=step_aux, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        internal_energy_pressure_work_change=step_pressure_work, &
        godunov_interval_history=step_godunov,godunov_replay_interval=step_replay_godunov, &
        godunov_replay_partial=partial_source_replay,godunov_source_recovery=step_source_recovery, &
        godunov_source_recovery_ready=source_prefix_ready, &
        spatial_flux_weight=spatial_flux_weight,provider_refused=step_refused, &
        source_funding_stock=funding_stock,source_funding_impulse=funding_receipt, &
        replay_interval=step_native_replay,freeze_native_external_sources=step_replaying_native, &
        native_spatial_weight=common_native_weight, &
        native_source_prefix=native_prefix_output,native_source_prefix_ready=native_prefix_ready, &
        native_source_replay=native_prefix_input, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      endif
      else
      if(present(predictor_builder)) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,redistribution=method, &
        gas_acceleration=gas_acceleration,evaluate_reaction=evaluate_reaction, &
        predictor_builder=predictor_adapter, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
        evaluate_diffusion=evaluate_diffusion,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        internal_energy=held_candidate,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=held_starts,ancestor_internal_energy_end=held_ends, &
        target_volume_fraction=redist_target,state_redist_max_order=redist_order, &
        state_redist_registers=stage_registers,held_state_redist_registers=held_stage_registers, &
        interval_history=step_history,energy_policy=energy_policy,energy_correction=correction_work, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        internal_energy_pressure_work_change=step_pressure_work, &
        godunov_interval_history=step_godunov,godunov_replay_interval=step_replay_godunov, &
        godunov_replay_partial=partial_source_replay,godunov_source_recovery=step_source_recovery, &
        godunov_source_recovery_ready=source_prefix_ready, &
        spatial_flux_weight=spatial_flux_weight,provider_refused=step_refused, &
        source_funding_stock=funding_stock,source_funding_impulse=funding_receipt, &
        replay_interval=step_native_replay,freeze_native_external_sources=step_replaying_native, &
        native_spatial_weight=common_native_weight, &
        native_source_prefix=native_prefix_output,native_source_prefix_ready=native_prefix_ready, &
        native_source_replay=native_prefix_input, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,redistribution=method, &
        gas_acceleration=gas_acceleration,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
        evaluate_diffusion=evaluate_diffusion,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        internal_energy=held_candidate,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=held_starts,ancestor_internal_energy_end=held_ends, &
        target_volume_fraction=redist_target,state_redist_max_order=redist_order, &
        state_redist_registers=stage_registers,held_state_redist_registers=held_stage_registers, &
        interval_history=step_history,energy_policy=energy_policy,energy_correction=correction_work, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        internal_energy_pressure_work_change=step_pressure_work, &
        godunov_interval_history=step_godunov,godunov_replay_interval=step_replay_godunov, &
        godunov_replay_partial=partial_source_replay,godunov_source_recovery=step_source_recovery, &
        godunov_source_recovery_ready=source_prefix_ready, &
        spatial_flux_weight=spatial_flux_weight,provider_refused=step_refused, &
        source_funding_stock=funding_stock,source_funding_impulse=funding_receipt, &
        replay_interval=step_native_replay,freeze_native_external_sources=step_replaying_native, &
        native_spatial_weight=common_native_weight, &
        native_source_prefix=native_prefix_output,native_source_prefix_ready=native_prefix_ready, &
        native_source_replay=native_prefix_input, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      endif
      endif
      else
      if(has_auxiliary) then
      if(present(predictor_builder)) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,redistribution=method, &
        gas_acceleration=gas_acceleration,evaluate_reaction=evaluate_reaction, &
        predictor_builder=predictor_adapter, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
        evaluate_diffusion=evaluate_diffusion,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        internal_energy=held_candidate,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=held_starts,ancestor_internal_energy_end=held_ends,interval_history=step_history, &
        energy_policy=energy_policy,energy_correction=correction_work, &
        auxiliary_source=source_adapter,auxiliary_change=step_aux, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        internal_energy_pressure_work_change=step_pressure_work, &
        godunov_interval_history=step_godunov,godunov_replay_interval=step_replay_godunov, &
        godunov_replay_partial=partial_source_replay,godunov_source_recovery=step_source_recovery, &
        godunov_source_recovery_ready=source_prefix_ready, &
        spatial_flux_weight=spatial_flux_weight,provider_refused=step_refused, &
        source_funding_stock=funding_stock,source_funding_impulse=funding_receipt, &
        replay_interval=step_native_replay,freeze_native_external_sources=step_replaying_native, &
        native_spatial_weight=common_native_weight, &
        native_source_prefix=native_prefix_output,native_source_prefix_ready=native_prefix_ready, &
        native_source_replay=native_prefix_input, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,redistribution=method, &
        gas_acceleration=gas_acceleration,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
        evaluate_diffusion=evaluate_diffusion,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        internal_energy=held_candidate,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=held_starts,ancestor_internal_energy_end=held_ends,interval_history=step_history, &
        energy_policy=energy_policy,energy_correction=correction_work, &
        auxiliary_source=source_adapter,auxiliary_change=step_aux, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        internal_energy_pressure_work_change=step_pressure_work, &
        godunov_interval_history=step_godunov,godunov_replay_interval=step_replay_godunov, &
        godunov_replay_partial=partial_source_replay,godunov_source_recovery=step_source_recovery, &
        godunov_source_recovery_ready=source_prefix_ready, &
        spatial_flux_weight=spatial_flux_weight,provider_refused=step_refused, &
        source_funding_stock=funding_stock,source_funding_impulse=funding_receipt, &
        replay_interval=step_native_replay,freeze_native_external_sources=step_replaying_native, &
        native_spatial_weight=common_native_weight, &
        native_source_prefix=native_prefix_output,native_source_prefix_ready=native_prefix_ready, &
        native_source_replay=native_prefix_input, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      endif
      else
      if(present(predictor_builder)) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,redistribution=method, &
        gas_acceleration=gas_acceleration,evaluate_reaction=evaluate_reaction, &
        predictor_builder=predictor_adapter, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
        evaluate_diffusion=evaluate_diffusion,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        internal_energy=held_candidate,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=held_starts,ancestor_internal_energy_end=held_ends,interval_history=step_history, &
        energy_policy=energy_policy,energy_correction=correction_work, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        internal_energy_pressure_work_change=step_pressure_work, &
        godunov_interval_history=step_godunov,godunov_replay_interval=step_replay_godunov, &
        godunov_replay_partial=partial_source_replay,godunov_source_recovery=step_source_recovery, &
        godunov_source_recovery_ready=source_prefix_ready, &
        spatial_flux_weight=spatial_flux_weight,provider_refused=step_refused, &
        source_funding_stock=funding_stock,source_funding_impulse=funding_receipt, &
        replay_interval=step_native_replay,freeze_native_external_sources=step_replaying_native, &
        native_spatial_weight=common_native_weight, &
        native_source_prefix=native_prefix_output,native_source_prefix_ready=native_prefix_ready, &
        native_source_replay=native_prefix_input, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,redistribution=method, &
        gas_acceleration=gas_acceleration,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
        evaluate_diffusion=evaluate_diffusion,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        internal_energy=held_candidate,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=held_starts,ancestor_internal_energy_end=held_ends,interval_history=step_history, &
        energy_policy=energy_policy,energy_correction=correction_work, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        internal_energy_pressure_work_change=step_pressure_work, &
        godunov_interval_history=step_godunov,godunov_replay_interval=step_replay_godunov, &
        godunov_replay_partial=partial_source_replay,godunov_source_recovery=step_source_recovery, &
        godunov_source_recovery_ready=source_prefix_ready, &
        spatial_flux_weight=spatial_flux_weight,provider_refused=step_refused, &
        source_funding_stock=funding_stock,source_funding_impulse=funding_receipt, &
        replay_interval=step_native_replay,freeze_native_external_sources=step_replaying_native, &
        native_spatial_weight=common_native_weight, &
        native_source_prefix=native_prefix_output,native_source_prefix_ready=native_prefix_ready, &
        native_source_replay=native_prefix_input, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      endif
      endif
      endif
      else if(prepared_redist) then
      if(present(predictor_builder)) then
      if(has_auxiliary) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,redistribution=method, &
        target_volume_fraction=redist_target,auxiliary_source=source_adapter,auxiliary_change=step_aux, &
        gas_acceleration=gas_acceleration,predictor_builder=predictor_adapter, &
        evaluate_diffusion=evaluate_diffusion, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=redist_order,state_redist_registers=stage_registers, &
      evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        interval_history=step_history, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,redistribution=method, &
        target_volume_fraction=redist_target,gas_acceleration=gas_acceleration,predictor_builder=predictor_adapter, &
        evaluate_diffusion=evaluate_diffusion, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=redist_order,state_redist_registers=stage_registers, &
      evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        interval_history=step_history, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      end if
      else
      if(has_auxiliary) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,redistribution=method, &
        target_volume_fraction=redist_target,auxiliary_source=source_adapter,auxiliary_change=step_aux, &
        gas_acceleration=gas_acceleration,evaluate_diffusion=evaluate_diffusion, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=redist_order,state_redist_registers=stage_registers, &
      evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        interval_history=step_history, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,redistribution=method, &
        target_volume_fraction=redist_target,gas_acceleration=gas_acceleration,evaluate_diffusion=evaluate_diffusion, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=redist_order,state_redist_registers=stage_registers, &
      evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        interval_history=step_history, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      end if
      end if
      else
      if(present(predictor_builder)) then
      if(has_auxiliary) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,method,transfer, &
        redist_target,source_adapter,step_aux,gas_acceleration=gas_acceleration,predictor_builder=predictor_adapter, &
        evaluate_diffusion=evaluate_diffusion, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=redist_order,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        interval_history=step_history, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,method,transfer, &
        redist_target,gas_acceleration=gas_acceleration,predictor_builder=predictor_adapter, &
        evaluate_diffusion=evaluate_diffusion, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=redist_order,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        interval_history=step_history, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      end if
      else
      if(has_auxiliary) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,method,transfer, &
        redist_target,source_adapter,step_aux,gas_acceleration=gas_acceleration,evaluate_diffusion=evaluate_diffusion, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=redist_order,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        interval_history=step_history, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,candidate,ct, &
        extended_geometry,interval,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces, &
        b,w,ch,ext,limiter,comm,accepted,ir,difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        wall_temperature,level,starts,ends,ts,te,fstart,fend,source0,source1,flux,op,np,method,transfer, &
        redist_target,gas_acceleration=gas_acceleration,evaluate_diffusion=evaluate_diffusion, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=redist_order,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        interval_history=step_history, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        evaluate_reaction_admission=evaluate_reaction_admission,reaction_admission_status=child_admission)
      admission_state=max(admission_state,child_admission)
      if(child_admission==reaction_admission_provider_refused) then
        godunov_provider_rejected=.true.;prepared_provider_rejected=.true.
      endif
      end if
      end if
      endif
      if(accepted.or..not.native_retry_enabled) exit
      if(.not.native_prefix_ready.or.step_refused.or.native_attempt==14) exit
      ! The common transaction has preserved candidate/T/H/I_R and receipts.
      ! Extend the actual source cache; never request a known packet again.
      if(allocated(native_prefix_input)) deallocate(native_prefix_input)
      allocate(native_prefix_input,source=native_prefix_output)
      if(.not.allocated(common_native_weight)) allocate(common_native_weight)
      common_native_weight=2._dp**(-native_attempt)
      if(native_attempt==13) common_native_weight=0
      if(present(native_common_retries)) native_common_retries=native_common_retries+1_int64
      enddo
      godunov_provider_rejected=godunov_provider_rejected.or.step_refused
      if(present(provider_refused)) provider_refused=godunov_provider_rejected
      if(.not.accepted) then
        if(godunov_capture) then
          if(source_prefix_ready) then
            call store_godunov_step(step_source_recovery,godunov_clock,godunov_index,time,accepted)
            if(.not.accepted) return
            if(godunov_automatic_repair) go to 880
          endif
        endif
        return
      endif
      if(godunov_capture) then
        call store_godunov_step(step_godunov,godunov_clock,godunov_index,time,accepted)
        if(.not.accepted) return
      endif
      if(allocated(step_pressure_work)) then
        pressure_work_change=pressure_work_change+step_pressure_work
        call agree(ieee_is_finite(pressure_work_change),comm,accepted)
        if(.not.accepted) return
      endif
      if(allocated(funding_context_work)) then
        call queue_eb_held_funding_receipt_3d(tree,geometry,nv,funding_local,level, &
          funding_context_work%pending,funding_receipt,funding_queue_next,accepted)
        call agree(accepted,comm,local_ok);if(.not.local_ok) return
        call move_alloc(funding_queue_next,funding_context_work%pending)
      endif
      if(collect_history) then
        history_count=history_count+1;history_index=history_count
        step_history%time=time
        call move_eb_held_funding_context_3d(clock_record%clock_source_funding_context, &
          step_history%clock_source_funding_context)
        call move_alloc(clock_record%clock_initial%fields,step_history%clock_initial%fields)
        call move_alloc(clock_record%clock_temperature%fields,step_history%clock_temperature%fields)
        call move_alloc(clock_record%clock_reaction_source%fields,step_history%clock_reaction_source%fields)
        call move_alloc(clock_record%clock_heat%fields,step_history%clock_heat%fields)
        call move_alloc(clock_record%clock_internal_energy%fields,step_history%clock_internal_energy%fields)
        step_history%clock_energy_correction=clock_record%clock_energy_correction
        call move_alloc(clock_record%clock_start,step_history%clock_start)
        call move_alloc(clock_record%clock_interval,step_history%clock_interval)
        call compute_reactive_eb_held_projection(0._dp,clock_record%held_clock_start_change,0._dp, &
          step_history%initial_internal_energy_change,step_history%held_clock_start_change,accepted)
        if(.not.accepted) return
        call move_alloc(clock_record%open_ancestor_registers,step_history%open_ancestor_registers)
        call move_alloc(clock_record%open_ancestor_held_registers,step_history%open_ancestor_held_registers)
        step_history%old_predictor=op;step_history%new_predictor=np
        call move_mol_interval(step_history,history(history_index))
      endif
      local_ok=.true.
      if(present(registered_flux)) then
        if(.not.allocated(registration)) then
          allocate(registration(size(tree%nodes)),stat=allocation_status)
          call agree(allocation_status==0,comm,accepted)
          if(.not.accepted) return
          local_ok=.true.
          do p=1,size(registration)
            if(owners(p)/=rank) cycle
            d=tree%nodes(p)%shape
            allocate(registration(p)%x(nv,0:d(1),d(2),d(3)),registration(p)%y(nv,d(1),0:d(2),d(3)), &
              registration(p)%z(nv,d(1),d(2),0:d(3)),source=0._dp,stat=allocation_status)
            if(allocation_status/=0) then
              local_ok=.false.;exit
            endif
          enddo
          call agree(local_ok,comm,accepted)
          if(.not.accepted) return
        endif
        do p=1,size(registration)
          if(owners(p)/=rank) cycle
          if(tree%nodes(p)%level/=level) then
            if(.not.allocated(step_registration(p)%x).and..not.allocated(step_registration(p)%y).and. &
              .not.allocated(step_registration(p)%z)) cycle
          endif
          if(.not.valid_eb_integrated_flux_3d(step_registration(p),geometry(p),nv)) then
            local_ok=.false.;exit
          endif
          registration(p)%x=registration(p)%x+step_registration(p)%x
          registration(p)%y=registration(p)%y+step_registration(p)%y
          registration(p)%z=registration(p)%z+step_registration(p)%z
          local_ok=local_ok.and.all(ieee_is_finite(registration(p)%x)).and. &
            all(ieee_is_finite(registration(p)%y)).and.all(ieee_is_finite(registration(p)%z))
        enddo
      endif
      call agree(local_ok,comm,accepted)
      if(.not.accepted) return
      if(prepared_redist) then
        local_ok=.true.
        do p=1,size(tree%nodes)
          if(owners(p)/=rank.or.tree%nodes(p)%level/=level) cycle
          d=tree%nodes(p)%shape
          if(level+1<levels) then
            register_part=stage_registers(p)
            if(allocated(register_part%fine_delta)) register_part%fine_delta=0
            call accumulate_state_redist_register_3d(interface_registers(level)%fields(p),register_part,nv,d+12,accepted)
            if(.not.accepted) then
              local_ok=.false.;exit
            endif
            if(present(internal_energy)) then
              register_part=held_stage_registers(p)
              if(allocated(register_part%fine_delta)) register_part%fine_delta=0
              call accumulate_state_redist_register_3d(held_interface_registers(level)%fields(p),register_part,1,d+12,accepted)
              if(.not.accepted) then
                local_ok=.false.;exit
              endif
            endif
          endif
          if(level>0) then
            register_part=stage_registers(p)
            if(allocated(register_part%coarse_delta)) register_part%coarse_delta=0
            call accumulate_state_redist_register_3d(interface_registers(level-1)%fields(p),register_part,nv,d+12,accepted)
            if(.not.accepted) then
              local_ok=.false.;exit
            endif
            if(present(internal_energy)) then
              register_part=held_stage_registers(p)
              if(allocated(register_part%coarse_delta)) register_part%coarse_delta=0
              call accumulate_state_redist_register_3d(held_interface_registers(level-1)%fields(p), &
                register_part,1,d+12,accepted)
              if(.not.accepted) then
                local_ok=.false.;exit
              endif
            endif
          endif
        enddo
        call agree(local_ok,comm,accepted)
        if(.not.accepted) return
      elseif(method/='none') then
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          target=tree%nodes(p)%level;bucket=min(level,target)
          pending(bucket)%fields(p)%values=pending(bucket)%fields(p)%values+transfer%fields(p)%values
        end do
      end if
      if(routing) then
        if(complete_routed_replay) then
          call copy_godunov_field_3d(step_replay_godunov%routed_auxiliary_impulse,off_transfer,local_ok)
          call agree(local_ok,comm,accepted);if(.not.accepted) return
        else
          off_transfer=off_old
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            if(old_auxiliary_only) then
              off_transfer%fields(p)%values=.5_dp*interval*off_old%fields(p)%values
            else
              off_transfer%fields(p)%values=.5_dp*interval*(off_old%fields(p)%values+off_new%fields(p)%values)
            endif
          enddo
        endif
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          target=tree%nodes(p)%level;bucket=min(level,target)
          auxiliary_pending(bucket)%fields(p)%values=auxiliary_pending(bucket)%fields(p)%values+ &
            off_transfer%fields(p)%values
        enddo
        if(godunov_capture) then
          call copy_godunov_field_3d(off_transfer,godunov_records(godunov_index)%routed_auxiliary_impulse,local_ok)
          call agree(local_ok,comm,accepted);if(.not.accepted) return
        endif
        if(collect_history) then
          history(history_index)%routed_source=off_transfer
          history(history_index)%routed_auxiliary=.true.
        endif
        call integrate_mpi_eb_tree_3d(tree,owners,geometry,off_transfer,off_inventory,comm,accepted)
        if(.not.accepted) return
        step_aux=step_aux+off_inventory
        if(godunov_capture) godunov_records(godunov_index)%auxiliary_change=step_aux
      end if
      if(has_auxiliary) then
        auxiliary=auxiliary+step_aux
        call agree(all(ieee_is_finite(auxiliary)),comm,accepted)
        if(.not.accepted) return
      end if
      counts(level+1)=counts(level+1)+1
      wall=wall+w;chemical=chemical+ch;external=external+ext;theta=min(theta,limiter)
      local_ok=.true.
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        if(tree%nodes(p)%level==level) then
          ends%fields(p)%values=candidate%fields(p)%values
          te%fields(p)%values=ct%fields(p)%values
          if(present(internal_energy)) held_ends%fields(p)=held_candidate%fields(p)
          source0%fields(p)%values=op%fields(p)%values;source1%fields(p)%values=np%fields(p)%values
        else
          ! Owner-local inactive clocks retain accepted donor face increments.
          call complete_interval_faces(flux(p),geometry(p),nv,local_ok)
          if(.not.local_ok) exit
          if(present(internal_energy)) then
            call complete_interval_faces(step_scalar_flux(p),geometry(p),1,local_ok)
            if(.not.local_ok) exit
          endif
        end if
      end do
      call agree(local_ok,comm,accepted)
      if(.not.accepted) return
      physical_flux=flux
      if(present(internal_energy)) scalar_flux=step_scalar_flux
      if(present(mol_mode)) then
        if(mol_mode.and..not.chemistry.and.iterations>1) then
          do p=1,size(flux)
            if(owners(p)/=rank) cycle
            if(tree%nodes(p)%level/=level) cycle
            flux(p)=step_registration(p)
          enddo
        endif
      endif
      if(allocated(funding_context_work)) then
        call set_eb_held_funding_clock_faces_3d(tree,geometry,nv,funding_local,level, &
          physical_flux,scalar_flux,funding_context_work,funding_next,accepted)
        call agree(accepted,comm,local_ok);if(.not.local_ok) return
        call move_eb_held_funding_context_3d(funding_next,funding_context_work)
      endif
      if(level+1<levels) then
        elapsed=0
        do step=1,ratios(level+2)
          sub_dt=interval/real(ratios(level+2),dp)
          if(step==ratios(level+2)) sub_dt=interval-elapsed
          if(sub_dt<=0.or.time+elapsed+sub_dt<=time+elapsed) return
          call advance_level(level+1,time+elapsed,sub_dt,step,ratios(level+2),fine_flux, &
            fine_physical_flux,fine_scalar_flux,accepted)
          if(.not.accepted) then
            if(godunov_automatic_repair) go to 880
            ! Every local interval must already have been captured, and this
            ! must be the last child in the enclosing clock window. Never
            ! retry a source or post callback that has reported failure.
            if(.not.replay_enabled.or. &
              (method/='flux_redist'.and..not.prepared_redist.and..not.present(internal_energy))) return
            if(step/=ratios(level+2).or.deferred_record(level+1)==0) return
            call replay_parent_interval(history_index,history_count,level,flux,accepted,scalar_flux)
            if(.not.accepted) then
              if(prepared_provider_rejected) return
              deferred_record(level)=history_index
              deferred_iteration(level)=iteration;deferred_cycle(level)=ncycle
              return
            endif
            call finish_deferred_closures(level,accepted)
            if(.not.accepted) return
            replayed_closure=.true.
            exit
          endif
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            ! Descendant clocks can also fund an exposed ancestor donor.
            ! Its returned faces are increments to this open native family.
            if(tree%nodes(p)%level<=level) then
              if(all(fine_flux(p)%x==0).and.all(fine_flux(p)%y==0).and.all(fine_flux(p)%z==0).and. &
                all(fine_physical_flux(p)%x==0).and.all(fine_physical_flux(p)%y==0).and. &
                all(fine_physical_flux(p)%z==0)) then
                if(.not.present(internal_energy)) cycle
                if(all(fine_scalar_flux(p)%x==0).and.all(fine_scalar_flux(p)%y==0).and. &
                  all(fine_scalar_flux(p)%z==0)) cycle
              endif
            endif
            flux(p)%x=flux(p)%x+fine_flux(p)%x
            flux(p)%y=flux(p)%y+fine_flux(p)%y
            flux(p)%z=flux(p)%z+fine_flux(p)%z
            physical_flux(p)%x=physical_flux(p)%x+fine_physical_flux(p)%x
            physical_flux(p)%y=physical_flux(p)%y+fine_physical_flux(p)%y
            physical_flux(p)%z=physical_flux(p)%z+fine_physical_flux(p)%z
            if(present(internal_energy)) then
              scalar_flux(p)%x=scalar_flux(p)%x+fine_scalar_flux(p)%x
              scalar_flux(p)%y=scalar_flux(p)%y+fine_scalar_flux(p)%y
              scalar_flux(p)%z=scalar_flux(p)%z+fine_scalar_flux(p)%z
            endif
          end do
          elapsed=elapsed+sub_dt
        end do
      end if
      if(allocated(funding_context_work)) then
        call take_eb_held_funding_bucket_3d(tree,geometry,nv,funding_local,level, &
          funding_context_work%pending,funding_payment,funding_queue_next,accepted)
        call agree(accepted,comm,local_ok);if(.not.local_ok) return
        call apply_native_funding_receipt(funding_payment,accepted);if(.not.accepted) return
        call move_alloc(funding_queue_next,funding_context_work%pending)
      endif
      if(replayed_closure) then
        if(allocated(funding_context_work)) paired_closure=.true.
        if(.not.prepared_redist.and.method=='flux_redist') then
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          candidate%fields(p)%values=candidate%fields(p)%values+pending(level)%fields(p)%values
          pending(level)%fields(p)%values=0
        enddo
        call recover_transfer(accepted)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        endif
        if(prepared_redist) then
          if(level+1<levels) then
            call close_prepared_level(candidate,ct,flux,scalar_flux,level,accepted,paired_closure, &
              scalar_change=clock_end_reset)
          else
            call close_prepared_level(candidate,ct,flux,scalar_flux,level,accepted,paired_closure, &
              scalar_change=clock_end_reset)
          endif
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted) return
        endif
        if(present(internal_energy)) then
          call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,accepted,frame,held_candidate)
        else
          call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,accepted,frame)
        endif
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
      else
      if(prepared_redist) then
        if(old_only) then
          ! Match serial: conserve the accepted two-half-step gas transfer
          ! while retaining the pinned old-half registration for the parent.
          if(present(internal_energy).and.allocated(correction_work)) then
            call close_prepared_level(candidate,ct,physical_flux,scalar_flux,level,accepted,paired_closure, &
              scalar_change=clock_end_reset)
          else if(level+1<levels) then
            call synchronize_mpi_eb_tree_register_flux_3d(species,tree,owners,geometry,register_geometry, &
              candidate,ct,physical_flux,interface_registers(level)%fields,level,comm,accepted,boundary_faces,frame)
          else
            call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,physical_flux,comm,accepted, &
              minimum_level=level,maximum_level=level,boundary_faces=boundary_faces,frame=frame)
          endif
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted) return
          call assemble_mpi_eb_tree_reflux_3d(species,tree,owners,geometry,candidate,ct,flux,registered_increment, &
            registered_replacement,comm,accepted,level,min(level+1,levels-1),boundary_faces,frame)
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted) return
          call move_alloc(registered_replacement,flux)
        else
          if(level+1<levels) then
            call close_prepared_level(candidate,ct,flux,scalar_flux,level,accepted,paired_closure, &
              scalar_change=clock_end_reset)
          else
            call close_prepared_level(candidate,ct,flux,scalar_flux,level,accepted,paired_closure, &
              scalar_change=clock_end_reset)
          endif
        endif
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted.and.replay_enabled.and.level+1<levels) then
          call replay_parent_interval(history_index,history_count,level,flux,accepted,scalar_flux)
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted) then
            if(prepared_provider_rejected) return
            deferred_record(level)=history_index
            deferred_iteration(level)=iteration;deferred_cycle(level)=ncycle
            return
          endif
          call finish_deferred_closures(level,accepted)
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted) return
          call close_prepared_level(candidate,ct,flux,scalar_flux,level,accepted,paired_closure, &
              scalar_change=clock_end_reset)
        endif
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
      else
      if(method/='none') then
        before_sync=candidate;before_t=ct
        if(replay_enabled.and.method=='flux_redist') before_pending=pending
      end if
      if(old_only) then
        call move_alloc(physical_flux,sync_flux)
      else
        call move_alloc(flux,sync_flux)
      endif
      if(present(internal_energy).and.allocated(correction_work).and.method=='none') then
        call close_held_prepared_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,held_candidate, &
          sync_flux,scalar_flux,held_policy,correction_work,level,comm,accepted,boundary_faces,frame, &
          register_geometry=extended_geometry,internal_energy_change=clock_end_reset, &
          reset_current_level=godunov_mol.and.chemistry)
        paired_closure=accepted
      else
      call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,sync_flux,comm,accepted, &
        minimum_level=level,maximum_level=min(level+1,levels-1),boundary_faces=boundary_faces,frame=frame)
      endif
      if(.not.accepted.and.godunov_automatic_repair) go to 880
      if(.not.accepted.and.level+1<levels) then
        ! Gas-only interface repair would leave the held family behind. Held
        ! chemistry instead retries the paired native interval transaction.
        if(.not.present(internal_energy)) then
          call repair_mpi_eb_tree_interface_budget_3d(species,tree,owners,geometry,candidate,ct,sync_flux,comm,accepted, &
            level,boundary_faces,frame)
        endif
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted.and.replay_enabled) then
          call replay_parent_interval(history_index,history_count,level,sync_flux,accepted,scalar_flux)
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted.and.prepared_provider_rejected) return
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted.and.(method=='flux_redist'.or.present(internal_energy))) then
            deferred_record(level)=history_index
            deferred_iteration(level)=iteration;deferred_cycle(level)=ncycle
          endif
          if(accepted.and.allocated(funding_context_work)) then
            replayed_closure=.true.;paired_closure=.true.
          endif
        endif
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        if(method/='none') then
          before_sync=candidate;before_t=ct
        endif
        if(replayed_closure.and.allocated(funding_context_work)) then
          ! The accepted native transaction has closed gas, scalar and queue.
          paired_closure=.true.
        else if(present(internal_energy).and.allocated(correction_work).and.method=='none') then
          call close_held_prepared_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,held_candidate, &
            sync_flux,scalar_flux,held_policy,correction_work,level,comm,accepted,boundary_faces,frame, &
            register_geometry=extended_geometry,internal_energy_change=clock_end_reset, &
            reset_current_level=godunov_mol.and.chemistry)
          paired_closure=accepted
        else
        call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,sync_flux,comm,accepted, &
          minimum_level=level,maximum_level=level+1,boundary_faces=boundary_faces,frame=frame)
        endif
      endif
      if(.not.accepted.and.godunov_automatic_repair) go to 880
      if(.not.accepted) return
      if(old_only) then
        call move_alloc(sync_flux,physical_flux)
      else
        call move_alloc(sync_flux,flux)
      endif
      if(method/='none') then
        correction=source0
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          correction%fields(p)%values=0
          target=tree%nodes(p)%level
          if(target<level.or.target>level+1) cycle
          correction%fields(p)%values=candidate%fields(p)%values-before_sync%fields(p)%values
        end do
        if(method=='state_redist') then
          call redistribute_mpi_eb_tree_state_3d(tree,owners,geometry,correction,nv,redistributed, &
            comm,accepted,redist_target,boundary_faces, &
      max_order=redist_order)
        else
          call redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree,owners,geometry,correction,nv,redistributed, &
            comm,accepted,boundary_faces)
        end if
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        if(collect_history) then
          history(history_index)%synchronization_transfer=source0
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            history(history_index)%synchronization_transfer%fields(p)%values=0
            if(tree%nodes(p)%level>=level) cycle
            history(history_index)%synchronization_transfer%fields(p)%values=redistributed%fields(p)%values
          enddo
        endif
        candidate=before_sync;ct=before_t
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          target=tree%nodes(p)%level
          if(target<level) then
            pending(target)%fields(p)%values=pending(target)%fields(p)%values+redistributed%fields(p)%values
          else
            candidate%fields(p)%values=candidate%fields(p)%values+redistributed%fields(p)%values+ &
              pending(level)%fields(p)%values
          end if
          pending(level)%fields(p)%values=0
        end do
        call recover_transfer(accepted)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted.and.replay_enabled.and.level+1<levels) then
          candidate=before_sync;ct=before_t
          pending=before_pending
          call replay_parent_interval(history_index,history_count,level,flux,accepted,scalar_flux)
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted) then
            if(prepared_provider_rejected) return
            deferred_record(level)=history_index
            deferred_iteration(level)=iteration;deferred_cycle(level)=ncycle
            return
          endif
          do p=1,size(tree%nodes)
            if(owners(p)/=rank) cycle
            candidate%fields(p)%values=candidate%fields(p)%values+pending(level)%fields(p)%values
            pending(level)%fields(p)%values=0
          enddo
          call recover_transfer(accepted)
        endif
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,accepted,frame)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
      end if
      endif
      endif ! normal synchronization or already closed ancestor replay
      if(present(internal_energy).and..not.paired_closure) then
        if(prepared_redist.and.level+1<levels) then
          call synchronize_mpi_eb_tree_scalar_register_flux_3d(species,tree,owners,geometry,register_geometry, &
            held_candidate,scalar_flux,held_interface_registers(level)%fields,level,comm,accepted,boundary_faces)
        else
          call synchronize_mpi_eb_tree_scalar_flux_3d(species,tree,owners,geometry,held_candidate,scalar_flux, &
            comm,accepted,minimum_level=level,maximum_level=min(level+1,levels-1),boundary_faces=boundary_faces)
        endif
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
      endif
      if(allocated(correction_work).and..not.paired_closure) then
        call reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,candidate,ct,held_candidate,reset_difference, &
          held_policy,level,reset_change,comm,accepted,frame,leaf_only=.true., &
          internal_energy_change=scalar_reset_change)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        clock_end_reset=clock_end_reset+scalar_reset_change
        correction_work=correction_work+reset_change
      endif
      if(routing) then
        if(held_auxiliary) then
          call apply_routed_held_source(candidate,ct,held_candidate,correction_work,auxiliary_pending(level),accepted, &
            scalar_change=scalar_reset_change)
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted) return
          clock_end_reset=clock_end_reset+scalar_reset_change
          do p=1,size(tree%nodes)
            if(owners(p)==rank) auxiliary_pending(level)%fields(p)%values=0
          enddo
        else
        do p=1,size(tree%nodes)
          if(owners(p)/=rank) cycle
          candidate%fields(p)%values=candidate%fields(p)%values+auxiliary_pending(level)%fields(p)%values
          auxiliary_pending(level)%fields(p)%values=0
        end do
        call recover_transfer(accepted)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,accepted,frame)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        endif
      endif
      if(chemistry.and.level+1<levels) then
        call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,ir,nv-4,comm,accepted,minimum_depth=level)
        if(.not.accepted) return
      endif
      if(present(heat_release).and.level+1<levels) then
        call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,heat,1,comm,accepted,minimum_depth=level)
        if(.not.accepted) return
      endif
      if(present(post_level)) then
        call post_level(level,time,interval,iteration,ncycle,comm,local_ok)
        call agree(local_ok,comm,accepted)
        if(.not.accepted) return
      end if
      if(present(mol_mode)) then
        if(.not.(mol_mode.and..not.chemistry.and.iterations>1)) physical_flux=flux
      else
        physical_flux=flux
      endif
      if(collect_history) then
        call copy_mol_level_flux(tree,owners,rank,level,physical_flux, &
          history(history_index)%synchronized_flux,comm,accepted)
        if(.not.accepted) return
      endif
      if(collect_history.and.present(internal_energy)) then
        if(replayed_closure.and.allocated(funding_context_work)) &
          clock_end_reset=history(history_index)%held_clock_end_change
        history(history_index)%held_clock_end_change=clock_end_reset
        history(history_index)%held_clock_events_closed=.true.
      endif
      if(godunov_capture) then
        godunov_records(godunov_index)%held_clock_end_change=clock_end_reset
        godunov_records(godunov_index)%held_clock_events_closed=.true.
        call close_godunov_clock(godunov_records(godunov_index),physical_flux,scalar_flux,accepted)
        if(.not.accepted) return
      endif
      if(allocated(funding_context_work).and..not.replayed_closure) then
        call close_eb_held_funding_clock_faces_3d(tree,geometry,nv,funding_local,level, &
          physical_flux,scalar_flux,funding_context_work,funding_next,accepted)
        call agree(accepted,comm,local_ok);if(.not.local_ok) return
        call move_eb_held_funding_context_3d(funding_next,funding_context_work)
      endif
      good=.true.
      return
880   continue
      call repair_godunov_parent(godunov_index,godunov_cursor,level,flux,physical_flux,scalar_flux,accepted)
      good=accepted
    end subroutine advance_level
    subroutine finish_deferred_closures(depth,good)
      integer, intent(in) :: depth
      logical, intent(out) :: good
      integer :: lev,r
      logical :: local_ok,valid
      good=.false.
      do lev=levels-1,depth+1,-1
        r=deferred_record(lev)
        if(r==0) cycle
        if(present(post_level)) then
          call post_level(lev,history(r)%time,history(r)%interval,deferred_iteration(lev), &
            deferred_cycle(lev),comm,local_ok)
          call agree(local_ok,comm,valid)
          if(.not.valid) return
        endif
        deferred_record(lev)=0
      enddo
      good=.true.
    end subroutine finish_deferred_closures

    subroutine apply_routed_held_source(gas,temp,stored,prefix,increment,good,scalar_change)
      type(amr_tree_state_3d), intent(inout) :: gas,temp,stored
      real(dp), optional, intent(inout) :: prefix
      type(amr_tree_state_3d), intent(in) :: increment
      logical, intent(out) :: good
      real(dp), optional, intent(inout) :: scalar_change
      real(dp) :: source_reset_change
      type(amr_tree_state_3d) :: difference
      real(dp) :: change
      logical :: valid
      source_reset_change=0
      call apply_eb_held_auxiliary_increment_3d(gas,stored,increment,valid)
      call agree(valid,comm,good);if(.not.good) return
      if(present(prefix)) then
      call reset_held_mpi_eb_tree_family_3d(species,tree,owners,geometry,gas,temp,stored,difference, &
        held_policy,-1,change,comm,good,frame,leaf_only=.true., &
          internal_energy_change=source_reset_change)
      if(.not.good) return
      prefix=prefix+change
      call agree(ieee_is_finite(prefix),comm,good);if(.not.good) return
      else
        call recover_routed_replay(gas,temp,good,stored)
        if(.not.good) return
      endif
      call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,gas,temp,comm,good,frame,stored)
      if(.not.good) return
      if(present(scalar_change)) scalar_change=source_reset_change
    end subroutine apply_routed_held_source

    subroutine recover_routed_replay(gas,temp,good,stored)
      type(amr_tree_state_3d), intent(inout) :: gas,temp
      type(amr_tree_state_3d), optional, intent(in) :: stored
      logical, intent(out) :: good
      real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed
      integer :: p,x,y,z,d(3)
      logical :: accepted,local_ok
      local_ok=.true.
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          if(any(gas%fields(p)%values(6:,x,y,z)<0)) then
            local_ok=.false.;cycle
          endif
          if(present(stored)) then
          call recover_reactive_eb_cell_3d(species,geometry(p),[x,y,z],gas%fields(p)%values(:,x,y,z), &
            temp%fields(p)%values(1,x,y,z),primitive,recovered,speed,accepted,frame,stored%fields(p)%values(1,x,y,z))
          else
          call recover_reactive_eb_cell_3d(species,geometry(p),[x,y,z],gas%fields(p)%values(:,x,y,z), &
            temp%fields(p)%values(1,x,y,z),primitive,recovered,speed,accepted,frame)
          endif
          if(.not.accepted) then
            local_ok=.false.;cycle
          endif
          temp%fields(p)%values(1,x,y,z)=recovered
        enddo
        enddo
        enddo
      enddo
      call agree(local_ok,comm,good)
    end subroutine recover_routed_replay

    subroutine recover_transfer(good)
      logical, intent(out) :: good
      real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed
      integer :: p,x,y,z,d(3)
      logical :: accepted,local_ok
      local_ok=.true.
      do p=1,size(tree%nodes)
        if(owners(p)/=rank) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
          if(any(candidate%fields(p)%values(6:,x,y,z)<0)) then
            local_ok=.false.;cycle
          end if
          call recover_reactive_eb_cell_3d(species,geometry(p),[x,y,z],candidate%fields(p)%values(:,x,y,z), &
            ct%fields(p)%values(1,x,y,z),primitive,recovered,speed,accepted,frame)
          if(.not.accepted) then
            local_ok=.false.;cycle
          end if
          ct%fields(p)%values(1,x,y,z)=recovered
        end do
        end do
        end do
      end do
      call agree(local_ok,comm,good)
    end subroutine recover_transfer
  end subroutine advance_mpi_eb_tree_sdc_subcycle_3d
end module mpi_eb_tree_full_3d_mod
